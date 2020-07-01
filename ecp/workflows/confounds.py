from nipype import Workflow

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from nipype.algorithms import confounds as nac

from niworkflows.interfaces.bids import DerivativesDataSink
from niworkflows.interfaces.confounds import ExpandModel, SpikeRegressors
from niworkflows.interfaces.images import SignalExtraction
from niworkflows.interfaces.masks import ROIsPlot
from niworkflows.interfaces.patches import (
    RobustACompCor as ACompCor,
    RobustTCompCor as TCompCor,
)
from niworkflows.interfaces.plotting import (
    CompCorVariancePlot, ConfoundsCorrelationPlot
)
from niworkflows.interfaces.utils import (
    AddTSVHeader, TSV2JSON, DictMerge
)

from fmriprep.interfaces import GatherConfounds

DEFAULT_MEMORY_MIN_GB = 0.01

def init_bold_confs_wf(
    mem_gb,
    metadata,
    regressors_all_comps,
    regressors_dvars_th,
    regressors_fd_th,
    name="bold_confs_wf",
):
    """
    This workflow calculates confounds for a BOLD series, and aggregates them
    into a :abbr:`TSV (tab-separated value)` file, for use as nuisance
    regressors in a :abbr:`GLM (general linear model)`.

    The following confounds are calculated, with column headings in parentheses:

    #. Region-wise average signal (``csf``, ``white_matter``, ``global_signal``)
    #. DVARS - original and standardized variants (``dvars``, ``std_dvars``)
    #. Framewise displacement, based on head-motion parameters
       (``framewise_displacement``)
    #. Temporal CompCor (``t_comp_cor_XX``)
    #. Anatomical CompCor (``a_comp_cor_XX``)
    #. Cosine basis set for high-pass filtering w/ 0.008 Hz cut-off
       (``cosine_XX``)
    #. Non-steady-state volumes (``non_steady_state_XX``)
    #. Estimated head-motion parameters, in mm and rad
       (``trans_x``, ``trans_y``, ``trans_z``, ``rot_x``, ``rot_y``, ``rot_z``)


    Prior to estimating aCompCor and tCompCor, non-steady-state volumes are
    censored and high-pass filtered using a :abbr:`DCT (discrete cosine
    transform)` basis.
    The cosine basis, as well as one regressor per censored volume, are included
    for convenience.

    .. workflow::
        :graph2use: orig
        :simple_form: yes

        from fmriprep.workflows.bold.confounds import init_bold_confs_wf
        wf = init_bold_confs_wf(
            mem_gb=1,
            metadata={},
            regressors_all_comps=False,
            regressors_dvars_th=1.5,
            regressors_fd_th=0.5,
        )

    **Parameters**

        mem_gb : float
            Size of BOLD file in GB - please note that this size
            should be calculated after resamplings that may extend
            the FoV
        metadata : dict
            BIDS metadata for BOLD file
        name : str
            Name of workflow (default: ``bold_confs_wf``)
        regressors_all_comps: bool
            Indicates whether CompCor decompositions should return all
            components instead of the minimal number of components necessary
            to explain 50 percent of the variance in the decomposition mask.
        regressors_dvars_th
            Criterion for flagging DVARS outliers
        regressors_fd_th
            Criterion for flagging framewise displacement outliers


    **Inputs**

        bold
            BOLD image, after the prescribed corrections (STC, HMC and SDC)
            when available.
        bold_mask
            BOLD series mask
        movpar_file
            SPM-formatted motion parameters file
        skip_vols
            number of non steady state volumes
        csf_mask
            csk mask in MNI 2mm space
        wm_mask
            wm mask in MNI 2mm space
        cortical_gm_mask
            gm mask in MNI 2mm space
        

    **Outputs**

        confounds_file
            TSV of all aggregated confounds
        confounds_metadata
            Confounds metadata dictionary.

    """
    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold', 'bold_mask', 'movpar_file', 'skip_vols', 'csf_mask', 
                'wm_mask', 'cortical_gm_mask']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['confounds_file', 'confounds_metadata']),
        name='outputnode')

    # create tcc mask: fslmaths cortical_gm_mask -dilD -mul -1 -add bold_mask -bin
    tcc_roi = pe.Node(fsl.utils.ImageMaths(op_string='-dilD -mul -1 -add',
                                           args='-bin'),
                       name='tcc_roi')

    # create acc mask fslmaths wm_mask -add csf_mask
    acc_roi = pe.Node(fsl.utils.ImageMaths(op_string='-add'),
                      name='acc_roi')

    # Ensure ROIs don't go off-limits (reduced FoV)
    csf_msk = pe.Node(niu.Function(function=_maskroi), name='csf_msk')
    wm_msk = pe.Node(niu.Function(function=_maskroi), name='wm_msk')
    acc_msk = pe.Node(niu.Function(function=_maskroi), name='acc_msk')
    tcc_msk = pe.Node(niu.Function(function=_maskroi), name='tcc_msk')

    # DVARS
    dvars = pe.Node(nac.ComputeDVARS(save_nstd=True, save_std=True, remove_zerovariance=True),
                    name="dvars", mem_gb=mem_gb)

    # Frame displacement
    fdisp = pe.Node(nac.FramewiseDisplacement(parameter_source="SPM"),
                    name="fdisp", mem_gb=mem_gb)

    # a/t-Compcor
    mrg_lbl_cc = pe.Node(niu.Merge(1), name='merge_rois_cc', run_without_submitting=True)

    tcompcor = pe.Node(
        TCompCor(components_file='tcompcor.tsv', header_prefix='t_comp_cor_', pre_filter='cosine',
                 save_pre_filter=True, save_metadata=True, percentile_threshold=.05,
                 failure_mode='NaN'),
        name="tcompcor", mem_gb=mem_gb)

    acompcor = pe.Node(
        ACompCor(components_file='acompcor.tsv', header_prefix='a_comp_cor_', pre_filter='cosine',
                 save_pre_filter=True, save_metadata=True, mask_names=['combined', 'CSF', 'WM'],
                 merge_method='none', failure_mode='NaN'),
        name="acompcor", mem_gb=mem_gb)

    # Set number of components
    if regressors_all_comps:
        acompcor.inputs.num_components = 'all'
        tcompcor.inputs.num_components = 'all'
    else:
        acompcor.inputs.variance_threshold = 0.5
        tcompcor.inputs.variance_threshold = 0.5

    # Set TR if present
    if 'RepetitionTime' in metadata:
        tcompcor.inputs.repetition_time = metadata['RepetitionTime']
        acompcor.inputs.repetition_time = metadata['RepetitionTime']

    # Global and segment regressors
    mrg_lbl = pe.Node(niu.Merge(3), name='merge_rois', run_without_submitting=True)
    signals = pe.Node(SignalExtraction(class_labels=["csf", "white_matter", "global_signal"]),
                      name="signals", mem_gb=mem_gb)

    # Arrange confounds
    add_dvars_header = pe.Node(
        AddTSVHeader(columns=["dvars"]),
        name="add_dvars_header", mem_gb=0.01, run_without_submitting=True)
    add_std_dvars_header = pe.Node(
        AddTSVHeader(columns=["std_dvars"]),
        name="add_std_dvars_header", mem_gb=0.01, run_without_submitting=True)
    add_motion_headers = pe.Node(
        AddTSVHeader(columns=["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]),
        name="add_motion_headers", mem_gb=0.01, run_without_submitting=True)
    concat = pe.Node(GatherConfounds(), name="concat", mem_gb=0.01, run_without_submitting=True)

    # CompCor metadata
    tcc_metadata_fmt = pe.Node(
        TSV2JSON(index_column='component', drop_columns=['mask'], output=None,
                 additional_metadata={'Method': 'tCompCor'}, enforce_case=True),
        name='tcc_metadata_fmt')
    acc_metadata_fmt = pe.Node(
        TSV2JSON(index_column='component', output=None,
                 additional_metadata={'Method': 'aCompCor'}, enforce_case=True),
        name='acc_metadata_fmt')
    mrg_conf_metadata = pe.Node(niu.Merge(2), name='merge_confound_metadata',
                                run_without_submitting=True)
    mrg_conf_metadata2 = pe.Node(DictMerge(), name='merge_confound_metadata2',
                                 run_without_submitting=True)

    # Expand model to include derivatives and quadratics
    model_expand = pe.Node(ExpandModel(
        model_formula='(dd1(rps + wm + csf + gsr))^^2 + others'),
        name='model_expansion')

    # Add spike regressors
    spike_regress = pe.Node(SpikeRegressors(
        fd_thresh=regressors_fd_th,
        dvars_thresh=regressors_dvars_th),
        name='spike_regressors')

    # Generate reportlet (ROIs)
    mrg_compcor = pe.Node(niu.Merge(2), name='merge_compcor', run_without_submitting=True)
    rois_plot = pe.Node(ROIsPlot(colors=['b', 'magenta'], generate_report=True),
                        name='rois_plot', mem_gb=mem_gb)

    ds_report_bold_rois = pe.Node(
        DerivativesDataSink(desc='rois', keep_dtype=True),
        name='ds_report_bold_rois', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)

    # Generate reportlet (CompCor)
    mrg_cc_metadata = pe.Node(niu.Merge(2), name='merge_compcor_metadata',
                              run_without_submitting=True)
    compcor_plot = pe.Node(
        CompCorVariancePlot(metadata_sources=['tCompCor', 'aCompCor']),
        name='compcor_plot')
    ds_report_compcor = pe.Node(
        DerivativesDataSink(desc='compcorvar', keep_dtype=True),
        name='ds_report_compcor', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)

    # Generate reportlet (Confound correlation)
    conf_corr_plot = pe.Node(
        ConfoundsCorrelationPlot(reference_column='global_signal', max_dim=70),
        name='conf_corr_plot')
    ds_report_conf_corr = pe.Node(
        DerivativesDataSink(desc='confoundcorr', keep_dtype=True),
        name='ds_report_conf_corr', run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB)

    workflow.connect([
        # generate tcc and acc rois
        (inputnode, tcc_roi, [('cortical_gm_mask', 'in_file'),
                              ('bold_mask', 'in_file2')]),
        (inputnode, acc_roi, [('wm_mask', 'in_file'),
                              ('csf_mask', 'in_file2')]),
        # Mask ROIs with bold_mask
        (inputnode, csf_msk, [('bold_mask', 'in_mask')]),
        (inputnode, wm_msk, [('bold_mask', 'in_mask')]),
        (inputnode, acc_msk, [('bold_mask', 'in_mask')]),
        (inputnode, tcc_msk, [('bold_mask', 'in_mask')]),
        # connect inputnode to each non-anatomical confound node
        (inputnode, dvars, [('bold', 'in_file'),
                            ('bold_mask', 'in_mask')]),
        (inputnode, fdisp, [('movpar_file', 'in_file')]),

        # tCompCor
        (inputnode, tcompcor, [('bold', 'realigned_file')]),
        (inputnode, tcompcor, [('skip_vols', 'ignore_initial_volumes')]),
        (tcc_roi, tcc_msk, [('out_file', 'roi_file')]),
        (tcc_msk, tcompcor, [('out', 'mask_files')]),

        # aCompCor
        (inputnode, acompcor, [('bold', 'realigned_file')]),
        (inputnode, acompcor, [('skip_vols', 'ignore_initial_volumes')]),
        (acc_roi, acc_msk, [('out_file', 'roi_file')]),
        (acc_msk, mrg_lbl_cc, [('out', 'in1')]),
        (inputnode, mrg_lbl_cc, [('csf_mask', 'in2')]),
        (inputnode, mrg_lbl_cc, [('wm_mask', 'in3')]),
        (mrg_lbl_cc, acompcor, [('out', 'mask_files')]),

        # Global signals extraction (constrained by anatomy)
        (inputnode, signals, [('bold', 'in_file')]),
        (inputnode, csf_msk, [('csf_mask', 'roi_file')]),
        (csf_msk, mrg_lbl, [('out', 'in1')]),
        (inputnode, wm_msk, [('wm_mask', 'roi_file')]),
        (wm_msk, mrg_lbl, [('out', 'in2')]),
        (inputnode, mrg_lbl, [('bold_mask', 'in3')]),
        (mrg_lbl, signals, [('out', 'label_files')]),

        # Collate computed confounds together
        (inputnode, add_motion_headers, [('movpar_file', 'in_file')]),
        (dvars, add_dvars_header, [('out_nstd', 'in_file')]),
        (dvars, add_std_dvars_header, [('out_std', 'in_file')]),
        (signals, concat, [('out_file', 'signals')]),
        (fdisp, concat, [('out_file', 'fd')]),
        (tcompcor, concat, [('components_file', 'tcompcor'),
                            ('pre_filter_file', 'cos_basis')]),
        (acompcor, concat, [('components_file', 'acompcor')]),
        (add_motion_headers, concat, [('out_file', 'motion')]),
        (add_dvars_header, concat, [('out_file', 'dvars')]),
        (add_std_dvars_header, concat, [('out_file', 'std_dvars')]),

        # Confounds metadata
        (tcompcor, tcc_metadata_fmt, [('metadata_file', 'in_file')]),
        (acompcor, acc_metadata_fmt, [('metadata_file', 'in_file')]),
        (tcc_metadata_fmt, mrg_conf_metadata, [('output', 'in1')]),
        (acc_metadata_fmt, mrg_conf_metadata, [('output', 'in2')]),
        (mrg_conf_metadata, mrg_conf_metadata2, [('out', 'in_dicts')]),

        # Expand the model with derivatives, quadratics, and spikes
        (concat, model_expand, [('confounds_file', 'confounds_file')]),
        (model_expand, spike_regress, [('confounds_file', 'confounds_file')]),

        # Set outputs
        (spike_regress, outputnode, [('confounds_file', 'confounds_file')]),
        (mrg_conf_metadata2, outputnode, [('out_dict', 'confounds_metadata')]),
        (inputnode, rois_plot, [('bold', 'in_file'),
                                ('bold_mask', 'in_mask')]),
        (tcompcor, mrg_compcor, [('high_variance_masks', 'in1')]),
        (acc_msk, mrg_compcor, [('out', 'in2')]),
        (mrg_compcor, rois_plot, [('out', 'in_rois')]),
        (rois_plot, ds_report_bold_rois, [('out_report', 'in_file')]),
        (tcompcor, mrg_cc_metadata, [('metadata_file', 'in1')]),
        (acompcor, mrg_cc_metadata, [('metadata_file', 'in2')]),
        (mrg_cc_metadata, compcor_plot, [('out', 'metadata_files')]),
        (compcor_plot, ds_report_compcor, [('out_file', 'in_file')]),
        (concat, conf_corr_plot, [('confounds_file', 'confounds_file')]),
        (conf_corr_plot, ds_report_conf_corr, [('out_file', 'in_file')]),
    ])

    return workflow

def _maskroi(in_mask, roi_file):
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    roi = nb.load(roi_file)
    roidata = roi.get_data().astype(np.uint8)
    msk = nb.load(in_mask).get_data().astype(bool)
    roidata[~msk] = 0
    roi.set_data_dtype(np.uint8)

    out = fname_presuffix(roi_file, suffix='_boldmsk')
    roi.__class__(roidata, roi.affine, roi.header).to_filename(out)
    return out
