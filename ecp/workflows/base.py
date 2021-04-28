import os

from niworkflows.interfaces import bids
from nipype import Workflow, Node, IdentityInterface, MapNode, Function

from nipype.interfaces.nipy.preprocess import Trim
from nipype.algorithms.confounds import NonSteadyStateDetector

from ecp import utils
from ecp.interfaces.paths import (
    PostFreeSurferFiles, HcpTaskVolumeFiles, HcpTaskCiftiFiles
)
from ecp.interfaces.afni import TProject
from ecp.interfaces.confounds import GetHcpMovement
from ecp.interfaces.workbench import (
    CiftiConvertToNifti, 
    CiftiConvertFromNifti, 
    CiftiParcellate, 
    CiftiCorrelation,
    CiftiMerge
)

from .anat import init_hcp_segment_anat_wf
from .confounds import init_bold_confs_wf, init_timeseries_wf
from .bold import init_generate_boldmask

def init_cleanprep_wf(
    data_dir,
    work_dir, 
    out_dir,
    subject,
    tasks,
    skip_begin,
    skip_end,
    skip_vols=None,
    name='cleanprep_wf',
):
    """
    Creates cleanprep workflow.

    Parameters
    ----------

    data_dir: str
        data directory holding subject folders
    work_dir: str
        working directory
    out_dir: str
        out directory. Final out directory is out_dir/cleanprep
    subject: str
        subject name; data_dir expected to have a subject directory
    tasks: list (str)
        the task names in HCP format
    skip_begin: list (int)
        this option is intended for pre-upgrade scans. It removes this many 
        volumes from the beginning of the HCP movement regressors file, because
        this number of volumes were removed from the scan after preprocessing
        (nifti and cifti match) but not the movement regresssors file.
    skip_end: list (int)
        this option is intended for pre-upgrade scans. It removes this many 
        volumes from the end of the HCP movement regressors file, because
        this number of volumes were removed from the scan after preprocessing
        (nifti and cifti match) but not the movement regresssors file.
    skip_vols: int
        count these beginning volumes as non steady state instead of relying
        on the algorithm

    Outputs
    -------

    cleanprep_wf: Workflow
        the cleanprep workflow

    
    """

    DerivativesDataSink = bids.DerivativesDataSink
    DerivativesDataSink.out_path_base = 'cleanprep'

    anat_name_template = os.path.join('anat', f'sub-{subject}_T1.nii.gz')

    reportlets_dir = os.path.join(work_dir, 'reportlets')

    cleanprep_wf = Workflow(name=name, base_dir=work_dir)
    
    anat_files = PostFreeSurferFiles(base_dir=data_dir,
                                     subject=subject).run().outputs

    hcp_segment_anat_wf = init_hcp_segment_anat_wf()
    inputnode = hcp_segment_anat_wf.inputs.inputnode
    inputnode.brainmask_fs = anat_files.brainmask_fs
    inputnode.l_atlasroi = anat_files.L_atlasroi_32k_fs_LR
    inputnode.l_midthickness = anat_files.L_midthickness_32k_fs_LR
    inputnode.l_white = anat_files.L_white_32k_fs_LR
    inputnode.l_pial = anat_files.L_pial_32k_fs_LR
    inputnode.r_atlasroi = anat_files.R_atlasroi_32k_fs_LR
    inputnode.r_midthickness = anat_files.R_midthickness_32k_fs_LR
    inputnode.r_white = anat_files.R_white_32k_fs_LR
    inputnode.r_pial = anat_files.R_pial_32k_fs_LR
    inputnode.wmparc = anat_files.wmparc
    inputnode.ROIs = anat_files.subcortical

    ds_csf_mask = Node(DerivativesDataSink(
        base_directory=out_dir, desc='csf', source_file=anat_name_template,
        space='mni', suffix='mask'), name='ds_csf_mask')
    ds_wm_mask = Node(DerivativesDataSink(
        base_directory=out_dir, desc='wm', source_file=anat_name_template,
        space='mni', suffix='mask'), name='ds_wm_mask')
    ds_cortical_gm_mask = Node(DerivativesDataSink(
        base_directory=out_dir, desc='cortgm', source_file=anat_name_template,
        space='mni', suffix='mask'), name='ds_coritical_gm_mask')

    cleanprep_wf.connect([
        (hcp_segment_anat_wf, ds_csf_mask, [('outputnode.csf_mask', 'in_file')]),
        (hcp_segment_anat_wf, ds_wm_mask, [('outputnode.wm_mask', 'in_file')]),
        (hcp_segment_anat_wf, ds_cortical_gm_mask, [('outputnode.cort_gm_mask', 'in_file')]),
    ])

    for task, task_skip_begin, task_skip_end in zip(tasks, 
                                                    skip_begin, 
                                                    skip_end):

        task_vol_files = HcpTaskVolumeFiles(
            mninonlinear=anat_files.mninonlinear,
            subject=subject, 
            task=task).run().outputs

        task_cifti_files = HcpTaskCiftiFiles(
            mninonlinear=anat_files.mninonlinear,
            subject=subject, 
            task=task).run().outputs

        func_name_template = os.path.join('func', utils.hcp_to_bids(task, subject))
        task_wf = Workflow(name=task + '_wf')

        input_node = Node(IdentityInterface(
            fields=['csf_mask', 'wm_mask', 'cortical_gm_mask', 'cifti']),
            name='inputnode')
        input_node.inputs.cifti = task_cifti_files.preproc

        movement = Node(
            GetHcpMovement(hcp_movement=task_vol_files.movement_regressors,
                           skip_begin=int(task_skip_begin),
                           skip_end=int(task_skip_end)),
            name='movement')
        
        bold_confs_wf = init_bold_confs_wf(
            mem_gb=1,
            metadata={},
            regressors_all_comps=False,
            regressors_dvars_th=1.5,
            regressors_fd_th=0.5)
        inputnode = bold_confs_wf.inputs.inputnode
        inputnode.bold = task_vol_files.preproc

        if skip_vols:
            inputnode.skip_vols = skip_vols
        else:
            calc_dummy_scans = Node(NonSteadyStateDetector(
                in_file=task_vol_files.preproc), name='calc_dummy_scans') 

        # bold mask
        generate_boldmask = init_generate_boldmask(task_vol_files.preproc)
        ds_boldmask = Node(DerivativesDataSink(
            base_directory=out_dir, 
            desc='confounds',
            source_file=func_name_template, 
            suffix='boldmask'),
            name='ds_boldmask')

        ds_confounds = Node(DerivativesDataSink(
            base_directory=out_dir, 
            desc='confounds', 
            source_file=func_name_template, 
            suffix='regressors'), 
            name='ds_confounds')

        ds_cifti = Node(DerivativesDataSink(
            base_directory=out_dir, 
            space='fsLR32k',
            source_file=func_name_template, 
            suffix='bold.dtseries'), name='ds_fakenifti')

        task_wf.connect([
            (input_node, bold_confs_wf, [('csf_mask', 'inputnode.csf_mask'),
                                         ('wm_mask', 'inputnode.wm_mask'),
                                         ('cortical_gm_mask', 'inputnode.cortical_gm_mask')]),
            (movement, bold_confs_wf, [('movement', 'inputnode.movpar_file')]),
            (generate_boldmask, bold_confs_wf, [('outputnode.bold_mask', 'inputnode.bold_mask')]),
            # derivatives
            (input_node, ds_cifti, [('cifti', 'in_file')]),
            (bold_confs_wf, ds_confounds, [('outputnode.confounds_file', 'in_file'),
                                           ('outputnode.confounds_metadata', 'meta_dict')]),
            (generate_boldmask, ds_boldmask, [('outputnode.bold_mask', 'in_file')]),
        ])
        if skip_vols is None:
            task_wf.connect([
                (calc_dummy_scans, bold_confs_wf,
                 [('n_volumes_to_discard', 'inputnode.skip_vols')]),
            ])

        cleanprep_wf.connect([
            (hcp_segment_anat_wf, task_wf,
             [('outputnode.csf_mask', 'inputnode.csf_mask'),
              ('outputnode.wm_mask', 'inputnode.wm_mask'),
              ('outputnode.cort_gm_mask', 'inputnode.cortical_gm_mask')]),
        ])

        for node in task_wf.list_node_names():
            if node.split('.')[-1].startswith('ds_report'):
                task_wf.get_node(node).inputs.base_directory = reportlets_dir
                task_wf.get_node(node).inputs.source_file = func_name_template

    return cleanprep_wf

def init_clean_wf(
    in_files,
    parcellation,
    out_dir,
    out_path_base,
    space,
    parcellation_space,
    source_file,
    work_dir,
    save_clean_dtseries=False,
    save_clean_ptseries=False,
    save_clean_pconn=False,
    save_clean_covariance=False,
    polort=-1,
    passband=None,
    stopband=None,
    name='clean_wf'
):
    """
    Generates a parcellated connectome.

    Example:
        wf = init_clean_wf(
            in_files = {'cifti': ['sub-001_task-rest_run-1_space-fsLR32k_bold.nii.gz',
                                  'sub-001_task-rest_run-2_space-fsLR32k_bold.nii.gz'],
                        'source_files: ['sub-001_task-rest_run-1_bold.nii.gz',
                                        'sub-001_tesk-rest_run-2_bold.nii.gz'],
                        'censor': ['sub-001_task-rest_run-1_censor.1D',
                                   'sub-001_task-rest_run-2_censor.1D'],
                        'ort': ['sub-001_task-rest_run-1_ort.1D',
                                'sub-001_task-rest_run-2_ort.1D'],
                        'dt': [0.8, 0.8],
                        'trim': [2, 3]},
            parcellation = 'glasser_parcellation.dlabel.nii',
            out_dir = '/path/to/out_dir',
            out_path_base = 'some_folder_name',
            space = 'fsLR32k',
            parcellation_space = 'glasser',
            source_file = 'sub-EC1008_task-rest_bold.nii.gz',
            work_dir = '/path/to/work_dir',
            save_clena_dtseries = False,
            save_clean_ptseries = False,
            save_clean_pconn = False,
            save_clean_covariance = False,
            polort=2,
            passband=[0.01, 0.1],
            stopband=None,
            name='really_cool_name'
        )

    Parameters
    ----------

    in_files : dict
        keys - cifti, source_files, censor, ort, dt
            cifti is a list of preprocessed cifti files
            source_files are used for naming output modified bold files
            censor is a list of afni censor 1D files or None
            ort is a list of ort 1D files or None
            dt is a list of floats for the repetition times
            trim removes these number of volumes from the begining before any
                preprocessing. Number of volumes after trimming must match
                censor and ort volumes. This number is zero-based. Use 0
                if you do not want to trim any beginning frames.
    parcellation : str
        path to cifti parcellation
    out_dir : str
        the output directory
    out_path_base : str
        the new directory for the output, to be crteate within out_dir
    space : str
        the space for the cift files
    parcellation_space : str    
        the space alias for the parcellation, so keep it terse
    source_file : str
        a filename for output naming purposes
    work_dir : str
        the working directory for the workflow
    save_clean_dtseries : bool
        save the clean dtseries for each in cifti
    save_clean_ptseries : bool
        save the clean ptseries for each in cifti
    save_clean_pconn : bool
        save the pconn for each in cifti
    save_clean_covariance : bool
        save the covariance for each in cifti
    polort : int
        remove polynomials up to and including degree polynomial
    passband : list
        fbot, ftop
    stopband : list
        sbot, stop
    name : str
        workflow name

    Returns
    -------

    workflow : nipype workflow
        Here are the workflow steps:
            1. convert cifti to nifti.
            2. Remove nuisance regressors with 3dTproject.
            3. Convert 'cleaned' nifti back to cifti.
            4. Parcellate cifti.
                a) Calculate task connectome.
            5. Merge all task ciftis.
            6. Calculate subject connectome.
    """

    cifti = in_files['cifti']
    source_files = in_files['source_files']
    censor = in_files['censor']
    ort = in_files['ort']
    dt = in_files['dt']
    trim = in_files['trim']

    DerivativesDataSink = bids.DerivativesDataSink
    DerivativesDataSink.out_path_base = out_path_base

    write_verbose = (ort is not None 
                     or polort > 0 
                     or passband is not None
                     or stopband is not None)

    run_tproject = (censor is not None
                    or write_verbose)

    tproject_iterfields = ['in_file', 'TR']
    if censor is not None:
        tproject_iterfields.append('censor')
    if ort is not None:
        tproject_iterfields.append('ort')

    # start workflow now
    clean_wf = Workflow(name=name, base_dir=work_dir)

    cifti_to_nifti = MapNode(CiftiConvertToNifti(
        out_file='fakenifti.nii.gz'),
        name='cifti_to_nifti', iterfield=['in_file'])
    cifti_to_nifti.inputs.in_file = cifti

    trim_begin = MapNode(Trim(),
        name='trim_begin', iterfield=['in_file', 'begin_index'])
    trim_begin.inputs.begin_index = trim

    tproject = MapNode(TProject(
        out_file='clean.nii.gz',
        polort=polort,
        verb=True), name='tproject', iterfield=tproject_iterfields)
    tproject.inputs.TR = dt
    if 'censor' in tproject_iterfields:
        tproject.inputs.censor = censor
    if 'ort' in tproject_iterfields:
        tproject.inputs.ort = ort

    nifti_to_cifti = MapNode(CiftiConvertFromNifti(
        out_file='converted.dtseries.nii'),
        name='nifti_to_cifti', 
        iterfield=['in_file', 'cifti_template', 'reset_timepoints'])
    nifti_to_cifti.inputs.cifti_template = cifti
    nifti_to_cifti.inputs.reset_timepoints = [(x, 0) for x in dt]

    cifti_parcellate = MapNode(CiftiParcellate(
        cifti_label=parcellation,
        direction='COLUMN',
        out_file='parcellated.ptseries.nii'),
        name='cifti_parcellate', iterfield=['in_file'])

    task_rvals = MapNode(CiftiCorrelation(
        out_file='task_rvals.pconn.nii'),
        name='task_rvals', iterfield=['in_file'])

    task_zvals = MapNode(CiftiCorrelation(
        out_file='task_zvals.pconn.nii',
        fisher_z=True),
        name='task_zvals', iterfield=['in_file'])

    task_cov = MapNode(CiftiCorrelation(
        out_file='task_cov.pconn.nii',
        covariance=True),
        name='task_cov', iterfield=['in_file'])

    task_merge = Node(CiftiMerge(
        out_file='merge.ptseries.nii'),
        name='task_merge')

    merge_rvals = Node(CiftiCorrelation(
        out_file='merge_rvals.pconn.nii'),
        name='merge_rvals')

    merge_zvals = Node(CiftiCorrelation(
        out_file='merge_zvals.pconn.nii',
        fisher_z=True),
        name='merge_zvals')

    merge_cov = Node(CiftiCorrelation(
        out_file='merge_cov.pconn.nii',
        covariance=True),
        name='merge_cov')

    # derivatives
    ds_clean_dtseries = MapNode(DerivativesDataSink(
        base_directory=out_dir,
        desc='clean',
        space=space,
        suffix='bold.dtseries'),
        iterfield=['in_file', 'source_file'],
        name='ds_clean_dtseries', 
        run_without_submitting=True)
    ds_clean_dtseries.inputs.source_file = source_files

    ds_task_ptseries = MapNode(DerivativesDataSink(
        base_directory=out_dir,
        desc='clean',
        space=parcellation_space,
        suffix='bold.ptseries'),
        iterfield=['in_file', 'source_file'],
        name='ds_task_ptseries', 
        run_without_submitting=True)
    ds_task_ptseries.inputs.source_file = source_files

    ds_task_zvals = MapNode(DerivativesDataSink(
        base_directory=out_dir,
        space=parcellation_space,
        suffix='zvals.pconn'),
        iterfield=['in_file', 'source_file'],
        name='ds_task_zvals', 
        run_without_submitting=True)
    ds_task_zvals.inputs.source_file = source_files
    
    ds_task_rvals = MapNode(DerivativesDataSink(
        base_directory=out_dir,
        space=parcellation_space,
        suffix='rvals.pconn'),
        iterfield=['in_file', 'source_file'],
        name='ds_task_rvals', 
        run_without_submitting=True)
    ds_task_rvals.inputs.source_file = source_files

    ds_task_cov = MapNode(DerivativesDataSink(
        base_directory=out_dir,
        space=parcellation_space,
        suffix='cov.pconn'),
        iterfield=['in_file', 'source_file'],
        name='ds_task_cov', 
        run_without_submitting=True)
    ds_task_cov.inputs.source_file = source_files

    ds_subject_zvals = Node(DerivativesDataSink(
        base_directory=out_dir,
        desc='concatenated',
        space=parcellation_space,
        suffix='zvals.pconn',
        source_file=source_file),
        name='ds_subject_zvals', run_without_submitting=True)

    ds_subject_rvals = Node(DerivativesDataSink(
        base_directory=out_dir,
        desc='concatenated',
        space=parcellation_space,
        suffix='rvals.pconn',
        source_file=source_file),
        name='ds_subject_rvals', run_without_submitting=True)

    ds_subject_cov = Node(DerivativesDataSink(
        base_directory=out_dir,
        desc='concatenated',
        space=parcellation_space,
        suffix='cov.pconn',
        source_file=source_file),
        name='ds_subject_cov', run_without_submitting=True)

    ds_ort = MapNode(DerivativesDataSink(
        base_directory=out_dir,
        desc='confounds',
        suffix='tproject'),
        iterfield=['in_file', 'source_file'],
        name='ds_ort',
        run_without_submittting=True)
    ds_ort.inputs.source_file = source_files

    ds_sval = MapNode(DerivativesDataSink(
        base_directory=out_dir,
        desc='sval',
        suffix='tproject'),
        iterfield=['in_file', 'source_file'],
        name='ds_sval',
        run_without_submitting=True)
    ds_sval.inputs.source_file = source_files

    ds_psinv = MapNode(DerivativesDataSink(
        base_directory=out_dir,
        desc='psinv',
        suffix='tproject'),
        iterfield=['in_file', 'source_file'],
        name='ds_psinv',
        run_without_submitting=True)
    ds_psinv.inputs.source_file = source_files

    if run_tproject:

        clean_wf.connect([
            (cifti_to_nifti, trim_begin, [('out_file', 'in_file')]),
            (trim_begin, tproject, [('out_file', 'in_file')]),
            (tproject, nifti_to_cifti, [('out_file', 'in_file')]),
            (nifti_to_cifti, cifti_parcellate, [('out_file', 'in_file')]),
        ])

        if write_verbose:
            clean_wf.connect([
                (tproject, ds_ort, [('matrix', 'in_file')]),
                (tproject, ds_sval, [('singular_values', 'in_file')]),
                (tproject, ds_psinv, [('pseudo_inv', 'in_file')]),
            ])
    else:

        clean_wf.connect([
            (cifti_to_nifti, trim_begin, [('out_file', 'in_file')]),
            (trim_begin, nifti_to_cifti, [('out_file', 'in_file')]),
            (nifti_to_cifti, cifti_parcellate, [('out_file', 'in_file')]),
        ])
        
    clean_wf.connect([
        (cifti_parcellate, task_merge, [('out_file', 'in_files')]),
        (task_merge, merge_rvals, [('out_file', 'in_file')]),
        (task_merge, merge_zvals, [('out_file', 'in_file')]),
        (task_merge, merge_cov, [('out_file', 'in_file')]),
        # derivatives
        (merge_rvals, ds_subject_rvals, [('out_file', 'in_file')]),
        (merge_zvals, ds_subject_zvals, [('out_file', 'in_file')]),
        (merge_cov, ds_subject_cov, [('out_file', 'in_file')]),
    ])

    if save_clean_dtseries:
        clean_wf.connect([
            (nifti_to_cifti, ds_clean_dtseries, [('out_file', 'in_file')]),
        ])

    if save_clean_ptseries:
        clean_wf.connect([
            (cifti_parcellate, ds_task_ptseries, [('out_file', 'in_file')]),
        ])

    if save_clean_pconn:
        clean_wf.connect([
            (cifti_parcellate, task_rvals, [('out_file', 'in_file')]),
            (cifti_parcellate, task_zvals, [('out_file', 'in_file')]),
            # derivatives
            (task_rvals, ds_task_rvals, [('out_file', 'in_file')]),
            (task_zvals, ds_task_zvals, [('out_file', 'in_file')]),
        ])

    if save_clean_covariance:
        clean_wf.connect([
            (cifti_parcellate, task_cov, [('out_file', 'in_file')]),
            # derivatives
            (task_cov, ds_task_cov, [('out_file', 'in_file')]),
        ])

    return clean_wf

def init_clean_cifti(
    out_dir,
    out_path_base,
    space,
    source_file,
    save_clean_dtseries,
    cifti,
    work_dir=None,
    censor=None,
    ort=None,
    polort=-1,
    passband=None,
    dt=None,
    name='clean_cifti_wf'
):
    """
    Removes 'nuisance' time series from each voxel in the input dataset.

    To skip TProject, use default TProject parameter values. If you want to run
    TProject and not remove anything, use a censor file with all 1's. When this
    is done, no matrix, singular values, or pseudo inverse will be written.

    Parameters
    ----------
    out_dir : str
        the output directory
    out_path_base : str
        the new directory for the output, to be created within out_dir
    space : str
        the space for the bold files
    source_file : str
        a filename for output naming purposes
    save_clean_dtseries : bool
        save clean dtseries
    work_dir : str
        the working directory for the workflow
    cifti : str
        path to cifti file
    censor : str
        path to 1D censor file
    ort : str
        path to 1D nuisance regressor file
    polort : int
        remove polynomials up to degree polort
    passband : list 
        fbot, ftop
    dt : float
        repetition time
    name : str
        the workflow name

    Returns
    -------

    workflow : nipype workflow
        Here are teh workflow steps:
            1. Convert cifti to nifti
            2. Remove nuisance regressors with 3dTproject
            3. Covnert 'cleaned' nifti back to cifti
    """

    DerivativesDataSink = BIDSDerivatives
    DerivativesDataSink.out_path_base = out_path_base

    write_verbose = (ort is not None 
                     and polort > 0 
                     and passband is not None)

    run_tproject = (censor is not None
                    and ort is not None
                    and polort > 0
                    and passband is not None)
                    

    # start workflow now
    workflow = Workflow(name=name, base_dir=work_dir)

    outputnode = Node(IdentityInterface(
        fields=['clean_cifti']), name='outputnode')

    cifti_to_nifti = Node(CiftiConvertToNifti(
        in_file=cifti, 
        out_file='fakenifti.nii.gz'), name='cifti_convert_to_nifti') 

    tproject = Node(TProject(
        verb=True,
        out_file='tproject.nii.gz',
        censor=censor,
        ort=ort,
        polort=polort,
        passband=passband,
        dt=dt,), name='tproject')

    nifti_to_cifti = Node(CiftiConvertFromNifti(
        cifti_template=cifti,
        reset_timepoints=(dt, 0),
        out_file='converted.nii.gz'), name='nifti_to_cifti')

    # derivatives
    ds_clean = Node(DerivativesDataSink(
        source_file=source_file,
        base_directory=out_dir,
        space=space,
        suffix='cleanbold'), name='ds_clean', run_without_submitting=True)

    ds_ort = Node(DerivativesDataSink(
        source_file=source_file,
        base_directory=out_dir,
        space=space,
        suffix='ort'), name='ds_ort', run_without_submitting=True)

    ds_sval = Node(DerivativesDataSink(
        source_file=source_file,
        base_directory=out_dir,
        space=space,
        suffix='sval'), name='ds_sval', run_without_submitting=True)

    ds_psinv = Node(DerivativesDataSink(
        source_file=source_file,
        base_directory=out_dir,
        space=space,
        suffix='psinv'), name='ds_psinv', run_without_submitting=True)

    if run_tproject:
        workflow.connect([
            (cifti_to_nifti, tproject, [('out_file', 'in_file')]),
            (tproject, nifti_to_cifti, [('out_file', 'in_file')]),
            (nifti_to_cifti, outputnode, [('out_file', 'clean_cifti')]),
        ])

        if save_clean_dtseries:
            workflow.connect([
                (tproject, ds_clean, [('out_file', 'in_file')])
            ])

        if write_verbose:
            workflow.connect([
                (tproject, ds_ort, [('matirx', 'in_file')]),
                (tproject, ds_sval, [('singular_values', 'in_file')]),
                (tproject, ds_psinv, [('pseudo_inv', 'in_file')]),
            ])
    else:
        workflow.connect([
            (cifti_to_nifti, nifti_to_cifti, [('in_file', 'out_file')]),
            (nifti_to_cifti, outputnode, [('out_file', 'clean_cifti')]),
        ])

        if save_clean_dtseries:
            workflow.connect([
                (nifti_to_cifti, ds_clean, [('out_file', 'in_file')])
            ])

    return workflow
                
def init_bidsify_hcp_wf(
    data_dir,
    work_dir, 
    out_dir,
    subject,
    tasks,
    skip_begin,
    skip_end,
    name='bidsify_hcp_wf',
):
    """
    Creates bidisfy_hcp workflow.

    Parameters
    ----------

    data_dir: str
        data directory holding subject folders
    work_dir: str
        working directory
    out_dir: str
        out directory. Final out directory is out_dir/bidsify_hcp
    subject: str
        subject name; data_dir expected to have a subject directory
    tasks: list (str)
        the task names in HCP format
    skip_begin: list (int)
        this option is intended for pre-upgrade scans. It removes this many 
        volumes from the beginning of the HCP movement regressors file, because
        this number of volumes were removed from the scan after preprocessing
        (nifti and cifti match) but not the movement regresssors file.
    skip_end: list (int)
        this option is intended for pre-upgrade scans. It removes this many 
        volumes from the end of the HCP movement regressors file, because
        this number of volumes were removed from the scan after preprocessing
        (nifti and cifti match) but not the movement regresssors file.

    Outputs
    -------

    bidsify_hcp_wf: Workflow
        the cleanprep workflow
    """

    DerivativesDataSink = bids.DerivativesDataSink
    DerivativesDataSink.out_path_base = 'bidsify_hcp'

    anat_name_template = os.path.join('anat', f'sub-{subject}_T1.nii.gz')

    bidsify_hcp_wf = Workflow(name=name, base_dir=work_dir)
    
    anat_files = PostFreeSurferFiles(base_dir=data_dir,
                                     subject=subject).run().outputs

    hcp_segment_anat_wf = init_hcp_segment_anat_wf()
    inputnode = hcp_segment_anat_wf.inputs.inputnode
    inputnode.brainmask_fs = anat_files.brainmask_fs
    inputnode.l_atlasroi = anat_files.L_atlasroi_32k_fs_LR
    inputnode.l_midthickness = anat_files.L_midthickness_32k_fs_LR
    inputnode.l_white = anat_files.L_white_32k_fs_LR
    inputnode.l_pial = anat_files.L_pial_32k_fs_LR
    inputnode.r_atlasroi = anat_files.R_atlasroi_32k_fs_LR
    inputnode.r_midthickness = anat_files.R_midthickness_32k_fs_LR
    inputnode.r_white = anat_files.R_white_32k_fs_LR
    inputnode.r_pial = anat_files.R_pial_32k_fs_LR
    inputnode.wmparc = anat_files.wmparc
    inputnode.ROIs = anat_files.subcortical

    ds_csf_mask = Node(DerivativesDataSink(
        base_directory=out_dir, desc='csf', source_file=anat_name_template,
        space='mni', suffix='mask'), name='ds_csf_mask')
    ds_wm_mask = Node(DerivativesDataSink(
        base_directory=out_dir, desc='wm', source_file=anat_name_template,
        space='mni', suffix='mask'), name='ds_wm_mask')
    ds_cortical_gm_mask = Node(DerivativesDataSink(
        base_directory=out_dir, desc='cortgm', source_file=anat_name_template,
        space='mni', suffix='mask'), name='ds_coritical_gm_mask')

    bidsify_hcp_wf.connect([
        (hcp_segment_anat_wf, ds_csf_mask, [('outputnode.csf_mask', 'in_file')]),
        (hcp_segment_anat_wf, ds_wm_mask, [('outputnode.wm_mask', 'in_file')]),
        (hcp_segment_anat_wf, ds_cortical_gm_mask, [('outputnode.cort_gm_mask', 'in_file')]),
    ])

    out_func_dir = os.path.join(out_dir, 
                                DerivativesDataSink.out_path_base,
                                f'sub-{subject}',
                                'func')

    for task, task_skip_begin, task_skip_end in zip(tasks, 
                                                    skip_begin, 
                                                    skip_end):

        out_vol = utils.hcp_to_bids(task, subject)
        entities = utils.get_entities(out_vol)
        out_vol = os.path.join(out_func_dir, out_vol)

        out_cifti = utils.generate_bold_name(
            subject, entities['task'], 'bold', 'dtseries.nii', dir=entities['dir'], 
            run=entities['run'], space='fsLR32k')
        out_cifti = os.path.join(out_func_dir, out_cifti)

        task_vol_files = HcpTaskVolumeFiles(
            mninonlinear=anat_files.mninonlinear,
            subject=subject, 
            task=task).run().outputs

        task_cifti_files = HcpTaskCiftiFiles(
            mninonlinear=anat_files.mninonlinear,
            subject=subject, 
            task=task).run().outputs

        func_name_template = os.path.join('func', utils.hcp_to_bids(task, subject))
        task_wf = Workflow(name=task + '_wf')

        movement = Node(
            GetHcpMovement(hcp_movement=task_vol_files.movement_regressors,
                           skip_begin=int(task_skip_begin),
                           skip_end=int(task_skip_end)),
            name='movement')
        
        link_vol = Node(Function(input_name=['originalfile', 'newfile'],
                                 function=_hardlink), name='link_vol')
        link_vol.inputs.originalfile = task_vol_files.preproc
        link_vol.inputs.newfile = out_vol

        link_cifti = Node(Function(input_name=['originalfile', 'newfile'],
                                   function=_hardlink), name='link_cifti')
        link_cifti.inputs.originalfile = task_cifti_files.preproc
        link_cifti.inputs.newfile = out_cifti

        generate_boldmask = init_generate_boldmask(task_vol_files.preproc)
        ds_boldmask = Node(DerivativesDataSink(
            base_directory=out_dir, 
            desc='confounds',
            source_file=func_name_template, 
            suffix='boldmask'),
            name='ds_boldmask')

        ds_movement = Node(DerivativesDataSink(
            base_directory=out_dir, 
            desc='movpar', 
            source_file=func_name_template, 
            suffix='timeseries'), 
            name='ds_movement')

        task_wf.connect([
            # derivatives
            (generate_boldmask, ds_boldmask, [('outputnode.bold_mask', 'in_file')]),
            (movement, ds_movement, [('movement', 'in_file')]),
        ])
        task_wf.add_nodes([link_vol, link_cifti])

        bidsify_hcp_wf.add_nodes([task_wf])

    return bidsify_hcp_wf

def init_multi_timeseries_wf(
    parameters,
    csf_mask,
    wm_mask,
    cortical_gm_mask,
    work_dir,
    out_dir,
    out_path_base,
    name='multi_timeseries_wf'
):
    """
    Creates nuissasnce timeseries for nifti files.

    Example:
        wf = init_timeseries_wf(
            parmaeters = [
                {
                    'bold': 'sub-EC1008_task-rest_dir-ap_run-1_bold.nii.gz',
                    'bold_mask': 'sub-EC1008_task-rest_dir-ap_run-1_desc-confounds_boldmask.nii.gz',
                    'movpar_file': 'sub-EC1008_task-rest_dir-ap_run-1_desc-movpar_timeseries.tsv',
                    'skip_vols': 1,
                    'wf_name': 'rest_ap_1',
                },
                {
                    'bold': 'sub-EC1008_task-rest_dir-pa_run-1_bold.nii.gz',
                    'bold_mask': 'sub-EC1008_task-rest_dir-pa_run-1_desc-confounds_boldmask.nii.gz',
                    'movpar_file': 'sub-EC1008_task-rest_dir-pa_run-1_desc-movpar_timeseries.tsv',
                    'skip_vols': 3,
                    'wf_name': 'rest_pa_1',
                },
            ]
            csf_mask = 'sub-EC1008_space-mni_desc-csf_mask.nii.gz',
            wm_mask = 'sub-EC1008_space-mni_desc-wm_mask.nii.gz',
            work_dir = '/path/to/work_dir',
            out_dir = '/path/to/out_dir',
            out_path_base = 'no_skip_vols',
            skip_vols = Nne
        )

    Parameters
    ----------

    parameters: list of dicts
        dict keys - bold, bold_mask, movpar_files, skip_vols, wf_name
            bold is a list of preprocessed bold files
            bold_mask is a list of bold masks
            movpar_files is a list of movement parameter files
            skip_vols is a list of skipped volumes
            wf_name a list of workflow names for each bold
    csf_mask: str
        path to subject csf mask
    wm_mask: str
        path to subject white matter mask
    work_dir: str
        path to working directory
    out_dir: str
        path to out directory
    out_path_base : str
        the new directory for the output, to be crteate within out_dir

    Outputs
    -------

    multi_timeseries_wf: Workflow
        the timeseries workflow
    """

    DerivativesDataSink = bids.DerivativesDataSink
    DerivativesDataSink.out_path_base = out_path_base

    multi_timeseries_wf = Workflow(name=name, base_dir=work_dir)

    anat = Node(IdentityInterface(
        fields=['csf_mask', 'wm_mask', 'cortical_gm_mask']), name='anat')
    anat.inputs.csf_mask = csf_mask
    anat.inputs.wm_mask = wm_mask
    anat.inputs.cortical_gm_mask = cortical_gm_mask

    for parameter in parameters:

        timeseries_wf = init_timeseries_wf(
            out_dir, 
            out_path_base, 
            parameter['source_file'],
            parameter['dt'],
            name=parameter['wf_name'])
        timeseries_wf.inputs.inputnode.bold_std = parameter['bold']
        timeseries_wf.inputs.inputnode.bold_mask_std = parameter['bold_mask']
        timeseries_wf.inputs.inputnode.movpar_file = parameter['movpar_file']
        timeseries_wf.inputs.inputnode.skip_vols = parameter['skip_vols']

        multi_timeseries_wf.connect([
            (anat, timeseries_wf, [('csf_mask', 'inputnode.csf_mask'),
                                   ('wm_mask', 'inputnode.wm_mask'),
                                   ('cortical_gm_mask', 'inputnode.cortical_gm_mask')])
        ])

    return multi_timeseries_wf


def _hardlink(originalfile, newfile):
    import os
    from shutil import copyfile

    out_dir = os.path.dirname(newfile)
    os.makedirs(out_dir, exist_ok=True)

    try:
        os.link(originalfile, newfile)
    except PermissionError:
        copyfile(originalfile, newfile)
        
    
