import os
import re

from nipype import Workflow, Node, IdentityInterface

from ..interfaces.paths import (
    PostFreeSurferFiles, HcpTaskVolumeFiles, HcpTaskCiftiFiles
)
from ..interfaces.confounds import GetHcpMovement, RegressorsTsvTo1D
from ..interfaces.workbench import CiftiConvertToNifti

from .anat import init_hcp_segment_anat_wf
from .confounds import init_bold_confs_wf

from niworkflows.interfaces import bids

class DerivativesDataSink(bids.DerivativesDataSink):
    out_path_base = 'cleanprep'

def init_cleanprep_wf(
    data_dir,
    work_dir, 
    out_dir,
    subject,
    tasks,
    skipped_vols,
):
    anat_name_template = f'sub-{subject}_T1.nii.gz'

    reportlets_dir = os.path.join(work_dir, 'reportlets')

    cleanprep_wf = Workflow(name='cleanprep_wf')
    cleanprep_wf.base_dir = os.path.join(work_dir, subject)
    
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
        base_directory=out_dir, desc='csf',source_file=anat_name_template,
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

    for task, task_skipped_vols in zip(tasks, skipped_vols):
        func_name_template = hcp_to_bids(task, subject)
        task_wf = Workflow(name=task + '_wf')

        input_node = Node(IdentityInterface(
            fields=['csf_mask', 'wm_mask', 'cortical_gm_mask']),
            name='inputnode')

        task_vol_files = HcpTaskVolumeFiles(
            mninonlinear=anat_files.mninonlinear,
            subject=subject, 
            task=task).run().outputs

        task_cifti_files = HcpTaskCiftiFiles(
            mninonlinear=anat_files.mninonlinear,
            subject=subject, 
            task=task).run().outputs

        movement = Node(
            GetHcpMovement(hcp_movement=task_vol_files.movement_regressors,
                           task_skipped_vols=skipped),
            name='movement')
        
        bold_confs_wf = init_bold_confs_wf(
            mem_gb=1,
            metadata={},
            regressors_all_comps=False,
            regressors_dvars_th=1.5,
            regressors_fd_th=0.5)
        inputnode = bold_confs_wf.inputs.inputnode
        inputnode.bold = task_vol_files.preproc
        inputnode.bold_mask = task_vol_files.mask
        inputnode.skip_vols = skip_vols

        cifti_convert_to_nifti = Node(CiftiConvertToNifti(
            in_file=task_cifti_files.preproc, out_file='fakenifti.nii.gz'), 
            name='cifti_convert_to_nifti') 

        ds_confounds = Node(DerivativesDataSink(
            base_directory=out_dir, desc='confounds', 
            source_file=func_name_template, suffix='regressors'), 
            name='ds_confounds')

        ds_fakenifti = Node(DerivativesDataSink(
            base_directory=out_dir, space='fsLR32k',
            source_file=func_name_template, suffix='fakenifti'),
            name='ds_fakenifti')

        task_wf.add_nodes([cifti_convert_to_nifti])
        task_wf.connect([
            (movement, bold_confs_wf, 
             [('movement', 'inputnode.movpar_file')]),
            (input_node, bold_confs_wf,
             [('csf_mask', 'inputnode.csf_mask'),
              ('wm_mask', 'inputnode.wm_mask'),
              ('cortical_gm_mask', 'inputnode.cortical_gm_mask')]),
             # derivatives
             (bold_confs_wf, ds_confounds,
              [('outputnode.confounds_file', 'in_file'),
               ('outputnode.confounds_metadata', 'meta_dict')]),
             (cifti_convert_to_nifti, ds_fakenifti,
              [('out_file', 'in_file')])
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

def hcp_to_bids(hcp_name, subject, ses=None):
    tokenized = hcp_name.split('_')
    task, run = re.match('([a-zA-Z]+)(\d+)', tokenized[1]).groups()
    task = task.lower()
    run = '{:02d}'.format(int(run))
    direction = tokenized[2].lower()

    if ses:
        template = 'sub-{subject}_ses-{ses}_task-{task}_dir-{direction}_run-{run}_bold.nii.gz'
    else:
        template = 'sub-{subject}_task-{task}_dir-{direction}_run-{run}_bold.nii.gz'

    return template.format(
        subject=subject, ses=ses, task=task, direction=direction, run=run)
    
    
            

            
            
    




                             
                     


        
                                            
                                    
