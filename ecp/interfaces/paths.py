from os.path import join as opj

from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec, File, Directory, isdefined,
    SimpleInterface
)

class PostFreeSurferFilesInputSpec(BaseInterfaceInputSpec):
    base_dir = Directory(
        desc='directory containing subject directories',
        mandatory=True,
        exists=True)
    subject = traits.Str(
        desc='subject id',
        mandatory=True)

class PostFreeSurferFilesOutputSpec(TraitedSpec):
    subject_dir = Directory(desc='subject directory', exists=True)
    mninonlinear = Directory(desc='MNINonLinear', exists=True)
    roi_folder = Directory(desc='ROIFolder', exists=True)
    down_sample_folder = Directory(desc='DownSampleFolder', exists=True)

    # volumes
    wmparc = File(desc='wmparc in 2 mm', exists=True)
    subcortical = File(desc='subcortical regions parcellated in 2 mm', 
                       exists=True)

    # surfaces in 32k_fs_LR
    L_atlasroi_32k_fs_LR = File(desc='atlasroi', exists=True)
    L_midthickness_32k_fs_LR = File(desc='midthickness', exists=True)
    L_pial_32k_fs_LR = File(desc='pial', exists=True)
    L_white_32k_fs_LR = File(desc='white', exits=True)

    R_atlasroi_32k_fs_LR = File(desc='atlasroi', exists=True)
    R_midthickness_32k_fs_LR = File(desc='midthickness', exists=True)
    R_pial_32k_fs_LR = File(desc='pial', exists=True)
    R_white_32k_fs_LR = File(desc='white', exits=True)
    

class PostFreeSurferFiles(SimpleInterface):
    input_spec = PostFreeSurferFilesInputSpec
    output_spec = PostFreeSurferFilesOutputSpec

    def _run_interface(self, runtime):
        self._results['subject_dir'] = opj(self.inputs.base_dir, self.inputs.subject)
        self._results['mninonlinear'] = opj(self._results['subject_dir'], 'MNINonLinear')
        self._results['roi_folder'] = opj(self._results['mninonlinear'], 'ROIs')
        self._results['down_sample_folder'] = opj(self._results['mninonlinear'], 'fsaverage_LR32k')

        self._results['wmparc'] = opj(self._results['roi_folder'], 'wmparc.2.nii.gz')
        self._results['subcortical'] = opj(self._results['roi_folder'], 'ROIs.2.nii.gz')

        surface_template = '{subject}.{hemi}.{structure}.32k_fs_LR.{gii_type}.gii'
        self._results['L_atlasroi_32k_fs_LR'] = opj(self._results['down_sample_folder'],
            surface_template.format(subject=self.inputs.subject,
                                    hemi='L',
                                    structure='atlasroi',
                                    gii_type='shape'))
        self._results['L_midthickness_32k_fs_LR'] = opj(self._results['down_sample_folder'],
            surface_template.format(subject=self.inputs.subject,
                                    hemi='L',
                                    structure='midthickness',
                                    gii_type='surf'))
        self._results['L_pial_32k_fs_LR'] = opj(self._results['down_sample_folder'],
            surface_template.format(subject=self.inputs.subject,
                                    hemi='L',
                                    structure='pial',
                                    gii_type='surf'))
        self._results['L_white_32k_fs_LR'] = opj(self._results['down_sample_folder'],
            surface_template.format(subject=self.inputs.subject,
                                    hemi='L',
                                    structure='white',
                                    gii_type='surf'))

        self._results['R_atlasroi_32k_fs_LR'] = opj(self._results['down_sample_folder'],
            surface_template.format(subject=self.inputs.subject,
                                    hemi='R',
                                    structure='atlasroi',
                                    gii_type='shape'))
        self._results['R_midthickness_32k_fs_LR'] = opj(self._results['down_sample_folder'],
            surface_template.format(subject=self.inputs.subject,
                                    hemi='R',
                                    structure='midthickness',
                                    gii_type='surf'))
        self._results['R_pial_32k_fs_LR'] = opj(self._results['down_sample_folder'],
            surface_template.format(subject=self.inputs.subject,
                                    hemi='R',
                                    structure='pial',
                                    gii_type='surf'))
        self._results['R_white_32k_fs_LR'] = opj(self._results['down_sample_folder'],
            surface_template.format(subject=self.inputs.subject,
                                    hemi='R',
                                    structure='white',
                                    gii_type='surf'))
        return runtime

class HcpTaskVolumeFilesInputSpec(BaseInterfaceInputSpec):
    mninonlinear = Directory(
        desc='MNINonLinear directory for subject',
        madatory=True,
        exists=True)
    subject = traits.Str(
        desc='subject id',
        mandatory=True)
    task = traits.Str(
        desc='task name',
        mandatory=True)
    
class HcpTaskVolumeFilesOutputSpec(TraitedSpec):
    results_dir = Directory(desc='Results', exists=True)
    task_dir = Directory(desc='results task directory', exists=True)

    # {trans_x, trans_y, tran_z, rot_x, rot_y, rot_z} trans=mm, rot=degrees
    movement_regressors = File(desc='movement regressors', exists=True)
    movement_regressors_dt = File(desc='derivative of movement regressors', exists=True)

    mask = File(desc='mask for volume in MNI 2mm space', exists=True)
    preproc = File(desc='preprocessed task volume in MNI space', exists=True)

class HcpTaskVolumeFiles(SimpleInterface):
    input_spec = HcpTaskVolumeFilesInputSpec
    output_spec = HcpTaskVolumeFilesOutputSpec

    def _run_interface(self, runtime):
        self._results['results_dir'] = opj(self.inputs.mninonlinear, 'Results')
        self._results['task_dir'] = opj(self._results['results_dir'], self.inputs.task)

        self._results['movement_regressors'] = opj(self._results['task_dir'],
                                                   'Movement_Regressors.txt')
        self._results['movement_regressors_dt'] = opj(self._results['task_dir'],
                                                      'Movement_Regressors_dt.txt')

        self._results['mask'] = opj(self._results['task_dir'],
                                    'brainmask_fs.2.nii.gz')
        self._results['preproc'] = opj(self._results['task_dir'],
                                       self.inputs.task + '.nii.gz')
        return runtime

class HcpTaskCiftiFilesInputSpec(BaseInterfaceSpec):
    mninonlinear = Directory(
        desc='MNINonLinear directory for subject',
        madatory=True,
        exists=True)
    subject = traits.Str(
        desc='subject id',
        mandatory=True)
    task = traits.Str(
        desc='task name',
        mandatory=True)

class HcpTaskCiftiFilesOutputSpec(TraitedSpec):
    results_dir = Directory(desc='Results', exists=True)
    task_dir = Directory(desc='results task directory', exists=True)

    # {trans_x, trans_y, tran_z, rot_x, rot_y, rot_z} trans=mm, rot=degrees
    preproc = File(desc'preprocessed task cifti file in fs_LR32k space', exist=True)

class HcpTaskCiftiFiles(SimpleInterface):
    input_spec = HcpTaskCiftiFilesInputSpec
    output_spec = HcpTaskCiftiFilesOutpuSpec

    def _run_interface(self, runtim):
        self._results['results_dir'] = opj(self.inputs.mninonlinear, 'Results')
        self._results['task_dir'] = opj(self._results['results_dir'], self.inputs.task)

        self._results['preproc'] = opj(self._results['task_dir'],
                                       self.inputs.task + '_Atlas.dtseries.nii')
