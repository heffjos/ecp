import os

from nipype.interfaces.base import (
    CommandLine, 
    traits, 
    CommandLineInputSpec, 
    isdefined,
    File, 
    TraitedSpec, 
    PackageInfo,
    InputMultiPath
)
from nipype.utils.filemanip import fname_presuffix

class Info(PackageInfo):
    """Handle workbench version information."""
    version_cmd = 'wb_command -version'

    @staticmethod
    def parse_version(raw_info):
        return raw_info

class WBCommand(CommandLine):
    """Base class for workbench commands."""
    input_spec = CommandLineInputSpec

    def _gen_fname(self,
                   basename,
                   cwd=None,
                   suffix=None,
                   change_ext=True,
                   ext=None):
        """Generate a filename based on the given parameters.

        The filename will take the form: cwd/basename<suffix><ext>.
        If change_ext is True, it will use the extentions specified in
        <instance>intputs.output_type.

        Parameters
        ----------
        basename : str
            Filename to base the new filename on.
        cwd : str
            Path to prefix to the new filename. (default is os.getcwd())
        suffix : str
            Suffix to add to the `basename`.  (defaults is '' )
        change_ext : bool
            Flag to change the filename extension to the FSL output type.
            (default True)

        Returns
        -------
        fname : str
            New filename based on given parameters.

        """

        if basename == '':
            msg = 'Unable to generate filename for command %s. ' % self.cmd
            msg += 'basename is not set!'
            raise ValueError(msg)
        if cwd is None:
            cwd = os.getcwd()
        if ext is None:
            ext = '.nii.gz'
        if change_ext:
            if suffix:
                suffix = ''.join((suffix, ext))
            else:
                suffix = ext
        if suffix is None:
            suffix = ''
        fname = fname_presuffix(
            basename, suffix=suffix, use_ext=False, newpath=cwd)
        return fname

class VolumeLabelImportInputSpec(CommandLineInputSpec):
    in_file = File(
        desc='the input volume file',
        argstr='%s',
        position=0,
        mandatory=True,
        exists=True,
        copyfile=False)
    label_list_file = File(
        desc='the text file containing the values and names for labels',
        argstr='%s',
        position=1,
        mandatory=True,
        exists=True,
        copyfile=False)
    out_file = File(
        desc='the output workbench label volume',
        argstr='%s',
        position=2,
        mandatory=True)
    discard_others = traits.Bool(
        desc='set any voxels not mentioned in the label list to the ??? label',
        argstr='-discard-others')
    unlabeled_value = traits.Int(
        desc='set the value that will be interpreted as unlabled (default is 0)',
        argstr='-datum %d')
    subvolume = traits.Str(
        desc='select a single subvolume to import by number or name',
        argstr='-subvolume %s')
    drop_unused_labels = traits.Bool(
        desc='remove any unused label values from the label table',
        argstr='-drop-unused-labels')

class VolumeLabelImportOutputSpec(TraitedSpec):
    out_file = File(desc='output file', exists=True)

class VolumeLabelImport(WBCommand):
    _cmd = 'wb_command -volume-label-import'
    input_spec = VolumeLabelImportInputSpec
    output_spec = VolumeLabelImportOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs

class MetricToVolumeMappingRCInputSpec(CommandLineInputSpec):
    in_file = File(
        desc='the input metric file',
        argstr='%s',
        position=0,
        mandatory=True,
        exists=True,
        copyfile=False)    
    surface = File(
        desc='the surface to use coordinates from',
        argstr='%s',
        position=1,
        mandatory=True,
        exists=True,
        copyfile=False)    
    volume_space = File(
        desc='a volume file in the desired output volume space',
        argstr='%s',
        position=2,
        mandatory=True,
        exists=True,
        copyfile=False)    
    out_file = File(
        desc='the output volume file',
        argstr='%s',
        position=3,
        mandatory=True)
    inner_surf = File(
        desc='inner surface of the ribbon',
        argstr='-ribbon-constrained %s',
        position=4,
        mandatory=True,
        exists=True,
        copyfile=False)
    outer_surf = File(
        desc='outer surface of the ribbon',
        argstr='%s',
        position=5,
        mandatory=True,
        exists=True,
        copyfile=False)
    voxel_subdiv = traits.Int(
        desc=('voxel division while estimating voxel weights\n'
              '<subdiv-num> number of subdivisions, default 3'),
        argstr='-voxel-subdiv %d')
    greedy = traits.Bool(
        desc=('instead of antialiasing partial-volumed voxels, put full\n'
              'metric values (legacy behavior)'),
        argstr='-greedy')
    thick_columns = traits.Bool(
        desc = 'use overlapping columns (legacy method)',
        argstr='-thick-columns')

class MetricToVolumeMappingRCOutputSpec(TraitedSpec):
    out_file = File(desc='output file', exists=True)

class MetricToVolumeMappingRC(CommandLine):
    _cmd = 'wb_command -metric-to-volume-mapping'
    input_spec = MetricToVolumeMappingRCInputSpec
    output_spec = MetricToVolumeMappingRCOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs

class CiftiConvertToNiftiInputSpec(CommandLineInputSpec):
    in_file = File(
        desc='the input cifti file',
        argstr='%s',
        position=0,
        mandatory=True,
        exists=True,
        copyfile=False)
    out_file = File(
        desc='the output nifit file',
        argstr='%s',
        position=1,
        mandatory=True)
    smaller_file = traits.Bool(
        desc='use bettter-fitting dimension lenghts',
        argstr='-smaller-file')
    smaller_dims=traits.Bool(
        desc='minimize the largest dimension, for tools taht do not like large indices',
        argstr='--smaller_dims')

class CiftiConvertToNiftiOutputSpec(TraitedSpec):
    out_file = File(desc='output file', exists=True)

class CiftiConvertToNifti(WBCommand):
    _cmd = 'wb_command -cifti-convert -to-nifti'
    input_spec = CiftiConvertToNiftiInputSpec
    output_spec = CiftiConvertToNiftiOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs

class CiftiConvertFromNiftiInputSpec(CommandLineInputSpec):
    in_file = File(
        desc='the input nifti file',
        argstr='%s',
        position=0,
        mandatory=True,
        exists=True,
        copyfile=False)
    cifti_template = File(    
        desc='a cifti file with th dimension(s) and mapping(s) that should be used',
        argstr='%s',
        position=1,
        mandatory=True,
        exists=True,
        copyfile=False)
    out_file = File(    
        desc='the output cifti file',
        argstr='%s',
        position=2,
        mandatory=True)
    reset_timepoints = traits.Tuple(
        traits.Float, 
        traits.Float,
        desc=('reset the mapping along rows to timepoints\n'
              'taking length from the nifti file\n'
              '<timestep> - the desired time between frames\n'
              '<timestart> - the desired tiem offset of the intitial frame'),
        argstr='-reset-timepoints %f %f')
    unit = traits.Str(
        desc='use a unit other than time (default SECOND)',
        requires=['reset_timepoints'],
        argstr='-unit %s')
    reset_scalars = traits.Bool(
        desc='reset mapping along rows to scalars, taking lenght from the nifti file',
        argstr='-reset-scalars')

class CiftiConvertFromNiftiOutputSpec(TraitedSpec):
    out_file = File(desc='output file', exists=True)

class CiftiConvertFromNifti(WBCommand):
    _cmd = 'wb_command -cifti-convert -from-nifti'
    input_spec = CiftiConvertFromNiftiInputSpec
    output_spec = CiftiConvertFromNiftiOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs

class CiftiParcellateInputSpec(CommandLineInputSpec):
    in_file = File(
        desc='the cifti file to parcellate',
        argstr='%s',
        position=0,
        mandatory=True,
        exists=True,
        copyfile=False)
    cifti_label = File(
        desc='a cifti label file to use for the parcellation',
        argstr='%s',
        position=1,
        mandatory=True,
        exists=True,
        copyfile=False)
    direction = traits.Str(
        desc='which mapping to parcellate (integer, ROW, or COLUMN)',
        argstr='%s',
        position=2,
        mandatory=True)
    out_file = File(
        desc='output cifti file',
        argstr='%s',
        position=3,
        mandatory=True)

class CiftiParcellateOutputSpec(TraitedSpec):
    out_file = File(desc='output cifti file', exists=True)

class CiftiParcellate(WBCommand):
    _cmd = 'wb_command -cifti-parcellate'
    input_spec = CiftiParcellateInputSpec
    output_spec = CiftiParcellateOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs

class CiftiCorrelationInputSpec(CommandLineInputSpec):
    in_file = File(
        desc='input cifti file',
        argstr='%s',
        position=0,
        mandatory=True,
        exists=True,
        copyfile=False)
    out_file = File(
        desc='output cifti file',
        argstr='%s',
        position=1,
        mandatory=True)
    fisher_z = traits.Bool(
        desc='apply fisher small z transform (ie, artanh) to correlation',
        argstr='-fisher-z',
        xor=['covariance'])
    covariance = traits.Bool(
        desc='compute covariance instead of correlation',
        argstr='-covariance',
        xor=['fisher_z'])

class CiftiCorrelationOutputSpec(TraitedSpec):
    out_file = File(desc='output cifti file', exist=True)

class CiftiCorrelation(WBCommand):
    _cmd = 'wb_command -cifti-correlation'
    input_spec = CiftiCorrelationInputSpec
    output_spec = CiftiCorrelationOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs

class CiftiCovarianceInputSpec(CommandLineInputSpec):
    in_file = File(
        desc='input cifti file',
        argstr='%s',
        position=0,
        mandatory=True,
        exists=True,
        copyfile=False)
    out_file = File(
        desc='output cifti file',
        argstr='%s',
        position=1,
        mandatory=True)
    

class CiftiMergeInputSpec(CommandLineInputSpec):
    in_files = InputMultiPath(
        File(exist=True),
        desc='input cifti files',
        argstr='-cifti %s',
        copyfile=False,
        sep=' -cifti ',
        position=1,
        mandatory=True)
    out_file = File(
        desc='output cifti file',
        argstr='%s',
        position=0,
        mandatory=True)

class CiftiMergeOutputSpec(TraitedSpec):
    out_file = File(desc='output cifti file', exist=True)

class CiftiMerge(WBCommand):
    _cmd = 'wb_command -cifti-merge'
    input_spec = CiftiMergeInputSpec
    output_spec = CiftiMergeOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs

def check_wb():
    ver = Info.version()
    if ver:
        return 0
    else:
        return 1

def no_wb():
    """Checks if wb_command is available"""
    if Info.version() is None:
        return True
    return False

