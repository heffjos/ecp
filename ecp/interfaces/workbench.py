from nipype.interfaces.base import (CommandLine, traits, CommandLineInputSpec, isdefined,
                                    File, TraitedSpec, PackageInfo)

# from nipype.interfaces.base import (
#     traits, TraitedSpec, BaseInterfaceInputSpec, File, Directory, isdefined,
#     SimpleInterface
# )

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

class VolumeLabelImport(CommandLine):
    _cmd = 'wb_command -volume-label-import'
    input_spec = VolumeLabelImportInputSpec
    output_spec = VolumeLabelImportOutputSpec

class MetricToVolumeMappingInputSpec(CommandLineInputSpec):
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
    nearest_vertex = traits.Float(
        desc=('use the value from the vertex closest to the voxel\n'
              'distance is how far from the surface to map values to voxels, in mm'),
        argstr='-nearest-vertex %f',
        xor=['ribbon_constrained'])
    ribbon_constrained = traits.Tuple(
        File(exists=True), File(exists=True),
        desc=('use ribbon constrained mapping algorithm\n'
              '<inner-surf> - the inner surface of the ribbon\n'
              '<outer-surf> - the outer surface of the ribbon'),
        argstr='-ribbon-constrained %s %s')
    voxel_subdiv = traits.Int(
        desc=('voxel division while estimating voxel weights\n'
              '<subdiv-num> number of subdivisions, default 3'),
        argstr='-voxel-subdiv %d',
        requires=['ribbon_constrained'])
    greedy = traits.Bool(
        desc=('instead of antialiasing partial-volumed voxels, put full\n'
              'metric values (legacy behavior)'),
        argstr='-greedy',
        requires=['ribbon_constrained'])
    thick_columns = traits.Bool(
        desc = 'use overlapping columns (legacy method)',
        argstr='-thick-columns',
        requires=['ribbon_constrained'])

class MetricToVolumeMappingOutputSpec(TraitedSpec):
    out_file = File(desc='output file', exists=True)

class MetricToVolumeMapping(CommandLine):
    _cmd = 'wb_command -metric-to-volume-mapping'
    input_spec = MetricToVolumeMappingInputSpec
    output_spec = MetricToVolumeMappingOutputSpec
  
