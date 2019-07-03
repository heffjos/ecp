from nipype.interfaces.base import traits, File
from nipype.interfaces.fsl.maths import MathsInput, MathsCommand
from nipype.interfaces.fsl.utils import ImageMathsInputSpec, ImageMaths

class Augmented1ImageMathsInputSpec(ImageMathsInputSpec):
    in_file3 = File(
        exists=True, 
        argstr='%s', 
        position=5)
    op_string2 = traits.Str(
        argstr='%s',
        position=4,
        desc='string defining second operation')
    op_string_end = traits.Str(
        argstr='%s',
        position=-3,
        desc='end operation')

class Augmented1ImageMaths(ImageMaths):
    input_spec = Augmented1ImageMathsInputSpec

class AddAddDilateInput(MathsInput):
    add_file1 = File(
        position=4,
        argstr='-add %s',
        exists=True,
        mandatory=True,
        desc='first image to add')
    add_file2 = File(
        position=5,
        argstr='-add %s -dilD -dilD',
        exists=True,
        mandatory=True,
        desc='second image to add')

class AddAddDilate(MathsCommand):
    """Use fslmaths to execute '%s -add %s -add %s -dilD -dilD' operation"""
    input_spec = AddAddDilateInput
    
