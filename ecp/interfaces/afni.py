import os

from nipype.interfaces.afni.base import AFNICommandOutputSpec
from nipype.interfaces.afni.preprocess import TProjectInputSpec, TProject

from nipype.interfaces.base import (
    CommandLineInputSpec,
    CommandLine,
    TraitedSpec,
    traits,
    isdefined,
    File,
    InputMultiPath,
    Undefined,
    Str,
    InputMultiObject,
)

class TProjectInputSpec(TProjectInputSpec):
    verb = traits.Bool(
        desc='write ort matrix, its singular values and psuedo-inverse',
        argstr='-verb',
        usedefault=True)

class TProjectOutputSpec(AFNICommandOutputSpec):
    matrix = File(desc='fixed ort matrix')
    pseudo_inv = File(desc='matrix pseudo inverse')
    singular_values = File(desc='matrix singular values')

class TProject(TProject):

    input_spec = TProjectInputSpec
    output_spec = TProjectOutputSpec

    def _list_outputs(self):
        outputs = super(TProject, self)._list_outputs()
        if self.inputs.verb:
            out_file = outputs['out_file']
            if out_file.endswith('BRIK'):
                loc = out_file.rfind('+')
                prefix = out_file[:loc]
            else:
                prefix = out_file

            outputs['matrix'] = os.path.abspath(prefix + '.ort.1D')
            outputs['pseudo_inv'] = os.path.abspath(prefix + '.ort_psinv.1D')
            outputs['singular_values'] = os.path.abspath(prefix + '.sval.1D')
        return outputs
            
    

