import os
import re
import numpy as np
import pandas as pd

from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec, File, Directory, isdefined,
    SimpleInterface
)

class GetHcpMovementInputSpec(BaseInterfaceInputSpec):
    hcp_movement = File(desc='HCP movement regressors', 
                        exists=True,
                        mandatory=True)
    skip_begin = traits.Int(desc='skip volumes in the beginning',
                               usedefault=True)
    skip_end = traits.Int(desc='skip volumes at the end',
                             usedefault=True)


class GetHcpMovementOutputSpec(TraitedSpec):
    movement = File(desc='movement regressors', exists=True)

class GetHcpMovement(SimpleInterface):
    """Grabs the 6 motion parmaters from the HCP motion file"""
    input_spec = GetHcpMovementInputSpec
    output_spec = GetHcpMovementOutputSpec

    def _run_interface(self, runtime):

        # for fmriprep
        # rotations are in radians
        # translations are in mm
        # same format as spm

        # for spm
        # rotations are in radians
        # translations are in mm
        # translations are first 3
        # rotations are next 3
        
        # for hcp
        # rotations are in degrees
        # translations are in mm
        # translations are first 3
        # rotations are next 3

        data = np.loadtxt(self.inputs.hcp_movement)

        begin = self.inputs.skip_begin
        end = data.shape[1] - self.inputs.skip_end

        param = data[begin:end, 0:6]
        param[:, 3:] = param[:, 3:] * 2 * np.pi / 360

        out_dir = runtime.cwd
        self._results['movement'] = os.path.join(out_dir, 'movement.tsv')

        np.savetxt(self._results['movement'], param, fmt='%0.6f',
                   delimiter='\t')


        return runtime

class CleaningRegressorsInputSpec(BaseInterfaceInputSpec):
    regressors_tsv = File(desc='fmriprep regressors tsv file', 
                          exists=True, 
                          mandatory=True)
    regressors_to_keep = traits.List(traits.Str(), 
                                     mandatory=True,
                                     desc='subset list of cleaning regressors')
    afni_censor = traits.Float(desc='framewise displacement threshold')

class CleaningRegressorsOutputSpec(TraitedSpec):
    regressors_tsv = File(desc='cleaning regressors', exists=True)
    regressors_1D = File(desc='afni formatted cleanring regressors', exists=True)
    censor_1D = File(desc='afni censor file derived from framewise displacement')

class CleaningRegressors(SimpleInterface):
    """Selects the cleaning regressors"""

    input_spec = CleaningRegressorsInputSpec
    output_spec = CleaningRegressorsInputSpec

    def _run_interface(self, runtime):
        out_dir = runtime.cwd

        data = pd.read_csv(self.inputs.regressors_tsv, sep='\t')
        keep_data = data[self.inputs.regressors_to_keep].fillna(0)

        regressors_fname = os.path.basename(self.inputs.regressors_tsv)
        no_ext_regressors_fname, _ = utils.splitext(regressors_fname)
        no_suffix = no_ext_regressors_fname.rsplit('_')[0]

        out_tsv = opj(out_dir, regressors_fname)
        out_1D = opj(out_dir, no_ext_regressors_fname + '.1D')

        keep_data.to_csv(out_tsv, sep='\t', na='n/a', index=False)
        keep_data.to_csv(out_1D, sep=' ', na='n/a', header=False, index=False) 

        self._results['regressors_tsv'] = out_tsv
        self._results['regressors_1D'] = regressors_1D

        if isdefined(self.inputs.afni_censor):
            data['afni_censor'] = data['framewise_displacement'] < self.inputs.afni_censor
            out_censor = opj(out_dir, no_suffx + '_censor.1D')
            data[['afni_censor']].to_csv(out_censor, sep=' ', na='n/a',
                                         header=False, index=False)
            self._results['censor_1d'] = out_censor
        else:
            self._results['censor_1d'] = ''

        return runtime
