import os
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

class GetHcpMovementOutputSpec(TraitedSpec):
    movement = File(desc='movement regressors', exists=True)

class GetHcpMovement(SimpleInterface):
    """Grabs the 6 motion parmaters from the HCP motion file"""
    input_spec = GetHcpMovementInputSpec
    output_spec = GetHcpMovementOutputSpec

    def _run_interface(self, runtime):

        data = np.loadtxt(self.inputs.hcp_movement)
        param = data[:, 0:6]
        param[:, 3:] = param[:, 3:] * 2 * np.pi / 360

        out_dir = runtime.cwd
        np.savetxt(os.path.join(out_dir, 'movement.tsv'), param, fmt='%0.6f',
                   delimiter='\t')

        self._results['movement'] = os.path.join(out_dir, 'movement.tsv')

        return runtime

class RegressorsTsvTo1DInputSpec(BaseInterfaceInputSpec):
    regressors = File(desc='regressors.tsv file', 
                      exists=True, 
                      mandatory=True)

class RegressorsTsvTo1DOutputSpec(TraitedSpec):
    movement = File(desc='movement parameters', exists=True)
    movement_sq = File(desc='movement paratmers squared', exists=True)
    movement_deriv = File(desc='movement derivatives', exists=True)
    movement_deriv_sq = File(desc='movement derivatives squared', exists=True)

    csf = File(desc='mean csf time course', exists=True)
    wm = File(desc='mean wm time course', exists=True)
    gs = File(desc='global signal time course', exists=True)

    csf_deriv = File(desc='mean csf derivative', exists=True)
    wm_deriv = File(desc='mean wm derivative', exists=True)
    gs_deriv = File(desc='mean global signal derivative', exists=True)

    csf_sq = File(desc='mean csf squared', exists=True)
    wm_sq = File(desc='mean wm squared', exists=True)
    gs_sq = File(desc='mean gs squared', exists=True)

    csf_deriv_sq = File(desc='mean csf derivative squared', exists=True)
    wm_deriv_sq = File(desc='mean wm derivative squared', exists=True)
    gs_sq = File(desc='mean gs derivative squared', exists=True)

    tcc = File(desc='temporal comp cor', exists=True)
    acc = File(desc='anatomical comp cor', exists=True)
    
    fd = File(desc='frame displacement', exists=True)
    std_dvars = File(desc='std dvars', exists=True)
    dvars = File(desc='dvars', exists=True)

class RegressorsTsvTo1D(SimpleInterface):
    """Converts the regressors.tsv file to 1D files for 3dTproject"""
    input_spec = RegressorsTsvTo1DInputSpec
    output_spec = RegressorsTsvTo1DOuptuSpec

    def _run_interface(self, runtime):
        out_dir = runtime.cwd
        data = pd.read_csv(self.inputs.regressors, sep='\t')

        return runtime
