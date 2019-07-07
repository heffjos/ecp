import os
import numpy as np

from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec, File, Directory, isdefined,
    SimpleInterface
)

class GetHcpMovementInputSpec(BaseInterfaceInputSpec):
    hcp_movement = File(desc='HCP movement regressors', exists=True)

class GetHcpMovementOutputSpec(TraitedSpec):
    movement = File(desc='movement regressors', exists=True)

class GetHcpMovement(SimpleInterface):
    """Converts HCP movement regressor file to ECP regressor files."""
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

