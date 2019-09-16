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
    skipped_vols = traits.Int(desc='skipped volumes',
                              usedefautl=True)


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
        param = data[self.inputs.skipped_vols:, 0:6]
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
    gs_deriv_sq = File(desc='mean gs derivative squared', exists=True)

    tcc = File(desc='temporal comp cor', exists=True)
    acc = File(desc='anatomical comp cor', exists=True)
    
    fd = File(desc='frame displacement', exists=True)
    std_dvars = File(desc='std dvars', exists=True)
    dvars = File(desc='dvars', exists=True)

class RegressorsTsvTo1D(SimpleInterface):
    """Converts the regressors.tsv file to 1D files for 3dTproject"""
    input_spec = RegressorsTsvTo1DInputSpec
    output_spec = RegressorsTsvTo1DOutputSpec

    def _run_interface(self, runtime):
        out_dir = runtime.cwd
        data = pd.read_csv(self.inputs.regressors, sep='\t')
        file_table = {
            'csf': 'csf',
            'wm': 'white_matter',
            'gs': 'global_signal',
            'csf_deriv': 'csf_derivative1',
            'wm_deriv': 'white_matter_derivative1',
            'gs_deriv': 'global_signal_derivative1',
            'csf_sq': 'csf_power2',
            'wm_sq': 'white_matter_power2',
            'gs_sq': 'global_singal_power2',
            'csf_deriv_sq': 'csf_derivative1_power2',
            'wm_deriv_sq': 'white_matter_derivative1_power2',
            'gs_deriv_sq': 'global_signal_derivative1_power2',
            'fd': 'framewise_displacement',
            'std_dvars': 'std_dvars',
            'dvars': 'dvars'
        }

        pd.to_csv(
            os.path.join(out_dir, 'movement.1D'),
            data[[x for x in data.columns if re.match('^(rot|trans)_[xyz]$', x)]],
            header=False, 
            index=False,
            sep=' ')
        self._results['movement'] = os.path.join(out_dir, 'movement.1D')

        pd.to_csv(
            os.path.join(out_dir, 'movement_sq.1D'),
            data[[x for x in data.columns if re.match('^(rot|trans)_[xyz]_power2$', x)]],
            header=False, 
            index=False,
            sep=' ')
        self._results['movement_sq'] = os.path.join(out_dir, 'movement_sq.1D')

        pd.to_csv(
            os.path.join(out_dir, 'movement_deriv.1D'),
            data[[x for x in data.columns if re.match('^(rot|trans)_[xyz]_derivative1$', x)]],
            header=False, 
            index=False,
            sep=' ')
        self._results['movement_deriv'] = os.path.join(out_dir, 'movement_deriv.1D')

        pd.to_csv(
            os.path.join(out_dir, 'movement_deriv_sq.1D'),
            data[[x for x in data.columns if re.match('^(rot|trans)_[xyz]_derivative1_power2$', x)]],
            header=False, 
            index=False,
            sep=' ')
        self._results['movement_deriv_sq'] = os.path.join(out_dir, 'movement_deriv_sq.1D')

        for file_name, column in file_table.items():
            out_file = os.path.join(out_dir, file_name + '.1D')
            pd.to_csv(out_file, data[column], header=False, index=False, sep=' ')
            self._results[file_name] = out_file

        pd.to_csv(
            os.path.join(out_dir, 'tcc.1D'),
            data[[x for x in data.columns if re.match('^t_comp_cor_\d+$', x)]],
            header=False,
            index=False,
            sep=' ')
        self._results['tcc'] = os.path.join(out_dir, 'tcc.1D')

        pd.to_csv(
            os.path.join(out_dir, 'acc.1D'),
            data[[x for x in data.columns if re.match('^a_comp_cor_\d+$', x)]],
            header=False,
            index=False,
            sep=' ')
        self._results['acc'] = os.path.join(out_dir, 'acc.1D')

        return runtime
