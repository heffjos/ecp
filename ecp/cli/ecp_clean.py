import os
import sys
import json

import numpy as np
import pandas as pd
import nibabel as nib

from pathlib import Path

from ecp import utils
from argparse import ArgumentParser

from ecp.interfaces.paths import (
    PostFreeSurferFiles, HcpTaskCiftiFiles, CleanPrepFiles
)

from os.path import join as opj
from ecp.workflows.base import init_clean_wf

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

import ecp.data as data

# skipped volumes are put in the regressors tsv and they are zero based
#   replace these with zeros before putting into 3dTproject

clean_regressors = [
    'csf',
    'csf_derivative1',
    'csf_power2',
    'csf_derivative1_power2',
    'white_matter',
    'white_matter_derivative1',
    'white_matter_derivative1_power2',
    'white_matter_power2',
    'global_signal',
    'global_signal_derivative1',
    'global_signal_power2',
    'global_signal_derivative1_power2',
    'std_dvars',
    'dvars',
    'framewise_displacement',]

# there are additionally: 
#   cosine\d\d
#   motion_outlier\d\d

alias_regressors = {
    'motion': ['trans_x', 
               'trans_y', 
               'trans_z', 
               'rot_x', 
               'rot_y', 
               'rot_z',],
    'motion_derivatives': ['trans_x_derivative1', 
                           'trans_y_derivative1', 
                           'trans_z_derivative1',
                           'rot_x_derivative1',
                           'rot_y_derivative1',
                           'rot_z_derivative1',],
    'motion_squared': ['trans_x_power2',
                       'trans_y_power2',
                       'trans_z_power2',
                       'rot_x_power2',
                       'rot_y_power2',
                       'rot_z_power2',],
    'motion_derivatives_squared': ['trans_x_derivative1_power2',
                                   'trans_y_derivative1_power2',
                                   'trans_z_derivative1_power2',
                                   'rot_x_derivative1_power2',
                                   'rot_y_derivative1_power2',
                                   'rot_z_derivative1_power2',],
}

parcellations = ['glasser']

def get_parser():
    """Define parse object"""

    parser = ArgumentParser(description='cleans prepped ECP resting state data')
    parser.add_argument('bidsify_dir', action='store', help='the bidsify directory')
    parser.add_argument('multitimeseries_dir', action='store', help='multitimeseries directory')
    parser.add_argument('work_dir', action='store', help='the working directory')
    parser.add_argument('out_dir', action='store', help='the output directory')
    parser.add_argument('participant', action='store', help='participant to be cleaned')
    parser.add_argument('clean_name', action='store', help='name for cleaning process')
    parser.add_argument('clean_desc', action='store', help='description entity value')

    parser.add_argument('--hcp_funcs', action='store', help='hcp func names', nargs='+',
                        required=True)

    parser.add_argument('--n-procs', action='store', type=int, default=8,
                        help='number of processors to use')

    # cleaning arguments
    parser.add_argument('--regressors', action='store', nargs='+',
                        help='regressors used for cleaning', 
                        choices=clean_regressors)
    parser.add_argument('--grouped-regressors', action='store', nargs='+',
                        help='adds multilple regressors; names are what you expect',
                        choices=list(alias_regressors.keys()))
    parser.add_argument('--fd-censor', action='store', type=float, 
                        help='censor volumes with frame displacement > FD_CENSOR in 3dTproject'
                             'sensible FD_CENSOR = 0.9, 0.5, or 0.2')
    parser.add_argument('--dvars-censor', action='store', type=float,
                        help='censor volumes with DVARS > DVARS_CENSOR in 3dTproject')
    parser.add_argument('--polort', type=int, default=2,
                        help='Remove polynomials up to and including degree POLORT')

    parser.add_argument('--passband', action='store', type=float, nargs=2, 
                        metavar=('FBOT', 'FTOP'),
                        help='Remove all frequencies EXCEPT those in the range FBOT..FTOP')
    parser.add_argument('--stopband', action='store', type=float, nargs=2,
                        metavar=('SBOT', 'STOP'),
                        help='Remove all frequencies in teh range SBOT..STOP')

    parser.add_argument('--parcellations', action='store', nargs='+',
                        help='cifti parcellation templates',
                        choices=parcellations)
    parser.add_argument('--remove-non-steady-state', action='store_true',
                        help='remove marked non steady state volumes before anything is done')

    acompcor = parser.add_mutually_exclusive_group()
    acompcor.add_argument('--n-acompcor-separate', action='store', type=int,
                          help='number of anatomical principle components (csf and wm separated)')
    acompcor.add_argument('--n-acompcor-combined', action='store', type=int,
                          help='number of anatomical principle components (csf and wm combined)')
    acompcor.add_argument('--pvar-acompcor-separate', action='store', type=float,
                          help='percent variance of anatomical principle components\n'
                               'the maximum value is 50 (csf and wm separated)')
    acompcor.add_argument('--pvar-acompcor-combined', action='store', type=float,
                          help='percent variance of anatomical principle components\n'
                               'the maximum value is 50 (csf and wm combined)')

    tcompcor = parser.add_mutually_exclusive_group()
    tcompcor.add_argument('--n-tcompcor', action='store', type=int,
                          help='number of tSTD principle components')
    tcompcor.add_argument('--pvar-tcompcor', action='store', type=float,
                          help='percent variance of tSTD principle components\n'
                               'the maximum value is 50')

    # save options
    parser.add_argument('--save-clean-dtseries', action='store_true',
                        help='save the clean dtseries for each cifti')
    parser.add_argument('--save-clean-ptseries', action='store_true',
                        help='save the clean ptseries for each cifti')
    parser.add_argument('--save-clean-pconn', action='store_true',
                        help='save the clean pconn for each in cifti')
    parser.add_argument('--save-clean-covariance', action='store_true',
                        help='save the clean covariance for each in cifti')

    parser.add_argument('--testing', action='store_true',
                        help='do everything excep run the workflow')
                        

    return parser

def setup_clean(args):
    """
    Setups a cleaning dictionary from the input arguments.

    Parameters
    ----------

    args: argparse args
        the argparse arguments from get_parser

    Returns
    -------

    setups: dict
        
        cifti         - cifti file
        source_file   - template file name used for naming purposes
        fd_censor     - censor frames with frame displacement greater than
                        this value
        dvars_censor  - censor frames with DVARS greater than this value
        confounds_tsv - confound file for functional output by cleanprep
        confounds     - confound columns names to create ort file
        dt            - functional repetition time
        skip_vols     - remove these volumes from beginning before anything
                        is done
         

    """
        
    bidsify_dir = Path(args.bidsify_dir).resolve()
    mts_dir = Path(args.multitimeseries_dir).resolve()
    work_dir = Path(args.work_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    participant = args.participant
    hcp_funcs = args.hcp_funcs
    clean_name = args.clean_name

    participant_bidsify_dir = bidsify_dir / f'sub-{participant}' / 'func'
    participant_mts_dir = mts_dir / f'sub-{participant}' / 'func'

    # check if participant was run previously
    if not participant_mts_dir.is_dir():
        raise Exception(f'Subject directory does not exist: {participant_mts_dir}')

    if args.pvar_acompcor_combined and args.pvar_acompcor_combined > 50:
        raise Exception('Maximum percent variance allowed for acompcor is 50'
                        'Your input: {}'.format(args.pvar_acompcor_combined))

    if args.pvar_acompcor_separate and args.pvar_acompcor_separate > 50:
        raise Exception('Maximum percent varaince allowed for acompcor is 50'
                        'Your input: {}'.format(args.pvar_acompcor_separate))

    if args.pvar_tcompcor and arg.pvar_tcompcor > 50:
        raise Exception('Maximum percent variace allowed for tcompcor is 50'
                        'Your input: {}'.format(args.pvar_tcompcor))

    common_regressors = args.regressors if args.regressors else []

    if args.grouped_regressors:
        for grouped_regressors in args.grouped_regressors:
            common_regressors.extend(alias_regressors[grouped_regressors])

    setups = {}
    for func in hcp_funcs:
        info = {}

        bold = utils.hcp_to_bids(func, participant)
        entities = utils.get_entities(bold)

        def _generate_bids_name(subject, suffix, ext, desc, space=None):
            return utils.generate_bold_name(
                subject, entities['task'], suffix, ext, dir=entities['dir'], 
                run=entities['run'], desc=desc, space=space)

        # make paths to relevant files
        confounds_tsv = _generate_bids_name(participant, 'timeseries', '.tsv', 'confounds')
        confounds_tsv = participant_mts_dir / confounds_tsv

        confounds_json = _generate_bids_name(participant, 'timeseries', '.json', 'confounds')
        confounds_json = participant_mts_dir / confounds_json

        cifti = _generate_bids_name(participant, 'bold', '.dtseries.nii', None, 'fsLR32k')
        cifti = participant_bidsify_dir / cifti
        
        confounds = pd.read_csv(confounds_tsv, sep='\t')
        with open(confounds_json) as in_json:
            confounds_json = json.load(in_json)

        # handle acompcor
        comp_num = lambda x: int(x[0].rsplit('_')[-1])

        combined = sorted([(x, y['CumulativeVarianceExplained'])
                           for x,y in confounds_json.items()
                           if y['Method'] == 'aCompCor'
                           and y['Mask'] == 'combined'
                           and y['Retained']], key=comp_num)

        csf = sorted([(x, y['CumulativeVarianceExplained'])
                      for x,y in confounds_json.items()
                      if y['Method'] == 'aCompCor'
                      and y['Mask'] == 'CSF'
                      and y['Retained']], key=comp_num)

        wm = sorted([(x, y['CumulativeVarianceExplained']) 
                    for x,y in confounds_json.items()
                    if y['Method'] == 'aCompCor'
                    and y['Mask'] == 'WM'
                    and y['Retained']], key=comp_num)
        
        if args.n_acompcor_combined:
            max_n = len(combined) 
        
            if args.n_acompcor_combined > max_n:
                raise Exception(
                    f'Input n-acompcor-combined ({args.n_acompcor_combined}) exceeds available components ({max_n})'
                )

            acompcor_regressors = [x for x,y in combined[:args.n_acompcor_combined]]

        elif args.n_acompcor_separate:
            max_n = min([len(csf), len(wm)])

            if args.n_acompcor_separate > max_n:
                raise Exception(
                    f'Input n-acompcor-separate ({args.n_acompcor_separate}) exceeds available components ({max_n})'
                )

            acompcor_regressors = ([x for x,y in csf[:args.n_acompcor_separate]] 
                                   + [x for x,y in wm[:args.n_acompcor_separate]])
        
        elif args.pvar_acompcor_combined:
            pvar_acompcor = args.pvar_acompcor_combined / 100
        
            acomp_cor_regressors = [x for x,y in combined if y < pvar_acompcor]
        elif args.pvar_acompcor_separate:
            pvar_acompcor = args.pvar_acompcor_separate / 100

            acompcor_regressors = ([x for x,y in wm if y < pvar_acompcor]
                                   + [x for x,y in csf if y < pvar_acompcor])
        else:
            acompcor_regressors = []
       
        # handle tcompcor
        temporal = sorted([(x, y['CumulativeVarianceExplained'])
                          for x,y in confounds_json.items()
                          if y['Method'] == 'tCompCor'
                          and y['Retained']], key=comp_num)

        if args.n_tcompcor:
            max_n = len(temporal) 
        
            if args.n_tcompcor > max_n:
                raise Exception(f'Input n-tcompcor ({args.n_tcompcor}) exceeds available components ({max_n})')
        
            tcompcor_regressors = [x for x,y in temporal[:args.n_tcompcor]]

        elif args.pvar_tcompcor:
            pvar_tcompcor = pvar_tcompcor / 100
        
            tcompcor_regressors = [x for x,y in temporal if y < pvar_tcompcor]
        else:
            tcompcor_regressors = []

        # handle source file - make it very basic
        source_file = _generate_bids_name(participant, 'bold', '.dtseries.nii', None)

        # get dt
        dt = nib.load(cifti).header.matrix[0].series_step

        # get skip skip_vols
        with open(confounds_tsv) as f:
            lines = f.readlines()
        header = lines[0].split()

        if args.remove_non_steady_state:
            skip_vols = len([x for x in header
                             if x.startswith('non_steady_state_outlier')])
        else:
            skip_vols = 0

        # assign to dictionary
        info['cifti'] = cifti
        info['source_file'] = source_file
        info['fd_censor'] = args.fd_censor
        info['dvars_censor'] = args.dvars_censor
        info['confounds_tsv'] = confounds_tsv
        info['confounds'] = (common_regressors 
                             + acompcor_regressors 
                             + tcompcor_regressors)
        info['dt'] = dt
        info['skip_vols'] = skip_vols

        setups[func] = info

    return setups

def create_confound_files(setup, args, out_dir):
    """
    Creates confound files to disk.

    Parameters
    ----------

    setup: dict
        image setup dictionary
    args: parse_args
        arguments from parse_args
    out_dir: Path
        output directory

    Returns
    -------

    setup: dict
        image setup dictionary
            ort_tsv confound file name in tsv format
            ort_1D confound file name in 1D format
            censor censor file name in 1D format            
    """

    regressors = pd.read_csv(setup['confounds_tsv'], sep='\t')
    regressors = regressors.iloc[setup['skip_vols']:]
    confounds = regressors[setup['confounds']].fillna(0)

    confounds_fname = os.path.basename(setup['confounds_tsv'])
    confounds_entities = utils.get_entities(confounds_fname)

    out_tsv = opj(out_dir, utils.generate_bold_name(
        suffix='timeseries', ext='.tsv', **confounds_entities))
    out_1D = opj(out_dir, utils.generate_bold_name(
        suffix='timeseries', ext='.1D', **confounds_entities))

    confounds.to_csv(out_tsv, sep='\t', na_rep='n/a', index=False)
    confounds.to_csv(out_1D, sep=' ', na_rep='n/a', header=False, index=False)

    censor = np.ones((regressors.shape[0], 1))
    if setup['fd_censor']:
        fd = regressors['framewise_displacement']
        censor[fd > setup['fd_censor']] = 0
    if setup['dvars_censor']:
        dvars = regressors['dvars']
        censor[dvars > setup['dvars_censor']] = 0

    out_censor = opj(out_dir, utils.generate_bold_name(
            suffix='censor', ext='.1D', **confounds_entities))
    np.savetxt(out_censor, censor, fmt='%d')

    setup['ort_tsv'] = out_tsv
    setup['ort_1D'] = out_1D
    setup['censor'] = out_censor

def run_clean_wf(args):

    setups = setup_clean(args)

    participant = args.participant
    work_dir = Path(args.work_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    clean_name = args.clean_name
    clean_desc = args.clean_desc
    participant = args.participant
    n_procs = args.n_procs
    testing = args.testing

    clean_out_dir = out_dir.joinpath(clean_name, f'sub-{participant}')
    clean_out_dir.mkdir(parents=True, exist_ok=True)

    for setup in setups.values():
        create_confound_files(setup, args, clean_out_dir)

    in_files = {
        'cifti': [],
        'source_files': [],
        'censor': [],
        'ort': [],
        'dt': [],
        'trim': []
    }
    for run, setup in setups.items():
        in_files['cifti'].append(setup['cifti'])
        in_files['source_files'].append(setup['source_file'])
        in_files['censor'].append(setup['censor'])
        in_files['ort'].append(setup['ort_1D'])
        in_files['dt'].append(setup['dt'])
        in_files['trim'].append(setup['skip_vols'])

    for parcellation in args.parcellations:

        if parcellation == 'glasser':
            with pkg_resources.path(data, 'glasser_conte.dlabel.nii') as tmp:
                parcellation_file = str(tmp)

        wf_name = f'{clean_name}_{participant}_{parcellation}_{args.clean_desc}_clean'
        clean_wf = init_clean_wf(
            in_files,
            parcellation_file,
            out_dir.as_posix(),
            clean_name,
            'fsLR32k',
            parcellation,
            f'sub-{participant}_task-rest_bold.nii.gz',
            work_dir.as_posix(),
            args.save_clean_dtseries,
            args.save_clean_ptseries,
            args.save_clean_pconn,
            args.save_clean_covariance,
            args.polort,
            args.passband,
            args.stopband,
            args.clean_desc,
            wf_name)

        if not args.testing:
            clean_wf.run(plugin='MultiProc', plugin_args={'n_procs': n_procs})

    return 0

def main():
    parser = get_parser()
    args = parser.parse_args()

    return run_clean_wf(args)        

if __name__ == "__main__":
    main()

# create these files:
#        regressors.tsv
#        regressors.1D (same as tsv above, just without header)
#        censor.1D
# These files will serve as documentation to know what was used during
# cleaning.
#
