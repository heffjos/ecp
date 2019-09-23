import os
import sys
import json
import pathlib

import pandas as pd

from ecp import utils
from argparse import ArgumentParser
from nipype.interfaces import afni

from ecp.interfaces.paths import (
    PostFreeSurferFiles, HcpTaskCiftiFiles, CleanPrepFiles
)

from os.path import join as opj
from ecp.interfaces import workbench as wb

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
    parser.add_argument('data_dir', action='store', help='the cleanprep directory')
    parser.add_argument('cleanprep_dir', action='store', help='the cleanprep directory')
    parser.add_argument('work_dir', action='store', help='the working directory')
    parser.add_argument('out_dir', action='store', help='the output directory')
    parser.add_argument('clean_name', action='store', help='name for cleaning process')
    parser.add_argument('spec_file', action='store', help='csv spec file')
    parser.add_argument('--participants', action='store', nargs='+',
                        help='participants to be prepped')
    parser.add_argument('--n-procs', action='store', type=int, default=None,
                        help='number of processors to use')

    parser.add_argument('--regressors', action='store', nargs='+',
                        help='regressors used for cleaning', 
                        choices=clean_regressors)
    parser.add_argument('--grouped-regressors', action='store', nargs='+',
                        help='adds multilple regressors; names are what you expect',
                        choices=list(alias_regressors.keys()))

    acompcor = parser.add_mutually_exclusive_group()
    acompcor.add_argument('--n-acompcor', action='store', type=int,
                          help='number of anatomical principle components')
    acompcor.add_argument('--pvar-acompcor', action='store', type=float,
                          help='percent variance of anatomical principle components\n'
                               'the maximum value is 50')

    tcompcor = parser.add_mutually_exclusive_group()
    tcompcor.add_argument('--n-tcompcor', action='store', type=int,
                          help='number of tSTD principle components')
    tcompcor.add_argument('--pvar-tcompcor', action='store', type=float,
                          help='percent variance of tSTD principle components\n'
                               'the maximum value is 50')

    censor = parser.add_mutually_exclusive_group()
    censor.add_argument('--fd-censor', action='store', type=float, 
                        help='censor volumes with frame displacement > FD_CENSOR in 3dTproject'
                             'sensible FD_CENSOR = 0.9, 0.5, or 0.2')
    censor.add_argument('--fmriprep-censor', action='store_true',
                        help='censor volumes in 3dTproject with the motion_outliers\n'
                             'generated by fmriprep')

    parser.add_argument('--fd-censor-method', action='store', default='KILL',
                        choices=['ZERO', 'KILL', 'NTRP'],
                        help='specifiy how censored time points are treated\n'
                             'only used if --fd-censor is specified')

    parser.add_argument('--polort', type=int, default=2,
                        help='Remove polynomials up to and including degree POLORT')

    parser.add_argument('--passband', action='store', type=float, nargs=2, 
                        metavar=('FBOT', 'FTOP'),
                        help='Remove all frequencies EXCEPT those in the range FBOT..FTOP')
    parser.add_argument('--stopband', action='store', type=float, nargs=2,
                        metavar=('SBOT', 'STOP'),
                        help='Remove all frequencies in teh range SBOT..STOP')

    parser.add_argument('--dt', action='store', type=float,
                        help='Use time step dt for the frequncey calculations')

    parser.add_argument('--parcellations', action='store', nargs='+',
                        help='cifti parcellation templates',
                        choices=parcellations)
    
    return parser

def run_clean_setup(args):

    data_dir = args.data_dir
    cleanprep_dir = args.cleanprep_dir
    work_dir = args.work_dir
    out_dir = args.out_dir
    clean_name = args.clean_name
    spec_file = args.spec_file

    participants = args.participants

    if not os.path.isdir(out_dir):
        raise Exception(f'Out directory does not exist: {out_dir}')

    spec = pd.read_csv(spec_file)
    spec = spec[spec['has_new_nifti'] & spec['has_cifti']]

    spec_participants = set(spec['subject'])
    extra_participants = set(participants).difference(spec_participants)
    if extra_participants:
        raise Exception('Input participants are not in spec_file: '
                        '{}'.format('\n'.join(extra_participants)))

    if args.pvar_acompcor and args.pvar_acompcor > 50:
        raise Exception('Maximum percent variance allowed for acompcor is 50'
                        'Your input: {}'.format(args.pvar_acompcor))

    if args.pvar_tcompcor and arg.pvar_tcompcor > 50:
        raise Exception('Maximum percent variace allowed for tcompcor is 50'
                        'Your input: {}'.format(args.pvar_tcompcor))

    common_regressors = args.regressors if args.regressors else []

    if args.grouped_regressors:
        for grouped_regressors in args.grouped_regressors:
            common_regressors.extend(alias_regressors[grouped_regressors])

    setups = {}
    for participant in participants:

        setups[participant] = {}
        print(f'Setting up participant: {participant}')

        participant_info = spec[spec['subject'] == participant]
        tasks = list(participant_info['task'])

        anat_files = PostFreeSurferFiles(base_dir=data_dir,
                                         subject=participant).run().outputs 
        
        for task in tasks:
            info = {}

            task_cifti_files = HcpTaskCiftiFiles(
                mninonlinear=anat_files.mninonlinear,
                subject=participant,
                task=task).run().outputs

            clean_prep_files = CleanPrepFiles(
                cleanprep_dir=cleanprep_dir,
                subject=participant,
                hcp_task=task).run().outputs

            confounds = pd.read_csv(clean_prep_files.confounds_tsv, sep='\t')
            with open(clean_prep_files.confounds_json) as in_json:
                confounds_json = json.load(in_json)

            if args.n_acompcor:
                max_n = max([int(x.rsplit('_', 1)[1]) 
                             for x in confounds.columns
                             if 'a_comp_cor' in x])

                if args.n_acompcor > max_n:
                    raise Exception(f'{participant} {task}: Input n-acompcor ({args.n_acompcor}) exceeds available components ({max_n})')

                acomp_cor_regressors = [f'a_comp_cor_{x:02d}' 
                                        for x in range(args.n_acompcor)]
            elif args.pvar_acompcor:
                pvar_acompcor = pvar_acompcor / 100

                acomp_cor_regressors = [x for x in confounds_json
                    if 'a_comp_cor' in x 
                    and confounds_json[x]['CumlativeVarianceExplained'] < pvar_acompcor]
            else:
                acomp_cor_regressors = []

            if args.n_tcompcor:
                max_n = max([int(x.rsplit('_', 1)[1]) 
                             for x in confounds.columns
                             if 't_comp_cor' in x])

                if args.n_tcompcor > max_n:
                    raise Exception(f'{participant} {task}: Input n-tcompcor ({args.n_tcompcor}) exceeds available components ({max_n})')

                tcomp_cor_regressors = [f't_comp_cor_{x:02d}' 
                                        for x in range(args.n_tcompcor)]
            elif args.pvar_tcompcor:
                pvar_tcompcor = pvar_tcompcor / 100

                tcomp_cor_regressors = [x for x in confounds_json
                    if 't_comp_cor' in x 
                    and confounds_json[x]['CumlativeVarianceExplained'] < pvar_tcompcor]
            else:
                tcomp_cor_regressors = []

            if args.fmriprep_censor:
                fmriprep_censor = [x for x in confounds.columns if 'motion_outlier' in x]
            else:
                fmriprep_censor = []
                                        
            info['cifti_ref'] = task_cifti_files.preproc
            info['confounds_tsv'] = clean_prep_files.confounds_tsv
            info['fakenifti'] = clean_prep_files.fakenifti
            info['to_keep'] = (common_regressors 
                               + acomp_cor_regressors 
                               + tcomp_cor_regressors
                               + fmriprep_censor)

            setups[participant][task] = info

    return setups

def run_clean(setups, args):

    data_dir = args.data_dir
    cleanprep_dir = args.cleanprep_dir
    work_dir = args.work_dir
    out_dir = args.out_dir
    clean_name = args.clean_name
    participants = args.participants

    clean_out_dir = opj(out_dir, 'clean', clean_name)
    pathlib.Path(clean_out_dir).mkdir(parents=True, exist_ok=True)

    for participant, tasks_info in setups.items():
        participant_dir = os.path.join(clean_out_dir, 'sub-' + participant)
        pathlib.Path(participant_dir).mkdir(exist_ok=True)

        for task, info in tasks_info.items():
            clean_task(participant, task, info, args, participant_dir)

def clean_task(participant, task, info, args, out_dir):
    """cleans the task using the info and args"""

    log_dir = os.path.join(out_dir, 'logs')
    pathlib.Path(log_dir).mkdir(exist_ok=True)

    regressors = pd.read_csv(info['confounds_tsv'], sep='\t')
    regressors = regressors[info['to_keep']].fillna(0)
    regressors_fname = os.path.basename(info['confounds_tsv'])
    regressors_entities = utils.get_entities(regressors_fname)

    out_tsv = opj(out_dir, utils.generate_bold_name(
        suffix='cleanregressors', ext='.tsv', **regressors_entities))
    out_1D = opj(out_dir, utils.generate_bold_name(
        suffix='cleanregressors', ext='.1D', **regressors_entities))

    regressors.to_csv(out_tsv, sep='\t', na_rep='n/a', index=False)
    regressors.to_csv(out_1D, sep=' ', na_rep='n/a', header=False, index=False)

    if args.fd_censor:
        regressors['censor'] = (regressors['framewise_displacement'] < args.fd_censor).astype(int)
        out_censor = opj(out_dir, utils.generate_bold_name(
            suffix='censor', ext='.1D', **regressors_entities))

        regressors[['censor']].to_csv(out_censor, sep=' ', na_rep='n/a', header=False, index=False)

    fakenifti_entities = utils.get_entities(info['fakenifti'])
    clean_nifti_file = opj(out_dir, utils.generate_bold_name(suffix='cleanfakenifti',
                                                             ext='.nii.gz',
                                                             **fakenifti_entities))

    # clean
    tproject = afni.TProject()
    tproject.inputs.in_file = info['fakenifti']
    tproject.inputs.out_file = clean_nifti_file
    
    if args.fd_censor:
        tproject.inputs.censortr = out_censor
        tproject.inputs.cenmode = args.fd_censor_method

    tproject.inputs.polort = args.polort

    if len(info['to_keep']) != 0:
        tproject.inputs.org = out_1D
    
    if args.passband:
        tproject.inputs.bandpass = tuple(args.passband)

    if args.stopband:
        tproject.inputs.stopband = tuple(args.stopband)

    if args.dt:
        tproject.inputs.TR = args.dt

    if args.n_procs:
        tproject.inputs.num_threads = args.n_procs

    tproject_run = tproject.run()
    log = opj(log_dir, utils.generate_bold_name(suffix='cleanfakenift',
                                                ext='.log',
                                                **fakenifti_entities))

    with open(log, 'w') as out_log:
        out_log.writelines(tproject_run.runtime.stderr)

    # convert back to cifti
    out_cifti = opj(out_dir, utils.generate_bold_name(suffix='cleanbold',
                                                      ext='.dtseries.nii',
                                                      **fakenifti_entities))

    cifticonvertfromnifti = wb.CiftiConvertFromNifti(
        in_file=tproject_run.outputs.out_file,
        cifti_template=info['cifti_ref'],
        out_file=out_cifti)

    cifticonvertfromnifti_results = cifticonvertfromnifti.run()

    # parcellate and correlate
    for parcellation in args.parcellations:

        if parcellation == 'glasser':
            with pkg_resources.path(data, 'glasser_conte.dlabel.nii') as tmp:
                parcellation_file = str(tmp)

        parcel_entities = fakenifti_entities.copy()
        parcel_entities['space'] = parcellation
        out_cifti = opj(out_dir, utils.generate_bold_name(suffix='cleanbold',
                                                          ext='.ptseries.nii',
                                                          **parcel_entities))
                
        cifti_parcellate = wb.CiftiParcellate(
            in_file=cifticonvertfromnifti_results.outputs.out_file,
            cifti_label=parcellation_file,
            direction='COLUMN',
            out_file=out_cifti)

        cifti_parcellate_results = cifti_parcellate.run()

        out_cifti = opj(out_dir, utils.generate_bold_name(suffix='cleanbold',
                                                          ext='.pconn.nii',
                                                          **parcel_entities))

        cifti_correlation = wb.CiftiCorrelation(
            in_file=cifti_parcellate_results.outputs.out_file,
            fisher_z=True,
            out_file=out_cifti)

        cifti_correlation_results = cifti_correlation.run()

def main(*args):
    parser = get_parser()
    args = parser.parse_args(*args)
    setups = run_clean_setup(args)
    run_clean(setups, args)

    return 0

if __name__ == "__main__":
    main()
            
    
    
    

    
    
    
        

    

    
        


    
    
            
            


        
    
    

# create these files:
#        regressors.tsv
#        regressors.1D (same as tsv above, just without header)
# These files will serve as documentation to know what was used during
# cleaning.
#
# run 3dTproject with inputs
# convert cleaned fake niftis to cifti
# parcellate the ciftis into the glasser parcellation

# Here is what the workflow would look like:
#   cleanprepfiles
#   get_regressors
#       ort
#       afni censor vector if specified
#   run 3dTproject
#   convert cleaned fake niftis to cifti
#   parcellate the ciftis into the glasser parcellation
#   do correlation

