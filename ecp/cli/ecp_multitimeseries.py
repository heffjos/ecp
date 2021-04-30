import pandas as pd
import nibabel as nib

from pathlib import Path
from argparse import ArgumentParser

from ecp.workflows.base import init_multi_timeseries_wf
from ecp.utils import get_entities, hcp_to_bids, generate_bold_name

def get_anat_mask(bidsify_dir, participant, anat):

    mask = f'sub-{participant}_space-mni_desc-{anat}_mask.nii.gz'
    mask = bidsify_dir / f'sub-{participant}' / 'anat' / mask

    return mask

def get_parser():
    """Define parse object"""
    parser = ArgumentParser(description='preps ECP resting data in HCP format'
                                        'for connectivity analysis')
    parser.add_argument('bidsify_dir', action='store', help='the bidsified hcp directory')
    parser.add_argument('work_dir', action='store', help='the working directory')
    parser.add_argument('out_dir', action='store', help='the output directory')
    parser.add_argument('spec_file', action='store', help='csv spec file')
    parser.add_argument('participant', action='store', help='prep participant')
    parser.add_argument('out_path_base', action='store', help='unique name; primarily reflects skip_vols')
    parser.add_argument('--skip_vols', action='store', type=int, 
                        help='remove these non-steady state volumes'
                             'use 0 to remove non'
                             'no input automatically detects how many volumes to remove')
    parser.add_argument('--nthreads', action='store', type=int, default=8,
                        help='number of threads')

    return parser

def run_multi_timeseries_wf(args):

    bidsify_dir = Path(args.bidsify_dir).resolve()
    work_dir = Path(args.work_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    spec_file = Path(args.spec_file).resolve().as_posix()

    participant = args.participant
    out_path_base = args.out_path_base
    skip_vols = args.skip_vols
    nthreads = args.nthreads

    out_dir = out_dir.joinpath('multitimeseries')
    work_dir = work_dir.joinpath(out_path_base)

    spec = pd.read_csv(spec_file, sep='\t')
    spec = spec.loc[spec.usable_data, :]

    if participant not in set(spec.subjects):
        raise Exception(f'Subject not found in spec file: {participant}')

    participant_info = spec[spec.subjects == participant]
    tasks = list(participant_info['func'])

    csf_mask = get_anat_mask(bidsify_dir, participant, 'csf')
    wm_mask = get_anat_mask(bidsify_dir, participant, 'wm')
    cortical_gm_mask = get_anat_mask(bidsify_dir, participant, 'cortgm')

    parameters = []
    func_dir = bidsify_dir / f'sub-{participant}' / 'func'
    for task in tasks:
        task_parameters = {}

        bold = hcp_to_bids(task, participant)
        entities = get_entities(bold)

        def _generate_bids_name(subject, suffix, ext, desc):
            return generate_bold_name(
                subject, entities['task'], suffix, ext, dir=entities['dir'], 
                run=entities['run'], desc=desc)

        task_parameters['bold'] = (func_dir / bold).as_posix()
        task_parameters['source_file'] = f'func/{bold}'
        task_parameters['dt'] = nib.load(task_parameters['bold']).header['pixdim'][4]
        task_parameters['bold_mask'] = (func_dir /
            _generate_bids_name(participant, 'boldmask', '.nii.gz', 'confounds')).as_posix()
        task_parameters['movpar_file'] = (func_dir / 
            _generate_bids_name(participant, 'timeseries', '.tsv', 'movpar')).as_posix()
        task_parameters['skip_vols'] = skip_vols
        task_parameters['wf_name'] = f'{task}_wf'

        parameters.append(task_parameters)

    multi_timeseries_wf = init_multi_timeseries_wf(
        parameters,
        csf_mask,
        wm_mask,
        cortical_gm_mask,
        work_dir,
        out_dir.as_posix(),
        out_path_base,
        name=f'{participant}_multi_timeseries')

    multi_timeseries_wf.run(plugin='MultiProc', plugin_args={'n_procs': nthreads})

    return 0

def main():
    parser = get_parser()
    args = parser.parse_args()

    return run_multi_timeseries_wf(args)

if __name__ == '__main__':
    main()
    
