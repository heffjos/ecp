import os
import sys

import pandas as pd

from argparse import ArgumentParser
from multiprocessing import cpu_count

from ecp.workflows.base import init_cleanprep_wf

def get_parser():
    """Define parse object"""
    parser = ArgumentParser(description='preps ECP resting data in HCP format'
                                        'for connectivity analysis')
    parser.add_argument('data_dir', action='store', help='the data directory')
    parser.add_argument('work_dir', action='store', help='the working directory')
    parser.add_argument('out_dir', action='store', help='the output directory')
    parser.add_argument('spec_file', action='store', help='csv spec file')
    parser.add_argument('--participants', action='store', nargs='+',
                        help='participants to be prepped')
    parser.add_argument('--n-procs', action='store', type=int, default=None,
                        help='number of processors to use')

    return parser

def run_cleanprep_wf(args):

    data_dir = args.data_dir
    work_dir = args.work_dir
    out_dir = args.out_dir
    spec_file = args.spec_file

    participants = args.participants
    n_procs = args.n_procs

    spec = pd.read_csv(spec_file)
    spec = spec[spec['has_new_nifti'] & spec['has_cifti']]

    spec_participants = set(spec['subject'])
    extra_participants = set(participants).difference(spec_participants)
    if extra_participants:
        raise Exception('Input participants are not in spec_file: '
                        '{}'.format('\n'.join(extra_participants)))

    n_procs = args.n_procs if args.n_procs else cpu_count()

    wfs = []
    for participant in participants:
        print(f'Building workflow for participant: {participant}')

        participant_info = spec[spec['subject'] == participant]
        tasks = list(participant_info['task'])
        skip_begin = list(participant_info['skip_begin'])
        skip_end = list(participant_info['skip_end'])

        wfs.append(init_cleanprep_wf(
            data_dir, work_dir, out_dir, participant, tasks, skip_begin, skip_end))

    for wf in wfs:
        # workflow.config['execution']['job_finished_timeout'] = 65
        wf.run(plugin='MultiProc', plugin_args={'n_procs' : n_procs,
                                                'maxtasksperchild': 1})

    return 0

def main():
    parser = get_parser()
    args = parser.parse_args()

    return run_cleanprep_wf(args)

if __name__ == '__main__':
    main()

        
