import os
import sys

import pandas as pd

from argparse import ArgumentParser
from multiprocessing import cpu_count

jheffernan = '/rcc/stor1/depts/neurology/users/jheffernan'
sys.path.insert(0, os.path.join(jheffernan, 'repositories', 'ecp'))

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
    parser.add_argument('--skip-vols', action='store', type=int, default=None,
                        help='the number of non-steady state volumes')
    parser.add_argument('--n-procs', action='store', type=int, default=None,
                        help='number of processors to use')

    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    data_dir = args.data_dir
    work_data = args.work_dir
    out_dir = args.out_dir
    spec_file = args.spec_file

    participants = args.participants
    skip_vols = args.skip_vols
    n_procs = args.n_procs

    spec = pd.read_csv(spec_file)
    spec = spec[spec['tsnr'].notnull()]

    spec_participants = set(spec['subject'])
    extra_participants = set(participants).difference(spec_participants)
    if extra_participants:
        raise Exception('Input participants are not in spec_file: '
                        '{}'.format('\n'.join(extra_participants)))

    n_procs = args.n_procs if args.n_procs else cpu_count()

    wfs = []
    for participant in participants:
        print(f'Building workflow for participant: {participant}')

        tasks = list((spec['subject'] == participant)['task'])

        wfs.append(init_cleanprep_wf(
            data_dir, work_dir, out_dir, participant, tasks, skip_vols))

    for wf in wfs:
        wf.run(plugin='MultiProc', plugin_args={'n_procs' : n_procs})

    return 0

if __name__ == '__main__':
    main()

        
