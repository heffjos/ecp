import os
import sys

from argparse import ArgumentParser
from multiprocessing import cpu_count

jheffernan = '/rcc/stor1/depts/neurology/users/jheffernan'
sys.path.insert(0, os.path.join(jheffernan, 'repositories', 'ecp'))

from ecp.workflows.base import init_cleanprep_wf

REST_TASKS=[
    'rfMRI_REST1_AP',
    'rfMRI_REST1_PA',
    'rfMRI_REST2_AP',
    'rfMRI_REST2_PA',
    'rfMRI_REST3_AP',
    'rfMRI_REST3_PA',
    'rfMRI_REST4_AP',
    'rfMRI_REST4_PA',
]

def get_parser():
    """Define parse object"""
    parser = ArgumentParser(description='preps ECP resting data in HCP format'
                                        'for connectivity analysis')
    parser.add_argument('data_dir', action='store', help='the data directory')
    parser.add_argument('work_dir', action='store', help='the working directory')
    parser.add_argument('out_dir', action='store', help='the output directory')
    parser.add_argument('--participants', action='store', nargs='+',
                        help='participants to be prepped')
    parser.add_argument('--tasks', action='store', nargs='+', default=REST_TASKS,
                        help='the task names to prep')
    parser.add_argument('--skip-vols', action='store', type=int, default=None,
                        help='the number of non-steady state volumes')
    parser.add_argument('--n-procs', action='store', type=int, default=None,
                        help='number of processors to use')

    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()

    n_procs = args.n_procs if args.n_procs else cpu_count()

    wfs = []
    for participant in args.participants:
        print('Building workflow for participant: {}'.format(participant))

        wfs.append(init_cleanprep_wf(
            args.data_dir, args.work_dir, args.out_dir, participant, args.tasks,
            args.skip_vols))

    for wf in wfs:
        wf.run(plugin='MultiProc', plugin_args={'n_procs' : n_procs})

    return 0

if __name__ == '__main__':
    main()

        
