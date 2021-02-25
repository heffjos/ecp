import os
import sys
import subprocess

import pandas as pd

from shlex import split
from pathlib import Path
from argparse import ArgumentParser

from nipype import Workflow
from ecp.workflows.base import init_cleanprep_wf

def get_parser():
    """Define parse object"""
    parser = ArgumentParser(description='preps ECP resting data in HCP format'
                                        'for connectivity analysis')
    parser.add_argument('data_dir', action='store', help='the data directory')
    parser.add_argument('work_dir', action='store', help='the working directory')
    parser.add_argument('out_dir', action='store', help='the output directory')
    parser.add_argument('spec_file', action='store', help='csv spec file')
    parser.add_argument('participant', action='store', help='prep participant')
    parser.add_argument('--n-procs', action='store', type=int, default=8,
                        help='number of processors to use')

    return parser

def run_cleanprep_wf(args):

    data_dir = Path(args.data_dir).resolve()
    work_dir = Path(args.work_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    spec_file = Path(args.spec_file).resolve().as_posix()

    participant = args.participant
    n_procs = args.n_procs

    spec = pd.read_csv(spec_file, sep='\t')
    spec = spec.loc[spec.usable_data, :]

    if participant not in set(spec.subjects):
        raise Exception(f'Subject not found in spec file: {participant}')

    subject_work_dir = work_dir.joinpath(f'{participant}_cleanprep')
    participant_info = spec[spec.subjects == participant]
    tasks = list(participant_info['func'])
    skip_begin = list(participant_info['movement_skip_begin'])
    skip_end = list(participant_info['movement_skip_end'])
    
    cleanprep_wf = init_cleanprep_wf(
        data_dir.as_posix(), 
        subject_work_dir.as_posix(), 
        out_dir.as_posix(), 
        participant, 
        tasks, 
        skip_begin, 
        skip_end)
    
    cleanprep_wf.run(plugin='MultiProc', plugin_args={'n_procs': n_procs})

    return 0

def main():
    parser = get_parser()
    args = parser.parse_args()

    return run_cleanprep_wf(args)

if __name__ == '__main__':
    main()

