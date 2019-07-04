import os
import sys

from argparse import ArgumentParser

jheffernan = '/rcc/stor1/depts/neurology/users/jheffernan'
sys.path.insert(0, os.path.join(jheffernan, 'repositories', 'ecp'))

from ecp.interfaces.paths import PostFreeSurferFiles
from ecp.workflows.anat import init_hcp_segment_anat_wf

def get_parser():
    """Define parse object"""
    parser = ArgumentParser(description='segments anatomical image HCP style')
    parser.add_argument('data_dir', action='store', help='the data directory')
    parser.add_argument('work_dir', action='store', help='the working directory')
    parser.add_argument('--participants', action='store', nargs='+',
                        help='participants anatomy image to be segmented')

    return parser

def run_ecp_segment_anat(data_dir, work_dir, participants):
    wfs = []
    for par in participants:
        par_work_dir = os.path.join(work_dir, par)

        par_files = PostFreeSurferFiles(base_dir=data_dir, 
                                        subject=par).run().outputs

        segment_wf = init_hcp_segment_anat_wf(out_dir=par_files.roi_folder)
        segment_wf.base_dir = par_work_dir
        inputnode = segment_wf.inputs.inputnode
        
        inputnode.wmparc = par_files.wmparc

        inputnode.l_atlasroi = par_files.L_atlasroi_32k_fs_LR
        inputnode.l_midthickness = par_files.L_midthickness_32k_fs_LR
        inputnode.l_white = par_files.L_white_32k_fs_LR
        inputnode.l_pial = par_files.L_pial_32k_fs_LR

        inputnode.r_atlasroi = par_files.R_atlasroi_32k_fs_LR
        inputnode.r_midthickness = par_files.R_midthickness_32k_fs_LR
        inputnode.r_white = par_files.R_white_32k_fs_LR
        inputnode.r_pial = par_files.R_pial_32k_fs_LR

        inputnode.ROIs = par_files.subcortical

        wfs.append(segment_wf)

    for wf in wfs:
        wf.run()

def main():
    parser = get_parser()
    args = parser.parse_args()
    run_ecp_segment_anat(args.data_dir, args.work_dir, args.participants)
    return 0

if __name__ == "__main__":
    main()
