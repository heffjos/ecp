from nipype import Workflow

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl

def init_generate_boldmask(preproc):
    """
    Generates a mask for preprocesed bold image.

    We need this because the original mask is created from the anatomy. 
    Sometimes the bold does not cover the brain as much as the T1, so
    there can be many voxels set to 0 for an entie timecourse.

    **Paramters**

        preproc: str
            hcp preprocessed bold
    
    **Inputs**

        None

    **Outputs**

        bold_mask
            BOLD series mask
    """

    workflow = Workflow(name='generate_bolmask_wf')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['bold_mask']), name='outputnode')

    mean_image = pe.Node(fsl.maths.MeanImage(
        in_file=preproc, out_file='out.nii.gz'), name='mean_image')
    thr_mean = pe.Node(fsl.maths.Threshold(thresh=0, out_file='out.nii.gz', args='-bin'),
        name='thr_mean')

    workflow.connect([
        (mean_image, thr_mean, [('out_file', 'in_file')]),
        (thr_mean, outputnode, [('out_file', 'bold_mask')]),
    ])

    return workflow
        
        

    

    

    
