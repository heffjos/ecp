import re
import os

from pathlib import Path

def get_entities(source_file):
    """Return the entities from a bids-like file name"""
    
    entities = {x: y for x, y in 
                [z.split('-') for z in re.findall('[a-zA-Z\d]+-[a-zA-Z\d]+', source_file)]}

    return entities

def generate_bold_name(sub, task, suffix, ext, 
                       session=None, acq=None, ce=None, 
                       dir=None, rec=None, run=None, echo=None,
                       space=None, desc=None, hemi=None):
    """generate bold file names for bids and derivatives"""

    subject = f'sub-{sub}'
    session = f'_ses-{session}' if session else ''
    task = f'_task-{task}'
    acq = f'_acq-{acq}' if acq else ''
    ce = f'_ce-{ce}' if ce else ''
    dir = f'_dir-{dir}' if dir else ''
    rec = f'_rec-{rec}' if rec else ''
    run = f'_run-{run}' if run else ''
    echo = f'_echo-{echo}' if echo else ''
    space = f'_space-{space}' if space else ''
    desc = f'_desc-{desc}' if desc else ''
    hemi = f'_hemi-{hemi}' if hemi else ''

    return (f'{subject}{session}{task}{acq}{ce}{dir}{rec}{run}{echo}{space}{desc}{hemi}'
            + '_' + suffix + ext)

def generate_bids_name(source_file, ext, desc=None, space=None, suffix=None):
    """Generates new bids file name based on source_file"""

    src_fname, _ = splitext(source_file)
    src_fname, dtype = src_fname.rsplit('_', 1)
    
    formatstr = '{bname}{space}{desc}{suffix}.{ext}'
    
    space = '_space-{}'.format(space) if space else ''
    desc = '_desc-{}'.format(desc) if desc else ''
    suffix = '_{}'.format(suffix) if suffix else ''
    
    out_file = formatstr.format(
        bname=src_fname,
        space=space,
        desc=desc,
        suffix=suffix,
        ext=ext,
    )

    return out_file

def splitext(fname):
    """Splits filename and extension (.gz safe)

    >>> splitext('some/file.nii.gz')
    ('file', '.nii.gz')
    >>> splitext('some/other/file.nii')
    ('file', '.nii')
    >>> splitext('otherext.tar.gz')
    ('otherext', '.tar.gz')
    >>> splitext('text.txt')
    ('text', '.txt')

    """
    basename = str(Path(fname).name)
    stem = Path(basename.rstrip('.gz')).stem
    return stem, basename[len(stem):]

def hcp_to_bids(hcp_name, subject, ses=None, suffix='bold'):
    """Convert hcp name to bids name"""

    tokenized = hcp_name.split('_')
    task, run = re.match('([a-zA-Z]+)(\d+)', tokenized[1]).groups()
    task = task.lower()
    direction = tokenized[2].lower()

    if ses:
        template = f'sub-{subject}_ses-{ses}_task-{task}_dir-{direction}_run-{run}_{suffix}.nii.gz'
    else:
        template = f'sub-{subject}_task-{task}_dir-{direction}_run-{run}_{suffix}.nii.gz'

    return template


def makedir(path):
    """make directory and not error out if it exists already"""

    try:
        os.mkdir(path, mode=0o775)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass
    

