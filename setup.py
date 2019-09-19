from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='ecp',
    version='0.0.1',
    author='Joe Heffernan',
    author_email='jheffernan@mcw.edu',
    license='MIT',
    description='A package to process ECP resting data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/heffjos/ecp',
    packages=find_packages(),
    entry_points={
        'gui_scripts': [
            'ecp_cleanprep=ecp.cli.ecp_cleanprep:main',
            'ecp_clean=ecp.cli.ecp_clean:main',
        ],
    },
    install_requires=[
        'numpy',
        'pandas',
        'importlib_resources',
        'fmriprep',
        'nipype',
        'nibabel',
    ],
    package_data={
        '': ['LICENSE', 'README.md'],
        'ecp': ['data/FreeSurferCSFRegLut.txt',
                'data/FreeSurferWMRegLut.txt',
                'data/glasser_conte.dlabel.nii',
                'data/fsl_identity.txt',
                'data/itk_identity.txt',],
    },
    include_package_data=True,
)

