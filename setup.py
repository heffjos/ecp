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
            'ecp_bidsify_hcp=ecp.cli.ecp_bidsify_hcp:main',
            'ecp_multitimeseries=ecp.cli.ecp_multitimeseries:main',
            'ecp_clean=ecp.cli.ecp_clean:main',
        ],
    },
    install_requires=[
        'numpy',
        'pandas',
        'importlib_resources',
        'fmriprep==20.0.5',
        'nipype',
        'nibabel',
    ],
    package_data={
        '': ['LICENSE', 'README.md'],
        'ecp': ['data/FreeSurferCSFRegLut.txt',
                'data/FreeSurferWMRegLut.txt',
                'data/glasser_conte.dlabel.nii',
                'data/fsl_identity.mat',
                'data/itk_identity.txt',],
    },
    include_package_data=True,
)

