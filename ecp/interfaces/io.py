import os

from nipype.interfaces.base import (
    traits, isdefined, Undefined,
    TraitedSpec, BaseInterfaceInputSpec, DynamicTraitedSpec,
    File, Directory, InputMultiObject, OutputMultiObject, Str,
    SimpleInterface,
)
from nipype import config
from nipype.utils.misc import str2bool
from nipype.utils.filemanip import ensure_list, copyfile
from nipype.interfaces.io import IOBase

# DataSink inputs
class DirectoryDataSinkInputSpec(DynamicTraitedSpec, BaseInterfaceInputSpec):
    '''
    '''

    # Init inputspec data attributes
    base_directory = Directory(
        desc='Path to the base directory for storing data.',
        madatory=True,
        exists=True)
    _outputs = traits.Dict(Str, value={}, usedefault=True)
    remove_dest_dir = traits.Bool(
        False, usedefault=True, desc='remove dest directory when copying dirs')

    # Set call-able inputs attributes
    def __setattr__(self, key, value):

        if key not in self.copyable_trait_names():
            if not isdefined(value):
                super(DirectoryDataSinkInputSpec, self).__setattr__(key, value)
            self._outputs[key] = value
        else:
            if key in self._outputs:
                self._outputs[key] = value
            super(DirectoryDataSinkInputSpec, self).__setattr__(key, value)


# DataSink outputs
class DirectoryDataSinkOutputSpec(TraitedSpec):

    # Init out file
    out_file = traits.Any(desc='datasink output')


# Custom DataSink class
class DirectoryDataSink(IOBase):
    """ Generic datasink module to store structured outputs in an existing
        directory

        Primarily for use within a workflow. This interface allows arbitrary
        creation of input attributes.
 
        the general form of the output is::

           'base_directory/filename'

           filename comesfrom the input to the connect statement.

        .. warning::

            This is not a thread-safe node because it can write to a common
            shared location. It will not complain when it overwrites a file.

        .. note::

            This interface **cannot** be used in a MapNode as the inputs are
            defined only when the connect statement is executed.

        Examples
        --------

        >>> ds = DataSink()
        >>> ds.inputs.base_directory = 'results_dir'
        >>> ds.inputs.structural = 'structural.nii'
        >>> setattr(ds.inputs, 'files1', ['cont1.nii', 'cont2.nii'])
        >>> setattr(ds.inputs, 'files2', ['cont1a.nii', 'cont2a.nii'])
        >>> ds.run()  # doctest: +SKIP

        To use DataSink in a MapNode, its inputs have to be defined at the
        time the interface is created.

        >>> ds = DataSink(infields=['contasts.@con'])
        >>> ds.inputs.base_directory = 'results_dir'
        >>> ds.inputs.structural = 'structural.nii'
        >>> setattr(ds.inputs, 'contrasts.@con', ['cont1.nii', 'cont2.nii'])
        >>> setattr(ds.inputs, 'contrasts.alt', ['cont1a.nii', 'cont2a.nii'])
        >>> ds.run()  # doctest: +SKIP

    """

    # Give obj .inputs and .outputs
    input_spec = DirectoryDataSinkInputSpec
    output_spec = DirectoryDataSinkOutputSpec

    # Initialization method to set up datasink
    def __init__(self, infields=None, force_run=True, **kwargs):
        """
        Parameters
        ----------
        infields : list of str
            Indicates the input fields to be dynamically created
        """

        super(DirectoryDataSink, self).__init__(**kwargs)
        undefined_traits = {}
        # used for mandatory inputs check
        self._infields = infields
        if infields:
            for key in infields:
                self.inputs.add_trait(key, traits.Any)
                self.inputs._outputs[key] = Undefined
                undefined_traits[key] = Undefined
        self.inputs.trait_set(trait_change_notify=False, **undefined_traits)
        if force_run:
            self._always_run = True

    # Get destination paths
    def _get_dst(self, src):
        # If path is directory with trailing os.path.sep,
        # then remove that for a more robust behavior
        src = src.rstrip(os.path.sep)
        path, fname = os.path.split(src)
        if fname:
            dst = fname
        else:
            dst = path.split(os.path.sep)[-1]
        if dst[0] == os.path.sep:
            dst = dst[1:]
        return dst

    # List outputs, main run routine
    def _list_outputs(self):
        """Execute this module.
        """

        # Init variables
        outputs = self.output_spec().get()
        out_files = []
        # Use hardlink
        use_hardlink = str2bool(
            config.get('execution', 'try_hard_link_datasink'))

        outdir = os.path.abspath(self.inputs.base_directory)

        # Iterate through outputs attributes {key : path(s)}
        for key, files in list(self.inputs._outputs.items()):
            if not isdefined(files):
                continue
            files = ensure_list(files)

            # flattening list
            if isinstance(files, list):
                if isinstance(files[0], list):
                    files = [item for sublist in files for item in sublist]

            # Iterate through passed-in source files
            for src in ensure_list(files):
                # Format src and dst files
                src = os.path.abspath(src)
                if not os.path.isfile(src):
                    src = os.path.join(src, '')
                dst = self._get_dst(src)
                dst = os.path.join(outdir, dst)

                # If src is a file, copy it to dst
                if os.path.isfile(src):
                    copyfile(
                        src,
                        dst,
                        copy=True,
                        hashmethod='content',
                        use_hardlink=use_hardlink)
                    out_files.append(dst)
                # If src is a directory, copy entire contents to dst dir
                elif os.path.isdir(src):
                    if os.path.exists(dst) and self.inputs.remove_dest_dir:
                        shutil.rmtree(dst)
                    copytree(src, dst)
                    out_files.append(dst)

        # Return outputs dictionary
        outputs['out_file'] = out_files

        return outputs
