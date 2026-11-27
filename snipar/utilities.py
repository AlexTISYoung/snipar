import errno
import hashlib
import os
from tempfile import TemporaryDirectory
from typing import Any, Dict, Optional, Union
from pathlib import Path
import numpy as np
import pandas as pd
from os import path
from typing import Tuple

from snipar.docgen import get_parser_doc

def make_id_dict(x,col=0):
    """
    Make a dictionary that maps from the values in the given column (col) to their row-index in the input array
    """
    if len(x.shape)>1:
        x = x[:,col]
    id_dict = {}
    for i in range(0,x.shape[0]):
        id_dict[x[i]] = i
    return id_dict

def convert_str_array(x):
    """
    Convert an ascii array to unicode array (UTF-8)
    """
    x = np.array(x)
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.decode('UTF-8') for y in x])
    return x_out.reshape(x_shape)

def encode_str_array(x):
    """
    Encode a unicode array as an ascii array
    """
    x = np.array(x)
    x_shape = x.shape
    x = x.flatten()
    x_out = np.array([y.encode('ascii') for y in x])
    return x_out.reshape(x_shape)

def parse_obsfiles(obsfiles, obsformat='bed', append = True, wildcard = '@', chromosomes=None):
    obs_files = []
    chroms = []
    if wildcard in obsfiles:
        bed_ixes = obsfiles.split(wildcard)
        if chromosomes is None:
            chromosomes = np.arange(1,23)
        for i in chromosomes:
            obsfile = bed_ixes[0]+str(i)+bed_ixes[1]+'.'+obsformat
            if path.exists(obsfile):
                if append:
                    obs_files.append(obsfile)
                else:
                    obs_files.append(obsfile[:-(len(obsformat)+1)])
                chroms.append(i)
        if len(obs_files)==0:
            raise(ValueError('Observed genotype files not found'))
        else:
            print(str(len(obs_files))+' files found')
    else:
            obsfile = obsfiles+'.'+obsformat
            if path.exists(obsfile):
                obs_files = [obsfile]
                chroms = [0]
            else:
                raise(ValueError(obsfile+' not found'))
    return np.array(obs_files), np.array(chroms,dtype=int)

def parse_filelist(obsfiles, impfiles, obsformat, chromosomes=None):
    obs_files = []
    imp_files = []
    chroms = []
    if '@' in obsfiles and impfiles:
        bed_ixes = obsfiles.split('@')
        imp_ixes = impfiles.split('@')
        if chromosomes is None:
            chromosomes = np.arange(1,23)
        for i in chromosomes:
            obsfile = bed_ixes[0]+str(i)+bed_ixes[1]+'.'+obsformat
            impfile = imp_ixes[0]+str(i)+imp_ixes[1]+'.hdf5'
            if path.exists(impfile) and path.exists(obsfile):
                obs_files.append(obsfile)
                imp_files.append(impfile)
                chroms.append(i)
        if len(imp_files)==0:
            raise(ValueError('Observed/imputed genotype files not found'))
        else:
            print(str(len(imp_files))+' matched observed and imputed genotype files found')
    else:
            obsfile = obsfiles+'.'+obsformat
            impfile = impfiles+'.hdf5'
            if path.exists(obsfile) and path.exists(impfile):
                obs_files = [obsfile]
                imp_files = [impfile]
                chroms = [0]
            else:
                if not path.exists(obsfile):
                    raise(ValueError(obsfile+' not found'))
                if not path.exists(impfile):
                    raise(ValueError(impfile+' not found'))
    return np.array(obs_files), np.array(imp_files), np.array(chroms,dtype=int)

def outfile_name(outprefix,outsuffix,chrom=None):
    if '@' in outprefix:
        if chrom is None:
            raise(ValueError('Must provide chromosome number with wildcard character'))
        outprefix = outprefix.split('@')
        outprefix = outprefix[0]+str(chrom)+outprefix[1]
        return outprefix+outsuffix
    elif chrom is not None:
        return outprefix+'_chr_'+str(chrom)+outsuffix
    else:
        return outprefix+outsuffix

def coord2linear(ind1: int, ind2: int) -> int:
    row_ind, col_ind = max(ind1, ind2), min(ind1, ind2)
    return int(row_ind * (row_ind + 1) / 2 + col_ind)


def linear2coord(ind: int) -> Tuple[int, int]:
    ind += 1
    ind1 = np.ceil((2 * ind + 0.25) ** 0.5 - 0.5)
    ind2 = ind - (ind1 - 1) * ind1 / 2
    return int(ind1 - 1), int(ind2 - 1)

def read_oxford_sample(filename):
    """Read individual IDs, preserving leading zeros and skipping the column-type row."""
    samples = pd.read_csv(filename, sep=r"\s+", dtype=str, skiprows=[1], keep_default_na=False)
    return samples["ID_2"].to_numpy()


def _legacy_bgen(filename, *, lazy=False, allow_complex=False, verbose=False):
    """Use bgen-reader's alternate index path when the adjacent index cannot be written."""
    import bgen_reader
    from bgen_reader._environment import BGEN_READER_CACHE_HOME

    source = Path(filename).resolve()
    suffix = '.complex.metadata2.mmm' if allow_complex else '.metadata2.mmm'
    if lazy:
        suffix = '.metafile'

    def read(index):
        if lazy:
            return bgen_reader.read_bgen(source, metafile_filepath=index, verbose=verbose)
        return bgen_reader.open_bgen(source, metadata_filepath=index, allow_complex=allow_complex, verbose=verbose)

    stat = source.stat()

    def current(index):
        try:
            return index.stat().st_mtime_ns >= stat.st_mtime_ns
        except FileNotFoundError:
            return False

    def open_at(index):
        if current(index):
            return read(index)
        # bgen-reader uses a fixed temporary filename. Publish a completed
        # index from a private directory when multiple workers build at once.
        with TemporaryDirectory(dir=index.parent) as directory:
            temporary = Path(directory) / 'index'
            if lazy:
                bgen_reader.create_metafile(source, temporary, verbose=verbose)
            else:
                read(temporary).close()
            if temporary.stat().st_mtime_ns < stat.st_mtime_ns:
                os.utime(temporary, ns=(stat.st_mtime_ns, stat.st_mtime_ns))
            try:
                temporary.replace(index)
            except PermissionError:
                # Windows can keep an index published by another worker mapped.
                if not current(index):
                    raise
        if lazy:
            # Mirror bgen-reader's invalidation after replacing a .metafile.
            from bgen_reader import _genotype, _metafile
            for module in (_genotype, _metafile):
                with module.lock:
                    module.cache.clear()
        return read(index)

    try:
        return open_at(source.with_name(source.name + suffix))
    except OSError as error:
        if error.errno not in (errno.EACCES, errno.EPERM, errno.EROFS):
            raise

    # Version the key so an updated input never replaces another reader's
    # mapped index. Full paths distinguish identically named chromosomes.
    identity = os.fsencode(source) + f'\0{stat.st_size}:{stat.st_mtime_ns}'.encode()
    key = hashlib.sha256(identity).hexdigest()
    return open_at(BGEN_READER_CACHE_HOME.resolve() / 'metafile' / (key + suffix))


def _checked_bgen_sample(filename, nsamples):
    samples = read_oxford_sample(filename)
    if len(samples) != nsamples:
        raise ValueError('sample file length and bgen file length do not match.')
    return samples


class _bgen:
    """Simple wrapper class of bgen_reader.open_bgen that controls access to sample ids."""
    def __init__(self, 
                 filepath: Union[str, Path],
                 samples_filepath: Optional[Union[str, Path]] = None,
                 allow_complex: bool = False,
                 verbose: bool = True, reader=None):
        self._bgen = reader if reader is not None else _legacy_bgen(filepath, verbose=verbose, allow_complex=allow_complex)
        try:
            self._samples = _checked_bgen_sample(samples_filepath, self._bgen.nsamples)
        except BaseException:
            self._bgen.close()
            raise
    
    @property
    def samples(self) -> np.ndarray:
        return self._samples
    
    def __getattr__(self, name: str) -> Any:
        """Access to attributes other than sample ids."""
        return getattr(self._bgen, name)
    
    def __enter__(self):
        return self

    def __exit__(self, *args):
        self._bgen.close()


def open_bgen(filename: str, verbose: bool = False):
    """Wrapper of bgen_reader.open_bgen that checks if sample ids make sense."""
    try:
        import bgen_reader
    except ModuleNotFoundError as exc:
        if exc.name != "bgen_reader":
            raise
        from snipar.bgen_io import BgenFile
        return BgenFile(filename)
    bgen = _legacy_bgen(filename, verbose=verbose)
    if bgen.samples[0] == 'sample_0':
        print('WARNING: Sample ids in bgen file are generic. Trying to read the corresponding .sample file ...')
        samples_filepath = Path(filename).with_suffix('.sample')
        bgen = _bgen(filename, samples_filepath=samples_filepath, reader=bgen)
    return bgen


def read_bgen(filename: str, verbose: bool = False) -> Dict:
    """Wrapper of bgen_reader.read_bgen that checks if sample ids make sense."""
    try:
        import bgen_reader
    except ModuleNotFoundError as exc:
        if exc.name != "bgen_reader":
            raise
        raise ImportError("The lazy read_bgen API requires snipar[legacy-bgen]. "
                          "Use open_bgen for array access.") from exc
    bgen = _legacy_bgen(filename, lazy=True, verbose=verbose)
    if bgen['samples'][0] == 'sample_0':
        print('WARNING: Sample ids in bgen file are generic. Trying to read the corresponding .sample file ...')
        samples_filepath = Path(filename).with_suffix('.sample')
        if not path.exists(samples_filepath):
            raise FileNotFoundError(f'{samples_filepath} does not exist.')
        bgen['samples'] = _checked_bgen_sample(samples_filepath, len(bgen['samples']))
    return bgen
