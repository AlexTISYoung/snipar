import numpy as np
from pysnptools.snpreader import SnpReader
from pysnptools.pstreader import _MergeCols

class _MergeSIDs(_MergeCols,SnpReader):
    def __init__(self,  *args, **kwargs):
        super(_MergeSIDs, self).__init__(*args, **kwargs)

    def _savez(self, cache_file):
        np.savez(cache_file,
                    _row=np.array(self._row,dtype='S'), _row_property=self._row_property,
                    _col=np.array(self._col,dtype='S'), _col_property=self._col_property,
                    sid_count_list=self.sid_count_list)

    def _load(self,cache_file):
            with np.load(cache_file,allow_pickle=True) as data:
                self._col = np.array(data['_col'],dtype='str')
                self._col_property = data['_col_property']
                self.sid_count_list = data['sid_count_list']
                assert ('_row' in data) == ('_row_property' in data)
                self._row = np.array(data['_row'],dtype='str')
                self._row_property = data['_row_property']
