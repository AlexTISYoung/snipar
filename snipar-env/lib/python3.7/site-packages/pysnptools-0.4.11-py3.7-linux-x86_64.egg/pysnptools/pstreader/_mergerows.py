#!!! Has not been tested

import logging
import os
import numpy as np
from pysnptools.pstreader import PstReader

class _MergeRows(PstReader): #!!!why does this start with _
    def __init__(self,reader_list,cache_file=None,skip_check=False):
        super(_MergeRows, self).__init__()
        assert len(reader_list) > 0, "Expect at least one reader"
        self.skip_check = skip_check

        self.reader_list = reader_list
        self._repr_string = "_MergeRows({0})".format(",".join([str(s) for s in reader_list]))

        if cache_file is not None:
            if not os.path.exists(cache_file):
                self._run_once()
                np.savez(cache_file, _row=self._row, _row_property=self._row_property, _row_count_list=self._row_count_list)
            else:
                with np.load(cache_file,allow_pickle=True) as data:
                    self._row = data['_row']
                    self._row_property = data['_row_property']
                    self._row_count_list = data['_row_count_list']
                    self._has_run_once = True


    def _run_once(self):
        if hasattr(self,'_has_run_once'):
            return
        self._has_run_once = True
        #Check that all iids are distinct and that all sids and pos are the same and in the same order

        row_list = []
        row_property_list = []
        row_set = set()
        col = self.reader_list[0].col
        col_property = self.reader_list[0].col_property
        for reader_index,reader in enumerate(self.reader_list):
            if reader_index % 10 == 0: logging.info("_MergeRows looking at reader #{0}: {1}".format(reader_index,reader))
            if not self.skip_check:
                assert np.array_equal(col,reader.col), "Expect columns to be the same across all files"
                np.testing.assert_equal(col_property,reader.col_property) #"Expect column_property to be the same across all files"
                size_before = len(row_set)
                row_set.update((tuple(item) for item in reader.row))
                assert len(row_set) == size_before + reader.row_count, "Expect rows to be distinct in all files"
            row_list.append(reader.row)
            row_property_list.append(reader.row_property)
        self._row = np.concatenate(row_list)
        self._row_property = np.concatenate(row_property_list)
        self._row_count_list = [len(row) for row in row_list]

    @property
    def row(self):
        self._run_once()
        return self._row

    @property
    def col(self):
        self._run_once()
        return self.reader_list[0].col

    @property
    def row_property(self):
        self._run_once()
        return self._row_property

    @property
    def col_property(self):
        self._run_once()
        return self.reader_list[0].col_property

    def __repr__(self): 
        #Don't need _run_once because based only on initial info
        return self._repr_string

    def copyinputs(self, copier):
        self._run_once()
        for reader in self.reader_list:
            copier.input(reader)

    def _create_reader_and_iid_index_list(self,iid_index):
        result = []
        start = 0
        for reader_index in range(len(self.reader_list)): #!!!this needs test coverage
            stop = start + self._row_count_list[reader_index]
            is_here = (iid_index >= start) * (iid_index < stop)
            if any(is_here):
                iid_index_rel = iid_index[is_here]-start
                result.append((reader_index,is_here,iid_index_rel))
            start = stop
        return result

    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
        #!!!tests to do: no iid's
        #!!!tests to do: no sid's
        #!!!test to do: from file 1, file2, and then file1 again

        iid_index = iid_index_or_none if iid_index_or_none is not None else np.arange(self.iid_count) #!!!might want to special case reading all
        sid_index_or_none_count = self.sid_count if sid_index_or_none is None else len(sid_index_or_none)
        reader_and_iid_index_list = self._create_reader_and_iid_index_list(iid_index)

        if len(reader_and_iid_index_list) == 0:
            return self.reader_list[0]._read(iid_index,sid_index_or_none,order,dtype,force_python_only, view_ok)
        elif len(reader_and_iid_index_list) == 1:
            reader_index,iid_index_in,iid_index_rel = reader_and_iid_index_list[0]
            reader = self.reader_list[reader_index]
            return reader._read(iid_index_rel,sid_index_or_none,order,dtype,force_python_only, view_ok)
        else:
            logging.info("Starting read from {0} subreaders".format(len(reader_and_iid_index_list)))
            if order == 'A' or order is None:
                order = 'F'
            val = np.empty((len(iid_index),sid_index_or_none_count),dtype=dtype,order=order)
            for reader_index,is_here,iid_index_rel in reader_and_iid_index_list:
                reader = self.reader_list[reader_index]
                if reader_index % 1 == 0: logging.info("Reading from #{0}: {1}".format(reader_index,reader))
                val[is_here,:] = reader._read(iid_index_rel,sid_index_or_none,order,dtype,force_python_only, view_ok=True)
            logging.info("Ended read from {0} subreaders".format(len(reader_and_iid_index_list)))
            return val

