from pathlib import Path

import gzip
import h5py
import numpy as np
from pysnptools.snpreader import Bed, SnpData

from snipar import ibd
from snipar.scripts.gts_matrix import gts_matrix_chromosome


def _observations(tmp_path):
    ids = np.array([f's{i:03}' for i in range(48)])
    values = np.random.default_rng(17).integers(0, 3, (len(ids), 512)).astype(float)
    values[10, 5] = np.nan
    stem = tmp_path / 'one.bed'
    Bed.write(str(stem), SnpData(iid=np.column_stack((ids, ids)),
                                sid=np.array([f'rs{i:03}' for i in range(512)]), val=values,
                                pos=np.column_stack((np.ones(512), np.arange(512) * .01,
                                                     np.arange(1, 513)))), count_A1=True)
    return stem, ids


def test_ibd_disk_backed_tiles_and_raw_bed_agree_with_in_memory(tmp_path, monkeypatch):
    bedfile, ids = _observations(tmp_path)
    pairs = ids.reshape(-1, 2)
    plan = ibd.plan_batches

    def capped(n_items, bytes_per_item, *args, **kwargs):
        if n_items == len(pairs):
            kwargs['batch_size'] = 3
        return plan(n_items, bytes_per_item, *args, **kwargs)

    monkeypatch.setattr(ibd, 'plan_batches', capped)
    versions = ('memory', 'disk', 'packed_bed')
    for version in versions:
        np.random.seed(72)
        monkeypatch.setattr(ibd, 'get_memory_budget', lambda: 1_000_000 if version == 'memory' else 5000)
        if version == 'packed_bed':
            # Enough disk to spill only the LD sample, forcing the raw BED path.
            monkeypatch.setattr(ibd, 'disk_usage', lambda path: type('Disk', (), {'free': 100_000})())
        ibd.infer_ibd_chr(pairs, .001, None, str(tmp_path / version), bedfile=bedfile,
                          chrom=1, ibdmatrix=True, ld_out=True, num_threads=2)
    for kind in ('ibd.segments', 'ibdmatrix', 'l2.ldscore'):
        contents = [gzip.decompress(Path(tmp_path / f'{v}_chr_1.{kind}.gz').read_bytes()) for v in versions]
        assert contents[0] == contents[1] == contents[2]


def test_genotype_export_streaming_keeps_rows_variants_and_filtering(tmp_path):
    bedfile, ids = _observations(tmp_path)
    children = ids.reshape(-1, 4)[:, :2].ravel()
    families = np.repeat(np.array([f'f{i:02}' for i in range(len(ids) // 4)]), 2)
    parents = np.repeat(ids.reshape(-1, 4)[:, 2:], 2, axis=0)
    pedigree = np.column_stack((families, children, parents))
    for label, batch in (('full', 512), ('streamed', 7)):
        gts_matrix_chromosome(0, pedigree, str(tmp_path / label), bedfile=bedfile,
                              max_missing=0, min_maf=.25, num_threads=2, batch_size=batch)
    with h5py.File(tmp_path / 'full.genoarray.hdf5') as full, h5py.File(tmp_path / 'streamed.genoarray.hdf5') as stream:
        assert 0 < stream['genotypes'].shape[2] < 512
        assert full.keys() == stream.keys()
        for key in full:
            np.testing.assert_array_equal(full[key][:], stream[key][:])
    assert not list(tmp_path.glob('.snipar-genotypes-*'))


def test_bed_ld_with_tile_halos_matches_full_calculation(tmp_path):
    from snipar.ld import _bounded_bed_ld
    from snipar.observed import open_observed

    bedfile, _ = _observations(tmp_path)
    columns = np.arange(0, 512, 3)
    genetic_map = columns * .01
    with open_observed(bedfile=bedfile, num_threads=2) as source:
        reference = _bounded_bed_ld(source, columns, genetic_map, 0.27)
        tiled = _bounded_bed_ld(source, columns, genetic_map, 0.27, batch_size=7)
    np.testing.assert_allclose(reference, tiled, atol=1e-12, rtol=1e-12)


def test_streaming_scores_keep_the_public_output_format(tmp_path):
    from snipar.pgs import pgarray

    ids = np.array(['001', '002'])
    families = np.array(['f1', 'f2'])
    parent_ids = np.array([['p1', 'm1'], ['p2', 'm2']])
    data = np.array([[.251, -.75, np.nan], [.42, -.33, .14]], dtype=np.float32)
    scores = pgarray(data, ids, sid=np.array(['proband', 'paternal', 'maternal']), fams=families,
                     ped=np.column_stack((families, ids, parent_ids)), par_status=np.array([[0, 1], [0, 0]]))
    expected = tmp_path / 'expected.txt'
    actual = tmp_path / 'actual.txt.gz'
    scores.write(expected)
    scores.write(actual, return_array=False)
    assert gzip.decompress(actual.read_bytes()) == expected.read_bytes()
