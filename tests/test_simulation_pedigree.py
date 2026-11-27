"""The parents recorded for each mating must carry their own phenotypes."""
import numpy as np
import pytest

from snipar.simulate import forward_sim, simulate_first_gen


@pytest.mark.parametrize('n_random, n_am', [(0, 1), (1, 2)])
def test_simulated_parent_phenotypes_follow_parent_ids(n_random, n_am):
    n_families = 512
    haps, maps, *_ = simulate_first_gen(n_families, 32, maf=0.3)
    result = forward_sim(haps, maps, n_random, n_am, unlinked=True,
                         n_causal=32, h2=0.5, v_indirect=0.5,
                         r_direct_indirect=0, r_y=0.9)
    header, *records = result[5]
    columns = {label: i for i, label in enumerate(header)}
    by_id = {row[columns['IID']]: row for row in records}
    for row in records[2 * n_families:]:
        for parent in ('FATHER', 'MOTHER'):
            linked = by_id[row[columns[parent + '_ID']]]
            for trait in ('PHENO', 'DIRECT', 'INDIRECT'):
                assert row[columns[parent + '_' + trait]] == linked[columns[trait]]

    pairs = np.array(records[-2 * n_families::2])
    parent_phenotypes = pairs[:, [columns['FATHER_PHENO'], columns['MOTHER_PHENO']]].astype(float)
    assert np.corrcoef(parent_phenotypes.T)[0, 1] > 0.8
