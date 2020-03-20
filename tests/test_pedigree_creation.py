import unittest
import subprocess
import pandas as pd
from sibreg.bin.impute_runner import create_pedigree, add_control
import networkx as nx



def node_match(node1, node2):
    return node1["observed"] == node2["observed"]

def create_graph_from_pedigree(pedigree):
    fam_tree = nx.DiGraph()
    observed = {}
    for index, row in pedigree.iterrows():
        iid = row["IID"]
        father = row["FATHER_ID"]
        mother = row["MOTHER_ID"]
        if iid not in fam_tree:
            fam_tree.add_node(iid)
        if father not in fam_tree:
            fam_tree.add_node(iid)
        if mother not in fam_tree:
            fam_tree.add_node(iid)
        observed[iid] = True
        observed[father] = observed.get(father, False)
        observed[mother] = observed.get(mother, False)
        fam_tree.add_edge(mother, iid)
        fam_tree.add_edge(father, iid)
    nx.set_node_attributes(fam_tree, observed, 'observed')

    return fam_tree


class TestPedigree(unittest.TestCase):
    def test_create_pedigree(self):
        result = create_pedigree("test_data/sample.king",
                                   "test_data/sample.agesex",
        ).sort_values(by=['FID', "IID"])
        expected = pd.read_csv("test_data/sample.ped", sep = " ").sort_values(by=['FID', "IID"])
        result_graph = create_graph_from_pedigree(result)
        expected_graph = create_graph_from_pedigree(expected)
        equality = nx.is_isomorphic(result_graph, expected_graph, node_match)
        self.assertTrue(equality)
