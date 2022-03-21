import unittest
import subprocess
import pandas as pd
from snipar.imputation.preprocess_data import create_pedigree, add_control
import networkx as nx
import os
import numpy as np
from snipar.tests.utils import *

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


class TestPedigree(SniparTest):
    def test_create_pedigree(self):
        result = create_pedigree(f"{tests_root}/test_data/pedigree_creation_sample.king",
                                   f"{tests_root}/test_data/pedigree_creation_sample.agesex",
        ).sort_values(by=['FID', "IID"])
        expected = pd.read_csv(f"{tests_root}/test_data/pedigree_creation_sample.ped", delim_whitespace=True).sort_values(by=['FID', "IID"])
        result_graph = create_graph_from_pedigree(result)
        expected_graph = create_graph_from_pedigree(expected)

        equality = nx.is_isomorphic(result_graph, expected_graph, node_match)
        self.assertTrue(equality)
    
    def test_add_control(self):
        pedigree = pd.read_csv(f"{tests_root}/test_data/pedigree_creation_sample.ped", delim_whitespace=True).sort_values(by=['FID', "IID"]).astype(str)
        pedigree["has_father"] = pedigree["FATHER_ID"].isin(pedigree["IID"])
        pedigree["has_mother"] = pedigree["MOTHER_ID"].isin(pedigree["IID"])
        controlled_pedigree = add_control(pedigree).astype(str)
        controlled_fams = [fid[3:] for fid in controlled_pedigree["FID"] if fid.startswith("_")]
        has_father = pedigree["FATHER_ID"].isin(pedigree["IID"])
        has_mother = pedigree["MOTHER_ID"].isin(pedigree["IID"])
        sibs = pedigree[has_father & has_mother].groupby(["FID", "FATHER_ID", "MOTHER_ID"]).agg({"IID":lambda x:len(list(x))})
        sibs = sibs.groupby("FID").agg({"IID":max}).reset_index()
        no_parent_control = sibs[sibs["IID"]>1]["FID"].values.tolist()
        one_parent_control = pedigree[pedigree["has_father"] & pedigree["has_mother"]]["FID"].values.tolist()
        expected_controlls = no_parent_control + one_parent_control
        self.assertEqual(set(controlled_fams), set(expected_controlls))
        controlled_pedigree = controlled_pedigree[controlled_pedigree["FID"].str.startswith("_")]
        controlled_pedigree["FID"] = controlled_pedigree["FID"].str[3:]
        merged = pd.merge(pedigree, controlled_pedigree, on=["FID", "IID"])
        self.assertEqual(merged.shape[0], controlled_pedigree.shape[0])
            
