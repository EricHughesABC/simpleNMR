# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:55:03 2021.

@author: ERIC
"""
from PyQt5.QtWidgets import QMessageBox

import pandas as pd
import numpy as np
from numpy import pi, sin, cos, exp
from scipy import stats
from scipy import fftpack
import os
import nmrglue as ng
import yaml
import re
import sys
import os
import json
from collections.abc import Iterable

import networkx as nx

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import PIL
from PIL import Image

excel_orig_df_columns_str = """{
    "h1": [
        "Name",
        "Shift",
        "Range",
        "H's",
        "Integral",
        "Class",
        "J's",
        "Method"
    ],
    "c13": [
        "ppm",
        "Intensity",
        "Width",
        "Area",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation"
    ],
    "pureshift": [
        "ppm",
        "Intensity",
        "Width",
        "Area",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation"
    ],
    "cosy": [
        "f2 (ppm)",
        "f1 (ppm)",
        "Intensity",
        "Width f2",
        "Width f1",
        "Volume",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation"
    ],
    "hsqc": [
        "f2 (ppm)",
        "f1 (ppm)",
        "Intensity",
        "Width f2",
        "Width f1",
        "Volume",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation"
    ],
    "hmbc": [        
        "f2 (ppm)",
        "f1 (ppm)",
        "Intensity",
        "Width f2",
        "Width f1",
        "Volume",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation"
    ]
}"""

new_dataframes = {}

excel_df_columns_str = """{
    "h1": [
        "ppm",
        "integral",
        "jCouplingClass",
        "jCouplingVals",
        "range",
        "label",
        "f1H_i",
        "f2H_i",
        "f1_i",
        "f2_i"
    ],
    "c13": [
        "ppm",
        "attached_protons",
        "ppmH1s",
        "max_bonds",
        "label",
        "f2C_i",
        "f1C_i",
        "f1_i",
        "f2_i"
    ],
    "pureshift": [
        "ppm",
        "intensity",
        "Width",
        "Area",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation"
    ],
    "cosy": [
        "f1_ppm",
        "f2_ppm",
        "intensity",
        "f1_i",
        "f1p_i",
        "f2_i",
        "f2p_i",
        "f1p_ppm",
        "f2p_ppm",
        "f1H_i",
        "f2H_i",
        "f1Cp_i",
        "f2Cp_i"
    ],
    "hsqc": [
        "f1_ppm",
        "f2_ppm",
        "intensity",
        "f1_i",
        "f2_i",
        "f2p_i",
        "f1C_i",
        "f2H_i",
        "f2Cp_i",
        "f2p_ppm"
    ],
    "hmbc": [
        "f1_ppm",
        "f2_ppm",
        "intensity",
        "f2p_ppm",
        "f1_i",
        "f2_i",
        "f2p_i",
        "f1C_i",
        "f2H_i",
        "f2Cp_i"
    ]
}"""

excel_df_columns = json.loads(excel_df_columns_str)
excel_orig_df_columns = json.loads(excel_orig_df_columns_str)


def read_in_cs_tables(
    h1: str, c13: str, scale_factor=6
) -> "list[pd.DataFrame, pd.DataFrame]":
    """

    Reads in json files for proton and chemical shift info and creates two
    pandas dataframes, one for proton and the other for carbon.


    Parameters.
    -----------
    h1 : str
        DESCRIPTION.
        filename of json file containing proton chemical shift information
    c13 : str
        DESCRIPTION.
        filename of json file containing proton chemical shift information
    scale_factor : TYPE, optional
        DESCRIPTION. The default is 6.

    Returns
    -------
    [pd.DataFrame, pd.DataFrame]

    """

    H1df = pd.read_json(h1)
    C13df = pd.read_json(c13)

    # create mean and sigma based on min max chemical shifts
    H1df["meanCS"] = (H1df.minCS + H1df.maxCS) / 2
    H1df["sigmaCS"] = (H1df.maxCS - H1df.minCS) / scale_factor

    C13df["meanCS"] = (C13df.minCS + C13df.maxCS) / 2
    C13df["sigmaCS"] = (C13df.maxCS - C13df.minCS) / scale_factor

    # create probability density functions for each chemical shitf group
    for i in H1df.index:
        H1df.loc[i, "norm"] = stats.norm(
            loc=H1df.loc[i, "meanCS"], scale=H1df.loc[i, "sigmaCS"]
        )

    for i in C13df.index:
        C13df.loc[i, "norm"] = stats.norm(
            loc=C13df.loc[i, "meanCS"], scale=C13df.loc[i, "sigmaCS"]
        )

    return H1df, C13df


def parse_argv(my_argv=None):

    if my_argv == None:
        my_argv = sys.argv

    datadirectories = None
    smilefilenames = None
    pngfilenames = None
    excelfilenames = None
    smiles_fn = None
    png_fn = None
    excel_fn = None

    datadirectories = [d for d in my_argv[1:] if os.path.isdir(d)]
    excelfilenames = [e for e in my_argv[1:] if e.endswith(".xlsx")]
    pngfilenames = [s for s in my_argv[1:] if s.endswith(".png")]
    smilefilenames = [s for s in my_argv[1:] if s.endswith(".smi")]
    xy3jsonfilenames = [s for s in my_argv[1:] if s == "xy3.json"]

    data_directory = "."
    if len(datadirectories) > 0:
        data_directory = datadirectories[0]

    if len(excelfilenames) > 0:
        if os.path.exists(os.path.join(data_directory, excelfilenames[0])):
            excel_fn = os.path.join(data_directory, excelfilenames[0])
        else:
            return {
                "data_directory": None,
                "excel_fn": None,
                "smiles_fn": None,
                "png_fn": None,
                "xy3_fn": None,
            }

    else:
        excelfilenames = [e for e in os.listdir(data_directory) if e.endswith(".xlsx")]
        if len(excelfilenames) > 0:
            excel_fn = os.path.join(data_directory, excelfilenames[0])
        else:
            # qt message box No excel file found in data directory

            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Information)
            msgBox.setText(
                "No excel file found in data directory\n{}".format(data_directory)
            )
            msgBox.setWindowTitle("Excel File Not Found")
            msgBox.setStandardButtons(QMessageBox.Ok)

            returnValue = msgBox.exec()

            return {
                "data_directory": os.getcwd(),
                "excel_fn": None,
                "smiles_fn": None,
                "png_fn": None,
                "xy3_fn": None,
            }
            # print("Program stopping")
            # sys.exit()

    if len(smilefilenames) > 0:
        if os.path.exists(os.path.join(data_directory, smilefilenames[0])):
            smiles_fn = os.path.join(data_directory, smilefilenames[0])
        else:
            smiles_fn = None
    else:
        smilefilenames = [s for s in os.listdir(data_directory) if s.endswith(".smi")]
        if len(smilefilenames) > 0:
            smiles_fn = os.path.join(data_directory, smilefilenames[0])
        else:
            smiles_fn = None

    if len(pngfilenames) > 0:
        if os.path.exists(os.path.join(data_directory, pngfilenames[0])):
            png_fn = os.path.join(data_directory, pngfilenames[0])
        else:
            png_fn = None
    else:
        pngfilenames = [s for s in os.listdir(data_directory) if s.endswith(".png")]
        if len(pngfilenames) > 0:
            png_fn = os.path.join(data_directory, pngfilenames[0])
        else:
            png_fn = None

    if len(xy3jsonfilenames) > 0:
        if os.path.exists(os.path.join(data_directory, xy3jsonfilenames[0])):
            xy3_fn = os.path.join(data_directory, xy3jsonfilenames[0])
        else:
            xy3_fn = None
    else:
        xy3jsonfilenames = [s for s in os.listdir(data_directory) if s == "xy3.json"]
        if len(xy3jsonfilenames) > 0:
            xy3_fn = os.path.join(data_directory, xy3jsonfilenames[0])
        else:
            xy3_fn = None

    print

    return {
        "data_directory": data_directory,
        "excel_fn": excel_fn,
        "smiles_fn": smiles_fn,
        "png_fn": png_fn,
        "xy3_fn": xy3_fn,
    }


def read_smiles_file(problemdata_info: dict) -> str:
    """Reads the smiles file and returns the smiles string

    Args:
        problemdata_info (dict): parsed arguments from command line

    Returns:
        str: smiles string or None
    """
    smiles_str = None
    if isinstance(problemdata_info["smiles_fn"], str):
        with open(problemdata_info["smiles_fn"], "r") as fp:
            smiles_str = fp.readline()

    return smiles_str


def read_png_file(problemdata_info: dict) -> PIL.Image.Image:
    """Reads the png file and returns the image object"""
    png = None
    if isinstance(problemdata_info["png_fn"], str):
        png = Image.open(problemdata_info["png_fn"])

    return png


def create_png_from_smiles(smiles_str: str) -> PIL.Image.Image:
    """Creates a png image from a smiles string via rdkit"""
    png = None
    molecule = Chem.MolFromSmiles(smiles_str)
    Draw.MolToFile(molecule, "molecule.png")
    png = Image.open("molecule.png")
    return png


def create_molecule_from_smiles(smiles_str: str) -> Chem.Mol:
    """Creates a RDKIT molecule from a smiles string via rdkit"""
    molecule = Chem.MolFromSmiles(smiles_str)
    molecule.Compute2DCoords()

    xy = [[xyz[0], xyz[1]] for xyz in molecule.GetConformer().GetPositions()]

    for a, (x, y) in zip(molecule.GetAtoms(), xy):
        a.SetDoubleProp("x", x)
        a.SetDoubleProp("y", y)

    # molecule.xy3 = {}
    # for i, n in enumerate(molecule.nodes()):
    #     molecule.xy3[n] = xy[i]
    return molecule


def read_xy3_jsonfile(problemdata_info: dict) -> dict:
    """Reads the xy3 json file and returns the json object"""
    xy3 = None
    print("problemdata_info[xy3_fn]", problemdata_info["xy3_fn"])
    if isinstance(problemdata_info["xy3_fn"], str):
        with open(problemdata_info["xy3_fn"], "r") as fp:
            xy3 = json.load(fp)

        for k, v in xy3.items():
            xy3[k] = np.asarray(v)

        print(xy3)

    return xy3


def center_molecule(nmrproblem, ax, pc_x, pc_y):

    xy3 = nmrproblem.xy3

    xy_np = np.array(list(xy3.values())).T
    xxx = xy_np[0]
    yyy = xy_np[1]

    delta_x = xxx.max() - xxx.min()
    delta_y = yyy.max() - yyy.min()

    xmin = xxx.min() - delta_x * pc_x
    xmax = xxx.max() + delta_x * pc_x
    ymin = yyy.min() - delta_y * pc_y
    ymax = yyy.max() + delta_y * pc_y

    # ax.set_xlim(xmax, xmin)
    # ax.set_ylim(ymin, ymax)
    return np.array([xmin, xmax, ymin, ymax])


def create_hmbc_edges_dict(nmrproblem):
    hmbc_edges = {}
    hmbc = nmrproblem.hmbc
    for i in hmbc.index:
        ci, cj = hmbc.loc[i, ["f1C_i", "f2Cp_i"]]
        if ci == None or cj == None:
            continue
        if ci in hmbc_edges:
            hmbc_edges[ci].add(cj)
        else:
            hmbc_edges[ci] = {cj}

        if cj in hmbc_edges:
            hmbc_edges[cj].add(ci)
        else:
            hmbc_edges[cj] = {ci}

    for c in nmrproblem.carbonAtoms:
        if c not in hmbc_edges:
            hmbc_edges[c] = set()

    return hmbc_edges


def build_model(nmrproblem):

    h1fn = r"csTables/h1_chemical_shift_table.jsn"
    c13fn = r"csTables/c13_chemical_shift_table.jsn"
    H1df_orig, C13df_orig = read_in_cs_tables(h1fn, c13fn)

    nmrproblem.H1df_orig = H1df_orig
    nmrproblem.C13df_orig = C13df_orig

    nmrproblem.calcProbDistFunctions(H1df_orig, C13df_orig)
    nmrproblem.identify1HC13peaks()
    nmrproblem.udic[0]["df"] = nmrproblem.H1df
    nmrproblem.udic[1]["df"] = nmrproblem.C13df
    nmrproblem.createInfoDataframes()
    nmrproblem.calculate1H13CSpectra1D()
    nmrproblem.save1DspecInfotoUdic()
    nmrproblem.create1H13Clabels(num_poss=3)


def build_molecule_graph_network(nmrproblem):
    """build  a graph of the molecule based on the COSY connections

    Args:
        nmrproblem (nmrProblem): the NMRproblem
    """

    df = nmrproblem.df
    catoms = nmrproblem.carbonAtoms
    hsqc = nmrproblem.hsqc
    c13 = nmrproblem.c13

    # create graph of molecule
    nmrproblem.molecule = nx.Graph()
    nmrproblem.molecule.add_nodes_from(catoms)

    molecule = nmrproblem.molecule

    # add some node information
    for n in catoms:
        molecule.nodes[n]["c13ppm"] = str(df.loc[n, "ppm"])
        molecule.nodes[n]["label"] = n

    # attach H1 chemical shifts to graph
    for n in molecule.nodes:
        # get list of labels of protons attached to carbon
        h1list = df.loc["hsqc", n]
        # convert to proton ppm
        molecule.nodes[n]["h1ppm"] = ",".join(
            [str(nmrproblem.H1labelH1ppm[hi]) for hi in h1list]
        )

    # add links to graph from cosy data
    for c in catoms:
        for l in df.loc["cosy", c]:
            # print(c,l)
            if (not molecule.has_edge(c, l)) and (c != l):
                molecule.add_edge(c, l)

    # add node color to graph
    # print(nmrproblem.h1)
    # print(nmrproblem.c13)
    # print(nmrproblem.hsqc)
    for i, n in enumerate(catoms):
        # print(n, i, c13.loc[i+1,"attached_protons"], nmrproblem.hmbc_edge_colors[c13.loc[i+1,"attached_protons"]])
        molecule.nodes[n]["node_color"] = nmrproblem.hmbc_edge_colors[
            c13.loc[i + 1, "attached_protons"]
        ]


def build_xy3_representation_of_molecule_from_smiles(nmrproblem):
    expected_molecule = nmrproblem.expected_molecule
    eigen_nodes = [
        a.GetSymbol() + str(a.GetIdx()) for a in expected_molecule.GetAtoms()
    ]
    eigen_carbons = [s for s in eigen_nodes if "C" in s]

    e_xy3 = {}
    xy3 = {}

    for n, (x, y, z) in zip(
        eigen_nodes, expected_molecule.GetConformer().GetPositions()
    ):

        if "C" in n:
            e_xy3[n] = np.array([x, y])

    for n1, n2 in zip(eigen_carbons, nmrproblem.carbonAtoms):
        # print("n1,n2", n1,n2)
        xy3[n2] = e_xy3[n1]

    nmrproblem.xy3 = xy3
    # print("nmrproblem.xy3", nmrproblem.xy3)


def build_xy3_representation_of_molecule(nmrproblem):
    """"""
    molecule = nmrproblem.molecule
    # create coodinates that look more molecule like
    atomicNumberfromSymbol = {"H": 1, "C": 6, "O": 8, "N": 7}

    for node in molecule.nodes:
        molecule.nodes[node]["atomic_num"] = atomicNumberfromSymbol[
            "".join([c for c in node if c.isalpha()])
        ]
        molecule.nodes[node][
            "chiral_tag"
        ] = rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED
        molecule.nodes[node]["formal_charge"] = 0
        molecule.nodes[node]["is_aromatic"] = False
        molecule.nodes[node]["hybridization"] = rdkit.Chem.rdchem.HybridizationType.SP3
        molecule.nodes[node]["num_explicit_hs"] = 0

    for e in molecule.edges:
        molecule.edges[e]["bond_type"] = rdkit.Chem.rdchem.BondType.SINGLE

    mmm = nx_to_mol(molecule)

    mmm.Compute2DCoords()

    xy = np.array([[xyz[0], xyz[1]] for xyz in mmm.GetConformer().GetPositions()])

    xy = xy.T
    xxx, yyy = xy

    yyy = yyy / yyy.max()
    yyy = yyy - yyy.mean()
    yyy = yyy / yyy.max()

    xxx = xxx / xxx.max()
    xxx = xxx - xxx.mean()
    xxx = xxx / xxx.max()

    xy = np.array([xxx, yyy])
    xy = xy.T

    # nmrproblem.xmin = -1.2
    # nmrproblem.xmax = 1.2

    # nmrproblem.ymin = -1.2
    # nmrproblem.ymax = 1.2

    print("nmrproblem.xy3", type(nmrproblem.xy3))

    if isinstance(nmrproblem.xy3, type(None)):
        xy3 = {}
        for i, n in enumerate(molecule.nodes()):
            xy3[n] = xy[i]

        nmrproblem.xy3 = xy3


def nx_to_mol(G: nx.Graph) -> rdkit.Chem.Mol:
    """convert a nx.Graph to a rdkit.Chem.Mol

    Args:
        G (nx.Graph): molecule netowrk graph derived from NMR data

    Returns:
        rdkit.Chem.Mol: rdkit chem Mol variable
    """

    atomicNumberfromSymbol = {"H": 1, "C": 6, "O": 8, "N": 7}

    mol = Chem.RWMol()
    atomic_nums = nx.get_node_attributes(G, "atomic_num")
    chiral_tags = nx.get_node_attributes(G, "chiral_tag")
    formal_charges = nx.get_node_attributes(G, "formal_charge")
    node_is_aromatics = nx.get_node_attributes(G, "is_aromatic")
    node_hybridizations = nx.get_node_attributes(G, "hybridization")
    num_explicit_hss = nx.get_node_attributes(G, "num_explicit_hs")
    node_to_idx = {}

    for node in G.nodes():
        a = Chem.Atom(atomic_nums[node])
        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(node_is_aromatics[node])
        a.SetHybridization(node_hybridizations[node])
        a.SetNumExplicitHs(num_explicit_hss[node])
        idx = mol.AddAtom(a)
        node_to_idx[node] = idx

    bond_types = nx.get_edge_attributes(G, "bond_type")

    for edge in G.edges():
        first, second = edge
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]
        bond_type = bond_types[first, second]
        mol.AddBond(ifirst, isecond, bond_type)

    Chem.SanitizeMol(mol)
    return mol


def center_molecule(nmrproblem, ax, pc_x, pc_y):

    xy3 = nmrproblem.xy3

    xy_np = np.array(list(xy3.values())).T
    xxx = xy_np[0]
    yyy = xy_np[1]

    delta_x = xxx.max() - xxx.min()
    delta_y = yyy.max() - yyy.min()

    xmin = xxx.min() - delta_x * pc_x
    xmax = xxx.max() + delta_x * pc_x
    ymin = yyy.min() - delta_y * pc_y
    ymax = yyy.max() + delta_y * pc_y

    # ax.set_xlim(xmax, xmin)
    # ax.set_ylim(ymin, ymax)
    return np.array([xmin, xmax, ymin, ymax])


def create_hmbc_graph_fragments(nmrproblem, hmbc_edges):
    hmbc_graphs = {}
    ntwk_labels = []

    catoms = nmrproblem.carbonAtoms
    df = nmrproblem.df
    xy3 = nmrproblem.xy3
    molecule = nmrproblem.molecule
    udic = nmrproblem.udic

    ret1 = None
    ret2 = None

    lineCollections = []
    for i, c in enumerate(catoms):

        if c not in hmbc_edges.keys():
            continue

        # create hmbc graph for node c and add  xy coodinates
        hmbc_graphs[c] = {}
        hmbc_graphs[c]["graph"] = nx.Graph()
        # print("[c] + list(hmbc_edges[c])\n", [c] + list(hmbc_edges[c]))
        hmbc_graphs[c]["xy"] = dict((k, xy3[k]) for k in [c] + list(hmbc_edges[c]) if k)
        hmbc_graphs[c]["colors"] = []

        # add nodes to hmbc graph
        hmbc_graphs[c]["graph"].add_nodes_from([c] + list(hmbc_edges[c]))

        # add edges
        for i, c1 in enumerate(hmbc_edges[c]):
            if c1:
                # print("c1", c1)
                hmbc_graphs[c]["graph"].add_edge(c, c1)
                hmbc_graphs[c]["colors"].append(nmrproblem.hmbc_edge_colors[i])

    return hmbc_graphs


class NMRproblem:
    def __init__(self, problemdata_info: dict, loadfromwhere=None, H1LarmorFreq=None):

        self.problemDirectoryPath = problemdata_info["data_directory"]
        self.problemDirectory = problemdata_info["data_directory"]
        self.rootDirectory, self.problemDirectory = os.path.split(
            problemdata_info["data_directory"]
        )

        self.hmbc_edge_colors = (
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        )

        self.xmin = -1.5
        self.ymin = -1.5
        self.xmax = 1.5
        self.ymax = 1.5

        self.data_complete = False

        self.pngFiles = []
        self.jupyterFiles = []
        self.fidFilesDirectories = []
        self.yamlFiles = []
        self.excelFiles = []
        self.csvFiles = []
        self.pklFiles = []

        self.df = pd.DataFrame()
        self.dfBackup = pd.DataFrame()
        self.dfColumns = None
        self.dfIndex = None

        self.moleculePNGpanel = None
        self.moleculeAtomsStr = ""

        self.protonAtoms = []
        self.carbonAtoms = []
        self.elements = {}
        self.iprobs = {}

        self.numProtonGroups = 0
        self.C = 0

        self.dbe = 0

        # pandas dataframes

        self.h1 = pd.DataFrame(columns=excel_df_columns["h1"])
        self.c13 = pd.DataFrame(columns=excel_df_columns["c13"])
        self.pureshift = pd.DataFrame(columns=excel_df_columns["pureshift"])
        self.cosy = pd.DataFrame(columns=excel_df_columns["cosy"])
        self.hsqc = pd.DataFrame(columns=excel_df_columns["hsqc"])
        self.hmbc = pd.DataFrame(columns=excel_df_columns["hmbc"])

        # graph coodinates of carbon skeleton molecule
        self.xy = None
        self.yx3 = None

        # png and smiles

        self.png = None
        self.smiles = None

        if isinstance(H1LarmorFreq, (int, float)):
            obs = H1LarmorFreq
        else:
            obs = 400.0

        self.udic = {
            "ndim": 2,
            0: {
                "obs": obs,
                "sw": 18 * obs,
                "dw": 1.0 / (18 * obs),
                "car": (((16.0 - (-2.0)) / 2) - 2.0) * obs,
                "size": int(32 * 1024),
                "label": "1H",
                "complex": True,
                "encoding": "direct",
                "time": False,
                "freq": True,
                "lb": 0.5,
            },
            1: {
                "obs": 100.0,
                "sw": 220 * obs / 4,
                "dw": 1.0 / (220 * obs / 4),
                "car": 220 * obs / 4.0 / 2.0,
                "size": int(1024 * 32),
                "label": "13C",
                "complex": True,
                "encoding": "direct",
                "time": False,
                "freq": True,
                "lb": 0.5,
            },
        }

        self.udic[0]["axis"] = ng.fileiobase.unit_conversion(
            self.udic[0]["size"],
            self.udic[0]["complex"],
            self.udic[0]["sw"],
            self.udic[0]["obs"],
            self.udic[0]["car"],
        )

        self.udic[1]["axis"] = ng.fileiobase.unit_conversion(
            self.udic[1]["size"],
            self.udic[1]["complex"],
            self.udic[1]["sw"],
            self.udic[1]["obs"],
            self.udic[1]["car"],
        )

        if isinstance(problemdata_info["excel_fn"], str):
            loadfromwhere = "xlsx"

        if isinstance(loadfromwhere, type(None)):

            if self.init_class_from_yml(self.problemDirectory):
                self.data_complete = True

            elif self.init_class_from_excel(self.problemDirectory):
                self.data_complete = True

        else:
            if loadfromwhere == "yml":
                if self.init_class_from_yml(self.problemDirectory):
                    self.data_complete = True
            elif loadfromwhere == "xlsx":
                if self.init_class_from_excel(
                    self.problemDirectory, problemdata_info["excel_fn"]
                ):
                    self.data_complete = True

        self.png = read_png_file(problemdata_info)
        self.smiles = read_smiles_file(problemdata_info)

        if isinstance(self.smiles, str):
            self.png = create_png_from_smiles(self.smiles)
            self.expected_molecule = create_molecule_from_smiles(self.smiles)
            # print("self.expected_molecule")
            # print(self.expected_molecule)

        self.xy3 = read_xy3_jsonfile(problemdata_info)

        # set xlims of 13C and 1H 1D spectra to +/- 10% of biggest and smallest ppm

        self.min_max_1D_ppm = []

        ppm_min = self.h1.ppm.min() - (self.h1.ppm.max() - self.h1.ppm.min()) / 10.0
        ppm_max = self.h1.ppm.max() + (self.h1.ppm.max() - self.h1.ppm.min()) / 10.0

        self.min_max_1D_ppm.append((ppm_max, ppm_min))

        ppm_min = self.c13.ppm.min() - (self.c13.ppm.max() - self.c13.ppm.min()) / 10.0
        ppm_max = self.c13.ppm.max() + (self.c13.ppm.max() - self.c13.ppm.min()) / 10.0

        self.min_max_1D_ppm.append((ppm_max, ppm_min))

    def init_class_from_excel(self, excelFileNameDirName: str, excel_fn=None):
        """
        read in class parameters from excel file found in problem directory if found

        Parameters
        ----------
        excelFileNameDirName : str
            DESCRIPTION. name and path to directory holding excel file

        Returns
        -------
        bool
            DESCRIPTION. Return True if excel found and processed

        """

        # Load the excel file if provided

        if isinstance(excel_fn, str):
            self.excelFiles = [excel_fn]

            # try to load excel file
            # display qt message box  if excel file not found or not readable
            # if excel file found, load it
            try:
                self.excelsheets = pd.read_excel(
                    self.excelFiles[0], sheet_name=None, index_col=0
                )
            except FileNotFoundError:
                # display qt message box  if excel file not found or not readable

                msgBox = QMessageBox(
                    "Excel file not found", "Excel file not found", QMessageBox.Ok
                )
                msgBox = QMessageBox()
                msgBox.setIcon(QMessageBox.Information)
                msgBox.setText("Excel File not Found\n{}".format(self.excelFiles[0]))
                msgBox.setWindowTitle("Excel File not Found")
                msgBox.setStandardButtons(QMessageBox.Ok)
                rtn = msgBox.exec_()
                return False
            except PermissionError:
                # display qt message box  if excel file not found or not readable

                msgBox = QMessageBox()
                msgBox.setIcon(QMessageBox.Information)
                msgBox.setText("Cannot Open\n{}".format(self.excelFiles[0]))
                msgBox.setWindowTitle("Excel File Access Error")
                msgBox.setStandardButtons(QMessageBox.Ok)
                rtn = msgBox.exec_()
                return False

        else:
            # Look for excel file in directory and use first one found

            if os.path.isdir(excelFileNameDirName):
                self.excelFiles = [
                    os.path.join(excelFileNameDirName, f)
                    for f in os.listdir(excelFileNameDirName)
                    if f.endswith("xlsx")
                ]
                if len(self.excelFiles) == 0:
                    return False
            elif excelFileNameDirName.endswith("xlsx") and os.path.exists(
                excelFileNameDirName
            ):
                self.excelFiles = [excelFileNameDirName]
            else:
                return False

            if not os.path.exists(self.excelFiles[0]):
                # print(self.excelFiles[0], "does not exist")
                return False

            self.excelsheets = pd.read_excel(
                self.excelFiles[0], sheet_name=None, index_col=0
            )

        # define a consistant set of column names for dataframes
        # For now we a re only working with excel sheets created from MestresNove
        for k, df in self.excelsheets.items():
            df.rename(
                columns={
                    "H's": "numProtons",
                    "Integral": "integral",
                    "J's": "jCouplingVals",
                    "Class": "jCouplingClass",
                    "Intensity": "intensity",
                    "Shift": "ppm",
                    "Range": "range",
                    "f2 (ppm)": "f2_ppm",
                    "f1 (ppm)": "f1_ppm",
                },
                inplace=True,
            )  # If they exist in dataframe!

        # replace any NaN values with 0
        for k, df in self.excelsheets.items():
            df.fillna(0, inplace=True)

        if "molecule" in self.excelsheets:
            self.molecule_defined = True
            self.moleculeAtomsStr = self.excelsheets["molecule"].molecule.values[0]
            self.smiles = self.excelsheets["molecule"].smiles.values[0]
            self.calculate_dbe()
        else:
            print("molecule not in sheet")
            self.dbe = 0
            self.elements = {"H": 0, "C": 0}
            self.molecule_defined = False

        # define short names for the dataframes
        self.h1_df = self.excelsheets["H1_1D"]
        self.cosy_df = self.excelsheets["COSY"]
        self.hsqc_df = self.excelsheets["HSQC"]
        self.hmbc_df = self.excelsheets["HMBC"]
        self.pureshift_df = self.excelsheets["H1_pureshift"]
        self.c13_df = self.excelsheets["C13_1D"]

        # sort h1, c13 and pureshift just in case they are out of order
        # reindex startting from 1

        self.h1_df = self.h1_df.sort_values("ppm", ascending=False, ignore_index=True)
        self.h1_df.index = self.h1_df.index + 1

        self.c13_df = self.c13_df.sort_values("ppm", ascending=False, ignore_index=True)
        self.c13_df.index = self.c13_df.index + 1

        self.pureshift_df = self.pureshift_df.sort_values(
            "ppm", ascending=False, ignore_index=True
        )
        self.pureshift_df.index = self.pureshift_df.index + 1

        # define short views of the dataframes and tidy up column names
        # attempt to remove solvent peaks
        self.hsqc = self.hsqc_df[self.hsqc_df.Type == "Compound"][
            ["f2_ppm", "f1_ppm", "intensity"]
        ].copy()
        self.c13 = self.c13_df[self.c13_df.Type == "Compound"][["ppm"]].copy()
        # define h1 datafram from h1_df instead of pureshift_df
        # self.h1 = self.pureshift_df[self.pureshift_df.Type == "Compound"][['ppm']].copy()
        self.h1 = self.h1_df[["ppm"]].copy()

        self.numCarbonGroups = self.c13.shape[0]
        self.numProtonGroups = self.h1.shape[0]

        self.symmetric_molecule = False
        if self.elements["C"] > 0 and self.elements["C"] > self.numCarbonGroups:
            self.symmetric_molecule = True

        self.c13["attached_protons"] = 0
        self.c13["ppmH1s"] = None

        self.hsqc["f1_i"] = 0
        self.hsqc["f2_i"] = 0
        self.hsqc["f2p_i"] = 0
        self.hsqc["f1C_i"] = 0
        self.hsqc["f2H_i"] = 0
        self.hsqc["f2Cp_i"] = 0
        self.hsqc["f2p_ppm"] = 0

        self.hmbc = self.hmbc_df[["f1_ppm", "f2_ppm", "intensity"]].copy()
        self.hmbc["f2p_ppm"] = 0
        self.hmbc["f1_i"] = 0
        self.hmbc["f2_i"] = 0
        self.hmbc["f2p_i"] = 0

        self.pureshift = self.pureshift_df.copy()
        self.cosy = self.cosy_df[["f1_ppm", "f2_ppm", "intensity"]].copy()
        self.cosy = self.cosy.assign(f1_i=lambda x: 0)
        self.cosy = self.cosy.assign(f1p_i=lambda x: 0)
        self.cosy = self.cosy.assign(f2_i=lambda x: 0)
        self.cosy = self.cosy.assign(f2p_i=lambda x: 0)

        self.cosy = self.cosy.assign(f1p_ppm=lambda x: np.nan)
        self.cosy = self.cosy.assign(f2p_ppm=lambda x: np.nan)

        # f1H_i	f2H_i	f1Cp_i	f2Cp_i
        self.cosy["f1H_i"] = ""
        self.cosy["f2H_i"] = ""
        self.cosy["f1Cp_i"] = ""
        self.cosy["f2Cp_i"] = ""

        self.c13["max_bonds"] = 4
        self.h1["integral"] = self.h1_df["numProtons"]
        self.h1["jCouplingClass"] = self.h1_df["jCouplingClass"]
        self.h1["jCouplingVals"] = self.h1_df["jCouplingVals"]
        self.h1["range"] = self.h1_df["range"]

        # tidy up chemical shift values by replacing cosy, hsqc and hmbc picked peaks with values from c13ppm and h1ppm dataframes

        # HMBC
        for i in self.hmbc.index:
            self.hmbc.loc[i, "f1_ppm"] = self.find_nearest(
                self.c13.ppm.tolist(), self.hmbc.loc[i, "f1_ppm"]
            )
            self.hmbc.loc[i, "f2_ppm"] = self.find_nearest(
                self.h1.ppm.tolist(), self.hmbc.loc[i, "f2_ppm"]
            )

        # HSQC
        for i in self.hsqc.index:
            self.hsqc.loc[i, "f1_ppm"] = self.find_nearest(
                self.c13.ppm.tolist(), self.hsqc.loc[i, "f1_ppm"]
            )
            self.hsqc.loc[i, "f2_ppm"] = self.find_nearest(
                self.h1.ppm.tolist(), self.hsqc.loc[i, "f2_ppm"]
            )

        # tidy up cosy H1 shifts
        for i in self.cosy.index:
            self.cosy.loc[i, "f1_ppm"] = self.find_nearest(
                self.h1.ppm.tolist(), self.cosy.loc[i, "f1_ppm"]
            )
            self.cosy.loc[i, "f2_ppm"] = self.find_nearest(
                self.h1.ppm.tolist(), self.cosy.loc[i, "f2_ppm"]
            )

        # add index columns to h1
        self.h1["label"] = ["H" + str(i) for i in self.h1.index]
        self.h1["f1H_i"] = ["H" + str(i) for i in self.h1.index]
        self.h1["f2H_i"] = ["H" + str(i) for i in self.h1.index]
        self.h1["f1_i"] = self.h1.index
        self.h1["f2_i"] = self.h1.index

        # add lookup dicts for dataframe h1
        self.H1ppmH1label = dict(zip(self.h1.ppm, self.h1.label))
        self.H1labelH1ppm = dict(zip(self.h1.label, self.h1.ppm))

        self.H1indexH1label = dict(zip(self.h1.index, self.h1.label))
        self.H1labelH1index = dict(zip(self.h1.label, self.h1.index))

        self.H1ppmH1index = dict(zip(self.h1.ppm, self.h1.index))
        self.H1indexH1ppm = dict(zip(self.h1.index, self.h1.ppm))

        # add index columns to c13
        self.c13["label"] = ["C" + str(i) for i in self.c13.index]
        self.c13["f2C_i"] = ["C" + str(i) for i in self.c13.index]
        self.c13["f1C_i"] = ["C" + str(i) for i in self.c13.index]
        self.c13["f1_i"] = self.c13.index
        self.c13["f2_i"] = self.c13.index

        # add lookup dicts for dataframe c13
        self.C13ppmC13label = dict(zip(self.c13.ppm, self.c13.label))
        self.C13labelC13ppm = dict(zip(self.c13.label, self.c13.ppm))

        self.C13indexC13label = dict(zip(self.c13.index, self.c13.label))
        self.C13labelC13index = dict(zip(self.c13.label, self.c13.index))

        self.C13ppmC13index = dict(zip(self.c13.ppm, self.c13.index))
        self.C13indexC13ppm = dict(zip(self.h1.index, self.c13.ppm))

        # open and read excel file into datframe and then process it
        # with open(self.yamlFiles[0], 'r') as fp:
        #    info = yaml.safe_load(fp)
        #    self.init_variables_from_dict(info)

        # add index columns to hsqc
        for i in self.hsqc.index:
            self.hsqc.loc[i, "f2_i"] = self.H1ppmH1index[self.hsqc.loc[i, "f2_ppm"]]
            self.hsqc.loc[i, "f2H_i"] = self.H1ppmH1label[self.hsqc.loc[i, "f2_ppm"]]
            self.hsqc.loc[i, "f1_i"] = self.C13ppmC13index[self.hsqc.loc[i, "f1_ppm"]]
            self.hsqc.loc[i, "f1C_i"] = self.C13ppmC13label[self.hsqc.loc[i, "f1_ppm"]]
            self.hsqc.loc[i, "f2Cp_i"] = self.C13ppmC13label[self.hsqc.loc[i, "f1_ppm"]]

        self.hsqc["f2p_i"] = self.hsqc["f1_i"]
        self.hsqc["f2p_ppm"] = self.hsqc["f1_ppm"]

        # add lookup dicts for hsqc
        self.hsqcH1ppmC13index = dict(zip(self.hsqc.f2_ppm, self.hsqc.f2p_i))
        self.hsqcH1ppmC13label = dict(zip(self.hsqc.f2_ppm, self.hsqc.f2Cp_i))
        self.hsqcH1ppmC13ppm = dict(zip(self.hsqc.f2_ppm, self.hsqc.f2p_ppm))

        self.hsqcH1indexC13index = dict(zip(self.hsqc.f2_i, self.hsqc.f2p_i))
        self.hsqcH1indexC13ppm = dict(zip(self.hsqc.f2_i, self.hsqc.f2p_ppm))
        self.hsqcH1indexC13label = dict(zip(self.hsqc.f2_i, self.hsqc.f2Cp_i))

        self.hsqcH1labelC13label = dict(zip(self.hsqc.f2H_i, self.hsqc.f1C_i))
        self.hsqcH1labelC13index = dict(zip(self.hsqc.f2H_i, self.hsqc.f1_i))
        self.hsqcH1labelC13ppm = dict(zip(self.hsqc.f2H_i, self.hsqc.f1_ppm))

        # fill in cosy dataframe
        for hppm in self.h1.ppm:
            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1_i"
            ] = self.H1ppmH1index.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2_i"
            ] = self.H1ppmH1index.get(hppm)

            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1p_i"
            ] = self.hsqcH1ppmC13index.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2p_i"
            ] = self.hsqcH1ppmC13index.get(hppm)

            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1p_ppm"
            ] = self.hsqcH1ppmC13ppm.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2p_ppm"
            ] = self.hsqcH1ppmC13ppm.get(hppm)

            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1H_i"
            ] = self.H1ppmH1label.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2H_i"
            ] = self.H1ppmH1label.get(hppm)

            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1Cp_i"
            ] = self.hsqcH1ppmC13label.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2Cp_i"
            ] = self.hsqcH1ppmC13label.get(hppm)

        # add index columns to hmbc
        self.hmbc["f1C_i"] = ""
        self.hmbc["f2H_i"] = ""
        self.hmbc["f2Cp_i"] = ""

        # fill in hmbc dataframe
        for i in self.hmbc.index:
            self.hmbc.loc[i, "f2p_ppm"] = self.hsqcH1ppmC13ppm.get(
                self.hmbc.loc[i, "f2_ppm"]
            )
            self.hmbc.loc[i, "f2p_i"] = self.hsqcH1ppmC13index.get(
                self.hmbc.loc[i, "f2_ppm"]
            )

            self.hmbc.loc[i, "f2_i"] = self.H1ppmH1index.get(self.hmbc.loc[i, "f2_ppm"])
            self.hmbc.loc[i, "f1_i"] = self.C13ppmC13index.get(
                self.hmbc.loc[i, "f1_ppm"]
            )

            self.hmbc.loc[i, "f1C_i"] = self.C13ppmC13label.get(
                self.hmbc.loc[i, "f1_ppm"]
            )
            self.hmbc.loc[i, "f2H_i"] = self.H1ppmH1label.get(
                self.hmbc.loc[i, "f2_ppm"]
            )
            self.hmbc.loc[i, "f2Cp_i"] = self.hsqcH1ppmC13label.get(
                self.hmbc.loc[i, "f2_ppm"]
            )

        self.create_new_nmrproblem_df()

        self.df.loc["ppm", self.protonAtoms] = self.h1.ppm.values
        self.df.loc["ppm", self.carbonAtoms] = self.c13.ppm.values

        self.df.loc[self.protonAtoms, "ppm"] = self.h1.ppm.values
        self.df.loc[self.carbonAtoms, "ppm"] = self.c13.ppm.values

        self.df.loc["integral", self.protonAtoms] = np.round(self.h1.integral.values)
        self.df.loc["J type", self.protonAtoms] = self.h1.jCouplingClass.values
        self.df.loc["J Hz", self.protonAtoms] = self.h1.jCouplingVals.values
        self.convertJHzToLists()

        self.df.loc["C13 hyb", self.protonAtoms + self.carbonAtoms] = 0
        self.df.loc["attached protons", self.protonAtoms + self.carbonAtoms] = 0

        self.df.loc["C13 hyb", self.protonAtoms] = np.round(self.h1.integral.values)
        self.df.loc["attached protons", self.protonAtoms] = np.round(
            self.h1.integral.values
        )

        self.df.loc["C13 hyb", self.carbonAtoms] = self.c13.attached_protons.values
        self.df.loc[
            "attached protons", self.carbonAtoms
        ] = self.c13.attached_protons.values

        self.updatecosygridfromExcel()

        # cosy = nmrproblem_excel.cosy
        # df = nmrproblem_excel.df
        for hi in self.cosy.f1H_i.unique():
            lll = self.cosy[self.cosy.f1H_i == hi]["f2H_i"].tolist()
            # print("lll", lll, hi)
            if hi in lll:
                lll.remove(hi)
            self.df.loc["cosy", hi] = lll

        for ci in self.cosy.f1Cp_i.unique():
            lll = self.cosy[self.cosy.f1Cp_i == ci]["f2Cp_i"].tolist()
            # remove None from list
            lll = [l for l in lll if l]
            if ci in lll:
                lll.remove(ci)
            # print("ci", ci, "lll", lll)
            # add ci only if not None
            if ci:
                self.df.loc["cosy", ci] = lll

        self.updateHSQCHMBCgridfromExcel()

        self.updatelistsfromgrid("hsqc", "o")
        self.updatelistsfromgrid("hmbc", "x")

        self.update_attachedprotons_c13hyb()

        return True

    def updatelistsfromgrid(self, exptname, marker="o"):
        hatoms = self.protonAtoms
        catoms = self.carbonAtoms
        hsqc = self.hsqc
        df = self.df

        H1toC13 = {}
        for hatom in hatoms:
            c_name = hsqc[hsqc["f2H_i"] == hatom]["f1C_i"].values
            if len(c_name == 1):
                H1toC13[hatom] = c_name[0]
            else:
                H1toC13[hatom] = hatom

        for atom in hatoms:
            df.loc[exptname, atom] = df.loc[catoms][
                df.loc[catoms, atom] == marker
            ].index.tolist()

        for atom in catoms:
            h1list = df.loc[hatoms][df.loc[hatoms, atom] == marker].index.tolist()
            # df.loc[exptname, atom] = [H1toC13[c] for c in h1list]
            df.loc[exptname, atom] = h1list

    def updateHSQCHMBCgridfromExcel(self):
        # hsqc = nmrproblem.hsqc
        # hmbc = nmrproblem.hmbc
        # df = nmrproblem.df
        # hatoms = nmrproblem.protonAtoms
        # catoms = nmrproblem.carbonAtoms

        for i in self.hsqc.index:
            f2H_i = self.hsqc.loc[i, "f2H_i"]
            f1C_i = self.hsqc.loc[i, "f1C_i"]

            self.df.loc[f1C_i, f2H_i] = "o"
            self.df.loc[f2H_i, f1C_i] = "o"

        for i in self.hmbc.index:
            f2H_i = self.hmbc.loc[i, "f2H_i"]
            f1C_i = self.hmbc.loc[i, "f1C_i"]

            self.df.loc[f1C_i, f2H_i] = "x"
            self.df.loc[f2H_i, f1C_i] = "x"

    def updatecosygridfromExcel(self):
        cosy = self.cosy
        df = self.df
        hatoms = self.protonAtoms
        catoms = self.carbonAtoms

        for i in cosy.index:
            f1H_i = cosy.loc[i, "f1H_i"]
            f2H_i = cosy.loc[i, "f2H_i"]
            if f1H_i != f2H_i:
                marker = "o"
            else:
                marker = "x"

            df.loc[f1H_i, f2H_i] = marker

    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def init_class_from_yml(self, ymlFileNameDirName: str):
        """
        read in class parameters from yaml found in problem directory if found

        Parameters
        ----------
        ymlFileNameDirName : str
            DESCRIPTION. name and path to directory holding yaml file

        Returns
        -------
        bool
            DESCRIPTION. Return True if yaml found and processed

        """

        # Look for yaml file in directory and use first one found
        if os.path.isdir(ymlFileNameDirName):
            self.yamlFiles = [
                os.path.join(ymlFileNameDirName, f)
                for f in os.listdir(ymlFileNameDirName)
                if f.endswith("yml")
            ]
            if len(self.yamlFiles) == 0:
                return False
        elif ymlFileNameDirName.endswith("yml") and os.path.exists(ymlFileNameDirName):
            self.yamlFiles = [ymlFileNameDirName]
        else:
            return False

        if not os.path.exists(self.yamlFiles[0]):
            return False

        # open and read yaml file into dictionary and then process it
        with open(self.yamlFiles[0], "r") as fp:
            info = yaml.safe_load(fp)
            self.init_variables_from_dict(info)

        return True

    def init_variables_from_dict(self, info: dict):
        """
        fill in class variables from dict produced by reading in yaml file.

        Parameters
        ----------
        info : dict
            DESCRIPTION. dictionary containing parameters that describe
            nmr problem.

        Returns
        -------
        None.

        """

        # columns and index  in best order for displaying df of parameters
        # describing problem
        self.dfColumns = info["dfColumns"]
        self.dfIndex = info["dfIndex"]

        # dataframe that has all the parameters that describe problem
        self.df = pd.DataFrame.from_dict(info["df"])
        if (len(self.dfColumns) > 0) and (len(self.dfIndex) > 0):
            self.df = self.df.loc[self.dfIndex, self.dfColumns]
        self.df_backup = self.df.copy()

        # parameters that describe the molecule in the problem and the number
        # of distinct peaks found in the carbon and proton NMR spectra
        self.moleculeAtomsStr = info.get("moleculeAtomsStr", self.moleculeAtomsStr)
        self.numProtonGroups = info.get("numProtonGroups", self.numProtonGroups)
        self.numCarbonGroups = info.get("numCarbonGroups", self.numCarbonGroups)

        # set of parameters that describe the spectrometer parameters
        # udic is based on nmrglue udic structure
        self.udic = info.get("udic", self.udic)

        # create nmr axes parameters from info in udic
        for i in range(self.udic["ndim"]):
            self.udic[i]["axis"] = ng.fileiobase.unit_conversion(
                self.udic[i]["size"],
                self.udic[i]["complex"],
                self.udic[i]["sw"],
                self.udic[i]["obs"],
                self.udic[i]["car"],
            )

        # create dbe and elements from moleculeAtomsStr
        self.calculate_dbe()

        # fill in protonAtoms and carbonAtoms list
        self.protonAtoms = ["H" + str(i + 1) for i in range(self.numProtonGroups)]
        self.carbonAtoms = ["C" + str(i + 1) for i in range(self.numCarbonGroups)]

        # put the list ofproton and carbon atoms in the udic for easier access
        self.udic[0]["atoms"] = self.protonAtoms
        self.udic[1]["atoms"] = self.carbonAtoms

    def calculate_dbe(self):
        """
        calculate DBE value for molecule and create a dictionary of elements
        and numbers found in molecule string

        Returns
        -------
        None.

        """
        # dbe_elements = ('C','H','N','F','Cl','Br')
        # match Element and number Cl, C3, O6, H
        aaa = re.findall(r"[A-Z][a-z]?\d?\d?\d?", self.moleculeAtomsStr)
        # match Element Cl, C, H, N
        eee = re.findall(r"[A-Z][a-z]?", self.moleculeAtomsStr)

        # create dictionary  of elements and number of elements
        self.elements = {}
        dbe_value = 0

        for e, a in zip(eee, aaa):
            if len(a) > len(e):
                num = a[len(e) :]
            else:
                num = "1"

            self.elements[e] = int(num)

        # calcluate DBE value formolecule
        if "C" in self.elements:
            dbe_value = self.elements["C"]
        if "N" in self.elements:
            dbe_value += self.elements["N"] / 2
        for e in ["H", "F", "Cl", "Br"]:
            if e in self.elements:
                dbe_value -= self.elements[e] / 2

        self.dbe = dbe_value + 1

    def dfToNumbers(self):
        """
        converts table contents from strings to floats and integers where appropriate

        Returns
        -------
        None.

        """
        self.df.loc["integral", self.protonAtoms + self.carbonAtoms] = self.df.loc[
            "integral", self.protonAtoms + self.carbonAtoms
        ].astype(int)
        self.df.loc["C13 hyb", self.protonAtoms + self.carbonAtoms] = self.df.loc[
            "C13 hyb", self.protonAtoms + self.carbonAtoms
        ].astype(int)
        self.df.loc["attached protons", self.carbonAtoms] = self.df.loc[
            "attached protons", self.carbonAtoms
        ].astype(int)
        self.df.loc["ppm", self.protonAtoms + self.carbonAtoms] = self.df.loc[
            "ppm", self.protonAtoms + self.carbonAtoms
        ].astype(float)
        self.df.loc[self.protonAtoms + self.carbonAtoms, "ppm"] = self.df.loc[
            self.protonAtoms + self.carbonAtoms, "ppm"
        ].astype(float)

    def updateDFtable(self, qdf: pd.DataFrame):
        """
        receives dataframe from widget and stores it in classs dataframe

        Parameters
        ----------
        qdf : pd.DataFrame
            DESCRIPTION. Dataframe that has been changed

        Returns
        -------
        bool
            DESCRIPTION. Returns True or False if problems in converting
            values in dataframes to the correct type ie floats, ints, etc

        """

        self.df_backup = self.df.copy()
        self.df = qdf.copy()

        df = self.df

        # copy ppm values in row to columns
        atoms = self.protonAtoms + self.carbonAtoms
        proton_atoms = self.protonAtoms
        carbon_atoms = self.carbonAtoms
        df.loc[atoms, "ppm"] = df.loc["ppm", atoms]

        # copy hmb and hsqc information over to carbon columns and proton rows
        for hi in proton_atoms:
            df.loc[hi, carbon_atoms] = df.loc[carbon_atoms, hi]

        hsqc = df.loc[proton_atoms, carbon_atoms]
        for hi in carbon_atoms:
            df.loc["hsqc", hi] = list(hsqc[hsqc[hi] == "o"].index)
            df.loc["hmbc", hi] = list(hsqc[hsqc[hi] == "x"].index)

        hsqc = df.loc[carbon_atoms, proton_atoms]
        for hi in proton_atoms:
            df.loc["hsqc", hi] = list(hsqc[hsqc[hi] == "o"].index)
            df.loc["hmbc", hi] = list(hsqc[hsqc[hi] == "x"].index)

        cosy = df.loc[proton_atoms]
        for hi in proton_atoms:
            df.loc["cosy", hi] = list(cosy[cosy[hi] == "o"].index)

        return True

        # # turn string values to ints, floats and lists
        # try:
        #     #self.dfToNumbers()
        #    # convertHSQCHMBCCOSYtoLists(self)
        #    # convertJHzToLists(self)

        #     # qdf = df

        #     return True
        # except:
        #     df = self.df_backup.copy()
        #     return False

    def createInfoDataframes(self):
        """
        Initialize info dataframes from main df dataframe

        Returns
        -------
        None.

        """

        self.udic[0]["info"] = self.df.loc[
            ["integral", "J type", "J Hz", "ppm", "cosy", "hsqc", "hmbc"],
            self.protonAtoms,
        ].T

        self.udic[0]["info"]["labels"] = self.udic[0]["info"].index

        self.udic[1]["info"] = self.df.loc[
            ["integral", "J type", "J Hz", "ppm", "cosy", "hsqc", "hmbc"],
            self.carbonAtoms,
        ].T

        self.udic[1]["info"]["labels"] = self.udic[1]["info"].index

    def update_molecule_ipywidgets(self, molecule_str: str, pGrps: int, cGrps: int):
        """
        updates the molecule and number of proton and xcarbon groups observed
        in NMR spectra.
        If any changed a new NMR problem dataframe is started from scratch

        Parameters
        ----------
        molecule_str : str
            DESCRIPTION. string containing the number and type of atoms
        pGrps : int
            DESCRIPTION. number of distinct proton signals in NMR spectrum
        cGrps : int
            DESCRIPTION. number of distinct carbon signals in NMR spectrum

        Returns
        -------
        None.

        """

        changed = False
        if self.moleculeAtomsStr != molecule_str:
            changed = True
        if self.numProtonGroups != pGrps:
            changed = True
        if self.numCarbonGroups != cGrps:
            changed = True
        self.moleculeAtomsStr = molecule_str
        self.numProtonGroups = pGrps
        self.numCarbonGroups = cGrps

        self.calculate_dbe()

        # print(self.elements)

        if changed:
            # delete old dataframe and then recreate it with new params
            self.create_new_nmrproblem_df()

    def update_molecule(self, molecule_str: str, pGrps: int, cGrps: int) -> bool:
        """
        updates the molecule and number of proton and xcarbon groups observed
        in NMR spectra.
        If any changed a new NMR problem dataframe is started from scratch

        Parameters
        ----------
        molecule_str : str
            DESCRIPTION. string containing the number and type of atoms
        pGrps : int
            DESCRIPTION. number of distinct proton signals in NMR spectrum
        cGrps : int
            DESCRIPTION. number of distinct carbon signals in NMR spectrum

        Returns
        -------
        bool
            DESCRIPTION. returns True if molecule or number of groups changed

        """

        changed = False
        if self.moleculeAtomsStr != molecule_str:
            changed = True
        if self.numProtonGroups != pGrps:
            changed = True
        if self.numCarbonGroups != cGrps:
            changed = True
        self.moleculeAtomsStr = molecule_str
        self.numProtonGroups = pGrps
        self.numCarbonGroups = cGrps

        self.calculate_dbe()

        # print(self.elements)

        return changed

    def create_new_nmrproblem_df(self):
        """
        Creates the dataframe template that holds the main info on the problem

        Returns
        -------
        None.

        """

        self.numProtonsInMolecule = self.elements["H"]
        self.numCarbonsInMolecule = self.elements["C"]

        self.protonAtoms = ["H" + str(i + 1) for i in range(self.numProtonGroups)]
        self.carbonAtoms = ["C" + str(i + 1) for i in range(self.numCarbonGroups)]

        self.dfIndex = (
            [
                "integral",
                "symmetry",
                "symmetry factor",
                "J type",
                "J Hz",
                "C13 hyb",
                "attached protons",
                "ppm",
            ]
            + self.protonAtoms[::-1]
            + self.carbonAtoms[::-1]
            + ["hsqc", "hmbc", "cosy"]
        )

        self.dfColumns = ["ppm"] + self.protonAtoms + self.carbonAtoms

        self.df = pd.DataFrame(index=self.dfIndex, columns=self.dfColumns)

        self.df = self.df.fillna("")

        # update df with default values
        self.df.loc["integral", self.protonAtoms + self.carbonAtoms] = [1,] * len(
            self.protonAtoms + self.carbonAtoms
        )
        self.df.loc["J type", self.protonAtoms + self.carbonAtoms] = ["s",] * len(
            self.protonAtoms + self.carbonAtoms
        )
        self.df.loc["J Hz", self.protonAtoms + self.carbonAtoms] = ["[0]",] * len(
            self.protonAtoms + self.carbonAtoms
        )
        self.df.loc["hsqc", self.protonAtoms + self.carbonAtoms] = ["[]",] * len(
            self.protonAtoms + self.carbonAtoms
        )
        self.df.loc["hmbc", self.protonAtoms + self.carbonAtoms] = ["[]",] * len(
            self.protonAtoms + self.carbonAtoms
        )
        self.df.loc["cosy", self.protonAtoms] = ["[]",] * len(self.protonAtoms)

        self.udic[0]["atoms"] = self.protonAtoms
        self.udic[1]["atoms"] = self.carbonAtoms

    def update_attachedprotons_c13hyb(self):

        for c in self.carbonAtoms:
            hsqc = self.df.loc["hsqc", c]
            attached_protons = 0
            for h in hsqc:
                attached_protons += int(self.df.loc["integral", h])

            self.df.loc["attached protons", c] = attached_protons
            self.df.loc["C13 hyb", c] = attached_protons

            for h in hsqc:
                self.df.loc["C13 hyb", h] = self.df.loc["C13 hyb", c]

    def extract_udic_base(self, udic):

        udicbasekeys = [
            "obs",
            "sw",
            "dw",
            "car",
            "size",
            "label",
            "complex",
            "encoding",
            "time",
            "lb",
        ]

        udicbase = {}
        udicbase["ndim"] = udic["ndim"]

        for i in range(udicbase["ndim"]):
            udicbase[i] = {}
            for k in udicbasekeys:
                udicbase[i][k] = udic[i][k]

        return udicbase

    def save_problem(self):

        tobesaved = {}

        if hasattr(self, "df"):
            if isinstance(self.df, pd.core.frame.DataFrame):
                tobesaved["df"] = self.df.to_dict()
                tobesaved["dfColumns"] = self.df.columns.tolist()
                tobesaved["dfIndex"] = self.df.index.tolist()

        if hasattr(self, "moleculeAtomsStr"):
            tobesaved["moleculeAtomsStr"] = self.moleculeAtomsStr

        if hasattr(self, "numProtonGroups"):
            tobesaved["numProtonGroups"] = self.numProtonGroups

        if hasattr(self, "numCarbonGroups"):
            tobesaved["numCarbonGroups"] = self.numCarbonGroups

        tobesaved["udic"] = self.extract_udic_base(self.udic)

        return tobesaved

    def create1H13Clabels(self, num_poss=3):

        udic = self.udic
        iprobs = self.iprobs
        df = self.df

        for i in range(udic["ndim"]):
            df = udic[i]["df"]
            udic[i]["poss_subst_list"] = []
            for k, n in enumerate(udic[i]["atoms"]):
                if num_poss == 1:
                    m = df.loc[iprobs[n][0], "sF_latex_matplotlib"]
                    sss = "{}: {}".format(n, m)
                    udic[i]["poss_subst_list"].append(sss)
                else:
                    output_list = [
                        df.loc[l, "sF_latex_matplotlib"] for l in iprobs[n][:2]
                    ]
                    sss = "{}$_{{{}}}$: ".format(n[0], n[1:])
                    sss += "\n".join(output_list)
                    udic[i]["poss_subst_list"].append(sss)

        for i in range(udic["ndim"]):

            h1_info_zipped = zip(
                udic[i]["info"].labels.apply(str).tolist(),
                udic[i]["info"].ppm.apply(str).tolist(),
                udic[i]["info"].integral.apply(str).tolist(),
            )

            labels1 = [
                "{}$_{{{}}}$ {} ppm\nIntegral: {}".format(s[0][0], s[0][1], s[1], s[2])
                for s in h1_info_zipped
            ]

            h1_info_zipped = zip(
                udic[i]["info"].labels.apply(str).tolist(),
                udic[i]["info"].ppm.apply(str).tolist(),
                udic[i]["info"]["J type"].apply(str).tolist(),
            )

            labels2 = [
                "{}$_{{{}}}$ {} ppm\nJ type: {}".format(s[0][0], s[0][1], s[1], s[2])
                for s in h1_info_zipped
            ]

            labels3 = udic[i]["poss_subst_list"]

            udic[i]["labels1_dict"] = {}
            for j, hi in enumerate(udic[i]["info"].index.tolist()):
                udic[i]["labels1_dict"][hi] = [labels3[j], labels2[j], labels1[j]]

    def calcProbDistFunctions(self, H1df_orig, C13df_orig):

        patoms = self.protonAtoms
        catoms = self.carbonAtoms

        # Figure out which dataframe has more rows
        if H1df_orig.index.size > C13df_orig.index.size:
            num_probs = H1df_orig.index.size
            iindex = H1df_orig.index
        else:
            num_probs = C13df_orig.index.size
            iindex = C13df_orig.index

        # create blank dataframe
        data = np.zeros((num_probs, len(patoms + catoms)))
        self.probDistFunctions = pd.DataFrame(
            data, index=iindex, columns=patoms + catoms
        )

        # Fill in probability density function table
        for H in patoms:
            for i in H1df_orig.index:
                ppm_val = float(self.df.loc[H, "ppm"])
                self.probDistFunctions.loc[i, H] = H1df_orig.loc[i, "norm"].pdf(ppm_val)

        for C in catoms:
            for i in C13df_orig.index:
                ppm_val = float(self.df.loc[C, "ppm"])
                self.probDistFunctions.loc[i, C] = C13df_orig.loc[i, "norm"].pdf(
                    ppm_val
                )

        self.H1df = pd.concat([H1df_orig, self.probDistFunctions[patoms]], axis=1)

        self.C13df = pd.concat([C13df_orig, self.probDistFunctions[catoms]], axis=1)

    def identify1HC13peaks(self):

        M_substituents = ["Li", "Na", "K", "Rb", "Cs", "Fr", "Si", "Al", "B"]

        elements = self.elements
        patoms = self.protonAtoms
        catoms = self.carbonAtoms
        df = self.df

        H1df = self.H1df
        C13df = self.C13df

        # reduce by DBE
        DBE = int(self.dbe)
        H1df = H1df[H1df["DBE"] <= DBE]
        C13df = C13df[C13df["DBE"] <= DBE]

        freeProtons = elements["H"] - df.loc["integral", patoms].sum()

        # If DBE equals 1, is it due to double bond oxygen, Nitrogen,
        # Halogen or alkene ?
        # carbonDBE = False
        oxygenDBE = False
        # halogenDBE = False
        # nitrogenDBE = False
        hydroxylsPresent = False
        numHydroxyls = 0
        freeOxygens = 0

        # Find out if there are hydroxyls
        freeProtons = elements["H"] - df.loc["integral", patoms].sum()
        if DBE == 0:
            if freeProtons > 0 and "O" in elements:
                hydroxylsPresent = True
                numHydroxyls = freeProtons
                freeOxygens = elements["O"] - freeProtons

        # Remove alkenes from table if DBE due to oxygen
        #
        # oxygenDBE = True
        # print(H1df.shape, C13df.shape)
        if ((DBE == 1) or (DBE >= 5 and elements["C"] > 6)) and (oxygenDBE):
            # remove alkene groups

            H1df = H1df[H1df.groupName != "alkene"]
            H1df = H1df[H1df.substituent != "alkene"]
            C13df = C13df[C13df.groupName != "alkene"]
            C13df = C13df[C13df.substituent != "alkene"]

        # remove choices where atoms not present, N, O, S, Halogens, metals
        # protons
        for e in ["O", "S", "N"]:
            if e not in elements:
                H1df = H1df[H1df[e] == 0]
            else:
                H1df = H1df[H1df[e] <= elements[e]]

        no_halides = True

        halide_elements = [e for e in elements.keys() if e in ["F", "Cl", "Br", "I"]]
        if len(halide_elements) > 0:
            no_halides = False

        if no_halides:
            H1df = H1df[H1df.substituent != "halide"]

        # remove metals from dataframe in none found
        no_Ms = True
        M_elements = [e for e in elements.keys() if e in M_substituents]
        if len(M_elements) > 0:
            no_Ms = False

        if no_Ms:
            H1df = H1df[H1df.substituent != "metal"]

        # carbons
        for e in ["O", "S", "N", "F", "Cl", "Br", "I"]:
            if e not in elements:
                C13df = C13df[C13df[e] == 0]
            else:
                C13df = C13df[C13df[e] <= elements[e]]

        # Now to reduce the choice further by using integrals, multiplicity.
        # This should now be done for individual peaks and the results kept separately
        # start with carbons
        self.iprobs = {}
        for c in catoms:
            kept_i = []
            for i in C13df.index:
                if int(df.loc["C13 hyb", c]) in C13df.loc[i, "attached_protons"]:
                    kept_i.append(i)
            self.iprobs[c] = kept_i

        # if attached protons / hybr == 3 remove other guesses
        # other than C13 methyl
        for c in catoms:
            if df.loc["C13 hyb", c] == 3:
                kept_i = []
                for i in self.iprobs[c]:
                    # if first pos in list is 3 then group in table is a methyl
                    if C13df.loc[i, "attached_protons"][0] == 3:
                        kept_i.append(i)
                self.iprobs[c] = kept_i

        # now reduce H1 options based on CH2, CH3, CH1
        # use the integral field in nmr properties table
        for h in patoms:
            kept_i = []
            for i in H1df.index:
                if df.loc["C13 hyb", h] in H1df.loc[i, "num_protons"]:
                    kept_i.append(i)
            self.iprobs[h] = kept_i

        # check to see if proton is attached to a carbon
        #
        for h in patoms:
            kept_i = []
            for i in self.iprobs[h]:
                if (int(df.loc["C13 hyb", h]) > 0) and (
                    H1df.loc[i, "carbon_attached"] == 1
                ):
                    kept_i.append(i)
                elif (int(df.loc["C13 hyb", h]) == 0) and (
                    H1df.loc[i, "carbon_attached"] == 0
                ):
                    kept_i.append(i)
            self.iprobs[h] = kept_i

        # sort the candidates, highes first
        for n, iii in self.iprobs.items():
            if n in patoms:
                sss = H1df.loc[iii, n].sort_values(ascending=False)
            elif n in catoms:
                sss = C13df.loc[iii, n].sort_values(ascending=False)
            self.iprobs[n] = sss.index.tolist()

    def convertJHzToLists(self):
        df = self.df
        patoms = self.protonAtoms
        catoms = self.carbonAtoms

        for i, j in enumerate(df.loc["J Hz", patoms]):
            if isinstance(j, str):
                df.loc["J Hz", patoms[i]] = [
                    float(k.strip()) for k in j.strip("][").split(",")
                ]
            elif isinstance(j, (float, int)):
                df.loc["J Hz", patoms[i]] = [
                    float(k.strip()) for k in str(j).strip("][").split(",")
                ]

        for i, j in enumerate(df.loc["J Hz", catoms]):
            if isinstance(j, str):
                df.loc["J Hz", catoms[i]] = [
                    float(k.strip()) for k in j.strip("][").split(",")
                ]
            elif isinstance(j, (float, int)):
                df.loc["J Hz", catoms[i]] = [
                    float(k.strip()) for k in str(j).strip("][").split(",")
                ]

    def calculate1H13CSpectra1D(self):

        couplings = {"s": 0, "d": 1, "t": 2, "q": 3, "Q": 4, "S": 6}

        udic = self.udic

        for e in [0, 1]:  # [proton, carbon]
            expt = udic[e]
            npts = expt["size"]
            lb = expt["lb"]

            fid = np.zeros(npts, dtype=np.complex128)
            dw = expt["dw"]
            ttt = np.linspace(0, dw * npts, npts)  # time array for calculating fid
            omega = expt["obs"]  # Larmor freq in MHz
            centre_freq = expt["car"]  # centre frequency in Hz

            expt["peak_ranges"] = {}

            for h1 in expt["atoms"]:
                iso_ppm = expt["info"].loc[h1, "ppm"]  # isotropic chemical shift in ppm
                integral = expt["info"].loc[h1, "integral"]
                jType = expt["info"].loc[
                    h1, "J type"
                ]  # coupling string "dd" doublet of doublets
                jHz = expt["info"].loc[h1, "J Hz"]  # list of J coupling values in Hz

                if not isinstance(jHz, Iterable):
                    jHz = [jHz]

                # calculate isotropic fid for indivudual resonances
                isofreq = (
                    iso_ppm * omega - centre_freq
                )  # isotropic chemical shift in Hz
                fid0 = integral * (
                    cos(-2.0 * pi * isofreq * ttt) + 1j * sin(-2.0 * pi * isofreq * ttt)
                )

                # add jcoupling modulation by iterating over string coupling values
                # for i, j in enumerate(jType):
                #     for k in range(couplings[j]):
                #         fid0 = fid0 * cos(pi * jHz[i] * ttt)

                for jc in jHz:
                    fid0 = fid0 * cos(pi * jc * ttt)
                fid += fid0

                # fft individual peaks to define peak limits
                # starting from the beginning and end of the spectrum
                # we move inwards if the value is not five times bigger
                # than the baseline value at the extremes
                fid0 = fid0 * exp(-pi * lb * ttt)
                spec = fftpack.fftshift(fftpack.fft(fid0)).real

                base_line_value = spec[-1]

                ileft = 0
                iright = -1
                while spec[ileft] < 5 * base_line_value:
                    ileft += 1
                while spec[iright] < 5 * base_line_value:
                    iright -= 1

                expt["peak_ranges"][h1] = [ileft, npts + iright]

            # fft complete fid and store it
            fid = fid * exp(-pi * lb * ttt)
            expt["spec"] = fftpack.fftshift(fftpack.fft(fid)).real

            expt["spec"] = 0.9 * (expt["spec"] / expt["spec"].max())

    def convertHSQCHMBCCOSYtoLists(self):
        #     global nmrproblem

        pcatoms = self.protonAtoms + self.carbonAtoms
        patoms = self.protonAtoms
        catoms = self.carbonAtoms
        df = self.df

        for i, j in enumerate(df.loc["hsqc", pcatoms]):
            if isinstance(j, str):
                df.loc["hsqc", pcatoms[i]] = [
                    float(k.strip()) for k in j.strip("][").split(",")
                ]

        for i, j in enumerate(df.loc["hmbc", pcatoms]):
            if isinstance(j, str):
                df.loc["hmbc", pcatoms[i]] = [
                    k.strip() for k in j.strip("][").split(",")
                ]

        for i, j in enumerate(df.loc["cosy", patoms]):
            if isinstance(j, str):
                df.loc["cosy", patoms[i]] = [
                    k.strip() for k in j.strip("][").split(",")
                ]

    def populateInfoDict(self):

        hinfo = self.udic[0]["info"]
        cinfo = self.udic[1]["info"]

        hinfo["peak_ranges_pts"] = hinfo["peak_ranges"].values()
        cinfo["peak_ranges_pts"] = cinfo["peak_ranges"].values()

        udic = self.udic

        # make sure ppm column values are floats
        for i in range(udic["ndim"]):
            udic[i]["info"].loc[udic[i]["atoms"], "ppm"] = udic[i]["info"].ppm.astype(
                float
            )

        # make sure integrals values are ints
        for i in range(udic["ndim"]):
            udic[i]["info"].loc[udic[i]["atoms"], "integral"] = udic[i][
                "info"
            ].integral.astype(int)

        # calculate peak ranges in ppm
        for i in range(udic["ndim"]):
            udic[i]["peak_ranges_ppm"] = {}

            for k, v in udic[i]["peak_ranges"].items():
                udic[i]["peak_ranges_ppm"][k] = [
                    udic[i]["axis"].ppm_scale()[udic[i]["peak_ranges"][k][0]],
                    udic[i]["axis"].ppm_scale()[udic[i]["peak_ranges"][k][1]],
                ]

            udic[i]["info"]["peak_ranges_ppm"] = udic[i]["peak_ranges_ppm"].values()

        # calculate peak positions in points and max height of peak
        for i in range(udic["ndim"]):
            udic[i]["info"]["pk_x"] = [
                -np.searchsorted(udic[i]["axis"].ppm_scale()[::-1], float(ppm))
                for ppm in udic[i]["info"].ppm
            ]

            udic[i]["info"]["pk_y"] = [
                udic[i]["spec"][p[0] : p[1]].max()
                for p in udic[i]["info"]["peak_ranges_pts"]
            ]

        # calculate peak widths in points starting from right hand side ie negative index
        for i in range(udic["ndim"]):
            for j in udic[i]["atoms"]:
                p = udic[i]["info"].loc[j, "peak_ranges_pts"]
                udic[i]["info"].loc[j, "pk_left"] = p[0] - udic[i]["size"]
                udic[i]["info"].loc[j, "pk_right"] = p[1] - udic[i]["size"]

    def save1DspecInfotoUdic(self):

        udic = self.udic

        udic[0]["info"]["peak_ranges_pts"] = udic[0]["peak_ranges"].values()
        udic[1]["info"]["peak_ranges_pts"] = udic[1]["peak_ranges"].values()

        for i in range(udic["ndim"]):
            udic[i]["peak_ranges_ppm"] = {}

            for k, v in udic[i]["peak_ranges"].items():
                udic[i]["peak_ranges_ppm"][k] = [
                    udic[i]["axis"].ppm_scale()[udic[i]["peak_ranges"][k][0]],
                    udic[i]["axis"].ppm_scale()[udic[i]["peak_ranges"][k][1]],
                ]

            udic[i]["info"]["peak_ranges_ppm"] = udic[i]["peak_ranges_ppm"].values()

        for i in range(udic["ndim"]):
            udic[i]["info"]["pk_x"] = [
                -np.searchsorted(udic[i]["axis"].ppm_scale()[::-1], float(ppm))
                for ppm in udic[i]["info"].ppm
            ]
            udic[i]["info"]["pk_y"] = [
                udic[i]["spec"][p[0] : p[1]].max()
                for p in udic[i]["info"]["peak_ranges_pts"]
            ]

        for i in range(udic["ndim"]):
            for j in udic[i]["atoms"]:
                p = udic[i]["info"].loc[j, "peak_ranges_pts"]
                udic[i]["info"].loc[j, "pk_left"] = p[0] - udic[i]["size"]
                udic[i]["info"].loc[j, "pk_right"] = p[1] - udic[i]["size"]


if __name__ == "__main__":

    h1 = r"csTables/h1_chemical_shift_table.jsn"
    c13 = r"csTables/c13_chemical_shift_table.jsn"

    H1df_orig, C13df_orig = read_in_cs_tables(h1, c13)

    nmrproblem = NMRproblem("exampleProblems\ch9_025", "yml")
