# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:55:03 2021.

@author: ERIC
"""

import os
import sys
from collections.abc import Iterable
import re
import json

import pandas as pd
# import openpyxl
import numpy as np
from numpy import pi, sin, cos, exp
from scipy import stats
from scipy import fftpack

import nmrglue as ng
import yaml

import networkx as nx

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

import PIL
from PIL import Image
from excelheaders import excel_orig_df_columns, excel_df_columns
# from .simpleNMR import JAVA_COMMAND
# from simpleNMR import JAVA_AVAILABLE

from PyQt5.QtWidgets import QMessageBox
import nmrmol

global JAVA_AVAILABLE
global JAVA_COMMAND

new_dataframes = {}


# return missing excelsheets names
def get_missing_sheets(excel_fname):

    # read in the excel file into a pandas dataframe

    excelsheets = pd.read_excel(
        excel_fname, sheet_name=None, index_col=0
    )
    # create a set of sheets keys that  are missing
    missing_sheets = set(excel_orig_df_columns.keys()) - set(excelsheets.keys())

    # determine if sheets present are empty
    empty_sheets = {key for key in excelsheets.keys() if excelsheets[key].empty}
    print("Empty sheets: ", empty_sheets)

    return missing_sheets.union(empty_sheets)



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

    h1_df = pd.read_json(h1)
    c13_df = pd.read_json(c13)

    # create mean and sigma based on min max chemical shifts
    h1_df["meanCS"] = (h1_df.minCS + h1_df.maxCS) / 2
    h1_df["sigmaCS"] = (h1_df.maxCS - h1_df.minCS) / scale_factor

    c13_df["meanCS"] = (c13_df.minCS + c13_df.maxCS) / 2
    c13_df["sigmaCS"] = (c13_df.maxCS - c13_df.minCS) / scale_factor

    # create probability density functions for each chemical shitf group
    for i in h1_df.index:
        h1_df.loc[i, "norm"] = stats.norm(
            loc=h1_df.loc[i, "meanCS"], scale=h1_df.loc[i, "sigmaCS"]
        )

    for i in c13_df.index:
        c13_df.loc[i, "norm"] = stats.norm(
            loc=c13_df.loc[i, "meanCS"], scale=c13_df.loc[i, "sigmaCS"]
        )

    return h1_df, c13_df


def parse_argv(my_argv=None,  qtstarted=True):
    """parse command line returnsa dictionary of arguments"""

    if my_argv is None:
        my_argv = sys.argv

    datadirectories = None
    smilefilenames = None
    pngfilenames = None
    excelfilenames = None
    smiles_fn = None
    png_fn = None
    excel_fn = None
    xy3_fn = None

    datadirectories = [d for d in my_argv[1:] if os.path.isdir(d)]
    excelfilenames = [e for e in my_argv[1:] if e.endswith(".xlsx")]
    pngfilenames = [s for s in my_argv[1:] if s.endswith(".png")]
    smilefilenames = [s for s in my_argv[1:] if (s.endswith(".smi") or s.endswith(".smiles"))]
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
                "xy3_fn": None
            }

    else:
        excelfilenames = [e for e in os.listdir(data_directory) if e.endswith(".xlsx")]
        if len(excelfilenames) > 0:
            excel_fn = os.path.join(data_directory, excelfilenames[0])
        else:
            # qt message box No excel file found in data directory
            if qtstarted:
                if len(my_argv) > 1:

                    msg_box = QMessageBox()
                    msg_box.setIcon(QMessageBox.Information)
                    msg_box.setText(f"No excel file found in data directory\n{data_directory}")
                    msg_box.setStandardButtons(QMessageBox.Ok)
                    msg_box.exec()
                    msg_box.setWindowTitle("Excel File Not Found")

            else:
                print(f"No excel file found in data directory\n{data_directory}")

            return {
                "data_directory": os.getcwd(),
                "excel_fn": None,
                "smiles_fn": None,
                "png_fn": None,
                "xy3_fn": None
            }

    if len(smilefilenames) > 0:
        if os.path.exists(os.path.join(data_directory, smilefilenames[0])):
            smiles_fn = os.path.join(data_directory, smilefilenames[0])
        else:
            smiles_fn = None
    else:
        smilefilenames = [s for s in os.listdir(data_directory) if (s.endswith(".smi") or s.endswith(".smiles"))]
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
        xy3jsonfilenames = [s for s in os.listdir(data_directory) if s=="xy3.json"]
        if len(xy3jsonfilenames) > 0:
            xy3_fn = os.path.join(data_directory, xy3jsonfilenames[0])
        else:
            xy3_fn = None

    return {
        "data_directory": data_directory,
        "excel_fn": excel_fn,
        "smiles_fn": smiles_fn,
        "png_fn": png_fn,
        "xy3_fn": xy3_fn
    }
     
def read_smiles_file(problemdata_info: dict)->str:
    """Reads the smiles file and returns the smiles string

    Args:
        problemdata_info (dict): parsed arguments from command line

    Returns:
        str: smiles string or None
    """
    smiles_str = None
    if isinstance(problemdata_info["smiles_fn"], str):
        with open(problemdata_info["smiles_fn"], "r", encoding="latin-1") as filepointer:
            smiles_str = filepointer.readline()

    return smiles_str


def read_png_file(problemdata_info: dict)->PIL.Image.Image:
    """Reads the png file and returns the image object"""
    png = None
    if isinstance(problemdata_info["png_fn"], str):
        png = Image.open(problemdata_info["png_fn"])

    return png


def create_png_from_smiles(smiles_str: str)->PIL.Image.Image:
    """Creates a png image from a smiles string via rdkit"""
    png = None
    molecule = Chem.MolFromSmiles(smiles_str)

    mol2 = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(mol2, randomSeed=3)
    molecule = Chem.RemoveHs(mol2)

    molecule.Compute2DCoords()
    
    # Draw.MolToFile(molecule, "molecule.png")
    # png = Image.open("molecule.png")
    png = Draw.MolToImage(molecule, size=(800,800))
    return png


def create_rdkit_molecule_from_smiles(smiles_str: str)->nmrmol.NMRmol:
    """Creates a RDKIT molecule from a smiles string via rdkit
       save the scaled coordinates of the atoms in the molecule
       min x and y =0, max x and y = 1"""

    
    # molecule = Chem.MolFromSmiles(smiles_str)
    print("create_rdkit_molecule_from_smiles")
    
    molecule = nmrmol.NMRmol.from_smiles(smiles_str)
    print(type(molecule))

    # attempt to fix coordinates of png
    mol2 = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(mol2, randomSeed=3)
    molecule = Chem.RemoveHs(mol2)

    molecule.Compute2DCoords()

    d2d = Draw.rdMolDraw2D.MolDraw2DSVG(800,800)
    d2d.DrawMolecule(molecule)
    d2d.FinishDrawing()

    for atom in molecule.GetAtoms():
        idx = atom.GetIdx()
        pt = d2d.GetDrawCoords(idx)
        print(idx, pt.x, pt.y)

        atom.SetDoubleProp('x', pt.x/800)
        atom.SetDoubleProp('y', pt.y/800)

    return nmrmol.NMRmol(molecule)


def return_carbon_xy3_positions(molecule: Chem.Mol)->dict:
    """Returns the xy3 positions of the carbon atoms in the molecule"""

    d2d = Draw.rdMolDraw2D.MolDraw2DSVG(800,800)
    d2d.DrawMolecule(molecule)
    d2d.FinishDrawing()
    xy3_positions = {}
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() == "C":
            idx = atom.GetIdx()
            point = d2d.GetDrawCoords(idx)
            xy3_positions['C'+str(idx+1)] = np.asarray([point.x/800., point.y/800.])
    return xy3_positions


def read_xy3_jsonfile(problemdata_info: dict)->dict:
    """Reads the xy3 json file and returns the json object"""
    xy3 = None
    if isinstance(problemdata_info["xy3_fn"], str):
        with open(problemdata_info["xy3_fn"], "r") as filepointer:
            xy3 = json.load(filepointer)

        for k, v in xy3.items():
            xy3[k] = np.asarray(v)

        return xy3
    else:
        return None


def create_hmbc_edges_dict(nmrproblem)->dict:
    """Creates a dictionary with the edges of the HMBC problem"""

    hmbc_edges = {}
    hmbc = nmrproblem.hmbc
    for i in hmbc.index:
        ci, cj = hmbc.loc[i, ["f1C_i", "f2Cp_i"]]
        if ci is None or cj is None:
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
    """Builds the model"""

    h1fn = r"csTables/h1_chemical_shift_table.jsn"
    c13fn = r"csTables/c13_chemical_shift_table.jsn"

    nmrproblem.H1df_orig, nmrproblem.C13df_orig = read_in_cs_tables(h1fn, c13fn)

    nmrproblem.calcProbDistFunctions(nmrproblem.H1df_orig, nmrproblem.C13df_orig)
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
    # hsqc = nmrproblem.hsqc
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
            if (not molecule.has_edge(c, l)) and (c != l):
                molecule.add_edge(c, l)

    # add node color to graph
    for i,n in enumerate(catoms):
        molecule.nodes[n]["node_color"] = nmrproblem.hmbc_edge_colors[c13.loc[i+1,"attached_protons"]]


def build_xy3_representation_of_molecule_from_smiles(nmrproblem):
    """Builds the xy3 representation of the molecule from the smiles string"""

    expected_molecule = nmrproblem.expected_molecule
    eigen_nodes = [a.GetSymbol() + str(a.GetIdx()) for a in expected_molecule.GetAtoms()]
    eigen_carbons = [s for s in eigen_nodes if 'C' in s]

    e_xy3 = {}
    xy3 = {}

    for n, (x,y,z) in zip(eigen_nodes, expected_molecule.GetConformer().GetPositions()):

        if 'C' in n:
            e_xy3[n]= np.array([x,y])

    for n1, n2 in zip(eigen_carbons, nmrproblem.carbonAtoms):
        xy3[n2] = e_xy3[n1]
    
    nmrproblem.xy3 = xy3





def build_xy3_representation_of_molecule(nmrproblem):
    """Builds the xy3 representation of the molecule"""
    molecule = nmrproblem.molecule
    # create coodinates that look more molecule like
    atomicNumberfromSymbol = {"H": 1, "Li": 3, "Be": 4, "B": 5, "C": 6,  "N": 7, "O": 8, "F": 9, "Na": 11, 
                              "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "K": 19, "Ca": 20,
                                "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
                                "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
                                "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
                                "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52,
                                "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
                                "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68,
                                "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76,
                                "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84,
                                "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92,
                                "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
                                "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107,
                                "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112, "Uut": 113, "Fl": 114,
                                "Uup": 115, "Lv": 116, "Uus": 117, "Uuo": 118}
               

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

    xy3 = return_carbon_xy3_positions(mmm)

    # print("xy3", xy3)

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


def create_hmbc_graph_fragments(nmrproblem, hmbc_edges: dict)-> dict:
    """create a graph of fragments from the hmbc_edges"""
    
    hmbc_graphs = {}
    # ntwk_labels = []

    catoms = nmrproblem.carbonAtoms
    # df = nmrproblem.df
    xy3 = nmrproblem.xy3
    # molecule = nmrproblem.molecule
    # udic = nmrproblem.udic

    # ret1 = None
    # ret2 = None

    lineCollections = []
    for i, c in enumerate(catoms):

        if c not in hmbc_edges.keys():
            continue

        # create hmbc graph for node c and add  xy coodinates
        hmbc_graphs[c] = {}
        hmbc_graphs[c]["graph"] = nx.Graph()
        hmbc_graphs[c]["xy"] = dict((k, xy3[k]) for k in [c] + list(hmbc_edges[c]) if k)
        hmbc_graphs[c]["colors"] = []

        # add nodes to hmbc graph
        hmbc_graphs[c]["graph"].add_nodes_from([c] + list(hmbc_edges[c]))

        # add edges
        for i, c1 in enumerate(hmbc_edges[c]):
            if c1:
                hmbc_graphs[c]["graph"].add_edge(c, c1)
                hmbc_graphs[c]["colors"].append(nmrproblem.hmbc_edge_colors[i])

    return hmbc_graphs


class NMRproblem:
    """NMR problem class"""

    def __init__(self, problemdata_info: dict, 
                       loadfromwhere=None, 
                       H1LarmorFreq=None, 
                       qtstarted=True, 
                       java_available=False, 
                       xy3_calc_method="xy3", 
                       java_command=None,
                       expts_available=[]):
        """initialize NMR problem class"""

        self.qtstarted = False
        self.c13_from_hsqc = False

        self.expts_available = expts_available
        self.java_available = java_available
        self.java_command = java_command
        self.xy3_calc_method = xy3_calc_method
        self.qtstarted = qtstarted
        self.H1LarmorFreq = H1LarmorFreq
        self.problemdata_info = problemdata_info
        self.loadfromwhere = loadfromwhere

        print("java available: ", self.java_available)
        print("xy3 calc method: ", self.xy3_calc_method)
        print("qtstarted: ", self.qtstarted)
        print("H1LarmorFreq: ", self.H1LarmorFreq)
        print("problemdata_info: ", self.problemdata_info)

        self.problemDirectoryPath = problemdata_info["data_directory"]
        self.problemDirectory = problemdata_info["data_directory"]
        self.rootDirectory, self.problemDirectory = os.path.split(
            problemdata_info["data_directory"]
        )

        self.hmbc_edge_colors = ('#1f77b4', '#ff7f0e', '#2ca02c',  '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                                  '#1f77b4', '#ff7f0e', '#2ca02c',  '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')

        self.xmin = -0.1
        self.ymin = -0.1
        self.xmax = 1.1
        self.ymax = 1.1
        
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

        self.molecule_defined = False
        self.expected_molecule = None  # rdkit molecule
        self.png = None                # image from png
        self.smiles = None
        self.smiles_defined = False
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
        self.xy3_from_json = False
        self.predict_c13ppm = True
        self.xy = None
        self.xy3 = None

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
                "car": (((16.0-(-2.0))/2)-2.0) * obs,
                "size": int(32 * 1024),
                "label": "1H",
                "complex": True,
                "encoding": "direct",
                "time": False,
                "freq": True,
                "lb": 0.5,
            },
            1: {
                "obs": obs / 4,
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
        if not isinstance(self.smiles, str):
            self.smiles = read_smiles_file(problemdata_info)
            if isinstance(self.smiles, str):
                self.smiles_defined = True
            else:
                self.smiles_defined = False

        if isinstance(self.smiles, str):
            self.smiles_defined = True
            self.png = create_png_from_smiles(self.smiles)
            self.expected_molecule = create_rdkit_molecule_from_smiles(self.smiles)



        # set xlims of 13C and 1H 1D spectra to +/- 10% of biggest and smallest ppm
        self.min_max_1D_ppm = []

        ppm_min = self.h1.ppm.min() - (self.h1.ppm.max() - self.h1.ppm.min()) * 0.3
        ppm_max = self.h1.ppm.max() + (self.h1.ppm.max() - self.h1.ppm.min()) * 0.3

        self.min_max_1D_ppm.append((ppm_max, ppm_min))

        ppm_min = self.c13.ppm.min() - (self.c13.ppm.max() - self.c13.ppm.min()) * 0.30
        ppm_max = self.c13.ppm.max() + (self.c13.ppm.max() - self.c13.ppm.min()) * 0.30

        self.min_max_1D_ppm.append((ppm_max, ppm_min))

        # check if c13 or h1 created from hsqc and hmbc
        # do this by checking for certain columns in c13 and h1 dataframes
        print("self.c13.columns", self.c13.columns)
        if 'CH3CH2' in self.c13.columns:
            self.c13_from_hsqc = True
        else:
            self.c13_from_hsqc = False

        print("*************************************")
        print("self.c13_from_hsqc", self.c13_from_hsqc)

        # print("self.expected_molecule\n", self.expected_molecule.molprops_df)

        if self.c13_from_hsqc:
            # see if any CH3 in expected molecule and then try to match them by the chemical shift
            # to the CH3CH2 column in c13
            # if there are matches then update the numProtons column to 3 where there is a match

            ch3_df = self.expected_molecule.molprops_df[self.expected_molecule.molprops_df.CH3]
            # drop duplicate rows based on c13ppm
            ch3_df = ch3_df.drop_duplicates(subset=['c13ppm'])
            
            # sort ch3_df by ppm, lowest to highest
            ch3_df = ch3_df.sort_values(by=['c13ppm'])
            ch3ch2_df = self.c13[self.c13.CH3CH2]

            print("\nch3_df\n", ch3_df)
            print("\nch3ch2_df\n", ch3ch2_df)

            for idx, ppm in zip(ch3_df.index, ch3_df.c13ppm):
                # find the closest match in the c13 dataframe
                # and update the numProtons column to 3
                closest_match = ch3ch2_df.iloc[(ch3ch2_df['ppm']-ppm).abs().argsort()[:1]]
                print("closest_match", closest_match)
                self.c13.loc[closest_match.index, 'numProtons'] = 3
                ch3ch2_df.drop(closest_match.index, inplace=True)

            # update corresponding integral in h1 frame
            ch3_df = self.c13[self.c13.numProtons == 3]
            for idx, ppm in zip(ch3_df.index, ch3_df.ppm):
                # find the closest match in the hsqc dataframe
                closest_match = self.hsqc.iloc[(self.hsqc['f1_ppm']-ppm).abs().argsort()[:1]]
                print("closest_match\n", closest_match) 
                # update the integral in the h1 dataframe
                self.h1.loc[closest_match.index, 'integral'] = 3

            print("self.h1\n", self.h1)

            # update h1 integrals for CH2 groups
            print("self.hsqc.columns", self.hsqc_df.columns)
            ch2_df = self.hsqc[self.hsqc_df.CH2]

            # count frequency of f1_ppm in ch2_df
            ch2_df_f1_ppm = ch2_df.f1_ppm.value_counts()
            print("ch2_df_f1_ppm\n", ch2_df_f1_ppm)
            # update corresponding integral in h1 frame
            for idx, f1_ppm, f2_ppm in zip(ch2_df.index, ch2_df.f1_ppm, ch2_df.f2_ppm):

                if ch2_df_f1_ppm.loc[f1_ppm] == 1:
                    # set the integral of f2_ppm in self.h1 to 2
                    self.h1.loc[self.h1.ppm == f2_ppm, 'integral'] = 2
                else:
                    # set the integral of f2_ppm in self.h1 to 1
                    self.h1.loc[self.h1.ppm == f2_ppm, 'integral'] = 1

            self.c13["attached_protons"] = self.c13["numProtons"]


    def tidy_up_excel_data(self):

        # define a consistant set of column names for dataframes
        # For now we are only working with excel sheets created from MestresNove
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

        # define short names for the dataframes
        if "molecule" in self.excelsheets:
            self.molecule_df = self.excelsheets["molecule"]
        else:
            print("No molecule sheet found")
            self.molecule_df = pd.DataFrame(columns=excel_orig_df_columns["molecule"])


        if "H1_1D" in self.excelsheets:  
            self.h1_df = self.excelsheets["H1_1D"]
        else:
            print("No H1_1D sheet found")
            self.h1_df = pd.DataFrame(columns=excel_orig_df_columns["H1_1D"])
            self.molecule_df
        if "COSY" in self.excelsheets:
            self.cosy_df = self.excelsheets["COSY"]
        else:
            print("No COSY sheet found")
            self.cosy_df = pd.DataFrame(columns=excel_orig_df_columns["COSY"])
        if "HSQC" in self.excelsheets:
            self.hsqc_df = self.excelsheets["HSQC"]
        else:
            print("No HSQC sheet found")
            self.hsqc_df = pd.DataFrame(columns=excel_orig_df_columns["HSQC"])
        if "HMBC" in self.excelsheets:
            self.hmbc_df = self.excelsheets["HMBC"]
        else:
            print("No HMBC sheet found")
            self.hmbc_df = pd.DataFrame(columns=excel_orig_df_columns["HMBC"])
        if "H1_pureshift" in self.excelsheets:
            self.pureshift_df = self.excelsheets["H1_pureshift"]
        else:
            print("No H1_pureshift sheet found")
            self.pureshift_df = pd.DataFrame(columns=excel_orig_df_columns["H1_pureshift"])
        if "C13_1D" in self.excelsheets:
            self.c13_df = self.excelsheets["C13_1D"]
        else:
            print("No C13_1D sheet found")
            self.c13_df = pd.DataFrame(columns=excel_orig_df_columns["C13_1D"])
        if "NOESY" in self.excelsheets:
            self.noesy_df = self.excelsheets["NOESY"]
        else:
            print("No NOESY sheet found")
            self.noesy_df = pd.DataFrame(columns=excel_orig_df_columns["NOESY"])

        # sanitize column names again just in case they were not in the original excel file
        for df in [self.h1_df, self.c13_df, self.pureshift_df, 
                    self.hsqc_df, self.hmbc_df, self.cosy_df,   self.molecule_df, self.noesy_df]:
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
            ) 

        


        # in case of nans values changed to 0 in h1_df.jCouplingVals
        # then set corresponding jCouplingClass to "s"
        self.h1_df.loc[self.h1_df.jCouplingVals == 0, "jCouplingClass"] = "s"

        # sort h1, c13 and pureshift jushapest in case they are out of order
        # reindex startting from 1

        self.h1_df = self.h1_df.sort_values('ppm', ascending=False, ignore_index=True)
        self.h1_df.index = self.h1_df.index + 1

        self.c13_df = self.c13_df.sort_values('ppm', ascending=False, ignore_index=True)
        self.c13_df.index = self.c13_df.index + 1

        self.pureshift_df = self.pureshift_df.sort_values('ppm', ascending=False, ignore_index=True)
        self.pureshift_df.index = self.pureshift_df.index + 1

    def init_h1_and_pureshift_from_c13_hsqc_hmbc_cosy(self):
        self.h1_df = pd.DataFrame(columns=['ppm', 'numProtons', 'integral', 'jCouplingVals', 'jCouplingClass', 'intensity', 'range', 'Type'])
        self.h1_df['ppm'] = self.hsqc_df['f2_ppm'].copy()
        self.h1_df['Type'] = self.hsqc_df['Type'].copy()

        # keep only rows where type is 'Compound'
        self.h1_df = self.h1_df[self.h1_df['Type'] == 'Compound']

        self.h1_df['numProtons'] = 1.0
        self.h1_df['integral'] = 1.0
        self.h1_df['jCouplingVals'] = 0.0
        self.h1_df['jCouplingClass'] = 's'
        self.h1_df['intensity'] = self.hsqc_df['intensity'].copy()
        self.h1_df['range'] = 0.0

        print("inside init_h1_and_pureshift_from_c13_hsqc_hmbc_cosy")
        print("self.h1_df\n", self.h1_df)

        # order h1_df by ppm in place highest to lowest, reset index to 1,2,3,4,5...
        self.h1_df = self.h1_df.sort_values('ppm', ascending=False, ignore_index=True)
        self.h1_df.index = self.h1_df.index + 1

        # create pureshift_df from h1_df
        self.pureshift_df = self.h1_df.copy()



    def find_and_group_CH2s(self, df1):

        # get a list of the CH2 resonances
        ch2_vals = df1[df1.CH2].f1_ppm.tolist()
        ch2_idx_vals = df1[df1.CH2].index.tolist()
        
        unique_ch2s = []
        unique_idxs = []

        #start with first hmbc value in the hmbc list
        # create a probability distribution around it and obtain the probability of all the other values to
        # to see if they are close to the first value.
        # all values that have a +ve probability are close to the first value
        # all values with a zero probability are not.
        # add the +ve to a saved list of lists "similar_hmbcs"
        # then remove them from the original list and repeat until original list length is zero

        while len(ch2_vals):
            # choose first from the list
            p0 = ch2_vals[0]
            # find list of hmbc values that are similar to the first in the list
            similar_ch2s = [p for p in ch2_vals if stats.norm.pdf(p, loc=p0, scale=0.01 ) > 0]
            similar_idxs = [i for i, p in zip(ch2_idx_vals, ch2_vals) if stats.norm.pdf(p, loc=p0, scale=0.01 ) > 0]
            # save the list
            unique_ch2s.append(similar_ch2s)
            unique_idxs.append(similar_idxs)

            # keep only hmbc values that were not similar 
            ch2_idx_vals = [i for i, p in zip(ch2_idx_vals, ch2_vals) if stats.norm.pdf(p, loc=p0, scale=0.01 ) == 0]
            ch2_vals = [p for p in ch2_vals if stats.norm.pdf(p, loc=p0, scale=0.01 ) == 0]
        

        # for i, h in zip(unique_idxs, unique_ch2s):
        #     print(np.mean(h), h, i)
            
        return unique_idxs, unique_ch2s



    def find_nearest(self, true_values:list, value:float):
        arraynp = np.asarray(true_values)
        idx = (np.abs(arraynp - value)).argmin()
        return arraynp[idx]


    def tidyup_ppm_values(self, df:pd.DataFrame, true_values:list, column_name: str):

        # make a copy of the column_name adding a suffix orig
        df[f"{column_name}_orig"] = df[column_name]
        
        # make a probability column to see how far replacement is from original
        df[f"{column_name}_prob"] = 0


        # create dataframe with ppm values and their nearest true value
        # dfnew = df.assign(nearest_true_value = df[column_name].apply(lambda x: find_nearest(true_values, x)))
        df[column_name] = df[column_name].apply(lambda x: self.find_nearest(true_values, x))
        
        #calculate probabilities
        for idx in df.index:
            df.loc[ idx, f"{column_name}_prob"] = stats.norm.pdf(df.loc[idx, column_name], loc=df.loc[idx, f"{column_name}_orig"], scale=0.005)

        return df


    def load_excel_file(self, excelFileNameDirName: str, excel_fn=None)->bool:
        # Load the excel file if provided

        print("self.qtstarted", self.qtstarted)

        if isinstance(excel_fn, str):
            self.excelFiles = [excel_fn]
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
                return False

        try:
            self.excelsheets = pd.read_excel(
                self.excelFiles[0], sheet_name=None, index_col=0
            )
            self.exceldf_attributes = {k: {'empty':df.empty, 'num_rows':df.shape[0]} 
                                            for k, df in self.excelsheets.items()}
            self.exceldf_attributes_df = pd.DataFrame(self.exceldf_attributes)
            
        except FileNotFoundError:
            # display qt message box  if excel file not found or not readable

            if self.qtstarted:
                msgBox = QMessageBox("Excel file not found", "Excel file not found", QMessageBox.Ok) 
                msgBox = QMessageBox()
                msgBox.setIcon(QMessageBox.Information)
                msgBox.setText("Excel File not Found\n{}".format(self.excelFiles[0]))
                msgBox.setWindowTitle("Excel File not Found")
                msgBox.setStandardButtons(QMessageBox.Ok)
                rtn = msgBox.exec_()
            else:
                print("Excel File not Found\n{}".format(self.excelFiles[0]))  
            return False
        except PermissionError:
            # display qt message box  if excel file not found or not readable

            if self.qtstarted:
                msgBox = QMessageBox()
                msgBox.setIcon(QMessageBox.Information)
                msgBox.setText("Cannot Open\n{}".format(self.excelFiles[0]))
                msgBox.setWindowTitle("Excel File Access Error")
                msgBox.setStandardButtons(QMessageBox.Ok)
                rtn = msgBox.exec_()
            else:
                print("Cannot Open\n{}".format(self.excelFiles[0]))
            return False
            
        return True

    def init_c13_from_hsqc_and_hmbc(self):

        # find CH2 groups in hsqc_df
        # CH2 groups where intensity is < 0.0

        self.hsqc_df['CH2'] = False
        self.hsqc_df.loc[self.hsqc_df.intensity < 0,"CH2"] = True

        unique_idxs, unique_ch2s = self.find_and_group_CH2s(self.hsqc_df)

        # replace hsqc CH2 values in f1_ppm of HSQC
        for idx, ch2 in zip(unique_idxs, unique_ch2s):
            self.hsqc_df.loc[idx, 'f1_ppm'] = np.mean(ch2)

        # replace values for f2_ppm in hmbc starting from f2_ppm hsqc
        self.hmbc_df = self.tidyup_ppm_values(self.hmbc_df, 
                            sorted(self.hsqc_df.f2_ppm.unique().tolist(), reverse=True), 'f2_ppm')

        # find all f1_ppm HMBC idx resonances that are not showing up in the HSQC f2_ppm
        iii = []
        for i in self.hmbc_df.index:
            prob_vals = []
            for c in self.hsqc_df.f1_ppm.unique():
                prob_vals.append(stats.norm.pdf(self.hmbc_df.loc[i,'f1_ppm'], loc=c, scale=0.01 ))
            if np.array(prob_vals).sum() == 0:  
                print(i, self.hmbc_df.loc[i,'f1_ppm'])
                iii.append(i)
            else:
                pass

        # keep only the unique hmbc resonances not in f1_ppm HSQC
        # get a list of the hmbc resonances
        hmbcs = self.hmbc_df.loc[iii, 'f1_ppm'].tolist()
        unique_hmbc = []

        #start with first hmbc value in the hmbc list
        # create a probability distribution around it and obtain the probability of all the other values to
        # to see if they are close to the first value.
        # all values that have a +ve probability are close to the first value
        # all values with a zero probability are not.
        # add the +ve to a saved list of lists "similar_hmbcs"
        # then remove them from the original list and repeat until original list length is zero

        while len(hmbcs):
            # choose first from the list
            p0 = hmbcs[0]
            # find list of hmbc values that are similar to the first in the list
            similar_hmbcs = [p for p in hmbcs if stats.norm.pdf(p, loc=p0, scale=0.01 ) > 0]
            # save the list
            unique_hmbc.append(similar_hmbcs)
            
            # keep only hmbc values that were not similar 
            hmbcs = [p for p in hmbcs if stats.norm.pdf(p, loc=p0, scale=0.01 ) == 0]

        for h in unique_hmbc:
            print(np.mean(h), h)

        # create an array of mean values for hmbc f1_ppm not found in hsqc f1_ppm
        mean_unique_hmbc_vals = [np.mean(h) for h in unique_hmbc]

        # tidyup f1_ppm values in hmbc that are not in f1_ppm HSQC
        hmbc_1  = self.tidyup_ppm_values(self.hmbc_df.loc[iii], mean_unique_hmbc_vals, 'f1_ppm')

        # tidyup f1_ppm values in hmbc that are in f1_ppm HSQC
        hmbc_2  = self.tidyup_ppm_values(self.hmbc_df.drop(iii), self.hsqc_df.f1_ppm.unique(), 'f1_ppm')

        # rejoin two parts of HMBC data
        self.hmbc_df = pd.concat([hmbc_1, hmbc_2])
        self.hmbc_df.sort_index(inplace=True)

        # add f2p_ppm column to HSQC and HMBC tables
        # f2p_ppm is C13 one to one relation between f1_ppm and f2_ppm in HSQC
        self.hmbc_df['f2p_ppm'] = 0.0
        for idx, f2ppmHSQC, f1ppmHSQC in zip(self.hsqc_df.index, self.hsqc_df.f2_ppm, self.hsqc_df.f1_ppm):
            self.hmbc_df.loc[self.hmbc_df.f2_ppm == f2ppmHSQC, 'f2p_ppm'] = f1ppmHSQC

        # return list of C13 values
        c13_list  = sorted(set(self.hmbc_df.f1_ppm).union(set(self.hmbc_df.f2p_ppm), set(self.hsqc_df.f1_ppm)), reverse=True)

        # create a dataframe of C13 values
        self.c13_df = pd.DataFrame(c13_list, columns=['ppm'])
        self.c13_df['Type'] = 'Compound'
        self.c13_df.index = self.c13_df.index+1
        
        # label c13 values as CH2 or not using the hsqc_df
        self.c13_df['CH2'] = False
        for idx, ppm in zip(self.c13_df.index, self.c13_df.ppm):
            if ppm in self.hsqc_df.f1_ppm.values:
                self.c13_df.loc[idx, 'CH2'] = self.hsqc_df[self.hsqc_df.f1_ppm == ppm]['CH2'].values[0]

        #  label c13 values as quartenary if not present in hsqc_df
        self.c13_df['Quartenary'] = False
        for idx, ppm in zip(self.c13_df.index, self.c13_df.ppm):
            if ppm not in self.hsqc_df.f1_ppm.values:
                self.c13_df.loc[idx, 'Quartenary'] = True

        # label c13 values as CH3CH2 if quaternary and CH2 both False
        self.c13_df['CH3CH2'] = False
        for idx, ppm in zip(self.c13_df.index, self.c13_df.ppm):
            if not self.c13_df.loc[idx, 'Quartenary'] and not self.c13_df.loc[idx, 'CH2']:
                self.c13_df.loc[idx, 'CH3CH2'] = True

        return self.c13_df

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


        if not self.load_excel_file(excelFileNameDirName, excel_fn):
            return False

        self.tidy_up_excel_data()

        if "molecule" in self.excelsheets:
            if "molecule" not in self.excelsheets["molecule"].columns:
                self.moleculeAtomsStr = ""
            else:
                self.moleculeAtomsStr = self.excelsheets["molecule"].molecule.values[0]
            self.molecule_defined = True
            print("self.moleculeAtomsStr", self.moleculeAtomsStr, type(self.moleculeAtomsStr))
            self.smiles = self.excelsheets["molecule"].smiles.values[0]
            if isinstance(self.smiles, str):
                self.smiles_defined = True
            else:
                self.smiles_defined = False
            print("self.smiles", self.smiles)
            if isinstance(self.moleculeAtomsStr, str) and len(self.moleculeAtomsStr) > 0:
                self.calculate_dbe()
            else: # define dbe, moleculeAtomStr and elements from smiles string
                self.expected_molecule = nmrmol.NMRmol.from_smiles(self.smiles)
                self.dbe = self.expected_molecule.dbe
                self.elements = self.expected_molecule.elements
                # reconstruct moleculeAtomsStr from elements
                self.moleculeAtomsStr = CalcMolFormula(self.expected_molecule)
                print("self.moleculeAtomsStr", self.moleculeAtomsStr)
        else:
            print("molecule not in sheet")
            self.dbe = 0
            self.elements = {'H': 0, 'C':0}
            self.molecule_defined = False


        print("self.expts_available", self.expts_available)
        # do basic checks on excels sheets and decide how to proceed
        if {"H1_pureshift", "HSQC"}.issubset(self.expts_available): 
            self.pureshift_data_present = True
        else:
            self.pureshift_data_present = False

        if {"C13_1D", "HSQC"}.issubset(self.expts_available):
            self.C13_data_present = True
        else:
            self.C13_data_present = False

        if {"H1_1D", "HSQC"}.issubset(self.expts_available):
            self.H1_data_present = True
        else:
            self.H1_data_present = False

        print("self.pureshift_data_present", self.pureshift_data_present)
        print("self.C13_data_present", self.C13_data_present)
        print("self.H1_data_present", self.H1_data_present)



        if (not self.pureshift_data_present and not self.H1_data_present and self.C13_data_present):

            print("##################################################")
            print("init h1 and pureshift from c13, hsqc, hmbc and cosy")

            # use f2_ppm values in hsqc_df to define h1_df:
            self.init_h1_and_pureshift_from_c13_hsqc_hmbc_cosy()

            print("self.h1_df", self.h1_df)
            self.c13_from_hsqc = True

        elif (not self.pureshift_data_present and not self.H1_data_present and not self.C13_data_present):
      
            self.c13_from_hsqc = True
            print("No c13_df, h1_df or pureshift_df")
            # use f2_ppm values in hsqc_df to define h1_df:
            self.c13_df = self.init_c13_from_hsqc_and_hmbc()
            self.init_h1_and_pureshift_from_c13_hsqc_hmbc_cosy()

        else:
            self.c13_from_hsqc = False
            print("*********************")
            print("h1_df not empty", self.h1_df.shape)
            print("pureshift_df not empty", self.pureshift_df.shape)

        # define short views of the dataframes and tidy up column names
        # attempt to remove solvent peaks
        print("self.hsqc_df.shape", self.hsqc_df.shape)
        print("self.hsqc_df", self.hsqc_df)
        # test if hsqc_df is empty
        if not self.hsqc_df.empty:
            self.hsqc = self.hsqc_df[self.hsqc_df.Type == "Compound"][['f2_ppm', 'f1_ppm', 'intensity']].copy()
        else:
            print("hsqc_df is empty")
            print("Program ending")
            return False
            # sys.exit()
            
        if self.c13_df.empty:
            return False
        else:
            self.c13 = self.c13_df[self.c13_df.Type == "Compound"][['ppm']].copy()
            self.c13['numProtons'] = 0
            self.c13['attached_protons'] = 0
            if {'CH2', 'CH3CH2', 'Quartenary'}.issubset(self.c13_df.columns):
                self.c13['CH2'] = self.c13_df[self.c13_df.Type == "Compound"]['CH2'].copy()
                self.c13['CH3CH2'] = self.c13_df[self.c13_df.Type == "Compound"]['CH3CH2'].copy()
                self.c13['Quartenary'] = self.c13_df[self.c13_df.Type == "Compound"]['Quartenary'].copy()

                # add two to numProtons for CH2
                self.c13.loc[self.c13.CH2, 'numProtons'] = 2
                # add one to numProtons for CH3CH2
                self.c13.loc[self.c13.CH3CH2, 'numProtons'] = 1

                print("self.expected_molecule\n", type(self.expected_molecule))

            self.c13.loc[:,'attached_protons'] = self.c13.numProtons


        if not self.pureshift_df.empty:
            self.h1 = self.pureshift_df[self.pureshift_df.Type == "Compound"][['ppm']].copy()
        elif not self.h1_df.empty:
            print("self.h1_df", self.h1_df.columns)
            self.h1 = self.h1_df[['ppm']].copy()
        else:
            print("No pureshift_df or h1_df")
            print("Program should not reach this point")
            return False

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

        # # replace any spaces in the column jCouplingVals with 0
        # print("self.h1_df", self.h1_df)
        # print("self.h1  before", self.h1)
        # for i in self.h1.index:
        #     if isinstance(self.h1.jCouplingVals[i], str):
        #         if self.h1.loc[i, "jCouplingVals"] in [" ", ""]:
        #             self.h1.loc[i, "jCouplingVals"] = 0
        #             self.h1.loc[i, "jCouplingClass"] = "s"
        # # self.h1["jCouplingVals"] = self.h1["jCouplingVals"].str.replace(" ", "0")
       


        # # set jCouplingClass to s if jCouplingVals is 0

        # print("self.h1 after", self.h1)


        # tidy up chemical shift values by replacing cosy, hsqc and hmbc picked peaks with values from c13ppm and h1ppm dataframes

        print(self.h1_df)

        print("self.h1.ppm.tolist()\n", self.h1.ppm.tolist())

        # HMBC
        self.hmbc = self.tidyup_ppm_values(self.hmbc, self.c13.ppm.tolist(), "f1_ppm")
        self.hmbc = self.tidyup_ppm_values(self.hmbc, self.h1.ppm.tolist(), "f2_ppm")
        # for i in self.hmbc.index:
        #     self.hmbc.loc[i, "f1_ppm"] = self.find_nearest(
        #         self.c13.ppm.tolist(), self.hmbc.loc[i, "f1_ppm"]
        #     )
        #     self.hmbc.loc[i, "f2_ppm"] = self.find_nearest(
        #         self.h1.ppm.tolist(), self.hmbc.loc[i, "f2_ppm"]
        #     )

        # HSQC
        self.hsqc = self.tidyup_ppm_values(self.hsqc, self.c13.ppm.tolist(), "f1_ppm")
        self.hsqc = self.tidyup_ppm_values(self.hsqc, self.h1.ppm.tolist(), "f2_ppm")
        
        # for i in self.hsqc.index:
        #     self.hsqc.loc[i, "f1_ppm"] = self.find_nearest(
        #         self.c13.ppm.tolist(), self.hsqc.loc[i, "f1_ppm"]
        #     )
        #     self.hsqc.loc[i, "f2_ppm"] = self.find_nearest(
        #         self.h1.ppm.tolist(), self.hsqc.loc[i, "f2_ppm"]
        #     )

        # tidy up cosy H1 shifts
        self.cosy = self.tidyup_ppm_values(self.cosy, self.h1.ppm.tolist(), "f1_ppm")
        self.cosy = self.tidyup_ppm_values(self.cosy, self.h1.ppm.tolist(), "f2_ppm")

        #check if any probability equals zero and remove the row
        # because it is likely that proton is not connected directly to carbon
        self.cosy.drop(self.cosy[self.cosy.f1_ppm_prob == 0].index, inplace=True)
        self.cosy.drop(self.cosy[self.cosy.f2_ppm_prob == 0].index, inplace=True)
        # for i in self.cosy.index:
        #     self.cosy.loc[i, "f1_ppm"] = self.find_nearest(
        #         self.h1.ppm.tolist(), self.cosy.loc[i, "f1_ppm"]
        #     )
        #     self.cosy.loc[i, "f2_ppm"] = self.find_nearest(
        #         self.h1.ppm.tolist(), self.cosy.loc[i, "f2_ppm"]
        #     )

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
            self.hmbc.loc[i, "f1_i"] = self.C13ppmC13index.get(self.hmbc.loc[i, "f1_ppm"])

            self.hmbc.loc[i, "f1C_i"] = self.C13ppmC13label.get(self.hmbc.loc[i, "f1_ppm"])
            self.hmbc.loc[i, "f2H_i"] = self.H1ppmH1label.get(self.hmbc.loc[i, "f2_ppm"])
            self.hmbc.loc[i, "f2Cp_i"] = self.hsqcH1ppmC13label.get(
                self.hmbc.loc[i, "f2_ppm"]
            )

        self.create_new_nmrproblem_df()

        print("self.protonAtoms", self.protonAtoms)
        print("self.carbonAtoms", self.carbonAtoms)

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

        print("self.c13.attached_protons.values", self.c13.attached_protons.values)
        print("self.h1.integral.values", self.h1.integral.values)
        self.df.loc["C13 hyb", self.carbonAtoms] = self.c13.attached_protons.values
        self.df.loc[
            "attached protons", self.carbonAtoms
        ] = self.c13.attached_protons.values

        print("self.df\n", self.df)

        self.updatecosygridfromExcel()

        # cosy = nmrproblem_excel.cosy
        # df = nmrproblem_excel.df
        for hi in self.cosy.f1H_i.unique():
            lll = self.cosy[self.cosy.f1H_i == hi]["f2H_i"].tolist()
            if hi in lll:
                lll.remove(hi)
            self.df.loc["cosy", hi] = lll

        for ci in self.cosy.f1Cp_i.unique():
            lll = self.cosy[self.cosy.f1Cp_i == ci]["f2Cp_i"].tolist()
            # remove None from list
            lll = [ l for l in lll if l]
            if ci in lll:
                lll.remove(ci)
            # add ci only if not None
            if ci:
                self.df.loc["cosy", ci] = lll

        self.updateHSQCHMBCgridfromExcel()

        self.updatelistsfromgrid("hsqc", "o")
        self.updatelistsfromgrid("hmbc", "x")

        self.update_attachedprotons_c13hyb()

        return True


    def build_xy3(self):

        #start off xy3 positions from random positions just in case no other method works

        self.init_xy3_from_random()

        if self.init_xy3_from_c13predictions():
            print("xy3 from c13 predictions")
            return True
        elif self.init_xy3_from_json():
            print("xy3 from json")
            return True
        else:
            print("xy3 from random positions")
            
            return self.init_xy3_from_random()

        print("xy3 not initiated")
        return False



    def init_xy3_from_random(self):
        """
        Initialise xy3 from random positions based on building a from cosy links
        """

        # if not self.xy3_calc_method == "random":
        #     return False

        build_xy3_representation_of_molecule(self)
        return True


    def init_xy3_from_json(self):
        if not self.xy3_calc_method == "xy3":
            return False

        self.xy3_from_json = True
        self.xy3 = read_xy3_jsonfile(self.problemdata_info)
        return isinstance(self.xy3, dict)


    def init_xy3_from_c13predictions(self):

        if  not isinstance(self.smiles, str):
            return False

        if not self.java_available:
            return False

        if not self.xy3_calc_method == "c13ppm":
            return False

        # just in case c13prediction fails initialise xy3 from random positions
        self.init_xy3_from_random()
        xy = np.array(list(self.xy3.values()))
        cxxx, cyyy = xy.T

        self.c13['x'] = cxxx
        self.c13['y'] = cyyy

        mol = create_rdkit_molecule_from_smiles(self.smiles)

        xy3 = return_carbon_xy3_positions(mol)
        # AllChem.Compute2DCoords(mol)

        xy = np.array(list(xy3.values()))
        cxxx, cyyy = xy.T

        molprops = [ [atom.GetIdx(), 
                      atom.GetNumImplicitHs(), 
                      atom.GetTotalNumHs(), 
                      atom.GetDegree(), 
                      atom.GetHybridization(), 
                      atom.GetIsAromatic()] for atom in mol.GetAtoms() if 'C' == atom.GetSymbol()]

        mol_df = pd.DataFrame(data=molprops, columns=['idx', 'implicitHs', 'totalNumHs', 'degree', 'hybridization', 'aromatic'])
        mol_df['x'] = cxxx
        mol_df['y'] = cyyy
        # mol_df['idx2'] = idx

        
        # mol_df['ppm'] = c13.sort_values(by=['C'], ignore_index=True)['predicted']

        with open("mol.mol", "w") as fp:
            fp.write(Chem.MolToMolBlock(mol))

        print("self.java_command", self.java_command)

        ret = os.system(self.java_command)

        if ret == 1:
            return False

        df = pd.read_csv("mol.csv")

        mol_df['C'] = df['C'].to_list()
        mol_df['predicted'] = df['mean'].to_list()

        print("mol_df", mol_df)

        mol_df = mol_df.sort_values(by=['predicted'], ignore_index=True, ascending=False)

        print("mol_df", mol_df)        

        self.c13['predicted'] = 0
        self.c13['C'] = 0
        # self.c13['x'] = 0
        # self.c13['y'] = 0

        print("mol_df", mol_df)

        # check if dataframes have the same number of rows
        # copy across mol information
        if self.c13.shape[0] == mol_df.shape[0]:
            self.c13['predicted'] = mol_df['predicted'].to_list()
            self.c13['C'] = mol_df['idx'].to_list()
            self.c13['x'] = mol_df['x'].to_list()
            self.c13['y'] = mol_df['y'].to_list()
        elif self.c13.shape[0] < df.shape[0]:
            # keep only unique rows based on mean ppm
            df_sym = mol_df.drop_duplicates(subset=['predicted'], ignore_index=True)
            print("df_sym", df_sym)
            print("=======================================")
            print(self.c13.shape[0], df_sym.shape[0])
            print("=======================================")
            if self.c13.shape[0] == df_sym.shape[0]:
                
                self.c13['predicted'] = df_sym['predicted'].to_list()
                self.c13['C'] = df_sym['idx'].to_list()
                print("idx", df_sym['idx'].to_list())
                self.c13['x'] = df_sym['x'].to_list()
                self.c13['y'] = df_sym['y'].to_list()

                mol_df = df_sym
        


        for numprotons in range(4):
            c13view = self.c13[self.c13.attached_protons==numprotons]
            mol_dfview = mol_df[mol_df.totalNumHs==numprotons].sort_values(by=['predicted'], ignore_index=True, ascending=False)

            print(numprotons, "c13view.shape[0]", c13view.shape[0], "mol_dfview.shape[0]", mol_dfview.shape[0])
            if c13view.shape[0]  == mol_dfview.shape[0]:
                for i1, i2 in zip(c13view.index, mol_dfview.index):
                    print(i1,i2)
                    self.c13.loc[i1, 'predicted'] = mol_dfview.loc[i2, 'predicted']
                    self.c13.loc[i1, 'C'] = mol_dfview.loc[i2, 'idx']
                    self.c13.loc[i1, 'x'] = mol_dfview.loc[i2, 'x']
                    self.c13.loc[i1, 'y'] = mol_dfview.loc[i2, 'y']
                    self.c13.loc[i1, 'attached_protons_predicted'] = mol_dfview.loc[i2,"totalNumHs"]
            else:
                print("shapes not equal")
       
        # initiate xy3

        print("self.c13", self.c13)
        self.xy3 = {}

        for i in self.c13.index:
            self.xy3[self.c13.loc[i,'label']] = [self.c13.loc[i,'x'], self.c13.loc[i,'y']]

        print(self.c13[['ppm', 'predicted', 'C', 'x', 'y']])

        return True





    # def init_class_from_excel(self, excelFileNameDirName: str, excel_fn=None):
    #     """
    #     read in class parameters from excel file found in problem directory if found

    #     Parameters
    #     ----------
    #     excelFileNameDirName : str
    #         DESCRIPTION. name and path to directory holding excel file

    #     Returns
    #     -------
    #     bool
    #         DESCRIPTION. Return True if excel found and processed

    #     """

    #     # Load the excel file if provided

    #     if isinstance(excel_fn, str):
    #         self.excelFiles = [excel_fn]
    #     else:
    #         # Look for excel file in directory and use first one found
    #         if os.path.isdir(excelFileNameDirName):
    #             self.excelFiles = [
    #                 os.path.join(excelFileNameDirName, f)
    #                 for f in os.listdir(excelFileNameDirName)
    #                 if f.endswith("xlsx")
    #             ]
    #             if len(self.excelFiles) == 0:
    #                 return False
    #         elif excelFileNameDirName.endswith("xlsx") and os.path.exists(
    #             excelFileNameDirName
    #         ):
    #             self.excelFiles = [excelFileNameDirName]
    #         else:
    #             return False

    #         if not os.path.exists(self.excelFiles[0]):
    #             return False

    #     try:
    #         self.excelsheets = pd.read_excel(
    #             self.excelFiles[0], sheet_name=None, index_col=0
    #         )
    #         self.exceldf_attributes = {k: {'empty':df.empty, 'num_rows':df.shape[0]} 
    #                                         for k, df in self.excelsheets.items()}
    #         self.exceldf_attributes_df = pd.DataFrame(self.exceldf_attributes)
            
    #     except FileNotFoundError:
    #         # display qt message box  if excel file not found or not readable

    #         if len(sys.argv == 1):
    #             print("Excel File not Found\n{}".format(self.excelFiles[0]))  
    #             return False
    #         if self.qtstarted:

    #             msgBox = QMessageBox("Excel file not found", "Excel file not found", QMessageBox.Ok) 
    #             msgBox = QMessageBox()
    #             msgBox.setIcon(QMessageBox.Information)
    #             msgBox.setText("Excel File not Found\n{}".format(self.excelFiles[0]))
    #             msgBox.setWindowTitle("Excel File not Found")
    #             msgBox.setStandardButtons(QMessageBox.Ok)
    #             rtn = msgBox.exec_()
    #         else:
    #             print("Excel File not Found\n{}".format(self.excelFiles[0]))  
    #         return False
    #     except PermissionError:
    #         # display qt message box  if excel file not found or not readable

    #         if self.qtstarted:
    #             msgBox = QMessageBox()
    #             msgBox.setIcon(QMessageBox.Information)
    #             msgBox.setText("Cannot Open\n{}".format(self.excelFiles[0]))
    #             msgBox.setWindowTitle("Excel File Access Error")
    #             msgBox.setStandardButtons(QMessageBox.Ok)
    #             rtn = msgBox.exec_()
    #         else:
    #             print("Cannot Open\n{}".format(self.excelFiles[0]))
    #         return False
            

    #     # define a consistant set of column names for dataframes
    #     # For now we are only working with excel sheets created from MestresNove
    #     for k, df in self.excelsheets.items():
    #         df.rename(
    #             columns={
    #                 "H's": "numProtons",
    #                 "Integral": "integral",
    #                 "J's": "jCouplingVals",
    #                 "Class": "jCouplingClass",
    #                 "Intensity": "intensity",
    #                 "Shift": "ppm",
    #                 "Range": "range",
    #                 "f2 (ppm)": "f2_ppm",
    #                 "f1 (ppm)": "f1_ppm",
    #             },
    #             inplace=True,
    #         )  # If they exist in dataframe!

    #     # replace any NaN values with 0
    #     print("self.excelsheets['H1_1D']", self.excelsheets['H1_1D'])
    #     for k, df in self.excelsheets.items():
    #         df.fillna(0, inplace=True)

        


    #     if "molecule" in self.excelsheets:
    #         self.molecule_defined = True
    #         self.moleculeAtomsStr = self.excelsheets["molecule"].molecule.values[0]
    #         print("self.moleculeAtomsStr", self.moleculeAtomsStr, len(self.moleculeAtomsStr))
    #         self.smiles = self.excelsheets["molecule"].smiles.values[0]
    #         if isinstance(self.smiles, str):
    #             self.smiles_defined = True
    #         else:
    #             self.smiles_defined = False
    #         print("self.smiles", self.smiles)
    #         if len(self.moleculeAtomsStr) > 0:
    #             self.calculate_dbe()
    #         else: # define dbe, moleculeAtomStr and elements from smiles string
    #             self.expected_molecule = nmrmol.NMRmol(self.smiles)
    #             self.dbe = self.expected_molecule.dbe
    #             self.elements = self.expected_molecule.elements
    #             # reconstruct moleculeAtomsStr from elements
    #             self.moleculeAtomsStr = ''
    #             for k, v in self.elements.items():
    #                 self.moleculeAtomsStr += k + str(v)
    #             print("self.moleculeAtomsStr", self.moleculeAtomsStr)
    #     else:
    #         print("molecule not in sheet")
    #         self.dbe = 0
    #         self.elements = {'H': 0, 'C':0}
    #         self.molecule_defined = False

    #     # define short names for the dataframes
    #     self.h1_df = self.excelsheets["H1_1D"]
    #     self.cosy_df = self.excelsheets["COSY"]
    #     self.hsqc_df = self.excelsheets["HSQC"]
    #     self.hmbc_df = self.excelsheets["HMBC"]
    #     self.pureshift_df = self.excelsheets["H1_pureshift"]
    #     self.c13_df = self.excelsheets["C13_1D"]

    #     # in case of nans values changed to 0 in h1_df.jCouplingVals
    #     # then set corresponding jCouplingClass to "s"
    #     self.h1_df.loc[self.h1_df.jCouplingVals == 0, "jCouplingClass"] = "s"

    #     # sort h1, c13 and pureshift jushapest in case they are out of order
    #     # reindex startting from 1

    #     self.h1_df = self.h1_df.sort_values('ppm', ascending=False, ignore_index=True)
    #     self.h1_df.index = self.h1_df.index + 1

    #     self.c13_df = self.c13_df.sort_values('ppm', ascending=False, ignore_index=True)
    #     self.c13_df.index = self.c13_df.index + 1

    #     self.pureshift_df = self.pureshift_df.sort_values('ppm', ascending=False, ignore_index=True)
    #     self.pureshift_df.index = self.pureshift_df.index + 1

    #     # define short views of the dataframes and tidy up column names
    #     # attempt to remove solvent peaks
    #     print("self.hsqc_df.shape", self.hsqc_df.shape)
    #     print("self.hsqc_df", self.hsqc_df)
    #     # test if hsqc_df is empty
    #     if not self.hsqc_df.empty:
    #         self.hsqc = self.hsqc_df[self.hsqc_df.Type == "Compound"][['f2_ppm', 'f1_ppm', 'intensity']].copy()
    #     else:
    #         print("hsqc_df is empty")
    #         print("Program ending")
    #         sys.exit()
            
    #     if not self.c13_df.empty:
    #         self.c13 = self.c13_df[self.c13_df.Type == "Compound"][['ppm']].copy()
    #     if not self.pureshift_df.empty:
    #         self.h1 = self.pureshift_df[self.pureshift_df.Type == "Compound"][['ppm']].copy()

    #     self.numCarbonGroups = self.c13.shape[0]
    #     self.numProtonGroups = self.h1.shape[0]

    #     self.symmetric_molecule = False
    #     if self.elements["C"] > 0 and self.elements["C"] > self.numCarbonGroups:
    #         self.symmetric_molecule = True


    #     self.c13["attached_protons"] = 0
    #     self.c13["ppmH1s"] = None

    #     self.hsqc["f1_i"] = 0
    #     self.hsqc["f2_i"] = 0
    #     self.hsqc["f2p_i"] = 0
    #     self.hsqc["f1C_i"] = 0
    #     self.hsqc["f2H_i"] = 0
    #     self.hsqc["f2Cp_i"] = 0
    #     self.hsqc["f2p_ppm"] = 0

    #     self.hmbc = self.hmbc_df[["f1_ppm", "f2_ppm", "intensity"]].copy()
    #     self.hmbc["f2p_ppm"] = 0
    #     self.hmbc["f1_i"] = 0
    #     self.hmbc["f2_i"] = 0
    #     self.hmbc["f2p_i"] = 0

    #     self.pureshift = self.pureshift_df.copy()
    #     self.cosy = self.cosy_df[["f1_ppm", "f2_ppm", "intensity"]].copy()
    #     self.cosy = self.cosy.assign(f1_i=lambda x: 0)
    #     self.cosy = self.cosy.assign(f1p_i=lambda x: 0)
    #     self.cosy = self.cosy.assign(f2_i=lambda x: 0)
    #     self.cosy = self.cosy.assign(f2p_i=lambda x: 0)

    #     self.cosy = self.cosy.assign(f1p_ppm=lambda x: np.nan)
    #     self.cosy = self.cosy.assign(f2p_ppm=lambda x: np.nan)

    #     # f1H_i	f2H_i	f1Cp_i	f2Cp_i
    #     self.cosy["f1H_i"] = ""
    #     self.cosy["f2H_i"] = ""
    #     self.cosy["f1Cp_i"] = ""
    #     self.cosy["f2Cp_i"] = ""

    #     self.c13["max_bonds"] = 4
    #     self.h1["integral"] = self.h1_df["numProtons"]
    #     self.h1["jCouplingClass"] = self.h1_df["jCouplingClass"]
    #     self.h1["jCouplingVals"] = self.h1_df["jCouplingVals"]
    #     self.h1["range"] = self.h1_df["range"]

    #     # # replace any spaces in the column jCouplingVals with 0
    #     # print("self.h1_df", self.h1_df)
    #     # print("self.h1  before", self.h1)
    #     # for i in self.h1.index:
    #     #     if isinstance(self.h1.jCouplingVals[i], str):
    #     #         if self.h1.loc[i, "jCouplingVals"] in [" ", ""]:
    #     #             self.h1.loc[i, "jCouplingVals"] = 0
    #     #             self.h1.loc[i, "jCouplingClass"] = "s"
    #     # # self.h1["jCouplingVals"] = self.h1["jCouplingVals"].str.replace(" ", "0")
       


    #     # # set jCouplingClass to s if jCouplingVals is 0

    #     # print("self.h1 after", self.h1)


    #     # tidy up chemical shift values by replacing cosy, hsqc and hmbc picked peaks with values from c13ppm and h1ppm dataframes

    #     # HMBC
    #     for i in self.hmbc.index:
    #         self.hmbc.loc[i, "f1_ppm"] = self.find_nearest(
    #             self.c13.ppm.tolist(), self.hmbc.loc[i, "f1_ppm"]
    #         )
    #         self.hmbc.loc[i, "f2_ppm"] = self.find_nearest(
    #             self.h1.ppm.tolist(), self.hmbc.loc[i, "f2_ppm"]
    #         )

    #     # HSQC
    #     for i in self.hsqc.index:
    #         self.hsqc.loc[i, "f1_ppm"] = self.find_nearest(
    #             self.c13.ppm.tolist(), self.hsqc.loc[i, "f1_ppm"]
    #         )
    #         self.hsqc.loc[i, "f2_ppm"] = self.find_nearest(
    #             self.h1.ppm.tolist(), self.hsqc.loc[i, "f2_ppm"]
    #         )

    #     # tidy up cosy H1 shifts
    #     for i in self.cosy.index:
    #         self.cosy.loc[i, "f1_ppm"] = self.find_nearest(
    #             self.h1.ppm.tolist(), self.cosy.loc[i, "f1_ppm"]
    #         )
    #         self.cosy.loc[i, "f2_ppm"] = self.find_nearest(
    #             self.h1.ppm.tolist(), self.cosy.loc[i, "f2_ppm"]
    #         )

    #     # add index columns to h1
    #     self.h1["label"] = ["H" + str(i) for i in self.h1.index]
    #     self.h1["f1H_i"] = ["H" + str(i) for i in self.h1.index]
    #     self.h1["f2H_i"] = ["H" + str(i) for i in self.h1.index]
    #     self.h1["f1_i"] = self.h1.index
    #     self.h1["f2_i"] = self.h1.index

    #     # add lookup dicts for dataframe h1
    #     self.H1ppmH1label = dict(zip(self.h1.ppm, self.h1.label))
    #     self.H1labelH1ppm = dict(zip(self.h1.label, self.h1.ppm))

    #     self.H1indexH1label = dict(zip(self.h1.index, self.h1.label))
    #     self.H1labelH1index = dict(zip(self.h1.label, self.h1.index))

    #     self.H1ppmH1index = dict(zip(self.h1.ppm, self.h1.index))
    #     self.H1indexH1ppm = dict(zip(self.h1.index, self.h1.ppm))

    #     # add index columns to c13
    #     self.c13["label"] = ["C" + str(i) for i in self.c13.index]
    #     self.c13["f2C_i"] = ["C" + str(i) for i in self.c13.index]
    #     self.c13["f1C_i"] = ["C" + str(i) for i in self.c13.index]
    #     self.c13["f1_i"] = self.c13.index
    #     self.c13["f2_i"] = self.c13.index

    #     # add lookup dicts for dataframe c13
    #     self.C13ppmC13label = dict(zip(self.c13.ppm, self.c13.label))
    #     self.C13labelC13ppm = dict(zip(self.c13.label, self.c13.ppm))

    #     self.C13indexC13label = dict(zip(self.c13.index, self.c13.label))
    #     self.C13labelC13index = dict(zip(self.c13.label, self.c13.index))

    #     self.C13ppmC13index = dict(zip(self.c13.ppm, self.c13.index))
    #     self.C13indexC13ppm = dict(zip(self.h1.index, self.c13.ppm))

    #     # open and read excel file into datframe and then process it
    #     # with open(self.yamlFiles[0], 'r') as fp:
    #     #    info = yaml.safe_load(fp)
    #     #    self.init_variables_from_dict(info)

    #     # add index columns to hsqc
    #     for i in self.hsqc.index:
    #         self.hsqc.loc[i, "f2_i"] = self.H1ppmH1index[self.hsqc.loc[i, "f2_ppm"]]
    #         self.hsqc.loc[i, "f2H_i"] = self.H1ppmH1label[self.hsqc.loc[i, "f2_ppm"]]
    #         self.hsqc.loc[i, "f1_i"] = self.C13ppmC13index[self.hsqc.loc[i, "f1_ppm"]]
    #         self.hsqc.loc[i, "f1C_i"] = self.C13ppmC13label[self.hsqc.loc[i, "f1_ppm"]]
    #         self.hsqc.loc[i, "f2Cp_i"] = self.C13ppmC13label[self.hsqc.loc[i, "f1_ppm"]]

    #     self.hsqc["f2p_i"] = self.hsqc["f1_i"]
    #     self.hsqc["f2p_ppm"] = self.hsqc["f1_ppm"]

    #     # add lookup dicts for hsqc
    #     self.hsqcH1ppmC13index = dict(zip(self.hsqc.f2_ppm, self.hsqc.f2p_i))
    #     self.hsqcH1ppmC13label = dict(zip(self.hsqc.f2_ppm, self.hsqc.f2Cp_i))
    #     self.hsqcH1ppmC13ppm = dict(zip(self.hsqc.f2_ppm, self.hsqc.f2p_ppm))

    #     self.hsqcH1indexC13index = dict(zip(self.hsqc.f2_i, self.hsqc.f2p_i))
    #     self.hsqcH1indexC13ppm = dict(zip(self.hsqc.f2_i, self.hsqc.f2p_ppm))
    #     self.hsqcH1indexC13label = dict(zip(self.hsqc.f2_i, self.hsqc.f2Cp_i))

    #     self.hsqcH1labelC13label = dict(zip(self.hsqc.f2H_i, self.hsqc.f1C_i))
    #     self.hsqcH1labelC13index = dict(zip(self.hsqc.f2H_i, self.hsqc.f1_i))
    #     self.hsqcH1labelC13ppm = dict(zip(self.hsqc.f2H_i, self.hsqc.f1_ppm))

    #     # fill in cosy dataframe
    #     for hppm in self.h1.ppm:
    #         self.cosy.loc[
    #             self.cosy[self.cosy.f1_ppm == hppm].index, "f1_i"
    #         ] = self.H1ppmH1index.get(hppm)
    #         self.cosy.loc[
    #             self.cosy[self.cosy.f2_ppm == hppm].index, "f2_i"
    #         ] = self.H1ppmH1index.get(hppm)

    #         self.cosy.loc[
    #             self.cosy[self.cosy.f1_ppm == hppm].index, "f1p_i"
    #         ] = self.hsqcH1ppmC13index.get(hppm)
    #         self.cosy.loc[
    #             self.cosy[self.cosy.f2_ppm == hppm].index, "f2p_i"
    #         ] = self.hsqcH1ppmC13index.get(hppm)

    #         self.cosy.loc[
    #             self.cosy[self.cosy.f1_ppm == hppm].index, "f1p_ppm"
    #         ] = self.hsqcH1ppmC13ppm.get(hppm)
    #         self.cosy.loc[
    #             self.cosy[self.cosy.f2_ppm == hppm].index, "f2p_ppm"
    #         ] = self.hsqcH1ppmC13ppm.get(hppm)

    #         self.cosy.loc[
    #             self.cosy[self.cosy.f1_ppm == hppm].index, "f1H_i"
    #         ] = self.H1ppmH1label.get(hppm)
    #         self.cosy.loc[
    #             self.cosy[self.cosy.f2_ppm == hppm].index, "f2H_i"
    #         ] = self.H1ppmH1label.get(hppm)

    #         self.cosy.loc[
    #             self.cosy[self.cosy.f1_ppm == hppm].index, "f1Cp_i"
    #         ] = self.hsqcH1ppmC13label.get(hppm)
    #         self.cosy.loc[
    #             self.cosy[self.cosy.f2_ppm == hppm].index, "f2Cp_i"
    #         ] = self.hsqcH1ppmC13label.get(hppm)

    #     # add index columns to hmbc
    #     self.hmbc["f1C_i"] = ""
    #     self.hmbc["f2H_i"] = ""
    #     self.hmbc["f2Cp_i"] = ""

    #     # fill in hmbc dataframe
    #     for i in self.hmbc.index:
    #         self.hmbc.loc[i, "f2p_ppm"] = self.hsqcH1ppmC13ppm.get(
    #             self.hmbc.loc[i, "f2_ppm"]
    #         )
    #         self.hmbc.loc[i, "f2p_i"] = self.hsqcH1ppmC13index.get(
    #             self.hmbc.loc[i, "f2_ppm"]
    #         )

    #         self.hmbc.loc[i, "f2_i"] = self.H1ppmH1index.get(self.hmbc.loc[i, "f2_ppm"])
    #         self.hmbc.loc[i, "f1_i"] = self.C13ppmC13index.get(self.hmbc.loc[i, "f1_ppm"])

    #         self.hmbc.loc[i, "f1C_i"] = self.C13ppmC13label.get(self.hmbc.loc[i, "f1_ppm"])
    #         self.hmbc.loc[i, "f2H_i"] = self.H1ppmH1label.get(self.hmbc.loc[i, "f2_ppm"])
    #         self.hmbc.loc[i, "f2Cp_i"] = self.hsqcH1ppmC13label.get(
    #             self.hmbc.loc[i, "f2_ppm"]
    #         )

    #     self.create_new_nmrproblem_df()

    #     self.df.loc["ppm", self.protonAtoms] = self.h1.ppm.values
    #     self.df.loc["ppm", self.carbonAtoms] = self.c13.ppm.values

    #     self.df.loc[self.protonAtoms, "ppm"] = self.h1.ppm.values
    #     self.df.loc[self.carbonAtoms, "ppm"] = self.c13.ppm.values

    #     self.df.loc["integral", self.protonAtoms] = np.round(self.h1.integral.values)
    #     self.df.loc["J type", self.protonAtoms] = self.h1.jCouplingClass.values
    #     self.df.loc["J Hz", self.protonAtoms] = self.h1.jCouplingVals.values
    #     self.convertJHzToLists()

    #     self.df.loc["C13 hyb", self.protonAtoms + self.carbonAtoms] = 0
    #     self.df.loc["attached protons", self.protonAtoms + self.carbonAtoms] = 0

    #     self.df.loc["C13 hyb", self.protonAtoms] = np.round(self.h1.integral.values)
    #     self.df.loc["attached protons", self.protonAtoms] = np.round(
    #         self.h1.integral.values
    #     )

    #     self.df.loc["C13 hyb", self.carbonAtoms] = self.c13.attached_protons.values
    #     self.df.loc[
    #         "attached protons", self.carbonAtoms
    #     ] = self.c13.attached_protons.values

    #     self.updatecosygridfromExcel()

    #     # cosy = nmrproblem_excel.cosy
    #     # df = nmrproblem_excel.df
    #     for hi in self.cosy.f1H_i.unique():
    #         lll = self.cosy[self.cosy.f1H_i == hi]["f2H_i"].tolist()
    #         if hi in lll:
    #             lll.remove(hi)
    #         self.df.loc["cosy", hi] = lll

    #     for ci in self.cosy.f1Cp_i.unique():
    #         lll = self.cosy[self.cosy.f1Cp_i == ci]["f2Cp_i"].tolist()
    #         # remove None from list
    #         lll = [ l for l in lll if l]
    #         if ci in lll:
    #             lll.remove(ci)
    #         # add ci only if not None
    #         if ci:
    #             self.df.loc["cosy", ci] = lll

    #     self.updateHSQCHMBCgridfromExcel()

    #     self.updatelistsfromgrid("hsqc", "o")
    #     self.updatelistsfromgrid("hmbc", "x")

    #     self.update_attachedprotons_c13hyb()

    #     return True

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

        # couplings = {"s": 0, "d": 1, "t": 2, "q": 3, "Q": 4, "S": 6}

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

                # calculate isotropic fid for individual resonances
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
