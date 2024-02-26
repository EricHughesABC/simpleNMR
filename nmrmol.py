# create nmrmol which is a subclass of rdkit.Chem.rdchem.Mol
import os
import platform
import numpy as np
import pandas as pd

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import PIL
import java

XYDIM = 800


class NMRmol(rdkit.Chem.rdchem.Mol):
    """NMRmol is a subclass of rdkit.Chem.rdchem.Mol"""

    def __init__(self, *args, **kwargs):

        """
        The default constructor.
        Note:
            This will be rarely used, as it can only create an empty molecule.
        Args:
            *args: Arguments to be passed to the rdkit Mol constructor.
            **kwargs: Arguments to be passed to the rdkit Mol constructor.
        """
        super().__init__(*args, **kwargs)

        # self.JAVA_AVAILABLE = java.JAVA_AVAILABLE
        # self.JAVA_COMMAND = java.JAVA_COMMAND

        # if platform.system() == "Darwin":
        #     self.MAC_OS = True
        #     self.WINDOWS_OS = False
        #     if not os.system(
        #         "jre/amazon-corretto-17.jdk/Contents/Home/bin/java --version"
        #     ):
        #         self.JAVA_AVAILABLE = True
        #         self.JAVA_COMMAND = "jre/amazon-corretto-17.jdk/Contents/Home/bin/java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"
        #         print("MAC Local JAVA is available")
        #     else:
        #         self.JAVA_AVAILABLE = False
        #         print("JAVA is not available")

        # elif platform.system() == "Windows":
        #     self.WINDOWS_OS = True
        #     self.MAC_OS = False
        #     # test if local windows ins installed
        #     if not os.system('"jre\\javawindows\\bin\\java -version"'):
        #         self.JAVA_AVAILABLE = True
        #         WINDOWS_OS = True
        #         self.JAVA_COMMAND = '"jre\\javawindows\\bin\\java -classpath predictorc.jar;cdk-2.7.1.jar;. NewTest mol.mol > mol.csv"'
        #         print("WINDOWS Local JAVA is available")
        #     else:
        #         self.JAVA_AVAILABLE = False
        #         print("JAVA is not available")

        self.dbe = self.calc_dbe()
        self.elements = self.init_elements_dict()
        self.num_carbon_atoms = self.elements.get("C", 0)
        self.num_hydrogen_atoms = self.elements.get("H", 0)
        self.png = self.create_png()
        self.smiles = self.create_smiles()
        self.has_symmetry = (
            len(self.GetSubstructMatches(self, uniquify=False, maxMatches=3)) > 1
        )

        # molprops = [
        #     [
        #         atom.GetIdx(),
        #         atom.GetNumImplicitHs(),
        #         atom.GetTotalNumHs(),
        #         atom.GetDegree(),
        #         atom.GetHybridization(),
        #         atom.GetIsAromatic(),
        #     ]
        #     for atom in self.GetAtoms()
        #     if "C" == atom.GetSymbol()
        # ]

        molprops = [
            [
                atom.GetIdx(),
                atom.GetNumImplicitHs(),
                atom.GetTotalNumHs(),
                atom.GetDegree(),
                atom.GetHybridization(),
                atom.GetIsAromatic(),
            ]
            for atom in self.GetAtoms()
            if atom.GetSymbol() == "C"
        ]

        self.molprops_df = pd.DataFrame(
            data=molprops,
            columns=[
                "idx",
                "implicitHs",
                "totalNumHs",
                "degree",
                "hybridization",
                "aromatic",
            ],
        )
        self.molprops_df = self.molprops_df.set_index(["idx"])

        print("\nstart calculating c13 chemical shifts\n")
        c13_chemical_shifts_df = self.calculated_c13_chemical_shifts()
        print("\nfinished calculating c13 chemical shifts\n")

        # add c13 chemical shifts to molprops_df if available
        if isinstance(c13_chemical_shifts_df, pd.DataFrame):
            self.molprops_df = self.molprops_df.join(c13_chemical_shifts_df, how="left")

        self.molprops_df["c13ppm"] = self.molprops_df["mean"]

        self.molprops_df["CH2"] = False
        # set CH2 to True if carbon has 2 protons attached
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 2, "CH2"] = True

        # define quaternary carbon atoms column for totalNumHs == 0
        self.molprops_df["quaternary"] = False
        self.molprops_df["CH0"] = False
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 0, "quaternary"] = True
        self.molprops_df["CH0"] = self.molprops_df["quaternary"]

        # define CH3 column for totalNumHs == 3
        self.molprops_df["CH3"] = False
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 3, "CH3"] = True

        # define CH1 column for totalNumHs == 1
        self.molprops_df["CH1"] = False
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 1, "CH1"] = True

        # define CH3CH1 column where CH3 or CH1 are True
        self.molprops_df["CH3CH1"] = False
        self.molprops_df["CH3CH1"] = self.molprops_df["CH3"] | self.molprops_df["CH1"]

        # calculate number of carbons with protons attached
        self.num_carbon_atoms_with_protons = self.molprops_df[
            self.molprops_df.totalNumHs > 0
        ].shape[0]

        # calculate number of carbons without protons attached
        self.num_quaternary_carbons = self.molprops_df[
            self.molprops_df.totalNumHs == 0
        ].shape[0]

        # calculate number of carbon with two protons attached
        self.num_ch2_carbon_atoms = self.molprops_df[
            self.molprops_df.totalNumHs == 2
        ].shape[0]

        # calculate number of carbon with three protons attached
        self.num_ch3_carbon_atoms = self.molprops_df[
            self.molprops_df.totalNumHs == 3
        ].shape[0]

        # calculate number of carbon with one proton  attached
        self.num_ch_carbon_atoms = self.molprops_df[
            self.molprops_df.totalNumHs == 1
        ].shape[0]

    # initialize class from smiles string
    @classmethod
    def from_smiles(cls, smiles: str):

        return cls(Chem.MolFromSmiles(smiles))

    def init_elements_dict(self):
        return (
            pd.DataFrame(
                [
                    [atom.GetIdx(), atom.GetSymbol()]
                    for atom in Chem.AddHs(self).GetAtoms()
                ],
                columns=["atom_index", "atom_symbol"],
            )["atom_symbol"]
            .value_counts()
            .to_dict()
        )

    # calculate DBE for molecule
    def calc_dbe(self) -> int:
        elements = self.init_elements_dict()
        if "C" in elements:
            dbe_value = elements["C"]
        if "N" in elements:
            dbe_value += elements["N"] / 2
        for e in ["H", "F", "Cl", "Br"]:
            if e in elements:
                dbe_value -= elements[e] / 2

        return dbe_value + 1

    def return_proton_groups(self) -> dict:
        # create pandas dataframe with carbon atom index and number of attached protons
        df = pd.DataFrame(
            [
                [atom.GetIdx(), atom.GetTotalNumHs()]
                for atom in self.GetAtoms()
                if atom.GetAtomicNum() == 6
            ],
            columns=["atom_index", "num_hydrogens"],
        )
        df = df["num_hydrogens"].value_counts().sort_index()
        return df.to_dict()

    # return dictionary of lists key is the number of protons attached to carbon, value is list of carbon atom indices
    def proton_groups(self) -> dict:
        proton_groups = {
            atom.GetTotalNumHs(): []
            for atom in self.GetAtoms()
            if atom.GetAtomicNum() == 6
        }
        for atom in self.GetAtoms():
            if atom.GetAtomicNum() == 6:
                proton_groups[atom.GetTotalNumHs()].append(atom.GetIdx())
        return proton_groups

    def calculated_c13_chemical_shifts(self) -> pd.DataFrame:
        return self.calc_c13_chemical_shifts_using_nmrshift2D()

    # return dictionary of dictionarys, first key is the number of protons, second key is carbon atom index, value is calculated C13 NMR chemical shift for molecule
    def c13_nmr_shifts(self) -> dict:
        c13_nmr_shifts = {
            k: {i: None for i in v} for k, v in self.proton_groups().items()
        }
        print("c13_nmr_shifts", c13_nmr_shifts)

        c13ppm_df = self.calc_c13_chemical_shifts_using_nmrshift2D()
        if isinstance(c13ppm_df, pd.DataFrame):
            # reset index to atom index
            print("c13ppm_df.index", c13ppm_df.index)
            print("c13ppm_df", c13ppm_df)
            for v in c13_nmr_shifts.values():
                for k2, v2 in v.items():
                    # find row in c13ppm_df with atom index k2
                    v[k2] = c13ppm_df.loc[k2, "mean"]

            print("c13_nmr_shifts", c13_nmr_shifts)
        return c13_nmr_shifts

    def calc_c13_chemical_shifts_using_nmrshift2D(self) -> pd.DataFrame:

        with open("mol.mol", "w") as fp:
            fp.write(Chem.MolToMolBlock(self))

        ret = os.system(self.JAVA_COMMAND)

        if ret == 1:
            print("NMRShift2D failed to calculate C13 chemical shifts")
            return False
        else:
            mol_df = pd.read_csv("mol.csv", index_col=0)
            mol_df.index = mol_df.index - 1
            return mol_df

    def num_protons_attached_to_carbons(self) -> int:
        return sum(k * v for k, v in self.return_proton_groups().items())

    def num_proton_groups(self) -> int:
        return sum(self.return_proton_groups().values())

    def aromatic_carbons(self) -> list:
        return [atom.GetIdx() for atom in self.GetAtoms() if atom.GetIsAromatic()]

    # return list idx of carbon atoms in molecule
    def carbon_idx(self) -> list:
        return [atom.GetIdx() for atom in self.GetAtoms() if atom.GetAtomicNum() == 6]

    def query_molprops_df(self, column_id: str, value) -> pd.DataFrame:
        return self.molprops_df[self.molprops_df[column_id] == value][
            ["implicitHs", "degree", "aromatic", "hybridization", "CH2", "c13ppm"]
        ]

    def create_smiles(self) -> str:
        """Creates a smiles string from a rdkit molecule"""
        return Chem.MolToSmiles(self)

    def create_png(self) -> PIL.Image.Image:
        """Creates a png image from a smiles string via rdkit"""
        png = None

        mol2 = Chem.AddHs(self)
        AllChem.EmbedMolecule(mol2, randomSeed=3)
        rdkit_molecule = Chem.RemoveHs(mol2)

        rdkit_molecule.Compute2DCoords()

        return Draw.MolToImage(rdkit_molecule, size=(XYDIM, XYDIM))


# class expectedMolecule:

#     if platform.system() == "Darwin":
#         MAC_OS = True
#         WINDOWS_OS = False
#         if not os.system(
#             "jre/amazon-corretto-17.jdk/Contents/Home/bin/java --version"
#         ):
#             JAVA_AVAILABLE = True
#             JAVA_COMMAND = "jre/amazon-corretto-17.jdk/Contents/Home/bin/java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"
#         else:
#             JAVA_AVAILABLE = False

#     elif platform.system() == "Windows":
#         WINDOWS_OS = True
#         MAC_OS = False
#         # test if local windows ins installed
#         if not os.system('"jre\\javawindows\\bin\\java -version"'):
#             JAVA_AVAILABLE = True
#             WINDOWS_OS = True
#             JAVA_COMMAND = '"jre\\javawindows\\bin\\java -classpath predictorc.jar;cdk-2.7.1.jar;. NewTest mol.mol > mol.csv"'
#         else:
#             JAVA_AVAILABLE = False

#     def __init__(self, smiles_str):
#         self.smiles_str = smiles_str
#         self.mol = Chem.MolFromSmiles(smiles_str)

#         self.dbe = self.calc_dbe()
#         self.elements = self.init_elements_dict()
#         self.num_carbon_atoms = self.elements.get("C", 0)
#         self.num_hydrogen_atoms = self.elements.get("H", 0)

#         molprops = [
#             [
#                 atom.GetIdx(),
#                 atom.GetNumImplicitHs(),
#                 atom.GetTotalNumHs(),
#                 atom.GetDegree(),
#                 atom.GetHybridization(),
#                 atom.GetIsAromatic(),
#             ]
#             for atom in self.mol.GetAtoms()
#             if atom.GetSymbol() == "C"
#         ]

#         self.molprops_df = pd.DataFrame(
#             data=molprops,
#             columns=[
#                 "idx",
#                 "implicitHs",
#                 "totalNumHs",
#                 "degree",
#                 "hybridization",
#                 "aromatic",
#             ],
#         )
#         self.molprops_df = self.molprops_df.set_index(["idx"])

#         # define quaternary carbon atoms column for totalNumHs == 0
#         self.molprops_df["quaternary"] = False
#         self.molprops_df["CH0"] = False
#         self.molprops_df.loc[self.molprops_df["totalNumHs"] == 0, "quaternary"] = True
#         self.molprops_df["CH0"] = self.molprops_df["quaternary"]

#         # define CH1 column for totalNumHs == 1
#         self.molprops_df["CH1"] = False
#         self.molprops_df.loc[self.molprops_df["totalNumHs"] == 1, "CH1"] = True

#         self.molprops_df["CH2"] = False
#         # set CH2 to True if carbon has 2 protons attached
#         self.molprops_df.loc[self.molprops_df["totalNumHs"] == 2, "CH2"] = True

#         # define CH3 column for totalNumHs == 3
#         self.molprops_df["CH3"] = False
#         self.molprops_df.loc[self.molprops_df["totalNumHs"] == 3, "CH3"] = True

#         # define CH3CH1 column where CH3 or CH1 are True
#         self.molprops_df["CH3CH1"] = False
#         self.molprops_df["CH3CH1"] = self.molprops_df["CH3"] | self.molprops_df["CH1"]

#         self.c13ppm = self.calculated_c13_chemical_shifts()
#         self.molprops_df["ppm"] = self.c13ppm["mean"]

#         # check if there are rings in the molecule
#         self.molprops_df["ring_idx"] = -1
#         self.molprops_df["ring_size"] = 0
#         ring_atoms = self.GetRingSystems()

#         # create sets of carbon atoms in rings store them in molprops_df column aromatic_rings
#         for ring_idx, ring  in enumerate(ring_atoms):
#             carbon_atoms_in_ring = [i for i in ring if self.mol.GetAtomWithIdx(i).GetSymbol() == "C"]
#             self.molprops_df.loc[carbon_atoms_in_ring, "ring_idx"] = ring_idx


#         # set the ring size for each carbon atom
#         for ring in ring_atoms:
#             for atom in ring:
#                 # check if atom is carbon
#                 if self.mol.GetAtomWithIdx(atom).GetSymbol() == "C":
#                     self.molprops_df.loc[atom, "ring_size"] = len(ring)

#         self.png = self.create_png()

#         self.has_symmetry = len(self.mol.GetSubstructMatches(self.mol, uniquify=False, maxMatches=3)) > 1
#         print("has symmetry", self.has_symmetry)

#         idx_list, xxx, yyy = self.calc_carbon_xy_positions_png(self.mol)
#         self.molprops_df["x"] = 0
#         self.molprops_df["y"] = 0
#         self.molprops_df.loc[idx_list, "x"] = xxx
#         self.molprops_df.loc[idx_list, "y"] = yyy


#         self.molprops_df["symmetry_idx1"] = -1
#         self.molprops_df["symmetry_idx2"] = -1

#         # find the symmetry atoms in each ring and non ring group and assign the symmetry indices to each
#         for ring_idx in self.molprops_df.ring_idx.unique():
#             ring_df = self.molprops_df[self.molprops_df.ring_idx == ring_idx]
#             for ppm in ring_df.ppm.unique():
#                 ring_df_ppm = ring_df[ring_df.ppm == ppm]
#                 if len(ring_df_ppm) == 2:
#                     self.molprops_df.loc[ring_df_ppm.index[0], "symmetry_idx1"] = ring_df_ppm.index[1]
#                     self.molprops_df.loc[ring_df_ppm.index[1], "symmetry_idx1"] = ring_df_ppm.index[0]
#                 elif len(ring_df_ppm) == 4 and ring_idx > -1:
#                     self.molprops_df.loc[ring_df_ppm.index[0], "symmetry_idx1"] = ring_df_ppm.index[3]
#                     self.molprops_df.loc[ring_df_ppm.index[1], "symmetry_idx1"] = ring_df_ppm.index[2]
#                     self.molprops_df.loc[ring_df_ppm.index[2], "symmetry_idx1"] = ring_df_ppm.index[1]
#                     self.molprops_df.loc[ring_df_ppm.index[3], "symmetry_idx1"] = ring_df_ppm.index[0]


#         # calculate number of carbons without protons attached
#         self.num_quaternary_carbons = self.molprops_df[self.molprops_df.quaternary].shape[0]

#         # calculate number of carbon with two protons attached
#         self.num_ch2_carbon_atoms = self.molprops_df[self.molprops_df.CH2].shape[0]

#         # calculate number of carbon with three protons attached
#         self.num_ch3_carbon_atoms = self.molprops_df[self.molprops_df.CH3].shape[0]

#         # calculate number of carbon with one proton  attached
#         self.num_ch_carbon_atoms = self.molprops_df[self.molprops_df.CH1].shape[0]

#         # calculate number of carbons with protons attached
#         self.num_carbon_atoms_with_protons = self.num_ch2_carbon_atoms  + self.num_ch3_carbon_atoms + self.num_ch_carbon_atoms

#         # check if there are aromatic rings that map symmetrically to each other
#         self.mapped_symmetric_aromatic_rings = self.map_symmetric_aromatic_rings(self.molprops_df)

#         # create reduced  molprops_df based on symmetry of NMR ppm values
#         self.sym_molprops_df = self.molprops_df.drop_duplicates(subset=['ppm', 'aromatic'], inplace=False)

#         # molecule has hose-code symmetry if there are less rows in sym_molprops_df than in molprops_df
#         self.has_hose_code_symmetry = self.sym_molprops_df.shape[0] < self.molprops_df.shape[0]

#         # calculate number of aromatic carbon atoms in sym_molprops_df
#         self.num_sym_aromatic_carbon_atoms = self.sym_molprops_df[self.sym_molprops_df.aromatic].shape[0]

#         # calculate the number of quaternary carbon atoms in sym_molprops_df
#         self.num_sym_quaternary_carbons = self.sym_molprops_df[self.sym_molprops_df.quaternary].shape[0]

#         # calculate the number of carbon atoms with two protons attached in sym_molprops_df
#         self.num_sym_ch2_carbon_atoms = self.sym_molprops_df[self.sym_molprops_df.CH2].shape[0]

#         # calculate the number of carbon atoms with three protons attached in sym_molprops_df
#         self.num_sym_ch3_carbon_atoms = self.sym_molprops_df[self.sym_molprops_df.CH3].shape[0]

#         # calculate the number of carbon atoms with one proton attached in sym_molprops_df
#         self.num_sym_ch_carbon_atoms = self.sym_molprops_df[self.sym_molprops_df.CH1].shape[0]

#         # calculate the number of carbon atoms with protons attached in sym_molprops_df
#         self.num_sym_carbon_atoms_with_protons = self.num_sym_ch2_carbon_atoms  + self.num_sym_ch3_carbon_atoms + self.num_sym_ch_carbon_atoms

#         # calculate number of carbon atoms in sym_molprops_df
#         self.num_sym_carbon_atoms = self.sym_molprops_df.shape[0]


#         # calculate number of aromatic rings
#         self.num_aromatic_rings = rdkit.Chem.rdMolDescriptors.CalcNumAromaticRings(self.mol)

#         # calculate number of all rings
#         self.num_rings = rdkit.Chem.rdMolDescriptors.CalcNumRings(self.mol)

#         # create a list of pairs of symmetry atoms if symmetry_idx1 is not -1
#         self.symmetry_pairs = []
#         for idx in self.molprops_df.index:
#             if self.molprops_df.loc[idx, "symmetry_idx1"] != -1:
#                 pair = {idx, self.molprops_df.loc[idx, "symmetry_idx1"]}
#                 # sort pair so that the lower index is first
#                 # pair.sort()
#                 # add pair to list if it is not already in the list
#                 if pair not in self.symmetry_pairs:
#                     self.symmetry_pairs.append(pair)

#     def map_symmetric_aromatic_rings(self, df):
#         # if there are no aromatic rings then return an empty list
#         if df.ring_idx.max() <= 0:
#             return []

#         # if there are more than one aromatic ring then map the aromatic rings to the symmetry pairs
#         # first get the symmetry pairs
#         symmetry_pairs = []
#         for x, y in np.array(np.triu_indices(df.ring_idx.max()+1, k=1)).T:
#             s1 = df.query("ring_idx == @x")["ppm"].to_list()
#             s2 = df.query("ring_idx == @y")["ppm"].to_list()
#             s1.sort()
#             s2.sort()
#             if s1 == s2:
#                 symmetry_pairs.append((x, y))

#         return symmetry_pairs


#     def GetRingSystems(self, includeSpiro=False):
#         ri = self.mol.GetRingInfo()
#         systems = []
#         for ring in ri.AtomRings():
#             ringAts = set(ring)
#             nSystems = []
#             for system in systems:
#                 nInCommon = len(ringAts.intersection(system))
#                 if nInCommon and (includeSpiro or nInCommon>1):
#                     ringAts = ringAts.union(system)
#                 else:
#                     nSystems.append(system)
#             nSystems.append(ringAts)
#             systems = nSystems
#         return systems

#     def _repr_png_(self):
#         return self.mol._repr_png_()

#     def _repr_svg_(self):
#         return self.mol._repr_svg_()

#     def __repr__(self):
#         return self.mol.__repr__()

#     def __str__(self):
#         return self.mol.__str__()

#     def GetAtoms(self):
#         return self.mol.GetAtoms()

#     def GetBonds(self):
#         return self.mol.GetBonds()

#     def GetAtomWithIdx(self, idx):
#         return self.mol.GetAtomWithIdx(idx)

#     def GetBondWithIdx(self, idx):
#         return self.mol.GetBondWithIdx(idx)

#     def GetNumAtoms(self):
#         return self.mol.GetNumAtoms()

#     def GetNumBonds(self):
#         return self.mol.GetNumBonds()

#     def GetAromaticAtoms(self):
#         return self.mol.GetAromaticAtoms()

#     def GetBondBetweenAtoms(self, i, j):
#         return self.mol.GetBondBetweenAtoms(i, j)

#     def GetRingInfo(self):
#         return self.mol.GetRingInfo()

#     def GetConformer(self):
#         return self.mol.GetConformer()

#     def init_elements_dict(self):
#         return (
#             pd.DataFrame(
#                 [
#                     [atom.GetIdx(), atom.GetSymbol()]
#                     for atom in Chem.AddHs(self.mol).GetAtoms()
#                 ],
#                 columns=["atom_index", "atom_symbol"],
#             )["atom_symbol"]
#             .value_counts()
#             .to_dict()
#         )

#     # calculate DBE for molecule
#     def calc_dbe(self) -> int:
#         elements = self.init_elements_dict()
#         if "C" in elements:
#             dbe_value = elements["C"]
#         if "N" in elements:
#             dbe_value += elements["N"] / 2
#         for e in ["H", "F", "Cl", "Br"]:
#             if e in elements:
#                 dbe_value -= elements[e] / 2

#         return dbe_value + 1

#     def calculated_c13_chemical_shifts(self) -> pd.DataFrame:
#         return self.calc_c13_chemical_shifts_using_nmrshift2D()

#     # return dictionary of dictionarys, first key is the number of protons, second key is carbon atom index, value is calculated C13 NMR chemical shift for molecule
#     def c13_nmr_shifts(self) -> dict:
#         c13_nmr_shifts = {
#             k: {i: None for i in v} for k, v in self.proton_groups().items()
#         }
#         print("c13_nmr_shifts", c13_nmr_shifts)

#         c13ppm_df = self.calc_c13_chemical_shifts_using_nmrshift2D()
#         if isinstance(c13ppm_df, pd.DataFrame):
#             # reset index to atom index
#             print("c13ppm_df.index", c13ppm_df.index)
#             print("c13ppm_df", c13ppm_df)
#             for v in c13_nmr_shifts.values():
#                 for k2, v2 in v.items():
#                     # find row in c13ppm_df with atom index k2
#                     v[k2] = c13ppm_df.loc[k2, "mean"]

#             print("c13_nmr_shifts", c13_nmr_shifts)
#         return c13_nmr_shifts

#     def calc_c13_chemical_shifts_using_nmrshift2D(self) -> pd.DataFrame:

#         with open("mol.mol", "w") as fp:
#             fp.write(Chem.MolToMolBlock(self.mol))

#         ret = os.system(expectedMolecule.JAVA_COMMAND)

#         if ret == 1:
#             print("NMRShift2D failed to calculate C13 chemical shifts")
#             return False
#         else:
#             mol_df = pd.read_csv("mol.csv", index_col=0)
#             mol_df.index = mol_df.index - 1
#             return mol_df

#     # return dictionary of lists key is the number of protons attached to carbon, value is list of carbon atom indices
#     def proton_groups(self) -> dict:
#         proton_groups = {
#             atom.GetTotalNumHs(): []
#             for atom in self.mol.GetAtoms()
#             if atom.GetAtomicNum() == 6
#         }
#         for atom in self.mol.GetAtoms():
#             if atom.GetAtomicNum() == 6:
#                 proton_groups[atom.GetTotalNumHs()].append(atom.GetIdx())
#         return proton_groups


#     def create_png(self) -> PIL.Image.Image:
#         """Creates a png image from a smiles string via rdkit"""
#         png = None

#         mol2 = Chem.AddHs(self.mol)
#         AllChem.EmbedMolecule(mol2, randomSeed=3)
#         rdkit_molecule = Chem.RemoveHs(mol2)

#         rdkit_molecule.Compute2DCoords()

#         return Draw.MolToImage(rdkit_molecule, size=(XYDIM, XYDIM))

#     def calc_carbon_xy_positions_png(self, rdkit_molecule: Chem.Mol) -> list:
#         """Returns the xy3 positions of the carbon atoms in the molecule"""

#         d2d = Draw.rdMolDraw2D.MolDraw2DSVG(XYDIM, XYDIM)
#         d2d.DrawMolecule(rdkit_molecule)
#         d2d.FinishDrawing()
#         idx_list = []
#         xxx = []
#         yyy = []
#         for atom in rdkit_molecule.GetAtoms():
#             if atom.GetSymbol() == "C":
#                 idx = atom.GetIdx()
#                 point = d2d.GetDrawCoords(idx)
#                 idx_list.append(idx)
#                 xxx.append(point.x / XYDIM)
#                 yyy.append(point.y / XYDIM)
#         return idx_list, xxx, yyy

if __name__ == "__main__":
    mol = NMRmol.from_smiles("CCC2Cc1ccccc1C2=O")
    print("mol.c13_nmr_shifts()\n", mol.c13_nmr_shifts())
    print("mol.num_carbon_atoms", mol.num_carbon_atoms)
    print("mol.num_carbon_atoms_with_protons", mol.num_carbon_atoms_with_protons)
    print("mol.num_quaternary_carbons", mol.num_quaternary_carbons)
    print("mol.num_ch_carbon_atoms", mol.num_ch_carbon_atoms)
    print("mol.num_ch2_carbon_atoms", mol.num_ch2_carbon_atoms)
    print("mol.num_ch3_carbon_atoms", mol.num_ch3_carbon_atoms)
