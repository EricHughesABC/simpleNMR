# create nmrmol which is a subclass of rdkit.Chem.rdchem.Mol
import os
import platform

import pandas as pd

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

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

        if platform.system() == "Darwin":
            self.MAC_OS = True
            self.WINDOWS_OS = False
            if not os.system("jre/amazon-corretto-17.jdk/Contents/Home/bin/java --version"):
                self.JAVA_AVAILABLE = True
                self.JAVA_COMMAND = "jre/amazon-corretto-17.jdk/Contents/Home/bin/java -classpath predictorc.jar:cdk-2.7.1.jar:. NewTest mol.mol > mol.csv"
                print("MAC Local JAVA is available")
            else:
                self.JAVA_AVAILABLE = False
                print("JAVA is not available")

        elif platform.system() == "Windows":
            self.WINDOWS_OS = True
            self.MAC_OS = False
            # test if local windows ins installed
            if not os.system('"jre\\javawindows\\bin\\java -version"'):
                self.JAVA_AVAILABLE = True
                WINDOWS_OS = True
                self.JAVA_COMMAND = '"jre\\javawindows\\bin\\java -classpath predictorc.jar;cdk-2.7.1.jar;. NewTest mol.mol > mol.csv"'
                print("WINDOWS Local JAVA is available")
            else:
                self.JAVA_AVAILABLE = False
                print("JAVA is not available")


        self.dbe = self.calc_dbe()
        self.elements = self.init_elements_dict()
        self.num_carbon_atoms = self.elements.get("C", 0)
        self.num_hydrogen_atoms = self.elements.get("H", 0)



        molprops = [ [atom.GetIdx(), 
                    atom.GetNumImplicitHs(), 
                    atom.GetTotalNumHs(), 
                    atom.GetDegree(), 
                    atom.GetHybridization(), 
                    atom.GetIsAromatic()] for atom in self.GetAtoms() if 'C' == atom.GetSymbol()]

        self.molprops_df = pd.DataFrame(data=molprops, columns=['idx', 'implicitHs', 'totalNumHs', 'degree', 'hybridization', 'aromatic'])
        self.molprops_df = self.molprops_df.set_index(['idx'])
        

        c13_chemical_shifts_df = self.calculated_c13_chemical_shifts()

        # add c13 chemical shifts to molprops_df if available
        if isinstance(c13_chemical_shifts_df, pd.DataFrame):
            self.molprops_df = self.molprops_df.join(c13_chemical_shifts_df, how='left')

        self.molprops_df["c13ppm"] = self.molprops_df["mean"]

        self.molprops_df["CH2"] = False
        # set CH2 to True if carbon has 2 protons attached
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 2, "CH2"] = True

        # define quaternary carbon atoms column for totalNumHs == 0
        self.molprops_df["quaternary"] = False
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 0, "quaternary"] = True

        # define CH3 column for totalNumHs == 3
        self.molprops_df["CH3"] = False
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 3, "CH3"] = True

        # define CH1 column for totalNumHs == 1
        self.molprops_df["CH1"] = False
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 1, "CH1"] = True

        # calculate number of carbons with protons attached
        self.num_carbon_atoms_with_protons = self.molprops_df[self.molprops_df.totalNumHs > 0].shape[0]

        # calculate number of carbons without protons attached
        self.num_quartenary_carbons = self.molprops_df[self.molprops_df.totalNumHs == 0].shape[0]

        # calculate number of carbon with two protons attached
        self.num_ch2_carbon_atoms = self.molprops_df[self.molprops_df.totalNumHs == 2].shape[0]

        # calculate number of carbon with three protons attached
        self.num_ch3_carbon_atoms = self.molprops_df[self.molprops_df.totalNumHs == 3].shape[0]

        # calculate number of carbon with one proton  attached
        self.num_ch_carbon_atoms = self.molprops_df[self.molprops_df.totalNumHs == 1].shape[0]




    #initialize class from smiles string
    @classmethod
    def from_smiles(cls, smiles: str):
        return cls(Chem.MolFromSmiles(smiles))

    def init_elements_dict(self):
        return pd.DataFrame([[atom.GetIdx(), atom.GetSymbol() ] 
                                for atom in Chem.AddHs(self).GetAtoms()], 
                                    columns=["atom_index", "atom_symbol"])["atom_symbol"].value_counts().to_dict()

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
        df = pd.DataFrame([[atom.GetIdx(), atom.GetTotalNumHs()] 
                               for atom in self.GetAtoms() if atom.GetAtomicNum() == 6], 
                                  columns=["atom_index", "num_hydrogens"])
        df = df["num_hydrogens"].value_counts().sort_index()
        return df.to_dict()

    # return dictionary of lists key is the number of protons attached to carbon, value is list of carbon atom indices
    def proton_groups(self) -> dict:
        proton_groups = {}
        for atom in self.GetAtoms():
            if atom.GetAtomicNum() == 6:
                proton_groups[atom.GetTotalNumHs()] = []
        for atom in self.GetAtoms():
            if atom.GetAtomicNum() == 6:
                proton_groups[atom.GetTotalNumHs()].append(atom.GetIdx())
        return proton_groups

    def calculated_c13_chemical_shifts(self) -> pd.DataFrame:
        return self.calc_c13_chemical_shifts_using_nmrshift2D()

    # return dictionary of dictionarys, first key is the number of protons, second key is carbon atom index, value is calculated C13 NMR chemical shift for molecule
    def c13_nmr_shifts(self) -> dict:
        c13_nmr_shifts = {}
        for k,v in self.proton_groups().items():
            c13_nmr_shifts[k] = {}
            for i in v:
                c13_nmr_shifts[k][i] = None

        print("c13_nmr_shifts", c13_nmr_shifts)

        c13ppm_df = self.calc_c13_chemical_shifts_using_nmrshift2D()




        if isinstance(c13ppm_df, pd.DataFrame):
            # reset index to atom index
            print("c13ppm_df.index", c13ppm_df.index)
            print("c13ppm_df", c13ppm_df)
            for k, v in c13_nmr_shifts.items():
                for k2, v2 in v.items():
                    # find row in c13ppm_df with atom index k2
                    v[k2] = c13ppm_df.loc[k2,"mean"]

            print("c13_nmr_shifts", c13_nmr_shifts)
        return c13_nmr_shifts


    def calc_c13_chemical_shifts_using_nmrshift2D(self)->pd.DataFrame:

        with open("mol.mol", "w") as fp:
            fp.write(Chem.MolToMolBlock(self))

        ret = os.system(self.JAVA_COMMAND)

        if ret == 1:
            print("NMRShift2D failed to calculate C13 chemical shifts")
            return False
        else:
            mol_df = pd.read_csv('mol.csv', index_col=0)
            mol_df.index = mol_df.index - 1
            return mol_df  



    def num_protons_attached_to_carbons(self) -> int:
        num_protons = sum([k*v for k,v in self.return_proton_groups().items()])
        return num_protons

    def num_proton_groups(self) -> int:
        return sum(self.return_proton_groups().values())

    def aromatic_carbons(self) -> list:
        aromatic_carbons = []
        for atom in self.GetAtoms():
            if atom.GetIsAromatic():
                aromatic_carbons.append(atom.GetIdx())
        return aromatic_carbons

    # return list idx of carbon atoms in molecule
    def carbon_idx(self) -> list:
        return [atom.GetIdx() for atom in self.GetAtoms() if atom.GetAtomicNum() == 6]


    def query_molprops_df(self, column_id: str, value) -> pd.DataFrame:
        return self.molprops_df[self.molprops_df[column_id] == value][['implicitHs', 
                                                                       'degree', 
                                                                       'aromatic', 
                                                                       'hybridization', 
                                                                       'CH2',
                                                                       'c13ppm']]


if __name__ == "__main__":
    mol = NMRmol.from_smiles("CCC2Cc1ccccc1C2=O")
    print("mol.c13_nmr_shifts()\n", mol.c13_nmr_shifts())
    print("mol.num_carbon_atoms", mol.num_carbon_atoms)
    print("mol.num_carbon_atoms_with_protons", mol.num_carbon_atoms_with_protons)
    print("mol.num_quartenary_carbons", mol.num_quartenary_carbons)
    print("mol.num_ch_carbon_atoms", mol.num_ch_carbon_atoms)
    print("mol.num_ch2_carbon_atoms", mol.num_ch2_carbon_atoms)
    print("mol.num_ch3_carbon_atoms", mol.num_ch3_carbon_atoms)



        