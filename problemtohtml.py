import jinja2

# import cairosvg
import webbrowser
import pandas as pd
import os
import re
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
# import Path
from pathlib import Path


class ProblemToHTML:
    def __init__(self, nmrproblem):
        self.nmrproblem = nmrproblem

        self.c13 = self.nmrproblem.c13
        self.hsqc = self.nmrproblem.hsqc
        self.h1 = self.nmrproblem.h1
        self.mol_df = self.nmrproblem.expected_molecule.molprops_df

        self.emol = Chem.Mol(self.nmrproblem.expected_molecule.mol)
        self.cmol = Chem.Mol(self.nmrproblem.expected_molecule.mol)
        self.pmol = Chem.Mol(self.nmrproblem.expected_molecule.mol)

        self.summary_df = pd.DataFrame(
            columns=[
                "predicted",
                "carbon",
                "proton",
                "integral",
                "J Class",
                "J Coupling",
            ]
        )

        # create html dir in problem dir
        # self.rootDirectory, self.problemDirectory
        # self.html_dir = os.path.join(self.nmrproblem.problemDirectoryPath, "html")
        # change to using Path
        
        self.html_dir = Path(self.nmrproblem.problemDirectoryPath, "html")

        # create png images and save them in the problem directory html directory

        # if not os.path.exists(self.html_dir):
        #     os.mkdir(self.html_dir)

        if not self.html_dir.exists():
            self.html_dir.mkdir()

        # add w3.css to html directory if it doesn't exist
        # w3css = os.path.join(self.html_dir, "w3.css")
        w3css = Path(self.html_dir, "w3.css")
        # if not os.path.exists(w3css):
        if not w3css.exists():
            import shutil

            shutil.copyfile(r"html/w3.css", w3css)

        # add notation to molecules
        # self._add_notation_to_mols()
        self.emol = self._annotate_rdkit_mol_with_ppm(
            self.emol, self.c13, self.hsqc, column_id="predicted", protons=False
        )
        self.cmol = self._annotate_rdkit_mol_with_ppm(
            self.cmol, self.c13, self.hsqc, column_id="ppm", protons=False
        )
        self.pmol = self._annotate_rdkit_mol_with_ppm(
            self.pmol, self.c13, self.hsqc, column_id="ppm", protons=True
        )

        # create svg strings for molecules
        self.emol_svg = self._create_svg_string(self.emol)
        self.cmol_svg = self._create_svg_string(self.cmol)
        self.pmol_svg = self._create_svg_string(self.pmol)

        # create summary table
        self.summary_df = pd.DataFrame(
            columns=[
               
                "carbon atom",
                 "atom idx",
                "predicted",
                "carbon",
                "proton",
                "predicted integral",
                "integral",
                "J Class",
                "J Coupling",
            ]
        )

        # create summary table
        self.summary_df = self._create_summary_table(
            self.summary_df, self.c13, self.hsqc, self.h1, self.mol_df
        )

        print(self.summary_df)

        allcarbons = self.c13.shape[0]
        allcarbonswithprotons = self.c13[self.c13.attached_protons > 0].shape[0]

        self.svg_str1 = self._add_atom_class_to_svg_string(
            self.emol_svg, self.summary_df, allcarbons, 0
        )
        self.svg_str2 = self._add_atom_class_to_svg_string(
            self.cmol_svg, self.summary_df, allcarbons, 0
        )
        self.svg_str3 = self._add_atom_class_to_svg_string(
            self.pmol_svg, self.summary_df, allcarbonswithprotons, 1
        )

    def _create_svg_string(self, mol, xdim=600, ydim=400):
        d2d = rdMolDraw2D.MolDraw2DSVG(xdim, ydim)
        d2d.DrawMolecule(mol)
        d2d.TagAtoms(mol)
        d2d.FinishDrawing()
        return d2d.GetDrawingText().replace("fill:#FFFFFF", "fill:none")

    def _annotate_rdkit_mol_with_ppm(
        self, mol, c13df, hsqcdf, column_id="ppm", protons=False
    ):

        if protons:
            for atom in mol.GetAtoms():

                if atom.GetSymbol() == "C":
                    idx = atom.GetIdx()
                    cppm = c13df.query("C == @idx")[column_id].values
                    if len(cppm) == 1:
                        cppm = cppm[0]
                    else:
                        continue

                    hppm = hsqcdf.query("f1_ppm == @cppm")["f2_ppm"].values
                    hppm_str = ", ".join([f"{x:.2f}" for x in hppm])
                    atom.SetProp("atomNote", hppm_str)
        else:
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "C":
                    idx = atom.GetIdx()
                    lbl = c13df.query("C == @idx")[column_id].values
                    if len(lbl) == 1:
                        atom.SetProp("atomNote", f"{lbl[0]:.1f}")

        return mol

    def _add_atom_class_to_svg_string(
        self, svg_str, summary_df, num_carbon_groups, num_protons=0
    ):
        """
        add atom class to svg string
        """
        # split svg string into lines
        svglines = svg_str.split("\n")

        # find the nodes/atoms and coordinates in the svg list
        circles = [
            dict(re.findall(r"(\w+)='([^']*)'", line))
            for line in svglines
            if "<circle" in line
        ]

        # find the anotations in the svg list
        # they are identified by the class attribute note
        # save the line number class and coordinates in a dictionary with the key the line number
        notes = {}
        for i, line in enumerate(svglines):
            if "note" in line:
                try:
                    # Extract the class and starting move coordinates using regular expressions
                    class_attr = re.search(r"class='([^']*)'", line).group(1)
                    start_coords = re.search(
                        r"M\s+([0-9.]+)\s+([0-9.]+)", line
                    ).groups()

                    notes[i] = {
                        "class": class_attr,
                        "x": float(start_coords[0]),
                        "y": float(start_coords[1]),
                    }
                except Exception as e:
                    print(e)

        # from  the circles extract the ones that have annotations based on the proton integrals
        # if proton ppm then integral must be greater than equal to 1

        pcircles = []
        pcircle_labels = []

        for circle in circles:
            # split class attribute on '-' keep only the last element
            atom_idx = int(circle["class"].split("-")[-1])

            # check if atom_idx is in the column atom idx of summary_df and print if predicted integral > 0
            if atom_idx in summary_df["atom idx"].values:
                if (
                    summary_df[summary_df["atom idx"] == atom_idx][
                        "predicted integral"
                    ].values[0]
                    >= num_protons
                ):
                    x, y = float(circle["cx"]), float(circle["cy"])
                    pcircles.append([x, y])
                    pcircle_labels.append(circle["class"])

        # set up kmeans clustering to associate the notes with the carbon atom
        features = [[v["x"], v["y"]] for v in notes.values()]

        # concatenate pcircles and features
        features = np.concatenate((pcircles, features))

        # standardize the features
        scaler = StandardScaler()
        scaled_features = scaler.fit_transform(features)

        # set up the kmeans clustering
        kmeans = KMeans(
            init="k-means++", n_clusters=num_carbon_groups, n_init=14, max_iter=300
        )

        # fit the kmeans clustering
        kmeans.fit(scaled_features)

        # start to assign the notes to the carbon atoms
        labels_df = pd.DataFrame(
            kmeans.labels_[num_carbon_groups:],
            columns=["labels"],
            index=list(notes.keys()),
        )
        labels_df["key"] = ""

        for line_no, atom_idx, k_idx in zip(
            list(notes.keys()), pcircle_labels, kmeans.labels_[:num_carbon_groups]
        ):
            labels_df.loc[labels_df.labels == k_idx, "key"] = f"{atom_idx.split()[-1]}"

        for idx in labels_df.index:
            svglines[idx] = svglines[idx].replace(
                "note", f"atom {labels_df.loc[idx, 'key']}"
            )

        return "\n".join(svglines)

    def _create_summary_table(self, summary_df, c13, hsqc, h1, mol_df):
        ### Create a summary table for report

        summary_df["predicted"] = c13["predicted"].values[::-1]
        summary_df["carbon"] = c13["ppm"].values[::-1]
        summary_df["atom idx"] = c13["C"].values[::-1]
        summary_df["carbon atom"] = list(c13.index.values[::-1])
        summary_df["proton"] = ""
        summary_df["predicted integral"] = ""
        summary_df["integral"] = ""
        summary_df["J Class"] = ""
        summary_df["J Coupling"] = ""

        # using hsqc dataframe add f2_ppm column to summary_df dataframe based on ppm column in c13 dataframe corresponding to f1_ppm column in hsqc dataframe
        for idx in summary_df.index:
            cppm = summary_df.loc[idx, "carbon"]
            protons = hsqc.query(f"f1_ppm == {cppm}")["f2_ppm"]

            if protons.empty:
                summary_df.loc[idx, "proton"] = ""
            else:
                summary_df.loc[idx, "proton"] = ",".join(
                    [f"{p:.2f}" for p in protons.values]
                )

        for idx in summary_df.index:
            cppm = summary_df.loc[idx, "carbon"]
            protons = hsqc.query(f"f1_ppm == {cppm}")["f2_ppm"]

            if protons.empty:
                summary_df.loc[idx, "integral"] = ""
            else:

                # find rows in h1 that have ppm values in ppm list
                ppm = protons.values.tolist()
                h1q = h1.query(f"ppm in {ppm}")

                if protons.empty:
                    summary_df.loc[idx, "integral"] = ""
                    summary_df.loc[idx, "J Class"] = ""
                    summary_df.loc[idx, "J Coupling"] = ""
                else:
                    summary_df.loc[idx, "integral"] = ", ".join(
                        [str(int(p)) for p in h1q.numProtons.values]
                    )
                    summary_df.loc[idx, "J Class"] = ", ".join(
                        [p for p in h1q.jCouplingClass.values]
                    )
                    summary_df.loc[idx, "J Coupling"] = ", ".join(
                        [f"[{p}]" for p in h1q.jCouplingVals.values if p != 0]
                    )

        # set carbon column to 1 decimal place
        summary_df["carbon"] = summary_df["carbon"].apply(lambda x: f"{x:.1f}")

        # add predicted integral from mol_df totalNumHs	column
        for idx in summary_df.index:
            cppm = summary_df.loc[idx, "predicted"]
            protons = mol_df.query(f"ppm == {cppm}")["totalNumHs"]
            summary_df.loc[idx, "predicted integral"] = protons.values[0]

        return summary_df

    def write_html_report(self):

        ### Create html page using jinja template

        summary_table = self.summary_df.to_html(table_id="summary_table", index=False)

        environment = jinja2.Environment()

        # read in template html file
        with open(r"html/highlight_table_template_002.html", "r") as f:
            html_template = f.read()

        template = environment.from_string(html_template)

        summary_html = template.render(
            filename=self.nmrproblem.problemDirectory,
            smiles=self.nmrproblem.smiles,
            svg1=self.svg_str1,
            svg2=self.svg_str2,
            svg3=self.svg_str3,
            table=summary_table,
        )

        # save summary html to file results_summary.html in directory html
        # html_file = os.path.join(
        #     self.html_dir, self.nmrproblem.problemDirectory + ".html"
        # )
        html_file = Path(self.html_dir, self.nmrproblem.problemDirectory + ".html")
        with open(html_file, "w") as f:
            f.write(summary_html)

        # open a web browser to display the html file

        # webbrowser.open("file://" + os.path.realpath(html_file))
        webbrowser.open("file://" + html_file.__str__())
        
