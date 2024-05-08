import networkx as nx
from networkx.readwrite import json_graph
import nx_pylab
from matplotlib import axes
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba


import seaborn as sns
from pathlib import Path
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import json

# import mplcursors
import numpy as np
import pandas as pd

import jinja2
import webbrowser

class ProblemToD3html:

    def __init__(self, problem):
        self.nmrproblem = problem

        self.c13 = self.nmrproblem.c13
        self.hsqc = self.nmrproblem.hsqc
        self.hmbc = self.nmrproblem.hmbc
        self.h1 = self.nmrproblem.h1
        self.cosy = self.nmrproblem.cosy
        self.mol_df = self.nmrproblem.expected_molecule.molprops_df

        self.mol = Chem.Mol(self.nmrproblem.expected_molecule.mol)

        self.html_dir = Path(self.nmrproblem.problemDirectoryPath, "html")

        if not self.html_dir.exists():
            self.html_dir.mkdir()

        self.molWidth = 1000
        self.molHeight = 600

        self.svgWidth = 1200
        self.svgHeight = 700

        # fix coordinates of molecule before creating png

        self.mol = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(self.mol, randomSeed=3)
        self.mol = Chem.RemoveHs(self.mol)

        self.mol.Compute2DCoords()

        self.svg_str, self.new_xy3 = self.create_svg_string(self.nmrproblem, 
                                                            self.mol, 
                                                            molWidth=self.molWidth, 
                                                            molHeight=self.molHeight, 
                                                            svgWidth=self.svgWidth, 
                                                            svgHeight=self.svgHeight)

        self.graph = self.create_network_graph(self.nmrproblem, self.mol, self.new_xy3)

        data = json_graph.node_link_data(self.graph)

        # tidy up the data

        links_str = "["
        for link in data["links"]:
            links_str += "\t{},\n".format(link)
        links_str += "]"
        links_str = links_str.replace("'", "\"")
        links_str = links_str.replace("True", "true")

        nodes_str = "["
        for node in data["nodes"]:
            nodes_str += "\t{},\n".format(node)
        nodes_str += "]"
        nodes_str = nodes_str.replace("'", "")
        node_str = nodes_str.replace("True", "true")

        # read in template html file
        with open(r"html/d3molplot_template.html", "r") as f:
            html_template = f.read()

        environment = jinja2.Environment()

        template = environment.from_string(html_template)

        d3plot_html = template.render(
            svg_container = self.svg_str,
            graph_edges = links_str,
            graph_nodes = nodes_str,
            translateX = int((self.svgWidth - self.molWidth) / 2),
            translateY = int((self.svgHeight - self.molHeight) / 2),
            title = self.nmrproblem.problemDirectory,
        )

        d3plot_html_path = Path(self.nmrproblem.problemDirectoryPath, 
                                "html", 
                                self.nmrproblem.problemDirectory + "_d3.html")
        
        with open(d3plot_html_path, "w") as f:
            f.write(d3plot_html)

        webbrowser.open("file://" + d3plot_html_path.__str__())

    def create_svg_string(self, 
                          nmrproblem, 
                          mol, 
                          molWidth=500, 
                          molHeight=500, 
                          svgWidth=700, 
                          svgHeight=700):

        c13 = nmrproblem.c13

        translateWidth = int((svgWidth - molWidth) / 2)
        translateHeight = int((svgHeight - molHeight) / 2)

        d2d = rdMolDraw2D.MolDraw2DSVG(molWidth, molHeight)
        d2d.DrawMolecule(mol)
        d2d.TagAtoms(mol)
        d2d.FinishDrawing()
        sss = d2d.GetDrawingText().replace(f"width='{molWidth}px' height='{molHeight}px'", f"width={molWidth} height={molHeight}")
        sss = sss.replace("fill:#FFFFFF", "fill:none").replace("<svg", '<svg class="center"')

        sss = sss.replace(f"<!-- END OF HEADER -->", f"<!-- END OF HEADER -->\n<g transform='translate({translateWidth}, {translateHeight})'>")
        sss = sss.replace("</svg>", "</g>\n</svg>")
        sss = sss.replace(f"width={molWidth} height={molHeight} viewBox='0 0 {molWidth} {molHeight}'", f"width={svgWidth} height={svgHeight} viewBox='0 0 {svgWidth} {svgHeight}'")

        idx_list = []
        xxx = []
        yyy = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "C":
                idx = atom.GetIdx()
                point = d2d.GetDrawCoords(idx)
                idx_list.append(idx)
                xxx.append(point.x / molWidth)
                yyy.append(point.y / molHeight)
        new_xy3 = {}
        for index, row in c13.iterrows():
            value_to_find = row['C']
            i = idx_list.index(value_to_find)
            new_xy3[index] = (xxx[i], yyy[i])
        return sss, new_xy3
    
    def create_network_graph(self, 
                             nmrproblem, 
                             mol, 
                             new_xy3):

        c13 = nmrproblem.c13
        hsqc = nmrproblem.hsqc
        hmbc = nmrproblem.hmbc
        h1 = nmrproblem.h1
        cosy = nmrproblem.cosy

        # Create an empty graph
        graph = nx.Graph()

        # Add nodes to the graph using the index values of the c13 dataframe
        graph.add_nodes_from([int(v) for v in c13.index.values])

        for index, row in cosy.iterrows():
            f1p_i = row['f1p_i']
            f2p_i = row['f2p_i']

            if f1p_i == f2p_i:
                continue

            # check if f1p_i and f2p_i are in c13 index if not skip
            if f1p_i not in c13.index or f2p_i not in c13.index:
                continue

            if not graph.has_edge(f1p_i, f2p_i):
                graph.add_edge(f1p_i, f2p_i)
            graph.get_edge_data(f1p_i, f2p_i)["cosy"] = True

        # Add the edges for the HMBC using columns f1_i and f2p_i and add the attribute "hmbc" 
        # if the edge is already present in the graph just add the attribute "hmbc" to the edge. Skip if f1_i  equals f2p_i
        for index, row in hmbc.iterrows():
            f1_i = row['f1_i']
            f2p_i = row['f2p_i']

            if f1_i == f2p_i:
                continue

            # check if f1_i and f2p_i are in c13 index if not skip
            if f1_i not in c13.index or f2p_i not in c13.index:
                continue

            if not graph.has_edge(f1_i, f2p_i):
                graph.add_edge(f1_i, f2p_i)
            graph.get_edge_data(f1_i,f2p_i)['hmbc'] = True

        # # Add  fake NOESY edeges to the graph
        # if not graph.has_edge(9, 11):
        #     graph.add_edge(9, 11)
        # graph.get_edge_data(9, 11)['noesy'] = True

        # add numprotons to the nodes
        for index, row in c13.iterrows():
            numprotons = row['attached_protons']
            graph.nodes[int(index)]['numProtons'] = numprotons

        # add c13 ppm to the nodes
        for index, row in c13.iterrows():
            ppm = row['ppm']
            graph.nodes[int(index)]['ppm'] = ppm

        # add the coodinates from new_xy3 to the nodes as x and y
        for key, value in new_xy3.items():
            graph.nodes[key]['x'] = value[0]
            graph.nodes[key]['y'] = value[1]

        for index, row in c13.iterrows():
            numprotons = row['attached_protons']
            if numprotons: # numprotons = 1, 2 or 3
                h1_index =  list(hsqc[hsqc['f1_i'] == int(index)]["f2_i"])
                # convert h1_index into a list of strings
                h1_index = [num for num in h1_index]

                graph.nodes[int(index)]["H1_ppm"] = list(h1.loc[h1_index, "ppm"].values)
                graph.nodes[int(index)]["jCouplingClass"] = ['"' + l + '"' for l in list(h1.loc[h1_index, "jCouplingClass"].values)]
                # exctact the jcoupling values as a pandas datframe
                df = h1.loc[h1_index, "jCouplingVals"]
                # if  numrows = 1 then convert to a  list
                if df.shape[0] == 1:
                    graph.nodes[int(index)]["jCouplingVals"] = [ '"' + str(l)  + '"' for l in list(h1.loc[h1_index, "jCouplingVals"].values)]
                else:
                    # loop through the rows of the df keeping the jCouplingVals as a list of lists 
                    graph.nodes[int(index)]["jCouplingVals"] = ['"' + str(row) + '"' for index, row in df.iteritems()]

            else:
                graph.nodes[int(index)]["H1_ppm"] = []
                graph.nodes[int(index)]["jCouplingClass"] = []
                graph.nodes[int(index)]["jCouplingVals"] = []

        return graph