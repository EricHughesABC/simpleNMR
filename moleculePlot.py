import sys

# from PyQt5.QtWidgets import QWidget, QApplication, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QMessageBox
# from PyQt5.QtWidgets import QApplication
# from PyQt5.QtWebEngineWidgets import QWebEngineView
# from PyQt5.QtCore import QUrl
# from PyQt5.QtCore import pyqtSlot

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import mplcursors

# import rdkit
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem import Draw

# import networkx as nx
import nx_pylab
import numpy as np

# import pandas as pd

# import PIL
# from PIL import Image

import nmrProblem


class MatplotlibMoleculePlot(Figure):
    def __init__(self, nmrproblem):

        self.nmrproblem = nmrproblem

        self.mol_ind = None
        self.hmbc_ind = None
        # global ptcoords
        # # global hmbc_vertices_moved
        # global mol_vertices_moved
        self.label_id = None
        self.nmrproblem = nmrproblem
        super(MatplotlibMoleculePlot, self).__init__(figsize=(4, 4), dpi=100)
        self.ax = self.add_subplot(label="molecule", gid="molecule_id")
        self.ax.tick_params(
            axis="both",
            which="both",
            bottom=False,
            left=False,
            labelbottom=False,
            labelleft=False,
        )
        
        self.ax.spines["top"].set_visible(False)
        self.ax.spines["bottom"].set_visible(False)
        self.ax.spines["left"].set_visible(False)
        self.ax.spines["right"].set_visible(False)

        self.xmin = -1.5
        self.ymin = -1.5
        self.xmax = 1.5
        self.ymax = 1.5

        if self.nmrproblem.data_complete:
            self.draw_molecule(self.nmrproblem, self.ax)

    def draw_molecule(self, nmrproblem, ax):
        self.nmrproblem = nmrproblem
        self.ax = ax

        self.nmrproblem.hmbc_graph_edges = nmrProblem.create_hmbc_edges_dict(
            self.nmrproblem
        )
        self.nmrproblem.hmbc_graphs = nmrProblem.create_hmbc_graph_fragments(
            self.nmrproblem, self.nmrproblem.hmbc_graph_edges
        )

        for k in self.nmrproblem.hmbc_graphs:
            if None in self.nmrproblem.hmbc_graphs[k]['graph'].nodes():
                self.nmrproblem.hmbc_graphs[k]['graph'].remove_node(None)

        self.mol_edges, self.mol_nodes, self.mol_labels = self.init_molecule_plots(
            self.ax, self.nmrproblem
        )
        self.mol_vertices_moved = self.init_mol_vertices_moved(self.mol_edges)

        self.hmbc_graph_plots = self.init_hmbc_graph_plots(
            self.ax, self.nmrproblem.hmbc_graphs
        )

        # self.xmin, self.xmax = self.ax.get_xlim()
        # self.ymin, self.ymax = self.ax.get_ylim()

        # self.xmin = self.xmin - 0.2 * np.abs(self.xmin)
        # self.xmax = self.xmax + 0.2 * self.xmax

        # ydiff = self.ymax - self.ymin
        # xdiff = self.xmax - self.xmin

        # print("ydiff", ydiff)

        # self.ymin = self.ymin - 5*ydiff
        # self.ymax = self.ymax + 5*ydiff



        self.ax.set_xlim(self.xmin, self.xmax)
        self.ax.set_ylim(self.ymin, self.ymax)

        # print("self.xmin, self.xmax, self.ymin, self.ymax")
        # print(self.xmin, self.xmax, self.ymin, self.ymax)

        print(nmrproblem.xy3)

        # self.ax.set_xlim(self.xmax, self.xmin)
        # self.ax.set_ylim(self.ymin, self.ymax)

        print("nmrproblem.png", type(nmrproblem.png))

        if not isinstance(nmrproblem.png, type(None)):
            self.bkgnd = self.ax.imshow(
                np.fliplr(nmrproblem.png),
                aspect="auto",
                extent=[1.0*self.xmax, 1.0*self.xmin, 1.5*self.ymin, 1.5*self.ymax],
                alpha=0.4,
            )

        self.canvas.mpl_connect("pick_event", self.onpick3)
        self.canvas.mpl_connect("motion_notify_event", self.motion_notify_callback)
        self.canvas.mpl_connect("button_release_event", self.button_release_callback)

    def onpick3(self, event):

        self.mol_ind = event.ind
        ptcoords = np.ma.compressed(self.mol_nodes.get_offsets()[self.mol_ind])

        # find label id
        self.label_id = list(self.mol_labels.keys())[self.mol_ind[0]]

        # for each carbon atom in moleule check to see if it is in
        # the associated hmbc network then set the vertices_moved flag
        # to true if the picked node is in the network
        # update edges, nodes and labels
        # update them even if they are not visible so that keep in sync
        # with the moved display
        for n in self.nmrproblem.molecule.nodes:
            hmbc_nodes = self.hmbc_graph_plots[n]["hmbc_nodes"]
            hmbc_edges = self.hmbc_graph_plots[n]["hmbc_edges"]
            hmbc_labels = self.hmbc_graph_plots[n]["hmbc_labels"]
            hmbc_vertices_moved = self.hmbc_graph_plots[n]["hmbc_vertices_moved"]

            # check if picked node is in hmbc network and then
            # if the picked coords match of the edges match set them to True
            # so that the values of the edges will be moved during the motion notify
            # event
            if self.label_id in hmbc_labels.keys():
                self.hmbc_ind = [list(hmbc_labels.keys()).index(self.label_id)]
                if not isinstance(hmbc_edges, list):
                    for i, e in enumerate(hmbc_edges.get_paths()):
                        for j, c in enumerate(e.vertices):
                            if c[0] == ptcoords[0] and c[1] == ptcoords[1]:
                                hmbc_vertices_moved[i][j] = True
                            else:
                                hmbc_vertices_moved[i][j] = False
            else:
                # if the hmbc fragment does not contain the picked node
                # set all the moved_vertices to False
                if not isinstance(hmbc_edges, list):
                    for i, e in enumerate(hmbc_edges.get_paths()):
                        for j, c in enumerate(e.vertices):
                            hmbc_vertices_moved[i][j] = False
                    self.hmbc_ind = None

        # do the same for the molecule network
        if self.label_id in self.mol_labels.keys():
            for i, e in enumerate(self.mol_edges.get_paths()):
                for j, c in enumerate(e.vertices):
                    if c[0] == ptcoords[0] and c[1] == ptcoords[1]:
                        self.mol_vertices_moved[i][j] = True
                    else:
                        self.mol_vertices_moved[i][j] = False
        else:
            for i, e in enumerate(self.mol_edges.get_paths()):
                for j, c in enumerate(e.vertices):
                    self.mol_vertices_moved[i][j] = False
            # mol_vertices_moved = None

        # print(vertices_moved)

    def motion_notify_callback(self, event):

        if self.mol_ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return

        self.node_moved = True
        x, y = event.xdata, event.ydata

        # move nodes, edges and labels associated with hmbc network for each carbon atom
        for n in self.nmrproblem.molecule.nodes:
            if n in self.hmbc_graph_plots:

                hmbc_nodes = self.hmbc_graph_plots[n]["hmbc_nodes"]
                hmbc_edges = self.hmbc_graph_plots[n]["hmbc_edges"]
                hmbc_labels = self.hmbc_graph_plots[n]["hmbc_labels"]
                hmbc_vertices_moved = self.hmbc_graph_plots[n]["hmbc_vertices_moved"]

                if self.label_id in hmbc_labels.keys():
                    self.hmbc_ind = [list(hmbc_labels.keys()).index(self.label_id)]

                    # readjust coords of nodes in hmbc network
                    xy = np.asarray(hmbc_nodes.get_offsets())
                    xy[self.hmbc_ind[0]] = np.array([x, y])
                    hmbc_nodes.set_offsets(xy)

                    # readjust coords of labels in hmbc network
                    hmbc_labels[self.label_id].set_x(x)
                    hmbc_labels[self.label_id].set_y(y)

                    if not isinstance(hmbc_edges, list):
                        verts = []
                        for i, e in enumerate(hmbc_edges.get_paths()):
                            v = []
                            for j, c in enumerate(e.vertices):
                                if hmbc_vertices_moved[i][j] == True:
                                    v.append([x, y])
                                else:
                                    v.append(c)

                            verts.append(v)
                        # readjust coords of edges in hmbc network
                        hmbc_edges.set_verts(verts)

        # move nodes, edges and labels associated with molecule
        xy = np.asarray(self.mol_nodes.get_offsets())
        xy[self.mol_ind] = np.array([x, y])

        # readjust coords of nodes in molecule network
        self.mol_nodes.set_offsets(xy)

        # readjust coords of labels in molecule network
        self.mol_labels[self.label_id].set_x(x)
        self.mol_labels[self.label_id].set_y(y)

        verts = []

        for i, e in enumerate(self.mol_edges.get_paths()):
            v = []
            for j, c in enumerate(e.vertices):
                if self.mol_vertices_moved[i][j] == True:
                    v.append([x, y])
                else:
                    v.append(c)

            verts.append(v)

        # readjust coords of edges in molecule network
        self.mol_edges.set_verts(verts)

        self.ax.figure.canvas.draw_idle()

        self.moved_node = self.mol_ind
        self.moved_x = x
        self.moved_y = y

    def button_release_callback(self, event):
        # global mol_ind
        # global hmbc_ind
        # global node_moved
        # global mol_nodes

        if event.button != 1:
            return
        self.mol_ind = None
        self.hmbc_ind = None
        self.node_moved = False
        # print("button release: ind", mol_ind, hmbc_ind)
        # print("mol_nodes.get_offsets()\n", dict(zip(nmrproblem.molecule.nodes,mol_nodes.get_offsets())))
        # print("nmrproblem.xy3\n", nmrproblem.xy3)
        # self.nmrproblem.xy3 = dict(zip(self.nmrproblem.molecule.nodes,self.mol_nodes.get_offsets()))
        # hmbc_graph_edges = nmrProblem.create_hmbc_edges_dict(nmrproblem)
        # self.hmbc_graphs = create_hmbc_graph_fragments(self.nmrproblem, hmbc_graph_edges)

        # print("hmbc_graphs\n", hmbc_graphs)

    def init_hmbc_graph_plots(self, ax, hmbc_graphs):
        hmbc_graph_plots = {}

        for n in self.nmrproblem.molecule.nodes:
            hmbc_graph_plots[n] = {}
        for n in self.nmrproblem.molecule.nodes:
            if n not in hmbc_graphs:
                continue
            hmbc_graph_plots[n]["hmbc_nodes"] = nx_pylab.draw_networkx_nodes(
                hmbc_graphs[n]["graph"],
                hmbc_graphs[n]["xy"],
                ax=ax,
                label=n + "_node",
                node_color="w",
                edgecolors=["r"] + hmbc_graphs[n]["colors"],
                node_size=500,
                picker=False,
            )

            hmbc_graph_plots[n]["hmbc_nodes"].set_visible(False)

            hmbc_graph_plots[n]["hmbc_edges"] = nx_pylab.draw_networkx_edges(
                hmbc_graphs[n]["graph"],
                hmbc_graphs[n]["xy"],
                edge_color=hmbc_graphs[n]["colors"],
                width=5,
                ax=ax,
                label=n + "_edge",
            )

            hmbc_graph_plots[n]["hmbc_labels"] = nx_pylab.draw_networkx_labels(
                hmbc_graphs[n]["graph"], hmbc_graphs[n]["xy"], ax=ax
            )

            for k, v in hmbc_graph_plots[n]["hmbc_labels"].items():
                v.set_visible(False)

            if not isinstance(hmbc_graph_plots[n]["hmbc_edges"], list):
                hmbc_graph_plots[n]["hmbc_edges"].set_visible(False)

                hmbc_graph_plots[n]["hmbc_vertices_moved"] = []

                for e in hmbc_graph_plots[n]["hmbc_edges"].get_paths():
                    vertices = []
                    for c in e.vertices:
                        vertices.append(False)
                    hmbc_graph_plots[n]["hmbc_vertices_moved"].append(vertices)

            else:
                hmbc_graph_plots[n]["hmbc_vertices_moved"] = [False]

        return hmbc_graph_plots

    def init_mol_vertices_moved(self, mol_edges):
        mol_vertices_moved = []
        for e in mol_edges.get_paths():
            vertices = []
            for c in e.vertices:
                vertices.append(False)
            mol_vertices_moved.append(vertices)

        return mol_vertices_moved

    def init_molecule_plots(self, ax, nmrproblem):

        mol_edges = nx_pylab.draw_networkx_edges(
            nmrproblem.molecule,
            nmrproblem.xy3,
            ax=ax,
            edge_color="r",
            width=2,
            label="mol_edges",
        )

        mol_nodes = nx_pylab.draw_networkx_nodes(
            nmrproblem.molecule,
            nmrproblem.xy3,
            node_color=[nmrproblem.molecule.nodes[n]['node_color'] for n in nmrproblem.molecule.nodes],
            edgecolors="r",
            linewidths=0.2,
            node_size=500,
            ax=ax,
            label="mol_nodes",
            picker=True,
            pickradius=5,
        )

        mol_labels = nx_pylab.draw_networkx_labels(
            nmrproblem.molecule, nmrproblem.xy3, ax=ax
        )

        return mol_edges, mol_nodes, mol_labels
