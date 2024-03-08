import sys
import os
import pathlib
import json
import webbrowser
import math

import PyQt5

from PyQt5 import QtCore
from PyQt5 import QtWidgets

from PyQt5.QtCore import QUrl
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtCore import Qt


from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5.QtWidgets import QLineEdit
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QSplitter
from PyQt5.QtWidgets import QMenu
from PyQt5.QtWidgets import QSpinBox
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QAction

from functools import partial

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)

from rdkit import Chem

# from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# import networkx as nx

# import nx_pylab
import numpy as np
import pandas as pd
from inspect import currentframe, getframeinfo

pd.options.mode.chained_assignment = "raise"
from PIL import Image
import platform

import nmrProblem

# from test_excel_simplenmr import TestExcelSimpleNMR
from qt5_tabs_001 import EditDataFrameDialog
from moleculePlot import MatplotlibMoleculePlot
from spectraPlot import MatplotlibH1C13Plot

from xy3_dialog import XY3dialog
from about_dialog import Aboutdialog
import problemtohtml
import java
import exampleProblems
from smiles_dialog import SmilesDialog

import simpleNMRutils
import mnovautils
from qtutils import warning_dialog


PLOTLINECOLORS = ("blue", "orange", "green", "red", "purple")
SCATTERFACECOLORS = ("blue", "orange", "green", "red", "purple")
SCATTEREDGECOLORS = ("blue", "orange", "green", "red", "purple")

PLOTLINECOLORS = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

SCATTERFACECOLORS = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

SCATTEREDGECOLORS = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


YELLOW = (1.0, 1.0, 0.0, 1.0)
RED = (1.0, 0.0, 0.0, 1.0)
WHITE = (1.0, 1.0, 1.0, 1.0)

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

# Test if we can run JAVA
# global JAVA_AVAILABLE
# global JAVA_COMMAND
# JAVA_AVAILABLE = False
# JAVA_COMMAND = "Undefined"

# Test what system we are running on
# global WINDOWS_OS
# global LINUX_OS
# global MAC_OS

# WINDOWS_OS = False
# LINUX_OS = False
# MAC_OS = False

# global XYDIM
# XYDIM = 800

# # set current working directory to the directory of this file
os.chdir(os.path.dirname(os.path.realpath(__file__)))


class MoleculePlotCanvas(FigureCanvasQTAgg):
    def __init__(self, fig, parent=None):
        super(MoleculePlotCanvas, self).__init__(fig)


class MainWidget(QMainWindow):
    def __init__(self, nmrproblem, parent=None):
        self.nmrproblem = nmrproblem
        self.hmbc_edge_colors = nmrproblem.hmbc_edge_colors
        self.old_node = ""
        self.new_node = ""
        self.selection = None

        super().__init__(parent)

        self.setGeometry(100, 100, 1200, 900)
        self.setWindowTitle(self.nmrproblem.problemDirectory)

        self._createActions()
        self._createMenuBar()
        # self._createToolBars()

        # Uncomment the call to ._createContextMenu() below to create a context
        # menu using menu policies. To test this out, you also need to
        # comment .contextMenuEvent() and uncomment ._createContextMenu()

        # self._createContextMenu()

        self._connectActions()
        self._createStatusBar()

        self.centralWidget = QWidget()
        # self.centralWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
        self.setCentralWidget(self.centralWidget)

        if self.nmrproblem.data_complete:

            print("line 142")
            print("nmrproblem.c13")
            print(nmrproblem.c13)
            print("nmrproblem.h1")
            print(self.nmrproblem.h1)
            print("nmrproblem.hsqc")
            print(self.nmrproblem.hsqc)

            nmrproblem.save_dataframes_in_nmrproblem_to_json()

            self.initiate_windows(self.nmrproblem)

    def initiate_windows(self, nmrproblem):

        del self.centralWidget

        self.centralWidget = QWidget()
        # self.centralWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
        self.setCentralWidget(self.centralWidget)

        self.nmrproblem = nmrproblem

        jsmewidget = QWidget()
        moleculewidget = QWidget()
        spectrawidget = QWidget()

        self.moleculePlot = MatplotlibMoleculePlot(self.nmrproblem)
        self.spectraPlot = MatplotlibH1C13Plot(self.nmrproblem)

        self.spectraCanvas = MoleculePlotCanvas(self.spectraPlot)
        self.spectraCanvas.setGeometry(0, 0, 600, 1100)
        self.moleculeCanvas = MoleculePlotCanvas(self.moleculePlot)
        self.moleculeCanvas.setGeometry(0, 0, 300, 300)

        self.moltoolbar = NavigationToolbar(self.moleculeCanvas, self)
        self.spctoolbar = NavigationToolbar(self.spectraCanvas, self)

        # self.webEngineView = QWebEngineView()
        # self.loadPage()

        # self.smilesInput = QLineEdit()
        # self.button = QPushButton("Display Smiles Molecule", self)

        self.sync_button = QPushButton("Synchronize with MesReNova", self)
        self.sync_button.clicked.connect(self.on_sync_click)

        # connect button to function on_click
        # self.button.clicked.connect(self.on_click)

        splitter1 = QSplitter(Qt.Vertical)
        splitter2 = QSplitter(Qt.Horizontal)

        hbox = QHBoxLayout()
        molvbox = QVBoxLayout()
        jsmevbox = QVBoxLayout()
        spcvbox = QVBoxLayout()

        spcvbox.addWidget(self.spctoolbar)
        spcvbox.addWidget(self.spectraCanvas)

        spectrawidget.setLayout(spcvbox)

        molvbox.addWidget(self.moltoolbar)
        molvbox.addWidget(self.moleculeCanvas)

        moleculewidget.setLayout(molvbox)

        jsmevbox.addWidget(self.sync_button)
        # jsmevbox.addWidget(self.smilesInput)
        # jsmevbox.addWidget(self.button)

        jsmewidget.setLayout(jsmevbox)

        splitter1.addWidget(moleculewidget)
        splitter1.addWidget(jsmewidget)

        splitter2.addWidget(splitter1)
        splitter2.addWidget(spectrawidget)

        hbox.addWidget(splitter2)
        # hbox.addLayout(molvbox)
        self.centralWidget.setLayout(hbox)
        splitter2.setSizes([600, 600])

        self.centralWidget.show()

        # add callbacks to the moleculeCanvas
        self.node_pick_ind = None
        self.node_hover_ind = None
        self.node_picked = False
        self.highlighted_peak_lbl = None
        self.old_lbl = None

        self.moleculePlot.canvas.mpl_connect(
            "button_release_event",
            lambda event: self.button_release_molecule(
                event, specplot=self.spectraPlot, molplot=self.moleculePlot
            ),
        )

        self.moleculePlot.canvas.mpl_connect(
            "motion_notify_event",
            lambda event: self.motion_notify_callback(
                event, specplot=self.spectraPlot, molplot=self.moleculePlot
            ),
        )

        self.moleculePlot.canvas.mpl_connect(
            "pick_event",
            lambda event: self.pick_molecule(
                event,
                event_name="pick_event",
                specplot=self.spectraPlot,
                molplot=self.moleculePlot,
            ),
        )

        self.spectraPlot.canvas.mpl_connect(
            "motion_notify_event",
            lambda event: self.hover_over_specplot(
                event, specplot=self.spectraPlot, molplot=self.moleculePlot
            ),
        )

        self.setWindowTitle(self.nmrproblem.problemDirectory)
        self.wcLabel.setText(
            f"Molecule: {self.nmrproblem.moleculeAtomsStr} DBE: {int(self.nmrproblem.dbe)}"
        )

        # self.spectraPlot.canvas.draw_idle()
        # self.moleculePlot.canvas.draw_idle()

    def hover_over_specplot(self, event, specplot, molplot):

        # just return if toolbar is active
        if (self.moltoolbar.mode != "") or (self.spctoolbar.mode != ""):
            return

        in_plot = []
        in_plot_label = []
        in_plot_index = []
        pos = None
        self.mol_nodes = molplot.mol_nodes
        for k, v in specplot.peak_overlays_dict.items():
            # in_c13plots, c13plots_index = v["highlight"].contains(event)
            in_c13plots, c13plots_index = v.contains(event)
            in_plot.append(in_c13plots)
            in_plot_label.append(k)
            in_plot_index.append(c13plots_index)

        if any(in_plot):
            lbl = in_plot_label[in_plot.index(True)]

            # if lbl != self.highlighted_peak_lbl:
            # highlight new peak
            specplot.peak_overlays_dict[lbl].set_visible(True)
            specplot.peak_overlays_dict[lbl].set_color(RED)
            specplot.peak_overlays_dict[lbl].set_linewidth(1.2)

            # annotate new peak
            if "H" in lbl:
                # set the annotation to the peak
                atom_index = int(lbl[1:])
                ppm = self.nmrproblem.h1.loc[atom_index, "ppm"]
                integral = float(self.nmrproblem.h1.loc[atom_index, "integral"])
                jcoupling = self.nmrproblem.h1.loc[atom_index, "jCouplingClass"]
                jcouplingvals = self.nmrproblem.h1.loc[atom_index, "jCouplingVals"]
                jcouplingvals = simpleNMRutils.stringify_vals(jcouplingvals)
                annot_text = f"{lbl}: {ppm:.2f} ppm\nInt:{integral:.1f}\nJ: {jcoupling}: {jcouplingvals} Hz"
                specplot.annot_H1.xy = (event.xdata, event.ydata)
                specplot.annot_H1.set_text(annot_text)
                specplot.annot_H1.set_visible(True)

            elif "C" in lbl:
                # set the annotation to the peak
                atom_index = int(lbl[1:])
                ppm = self.nmrproblem.c13.loc[atom_index, "ppm"]
                annot_text = f"{lbl}: {ppm:.1f} ppm"
                x = event.xdata
                y = event.ydata
                specplot.annot_C13.set_text(annot_text)
                specplot.annot_C13.xy = (x, y)
                specplot.annot_C13.set_visible(True)

            self.highlighted_peak_lbl = lbl

            if "H" in lbl:
                try:
                    clbl = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f2H_i == lbl][
                        "f1C_i"
                    ].values[0]
                except IndexError:
                    print(f"no HSQC peak found for this H1 peak {lbl}")
                    clbl = ""

            else:
                clbl = lbl

            if "C" not in clbl:
                specplot.canvas.draw_idle()
                return

            self.node_hover_lbl = clbl

            # update molplot highlights
            self.update_molplot_highlights(molplot, specplot, clbl)

            # uddate specplot canvas
            specplot.canvas.draw_idle()

            # uddate title in molplot
            c13_ind = int(clbl[1:])
            ppm_val = self.nmrproblem.c13.loc[c13_ind]["ppm"]
            attached_protons = self.nmrproblem.c13.loc[c13_ind]["attached_protons"]
            title_str = (
                f"{clbl}: {ppm_val:.1f} ppm, attached protons: {attached_protons}"
            )
            self.moleculePlot.ax.set_title(title_str)
        else:
            # unhilight old peak
            if self.highlighted_peak_lbl is not None:
                self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
                self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

                self.moleculePlot.hide_hmbc_graph_networks()

                specplot.reset_peak_overlays_eeh()

                specplot.hide_annotation(specplot.annot_C13)
                specplot.hide_annotation(specplot.annot_H1)

                # unhighlight distributions
                specplot.reset_distributions_eeh()

                # uddate specplot canvas
                specplot.canvas.draw_idle()
                molplot.canvas.draw_idle()

    def update_molplot_highlights(self, molplot, specplot, lbl):

        molplot.mol_nodes.node_highlighted = True

        ind = int(lbl[1:]) - 1

        scatter_facecolors_highlight = molplot.mol_nodes.get_facecolors()
        scatter_edgecolors_highlight = molplot.mol_nodes.get_edgecolors()

        scatter_facecolors_highlight[ind] = WHITE
        scatter_edgecolors_highlight[ind] = RED
        molplot.mol_nodes.set_fc(scatter_facecolors_highlight)

        if lbl in self.nmrproblem.hmbc_graphs.keys():
            # self.hide_hmbc_graph_networks()
            self.moleculePlot.draw_hmbc_graph_network(lbl)
            # higlight C13 hmbc peaks
            specplot.highlight_hmbc_C13_peaks(lbl)
            # self.redraw_axes()
            # highlight H1 hmbc peaks
            specplot.highlight_H1_HMBC_peaks(lbl)

        # highlight C13 peak in graph x1
        specplot.highlight_C13_peak(lbl)

        # higlight H1 peaks in graph x1
        specplot.highlight_H1_peaks_from_highlighted_carbon_atom(lbl)
        # self.c13_plots[lbl]['highlight'].set_visible(True)

        # annotate C13 peak in graph x1
        specplot.display_annotation_C13_from_molplot(lbl, specplot.annot_C13)

        # annotate H1 peaks in graph x1
        specplot.display_annotation_H1_from_molplot(
            lbl, specplot.annot_H1, self.nmrproblem
        )

        # annotate distributions
        hpks = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f1C_i == lbl]["f2H_i"].values
        cpks = [lbl]
        self.display_distributions(cpks, hpks, specplot)

        molplot.canvas.draw_idle()
        specplot.canvas.draw_idle()

    def display_distributions(self, C_pks, H_pks, specplot):
        # add highlighted distributions
        colors = ["b", "g", "r", "c", "m", "y", "k"]

        # used to dsplay legends of highlighted distributions
        plinesH = []
        plabelsH = []
        plinesC = []
        plabelsC = []
        ppmH = []
        ppmC = []

        # set visible
        # circle through colors
        # create proton distribution legends
        for pk in H_pks:
            numplots = len(specplot.h1c13distlist[0][pk]) - 1
            for i, aa in enumerate(specplot.h1c13distlist[0][pk]):
                ppmH.append(self.nmrproblem.udic[0]["info"].loc[pk, "ppm"])
                # sel.extras.append(self.cursor.add_highlight(aa))
                aa.set_visible(True)
                aa.set_linewidth(0.75)
                aa.set_color(colors[i % 7])
                # do not add legend info if plot is just single line showing where peak pos is
                if i < numplots:
                    plabelsH.append(aa.get_label())
                    plinesH.append(aa)

        # create carbon distribution legends
        for pk in C_pks:
            numplots = len(specplot.h1c13distlist[1][pk]) - 1
            for i, aa in enumerate(specplot.h1c13distlist[1][pk]):
                ppmC.append(self.nmrproblem.udic[1]["info"].loc[pk, "ppm"])
                # sel.extras.append(self.cursor.add_highlight(aa))
                aa.set_visible(True)
                aa.set_linewidth(0.75)
                aa.set_color(colors[i % 7])
                if i < numplots:
                    plabelsC.append(aa.get_label())
                    plinesC.append(aa)

        # adjust x-axis width H1 and C13 of distribution plots
        # and add legend information
        if len(ppmH) > 0:
            # calculate average position of peaks which will be used to adjust the x axis ppm range
            ppmmmH = np.mean(ppmH)
            specplot.h1dist_ax.set_xlim(ppmmmH + 2, ppmmmH - 2)
            specplot.h1dist_ax.legend(plinesH, plabelsH)
        if len(ppmC) > 0:
            ppmmmC = np.mean(ppmC)
            specplot.c13dist_ax.set_xlim(ppmmmC + 50, ppmmmC - 50)
            specplot.c13dist_ax.legend(plinesC, plabelsC)

    def button_release_molecule(self, event, molplot, specplot):

        if (self.moltoolbar.mode != "") or (self.spctoolbar.mode != ""):
            return

        self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
        self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

        self.moleculePlot.hide_hmbc_graph_networks()

        specplot.reset_peak_overlays_eeh()

        specplot.hide_annotation(specplot.annot_C13)
        specplot.hide_annotation(specplot.annot_H1)

        # unhighlight distributions
        specplot.reset_distributions_eeh()

        molplot.canvas.draw_idle()
        specplot.canvas.draw_idle()

        self.node_pick_ind = None
        self.node_hover_ind = None
        self.node_picked = False
        self.node_moved = False

    def motion_notify_callback(self, event, specplot, molplot):

        if (self.moltoolbar.mode != "") or (self.spctoolbar.mode != ""):
            return

        self.hover_over_molecule(
            event, event_name="motion_notify_event", molplot=molplot, specplot=specplot
        )

        if not hasattr(self, "label_id"):
            # print("not hasattr(self, 'label_id')")
            return
        # else:
        # print("self.label_id: ", self.label_id)

        in_node, _ = self.mol_nodes.contains(event)
        if (not in_node) and (event.button != 1):

            self.mol_nodes.node_highlighted = False

            self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
            self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

            self.moleculePlot.hide_hmbc_graph_networks()

            specplot.reset_peak_overlays_eeh()

            specplot.hide_annotation(specplot.annot_C13)
            specplot.hide_annotation(specplot.annot_H1)

            # unhighlight distributions
            specplot.reset_distributions_eeh()

            molplot.canvas.draw_idle()
            specplot.canvas.draw_idle()
            return

        if event.button != 1:
            # print("event.button != 1, return")
            return

        self.node_moved = True
        x, y = event.xdata, event.ydata

        # if x or y is None: return mouse beyond plotting area
        if x is None:
            return
        if y is None:
            return

        self.hmbc_graph_plots = self.moleculePlot.hmbc_graph_plots
        self.mol_nodes = self.moleculePlot.mol_nodes
        self.mol_edges = self.moleculePlot.mol_edges
        self.mol_labels = self.moleculePlot.mol_labels
        self.mol_ind = self.node_pick_ind

        # move nodes, edges and labels associated with hmbc network for each carbon atom
        for n in self.nmrproblem.nx_graph_molecule.nodes:
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

        if not isinstance(self.mol_edges, list):

            for i, e in enumerate(self.mol_edges.get_paths()):
                v = []
                for j, c in enumerate(e.vertices):
                    if self.moleculePlot.mol_vertices_moved[i][j] == True:
                        v.append([x, y])
                    else:
                        v.append(c)

                verts.append(v)

            # readjust coords of edges in molecule network
            self.mol_edges.set_verts(verts)

        self.moleculePlot.ax.figure.canvas.draw_idle()

        self.moved_node = self.mol_ind
        self.moved_x = x
        self.moved_y = y

    def hover_over_molecule(self, event, event_name, molplot, specplot):

        self.mol_nodes = molplot.mol_nodes
        self.mol_labels = molplot.mol_labels

        in_node, node_index = self.mol_nodes.contains(event)

        if self.node_picked:
            self.mol_nodes.node_highlighted = False
            return

        if in_node:

            self.mol_nodes.node_highlighted = True

            ind = node_index["ind"][0]
            lbl = f"{self.mol_nodes.my_labels[ind]}"
            x, y = self.mol_nodes.get_offsets()[ind]
            self.node_hover_ind = ind
            self.node_hover_lbl = lbl
            self.node_hover_x = x
            self.node_hover_y = y

            if lbl != self.old_lbl:
                self.old_lbl = lbl
                self.mol_nodes.node_highlighted = False

                self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
                self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

                self.moleculePlot.hide_hmbc_graph_networks()

                specplot.reset_peak_overlays_eeh()

                specplot.hide_annotation(specplot.annot_C13)
                specplot.hide_annotation(specplot.annot_H1)

                # unhighlight distributions
                specplot.reset_distributions_eeh()

                molplot.canvas.draw_idle()
                specplot.canvas.draw_idle()

            c13_ind = int(lbl[1:])
            ppm_val = self.nmrproblem.c13.loc[c13_ind]["ppm"]
            attached_protons = self.nmrproblem.c13.loc[c13_ind]["attached_protons"]
            title_str = (
                f"{lbl}: {ppm_val:.1f} ppm, attached protons: {attached_protons}"
            )
            self.moleculePlot.ax.set_title(title_str)

            scatter_facecolors_highlight = self.mol_nodes.get_facecolors()
            scatter_edgecolors_highlight = self.mol_nodes.get_edgecolors()

            scatter_facecolors_highlight[ind] = WHITE
            scatter_edgecolors_highlight[ind] = RED
            self.mol_nodes.set_fc(scatter_facecolors_highlight)

            if lbl in self.nmrproblem.hmbc_graphs.keys():
                # self.hide_hmbc_graph_networks()
                self.moleculePlot.draw_hmbc_graph_network(lbl)
                # higlight C13 hmbc peaks
                specplot.highlight_hmbc_C13_peaks(lbl)
                # self.redraw_axes()
                # highlight H1 hmbc peaks
                specplot.highlight_H1_HMBC_peaks(lbl)

            # highlight C13 peak in graph x1
            specplot.highlight_C13_peak(lbl)

            # higlight H1 peaks in graph x1
            specplot.highlight_H1_peaks_from_highlighted_carbon_atom(lbl)
            # self.c13_plots[lbl]['highlight'].set_visible(True)

            # annotate C13 peak in graph x1
            specplot.display_annotation_C13_from_molplot(lbl, specplot.annot_C13)

            # annotate H1 peaks in graph x1
            specplot.display_annotation_H1_from_molplot(
                lbl, specplot.annot_H1, self.nmrproblem
            )

            # annotate distributions
            hpks = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f1C_i == lbl][
                "f2H_i"
            ].values
            cpks = [lbl]
            self.display_distributions(cpks, hpks, specplot)

            molplot.canvas.draw_idle()
            specplot.canvas.draw_idle()
        else:
            # if self.mol_nodes.node_highlighted:
            self.mol_nodes.node_highlighted = False

            self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
            self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

            self.moleculePlot.hide_hmbc_graph_networks()

            specplot.reset_peak_overlays_eeh()

            specplot.hide_annotation(specplot.annot_C13)
            specplot.hide_annotation(specplot.annot_H1)

            # unhighlight distributions
            specplot.reset_distributions_eeh()

            molplot.canvas.draw_idle()
            specplot.canvas.draw_idle()

    def pick_molecule(self, event, event_name, specplot, molplot):
        global JAVA_AVAILABLE
        in_node, node_index = self.moleculePlot.mol_nodes.contains(event.mouseevent)
        print("in_node: ", in_node)
        print("node_index: ", node_index)   
        print("event_name: ", event_name)
        print("event.mouseevent: ", event.mouseevent)
        print(dir(event))

        if in_node:
            self.node_picked = True

        # self.mol_ind = event.ind
        # self.node_pick_ind = event.ind
        self.node_pick_ind = node_index["ind"][0]
        ptcoords = np.ma.compressed(
            self.moleculePlot.mol_nodes.get_offsets()[self.node_pick_ind]
        )

        # find label id
        self.label_id = list(self.moleculePlot.mol_labels.keys())[self.node_pick_ind]
        print("self.label_id: ", self.label_id, self.node_pick_ind, f"C{self.node_pick_ind+1}")

        # hide all highlights before dragging node
        if self.mol_nodes.node_highlighted:

            self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
            self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

            self.moleculePlot.hide_hmbc_graph_networks()

            specplot.reset_peak_overlays_eeh()

            specplot.hide_annotation(specplot.annot_C13)
            specplot.hide_annotation(specplot.annot_H1)

            # unhighlight distributions
            specplot.reset_distributions_eeh()
            specplot.canvas.draw_idle()
            molplot.canvas.draw_idle()

        # for each carbon atom in moleule check to see if it is in
        # the associated hmbc network then set the vertices_moved flag
        # to true if the picked node is in the network
        # update edges, nodes and labels
        # update them even if they are not visible so that keep in sync
        # with the moved display
        for n in self.nmrproblem.nx_graph_molecule.nodes:
            hmbc_nodes = self.moleculePlot.hmbc_graph_plots[n]["hmbc_nodes"]
            hmbc_edges = self.moleculePlot.hmbc_graph_plots[n]["hmbc_edges"]
            hmbc_labels = self.moleculePlot.hmbc_graph_plots[n]["hmbc_labels"]
            hmbc_vertices_moved = self.moleculePlot.hmbc_graph_plots[n][
                "hmbc_vertices_moved"
            ]

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
        if self.label_id in self.moleculePlot.mol_labels.keys():
            if not isinstance(self.moleculePlot.mol_edges, list):
                for i, e in enumerate(self.moleculePlot.mol_edges.get_paths()):
                    for j, c in enumerate(e.vertices):
                        if c[0] == ptcoords[0] and c[1] == ptcoords[1]:
                            self.moleculePlot.mol_vertices_moved[i][j] = True
                        else:
                            self.moleculePlot.mol_vertices_moved[i][j] = False
        else:
            if not isinstance(self.moleculePlot.mol_edges, list):
                for i, e in enumerate(self.moleculePlot.mol_edges.get_paths()):
                    for j, c in enumerate(e.vertices):
                        self.moleculePlot.mol_vertices_moved[i][j] = False

    def _createMenuBar(self):
        menuBar = self.menuBar()
        # File menu
        fileMenu = QMenu("&File", self)
        menuBar.addMenu(fileMenu)
        # fileMenu.addAction(self.newAction)
        fileMenu.addAction(self.newFromMresNovaAction)
        # fileMenu.addAction(self.newFromTopspinAction)
        fileMenu.addAction(self.openAction)
        # Open Recent submenu
        # self.openRecentMenu = fileMenu.addMenu("Open Recent")
        fileMenu.addAction(self.saveAction)
        # Separator
        fileMenu.addSeparator()
        fileMenu.addAction(self.exitAction)

        # # exampleProblems menu
        # exampleProblemsMenu = menuBar.addMenu("&Example Problems")

        examples_menu = menuBar.addMenu("Examples")
        self.initExampleMenu(examples_menu)

        # for  action in self.exampleProblemsActions:
        #     exampleProblemsMenu.addAction(action)

        # Edit menu
        # editMenu = menuBar.addMenu("&Edit")
        # editMenu.addAction(self.copyAction)
        # editMenu.addAction(self.pasteAction)
        # editMenu.addAction(self.cutAction)

        # Separator
        # editMenu.addSeparator()
        # Find and Replace submenu
        # findMenu = editMenu.addMenu("Find and Replace")
        # findMenu.addAction("Find...")
        # findMenu.addAction("Replace...")

        # Help menu
        # helpMenu = menuBar.addMenu(QIcon(":help-content.svg"), "&Help")
        helpMenu = menuBar.addMenu("&Help")
        helpMenu.addAction(self.helpContentAction)
        helpMenu.addAction(self.aboutAction)

    def initExampleMenu(self, menu):

        self.exampleProblems = exampleProblems.exampleproblems_names

        for actionName in self.exampleProblems:
            qa = QAction(actionName, self)
            menu.addAction(qa)
            qa.triggered.connect(self.exampleProblem_selected)

    def exampleProblem_selected(self):
        action = self.sender()
        exampleproblem_name = action.text()

        self.excel_path, self.tmpdir = exampleProblems.create_tmp_excel_file(
            exampleproblem_name, exampleProblems.exampleproblems_df[exampleproblem_name]
        )
        data_info = nmrProblem.parse_argv(
            my_argv=[
                "simplepy",
                self.tmpdir.name,
            ]
        )

        xy3_dlg = XY3dialog(
            java_available=java.JAVA_AVAILABLE,
            sheets_missing=nmrProblem.get_missing_sheets(data_info["excel_fn"], True),
        )

        if xy3_dlg.exec():
            xy3_calc_method, expts_available = xy3_dlg.get_method()
            # add "molecule" to expts_available
            expts_available.append("molecule")

            self.nmrproblem = nmrProblem.NMRproblem(
                data_info,
                java_available=java.JAVA_AVAILABLE,
                xy3_calc_method=xy3_calc_method,
                java_command=java.JAVA_COMMAND,
                expts_available=expts_available,
            )

            if self.nmrproblem.data_complete:

                print("line 918")
                print("nmrproblem.c13")
                print(self.nmrproblem.c13)
                print("nmrproblem.h1")
                print(self.nmrproblem.h1)
                print("nmrproblem.hsqc")
                print(self.nmrproblem.hsqc)
                self.nmrproblem.initiate_df_data_complete()
                self.initiate_windows(self.nmrproblem)
        else:
            print("xy3_dlg.exec() failed")

    def _createToolBars(self):
        # File toolbar
        fileToolBar = self.addToolBar("File")
        fileToolBar.setMovable(False)
        # fileToolBar.addAction(self.newAction)
        fileToolBar.addAction(self.newFromMresNovaAction)
        # fileToolBar.addAction(self.newFromTopspinAction)
        fileToolBar.addAction(self.openAction)
        fileToolBar.addAction(self.saveAction)
        # Edit toolbar
        # editToolBar = QToolBar("Edit", self)
        # self.addToolBar(editToolBar)
        # editToolBar.addAction(self.copyAction)
        # editToolBar.addAction(self.pasteAction)
        # editToolBar.addAction(self.cutAction)
        # Widgets
        self.fontSizeSpinBox = QSpinBox()
        self.fontSizeSpinBox.setFocusPolicy(Qt.NoFocus)
        # editToolBar.addWidget(self.fontSizeSpinBox)

    def _createStatusBar(self):
        self.statusbar = self.statusBar()
        # Temporary message
        self.statusbar.showMessage("Ready", 3000)
        # Permanent widget
        self.wcLabel = QLabel(
            f"Molecule: {self.nmrproblem.moleculeAtomsStr} DBE: {int(self.nmrproblem.dbe)}"
        )
        self.statusbar.addPermanentWidget(self.wcLabel)

    def _createActions(self):
        # File actions
        self.newAction = QAction(self)
        self.newAction.setText("&New")

        self.newFromMresNovaAction = QAction(self)
        self.newFromMresNovaAction.setText("New from MNova")

        # self.newFromTopspinAction = QAction(self)
        # self.newFromTopspinAction.setText("New from Topspin")

        # self.newAction.setIcon(QIcon(":file-new.svg"))
        # self.openAction = QAction(QIcon(":file-open.svg"), "&Open...", self)
        # self.saveAction = QAction(QIcon(":file-save.svg"), "&Save", self)
        self.openAction = QAction("&Open...", self)
        self.saveAction = QAction("&Save", self)
        self.exitAction = QAction("&Exit", self)
        # String-based key sequences
        # self.newAction.setShortcut("Ctrl+N")
        self.openAction.setShortcut("Ctrl+O")
        self.saveAction.setShortcut("Ctrl+S")
        # Help tips
        # newTip = "Create a new file"
        # self.newAction.setStatusTip(newTip)
        # self.newAction.setToolTip(newTip)
        # self.newAction.setWhatsThis("Create a new and empty text file")

        newMNovaTip = "Create a new problem from MNova"
        self.newFromMresNovaAction.setStatusTip(newMNovaTip)
        self.newFromMresNovaAction.setToolTip(newMNovaTip)
        self.newFromMresNovaAction.setWhatsThis("Create a new problem from MNova data")

        # newTopspinTip = "Create a new problem from topSpin"
        # self.newFromTopspinAction.setStatusTip(newTopspinTip)
        # self.newFromTopspinAction.setToolTip(newTopspinTip)
        # self.newFromTopspinAction.setWhatsThis("Create a new problem from Topspin data")

        # Edit actions
        # self.copyAction = QAction(QIcon(":edit-copy.svg"), "&Copy", self)
        # self.pasteAction = QAction(QIcon(":edit-paste.svg"), "&Paste", self)
        # self.cutAction = QAction(QIcon(":edit-cut.svg"), "C&ut", self)
        # self.copyAction = QAction("&Copy", self)
        # self.pasteAction = QAction("&Paste", self)
        # self.cutAction = QAction("C&ut", self)

        # Standard key sequence
        # self.copyAction.setShortcut(QKeySequence.Copy)
        # self.pasteAction.setShortcut(QKeySequence.Paste)
        # self.cutAction.setShortcut(QKeySequence.Cut)

        # Help actions
        self.helpContentAction = QAction("&Help Content...", self)
        self.aboutAction = QAction("&About...", self)

    # Uncomment this method to create a context menu using menu policies
    # def _createContextMenu(self):
    #     # Setting contextMenuPolicy
    #     self.centralWidget.setContextMenuPolicy(Qt.ActionsContextMenu)
    #     # Populating the widget with actions
    #     self.centralWidget.addAction(self.newAction)
    #     self.centralWidget.addAction(self.openAction)
    #     self.centralWidget.addAction(self.saveAction)
    #     self.centralWidget.addAction(selfnmrproblem.dbe.copyAction)
    #     self.centralWidget.addAction(self.pasteAction)
    #     self.centralWidget.addAction(self.cutAction)

    def contextMenuEvent(self, event):
        # Context menu
        menu = QMenu(self.centralWidget)
        # Populating the menu with actions
        # menu.addAction(self.newAction)
        menu.addAction(self.newFromMresNovaAction)
        # menu.addAction(self.newFromTopspinAction)
        menu.addAction(self.openAction)
        menu.addAction(self.saveAction)

        # Separator
        # separator = QAction(self)
        # separator.setSeparator(True)
        # menu.addAction(separator)
        # menu.addAction(self.copyAction)
        # menu.addAction(self.pasteAction)
        # menu.addAction(self.cutAction)

        # Launching the menu
        menu.exec(event.globalPos())

    def _connectActions(self):
        # Connect File actions
        # self.newAction.triggered.connect(self.newFile)
        self.newFromMresNovaAction.triggered.connect(self.newFromMresNova)
        # self.newFromTopspinAction.triggered.connect(self.newFromTopspin)
        # self.newAction.triggered.connect(self.newFile)
        self.openAction.triggered.connect(self.openFile)
        self.saveAction.triggered.connect(self.saveFile)
        self.exitAction.triggered.connect(self.close)

        # Connect Edit actions
        # self.copyAction.triggered.connect(self.copyContent)
        # self.pasteAction.triggered.connect(self.pasteContent)
        # self.cutAction.triggered.connect(self.cutContent)

        # Connect Help actions
        self.helpContentAction.triggered.connect(self.helpContent)
        self.aboutAction.triggered.connect(self.about)

        # Connect Open Recent to dynamically populate it
        # self.openRecentMenu.aboutToShow.connect(self.populateOpenRecent)

    # Slots
    def newFile(self):
        pass
        # Logic for creating a new file goes here...
        # self.csideWidget.setText("<b>File > New</b> clicked")

    def newFromMresNova(self):

        jsonfilepathlist = QFileDialog.getOpenFileName(
            self, "Select MResNova JSON file", "", "JSON (*_mresnova.json)"
        )

        if jsonfilepathlist[0] == "":
            return
        jsonfilepath = pathlib.Path(jsonfilepathlist[0])

        # check if jsoonfilepath exists
        if not jsonfilepath.exists():
            return

        # read in json file
        data = mnovautils.read_in_mesrenova_json(jsonfilepath)
        data = mnovautils.check_for_multiple_HSQC_expts(data)

        # check if json file contains a molecule ask user to add one if not
        smiles_str = ""
        if "smiles" not in data.keys():
            # read in excel file to get smiles string before over writing it
            excel_path = jsonfilepath.parent / (jsonfilepath.parent.name + ".xlsx")
            if excel_path.exists():
                try:
                    df = pd.read_excel(excel_path, sheet_name="molecule")
                    smiles_str = df["smiles"][0]
                except ValueError as ve:
                    smiles_str = ""
                    print("No sheet named molecule found in excel file")
                    warning_dialog("No sheet named molecule found in excel file", "No sheet named molecule found in excel file", self.nmrproblem.qtstarted)

            else:
                smiles_str = ""

            try:
                mol = Chem.MolFromSmiles("")
                smilesPNG = Draw.MolToImage(mol)
                print("smilesPNG:", smilesPNG)
            except Exception as e:
                print("Exception:", e)
                warning_dialog("Exception:", e, self.nmrproblem.qtstarted)
                smilesPNG = None

            dlg = SmilesDialog(smiles_str, smilesPNG)

            # dlg.setWindowTitle("HELLO!")

            # table_widget = MyTabWidget(dlg, self.nmrproblem)
            if dlg.exec():
                smile_str = dlg.get_smiles_string()
                if Chem.MolFromSmiles(smile_str):
                    warning_dialog(
                        "Invalid smiles string",
                        "Invalid smiles string",
                        self.nmrproblem.qtstarted,
                    )
                    return
                else:
                    data["smiles"] = smile_str
                # print(type(dlg.table_widget))
            else:
                print("Cancel!")
                return

        dataframes = mnovautils.create_dataframes_from_mresnova_json(data)

        # save dataframes to excel file

        # create a new excel file based on the json file parent folder name
        excel_fn = jsonfilepath.parent.name + ".xlsx"

        # create a a pathlib Path object for the new excel file and the parent folder
        folderpath = jsonfilepath.parent
        excel_path = folderpath / excel_fn

        if excel_path.exists():
            if not mnovautils.make_excel_backup(excel_path):
                return
        if not mnovautils.write_excel_file_from_mresnova_df_dict(
            dataframes, excel_path, True
        ):
            return

        if excel_path.parent.exists():
            data_info = nmrProblem.parse_argv(
                [excel_fn, str(folderpath), excel_fn], True
            )
            xy3_dlg = XY3dialog(
                java_available=java.JAVA_AVAILABLE,
                sheets_missing=nmrProblem.get_missing_sheets(
                    data_info["excel_fn"], True
                ),
            )
            if xy3_dlg.exec():
                xy3_calc_method, expts_available = xy3_dlg.get_method()
                # add "molecule" to expts_available
                expts_available.append("molecule")

                self.nmrproblem = nmrProblem.NMRproblem(
                    # nmrproblem = TestExcelSimpleNMR(
                    data_info,
                    java_available=java.JAVA_AVAILABLE,
                    xy3_calc_method=xy3_calc_method,
                    java_command=java.JAVA_COMMAND,
                    expts_available=expts_available,
                )
            else:
                return

            if self.nmrproblem.data_complete:
                print("line 1180")
                print("nmrproblem.h1")
                print(nmrproblem.h1)
                print("nmrproblem.hsqc")
                print(nmrproblem.hsqc)
                self.nmrproblem.initiate_df_data_complete()
                self.initiate_windows(self.nmrproblem)

    def newFromTopspin(self):
        print("New Topspin")

        # Logic for creating a new file goes here...
        # self.csideWidget.setText("<b>File > New</b> clicked")

    def openFile(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(
            self,
            "QFileDialog.getOpenFileName()",
            "",
            "All Files (*);;Python Files (*.py)",
            options=options,
        )
        if fileName:
            workingdir, fn = os.path.split(fileName)
            data_info = nmrProblem.parse_argv([fn, workingdir, fn], True)
            xy3_dlg = XY3dialog(
                java_available=java.JAVA_AVAILABLE,
                sheets_missing=nmrProblem.get_missing_sheets(
                    data_info["excel_fn"], True
                ),
            )
            if xy3_dlg.exec():
                xy3_calc_method, expts_available = xy3_dlg.get_method()
                # add "molecule" to expts_available
                expts_available.append("molecule")

                self.nmrproblem = nmrProblem.NMRproblem(
                    # nmrproblem = TestExcelSimpleNMR(
                    data_info,
                    java_available=java.JAVA_AVAILABLE,
                    xy3_calc_method=xy3_calc_method,
                    java_command=java.JAVA_COMMAND,
                    expts_available=expts_available,
                )
            else:
                return

            if self.nmrproblem.data_complete:

                self.nmrproblem.define_hsqc_f2integral()
                self.nmrproblem.initiate_df_data_complete()

                print("nmrproblem.c13")
                print(self.nmrproblem.c13[['ppm', 'numProtons', 'f1_i', 'f2_i', 'f1C_i', 'f2C_i']])
                print(self.nmrproblem.c13.columns)

                print("line 1234")
                print("nmrproblem.h1")
                print(self.nmrproblem.h1)
                print("nmrproblem.hsqc")
                print(self.nmrproblem.hsqc)
                self.initiate_windows(self.nmrproblem)
        # Logic for opening an existing file goes here...
        # self.centralWidget.setText("<b>File > Open...</b> clicked")

    def saveFile(self):
        for n in self.nmrproblem.nx_graph_molecule.nodes:
            self.nmrproblem.xy3[n] = self.moleculePlot.hmbc_graph_plots[n][
                "hmbc_nodes"
            ].get_offsets()[0]

            for k, v in self.nmrproblem.xy3.items():
                if isinstance(v, np.ndarray):
                    self.nmrproblem.xy3[k] = v.tolist()

        # taken outside of the loop to avoid overwriting xy3
        with open(
            os.path.join(self.nmrproblem.problemDirectoryPath, "xy3.json"), "w"
        ) as fp:
            json.dump(self.nmrproblem.xy3, fp, indent=2)

        # update C13 info based on xy3

        # using coordinates of xy3, match them to the coordinates in expected_molecule.molprops_df
        for k, [x, y] in self.nmrproblem.xy3.items():
            idx_c13 = int(k[1:])
            # find the closest row in expected_molecule.molprops_df based on x,y
            idx = self.nmrproblem.expected_molecule.molprops_df.apply(
                lambda row: math.sqrt((row["x"] - x) ** 2 + (row["y"] - y) ** 2), axis=1
            ).idxmin()

            self.nmrproblem.c13.loc[idx_c13, "C"] = idx
            self.nmrproblem.c13.loc[
                idx_c13, "predicted"
            ] = self.nmrproblem.expected_molecule.molprops_df.loc[idx, "ppm"]
            self.nmrproblem.c13.loc[
                idx_c13, "x"
            ] = self.nmrproblem.expected_molecule.molprops_df.loc[idx, "x"]
            self.nmrproblem.c13.loc[
                idx_c13, "y"
            ] = self.nmrproblem.expected_molecule.molprops_df.loc[idx, "y"]

        # create a html results file
        tohtml = problemtohtml.ProblemToHTML(self.nmrproblem)
        tohtml.write_html_report()

        # Logic for saving a file goes here...
        # self.centralWidget.setText("<b>File > Save</b> clicked")

    # def copyContent(self):
    #     pass
    #     # Logic for copying content goes here...
    #     # self.centralWidget.setText("<b>Edit > Copy</b> clicked")

    # def pasteContent(self):
    #     pass
    #     # Logic for pasting content goes here...
    #     # self.centralWidget.setText("<b>Edit > Paste</b> clicked")

    # def cutContent(self):
    #     pass
    #     # Logic for cutting content goes here...
    #     # self.centralWidget.setText("<b>Edit > Cut</b> clicked")

    def helpContent(self):
        webbrowser.open("https://github.com/EricHughesABC/simpleNMR")
        # webbrowser.open("C:\\Users\\vsmw51\\OneDrive - Durham University\\projects\\programming\\2022\\brukerWebinars\\docs\\source\\brukerWebinars\\test_template_rendered.html")
        # Logic for launching help goes here...
        # self.centralWidget.setText("<b>Help > Help Content...</b> clicked")

    def about(self):
        dlg = Aboutdialog(self)
        dlg.exec()
        # Logic for showing an about dialog content goes here...
        # self.centralWidget.setText("<b>Help > About...</b> clicked")

    # def populateOpenRecent(self):
    #     # Step 1. Remove the old options from the menu
    #     self.openRecentMenu.clear()
    #     # Step 2. Dynamically create the actions
    #     actions = []
    #     filenames = [f"File-{n}" for n in range(5)]
    #     for filename in filenames:
    #         action = QAction(filename, self)
    #         action.triggered.connect(partial(self.openRecentFile, filename))
    #         actions.append(action)
    #     # Step 3. Add the actions to the menu
    #     self.openRecentMenu.addActions(actions)

    # def openRecentFile(self, filename):
    #     pass
    #     # Logic for opening a recent file goes here...
    #     # self.centralWidget.setText(f"<b>{filename}</b> opened")

    # def getWordCount(self):
    #     # Logic for computing the word count goes here...
    #     return 42

    #     # self.show()

    def axes_leave_callback(self, event):

        if self.selection == None:
            return
        else:

            if "molecule" in self.selection.artist.axes.get_gid():
                lbl = "C" + str(int(self.selection.index) + 1)

                self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_nodes"].set_visible(False)
                self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_edges"].set_visible(False)
                for k, v in self.moleculePlot.hmbc_graph_plots[lbl][
                    "hmbc_labels"
                ].items():
                    v.set_visible(False)

                self.moleculePlot.ax.figure.canvas.draw_idle()
                self.spectraPlot.c13spec_ax.figure.canvas.draw_idle()
                self.spectraPlot.h1spec_ax.figure.canvas.draw_idle()

    # def loadPage(self):

    #     url = QUrl.fromUserInput("https://jsme-editor.github.io/dist/JSME_minimal.html")

    #     if url.isValid():
    #         self.webEngineView.load(url)

    # def warning_dialog(self, return_message, title_message, qstarted):
    #     if qstarted:
    #         # create qt5 warning dialog with return_message
    #         msg_box = QMessageBox()
    #         msg_box.setIcon(QMessageBox.Warning)
    #         msg_box.setText(return_message)
    #         msg_box.setStandardButtons(QMessageBox.Ok)
    #         msg_box.exec()
    #         msg_box.setWindowTitle(title_message)
    #     else:
    #         print(return_message)

    @pyqtSlot()
    def on_sync_click(self):
        mnova_dir = pathlib.Path(self.nmrproblem.problemDirectoryPath)
        mnova_fn = mnova_dir / (mnova_dir.name + "_mresnova.json")

        if not mnova_fn.is_file():
            warning_dialog(
                f"{mnova_fn} does not exist",
                "MNOVA json file does not exist",
                self.nmrproblem.qtstarted,
            )
            return

        # read in json file
        jsonfilepath = pathlib.Path(mnova_fn)

        # read in json file
        data = mnovautils.read_in_mesrenova_json(jsonfilepath)
        data = mnovautils.check_for_multiple_HSQC_expts(data)

        # check that the json file contains a smiles string
        if "smiles" not in data.keys():
            warning_dialog(
                "smiles not in json file",
                "No smiles in json file",
                self.nmrproblem.qtstarted,
            )
            # allow user to add smiles string to dataframe or  use old smiles string from nmrproblem if present
            # open smiles dialog
            if self.nmrproblem.smiles != "":
                smilesStr = self.nmrproblem.smiles
            else:
                smilesStr = ""

            mol = Chem.MolFromSmiles(smilesStr)
            smilesPNG = Draw.MolToImage(mol)
            dlg = SmilesDialog(smilesStr, smilesPNG)

            # dlg.setWindowTitle("HELLO!")

            # table_widget = MyTabWidget(dlg, self.nmrproblem)
            if dlg.exec():
                smilesStr = dlg.get_smiles_string()
                # check if smiles string is valid
                mol = Chem.MolFromSmiles(smilesStr)
                if mol is None:
                    warning_dialog(
                        "smiles string is invalid",
                        "smiles string is invalid",
                        self.nmrproblem.qtstarted,
                    )
                    return
                else:
                    # add smiles string to molecule dataframe
                    data["smiles"] = smilesStr
            else:
                print("Cancel!")
                return

        # check that the smiles string is the same as the one in the nmrproblem
        if data["smiles"] != self.nmrproblem.smiles:
            warning_dialog(
                "smiles string mismatch with json file",
                "smiles string mismatch with json file",
                self.nmrproblem.qtstarted,
            )
            return

        # check that HSQC data is in the json file
        if "HSQC" not in data.keys():
            warning_dialog(
                "HSQC data not in json file",
                "HSQC data not in json file",
                self.nmrproblem.qtstarted,
            )
            return

        dataframes = mnovautils.create_dataframes_from_mresnova_json(data)

        # check that the HSQC dataframe is not empty
        if dataframes["HSQC"].empty:
            warning_dialog(
                "HSQC dataframe is empty",
                "HSQC dataframe is empty",
                self.nmrproblem.qtstarted,
            )
            return

        # check that the molecule dataframe is not empty
        if dataframes["molecule"].empty:
            warning_dialog(
                "molecule dataframe is empty",
                "molecule dataframe is empty",
                self.nmrproblem.qtstarted,
            )
            # allow user to add smiles string to dataframe or  or use old smiles string from nmrproblem if present
            # open smiles dialog
            if self.nmrproblem.smiles != "":
                smilesStr = self.nmrproblem.smiles
            else:
                smilesStr = ""

            mol = Chem.MolFromSmiles(smilesStr)
            smilesPNG = Draw.MolToImage(mol)
            dlg = SmilesDialog(smilesStr, smilesPNG)

            # dlg.setWindowTitle("HELLO!")

            # table_widget = MyTabWidget(dlg, self.nmrproblem)
            if dlg.exec():
                smilesStr = dlg.get_smiles_string()
                # check if smiles string is valid
                mol = Chem.MolFromSmiles(smilesStr)
                if mol is None:
                    warning_dialog(
                        "smiles string is invalid",
                        "smiles string is invalid",
                        self.nmrproblem.qtstarted,
                    )
                    return
                else:
                    # add smiles string to molecule dataframe
                    data["smiles"] = smilesStr
                # print(type(dlg.table_widget))
            else:
                print("Cancel!")
                return

        # save dataframes to excel file

        # create a new excel file based on the json file parent folder name
        excel_fn = jsonfilepath.parent.name + ".xlsx"

        # create a a pathlib Path object for the new excel file and the parent folder
        folderpath = jsonfilepath.parent
        excel_path = folderpath / excel_fn
        if excel_path.exists():
            if not mnovautils.make_excel_backup(excel_path):
                return
        # save the dataframes to the excel file
        if not mnovautils.write_excel_file_from_mresnova_df_dict(
            dataframes, excel_path
        ):
            warning_dialog(
                "could not save excel file",
                "could not save excel file",
                self.nmrproblem.qtstarted,
            )
            return
        # with pd.ExcelWriter(excel_path) as writer:
        #     for k, v in dataframes.items():
        #         v.to_excel(writer, sheet_name=k)

        # sync nmrproblem from newly saved excel file
        # rtn_ok, msg = self.nmrproblem.sync_class_from_mnova_excel(str(folderpath))
        rtn_ok, msg = self.nmrproblem.init_class_from_excel(str(folderpath))
        if rtn_ok:
            self.nmrproblem.data_complete = True
        else:
            self.nmrproblem.data_complete = False
            return

        self.nmrproblem.attempt_to_assign_integrals_to_nmrproblem()

        if self.nmrproblem.data_complete:

            # set xlims of 13C and 1H 1D spectra to +/- 10% of biggest and smallest ppm
            self.nmrproblem.min_max_1D_ppm = []

            h1 = self.nmrproblem.h1

            ppm_min = h1.ppm.min() - (h1.ppm.max() - h1.ppm.min()) * 0.3
            ppm_max = h1.ppm.max() + (h1.ppm.max() - h1.ppm.min()) * 0.3

            self.nmrproblem.min_max_1D_ppm.append((ppm_max, ppm_min))
            c13 = self.nmrproblem.c13

            ppm_min = c13.ppm.min() - (c13.ppm.max() - c13.ppm.min()) * 0.30
            ppm_max = c13.ppm.max() + (c13.ppm.max() - c13.ppm.min()) * 0.30

            self.nmrproblem.min_max_1D_ppm.append((ppm_max, ppm_min))

            self.nmrproblem.initiate_df_data_complete()

            print("nmrproblem.c13")
            print(self.nmrproblem.c13)
            print("line 1561")
            print("nmrproblem.h1")
            print(self.nmrproblem.h1)
            print("nmrproblem.hsqc")
            print(self.nmrproblem.hsqc)
            self.initiate_windows(self.nmrproblem)

    def highlight_hmbc_peaks(self, lbl, sel):
        if lbl[0] == "H":
            if lbl in self.nmrproblem.hsqcH1labelC13label.keys():
                lbl = self.nmrproblem.hsqcH1labelC13label[lbl]

        if lbl in self.nmrproblem.hmbc_graph_edges.keys():
            for i, ci in enumerate(self.nmrproblem.hmbc_graph_edges[lbl]):
                sel.extras.append(
                    self.cursor.add_highlight(self.spectraPlot.peak_overlays_dict[ci])
                )
                sel.extras[-1].set_linewidth(1)
                sel.extras[-1].set_color(self.hmbc_edge_colors[i])

            # highlight corresponding HMBC connected peaks in H1 spectrum
            for i, ci in enumerate(self.nmrproblem.hmbc_graph_edges[lbl]):
                hmbc_h1s = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f1C_i == ci][
                    "f2H_i"
                ].tolist()
                for j, hi in enumerate(hmbc_h1s):
                    sel.extras.append(
                        self.cursor.add_highlight(
                            self.spectraPlot.peak_overlays_dict[hi]
                        )
                    )
                    sel.extras[-1].set_linewidth(1)
                    sel.extras[-1].set_color(self.hmbc_edge_colors[i])

    def redraw_axes(self):
        self.spectraPlot.c13spec_ax.figure.canvas.draw_idle()
        self.spectraPlot.h1spec_ax.figure.canvas.draw_idle()
        self.spectraPlot.c13dist_ax.figure.canvas.draw_idle()
        self.spectraPlot.h1dist_ax.figure.canvas.draw_idle()
        self.moleculePlot.ax.figure.canvas.draw_idle()

    def draw_hmbc_graph_network(self, lbl):
        if lbl in self.nmrproblem.hmbc_graphs.keys():
            self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_nodes"].set_visible(True)
            if not isinstance(
                self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_edges"], list
            ):
                self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_edges"].set_visible(True)
            for k, v in self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_labels"].items():
                v.set_visible(True)
            self.old_node = lbl

    def hide_hmbc_graph_networks(self, lbls=None):

        if isinstance(lbls, str):
            if lbls[0] == "H":
                return
            lbls_list = [lbls]
        else:
            lbls_list = self.nmrproblem.hmbc_graphs.keys()

        for lbl in lbls_list:
            self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_nodes"].set_visible(False)
            if not isinstance(
                self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_edges"], list
            ):
                self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_edges"].set_visible(False)
            for v in self.moleculePlot.hmbc_graph_plots[lbl]["hmbc_labels"].values():
                v.set_visible(False)

        self.old_node = False


if __name__ == "__main__":

    # set current working directory to the directory of this file
    os.chdir(os.path.dirname(os.path.realpath(__file__)))

    app = QApplication(sys.argv)

    xy3_calc_method_str = "xy3"

    # check if  command line arguments are given
    if len(sys.argv) > 1:

        data_info = nmrProblem.parse_argv()
        xy3_dlg = XY3dialog(
            java_available=java.JAVA_AVAILABLE,
            sheets_missing=nmrProblem.get_missing_sheets(data_info["excel_fn"], True),
        )

        if xy3_dlg.exec_():
            xy3_calc_method_str, expts_available = xy3_dlg.get_method()
            # add "molecule" to the list of available experiments
            expts_available.append("molecule")

            # data_info = nmrProblem.parse_argv()
            nmrproblem = nmrProblem.NMRproblem(
                # nmrproblem = TestExcelSimpleNMR(
                data_info,
                java_available=java.JAVA_AVAILABLE,
                xy3_calc_method=xy3_calc_method_str,
                java_command=java.JAVA_COMMAND,
                expts_available=expts_available,
            )

        else:
            # start the GUI with empty data
            data_info = nmrProblem.parse_argv(my_argv=[""])
            nmrproblem = nmrProblem.NMRproblem(
                # nmrproblem = TestExcelSimpleNMR(
                data_info,
                java_available=java.JAVA_AVAILABLE,
                xy3_calc_method=xy3_calc_method_str,
                java_command=java.JAVA_COMMAND,
                expts_available=[],
            )

    else:
        data_info = nmrProblem.parse_argv()
        nmrproblem = nmrProblem.NMRproblem(
            # nmrproblem = TestExcelSimpleNMR(
            data_info,
            java_available=java.JAVA_AVAILABLE,
            xy3_calc_method=xy3_calc_method_str,
            java_command=java.JAVA_COMMAND,
            expts_available=[],
        )

    if nmrproblem.data_complete:

        print("line 1692")
        print("nmrproblem.h1")
        print(nmrproblem.h1)
        print("nmrproblem.hsqc")
        print(nmrproblem.hsqc)
        if (nmrproblem.c13_from_hsqc) or (nmrproblem.h1_df_integral_added):
            # nmrproblem.df.loc["integral", nmrproblem.carbonAtoms] = nmrproblem.c13.attached_protons.values
            nmrproblem.df.loc[
                "attached protons", nmrproblem.carbonAtoms
            ] = nmrproblem.c13.attached_protons.values
            nmrproblem.df.loc[
                "C13 hyb", nmrproblem.carbonAtoms
            ] = nmrproblem.c13.attached_protons.values
            # do the same for nmrproblem,protonAtoms
            nmrproblem.df.loc[
                "integral", nmrproblem.protonAtoms
            ] = nmrproblem.h1.integral.values
            nmrproblem.df.loc[
                "attached protons", nmrproblem.protonAtoms
            ] = nmrproblem.h1.integral.values
            nmrproblem.df.loc[
                "C13 hyb", nmrproblem.protonAtoms
            ] = nmrproblem.h1.integral.values
        else:
            nmrproblem.define_hsqc_f2integral()
            nmrproblem.define_c13_attached_protons()

        nmrProblem.build_model(nmrproblem)
        nmrProblem.build_molecule_graph_network(nmrproblem)
        # nmrProblem.build_xy3_representation_of_molecule(nmrproblem)
        nmrproblem.build_xy3()
        # save dataframes to json
        nmrproblem.save_dataframes_in_nmrproblem_to_json()

    ex = MainWidget(nmrproblem)

    ex.show()
    sys.exit(app.exec_())
