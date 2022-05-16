import sys
from PyQt5.QtWidgets import (QWidget,
                             QMainWindow,
                             QApplication, 
                             QVBoxLayout, 
                             QHBoxLayout, 
                             QLineEdit, 
                             QPushButton, 
                             QMessageBox, 
                             QSplitter,
                             QMenu,
                             QToolBar,
                             QSpinBox,
                             QLabel,
                             QAction)

from PyQt5.QtGui import QIcon, QKeySequence
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtCore import Qt

from functools import partial

# from PyQt5.Qt import Vertical

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import mplcursors

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import networkx as nx
import nx_pylab
import numpy as np
import pandas as pd

import PIL
from PIL import Image

import nmrProblem
import qt5_tabs_001
from moleculePlot import MatplotlibMoleculePlot
from spectraPlot import MatplotlibH1C13Plot

from numbers import Number

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
    hmbc_edge_colors = ["b", "g", "c", "y", "m"]
    for i, c in enumerate(catoms):

        if c not in hmbc_edges.keys():
            continue

        # create hmbc graph for node c and add  xy coodinates
        hmbc_graphs[c] = {}
        hmbc_graphs[c]["graph"] = nx.Graph()
        hmbc_graphs[c]["xy"] = dict(
            (k, xy3[k]) for k in [c] + list(hmbc_edges[c])
        )
        hmbc_graphs[c]["colors"] = []

        # add nodes to hmbc graph
        hmbc_graphs[c]["graph"].add_nodes_from([c] + list(hmbc_edges[c]))

        # add edges
        for i, c1 in enumerate(hmbc_edges[c]):
            hmbc_graphs[c]["graph"].add_edge(c, c1)
            hmbc_graphs[c]["colors"].append(hmbc_edge_colors[i])

    return hmbc_graphs




class MoleculePlotCanvas(FigureCanvasQTAgg):

    def __init__(self, fig, parent=None):
        # self.molecule_fig = Figure(figsize=(width, height), dpi=dpi)
        super(MoleculePlotCanvas, self).__init__(fig)
        # gs = GridSpec(1, 1, figure=self.fig)
        # self.ax = self.fig.add_subplot(gs[:, :])
        # self.molecule_ax = self.molecule_fig.add_subplot()
        # self.molecule_ax.plot([1,2,3], [1,2,3])

class MainWidget(QMainWindow):


          

    def __init__(self,  nmrproblem, parent=None):
        self.nmrproblem = nmrproblem
        self.hmbc_edge_colors = ("b", "g", "c", "y", "m")
        self.old_node = ""
        self.new_node = ""
        self.selection = None

        super().__init__(parent)

        self.setGeometry(100, 100, 1200, 900)
        self.setWindowTitle('SimpleNMR')

        self.centralWidget = QWidget()
        # self.centralWidget.setAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
        self.setCentralWidget(self.centralWidget)

        self._createActions()
        self._createMenuBar()
        # self._createToolBars()

        # Uncomment the call to ._createContextMenu() below to create a context
        # menu using menu policies. To test this out, you also need to
        # comment .contextMenuEvent() and uncomment ._createContextMenu()

        # self._createContextMenu()

        self._connectActions()
        self._createStatusBar()

        if self.nmrproblem.data_complete:


            jsmewidget = QWidget()
            moleculewidget = QWidget()
            spectrawidget = QWidget()

            self.moleculePlot = MatplotlibMoleculePlot(self.nmrproblem)
            self.spectraPlot = MatplotlibH1C13Plot(self.nmrproblem)

            self.spectraCanvas =  MoleculePlotCanvas(self.spectraPlot)
            self.moleculeCanvas =  MoleculePlotCanvas(self.moleculePlot)
            self.moleculeCanvas.setGeometry(0,0,400,400)

            moltoolbar = NavigationToolbar(self.moleculeCanvas, self)
            spctoolbar = NavigationToolbar(self.spectraCanvas, self)


            self.webEngineView = QWebEngineView()
            self.loadPage()

            self.smilesInput = QLineEdit()
            self.button = QPushButton('Display Smiles Molecule', self)

            # connect button to function on_click
            self.button.clicked.connect(self.on_click)


            splitter1 = QSplitter(Qt.Vertical)
            splitter2 = QSplitter(Qt.Horizontal)
        

            hbox = QHBoxLayout()
            molvbox = QVBoxLayout()
            jsmevbox = QVBoxLayout()
            spcvbox = QVBoxLayout()

            spcvbox.addWidget(spctoolbar)
            spcvbox.addWidget(self.spectraCanvas)

            spectrawidget.setLayout(spcvbox)

            molvbox.addWidget(moltoolbar)
            molvbox.addWidget(self.moleculeCanvas)

            moleculewidget.setLayout(molvbox)
            
            jsmevbox.addWidget(self.webEngineView)
            jsmevbox.addWidget(self.smilesInput)
            jsmevbox.addWidget(self.button)

            jsmewidget.setLayout(jsmevbox)

            splitter1.addWidget(moleculewidget)
            splitter1.addWidget(jsmewidget)

            splitter2.addWidget(spectrawidget)
            splitter2.addWidget(splitter1)



            hbox.addWidget(splitter2)
            # hbox.addLayout(molvbox)
            self.centralWidget.setLayout(hbox)

            # self.moleculeCanvas.mpl_connect("axes_leave_event", 
            #                                          self.axes_leave_callback)
            # self.spectraCanvas.mpl_connect("axes_leave_event", 
            #                                          self.axes_leave_callback)

            # hover=mplcursors.HoverMode.Transient
            self.cursor = mplcursors.cursor( [self.moleculePlot.mol_nodes] + 
                                        self.spectraPlot.peak_overlays[0] + 
                                        self.spectraPlot.peak_overlays[1], 
                                        hover=mplcursors.HoverMode.Transient, 
                                        highlight=True)


            self.setup_mplcursors()

    def _createMenuBar(self):
        menuBar = self.menuBar()
        # File menu
        fileMenu = QMenu("&File", self)
        menuBar.addMenu(fileMenu)
        fileMenu.addAction(self.newAction)
        fileMenu.addAction(self.openAction)
        # Open Recent submenu
        self.openRecentMenu = fileMenu.addMenu("Open Recent")
        fileMenu.addAction(self.saveAction)
        # Separator
        fileMenu.addSeparator()
        fileMenu.addAction(self.exitAction)
        # Edit menu
        editMenu = menuBar.addMenu("&Edit")
        editMenu.addAction(self.copyAction)
        editMenu.addAction(self.pasteAction)
        editMenu.addAction(self.cutAction)
        # Separator
        editMenu.addSeparator()
        # Find and Replace submenu
        findMenu = editMenu.addMenu("Find and Replace")
        findMenu.addAction("Find...")
        findMenu.addAction("Replace...")
        # Help menu
        helpMenu = menuBar.addMenu(QIcon(":help-content.svg"), "&Help")
        helpMenu.addAction(self.helpContentAction)
        helpMenu.addAction(self.aboutAction)

    def _createToolBars(self):
        # File toolbar
        fileToolBar = self.addToolBar("File")
        fileToolBar.setMovable(False)
        fileToolBar.addAction(self.newAction)
        fileToolBar.addAction(self.openAction)
        fileToolBar.addAction(self.saveAction)
        # Edit toolbar
        editToolBar = QToolBar("Edit", self)
        self.addToolBar(editToolBar)
        editToolBar.addAction(self.copyAction)
        editToolBar.addAction(self.pasteAction)
        editToolBar.addAction(self.cutAction)
        # Widgets
        self.fontSizeSpinBox = QSpinBox()
        self.fontSizeSpinBox.setFocusPolicy(Qt.NoFocus)
        editToolBar.addWidget(self.fontSizeSpinBox)

    def _createStatusBar(self):
        self.statusbar = self.statusBar()
        # Temporary message
        self.statusbar.showMessage("Ready", 3000)
        # Permanent widget
        self.wcLabel = QLabel(f"{self.getWordCount()} Words")
        self.statusbar.addPermanentWidget(self.wcLabel)

    def _createActions(self):
        # File actions
        self.newAction = QAction(self)
        self.newAction.setText("&New")
        self.newAction.setIcon(QIcon(":file-new.svg"))
        self.openAction = QAction(QIcon(":file-open.svg"), "&Open...", self)
        self.saveAction = QAction(QIcon(":file-save.svg"), "&Save", self)
        self.exitAction = QAction("&Exit", self)
        # String-based key sequences
        self.newAction.setShortcut("Ctrl+N")
        self.openAction.setShortcut("Ctrl+O")
        self.saveAction.setShortcut("Ctrl+S")
        # Help tips
        newTip = "Create a new file"
        self.newAction.setStatusTip(newTip)
        self.newAction.setToolTip(newTip)
        self.newAction.setWhatsThis("Create a new and empty text file")
        # Edit actions
        self.copyAction = QAction(QIcon(":edit-copy.svg"), "&Copy", self)
        self.pasteAction = QAction(QIcon(":edit-paste.svg"), "&Paste", self)
        self.cutAction = QAction(QIcon(":edit-cut.svg"), "C&ut", self)
        # Standard key sequence
        self.copyAction.setShortcut(QKeySequence.Copy)
        self.pasteAction.setShortcut(QKeySequence.Paste)
        self.cutAction.setShortcut(QKeySequence.Cut)
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
    #     self.centralWidget.addAction(self.copyAction)
    #     self.centralWidget.addAction(self.pasteAction)
    #     self.centralWidget.addAction(self.cutAction)

    def contextMenuEvent(self, event):
        # Context menu
        menu = QMenu(self.centralWidget)
        # Populating the menu with actions
        menu.addAction(self.newAction)
        menu.addAction(self.openAction)
        menu.addAction(self.saveAction)
        # Separator
        separator = QAction(self)
        separator.setSeparator(True)
        menu.addAction(separator)
        menu.addAction(self.copyAction)
        menu.addAction(self.pasteAction)
        menu.addAction(self.cutAction)
        # Launching the menu
        menu.exec(event.globalPos())

    def _connectActions(self):
        # Connect File actions
        self.newAction.triggered.connect(self.newFile)
        self.openAction.triggered.connect(self.openFile)
        self.saveAction.triggered.connect(self.saveFile)
        self.exitAction.triggered.connect(self.close)
        # Connect Edit actions
        self.copyAction.triggered.connect(self.copyContent)
        self.pasteAction.triggered.connect(self.pasteContent)
        self.cutAction.triggered.connect(self.cutContent)
        # Connect Help actions
        self.helpContentAction.triggered.connect(self.helpContent)
        self.aboutAction.triggered.connect(self.about)
        # Connect Open Recent to dynamically populate it
        self.openRecentMenu.aboutToShow.connect(self.populateOpenRecent)

    # Slots
    def newFile(self):
        pass
        # Logic for creating a new file goes here...
        # self.csideWidget.setText("<b>File > New</b> clicked")

    def openFile(self):
        pass
        # Logic for opening an existing file goes here...
        # self.centralWidget.setText("<b>File > Open...</b> clicked")

    def saveFile(self):
        pass
        # Logic for saving a file goes here...
        # self.centralWidget.setText("<b>File > Save</b> clicked")

    def copyContent(self):
        pass
        # Logic for copying content goes here...
        # self.centralWidget.setText("<b>Edit > Copy</b> clicked")

    def pasteContent(self):
        pass
        # Logic for pasting content goes here...
        # self.centralWidget.setText("<b>Edit > Paste</b> clicked")

    def cutContent(self):
        pass
        # Logic for cutting content goes here...
        # self.centralWidget.setText("<b>Edit > Cut</b> clicked")

    def helpContent(self):
        pass
        # Logic for launching help goes here...
        # self.centralWidget.setText("<b>Help > Help Content...</b> clicked")

    def about(self):
        pass
        # Logic for showing an about dialog content goes here...
        # self.centralWidget.setText("<b>Help > About...</b> clicked")

    def populateOpenRecent(self):
        # Step 1. Remove the old options from the menu
        self.openRecentMenu.clear()
        # Step 2. Dynamically create the actions
        actions = []
        filenames = [f"File-{n}" for n in range(5)]
        for filename in filenames:
            action = QAction(filename, self)
            action.triggered.connect(partial(self.openRecentFile, filename))
            actions.append(action)
        # Step 3. Add the actions to the menu
        self.openRecentMenu.addActions(actions)

    def openRecentFile(self, filename):
        pass
        # Logic for opening a recent file goes here...
        # self.centralWidget.setText(f"<b>{filename}</b> opened")

    def getWordCount(self):
        # Logic for computing the word count goes here...
        return 42

        
        # self.show()


    def axes_leave_callback(self, event):
        print("axes_leave_event", event) 


        if self.selection == None:
            print("return None")
            return
        else:
            print("self.selection.artist.axes.get_gid()\n\t", 
                    self.selection.artist.axes.get_gid())

            print("self.selection.index\n\t", self.selection.index)
            

            if "molecule" in self.selection.artist.axes.get_gid():
                lbl = 'C' + str(int(self.selection.index) + 1)
                print("lbl", lbl)

                self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_nodes'].set_visible(False)
                self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_edges'].set_visible(False)
                for k, v  in self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_labels'].items():                                    
                    v.set_visible(False)

                self.moleculePlot.ax.figure.canvas.draw_idle()
                self.spectraPlot.c13spec_ax.figure.canvas.draw_idle()
                self.spectraPlot.h1spec_ax.figure.canvas.draw_idle()
            else:
                print("sel.index\n\t", self.selection.index)

    def loadPage(self):

        url = QUrl.fromUserInput("https://jsme-editor.github.io/dist/JSME_minimal.html")

        if url.isValid():
            self.webEngineView.load(url)


    def setup_mplcursors(self):

        @self.cursor.connect("remove")
        def on_remove(sel):
            if "molecule" in sel.artist.axes.get_gid():
                lbl = 'C' + str(int(sel.index) + 1)
            else:
                lbl = sel.artist.get_label()

            self.hide_hmbc_graph_networks(lbl)
            self.redraw_axes()
            self.selection = None

        @self.cursor.connect("add")
        def on_add(sel):
            self.selection = sel

            udic = self.nmrproblem.udic
            catoms = self.nmrproblem.carbonAtoms
            hatoms = self.nmrproblem.protonAtoms

            lbl = 'C' + str(int(sel.index) + 1)
            print("lbl", lbl)
            # print(dir(sel))
            # print("sel.artist.axes.get_label()", sel.artist.axes.get_label())
            # print("sel.artist.axes.get_gid()", sel.artist.axes.get_gid())
            # print("sel.target\n", sel.target)

            if "molecule" in sel.artist.axes.get_gid():
                lbl = 'C' + str(int(sel.index) + 1)
                print("lbl molecule", lbl)




                sel.annotation.set_text(lbl)
                sel.extras[0].set_edgecolor("r")
                sel.extras[0].set_linewidth(1)
                sel.extras[0].set_facecolor("w")

                # set annotation text
                c13pmm_text, _ = udic[1]["labels1_dict"][lbl][1].split("\n")
                _, c13groups_text = udic[1]["labels1_dict"][lbl][0].split(":")
                sel.annotation.set_text(c13pmm_text + "\n\n" + c13groups_text)

                            # cursor.add_highlight(self.spectraPlot.peak_overlays_dict[ci])
                # sel.extras.append(cursor.add_highlight(pairs[sel.artist]))
                print("lbl", lbl)
                sel.extras.append(self.cursor.add_highlight(self.spectraPlot.peak_overlays_dict[lbl]))
                sel.extras[-1].set_visible(True)
                sel.extras[-1].set_linewidth(0.75)
                sel.extras[-1].set_color("red")

                highlighted_H1_pks = udic[1]["info"].loc[lbl, "hsqc"]

                # highlight corresponding H1 HSQC peaks ins 1D proton Spectrum
                for hpk in highlighted_H1_pks:
                    sel.extras.append(
                        self.cursor.add_highlight(self.spectraPlot.peak_overlays_dict[hpk])
                    )
                    sel.extras[-1].set_linewidth(0.75)
                    sel.extras[-1].set_color("red")

                # highlight corresponding HMBC connected peaks in C13 spectrum
                self.highlight_hmbc_peaks(lbl, sel)

                if lbl in self.nmrproblem.hmbc_graphs.keys():
                    self.hide_hmbc_graph_networks()
                    self.draw_hmbc_graph_network(lbl)
                    self.redraw_axes()

            else:
                lbl = sel.artist.get_label()
                print("lbl", lbl)

                self.new_node = lbl
                print("lbl =", lbl)
                if str(lbl) in hatoms:
                    ii = 0
                else:
                    ii = 1

                print("ii", ii)
                x, y = sel.target

                # use artist to labal to find out peak id, H1, H2, ... or C1, C2

                selected_pk = [str(lbl)]
                highlighted_pks = udic[ii]["info"].loc[selected_pk[0], "hsqc"]
                print("highlighted peaks", highlighted_pks)

                if "H" in selected_pk[0]:
                    H_pks = selected_pk
                    C_pks = highlighted_pks
                else:
                    H_pks = highlighted_pks
                    C_pks = selected_pk

                # set selected peak to red with a linewidth of 0.75
                sel.extras[0].set_linewidth(0.75)
                sel.extras[0].set_color("red")

                ii2 = 0
                if ii == 0:
                    ii2 = 1

                # change selection text depending on what hight peak was picked
                # top, middle, bottom
                # pk_y is an index position from 0,1,2
                pk_y = 2 - int(y * 100 / udic[ii]["info"].loc[lbl, "pk_y"]) // 35

                # set annotation text
                # sel.annotation.set_text(udic[ii]['labels1_dict'][selected_pk[0]][pk_y])
                sel.annotation.set_text(udic[ii]["labels1_dict"][lbl][pk_y])
                print(ii, lbl, pk_y, udic[ii]["labels1_dict"][lbl][pk_y])

                # highlight coresponding proton or carbon peaks
                # from hsqc data find corresponding peaks
                highlighted_pks = udic[ii]["info"].loc[selected_pk[0], "hsqc"]
                for hpk in highlighted_pks:
                    sel.extras.append(
                        self.cursor.add_highlight(self.spectraPlot.peak_overlays_dict[hpk])
                    )
                    sel.extras[-1].set_linewidth(0.75)
                    sel.extras[-1].set_color("red")

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
                    numplots = len(self.spectraPlot.h1c13distlist[0][pk]) - 1
                    for i, aa in enumerate(self.spectraPlot.h1c13distlist[0][pk]):
                        ppmH.append(self.nmrproblem.udic[0]["info"].loc[pk, "ppm"])
                        sel.extras.append(self.cursor.add_highlight(aa))
                        sel.extras[-1].set_visible(True)
                        sel.extras[-1].set_linewidth(0.75)
                        sel.extras[-1].set_color(colors[i % 7])
                        # do not add legend info if plot is just single line showing where peak pos is
                        if i < numplots:
                            plabelsH.append(aa.get_label())
                            plinesH.append(sel.extras[-1])

                # create carbon distribution legends
                for pk in C_pks:
                    numplots = len(self.spectraPlot.h1c13distlist[1][pk]) - 1
                    for i, aa in enumerate(self.spectraPlot.h1c13distlist[1][pk]):
                        ppmC.append(self.nmrproblem.udic[1]["info"].loc[pk, "ppm"])
                        sel.extras.append(self.cursor.add_highlight(aa))
                        sel.extras[-1].set_visible(True)
                        sel.extras[-1].set_linewidth(0.75)
                        sel.extras[-1].set_color(colors[i % 7])
                        if i < numplots:
                            plabelsC.append(aa.get_label())
                            plinesC.append(sel.extras[-1])

                # adjust x-axis width H1 and C13 of distribution plots
                # and add legend information
                if len(ppmH) > 0:
                    # calculate average position of peaks which will be used to adjust the x axis ppm range
                    ppmmmH = np.mean(ppmH)
                    self.spectraPlot.h1dist_ax.set_xlim(ppmmmH + 2, ppmmmH - 2)
                    self.spectraPlot.h1dist_ax.legend(plinesH, plabelsH)
                if len(ppmC) > 0:
                    ppmmmC = np.mean(ppmC)
                    self.spectraPlot.c13dist_ax.set_xlim(ppmmmC + 50, ppmmmC - 50)
                    self.spectraPlot.c13dist_ax.legend(plinesC, plabelsC)

                # update molecule plot to highlight HMBC graph fragment

                # highlight corresponding HMBC connected peaks in C13 spectrum
                self.highlight_hmbc_peaks(lbl, sel)

                # if peak selected from H1 spectrum obtained corresponding carbon label
                if lbl not in catoms:
                    print(lbl)
                    if lbl in self.nmrproblem.hsqcH1labelC13label.keys():
                        lbl = self.nmrproblem.hsqcH1labelC13label[lbl]

                self.hide_hmbc_graph_networks()
                self.draw_hmbc_graph_network(lbl)
                self.redraw_axes()






    @pyqtSlot()
    def on_click(self):
        smilesText = self.smilesInput.text()
        QMessageBox.question(self, 'Smiles', "Molecule: " + smilesText, QMessageBox.Ok, QMessageBox.Ok)
        # self.mplplot.molecule_ax.clear()
        # self.mplplot.molecule_ax.set_title(smilesText)
        
        # print(type(self.mplplot.ax))
        # ax = self.mplplot.pltlines.axes
        # ax.text(1, 1, smilesText, ha='left', rotation=15, wrap=True)

        molecule = Chem.MolFromSmiles(smilesText)
        Draw.MolToFile(molecule,'molecule.png')
        self.nmrproblem.png  = Image.open("molecule.png").transpose(PIL.Image.FLIP_LEFT_RIGHT)

        self.moleculePlot.bkgnd.set_data(self.nmrproblem.png)
        # self.mplplot.ax.imshow(molecule_image, alpha=0.5)  

        # self.mplplot.ax.set_xlim(350,-50)
        # self.mplplot.ax.set_ylim(350, -50)  

        # self.mplplot.draw()

        # self.mplplot.update_molecule_plot()



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
                hmbc_h1s = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f1C_i == ci]["f2H_i"].tolist()
                for j, hi in enumerate(hmbc_h1s):
                    sel.extras.append(
                        self.cursor.add_highlight(self.spectraPlot.peak_overlays_dict[hi])
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
            self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_nodes'].set_visible(True)
            if not isinstance(self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_edges'], list):
                self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_edges'].set_visible(True)
            for k, v  in self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_labels'].items():                                    
                v.set_visible(True)
            self.old_node = lbl

    def hide_hmbc_graph_networks(self,lbls=None):
        
        if isinstance(lbls,str):
            if lbls[0] == 'H':
                return
            lbls_list = [lbls]
        else:
            lbls_list = self.nmrproblem.hmbc_graphs.keys()

        for lbl in lbls_list:
            self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_nodes'].set_visible(False)
            if not isinstance(self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_edges'], list):
                self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_edges'].set_visible(False)
            for v in self.moleculePlot.hmbc_graph_plots[lbl]['hmbc_labels'].values():                                    
                v.set_visible(False) 

        self.old_node = False      


def main(nmrproblem):

    app = QApplication(sys.argv)
    ex = MainWidget(nmrproblem)
    ex.show()
    sys.exit(app.exec_())


if __name__ == '__main__':

    data_info = nmrProblem.parse_argv()


    print(data_info)

    nmrproblem = nmrProblem.NMRproblem(data_info)


    if nmrproblem.data_complete:

        nmrProblem.build_model(nmrproblem)
        nmrProblem.build_molecule_graph_network(nmrproblem)
        nmrProblem.build_xy3_representation_of_molecule(nmrproblem)


    main(nmrproblem)