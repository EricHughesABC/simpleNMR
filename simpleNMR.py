import sys
import os
import json

from PyQt5.QtWidgets import (
    QWidget,
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
    QAction,
    QFileDialog
)

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
# from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
# import matplotlib.pyplot as plt
import mplcursors

import rdkit
from rdkit import Chem
# from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import networkx as nx
# import nx_pylab
import numpy as np
import pandas as pd

import PIL
from PIL import Image

import nmrProblem
import qt5_tabs_001 
from qt5_tabs_001 import EditDataFrameDialog
from moleculePlot import MatplotlibMoleculePlot
from spectraPlot import MatplotlibH1C13Plot

from numbers import Number


# def create_hmbc_graph_fragments(nmrproblem, hmbc_edges):
#     hmbc_graphs = {}
#     ntwk_labels = []

#     catoms = nmrproblem.carbonAtoms
#     df = nmrproblem.df
#     xy3 = nmrproblem.xy3
#     molecule = nmrproblem.molecule
#     udic = nmrproblem.udic

#     ret1 = None
#     ret2 = None

#     lineCollections = []
#     hmbc_edge_colors = nmrproblem.hmbc_edge_colors
    
#     for i, c in enumerate(catoms):

#         if c not in hmbc_edges.keys():
#             continue

#         # create hmbc graph for node c and add  xy coodinates
#         hmbc_graphs[c] = {}
#         hmbc_graphs[c]["graph"] = nx.Graph()
#         hmbc_graphs[c]["xy"] = dict((k, xy3[k]) for k in [c] + list(hmbc_edges[c]))
#         hmbc_graphs[c]["colors"] = []

#         # add nodes to hmbc graph
#         hmbc_graphs[c]["graph"].add_nodes_from([c] + list(hmbc_edges[c]))

#         # add edges
#         for i, c1 in enumerate(hmbc_edges[c]):
#             hmbc_graphs[c]["graph"].add_edge(c, c1)
#             hmbc_graphs[c]["colors"].append(hmbc_edge_colors[i])

#     return hmbc_graphs


class MoleculePlotCanvas(FigureCanvasQTAgg):
    def __init__(self, fig, parent=None):
        # self.molecule_fig = Figure(figsize=(width, height), dpi=dpi)
        super(MoleculePlotCanvas, self).__init__(fig)
        # gs = GridSpec(1, 1, figure=self.fig)
        # self.ax = self.fig.add_subplot(gs[:, :])
        # self.molecule_ax = self.molecule_fig.add_subplot()
        # self.molecule_ax.plot([1,2,3], [1,2,3])


class MainWidget(QMainWindow):
    def __init__(self, nmrproblem, parent=None):
        self.nmrproblem = nmrproblem
        self.hmbc_edge_colors = nmrproblem.hmbc_edge_colors
        self.old_node = ""
        self.new_node = ""
        self.selection = None

        super().__init__(parent)

        self.setGeometry(100, 100, 1200, 900)
        self.setWindowTitle("SimpleNMR")



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
        self.spectraCanvas.setGeometry(0,0,600,1100)
        self.moleculeCanvas = MoleculePlotCanvas(self.moleculePlot)
        self.moleculeCanvas.setGeometry(0, 0, 300, 300)

        moltoolbar = NavigationToolbar(self.moleculeCanvas, self)
        spctoolbar = NavigationToolbar(self.spectraCanvas, self)

        # self.webEngineView = QWebEngineView()
        # self.loadPage()

        self.smilesInput = QLineEdit()
        self.button = QPushButton("Display Smiles Molecule", self)

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

        # jsmevbox.addWidget(self.webEngineView)
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
        splitter2.setSizes([700,400])

        # hover=mplcursors.HoverMode.Transient
        self.cursor = mplcursors.cursor(
            [self.moleculePlot.mol_nodes]
            + self.spectraPlot.peak_overlays[0]
            + self.spectraPlot.peak_overlays[1],
            hover=mplcursors.HoverMode.Transient,
            highlight=True,
        )

        self.setup_mplcursors()

        self.centralWidget.show()

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
        helpMenu = menuBar.addMenu( "&Help")
        helpMenu.addAction(self.helpContentAction)
        helpMenu.addAction(self.aboutAction)

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
        editToolBar.addWidget(self.fontSizeSpinBox)

    def _createStatusBar(self):
        self.statusbar = self.statusBar()
        # Temporary message
        self.statusbar.showMessage("Ready", 3000)
        # Permanent widget
        self.wcLabel = QLabel(f"Molecule: {self.nmrproblem.moleculeAtomsStr} DBE: {int(self.nmrproblem.dbe)}")
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
    #     self.centralWidget.addAction(self.copyAction)
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
        print("New MresNova")
        folderpath = QFileDialog.getExistingDirectory(self, 'Select/Create New Problem Folder')
        _, excel_fn = os.path.split(folderpath)
        excel_fn = excel_fn + ".xlsx"
        print("folderpath", folderpath, excel_fn)

        dlg = EditDataFrameDialog(self.nmrproblem)

        if dlg.exec():
            print("Success!")
            writer = pd.ExcelWriter(os.path.join(folderpath, excel_fn), engine='xlsxwriter')
            for sheetname, df in nmrProblem.new_dataframes.items():


                print(sheetname)
                print(df)

                df.to_excel(writer, sheet_name=sheetname)
            # print(type(dlg.table_widget))
            writer.save()


            workingdir, fn = os.path.split(folderpath)
            data_info = nmrProblem.parse_argv([fn, workingdir, fn])
            nmrproblem = nmrProblem.NMRproblem(data_info)

            if nmrproblem.data_complete:

                define_hsqc_f2integral(nmrproblem)
                define_c13_attached_protons(nmrproblem)

                nmrProblem.build_model(nmrproblem)
                nmrProblem.build_molecule_graph_network(nmrproblem)
                nmrProblem.build_xy3_representation_of_molecule(nmrproblem)

                self.initiate_windows(nmrproblem)
        else:
            print("Cancel!")

        # Logic for creating a new file goes here...
        # self.csideWidget.setText("<b>File > New</b> clicked")

    def newFromTopspin(self):
        print("New Topspin")
        
        # Logic for creating a new file goes here...
        # self.csideWidget.setText("<b>File > New</b> clicked")

    def openFile(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            workingdir, fn = os.path.split(fileName)
            data_info = nmrProblem.parse_argv([fn, workingdir, fn])
            nmrproblem = nmrProblem.NMRproblem(data_info)

            if nmrproblem.data_complete:

                define_hsqc_f2integral(nmrproblem)
                define_c13_attached_protons(nmrproblem)

                nmrProblem.build_model(nmrproblem)
                nmrProblem.build_molecule_graph_network(nmrproblem)
                nmrProblem.build_xy3_representation_of_molecule(nmrproblem)

                self.initiate_windows(nmrproblem)
        # Logic for opening an existing file goes here...
        # self.centralWidget.setText("<b>File > Open...</b> clicked")

    def saveFile(self):
        for n in self.nmrproblem.molecule.nodes:
            print(n, self.moleculePlot.hmbc_graph_plots[n]["hmbc_nodes"].get_offsets()[0])
            self.nmrproblem.xy3[n] = self.moleculePlot.hmbc_graph_plots[n]["hmbc_nodes"].get_offsets()[0]

        if hasattr(self.nmrproblem, 'xy'):
            print( "xy\n", self.nmrproblem.xy)
        if hasattr(self.nmrproblem, 'xy3'):
            print( "xy3\n", self.nmrproblem.xy3)

            for k,v in self.nmrproblem.xy3.items():
                if isinstance(v, np.ndarray):
                    self.nmrproblem.xy3[k] = v.tolist()

            print(os.path.join(self.nmrproblem.problemDirectoryPath, "xy3.json"))
            with open(os.path.join(self.nmrproblem.problemDirectoryPath, "xy3.json"), "w") as fp:
                json.dump(self.nmrproblem.xy3, fp, indent=2)
                print("saving xy3.json")

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

    def loadPage(self):

        url = QUrl.fromUserInput("https://jsme-editor.github.io/dist/JSME_minimal.html")

        if url.isValid():
            self.webEngineView.load(url)

    def setup_mplcursors(self):
        @self.cursor.connect("remove")
        def on_remove(sel):
            if "molecule" in sel.artist.axes.get_gid():
                lbl = "C" + str(int(sel.index) + 1)
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

            lbl = "C" + str(int(sel.index) + 1)
            # print(dir(sel))
            # print("sel.artist.axes.get_label()", sel.artist.axes.get_label())
            # print("sel.artist.axes.get_gid()", sel.artist.axes.get_gid())
            # print("sel.target\n", sel.target)

            if "molecule" in sel.artist.axes.get_gid():
                lbl = "C" + str(int(sel.index) + 1)

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
                sel.extras.append(
                    self.cursor.add_highlight(self.spectraPlot.peak_overlays_dict[lbl])
                )
                sel.extras[-1].set_visible(True)
                sel.extras[-1].set_linewidth(0.75)
                sel.extras[-1].set_color("red")

                highlighted_H1_pks = udic[1]["info"].loc[lbl, "hsqc"]

                # highlight corresponding H1 HSQC peaks ins 1D proton Spectrum
                for hpk in highlighted_H1_pks:
                    sel.extras.append(
                        self.cursor.add_highlight(
                            self.spectraPlot.peak_overlays_dict[hpk]
                        )
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
                self.new_node = lbl
                if str(lbl) in hatoms:
                    ii = 0
                else:
                    ii = 1

                x, y = sel.target

                # use artist to labal to find out peak id, H1, H2, ... or C1, C2

                selected_pk = [str(lbl)]
                highlighted_pks = udic[ii]["info"].loc[selected_pk[0], "hsqc"]

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

                # highlight coresponding proton or carbon peaks
                # from hsqc data find corresponding peaks
                highlighted_pks = udic[ii]["info"].loc[selected_pk[0], "hsqc"]
                for hpk in highlighted_pks:
                    sel.extras.append(
                        self.cursor.add_highlight(
                            self.spectraPlot.peak_overlays_dict[hpk]
                        )
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
                    if lbl in self.nmrproblem.hsqcH1labelC13label.keys():
                        lbl = self.nmrproblem.hsqcH1labelC13label[lbl]

                self.hide_hmbc_graph_networks()
                self.draw_hmbc_graph_network(lbl)
                self.redraw_axes()

    @pyqtSlot()
    def on_click(self):
        smilesText = self.smilesInput.text()
        QMessageBox.question(
            self, "Smiles", "Molecule: " + smilesText, QMessageBox.Ok, QMessageBox.Ok
        )


        with open(os.path.join(self.nmrproblem.problemDirectoryPath, self.nmrproblem.problemDirectory+".smi"), "w") as fp:
            fp.write(smilesText)

        # save smiles string to nmrproblem
        self.nmrproblem.smiles = smilesText

        molecule = Chem.MolFromSmiles(smilesText)
        Draw.MolToFile(molecule, "molecule.png")
        self.nmrproblem.png = Image.open("molecule.png").transpose(
            PIL.Image.FLIP_LEFT_RIGHT
        )

        if hasattr(self.moleculePlot, "bkgnd"):
            self.moleculePlot.bkgnd.set_data(self.nmrproblem.png)
        else:
            self.moleculePlot.bkgnd = self.moleculePlot.ax.imshow(
                np.fliplr(self.nmrproblem.png),
                aspect="auto",
                extent=[1.0*self.moleculePlot.xmax, 1.0*self.moleculePlot.xmin, 1.5*self.moleculePlot.ymin, 1.5*self.moleculePlot.ymax],
                alpha=0.4,
            )


    def highlight_hmbc_peaks(self, lbl, sel):
        print("lbl", lbl)
        if lbl[0] == "H":
            if lbl in self.nmrproblem.hsqcH1labelC13label.keys():
                lbl = self.nmrproblem.hsqcH1labelC13label[lbl]

        if lbl in self.nmrproblem.hmbc_graph_edges.keys():
            for i, ci in enumerate(self.nmrproblem.hmbc_graph_edges[lbl]):
                print('i, ci', i, ci, len(self.hmbc_edge_colors))
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


def define_hsqc_f2integral(nmrproblem):

    h1 = nmrproblem.h1
    hsqc = nmrproblem.hsqc

    for i in h1.index:
        if i in hsqc.index:
            print(i, np.round(h1.loc[hsqc.loc[i, "f2_i"], "integral"]))
            hsqc.loc[i, "f2_integral"] = int(
                np.round(h1.loc[hsqc.loc[i, "f2_i"], "integral"])
            )

def define_c13_attached_protons(nmrproblem):
    c13 = nmrproblem.c13
    hsqc = nmrproblem.hsqc

    c13["attached_protons"] = 0
    
    for i in c13.ppm.index:
        dddf = hsqc[hsqc.f1_ppm == c13.loc[i, "ppm"]]

        if dddf.shape[0]:
            c13.loc[i, "attached_protons"] = int(dddf.f2_integral.sum())

# def main(nmrproblem):

#     ex = MainWidget(nmrproblem)
#     ex.show()
#     sys.exit(app.exec_())


if __name__ == "__main__":

    app = QApplication(sys.argv)


    data_info = nmrProblem.parse_argv()
    print("data_info\n", data_info)
    nmrproblem = nmrProblem.NMRproblem(data_info)



    if nmrproblem.data_complete:
        define_hsqc_f2integral(nmrproblem)
        define_c13_attached_protons(nmrproblem)
        nmrProblem.build_model(nmrproblem)
        nmrProblem.build_molecule_graph_network(nmrproblem)
        nmrProblem.build_xy3_representation_of_molecule(nmrproblem)
        # nmrProblem.build_xy3_representation_of_molecule_from_smiles(nmrproblem)


    print(nmrproblem.c13)

    # main(nmrproblem)
    ex = MainWidget(nmrproblem)

    ex.show()
    sys.exit(app.exec_())
