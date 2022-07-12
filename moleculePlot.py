import sys

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import nx_pylab
import numpy as np

import nmrProblem

from PyQt5.QtWidgets import (
    QWidget,
    QMainWindow,
    QApplication,
    QHBoxLayout,
)

from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import pyqtSlot

from spectraPlot import MatplotlibH1C13Plot

PLOTLINECOLORS = ('blue', 'orange', 'green', 'red', 'purple')
SCATTERFACECOLORS = ('blue', 'orange', 'green', 'red', 'purple')
SCATTEREDGECOLORS = ('blue', 'orange', 'green', 'red', 'purple')

YELLOW = (1., 1., 0., 1.)
RED    = (1., 0., 0., 1.)
WHITE  = (1., 1., 1., 1.)




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
            if None in self.nmrproblem.hmbc_graphs[k]["graph"].nodes():
                self.nmrproblem.hmbc_graphs[k]["graph"].remove_node(None)

        self.mol_edges, self.mol_nodes, self.mol_labels = self.init_moleculePlots(
            self.ax, self.nmrproblem
        )
        self.mol_vertices_moved = self.init_mol_vertices_moved(self.mol_edges)

        self.hmbc_graph_plots = self.init_hmbc_graph_plots(
            self.ax, self.nmrproblem.hmbc_graphs
        )

        self.ax.set_xlim(self.xmin, self.xmax)
        self.ax.set_ylim(self.ymin, self.ymax)

        if not isinstance(nmrproblem.png, type(None)):
            self.bkgnd = self.ax.imshow(
                np.fliplr(nmrproblem.png),
                aspect="auto",
                extent=[
                    1.0 * self.xmax,
                    1.0 * self.xmin,
                    1.5 * self.ymin,
                    1.5 * self.ymax,
                ],
                alpha=0.4,
            )


    def draw_hmbc_graph_network(self, lbl:str):
        if lbl in self.nmrproblem.hmbc_graphs.keys():
            self.hmbc_graph_plots[lbl]["hmbc_nodes"].set_visible(True)
            if not isinstance(
                self.hmbc_graph_plots[lbl]["hmbc_edges"], list
            ):
                self.hmbc_graph_plots[lbl]["hmbc_edges"].set_visible(True)
            for k, v in self.hmbc_graph_plots[lbl]["hmbc_labels"].items():
                v.set_visible(True)
            self.old_node = lbl

    def hide_hmbc_graph_networks(self, lbls:str=None):


        if isinstance(lbls, str):
            if lbls[0] == "H":
                return
            lbls_list = [lbls]
        else:
            lbls_list = self.nmrproblem.hmbc_graphs.keys()

        for lbl in lbls_list:
            self.hmbc_graph_plots[lbl]["hmbc_nodes"].set_visible(False)
            if not isinstance(
                self.hmbc_graph_plots[lbl]["hmbc_edges"], list
            ):
                self.hmbc_graph_plots[lbl]["hmbc_edges"].set_visible(False)
            for v in self.hmbc_graph_plots[lbl]["hmbc_labels"].values():
                v.set_visible(False)

        self.old_node = False        


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

    def init_moleculePlots(self, ax, nmrproblem):

        mol_edges = nx_pylab.draw_networkx_edges(
            nmrproblem.molecule,
            nmrproblem.xy3,
            ax=ax,
            edge_color='r',
            width=3,
            label="mol_edges",
        )

        mol_nodes = nx_pylab.draw_networkx_nodes(
            nmrproblem.molecule,
            nmrproblem.xy3,
            node_color=[
                nmrproblem.molecule.nodes[n]["node_color"]
                for n in nmrproblem.molecule.nodes
            ],
            edgecolors=[
                nmrproblem.molecule.nodes[n]["node_color"]
                for n in nmrproblem.molecule.nodes
            ],
            linewidths=0.2,
            node_size=500,
            ax=ax,
            label="mol_nodes",
            picker=True,
            pickradius=5,
        )
                # scatterplt = self.mol_nodes
        mol_nodes.scatter_facecolors_rgba = mol_nodes.get_facecolors()
        mol_nodes.scatter_edgecolors_rgba = mol_nodes.get_edgecolors()
        mol_nodes.my_labels = [f"C{n+1}" for n in range(len(nmrproblem.molecule.nodes))]
        mol_nodes.node_highlighted = False

        mol_labels = nx_pylab.draw_networkx_labels(
            nmrproblem.molecule, nmrproblem.xy3, ax=ax
        )

        return mol_edges, mol_nodes, mol_labels


def define_hsqc_f2integral(nmrproblem):

    h1 = nmrproblem.h1
    hsqc = nmrproblem.hsqc

    for i in h1.index:
        if i in hsqc.index:
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


if __name__ == "__main__":

    class MoleculePlotCanvas(FigureCanvasQTAgg):
        def __init__(self, fig, parent=None):
            # self.molecule_fig = Figure(figsize=(width, height), dpi=dpi)
            super(MoleculePlotCanvas, self).__init__(fig)


    class MainWidget(QMainWindow):

        def __init__(self, nmrproblem, parent=None):
            self.nmrproblem = nmrproblem
            self.centralWidget = QWidget()
            super().__init__(parent)

            self.setGeometry(100, 100, 1200, 900)
            self.setWindowTitle("moleculePlot")

            self.moleculePlot = MatplotlibMoleculePlot(nmrproblem)
            self.moleculeCanvas = MoleculePlotCanvas(self.moleculePlot)

            self.spectraPlot = MatplotlibH1C13Plot(self.nmrproblem)
            self.spectraCanvas = MoleculePlotCanvas(self.spectraPlot)

            # add callbacks to the moleculeCanvas
            self.node_pick_ind = None
            self.node_hover_ind = None
            self.node_picked = False 
            self.highlighted_peak_lbl = None

            self.moleculePlot.canvas.mpl_connect("button_release_event", 
                lambda event: self.button_release_molecule(event, 
                                                        event_name="button_release_event"))

            self.moleculePlot.canvas.mpl_connect("motion_notify_event", 
                lambda event: self.motion_notify_callback(event, 
                                                        specplot=self.spectraPlot, 
                                                        molplot=self.moleculePlot))

            self.moleculePlot.canvas.mpl_connect("pick_event", 
                lambda event: self.pick_molecule(event, 
                                                event_name="pick_event", 
                                                specplot=self.spectraPlot, 
                                                molplot=self.moleculePlot))
            
            self.spectraPlot.canvas.mpl_connect("motion_notify_event", 
                lambda event: self.hover_over_specplot(event, 
                                                    specplot=self.spectraPlot, 
                                                    molplot=self.moleculePlot))

            hbox = QHBoxLayout()
            hbox.addWidget(self.moleculeCanvas)
            hbox.addWidget(self.spectraCanvas)

            self.centralWidget.setLayout(hbox)
            self.centralWidget.show()

        def hover_over_specplot(self, event, specplot, molplot):
            in_plot = []
            in_plot_label = []
            in_plot_index = []
            pos = None
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
                specplot.peak_overlays_dict[lbl].set_linewidth(0.75)

                # annotate new peak
                if 'H' in lbl:
                    # set the annotation to the peak
                    atom_index = int(lbl[1:])
                    ppm = self.nmrproblem.h1.loc[atom_index, "ppm"] 
                    integral = self.nmrproblem.h1.loc[atom_index, "integral"]
                    jcoupling = self.nmrproblem.h1.loc[atom_index, "jCouplingClass"]
                    annot_text = f"{lbl}: {ppm:.2f} ppm\nInt:{integral}\nJ: {jcoupling}"
                    specplot.annot_H1.xy = (event.xdata, event.ydata)
                    specplot.annot_H1.set_text(annot_text)
                    specplot.annot_H1.set_visible(True)

                elif 'C' in lbl:
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
                    clbl = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f2H_i==lbl]["f1C_i"].values[0]

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

                #uddate title in molplot
                c13_ind = int(clbl[1:])
                ppm_val = self.nmrproblem.c13.loc[c13_ind]["ppm"]
                attached_protons = self.nmrproblem.c13.loc[c13_ind]["attached_protons"]
                title_str = f"{clbl}: {ppm_val:.1f} ppm, attached protons: {attached_protons}"
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

                    #unhighlight distributions
                    specplot.reset_distributions_eeh()

                    # uddate specplot canvas
                    specplot.canvas.draw_idle()
                    molplot.canvas.draw_idle()


        def update_molplot_highlights(self, molplot, specplot, lbl):

            molplot.mol_nodes.node_highlighted = True

            ind = int(lbl[1:])-1

            scatter_facecolors_highlight = molplot.mol_nodes.get_facecolors()
            scatter_edgecolors_highlight  = molplot.mol_nodes.get_edgecolors()

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
            specplot.display_annotation_H1_from_molplot(lbl, specplot.annot_H1)

            # annotate distributions
            hpks = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f1C_i==lbl]["f2H_i"].values
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


        def button_release_molecule(self, event, **event_argv):
            self.node_pick_ind = None
            self.node_hover_ind = None
            self.node_picked = False 


        def motion_notify_callback(self, event, specplot, molplot):

            self.hover_over_molecule(event, event_name="motion_notify_event", molplot=molplot, specplot=specplot)

            if self.node_pick_ind is None:
                return
            if event.inaxes is None:
                return
            if event.button != 1:
                return

            self.node_moved = True
            x, y = event.xdata, event.ydata

            self.hmbc_graph_plots = self.moleculePlot.hmbc_graph_plots
            self.mol_nodes = self.moleculePlot.mol_nodes
            self.mol_edges = self.moleculePlot.mol_edges
            self.mol_labels = self.moleculePlot.mol_labels
            self.mol_ind = self.node_pick_ind

            
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
                return

            if in_node:
                self.mol_nodes.node_highlighted = True

                ind = node_index['ind'][0]
                lbl = f"{self.mol_nodes.my_labels[ind]}"
                x, y = self.mol_nodes.get_offsets()[ind]
                self.node_hover_ind = ind
                self.node_hover_lbl = lbl
                self.node_hover_x = x
                self.node_hover_y = y
                
                c13_ind = int(lbl[1:])
                ppm_val = self.nmrproblem.c13.loc[c13_ind]["ppm"]
                attached_protons = self.nmrproblem.c13.loc[c13_ind]["attached_protons"]
                title_str = f"{lbl}: {ppm_val:.1f} ppm, attached protons: {attached_protons}"
                self.moleculePlot.ax.set_title(title_str)
            
                scatter_facecolors_highlight = self.mol_nodes.get_facecolors()
                scatter_edgecolors_highlight  = self.mol_nodes.get_edgecolors()
                
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
                specplot.display_annotation_H1_from_molplot(lbl, specplot.annot_H1)

                # annotate distributions
                hpks = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f1C_i==lbl]["f2H_i"].values
                cpks = [lbl]
                self.display_distributions(cpks, hpks, specplot)

                molplot.canvas.draw_idle()
                specplot.canvas.draw_idle()
            else:
                if self.mol_nodes.node_highlighted:
                    self.mol_nodes.node_highlighted = False

                    self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
                    self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

                    self.moleculePlot.hide_hmbc_graph_networks()

                    specplot.reset_peak_overlays_eeh()

                    specplot.hide_annotation(specplot.annot_C13)
                    specplot.hide_annotation(specplot.annot_H1)

                    #unhighlight distributions
                    specplot.reset_distributions_eeh()

                    molplot.canvas.draw_idle()
                    specplot.canvas.draw_idle()


        def pick_molecule(self, event, event_name, specplot, molplot):
            in_node, node_index = self.moleculePlot.mol_nodes.contains(event.mouseevent)

            if in_node:
                self.node_picked = True

            # self.mol_ind = event.ind
            # self.node_pick_ind = event.ind
            self.node_pick_ind = node_index['ind'][0]
            ptcoords = np.ma.compressed(self.moleculePlot.mol_nodes.get_offsets()[self.node_pick_ind])

            # find label id
            self.label_id = list(self.moleculePlot.mol_labels.keys())[self.node_pick_ind]

            # hide all highlights before dragging node
            if self.mol_nodes.node_highlighted:
                self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
                self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

                self.moleculePlot.hide_hmbc_graph_networks()

                specplot.reset_peak_overlays_eeh()

                specplot.hide_annotation(specplot.annot_C13)
                specplot.hide_annotation(specplot.annot_H1)

                #unhighlight distributions
                specplot.reset_distributions_eeh()
                specplot.canvas.draw_idle()
                molplot.canvas.draw_idle()

            # for each carbon atom in moleule check to see if it is in
            # the associated hmbc network then set the vertices_moved flag
            # to true if the picked node is in the network
            # update edges, nodes and labels
            # update them even if they are not visible so that keep in sync
            # with the moved display
            for n in self.nmrproblem.molecule.nodes:
                hmbc_nodes = self.moleculePlot.hmbc_graph_plots[n]["hmbc_nodes"]
                hmbc_edges = self.moleculePlot.hmbc_graph_plots[n]["hmbc_edges"]
                hmbc_labels = self.moleculePlot.hmbc_graph_plots[n]["hmbc_labels"]
                hmbc_vertices_moved = self.moleculePlot.hmbc_graph_plots[n]["hmbc_vertices_moved"]

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
                for i, e in enumerate(self.moleculePlot.mol_edges.get_paths()):
                    for j, c in enumerate(e.vertices):
                        if c[0] == ptcoords[0] and c[1] == ptcoords[1]:
                            self.moleculePlot.mol_vertices_moved[i][j] = True
                        else:
                            self.moleculePlot.mol_vertices_moved[i][j] = False
            else:
                for i, e in enumerate(self.moleculePlot.mol_edges.get_paths()):
                    for j, c in enumerate(e.vertices):
                        self.moleculePlot.mol_vertices_moved[i][j] = False


    app = QApplication(sys.argv)

    data_info = nmrProblem.parse_argv()
    nmrproblem = nmrProblem.NMRproblem(data_info)



    if nmrproblem.data_complete:
        define_hsqc_f2integral(nmrproblem)
        define_c13_attached_protons(nmrproblem)
        nmrProblem.build_model(nmrproblem)
        nmrProblem.build_molecule_graph_network(nmrproblem)
        nmrProblem.build_xy3_representation_of_molecule(nmrproblem)

        ex = MainWidget(nmrproblem)

        # moleculeCanvas.show()
        sys.exit(app.exec_())
