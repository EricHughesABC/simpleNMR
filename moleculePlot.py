"""moleculePlot.py"""
import sys

# from matplotlib.backends.backend_qt5agg import (
#     FigureCanvasQTAgg,
#     NavigationToolbar2QT as NavigationToolbar,
# )

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

from matplotlib.figure import Figure
from matplotlib.lines import Line2D

# import matplotlib.pyplot as plt


import numpy as np


from PyQt5.QtWidgets import (
    QWidget,
    QMainWindow,
    QApplication,
    QHBoxLayout,
)

# from PyQt5.QtCore import pyqtSlot
import nmrProblem
import nx_pylab
import simpleNMRutils
from spectraPlot import MatplotlibH1C13Plot

PLOTLINECOLORS = ("blue", "orange", "green", "red", "purple")
SCATTERFACECOLORS = ("blue", "orange", "green", "red", "purple")
SCATTEREDGECOLORS = ("blue", "orange", "green", "red", "purple")

YELLOW = (1.0, 1.0, 0.0, 1.0)
RED = (1.0, 0.0, 0.0, 1.0)
WHITE = (1.0, 1.0, 1.0, 1.0)


class MatplotlibMoleculePlot(Figure):
    """class to plot molecule in matplotlib figure enbedded in pyqt5 window"""

    def __init__(self, nmrprblm):

        self.nmrproblem = nmrprblm

        self.mol_ind = None
        self.hmbc_ind = None
        self.old_node = False
        self.mol_nodes = None
        self.mol_labels = None

        # global ptcoords
        # # global hmbc_vertices_moved
        # global mol_vertices_moved
        self.label_id = None
        # super(MatplotlibMoleculePlot, self).__init__(constrained_layout=True, figsize=(4, 4), dpi=100)
        super(MatplotlibMoleculePlot, self).__init__(constrained_layout=True)
        self.ax = self.add_subplot(label="molecule", gid="molecule_id")

        # set backfround of plot to light red
        # self.ax.set_facecolor(YELLOW)
        # self.tight_layout()

        # self.ax.set_facecolor((0.9, 0.9, 0.9))
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

        self.xmin = -0.1
        self.ymin = -0.1
        self.xmax = 1.1
        self.ymax = 1.1

        custom_lines = [
            Line2D([0], [0], color="blue", lw=4),
            Line2D([0], [0], color="orange", lw=4),
            Line2D([0], [0], color="green", lw=4),
            Line2D([0], [0], color="purple", lw=4),
        ]

        self.ax.legend(custom_lines, ["- C -", "- CH", "- CH$_2$", "- CH$_3$"])

        if self.nmrproblem.data_complete:
            self.draw_molecule(self.nmrproblem, self.ax)

    def draw_molecule(self, nmrprblm, ax):
        """draw molecule in matplotlib figure"""
        self.nmrproblem = nmrprblm
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
        print("mol_edges", type(self.mol_edges))
        self.mol_vertices_moved = self.init_mol_vertices_moved(self.mol_edges)

        self.hmbc_graph_plots = self.init_hmbc_graph_plots(
            self.ax, self.nmrproblem.hmbc_graphs
        )

        if not isinstance(self.nmrproblem.png, type(None)):
            self.bkgnd = self.ax.imshow(
                self.nmrproblem.png,
                aspect="auto",
                extent=[0, 1, 1, 0],
                alpha=0.6,
            )

        self.ax.set_xlim(-0.1, 1.1)
        self.ax.set_ylim(1.1, -0.1)

    def draw_hmbc_graph_network(self, lbl: str):
        """draw hmbc graph network"""
        if lbl in self.nmrproblem.hmbc_graphs.keys():
            self.hmbc_graph_plots[lbl]["hmbc_nodes"].set_visible(True)
            if not isinstance(self.hmbc_graph_plots[lbl]["hmbc_edges"], list):
                self.hmbc_graph_plots[lbl]["hmbc_edges"].set_visible(True)

            for hmbc_graph in self.hmbc_graph_plots[lbl]["hmbc_labels"].values():
                hmbc_graph.set_visible(True)
            self.old_node = lbl

    def hide_hmbc_graph_networks(self, lbls: str = None):
        """hide hmbc graph networks"""
        if isinstance(lbls, str):
            if lbls[0] == "H":
                return
            lbls_list = [lbls]
        else:
            lbls_list = self.nmrproblem.hmbc_graphs.keys()

        for lbl in lbls_list:
            self.hmbc_graph_plots[lbl]["hmbc_nodes"].set_visible(False)
            if not isinstance(self.hmbc_graph_plots[lbl]["hmbc_edges"], list):
                self.hmbc_graph_plots[lbl]["hmbc_edges"].set_visible(False)
            for hmbc_graph in self.hmbc_graph_plots[lbl]["hmbc_labels"].values():
                hmbc_graph.set_visible(False)

        self.old_node = False

    def init_hmbc_graph_plots(self, ax, hmbc_graphs):
        """initialize hmbc graph plots"""
        hmbc_graph_plots = {}

        for node in self.nmrproblem.nx_graph_molecule.nodes:
            hmbc_graph_plots[node] = {}
        for node in self.nmrproblem.nx_graph_molecule.nodes:
            if node not in hmbc_graphs:
                continue
            hmbc_graph_plots[node]["hmbc_nodes"] = nx_pylab.draw_networkx_nodes(
                hmbc_graphs[node]["graph"],
                hmbc_graphs[node]["xy"],
                ax=ax,
                label=node + "_node",
                node_color="w",
                edgecolors=["r"] + hmbc_graphs[node]["colors"],
                node_size=500,
                linewidths=4,
                picker=False,
            )

            hmbc_graph_plots[node]["hmbc_nodes"].set_visible(False)

            hmbc_graph_plots[node]["hmbc_edges"] = nx_pylab.draw_networkx_edges(
                hmbc_graphs[node]["graph"],
                hmbc_graphs[node]["xy"],
                # edge_color=hmbc_graphs[node]["colors"],
                edge_color="grey",
                width=5,
                ax=ax,
                label=node + "_edge",
            )

            hmbc_graph_plots[node]["hmbc_labels"] = nx_pylab.draw_networkx_labels(
                hmbc_graphs[node]["graph"], hmbc_graphs[node]["xy"], ax=ax
            )

            for k, vertex in hmbc_graph_plots[node]["hmbc_labels"].items():
                vertex.set_visible(False)

            if not isinstance(hmbc_graph_plots[node]["hmbc_edges"], list):
                hmbc_graph_plots[node]["hmbc_edges"].set_visible(False)

                hmbc_graph_plots[node]["hmbc_vertices_moved"] = []

                for edge_path in hmbc_graph_plots[node]["hmbc_edges"].get_paths():
                    vertices = []
                    for vertex in edge_path.vertices:
                        vertices.append(False)
                    hmbc_graph_plots[node]["hmbc_vertices_moved"].append(vertices)

            else:
                hmbc_graph_plots[node]["hmbc_vertices_moved"] = [False]

        return hmbc_graph_plots

    def init_mol_vertices_moved(self, mol_edges):
        """initialize molecule vertices moved"""
        mol_vertices_moved = []
        if not isinstance(mol_edges, list):
            for e in mol_edges.get_paths():
                vertices = []
                for c in e.vertices:
                    vertices.append(False)
                mol_vertices_moved.append(vertices)

        return mol_vertices_moved

    def init_moleculePlots(self, ax, nmrprblm):
        """initialize molecule plots"""

        mol_edges = nx_pylab.draw_networkx_edges(
            nmrprblm.nx_graph_molecule,
            nmrprblm.xy3,
            ax=ax,
            edge_color="r",
            width=3,
            label="mol_edges",
        )

        mol_nodes = nx_pylab.draw_networkx_nodes(
            nmrprblm.nx_graph_molecule,
            nmrprblm.xy3,
            node_color=[
                nmrprblm.nx_graph_molecule.nodes[node]["node_color"]
                for node in nmrprblm.nx_graph_molecule.nodes
            ],
            edgecolors=[
                nmrprblm.nx_graph_molecule.nodes[node]["node_color"]
                for node in nmrprblm.nx_graph_molecule.nodes
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
        mol_nodes.my_labels = [f"C{n+1}" for n in range(len(nmrprblm.nx_graph_molecule.nodes))]
        mol_nodes.node_highlighted = False

        mol_labels = nx_pylab.draw_networkx_labels(
            nmrprblm.nx_graph_molecule, nmrprblm.xy3, ax=ax
        )

        return mol_edges, mol_nodes, mol_labels


def define_hsqc_f2integral(nmrprblm):
    """define hsqc f2 integral"""

    h1 = nmrprblm.h1
    hsqc = nmrprblm.hsqc

    for i in h1.index:
        if i in hsqc.index:
            hsqc.loc[i, "f2_integral"] = int(
                np.round(h1.loc[hsqc.loc[i, "f2_i"], "integral"])
            )


def define_c13_attached_protons(nmrprblm):
    """define c13 attached protons"""
    c13 = nmrprblm.c13
    hsqc = nmrprblm.hsqc

    c13["attached_protons"] = 0

    for i in c13.ppm.index:
        dddf = hsqc[hsqc.f1_ppm == c13.loc[i, "ppm"]]

        if dddf.shape[0]:
            c13.loc[i, "attached_protons"] = int(dddf.f2_integral.sum())


if __name__ == "__main__":

    class MoleculePlotCanvas(FigureCanvasQTAgg):
        """MoleculePlotCanvas"""

        def __init__(self, fig, parent=None):
            # self.molecule_fig = Figure(figsize=(width, height), dpi=dpi)
            super(MoleculePlotCanvas, self).__init__(fig)

    class MainWidget(QMainWindow):
        """MainWidget"""

        def __init__(self, nmrproblem, parent=None):
            self.nmrproblem = nmrproblem
            self.centralWidget = QWidget()
            super().__init__(parent)

            self.setGeometry(100, 100, 1200, 900)
            self.setWindowTitle("moleculePlot")

            self.molecule_plot = MatplotlibMoleculePlot(nmrproblem)
            self.molecule_canvas = MoleculePlotCanvas(self.molecule_plot)

            self.spectra_plot = MatplotlibH1C13Plot(self.nmrproblem)
            self.spectra_canvas = MoleculePlotCanvas(self.spectra_plot)

            # add callbacks to the moleculeCanvas
            self.node_pick_ind = None
            self.node_hover_ind = None
            self.node_picked = False
            self.highlighted_peak_lbl = None
            self.node_hover_lbl = None
            self.node_hover_x = None
            self.node_hover_y = None

            self.molecule_plot.canvas.mpl_connect(
                "button_release_event",
                lambda event: self.button_release_molecule(
                    event, event_name="button_release_event"
                ),
            )

            self.molecule_plot.canvas.mpl_connect(
                "motion_notify_event",
                lambda event: self.motion_notify_callback(
                    event, specplot=self.spectra_plot, molplot=self.molecule_plot
                ),
            )

            self.molecule_plot.canvas.mpl_connect(
                "pick_event",
                lambda event: self.pick_molecule(
                    event,
                    event_name="pick_event",
                    specplot=self.spectra_plot,
                    molplot=self.molecule_plot,
                ),
            )

            self.spectra_plot.canvas.mpl_connect(
                "motion_notify_event",
                lambda event: self.hover_over_specplot(
                    event, specplot=self.spectra_plot, molplot=self.molecule_plot
                ),
            )

            hbox = QHBoxLayout()
            hbox.addWidget(self.molecule_canvas)
            hbox.addWidget(self.spectra_canvas)

            self.centralWidget.setLayout(hbox)
            self.centralWidget.show()

        def hover_over_specplot(self, event, specplot, molplot):
            """hover over specplot"""
            in_plot = []
            in_plot_label = []
            in_plot_index = []
            # pos = None
            for k, peak_overlay in specplot.peak_overlays_dict.items():
                # in_c13plots, c13plots_index = v["highlight"].contains(event)
                in_c13plots, c13plots_index = peak_overlay.contains(event)
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
                if "H" in lbl:
                    # set the annotation to the peak
                    atom_index = int(lbl[1:])
                    ppm = self.nmrproblem.h1.loc[atom_index, "ppm"]
                    integral = float(self.nmrproblem.h1.loc[atom_index, "integral"])
                    jcoupling = self.nmrproblem.h1.loc[atom_index, "jCouplingClass"]
                    jcouplingvals = self.nmrproblem.h1.loc[atom_index, "jCouplingVals"]
                    annot_text = f"{lbl}: {ppm:.2f} ppm\nInt:{integral:.1f}\nJ: {jcoupling}: {jcouplingvals}"
                    print("annot_text", annot_text)
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
                    clbl = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f2H_i == lbl][
                        "f1C_i"
                    ].values[0]

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
                self.molecule_plot.ax.set_title(title_str)
            else:
                # unhilight old peak
                if self.highlighted_peak_lbl is not None:
                    self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
                    self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

                    self.molecule_plot.hide_hmbc_graph_networks()

                    specplot.reset_peak_overlays_eeh()

                    specplot.hide_annotation(specplot.annot_C13)
                    specplot.hide_annotation(specplot.annot_H1)

                    # unhighlight distributions
                    specplot.reset_distributions_eeh()

                    # uddate specplot canvas
                    specplot.canvas.draw_idle()
                    molplot.canvas.draw_idle()

        def update_molplot_highlights(self, molplot, specplot, lbl):
            """update molplot highlights"""

            molplot.mol_nodes.node_highlighted = True

            ind = int(lbl[1:]) - 1

            scatter_facecolors_highlight = molplot.mol_nodes.get_facecolors()
            scatter_edgecolors_highlight = molplot.mol_nodes.get_edgecolors()

            scatter_facecolors_highlight[ind] = WHITE
            scatter_edgecolors_highlight[ind] = RED
            molplot.mol_nodes.set_fc(scatter_facecolors_highlight)

            if lbl in self.nmrproblem.hmbc_graphs.keys():
                # self.hide_hmbc_graph_networks()
                self.molecule_plot.draw_hmbc_graph_network(lbl)
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
            specplot.display_annotation_H1_from_molplot(lbl, specplot.annot_H1, self.nmrproblem)

            # annotate distributions
            hpks = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f1C_i == lbl][
                "f2H_i"
            ].values
            cpks = [lbl]
            self.display_distributions(cpks, hpks, specplot)

            molplot.canvas.draw_idle()
            specplot.canvas.draw_idle()

        def display_distributions(self, carbon_pks, hydrogen_pks, specplot):
            """display distributions"""
            # add highlighted distributions
            colors = ["b", "g", "r", "c", "m", "y", "k"]

            # used to dsplay legends of highlighted distributions
            plines_h = []
            plabels_h = []
            plines_c = []
            plabels_c = []
            ppm_h = []
            ppm_c = []

            # set visible
            # circle through colors
            # create proton distribution legends
            for peak in hydrogen_pks:
                numplots = len(specplot.h1c13distlist[0][peak]) - 1
                for i, mpl_component in enumerate(specplot.h1c13distlist[0][peak]):
                    ppm_h.append(self.nmrproblem.udic[0]["info"].loc[peak, "ppm"])
                    # sel.extras.append(self.cursor.add_highlight(aa))
                    mpl_component.set_visible(True)
                    mpl_component.set_linewidth(0.75)
                    mpl_component.set_color(colors[i % 7])
                    # do not add legend info if plot is just single line showing where peak pos is
                    if i < numplots:
                        plabels_h.append(mpl_component.get_label())
                        plines_h.append(mpl_component)

            # create carbon distribution legends
            for peak in carbon_pks:
                numplots = len(specplot.h1c13distlist[1][peak]) - 1
                for i, mpl_component in enumerate(specplot.h1c13distlist[1][peak]):
                    ppm_c.append(self.nmrproblem.udic[1]["info"].loc[peak, "ppm"])
                    # sel.extras.append(self.cursor.add_highlight(aa))
                    mpl_component.set_visible(True)
                    mpl_component.set_linewidth(0.75)
                    mpl_component.set_color(colors[i % 7])
                    if i < numplots:
                        plabels_c.append(mpl_component.get_label())
                        plines_c.append(mpl_component)

            # adjust x-axis width H1 and C13 of distribution plots
            # and add legend information
            if len(ppm_h) > 0:
                # calculate average position of peaks which will be used to adjust the x axis ppm range
                ppmmm_h = np.mean(ppm_h)
                specplot.h1dist_ax.set_xlim(ppmmm_h + 2, ppmmm_h - 2)
                specplot.h1dist_ax.legend(plines_h, plabels_h)
            if len(ppm_c) > 0:
                ppmmm_c = np.mean(ppm_c)
                specplot.c13dist_ax.set_xlim(ppmmm_c + 50, ppmmm_c - 50)
                specplot.c13dist_ax.legend(plines_c, plabels_c)

        def button_release_molecule(self, event, **event_argv):
            """button release event for molecule plot"""
            self.node_pick_ind = None
            self.node_hover_ind = None
            self.node_picked = False

        def motion_notify_callback(self, event, specplot, molplot):
            """motion notify event for molecule plot"""
            self.hover_over_molecule(
                event,
                event_name="motion_notify_event",
                molplot=molplot,
                specplot=specplot,
            )

            if self.node_pick_ind is None:
                return
            if event.inaxes is None:
                return
            if event.button != 1:
                return

            self.node_moved = True
            x, y = event.xdata, event.ydata

            self.hmbc_graph_plots = self.molecule_plot.hmbc_graph_plots
            self.mol_nodes = self.molecule_plot.mol_nodes
            self.mol_edges = self.molecule_plot.mol_edges
            self.mol_labels = self.molecule_plot.mol_labels
            self.mol_ind = self.node_pick_ind

            # move nodes, edges and labels associated with hmbc network for each carbon atom
            for node in self.nmrproblem.nx_graph_molecule.nodes:
                if node in self.hmbc_graph_plots:

                    hmbc_nodes = self.hmbc_graph_plots[node]["hmbc_nodes"]
                    hmbc_edges = self.hmbc_graph_plots[node]["hmbc_edges"]
                    hmbc_labels = self.hmbc_graph_plots[node]["hmbc_labels"]
                    hmbc_vertices_moved = self.hmbc_graph_plots[node][
                        "hmbc_vertices_moved"
                    ]

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
                            for i, edge_paths in enumerate(hmbc_edges.get_paths()):
                                vertices = []
                                for j, vertex in enumerate(edge_paths.vertices):
                                    if hmbc_vertices_moved[i][j] is True:
                                        vertices.append([x, y])
                                    else:
                                        vertices.append(vertex)

                                verts.append(vertices)
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

            for i, edge_path in enumerate(self.mol_edges.get_paths()):
                vertices = []
                for j, vertex in enumerate(edge_path.vertices):
                    if self.molecule_plot.mol_vertices_moved[i][j] is True:
                        vertices.append([x, y])
                    else:
                        vertices.append(vertex)

                verts.append(vertices)

            # readjust coords of edges in molecule network
            self.mol_edges.set_verts(verts)

            self.molecule_plot.ax.figure.canvas.draw_idle()

            self.moved_node = self.mol_ind
            self.moved_x = x
            self.moved_y = y

        def hover_over_molecule(self, event, event_name, molplot, specplot):
            """hover over event for molecule plot"""
            self.mol_nodes = molplot.mol_nodes
            self.mol_labels = molplot.mol_labels

            in_node, node_index = self.mol_nodes.contains(event)

            if self.node_picked:
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

                c13_ind = int(lbl[1:])
                ppm_val = self.nmrproblem.c13.loc[c13_ind]["ppm"]
                attached_protons = self.nmrproblem.c13.loc[c13_ind]["attached_protons"]
                title_str = (
                    f"{lbl}: {ppm_val:.1f} ppm, attached protons: {attached_protons}"
                )
                self.molecule_plot.ax.set_title(title_str)

                scatter_facecolors_highlight = self.mol_nodes.get_facecolors()
                scatter_edgecolors_highlight = self.mol_nodes.get_edgecolors()

                scatter_facecolors_highlight[ind] = WHITE
                scatter_edgecolors_highlight[ind] = RED
                self.mol_nodes.set_fc(scatter_facecolors_highlight)

                if lbl in self.nmrproblem.hmbc_graphs.keys():
                    # self.hide_hmbc_graph_networks()
                    self.molecule_plot.draw_hmbc_graph_network(lbl)
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
                specplot.display_annotation_H1_from_molplot(lbl, specplot.annot_H1, self.nmrproblem)

                # annotate distributions
                hpks = self.nmrproblem.hsqc[self.nmrproblem.hsqc.f1C_i == lbl][
                    "f2H_i"
                ].values
                cpks = [lbl]
                self.display_distributions(cpks, hpks, specplot)

                molplot.canvas.draw_idle()
                specplot.canvas.draw_idle()
            else:
                if self.mol_nodes.node_highlighted:
                    self.mol_nodes.node_highlighted = False

                    self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
                    self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

                    self.molecule_plot.hide_hmbc_graph_networks()

                    specplot.reset_peak_overlays_eeh()

                    specplot.hide_annotation(specplot.annot_C13)
                    specplot.hide_annotation(specplot.annot_H1)

                    # unhighlight distributions
                    specplot.reset_distributions_eeh()

                    molplot.canvas.draw_idle()
                    specplot.canvas.draw_idle()

        def pick_molecule(self, event, event_name, specplot, molplot):
            """pick event for molecule plot"""
            in_node, node_index = self.molecule_plot.mol_nodes.contains(
                event.mouseevent
            )

            if in_node:
                self.node_picked = True

            # self.mol_ind = event.ind
            # self.node_pick_ind = event.ind
            self.node_pick_ind = node_index["ind"][0]
            ptcoords = np.ma.compressed(
                self.molecule_plot.mol_nodes.get_offsets()[self.node_pick_ind]
            )

            # find label id
            self.label_id = list(self.molecule_plot.mol_labels.keys())[
                self.node_pick_ind
            ]

            # hide all highlights before dragging node
            if self.mol_nodes.node_highlighted:
                self.mol_nodes.set_fc(self.mol_nodes.scatter_facecolors_rgba)
                self.mol_nodes.set_ec(self.mol_nodes.scatter_edgecolors_rgba)

                self.molecule_plot.hide_hmbc_graph_networks()

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
            for node in self.nmrproblem.nx_graph_molecule.nodes:
                # hmbc_nodes = self.molecule_plot.hmbc_graph_plots[n]["hmbc_nodes"]
                hmbc_edges = self.molecule_plot.hmbc_graph_plots[node]["hmbc_edges"]
                hmbc_labels = self.molecule_plot.hmbc_graph_plots[node]["hmbc_labels"]
                hmbc_vertices_moved = self.molecule_plot.hmbc_graph_plots[node][
                    "hmbc_vertices_moved"
                ]

                # check if picked node is in hmbc network and then
                # if the picked coords match of the edges match set them to True
                # so that the values of the edges will be moved during the motion notify
                # event
                if self.label_id in hmbc_labels.keys():
                    self.hmbc_ind = [list(hmbc_labels.keys()).index(self.label_id)]
                    if not isinstance(hmbc_edges, list):
                        for i, edge_path in enumerate(hmbc_edges.get_paths()):
                            for j, vertex in enumerate(edge_path.vertices):
                                if (
                                    vertex[0] == ptcoords[0]
                                    and vertex[1] == ptcoords[1]
                                ):
                                    hmbc_vertices_moved[i][j] = True
                                else:
                                    hmbc_vertices_moved[i][j] = False
                else:
                    # if the hmbc fragment does not contain the picked node
                    # set all the moved_vertices to False
                    if not isinstance(hmbc_edges, list):
                        for i, edge_path in enumerate(hmbc_edges.get_paths()):
                            for j, vertex in enumerate(edge_path.vertices):
                                hmbc_vertices_moved[i][j] = False
                        self.hmbc_ind = None

            # do the same for the molecule network
            if self.label_id in self.molecule_plot.mol_labels.keys():
                # for i, e in enumerate(self.molecule_plot.mol_edges.get_paths()):svg
                for i, edge_path in enumerate(self.molecule_plot.mol_edges.get_paths()):
                    for j, vertex in enumerate(edge_path.vertices):
                        if vertex[0] == ptcoords[0] and vertex[1] == ptcoords[1]:
                            self.molecule_plot.mol_vertices_moved[i][j] = True
                        else:
                            self.molecule_plot.mol_vertices_moved[i][j] = False
            else:
                for i, edge_path in enumerate(self.molecule_plot.mol_edges.get_paths()):
                    for j, vertex in enumerate(edge_path.vertices):
                        self.molecule_plot.mol_vertices_moved[i][j] = False

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
