import sys

# 0345 050 4585

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

# import nx_pylab
import numpy as np
import pandas as pd

import nmrProblem


def createH1C13interactivePlot(nmrproblem, h1c13distlist, ax0):

    # w1 = widgets.Output()

    udic = nmrproblem.udic
    if "info" not in udic[0]:
        print("info not in udic[0]")
        return
    peak_overlays = []
    peak_overlays_dict = {}

    # for proton and carbon spectra
    for i in range(udic["ndim"]):
        peak_overlays1 = []
        for Hi in udic[i]["info"].index:
            # print("Hi", Hi)

            il = int(udic[i]["info"].loc[Hi, "pk_left"])
            ir = int(udic[i]["info"].loc[Hi, "pk_right"])

            (pk,) = ax0[1 - i].plot(
                udic[i]["axis"].ppm_scale()[il:ir],
                udic[i]["spec"][il:ir],
                lw=0.5,
                c="black",
                label=Hi,
                gid=Hi,
            )

            peak_overlays1.append(pk)
            peak_overlays_dict[Hi] = pk

        peak_overlays.append(peak_overlays1)

    return peak_overlays_dict, peak_overlays


class MatplotlibH1C13Plot(Figure):
    def __init__(self, nmrproblem):

        self.nmrproblem = nmrproblem

        super(MatplotlibH1C13Plot, self).__init__(
            constrained_layout=True, figsize=(4, 4), dpi=100
        )
        gs = GridSpec(2, 6, figure=self)

        self.c13spec_ax = self.add_subplot(
            gs[0, :4], label="C13_1Dspect", gid="C13_1Dspect_id"
        )  # carbon 1D spectrum
        self.c13dist_ax = self.add_subplot(
            gs[0, 4:], label="C13_1Ddist", gid="C13_1Ddist_id"
        )  # carbon distribution
        self.h1spec_ax = self.add_subplot(
            gs[-1, :4], label="H1_1Dspect", gid="H1_1Dspect_id"
        )  # proton 1D spectrum
        self.h1dist_ax = self.add_subplot(
            gs[-1, 4:], label="H1_1Ddist", gid="H1_1Ddist_id"
        )  # proton ppm distribution


        if self.nmrproblem.data_complete:
            self.draw_spectra(self.nmrproblem, self.c13spec_ax, self.c13dist_ax, self.h1spec_ax, self.h1dist_ax)

    def draw_spectra(self, nmrproblem, c13spec_ax, c13dist_ax, h1spec_ax, h1dist_ax):

        self.nmrproblem = nmrproblem 
        self.c13spec_ax = c13spec_ax
        self.c13dist_ax = c13dist_ax
        self.h1spec_ax  = h1spec_ax
        self.h1dist_ax  = h1dist_ax

        self.nmrproblem.peak_overlays_data = self.create1H13C1DSpectraOverlayData(
            self.nmrproblem
        )
        self.nmrproblem.spectra1D = self.create1H13C1DspectraData(self.nmrproblem)
        self.nmrproblem.distribution_data = self.createH1C13PlotDistributionsData(
            self.nmrproblem, 3
        )

        self.display1H13C1DmatplotlibSpectra(
            nmrproblem, [self.h1spec_ax, self.c13spec_ax]
        )
        (
            self.peak_overlays_dict,
            self.peak_overlays,
        ) = self.createH1C13matplotlibOverlaysPlot(
            self.nmrproblem, [self.h1spec_ax, self.c13spec_ax]
        )
        # self.display1H13C1Dspectra([self.c13spec_ax, self.h1spec_ax], nmrproblem)

        self.h1distdict, self.c13distdict = self.plotDistributions(
            self.nmrproblem, [self.h1dist_ax, self.c13dist_ax]
        )

        # self.c13distdict = self.plotC13Distributions(self.c13dist_ax, 3, self.nmrproblem)
        # self.h1distdict = self.plotH1Distributions(self.h1dist_ax, 3, self.nmrproblem)

        self.h1c13distlist = [self.h1distdict, self.c13distdict]

    def createH1C13PlotDistributionsData(self, nmrproblem, numCandidates):

        C13_ppm_axis = np.linspace(-30, 250, 500)
        H1_ppm_axis = np.linspace(-2, 16, 500)
        catoms = nmrproblem.carbonAtoms
        hatoms = nmrproblem.protonAtoms
        atoms = [hatoms, catoms]
        ppm_axis = [H1_ppm_axis, C13_ppm_axis]
        iprobs = nmrproblem.iprobs
        df = nmrproblem.df
        C13df = nmrproblem.udic[1]["df"]
        H1df = nmrproblem.udic[0]["df"]
        H1C13df = [H1df, C13df]

        distributions = {}
        for i in range(2):
            for k, ci in enumerate(atoms[i]):
                distributions[ci] = {}
                for j in iprobs[ci][:numCandidates]:

                    distributions[ci][j] = pd.DataFrame(
                        {
                            "xxx": ppm_axis[i],
                            "yyy": H1C13df[i].loc[j, "norm"].pdf(ppm_axis[i]),
                        }
                    )

                    distributions[ci][j]["label"] = H1C13df[i].loc[
                        j, "sF_latex_matplotlib"
                    ]
                    distributions[ci][j]["vline"] = float(df.loc["ppm", ci])

        return distributions

    def create1H13C1DspectraData(self, nmrproblem: nmrProblem.NMRproblem):

        spectra1D = {}

        h1c13 = ["proton1Dspectrum", "carbon1Dspectrum"]

        udic = nmrproblem.udic

        for i in range(2):
            xxx = udic[i]["axis"].ppm_scale()
            yyy = udic[i]["spec"]

            iii = ((np.roll(yyy, 1) - yyy) ** 2) > 1e-8
            iii[0] = 1.0
            iii[-1] = 1.0

            df = pd.DataFrame(
                data=np.array([xxx, yyy, iii]).transpose(),
                columns=["xxx", "yyy", "iii"],
            )

            df = df[df["iii"] == 1.0][["xxx", "yyy"]]

            spectra1D[h1c13[i]] = df

        return spectra1D

    def create1H13C1DSpectraOverlayData(self, nmrproblem):

        # w1 = widgets.Output()

        udic = nmrproblem.udic
        if "info" not in udic[0]:
            print("info not in udic[0]")
            return
        peak_overlays = []
        peak_overlays_data = {}

        # for proton and carbon spectra
        for i in range(udic["ndim"]):
            peak_overlays1 = []
            for Hi in udic[i]["info"].index:
                # print("Hi", Hi)

                il = int(udic[i]["info"].loc[Hi, "pk_left"])
                ir = int(udic[i]["info"].loc[Hi, "pk_right"])

                peak_overlays_data[Hi] = pd.DataFrame(
                    {
                        "xxx": udic[i]["axis"].ppm_scale()[il:ir],
                        "yyy": udic[i]["spec"][il:ir],
                    }
                )

        return peak_overlays_data

    def createH1C13interactivePlot(self, nmrproblem, h1c13distlist, ax0):

        # w1 = widgets.Output()

        udic = nmrproblem.udic
        if "info" not in udic[0]:
            print("info not in udic[0]")
            return
        peak_overlays = []
        peak_overlays_dict = {}

        # for proton and carbon spectra
        for i in range(udic["ndim"]):
            peak_overlays1 = []
            for Hi in udic[i]["info"].index:
                # print("Hi", Hi)

                il = int(udic[i]["info"].loc[Hi, "pk_left"])
                ir = int(udic[i]["info"].loc[Hi, "pk_right"])

                (pk,) = ax0[1 - i].plot(
                    udic[i]["axis"].ppm_scale()[il:ir],
                    udic[i]["spec"][il:ir],
                    lw=0.5,
                    c="black",
                    label=Hi,
                    gid=Hi,
                )

                peak_overlays1.append(pk)
                peak_overlays_dict[Hi] = pk

            peak_overlays.append(peak_overlays1)

        return peak_overlays_dict, peak_overlays

    def createH1C13matplotlibOverlaysPlot(self, nmrproblem, ax0):

        # w1 = widgets.Output()

        atoms = [nmrproblem.protonAtoms, nmrproblem.carbonAtoms]

        peak_overlays = []
        peak_overlays_dict = {}

        # for proton and carbon spectra
        for i, ax in enumerate(ax0):
            peak_overlays1 = []
            for atom_id in atoms[i]:

                (pk,) = ax.plot(
                    nmrproblem.peak_overlays_data[atom_id]["xxx"],
                    nmrproblem.peak_overlays_data[atom_id]["yyy"],
                    lw=0.5,
                    c="black",
                    label=atom_id,
                    gid=atom_id,
                )

                peak_overlays1.append(pk)
                peak_overlays_dict[atom_id] = pk

            peak_overlays.append(peak_overlays1)

        return peak_overlays_dict, peak_overlays

    def display1H13C1DmatplotlibSpectra(self, nmrproblem: nmrProblem.NMRproblem, ax0):

        xxx_labels = ["$^{1}$H [ppm]", "$^{13}$C [ppm]"]
        gid = ["h1ppm", "c13ppm"]
        spectra_id = list(nmrproblem.spectra1D.keys())

        for i, ax in enumerate(ax0):
            ax.plot(
                nmrproblem.spectra1D[spectra_id[i]]["xxx"],
                nmrproblem.spectra1D[spectra_id[i]]["yyy"],
                color="black",
                lw=0.5,
            )

            ax.set_xlim(ax.get_xlim()[::-1])

            ax.spines["top"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_xlabel(xxx_labels[i], fontsize=10, gid=gid[i])
            ax.set_yticks([])

    # def display1H13C1Dspectra(self, ax0, nmrproblem: nmrProblem.NMRproblem):
    #     """_summary_

    #     Args:
    #         ax (matplotlib.axes.Axes]): _description_
    #         nmrproblem (nmrProblem.NMRproblem): _description_
    #     """

    #     udic = nmrproblem.udic

    #     if (1 in udic) and (0 in udic):
    #         if (
    #             ("axis" in udic[1])
    #             and ("spec" in udic[1])
    #             and ("axis" in udic[0])
    #             and ("spec" in udic[0])
    #         ):
    #             xxx = udic[1]["axis"].ppm_scale()
    #             yyy = udic[1]["spec"]

    #             iii = ((np.roll(yyy, 1) - yyy) ** 2) > 1e-8
    #             iii[0] = 1.0
    #             iii[-1] = 1.0

    #             df = pd.DataFrame(
    #                 data=np.array([xxx, yyy, iii]).transpose(),
    #                 columns=["ppm", "yyy", "iii"],
    #             )

    #             full_spectrum1 = ax0[0].plot(
    #                 df[df["iii"] == 1.0]["ppm"],
    #                 df[df["iii"] == 1.0]["yyy"],
    #                 color="black",
    #                 lw=0.5,
    #             )
    #             ax0[0].set_xlim(udic[1]["axis"].ppm_limits())

    #             xxx = udic[0]["axis"].ppm_scale()
    #             yyy = udic[0]["spec"]

    #             iii = ((np.roll(yyy, 1) - yyy) ** 2) > 1e-8
    #             iii[0] = 1.0
    #             iii[-1] = 1.0

    #             df = pd.DataFrame(
    #                 data=np.array([xxx, yyy, iii]).transpose(),
    #                 columns=["ppm", "yyy", "iii"],
    #             )
    #             full_spectrum2 = ax0[1].plot(
    #                 df[df["iii"] == 1.0]["ppm"],
    #                 df[df["iii"] == 1.0]["yyy"],
    #                 color="black",
    #                 lw=0.5,
    #             )

    #             ax0[1].set_xlim(udic[0]["axis"].ppm_limits())

    #     # ax.set_xlim(10,-1)
    #     ax0[0].spines["top"].set_visible(False)
    #     ax0[0].spines["left"].set_visible(False)
    #     ax0[0].spines["right"].set_visible(False)
    #     ax0[0].set_xlabel("$^{13}$C [ppm]", fontsize=10, gid="c13ppm")
    #     ax0[0].set_yticks([])

    #     ax0[1].spines["top"].set_visible(False)
    #     ax0[1].spines["left"].set_visible(False)
    #     ax0[1].spines["right"].set_visible(False)
    #     ax0[1].set_xlabel("$^{1}$H [ppm]", fontsize=10, gid="h1ppm")
    #     ax0[1].set_yticks([])

    def plotDistributions(self, nmrproblem, ax0):

        atoms = [nmrproblem.protonAtoms, nmrproblem.carbonAtoms]
        xxx_labels = ["$^{1}$H [ppm]", "$^{13}$C [ppm]"]

        distdict = {}
        # for proton and carbon spectra
        for i, ax in enumerate(ax0):
            distdict[i] = {}
            ax.spines["top"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_xlabel(xxx_labels[i], fontsize=10, gid="ppm")
            ax.set_yticks([])

            for atom_id in atoms[i]:
                distlist = []

                for j, d in nmrproblem.distribution_data[atom_id].items():
                    if isinstance(j, (int, float)):
                        (distr,) = ax.plot(
                            d["xxx"], d["yyy"], "-", label=d.loc[0, "label"]
                        )

                        distr.set_visible(False)
                        distlist.append(distr)

                dline = ax.axvline(
                    nmrproblem.distribution_data[atom_id][j].loc[0, "vline"], c="r"
                )
                dline.set_visible(False)
                distlist.append(dline)

                distdict[i][atom_id] = distlist

            ax.set_xlim(ax.get_xlim()[::-1])

        return distdict[0], distdict[1]

        # plot top three candidates for each carbon present

        C13_ppm_axis = np.linspace(-30, 250, 500)
        catoms = nmrproblem.carbonAtoms
        iprobs = nmrproblem.iprobs
        df = nmrproblem.df
        C13df = nmrproblem.udic[1]["df"]

        c13distdict = {}

        for k, ci in enumerate(catoms):
            distlist = []
            for i in iprobs[ci][:numCandidates]:
                (c13distr,) = ax.plot(
                    C13_ppm_axis,
                    C13df.loc[i, "norm"].pdf(C13_ppm_axis),
                    label=C13df.loc[i, "sF_latex_matplotlib"],
                )

                c13distr.set_visible(False)
                distlist.append(c13distr)

            c13line = ax.axvline(float(df.loc["ppm", ci]))
            c13line.set_visible(False)
            distlist.append(c13line)

            c13distdict[ci] = distlist

        ax.set_xlabel("$^{13}$C [ppm]")
        ax.set_xlim(260, -40)
        return c13distdict

    def plotC13Distributions(self, ax, numCandidates, nmrproblem):

        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel("ppm", fontsize=10, gid="ppm")
        ax.set_yticks([])

        # plot top three candidates for each carbon present

        C13_ppm_axis = np.linspace(-30, 250, 500)
        catoms = nmrproblem.carbonAtoms
        iprobs = nmrproblem.iprobs
        df = nmrproblem.df
        C13df = nmrproblem.udic[1]["df"]

        c13distdict = {}

        for k, ci in enumerate(catoms):
            distlist = []
            for i in iprobs[ci][:numCandidates]:
                (c13distr,) = ax.plot(
                    C13_ppm_axis,
                    C13df.loc[i, "norm"].pdf(C13_ppm_axis),
                    label=C13df.loc[i, "sF_latex_matplotlib"],
                )

                c13distr.set_visible(False)
                distlist.append(c13distr)

            c13line = ax.axvline(float(df.loc["ppm", ci]))
            c13line.set_visible(False)
            distlist.append(c13line)

            c13distdict[ci] = distlist

        ax.set_xlabel("$^{13}$C [ppm]")
        ax.set_xlim(260, -40)
        return c13distdict

    def plotH1Distributions(self, ax, numCandidates, nmrproblem):

        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel("ppm", fontsize=10, gid="ppm")
        ax.set_yticks([])

        # plot top three candidates for each carbon present

        H1_ppm_axis = np.linspace(-5, 16, 500)
        patoms = nmrproblem.protonAtoms
        H1df = nmrproblem.udic[0]["df"]
        iprobs = nmrproblem.iprobs
        df = nmrproblem.df

        h1distdict = {}

        for k, hi in enumerate(patoms):
            distlist = []
            for i in iprobs[hi][:numCandidates]:
                (h1distr,) = ax.plot(
                    H1_ppm_axis,
                    H1df.loc[i, "norm"].pdf(H1_ppm_axis),
                    label=H1df.loc[i, "sF_latex_matplotlib"],
                )

                h1distr.set_visible(False)
                distlist.append(h1distr)

            h1line = ax.axvline(float(df.loc["ppm", hi]))
            h1line.set_visible(False)
            distlist.append(h1line)

            h1distdict[hi] = distlist

        ax.set_xlabel("$^{1}$H [ppm]")
        ax.set_xlim(12, -2)

        return h1distdict
