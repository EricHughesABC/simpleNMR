import sys
from PyQt5.QtWidgets import (
    QMainWindow,
    QApplication,
    QPushButton,
    QRadioButton,
    QVBoxLayout,
    QHBoxLayout,
    QDialog,
    QDialogButtonBox,
    QLabel,
    QCheckBox,
)


class XY3dialog(QDialog):
    def __init__(self, parent=None, sheets_missing=[], java_available=False):
        self.java_available = java_available
        self.sheets_missing = sheets_missing
        super().__init__(parent)
        self.setWindowTitle("Solve Problem with ...")
        self.setFixedSize(400, 400)
        self.create_widgets()
        self.create_layout()
        self.create_connections()

    # def create_widgets(self):
    #     # QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
    #     QBtn = QDialogButtonBox.Ok

    #     self.buttonBox = QDialogButtonBox(QBtn)

    #     self.use_xy3 = QRadioButton("Use saved carbon coordinates")
    #     self.use_xy3.setChecked(True)
    #     self.use_xy3.method = "xy3"

    #     self.use_c13ppm_predictor = QRadioButton(
    #         "Use NMRshiftDB2 C13 PPM Predictor for coordinates"
    #     )
    #     if not self.java_available:
    #         self.use_c13ppm_predictor.setCheckable(False)
    #         font = self.use_c13ppm_predictor.font()
    #         font.setItalic(True)
    #         font.setStrikeOut(True)
    #         self.use_c13ppm_predictor.setFont(font)
    #     self.use_c13ppm_predictor.method = "c13ppm"
    #     self.use_random_positions = QRadioButton(
    #         "Use Random Positions for initial carbon coordinates"
    #     )
    #     self.use_random_positions.method = "random"

    def create_widgets(self):
        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel

        self.buttonBox = QDialogButtonBox(QBtn)

        self.checkboxes = {}

        self.checkboxes["H1_1D"] = QCheckBox("H1 1D")
        self.checkboxes["C13_1D"] = QCheckBox("C13 1D")
        self.checkboxes["H1_pureshift"] = QCheckBox("H1 Pure Shift")
        self.checkboxes["HSQC"] = QCheckBox("HSQC")
        self.checkboxes["HMBC"] = QCheckBox("HMBC")
        self.checkboxes["COSY"] = QCheckBox("COSY")
        self.checkboxes["NOESY"] = QCheckBox("NOESY")

        # create an error  message text widget
        self.error_message = QLabel("OK")
        self.error_message.setStyleSheet("color: green")

        # set check all the boxes as default
        for key in self.checkboxes.keys():
            self.checkboxes[key].setChecked(True)

        # disable the boxes where the sheets are missing
        print("sheets missing: ", self.sheets_missing)
        
        for key in self.sheets_missing:
            self.checkboxes[key].setEnabled(False)
            self.checkboxes[key].setChecked(False)

        # disable and uncheck NOESY by default
        self.checkboxes["NOESY"].setChecked(False)
        self.checkboxes["NOESY"].setEnabled(False)

        # Display warning if HSQC is one of the missing sheets
        if "HSQC" in self.sheets_missing:
            self.checkboxes["HSQC"].setText(
                "WARNING: HSQC data is missing! Program requires HSQC data to run."
            )
            # set the HSQC checkbox text to red
            self.checkboxes["HSQC"].setStyleSheet("color: red")


            # set the okay button to disabled
            self.buttonBox.button(QDialogButtonBox.Ok).setEnabled(False)

        # do not allow editing of HSQC, COSY and HMBC by default
        self.checkboxes["HSQC"].setEnabled(False)

        # add radio buttons for calculating xy3 coordinates of the molecule

        self.use_xy3 = QRadioButton("Use saved carbon coordinates")

        self.use_xy3.method = "xy3"

        self.use_c13ppm_predictor = QRadioButton(
            "Use NMRshiftDB2 C13 PPM Predictor for coordinates"
        )
        if self.java_available:
            # set the radio button to checked
            self.use_c13ppm_predictor.setChecked(True)
        else:
            self.use_c13ppm_predictor.setCheckable(False)
            font = self.use_c13ppm_predictor.font()
            font.setItalic(True)
            font.setStrikeOut(True)
            self.use_c13ppm_predictor.setFont(font)
            self.use_xy3.setChecked(True)

        self.use_c13ppm_predictor.method = "c13ppm"
        self.use_random_positions = QRadioButton(
            "Use Random Positions for initial carbon coordinates"
        )
        self.use_random_positions.method = "random"

    # def create_layout(self):
    #     layout = QVBoxLayout()
    #     layout.addWidget(self.use_c13ppm_predictor)
    #     layout.addWidget(self.use_xy3)
    #     layout.addWidget(self.use_random_positions)
    #     layout.addWidget(self.buttonBox)
    #     self.setLayout(layout)

    def create_layout(self):

        vlayout = QVBoxLayout()
        layout_H1 = QHBoxLayout()
        layout_H2 = QHBoxLayout()
        layout_H3 = QHBoxLayout()
        for key in ["H1_1D", "C13_1D", "H1_pureshift"]:
            layout_H1.addWidget(self.checkboxes[key])
        for key in ["COSY", "HMBC", "NOESY"]:
            layout_H2.addWidget(self.checkboxes[key])
        for key in ["HSQC"]:
            layout_H3.addWidget(self.checkboxes[key])

        # add an explanatory text for the checkboxes
        vlayout.addWidget(QLabel("Select NMR spectra to use for solving the problem:"))

        # add horizontal layouts to vertical layout
        vlayout.addLayout(layout_H1)
        vlayout.addLayout(layout_H2)
        vlayout.addLayout(layout_H3)

        # add radio buttons to vertical layout
        # add a surrounding line around the radio buttons
        vlayout.addWidget(QLabel("Select method for calculating carbon coordinates:"))

        vlayout.addWidget(self.use_xy3)
        vlayout.addWidget(self.use_c13ppm_predictor)
        vlayout.addWidget(self.use_random_positions)

        # vlayout.addWidget(self.error_message)
        vlayout.addWidget(self.buttonBox)

        self.setLayout(vlayout)

    # def create_connections(self):
    #     self.buttonBox.accepted.connect(self.accept)
    #     # self.buttonBox.rejected.connect(self.reject)
    def create_connections(self):
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

    # def get_method(self):

    #     if self.use_c13ppm_predictor.isChecked():
    #         return self.use_c13ppm_predictor.method
    #     elif self.use_xy3.isChecked():
    #         return self.use_xy3.method
    #     elif self.use_random_positions.isChecked():
    #         return self.use_random_positions.method
    def get_method(self):
        kys = []

        for key in self.checkboxes.keys():
            if self.checkboxes[key].isChecked():
                kys.append(key)

        method = None
        if self.use_c13ppm_predictor.isChecked():
            method = self.use_c13ppm_predictor.method
        elif self.use_xy3.isChecked():
            method = self.use_xy3.method
        elif self.use_random_positions.isChecked():
            method = self.use_random_positions.method

        return method, kys


if __name__ == "__main__":

    class MainWindow(QMainWindow):
        def __init__(self):
            super().__init__()

            self.setWindowTitle("My App")

            button = QPushButton("Press me for a dialog!")
            button.clicked.connect(self.button_clicked)
            self.setCentralWidget(button)

        def button_clicked(self, s):
            # print("click", s)

            dlg = XY3dialog(self)

            # dlg.setWindowTitle("HELLO!")

            # table_widget = MyTabWidget(dlg, self.nmrproblem)
            if dlg.exec():
                print("Success!", dlg.get_method())
                # print(type(dlg.table_widget))
            else:
                print("Cancel!")

    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec()
