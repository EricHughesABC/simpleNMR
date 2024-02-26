import sys
import pathlib
from PyQt5.QtWidgets import (
    QMainWindow,
    QApplication,
    QPushButton,
    QVBoxLayout,
    QDialog,
    QDialogButtonBox,
    QLabel,
    QLineEdit,
)


import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


class SmilesDialog(QDialog):
    def __init__(self, smilesStr="", smilesPNG=None):
        self.smilesStr = smilesStr
        super().__init__()
        self.setWindowTitle("Use old smiles string ...")
        self.setFixedSize(350, 300)
        self.create_widgets(self.smilesStr, smilesPNG)
        self.create_layout()
        self.create_connections()

    def create_widgets(self, smilesStr, smilesPNG):
        # create a dialog box with a text entry field
        # a widget to dsisplay a png image
        # a button to accept the entry or cancel

        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel

        self.buttonBox = QDialogButtonBox(QBtn)

        self.smiles_entry = QLineEdit()
        self.smiles_entry.setText(smilesStr)
        self.smiles_entry.setPlaceholderText("Enter SMILES string")
        self.smiles_entry.setFixedWidth(300)
        self.smiles_png = QLabel()
        if smilesPNG:
            # self.smiles_png.setPixmap(smilesPNG.toqpixmap())
            self.smiles_png.setPixmap(smilesPNG.toqpixmap().scaled(300, 200))

    def create_layout(self):
        # layout = QVBoxLayout()
        # layout.addWidget(self.use_xy3)
        # layout.addWidget(self.use_c13ppm_predictor)
        # layout.addWidget(self.use_random_positions)
        # layout.addWidget(self.buttonBox)
        # self.setLayout(layout)

        layout = QVBoxLayout()
        layout.addWidget(self.smiles_entry)
        layout.addWidget(self.smiles_png)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def create_connections(self):
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

    def get_smiles_string(self):
        return self.smiles_entry.text()


if __name__ == "__main__":

    import numpy as np
    import PIL

    from PyQt5.QtWidgets import QMainWindow, QApplication

    # from PyQt5.QtGui import QPixmap, QImage

    class MainWindow(QMainWindow):
        def __init__(self):
            super().__init__()

            self.setWindowTitle("My App")

            button = QPushButton("Press me for a dialog!")
            button.clicked.connect(self.button_clicked)
            self.setCentralWidget(button)

        def button_clicked(self):
            # print("click", s)

            smilesStr = "C1=CC=C(C=C1)C(=O)O"
            mol = Chem.MolFromSmiles(smilesStr)
            smilesPNG = Draw.MolToImage(mol)
            dlg = SmilesDialog(smilesStr, smilesPNG)

            # dlg.setWindowTitle("HELLO!")

            # table_widget = MyTabWidget(dlg, self.nmrproblem)
            if dlg.exec():
                print("Success!", dlg.get_smiles_string())
                # print(type(dlg.table_widget))
            else:
                print("Cancel!")

    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec()
