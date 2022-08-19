import sys
from PyQt5.QtWidgets import (
    QMainWindow,
    QApplication,
    QPushButton,
    QRadioButton,
    QVBoxLayout,
    QDialog,
    QDialogButtonBox,
)

class XY3dialog(QDialog):
    def __init__(self, parent=None, java_available=False):
        self.java_available = java_available
        super(XY3dialog, self).__init__(parent)
        self.setWindowTitle("Carbon Coordinates Methods")
        self.setFixedSize(300, 200)
        self.create_widgets()
        self.create_layout()
        self.create_connections()

    def create_widgets(self):
        # QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        QBtn = QDialogButtonBox.Ok 

        self.buttonBox = QDialogButtonBox(QBtn)


        self.use_xy3 = QRadioButton("Use saved carbon coordinates")
        self.use_xy3.setChecked(True)
        self.use_xy3.method = "xy3"
        
        self.use_c13ppm_predictor = QRadioButton("Use NMRshiftDB2 C13 PPM Predictor for coordinates")
        if not self.java_available:
            self.use_c13ppm_predictor.setCheckable(False)
            font = self.use_c13ppm_predictor.font()
            font.setItalic(True)
            font.setStrikeOut(True)
            self.use_c13ppm_predictor.setFont(font)
        self.use_c13ppm_predictor.method = "c13ppm"
        self.use_random_positions = QRadioButton("Use Random Positions for initial carbon coordinates")
        self.use_random_positions.method = "random"


    def create_layout(self):
        layout = QVBoxLayout()
        layout.addWidget(self.use_c13ppm_predictor)
        layout.addWidget(self.use_xy3)
        layout.addWidget(self.use_random_positions)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def create_connections(self):
        self.buttonBox.accepted.connect(self.accept)
        # self.buttonBox.rejected.connect(self.reject)

    def get_method(self):
        if self.use_c13ppm_predictor.isChecked():
            return self.use_c13ppm_predictor.method
        elif self.use_xy3.isChecked():
            return self.use_xy3.method
        elif self.use_random_positions.isChecked():
            return self.use_random_positions.method

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
