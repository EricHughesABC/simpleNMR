import sys
from PyQt5.QtWidgets import (
    QMainWindow,
    QApplication,
    QPushButton,
    QVBoxLayout,
    QDialog,
    QLabel,
    QDialogButtonBox,
)


class Aboutdialog(QDialog):
    def __init__(self, parent=None):
        super(Aboutdialog, self).__init__(parent)
        self.setWindowTitle("About simpleNMR")
        self.setFixedSize(300, 200)
        self.create_widgets()
        # self.create_layout()
        self.create_connections()

    def create_widgets(self):
        # QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        QBtn = QDialogButtonBox.Ok

        # add message about the program
        self.message = QLabel(
            "simpleNMR version 0.0.15  08March2024\n\nWritten and designed by\n\n Eric Hughes and Alan Kenwright\n\nMIT License"
        )
        self.message.setWordWrap(True)

        self.buttonBox = QDialogButtonBox(QBtn)

        # create vertical layout box and add widgets
        layout = QVBoxLayout()
        layout.addWidget(self.message)
        layout.addWidget(self.buttonBox)

        # set layout
        self.setLayout(layout)

    def create_connections(self):
        self.buttonBox.accepted.connect(self.accept)


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

            dlg = Aboutdialog(self)

            # dlg.setWindowTitle("HELLO!")

            # table_widget = MyTabWidget(dlg, self.nmrproblem)
            dlg.exec()

    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec()
