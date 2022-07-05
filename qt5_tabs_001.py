import sys
from PyQt5.QtWidgets import (QMainWindow, 
                                QApplication, 
                                QPushButton, 
                                QWidget, 
                                QAction, 
                                QTabWidget,
                                QVBoxLayout, 
                                QTableWidget, 
                                QHeaderView, 
                                QTableWidgetItem,
                                QTableView,
                                QDialog,
                                QDialogButtonBox,
                                QLabel)

from PyQt5.Qt import QApplication, QClipboard

from PyQt5.QtGui import QIcon, QKeySequence
from PyQt5.QtCore import pyqtSlot, QAbstractTableModel, Qt
from PyQt5.QtWebEngineWidgets import QWebEngineView

import pandas as pd
import io

import nmrProblem

from excelheaders import excel_orig_df_columns

# excel_orig_df_columns = {
#     "molecule": [
#         "moleclule", 
#         "smile"
#         ],
#     "H1_1D": [
#         "Name",
#         "Shift",
#         "Range",
#         "H's",
#         "Integral",
#         "Class",
#         "J's",
#         "Method"
#     ],
#     "C13_1D": [
#         "ppm",
#         "Intensity",
#         "Width",
#         "Area",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation"
#     ],
#     "H1_pureshift": [
#         "ppm",
#         "Intensity",
#         "Width",
#         "Area",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation"
#     ],
#     "COSY": [
#         "f2 (ppm)",
#         "f1 (ppm)",
#         "Intensity",
#         "Width f2",
#         "Width f1",
#         "Volume",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation"
#     ],
#     "HSQC": [
#         "f2 (ppm)",
#         "f1 (ppm)",
#         "Intensity",
#         "Width f2",
#         "Width f1",
#         "Volume",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation"
#     ],
#     "HMBC": [        
#         "f2 (ppm)",
#         "f1 (ppm)",
#         "Intensity",
#         "Width f2",
#         "Width f1",
#         "Volume",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation"
#     ]
# }


class TableModel(QAbstractTableModel):

    def __init__(self, dataframe: pd.DataFrame):
        super(TableModel, self).__init__()
        self._df = dataframe

    def data(self, index, role):
        if role == Qt.DisplayRole:
            value = self._df.iloc[index.row(), index.column()]
            return str(value)

    def rowCount(self, index):
        return self._df.shape[0]

    def columnCount(self, index):
        return self._df.shape[1]

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._df.columns[section])

            if orientation == Qt.Vertical:
                return str(self._df.index[section])

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole or role == Qt.EditRole:
                return str(self._df.iloc[index.row(), index.column()])
        return None

    def setData(self, index, value, role):
        if role == Qt.EditRole:
            self._df.iloc[index.row(), index.column()] = value
            print(self._df)
            return True

    def flags(self, index):
        return Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsEditable

   

class tablePane(QTableView):

    def __init__(self, parent, tab_title: str, pd_df: pd.DataFrame):
        # super(QTableView, self).__init__(parent)
        super().__init__(parent)
        self.df = pd_df

        self.model = TableModel(self.df)
        self.setModel(self.model)

        self.setAlternatingRowColors(True)
        parent.addTab(self, tab_title)


    def keyPressEvent(self, event):
    #     # clipboard = QApplication.clipboard()
    #     # if event.matches(QKeySequence.Copy):
    #     #     print('Ctrl + C')
    #     #     clipboard.setText("some text")
        if event.matches(QKeySequence.Paste):
            print(QApplication.clipboard().text())
            print('Ctrl + V')

            ctext = QApplication.clipboard().text()
            if "\t" in ctext:
                print("tabs in ctext")
                if ctext.endswith("\t"):
                    print("remove extra tab")
                    ctext = ctext.rstrip("\t")
                else:
                    print("No extra tab at end")
                print(ctext)
                iotext = io.StringIO(ctext)
                # self.b.insertPlainText(text + '\n')
                print("\nPandas dataframe\n")

                self.df = pd.read_csv(iotext, sep="\t", index_col=0, dtype=str)
                print(self.df)

                self.model = TableModel(self.df)
                self.setModel(self.model)
                self.setAlternatingRowColors(True)


        QTableView.keyPressEvent(self, event)
        


class MyTabWidget(QTabWidget):
    
    def __init__(self, parent, tabTitles_dataframes):
        super(QTabWidget, self).__init__(parent)
        self.layout = QVBoxLayout()
        self.resize(300,200)

        self.tables = {}
        
        for title, table in tabTitles_dataframes.items():
            self.tables[title] = tablePane(self, title, table)


        
    # @pyqtSlot()
    # def on_click(self):
    #     print("Button Clicked\n")
    #     for currentQTabWidgetItem in self.tableWidget.selectedItems():
    #         print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())

class EditDataFrameDialog(QDialog):
    def __init__(self, nmrproblem):
        super().__init__()

        self.nmrproblem = nmrproblem



        # tabtitles_dataframes = {"Molecule": nmrproblem.molecule_df,
        #                         "H1": nmrproblem.h1_df,
        #                         "C13": nmrproblem.c13_df}

        tabtitles_dataframes = {}
        for t, v in excel_orig_df_columns.items():
            tabtitles_dataframes[t] = pd.DataFrame( columns=v)
            if t == "molecule":
                tabtitles_dataframes[t].loc[1] = [""] * len(v)

        self.table_widget = MyTabWidget(self, tabtitles_dataframes)

        self.setWindowTitle("MestreNova - Edit DataFrame")
        self.setGeometry(100, 100, 600, 400)

        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel

        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.on_accept)
        self.buttonBox.rejected.connect(self.reject)

        self.layout = QVBoxLayout()
        message = QLabel("Something happened, is that OK?")
        self.layout.addWidget(self.table_widget)
        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)

    def on_accept(self):
        print(list(self.table_widget.tables.keys()))

        # copy pandas df tables to nmrproblem

        for k in excel_orig_df_columns.keys():
            print(k)

            nmrProblem.new_dataframes[k] = self.table_widget.tables[k].df.copy()



        self.accept()

if __name__ == '__main__':

    class MainWindow(QMainWindow):
        def __init__(self):
            super().__init__()

            self.setWindowTitle("My App")

            self.nmrproblem = nmrProblem.NMRproblem(nmrProblem.parse_argv())

            button = QPushButton("Press me for a dialog!")
            button.clicked.connect(self.button_clicked)
            self.setCentralWidget(button)

            

        def button_clicked(self, s):
            print("click", s)

            dlg = EditDataFrameDialog(self.nmrproblem)

            # dlg.setWindowTitle("HELLO!")

            # table_widget = MyTabWidget(dlg, self.nmrproblem)
            if dlg.exec():
                print("Success!")
                # print(type(dlg.table_widget))
            else:
                print("Cancel!")


    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    app.exec()
