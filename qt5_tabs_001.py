import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget,QVBoxLayout, QTableWidget, QHeaderView, QTableWidgetItem
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWebEngineWidgets import QWebEngineView

import nmrProblem



class App(QMainWindow):

    def __init__(self):
        super().__init__()

        self.nmrproblem = nmrProblem.NMRproblem(nmrProblem.parse_argv())
        self.title = 'PyQt5 tabs - pythonspot.com'
        self.left = 100
        self.top = 100
        self.width = 900
        self.height = 600
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        
        self.table_widget = MyTabWidget(self, self.nmrproblem)
        # self.setCentralWidget(self.table_widget)

        layout = QVBoxLayout()
        layout.addWidget(self.table_widget)
        # layout.addWidget(sc)

        # Create a placeholder widget to hold our toolbar and canvas.
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)
        
        self.show()
    

class tablePane(QTableWidget):

    def __init__(self, parent, title, pd_df):
        super(QTableWidget, self).__init__(parent)
        self.df = pd_df

        # set table dimension
        nRows, nColumns = self.df.shape
        self.setColumnCount(nColumns)
        self.setRowCount(nRows)

        self.setHorizontalHeaderLabels(self.df.columns)
        self.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        # data insertion
        for i in range(self.rowCount()):
            for j in range(self.columnCount()):
                self.setItem(i, j, QTableWidgetItem(str(self.df.iloc[i, j])))

        parent.addTab(self, title)

class MyTabWidget(QTabWidget):
    
    def __init__(self, parent, nmrproblem):
        super(QTabWidget, self).__init__(parent)
        self.layout = QVBoxLayout()
        self.resize(300,200)
        
        # Initialize tab screen
        tab1 = tablePane(self, "H1", nmrproblem.h1)
        tab2 = tablePane(self, "C13", nmrproblem.c13)
        tab3 = tablePane(self, "COSY", nmrproblem.cosy)
        tab4 = tablePane(self, "HSQC", nmrproblem.hsqc)
        tab5 = tablePane(self, "HMBC", nmrproblem.hmbc)
        tab6 = tablePane(self, "Pureshift", nmrproblem.pureshift)
        tab7 = tablePane(self, "Summary", nmrproblem.df)
        
        
        # Add tabs
        # self.addTab(tab1,"Tab 1")
        # self.addTab(tab2,"Tab 2")
        
        # Create first tab
        # tab1.layout = QVBoxLayout(self)
        # pushButton1 = QPushButton("PyQt5 button")
        # tab1.layout.addWidget(pushButton1)
        # tab1.setLayout(tab1.layout)
        
        # Add tabs to widget
        # self.layout.addWidget(self.tabs)
        # self.setLayout(self.layout)
        
    # @pyqtSlot()
    # def on_click(self):
    #     print("Button Clicked\n")
    #     for currentQTabWidgetItem in self.tableWidget.selectedItems():
    #         print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
