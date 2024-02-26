
from PyQt5.QtWidgets import QMessageBox

def warning_dialog(return_message, title_message, qstarted=True):
    if qstarted:
        # create qt5 warning dialog with return_message
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Warning)
        msg_box.setText(return_message)
        msg_box.setStandardButtons(QMessageBox.Ok)
        msg_box.exec()
        msg_box.setWindowTitle(title_message)
    else:
        print(return_message)