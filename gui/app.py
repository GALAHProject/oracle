
from __future__ import with_statement

import numpy as np
import sys

from PySide import QtCore, QtGui

from equalibria import Ui_MainWindow

class DesignerMainWindow(QtGui.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):

        super(DesignerMainWindow, self).__init__(parent)
        self.setupUi(self)
        foo = self.tabWidget.setContentsMargins(0,0,0,0)
        print(foo)

    def calculate_state(self, info=None):
        print("Calculating state with {}".format(info))


app = QtGui.QApplication(sys.argv)
dmw = DesignerMainWindow()
dmw.show()
sys.exit(app.exec_())
