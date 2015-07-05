
from __future__ import with_statement

import numpy as np
import sys

from PySide import QtCore, QtGui

from equalibrium import Ui_MainWindow


from oracle.transitions import AtomicTransition




class AtomicTransitionsTableModel(QtCore.QAbstractTableModel):

    def __init__(self, parent, data, header, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)

        self._data = data
        self._header = header

    def rowCount(self, parent):
        return len(self._data)

    def columnCount(self, parent):
        return len(self._header)

    def data(self, index, role):
        if not index.isValid(): return
        if role == QtCore.Qt.DisplayRole:
            if index.column() == 0: return
            return getattr(self._data[index.row()], self._header[index.column()])

        elif role == QtCore.Qt.CheckStateRole and index.column() == 0:
            return self._data[index.row()].acceptable

        return None



    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal \
        and role == QtCore.Qt.DisplayRole:
            return self._header[col]
        return None


    def sort(self, column_index, order):
        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        self._data = sorted(self._data,
            key=lambda t: getattr(t, self._header[column_index]))
        if order == QtCore.Qt.DescendingOrder:
            self._data.reverse()
        self.emit(QtCore.SIGNAL("layoutChanged()"))


    def flags(self, index):
        if not index.isValid(): return
        if index.column() == 0:
            return QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEditable|QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsUserCheckable


    def setData(self, index, value, role):
        if not index.isValid() or role != QtCore.Qt.CheckStateRole: return False
        self._data[index.row()].acceptable = value
        self.dataChanged.emit(index, index)
        return True


class DesignerMainWindow(QtGui.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):

        super(DesignerMainWindow, self).__init__(parent)
        self.setupUi(self)

        self._initialise_plots()

        data = [
            AtomicTransition(wavelength=0, species=26., e_low=10, log_gf=-1),
            AtomicTransition(wavelength=1, species=23., e_low=0, log_gf=-1),
            AtomicTransition(wavelength=2, species=2., e_low=4, log_gf=-1),
        ]
        data[0].acceptable = True
        data[1].acceptable = False
        data[2].acceptable = True


        table_model = AtomicTransitionsTableModel(self, data, ["acceptable", 
            "wavelength",
            "species", "e_low", "log_gf"])

        self.table_atomic_transitions.setModel(table_model)
        self.table_atomic_transitions.resizeColumnsToContents()
        self.table_atomic_transitions.setSortingEnabled(True)


        #self.table_atomic_transitions.cellChanged.connect(self._func)


    def _initialise_plots(self):
        self._initialise_plot_state()
        self._initialise_plot_atomic_transition()


    def _initialise_plot_state(self):
        """
        Initialise the plot state figure.
        """

        # Create the figures
        ax_excitation = self.plot_state.canvas.fig.add_subplot(211)
        ax_line_strength = self.plot_state.canvas.fig.add_subplot(212,
            sharey=ax_excitation)
        
        ax_excitation.axhline(0, c="k")
        ax_excitation.set_ylim(-1, 1)
        ax_excitation.set_xlim(0, 1)

        ax_line_strength.axhline(0, c="k")
        ax_line_strength.set_xlim(0, 1)

        # Labels.
        ax_excitation.set_xlabel("$E_{\\rm low}$ (eV)")
        ax_excitation.set_ylabel("[X/H]")

        ax_line_strength.set_xlabel("log(EW/$\lambda$)")
        ax_line_strength.set_ylabel("[X/H]")

        self.plot_state.canvas.fig.tight_layout()
        self.plot_state.canvas.draw()
    

    def _initialise_plot_atomic_transition(self):
        """
        Initialise the plot for the selected atomic transition.
        """

        ax = self.plot_atomic_transition.canvas.fig.add_subplot(111)

        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Flux")

        self.plot_atomic_transition.canvas.fig.tight_layout()
        self.plot_atomic_transition.canvas.draw()


    #def calculate_state(self, info=None):
    #    print("Calculating state with {}".format(info))


app = QtGui.QApplication(sys.argv)
dmw = DesignerMainWindow()
dmw.show()
sys.exit(app.exec_())
