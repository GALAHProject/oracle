# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'equalibrium.ui'
#
# Created: Sat Jul  4 23:50:55 2015
#      by: pyside-uic 0.2.14 running on PySide 1.2.1
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui
from mplwidget import MPLWidget

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 480)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayoutWidget_2 = QtGui.QWidget(self.centralwidget)
        self.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(0, 0, 872, 341))
        self.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.horizontalLayoutWidget_2)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label_3 = QtGui.QLabel(self.horizontalLayoutWidget_2)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 2, 0, 1, 1)
        self.label_4 = QtGui.QLabel(self.horizontalLayoutWidget_2)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 3, 0, 1, 1)
        self.label = QtGui.QLabel(self.horizontalLayoutWidget_2)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label_2 = QtGui.QLabel(self.horizontalLayoutWidget_2)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.sp_teff = QtGui.QSpinBox(self.horizontalLayoutWidget_2)
        self.sp_teff.setMinimum(3500)
        self.sp_teff.setMaximum(8000)
        self.sp_teff.setSingleStep(50)
        self.sp_teff.setProperty("value", 5500)
        self.sp_teff.setObjectName("sp_teff")
        self.gridLayout.addWidget(self.sp_teff, 0, 1, 1, 1)
        self.sp_logg = QtGui.QDoubleSpinBox(self.horizontalLayoutWidget_2)
        self.sp_logg.setDecimals(3)
        self.sp_logg.setMaximum(5.5)
        self.sp_logg.setSingleStep(0.1)
        self.sp_logg.setProperty("value", 4.5)
        self.sp_logg.setObjectName("sp_logg")
        self.gridLayout.addWidget(self.sp_logg, 1, 1, 1, 1)
        self.sp_metallicity = QtGui.QDoubleSpinBox(self.horizontalLayoutWidget_2)
        self.sp_metallicity.setDecimals(3)
        self.sp_metallicity.setMinimum(-5.0)
        self.sp_metallicity.setMaximum(0.5)
        self.sp_metallicity.setSingleStep(0.1)
        self.sp_metallicity.setObjectName("sp_metallicity")
        self.gridLayout.addWidget(self.sp_metallicity, 2, 1, 1, 1)
        self.sp_xi = QtGui.QDoubleSpinBox(self.horizontalLayoutWidget_2)
        self.sp_xi.setDecimals(3)
        self.sp_xi.setMaximum(5.0)
        self.sp_xi.setSingleStep(0.1)
        self.sp_xi.setProperty("value", 1.0)
        self.sp_xi.setObjectName("sp_xi")
        self.gridLayout.addWidget(self.sp_xi, 3, 1, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.btn_edit_model = QtGui.QPushButton(self.horizontalLayoutWidget_2)
        self.btn_edit_model.setObjectName("btn_edit_model")
        self.verticalLayout.addWidget(self.btn_edit_model)
        self.btn_calculate_state = QtGui.QPushButton(self.horizontalLayoutWidget_2)
        self.btn_calculate_state.setObjectName("btn_calculate_state")
        self.verticalLayout.addWidget(self.btn_calculate_state)
        self.line = QtGui.QFrame(self.horizontalLayoutWidget_2)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout.addWidget(self.line)
        self.btn_solve_stellar_parameters = QtGui.QPushButton(self.horizontalLayoutWidget_2)
        self.btn_solve_stellar_parameters.setObjectName("btn_solve_stellar_parameters")
        self.verticalLayout.addWidget(self.btn_solve_stellar_parameters)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.plot_state = MPLWidget(self.horizontalLayoutWidget_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plot_state.sizePolicy().hasHeightForWidth())
        self.plot_state.setSizePolicy(sizePolicy)
        self.plot_state.setMinimumSize(QtCore.QSize(400, 0))
        self.plot_state.setObjectName("plot_state")
        self.horizontalLayout.addWidget(self.plot_state)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.plot_atomic_transition = MPLWidget(self.horizontalLayoutWidget_2)
        self.plot_atomic_transition.setMinimumSize(QtCore.QSize(140, 140))
        self.plot_atomic_transition.setObjectName("plot_atomic_transition")
        self.verticalLayout_2.addWidget(self.plot_atomic_transition)
        self.table_atomic_transitions = QtGui.QTableView(self.horizontalLayoutWidget_2)
        self.table_atomic_transitions.setObjectName("table_atomic_transitions")
        self.verticalLayout_2.addWidget(self.table_atomic_transitions)
        self.horizontalLayout.addLayout(self.verticalLayout_2)
        self.horizontalLayout_2.addLayout(self.horizontalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar()
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "Metallicity", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("MainWindow", "Microturbulence (km/s)", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Effective temperature (K)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "Surface gravity", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_edit_model.setText(QtGui.QApplication.translate("MainWindow", "Edit model configuration..", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_calculate_state.setText(QtGui.QApplication.translate("MainWindow", "Calculate state", None, QtGui.QApplication.UnicodeUTF8))
        self.btn_solve_stellar_parameters.setText(QtGui.QApplication.translate("MainWindow", "Solve stellar parameters", None, QtGui.QApplication.UnicodeUTF8))
