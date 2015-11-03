# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'imContrastDialog1.ui'
#
# Created: Fri Jul 31 10:36:07 2015
#      by: PyQt4 UI code generator 4.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(425, 438)
        Dialog.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.gridLayout = QtGui.QGridLayout(Dialog)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.buttonBoxDiag = QtGui.QDialogButtonBox(Dialog)
        self.buttonBoxDiag.setMaximumSize(QtCore.QSize(120, 16777215))
        self.buttonBoxDiag.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBoxDiag.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBoxDiag.setObjectName(_fromUtf8("buttonBoxDiag"))
        self.gridLayout.addWidget(self.buttonBoxDiag, 2, 4, 1, 1)
        self.stackedWidget = QtGui.QStackedWidget(Dialog)
        self.stackedWidget.setObjectName(_fromUtf8("stackedWidget"))
        self.stackedWidgetPage1 = QtGui.QWidget()
        self.stackedWidgetPage1.setObjectName(_fromUtf8("stackedWidgetPage1"))
        self.gridLayout_2 = QtGui.QGridLayout(self.stackedWidgetPage1)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.hSliderContrast = QtGui.QSlider(self.stackedWidgetPage1)
        self.hSliderContrast.setProperty("value", 49)
        self.hSliderContrast.setOrientation(QtCore.Qt.Horizontal)
        self.hSliderContrast.setObjectName(_fromUtf8("hSliderContrast"))
        self.gridLayout_2.addWidget(self.hSliderContrast, 5, 1, 1, 6)
        self.label = QtGui.QLabel(self.stackedWidgetPage1)
        self.label.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout_2.addWidget(self.label, 4, 0, 1, 1)
        self.frame = QtGui.QFrame(self.stackedWidgetPage1)
        self.frame.setFrameShape(QtGui.QFrame.NoFrame)
        self.frame.setFrameShadow(QtGui.QFrame.Plain)
        self.frame.setObjectName(_fromUtf8("frame"))
        self.gridLayout_2.addWidget(self.frame, 3, 0, 1, 7)
        self.hSliderIntensity = QtGui.QSlider(self.stackedWidgetPage1)
        self.hSliderIntensity.setProperty("value", 49)
        self.hSliderIntensity.setOrientation(QtCore.Qt.Horizontal)
        self.hSliderIntensity.setObjectName(_fromUtf8("hSliderIntensity"))
        self.gridLayout_2.addWidget(self.hSliderIntensity, 4, 1, 1, 6)
        self.doubleSpinBox_range = QtGui.QDoubleSpinBox(self.stackedWidgetPage1)
        self.doubleSpinBox_range.setMinimumSize(QtCore.QSize(100, 0))
        self.doubleSpinBox_range.setMaximumSize(QtCore.QSize(100, 16777215))
        self.doubleSpinBox_range.setPrefix(_fromUtf8(""))
        self.doubleSpinBox_range.setObjectName(_fromUtf8("doubleSpinBox_range"))
        self.gridLayout_2.addWidget(self.doubleSpinBox_range, 9, 4, 1, 1)
        self.label_3 = QtGui.QLabel(self.stackedWidgetPage1)
        self.label_3.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_2.addWidget(self.label_3, 7, 0, 1, 1)
        self.label_2 = QtGui.QLabel(self.stackedWidgetPage1)
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_2.addWidget(self.label_2, 5, 0, 1, 1)
        self.doubleSpinBox_min = QtGui.QDoubleSpinBox(self.stackedWidgetPage1)
        self.doubleSpinBox_min.setMinimumSize(QtCore.QSize(100, 0))
        self.doubleSpinBox_min.setMaximumSize(QtCore.QSize(100, 16777215))
        self.doubleSpinBox_min.setPrefix(_fromUtf8(""))
        self.doubleSpinBox_min.setObjectName(_fromUtf8("doubleSpinBox_min"))
        self.gridLayout_2.addWidget(self.doubleSpinBox_min, 7, 1, 1, 1)
        self.doubleSpinBox_max = QtGui.QDoubleSpinBox(self.stackedWidgetPage1)
        self.doubleSpinBox_max.setMinimumSize(QtCore.QSize(100, 0))
        self.doubleSpinBox_max.setMaximumSize(QtCore.QSize(100, 16777215))
        self.doubleSpinBox_max.setPrefix(_fromUtf8(""))
        self.doubleSpinBox_max.setObjectName(_fromUtf8("doubleSpinBox_max"))
        self.gridLayout_2.addWidget(self.doubleSpinBox_max, 9, 1, 1, 1)
        self.label_5 = QtGui.QLabel(self.stackedWidgetPage1)
        self.label_5.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_2.addWidget(self.label_5, 9, 0, 1, 1)
        self.label_4 = QtGui.QLabel(self.stackedWidgetPage1)
        self.label_4.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_2.addWidget(self.label_4, 7, 3, 1, 1)
        self.label_6 = QtGui.QLabel(self.stackedWidgetPage1)
        self.label_6.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout_2.addWidget(self.label_6, 9, 3, 1, 1)
        self.doubleSpinBox_mid = QtGui.QDoubleSpinBox(self.stackedWidgetPage1)
        self.doubleSpinBox_mid.setMinimumSize(QtCore.QSize(100, 0))
        self.doubleSpinBox_mid.setMaximumSize(QtCore.QSize(100, 16777215))
        self.doubleSpinBox_mid.setPrefix(_fromUtf8(""))
        self.doubleSpinBox_mid.setObjectName(_fromUtf8("doubleSpinBox_mid"))
        self.gridLayout_2.addWidget(self.doubleSpinBox_mid, 7, 4, 1, 1)
        self.radio_right = QtGui.QRadioButton(self.stackedWidgetPage1)
        self.radio_right.setMaximumSize(QtCore.QSize(60, 16777215))
        self.radio_right.setObjectName(_fromUtf8("radio_right"))
        self.gridLayout_2.addWidget(self.radio_right, 9, 6, 1, 1)
        self.radio_left = QtGui.QRadioButton(self.stackedWidgetPage1)
        self.radio_left.setMaximumSize(QtCore.QSize(60, 16777215))
        self.radio_left.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.radio_left.setChecked(True)
        self.radio_left.setObjectName(_fromUtf8("radio_left"))
        self.gridLayout_2.addWidget(self.radio_left, 7, 6, 1, 1)
        self.stackedWidget.addWidget(self.stackedWidgetPage1)
        self.stackedWidgetPage2 = QtGui.QWidget()
        self.stackedWidgetPage2.setObjectName(_fromUtf8("stackedWidgetPage2"))
        self.stackedWidget.addWidget(self.stackedWidgetPage2)
        self.gridLayout.addWidget(self.stackedWidget, 0, 1, 1, 4)
        self.pushButtonAuto = QtGui.QPushButton(Dialog)
        self.pushButtonAuto.setMaximumSize(QtCore.QSize(50, 16777215))
        self.pushButtonAuto.setObjectName(_fromUtf8("pushButtonAuto"))
        self.gridLayout.addWidget(self.pushButtonAuto, 2, 2, 1, 1)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 2, 3, 1, 1)

        self.retranslateUi(Dialog)
        self.stackedWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label.setText(_translate("Dialog", "Intensity", None))
        self.label_3.setText(_translate("Dialog", "min", None))
        self.label_2.setText(_translate("Dialog", "Contrast", None))
        self.label_5.setText(_translate("Dialog", "max", None))
        self.label_4.setText(_translate("Dialog", "mid", None))
        self.label_6.setText(_translate("Dialog", "range", None))
        self.radio_right.setText(_translate("Dialog", "Right", None))
        self.radio_left.setText(_translate("Dialog", "Left", None))
        self.pushButtonAuto.setText(_translate("Dialog", "Auto", None))

