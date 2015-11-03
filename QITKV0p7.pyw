# -*- coding: utf-8 -*-
"""
Created on Thu Apr 02 18:01:03 2015

@author: liw5
"""


from __future__ import print_function
import sys
#import os
from PyQt4 import QtCore, QtGui
#from PyQt4.QtGui import QApplication, QCursor
#from PyQt4.QtCore import *
#from PyQt4.QtGui import *

#import SimpleITK as sitk
#import vtk
#from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
#from scipy import ndimage
import numpy as np
#import random
#from numpy import arange, sin, pi
#import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.widgets import Cursor
#from matplotlib.backends import qt4_compat
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
#from matplotlib.backend_bases import key_press_handler

import nibabel as nib
from unwrap import unwrap
from skimage.restoration import unwrap_phase

import time
import LaplacianUnwrap as LPuwp
#from numpy.random import normal
    
from QITKV0p9 import Ui_MainWindow
from imContrastDialog2 import Ui_Dialog


sys.stderr = sys.stdout


if sys.version_info[:2] < (2, 5):
    def partial(func, arg):
        def callme():
            return func(arg)
        return callme
else:
    pass    
    #from functools import partial

class MyDialog(QtGui.QDialog,Ui_Dialog):
    _window_list = []
    def __init__(self, parent=None):
        super(MyDialog, self).__init__(parent)
        MyDialog._window_list.append(self)
        self.setupUi(self)
    
class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
    _window_list = []
    def __init__(self):
        super(MainWindow, self).__init__()
        MainWindow._window_list.append(self)
        self.setupUi(self)
        self.statusBar().showMessage('Ready' )
        self.stackedWidget_image.setContentsMargins(0,0,0,0)
        self.treeWidget_parameters.setContentsMargins(0,0,0,0)
        self.stackedWidget_control.setContentsMargins(0,0,0,0)
        self.menuFile.addAction('&Exit',self.fileQuit)
        self.menuHelp.addAction('&About', self.about)
        self.menuOptions.addAction('&ROI',self.on_ROI_view)
        self.toolBar.setToolButtonStyle(QtCore.Qt.ToolButtonIconOnly)
        self.actionOpen.triggered.connect(self.on_actionOpen)
        self.actionSave.triggered.connect(self.on_actionSave)
        self.actionZoom_in.triggered.connect(self.on_Zoom_in)
        self.actionZoom_out.triggered.connect(self.on_Zoomout)
        self.actionAdjust_Contrast.triggered.connect(self.on_AdjustContrast)
        self.imagepage1_panelselect=1
        self.horizontalScrollBar_z_2panel.valueChanged[int].connect(self.VScrollBarChangeValue)
        self.spinBox_x_3panel.valueChanged[str].connect(self.onSpinBoxX_Changed)      
        self.spinBox_y_3panel.valueChanged[str].connect(self.onSpinBoxY_Changed)      
        self.spinBox_z_3panel.valueChanged[str].connect(self.onSpinBoxZ_Changed)   
        self.spinBox_v_3panel.valueChanged[str].connect(self.onSpinBoxV_Changed)   
        self.spinBox_v_2panel.valueChanged[str].connect(self.onspinBox_v_2panel_Changed)   

        self.pushButton_opened.clicked.connect(self.on_Opened)
        self.treeWidget.clicked.connect(self.on_treeclicked)
        self.treeWidget.currentItemChanged.connect(self.on_treeclicked)
        
        self.ProcessPhase.clicked.connect(self.on_ProcessPhase)
        self.data = self.image_init()
        self.create_imagepage1()
        self.imagepage1_draw()
        self.create_imagepage2()
        self.imagepage2_draw()
        self.actionOpen.setChecked(True)
        self.action3D_view_2.triggered.connect(self.On_3D_view_2)   
        self.actionAxial.triggered.connect(self.On_Axial)   
        self.actionCoronal.triggered.connect(self.On_Coronal)   
        self.actionSaggital.triggered.connect(self.On_Saggital)
        self.treeWidget.setColumnWidth(0,130)
        self.treeWidget_parameters.setColumnWidth(0,140)


        self.hist_fig = Figure((2.0, 1.0), dpi=100,facecolor='none')
        self.hist_Diag = FigureCanvas(self.hist_fig)
        self.hist_Diag.setParent(self.frame_hist)
        self.hist_Diag.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.hist_Diag.setFocus()
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.hist_Diag)  
        self.frame_hist.setLayout(vbox)
        self.diag_axes1 = self.hist_fig.add_subplot(111,axis_bgcolor='none',alpha=0)
        self.diag_axes1.axes.get_yaxis().set_visible(False)
        self.diag_axes1.spines['top'].set_visible(False)
        self.diag_axes1.spines['right'].set_visible(False)
        self.diag_axes1.spines['left'].set_visible(False)
        self.diag_axes1.xaxis.set_ticks_position('bottom')
        self.hist_fig.subplots_adjust(top=0.95, bottom=0.2, left=0.05, right=0.95)
        #self.histogram1,self.diagline1,=self.diag_axes1.plot([0,1],[0 ,1],'g',
        #                                      [0,1],[0 ,1], 'r')

        self.histogram1,=self.diag_axes1.plot([0,1],[0 ,0],'g')
        self.diagline1, =self.diag_axes1.plot([0,1],[0 ,0], 'r')
        self.diagline2, =self.diag_axes1.plot([0,1],[0 ,0], 'b')
        self.diagline3, =self.diag_axes1.plot([0,1],[0 ,0], 'b')
        
    def On_3D_view_2(self):
        self.actionSaggital.setChecked(False)
        self.actionAxial.setChecked(False)
        self.actionCoronal.setChecked(False)
        self.actionSave.setChecked(False)
        self.actionOpen.setChecked(False)
        self.stackedWidget_image.setCurrentIndex(0)
        self.imagepage1_init() 
        
    def On_Axial(self):
        self.actionSaggital.setChecked(False)
        self.action3D_view_2.setChecked(False)
        self.actionCoronal.setChecked(False)      
        self.actionSave.setChecked(False)
        self.actionOpen.setChecked(False)
        self.stackedWidget_image.setCurrentIndex(1)
        self.horizontalScrollBar_z_2panel.setRange(0,self.imagesize[2]-1)
        self.horizontalScrollBar_z_2panel.setValue(int(self.imagesize[2]/2))
        self.on_imagepage2_init()  
        
    def On_Coronal(self):
        self.actionSaggital.setChecked(False)
        self.actionAxial.setChecked(False)
        self.action3D_view_2.setChecked(False)
        self.actionSave.setChecked(False)
        self.actionOpen.setChecked(False)
        self.stackedWidget_image.setCurrentIndex(1)
        self.horizontalScrollBar_z_2panel.setRange(0,self.imagesize[1]-1)
        self.horizontalScrollBar_z_2panel.setValue(int(self.imagesize[1]/2))
        self.on_imagepage2_init()  
        
    def On_Saggital(self):
        self.action3D_view_2.setChecked(False)
        self.actionAxial.setChecked(False)
        self.actionCoronal.setChecked(False)
        self.actionSave.setChecked(False)
        self.actionOpen.setChecked(False)
        self.stackedWidget_image.setCurrentIndex(1)
        self.horizontalScrollBar_z_2panel.setRange(0,self.imagesize[0]-1)
        self.horizontalScrollBar_z_2panel.setValue(int(self.imagesize[0]/2))
        self.on_imagepage2_init()  
        
    def on_actionOpen(self):
        self.actionSave.setChecked(False)
        self.actionSaggital.setChecked(False)
        self.actionAxial.setChecked(False)
        self.actionCoronal.setChecked(False)
        self.action3D_view_2.setChecked(False)
        
    def on_actionSave(self):
        self.actionOpen.setChecked(False)
        self.action3D_view_2.setChecked(False)
        self.actionSaggital.setChecked(False)
        self.actionAxial.setChecked(False)
        self.actionCoronal.setChecked(False)

    def on_AdjustContrast(self):
        v_number=self.spinBox_v_3panel.value()   
        if self.actionAdjust_Contrast.isChecked():
            self.stackedWidget_control.setCurrentIndex(1)

            currentvalue=self.spinBox_z_3panel.value()
            if self.imagepage1_panelselect==1:
                imgdatashow1 = self.imagdata1[:,:,currentvalue,v_number]
            elif   self.imagepage1_panelselect==2:         
                imgdatashow1 = self.imagdata2[:,:,currentvalue,v_number]
            imageshape=imgdatashow1.shape
            imgdatashow1=np.reshape(imgdatashow1,imageshape[0]*imageshape[1])
            imgdatashow1=imgdatashow1[imgdatashow1!=0]
            self.minval=np.amin(imgdatashow1)
            self.maxval=np.amax(imgdatashow1)
            mean1=(self.minval+self.maxval)/2.
            range1=(self.maxval-self.minval)/2.
    
            self.doubleSpinBox_min.setMaximum(self.maxval+range1)
            self.doubleSpinBox_min.setMinimum(self.minval-range1)
            self.doubleSpinBox_min.setValue(self.minval)
            
            self.doubleSpinBox_max.setMaximum(self.maxval+range1)
            self.doubleSpinBox_max.setMinimum(self.minval-range1)
            self.doubleSpinBox_max.setValue(self.maxval)
            
            self.doubleSpinBox_mid.setMaximum(self.maxval)
            self.doubleSpinBox_mid.setMinimum(self.minval)
            self.doubleSpinBox_mid.setValue(mean1)
            
            self.doubleSpinBox_range.setMaximum(range1*2)
            self.doubleSpinBox_range.setMinimum(0)
            self.doubleSpinBox_range.setValue(range1)
            
            hist, bin_edges = np.histogram(imgdatashow1, 100)
            bin_edges1=(bin_edges[1:]+bin_edges[:-1])/2.0
            histmax=hist.max()
            self.histogram1.set_data(bin_edges1,hist)
            self.diagline1.set_data([mean1,mean1],[0,histmax*1.2])
            self.diagline2.set_data([mean1-range1,mean1-range1],[0,histmax*1.1])
            self.diagline3.set_data([mean1+range1,mean1+range1],[0,histmax*1.1])
            self.diagline1.set_linewidth(2)
            self.diagline2.set_linewidth(2)
            self.diagline3.set_linewidth(2)
            self.diag_axes1.get_yaxis().set_visible(False)
            self.diag_axes1.relim()
            self.diag_axes1.autoscale_view()
            self.hist_Diag.draw()
            
            self.hSliderIntensity.valueChanged.connect(self.contrastAdjSlider)
            self.hSliderContrast.valueChanged.connect(self.contrastAdjSlider)        
        else:
            self.stackedWidget_control.setCurrentIndex(0)
        
        return
#        self.MyDialog=MyDialog()
#        self.diag_fig = Figure((2.0, 1.0), dpi=100,facecolor='none')
#        self.canvas_Diag = FigureCanvas(self.diag_fig)
#        self.canvas_Diag.setParent(self.MyDialog.frame)
        
        
    def contrastAdjSlider(self):
        RangeValue=self.hSliderContrast.value()
        MiddleValue = self.hSliderIntensity.value()
        Mid_adj=self.minval*(1-MiddleValue/100.)+self.maxval*(MiddleValue/100.)
        self.doubleSpinBox_mid.setValue(Mid_adj)
        range1=(self.maxval-self.minval)/2.
        range_adj=range1*(RangeValue+1)/51.
        min_adj = Mid_adj-range_adj
        max_adj = Mid_adj+range_adj
        self.doubleSpinBox_range.setValue(range_adj)
        self.doubleSpinBox_min.setValue(min_adj)
        self.doubleSpinBox_max.setValue(max_adj)
        self.diagline1.set_xdata([Mid_adj, Mid_adj])
        self.diagline2.set_xdata([min_adj, min_adj])
        self.diagline3.set_xdata([max_adj, max_adj])
        self.imagepage1_leftimage.set_clim(min_adj,max_adj)
        self.imagepage1_mid_image.set_clim(min_adj,max_adj)
        self.imagepage1_rightimage.set_clim(min_adj,max_adj)
        self.diag_axes1.relim()
        self.diag_axes1.autoscale_view()
        self.imagepage1_canvas.draw()
        self.hist_Diag.draw()
        self.statusBar().showMessage(np.array_str(min_adj))

    
    def contrastAdjYes(self):
        self.statusBar().showMessage('Contrast Adjusted' )
        self.MyDialog.close()
        
    def contrastAdjNO(self):
        pass

    def on_treeclicked(self):
        row = self.treeWidget.selectionModel().currentIndex().row()
        if self.actionOpen.isChecked():    
            if row in ([0,1,4,6,7]) :  
                fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file','')
                img = nib.load(fname)
                imagedata= img.get_data()
                if len(img.shape)==3:
                     imagedata=imagedata.reshape((img.shape[0],img.shape[1],img.shape[2],1)) 
                root1 = self.treeWidget.invisibleRootItem()
                Str1=np.array_str(np.array(img.shape))
                if row == 0:
                    self.imageheader=img.header
                    self.affine= img.affine
                    self.voxelsize=img.header.get_zooms()
                    self.imagesize=img.shape
                    self.spinBox_x_3panel.setMinimum(0)
                    self.spinBox_x_3panel.setMaximum(self.imagesize[0]-1)
                    self.spinBox_y_3panel.setMinimum(0)
                    self.spinBox_y_3panel.setMaximum(self.imagesize[1]-1)
                    self.spinBox_z_3panel.setMinimum(0)
                    self.spinBox_z_3panel.setMaximum(self.imagesize[2]-1)
                    if len(img.shape)==4:
                        self.spinBox_v_3panel.setMaximum(self.imagesize[3]-1)
                        self.spinBox_v_2panel.setMaximum(self.imagesize[3]-1)
                    else:
                        self.spinBox_v_3panel.setMaximum(0)
                    
                    self.spinBox_x_3panel.setValue(self.imagesize[0]/2)                  
                    self.spinBox_y_3panel.setValue(self.imagesize[1]/2)                  
                    self.spinBox_z_3panel.setValue(self.imagesize[2]/2)         
                    
                    self.stackedWidget_image.setCurrentIndex(0)
                    self.horizontalScrollBar_z_2panel.setRange(0,self.imagesize[2]-1)
                    self.horizontalScrollBar_z_2panel.setValue(int(self.imagesize[2]/2))
                    self.lineEdit_z_2panel.setText(str(int(self.imagesize[2]/2)))
                    self.treeWidget_parameters.topLevelItem(0).child(1).setText(1,str(self.voxelsize))  
                    
                    self.imagepage1_panelselect=1
                    self.magdata = imagedata/np.mean(imagedata)
                    self.phasedata = np.zeros(self.magdata.shape)
                    self.imagdata1 = self.magdata
                    self.imagdata2 = self.phasedata
                    root1.child(0).setText(1,Str1) 
                    self.maskdata = self.magdata>0
                    root1.child(4).setText(1,Str1)
                
                    
                elif row ==1:
                    imagedata=imagedata/(imagedata.max()-imagedata.min())*2*np.pi
                    self.phasedata = imagedata
                    self.imagdata2 = self.phasedata
                    self.imagepage1_panelselect=2
                    root1.child(1).setText(1,Str1) 
                elif row ==4:
                    self.maskdata = imagedata
                    self.imagdata2 = self.maskdata
                    self.imagepage1_panelselect=2
                    root1.child(4).setText(1,Str1)  
                self.imagepage1_init()    
                self.statusBar().showMessage(fname + ' ' + str(self.imagesize))
            
        elif self.actionSave.isChecked():       # writting data
            if row ==0:
                filename = QtGui.QFileDialog.getSaveFileName(self,"saveFile","Mag",filter ="nifti (*.nii.gz)") 
                array_img = nib.Nifti1Image(self.magdata , self.affine)
                nib.save(array_img, filename)
            elif row ==1:
                filename = QtGui.QFileDialog.getSaveFileName(self,"saveFile","Phase",filter ="nifti (*.nii.gz)") 
                array_img = nib.Nifti1Image(self.phasedata , self.affine)
                nib.save(array_img, filename)

            elif row ==4:
                filename = QtGui.QFileDialog.getSaveFileName(self,"saveFile","Mask",filter ="nifti (*.nii.gz)") 
                array_img = nib.Nifti1Image(self.maskdata , self.affine)
                nib.save(array_img, filename)
                
            elif row ==6:
                filename = QtGui.QFileDialog.getSaveFileName(self,"saveFile","Phase_Path",filter ="nifti (*.nii.gz)") 
                array_img = nib.Nifti1Image(self.Phase_Path , self.affine)
                nib.save(array_img, filename)
            elif row ==7:
                filename = QtGui.QFileDialog.getSaveFileName(self,"saveFile","Phase_LP",filter ="nifti (*.nii.gz)") 
                array_img = nib.Nifti1Image(self.Phase_LP , self.affine)
                nib.save(array_img, filename)
            pass
        elif self.action3D_view_2.isChecked():        # viewing data
            if row ==0:
                self.imagepage1_panelselect = 1
                self.imagepage1_init()  
            if row ==1:
                self.imagepage1_panelselect = 2
                self.imagdata2 = self.phasedata
                self.imagepage1_init()  
                
            if row ==4:
                self.imagepage1_panelselect = 2
                self.imagdata2 = self.maskdata
                self.imagepage1_init()  
                
            if row ==6:
                self.imagepage1_panelselect = 2
                self.imagdata2 = self.Phase_Path
                self.imagepage1_init()  
            if row ==7:
                self.imagepage1_panelselect = 2
                self.imagdata2 = self.Phase_LP
                self.imagepage1_init()  
                pass
        else:        # viewing data
            if row ==0:
                self.imagepage1_panelselect = 2
                if self.radioButton_left.isChecked():
                    self.imagdata1 = self.magdata
                else:
                    self.imagdata2 = self.magdata    
                self.on_imagepage2_init()  
            if row ==1:
                self.imagepage1_panelselect = 2
                if self.radioButton_left.isChecked():
                    self.imagdata1 = self.phasedata
                else:
                    self.imagdata2 = self.phasedata
                    
                self.on_imagepage2_init()  
            if row ==4:
                self.imagepage1_panelselect = 2
                if self.radioButton_left.isChecked():
                    self.imagdata1 = self.maskdata
                else:
                    self.imagdata2 = self.maskdata
                self.on_imagepage2_init()  
                
            if row ==6:
                self.imagepage1_panelselect = 2
                if self.radioButton_left.isChecked():
                    self.imagdata1 = self.Phase_Path
                else:
                    self.imagdata2 = self.Phase_Path
                self.on_imagepage2_init()  
            if row ==7:
                self.imagepage1_panelselect = 2
                if self.radioButton_left.isChecked():
                    self.imagdata1 = self.Phase_LP
                else:
                    self.imagdata2 = self.Phase_LP
                self.on_imagepage2_init()  
                pass
    def on_itemchanged(self):
        row = self.treeWidget.selectionModel().currentIndex().row()
        print('this is good')        
        print(row)
        pass
    
    def on_ProcessPhase(self):
        #maskdata = 1-self.maskdata
        #phi_wrapped_masked = np.ma.array(self.phasedata, dtype = np.float32, mask = maskdata)
#        start = timeit.timeit()
#        phi_unwrapped_masked = unwrap(phi_wrapped_masked)
#        end = timeit.timeit()
        self.actionOpen.setChecked(False)
        self.action3D_view_2.setChecked(True)
        root = self.treeWidget_parameters.invisibleRootItem()
        #child_count = root.childCount()
        self.TE1 = self.treeWidget_parameters.topLevelItem(0).child(0).text(1)
        #self.voxelsize =  self.treeWidget_parameters.topLevelItem(0).child(1).text(1)        
        self.deltaTE   =  self.treeWidget_parameters.topLevelItem(0).child(2).text(1)        
        self.B0 =         self.treeWidget_parameters.topLevelItem(0).child(3).text(1)        
        self.B0_dir    =  self.treeWidget_parameters.topLevelItem(0).child(4).text(1)
        self.VSHARPMaxR =  self.treeWidget_parameters.topLevelItem(1).child(0).text(1)        
        self.iHARPnIter =  self.treeWidget_parameters.topLevelItem(2).child(0).text(1)        
        self.VSHARP    =  root.child(1).text(1)        
        self.iHARP     =  root.child(2).text(1)     

        self.statusBar().showMessage('start path unwrapping, please wait')
        start = time.time()
        #phi_unwrapped_masked = unwrap(phi_wrapped_masked)
        self.Phase_Path=np.zeros(self.imagesize)
        for ii in range(self.imagesize[3]):
            #self.Phase_Path[:,:,:,ii] = unwrap(self.phasedata[:,:,:,ii])
            self.Phase_Path[:,:,:,ii] = unwrap_phase(self.phasedata[:,:,:,ii])
            #self.statusBar().showMessage('Path-Based phase Unwrapping Echo ' + str(ii))
            print('Path-Based phase Unwrapping Echo ' + str(ii))
            
        end = time.time()
        elapsed_Time1=end - start
        self.statusBar().showMessage('start LP unwrapping, please wait')
        start = time.time()
        self.Phase_LP=np.zeros(self.imagesize)
        #for ii in range(self.imagesize[3]):
        #    self.Phase_LP[:,:,:,ii], dim =  LPuwp.uwnrap_LP(self.phasedata[:,:,:,ii],self.voxelsize)
        #    self.statusBar().showMessage('Laplacian-based phase Unwrapping Echo ' + str(ii))
        #self.Phase_LP=self.Phase_LP*self.maskdata
        
        end = time.time()
        elapsed_Time2=end - start
        
        self.imagdata2 = self.Phase_LP
        self.imagepage1_panelselect=2
        self.imagepage1_init()    
        self.statusBar().showMessage('Path-based' + str(round(elapsed_Time1,1)) +' Laplacian:' +str(round(elapsed_Time2,1))  + ' sec')

        root1 = self.treeWidget.invisibleRootItem()
        Str1=np.array_str(np.array(self.imagesize))
        root1.child(6).setText(1,Str1)                     
        root1.child(7).setText(1,Str1)                     

    def on_Opened(self):
        self.stackedWidget_control.setCurrentIndex(0)
        self.actionAdjust_Contrast.setChecked(False)
        self.statusBar().showMessage('Contrast Adjusted' )

    def on_Zoom_in(self):
        self.statusBar().showMessage('Zoom in' )
        self.actionZoom_out.setChecked(False)
        
    def on_Zoomout(self):
        self.statusBar().showMessage('Zoom out' )
        self.actionZoom_in.setChecked(False)
       
    def image_init(self):
        return np.arange(24).reshape([6, 4]).copy()
           
    def create_imagepage2(self):
        self.imagepage2_fig = Figure((5.0, 4.0), dpi=100,facecolor='none')
        self.canvas = FigureCanvas(self.imagepage2_fig)
        self.canvas.setParent(self.frame_2panels)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.setFocus()
        self.canvas.mpl_connect('key_press_event', self.on_key_press)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)  
        self.frame_2panels.setLayout(vbox)

    def imagepage2_draw(self):
        self.imagepage2_fig.clear()
        self.imagepage2_axes1 = self.imagepage2_fig.add_subplot(121,axis_bgcolor='red',alpha=0)
        self.page2_axes2 = self.imagepage2_fig.add_subplot(122)
        self.imagepage2_fig.subplots_adjust(top=0.999, bottom=0.01, left=0.001, right=0.999,wspace=0.005,hspace=0.005)
        self.leftimage=self.imagepage2_axes1.imshow(self.data, interpolation='nearest', aspect='auto', cmap='gray',origin="lower")
        self.rightimage=self.page2_axes2.imshow(self.data, interpolation='nearest', aspect='auto', cmap='gray',origin="lower")
        self.imagepage2_axes1.set_axis_off()
        self.page2_axes2.set_axis_off()

        self.leftimage.set_extent((-0.5,20-.5,20-.5,-.5))
        self.imagepage2_axes1.set_xlim((-0.5,20-.5))
        self.imagepage2_axes1.set_ylim((20-.5,-.5))
        self.imagepage2_axes1.get_xaxis().set_visible(False)
        self.imagepage2_axes1.get_yaxis().set_visible(False)
        self.canvas.mpl_connect('scroll_event', self.on_page2mousescroll)
        self.canvas.mpl_connect('button_press_event', self.on_page2buttonpress)
        self.canvas.mpl_connect('motion_notify_event', self.on_page2_mousemove)
        self.canvas.mpl_connect('button_release_event', self.on_page2buttonrelease)
        self.canvas.draw()
        
    def on_page2buttonpress(self, event):
        self.ImagePageMousePressed=1
        self.ImagePageMouseLoc=np.array([event.xdata,event.ydata])
        yy=self.imagepage2_axes1.get_ylim()
        xx=self.imagepage2_axes1.get_xlim()
        self.ImagePagexylims=np.array([xx[1]-xx[0],yy[1]-yy[0]])

    def on_page2_mousemove(self, event):
        pass
        #if self.ImagePageMousePressed:
            #ImagePageMouseLoc1=np.array([(event.xdata),(event.ydata)])
            #print((ImagePageMouseLoc1-self.ImagePageMouseLoc)/self.ImagePagexylims)
        
    def on_page2buttonrelease(self, event):
        self.ImagePageMousePressed=0
        
    def on_page2mousescroll(self,event):
        
        delta=self.spinBox_P2_scroll.value()
        if event.button == 'up':
            currentvalue=self.horizontalScrollBar_z_2panel.value()
            self.horizontalScrollBar_z_2panel.setValue(currentvalue+delta)
        else:
            currentvalue=self.horizontalScrollBar_z_2panel.value()
            self.horizontalScrollBar_z_2panel.setValue(currentvalue-delta)

    def create_imagepage1(self):
        self.imagepage1_fig = Figure((5.0, 4.0), dpi=100,facecolor='none')
        self.imagepage1_canvas = FigureCanvas(self.imagepage1_fig)
        self.imagepage1_canvas.setParent(self.frame_3panels)
        self.imagepage1_canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.imagepage1_canvas.setFocus()
        
        vbox =QtGui.QVBoxLayout()
        vbox.setContentsMargins(0,0,0,0)
        vbox.addWidget(self.imagepage1_canvas)  
        self.frame_3panels.setLayout(vbox)

    def imagepage1_draw(self):
        self.imagepage1_fig.clear()
        self.imagepage1_axes1 = self.imagepage1_fig.add_subplot(131,axis_bgcolor='red',alpha=0)
        self.imagepage1_axes2 = self.imagepage1_fig.add_subplot(132)
        self.imagepage1_axes3 = self.imagepage1_fig.add_subplot(133)
        self.imagepage1_fig.subplots_adjust(top=0.99, bottom=0.01, left=0.001, right=0.999,wspace=0.005,hspace=0.005)
        self.imagepage1_leftimage =self.imagepage1_axes1.imshow(self.data, interpolation='nearest', aspect='auto', cmap='gray',origin="lower")
        self.imagepage1_mid_image =self.imagepage1_axes2.imshow(self.data, interpolation='nearest', aspect='auto', cmap='gray',origin="lower")
        self.imagepage1_rightimage=self.imagepage1_axes3.imshow(self.data, interpolation='nearest', aspect='auto', cmap='gray',origin="lower")
        self.imagepage1_axes1.set_axis_off()
        self.imagepage1_axes2.set_axis_off()
        self.imagepage1_axes3.set_axis_off()
        self.imagepage1_leftimage.set_extent((-0.5,20-.5,20-.5,-.5))
        self.imagepage1_axes1.set_xlim((-0.5,20-.5))
        self.imagepage1_axes1.set_ylim((20-.5,-.5))
        self.imagepage1_axes1.get_xaxis().set_visible(False)
        self.imagepage1_axes1.get_yaxis().set_visible(False)

        self.imagepage1_canvas.mpl_connect('button_press_event', self.on_page1_mousepress)
        self.imagepage1_canvas.mpl_connect('motion_notify_event', self.on_page1_mousemove)
        self.imagepage1_canvas.mpl_connect('button_release_event', self.on_page1_mouserelease)
        self.imagepage1_canvas.mpl_connect('scroll_event', self.on_mousescroll)

#        self.cursor = Cursor(self.imagepage1_axes1, color='yellow', linewidth=0.5 )
#        self.cursor.horizOn = True
#        self.cursor.vertOn = True
        self.imagepage1_canvas.draw()

    def on_key_press(self, event):
        if event.key == 'up':
            currentvalue=self.horizontalScrollBar_z_2panel.value()
            self.horizontalScrollBar_z_2panel.setValue(currentvalue+1)
        elif event.key == 'down':
            currentvalue=self.horizontalScrollBar_z_2panel.value()
            self.horizontalScrollBar_z_2panel.setValue(currentvalue-1)

    def fileQuit(self):
        self.close()
        
    def on_2panel_view(self):
        self.ControlPanel.show()
        self.statusBar().showMessage('Matplot viewer' )
        self.stackedWidget_image.setCurrentIndex(0)
        self.stackedWidget_control.setCurrentIndex(0)
    
    def on_3panel_view(self):
        self.ControlPanel.show()
        self.statusBar().showMessage('VTK viewer' )
        self.stackedWidget_image.setCurrentIndex(1)
        self.stackedWidget_control.setCurrentIndex(1)
    
    def on_page1_mousepress(self, event):
        v_number=self.spinBox_v_3panel.value()   
        if event.inaxes==self.imagepage1_axes1:
            if self.imagepage1_panelselect==1:
                imgdatashow1 = self.imagdata1[int(event.xdata),:,:,v_number]
                imgdatashow2 = self.imagdata1[:,int(event.ydata),:,v_number]
            elif self.imagepage1_panelselect==2:
                imgdatashow1 = self.imagdata2[int(event.xdata),:,:,v_number]
                imgdatashow2 = self.imagdata2[:,int(event.ydata),:,v_number]
            self.imagepage1_rightimage.set_data(imgdatashow1.T)
            self.imagepage1_mid_image.set_data(imgdatashow2.T)
            self.spinBox_x_3panel.setValue(int(event.xdata))
            self.spinBox_y_3panel.setValue(int(event.ydata))
            
        elif event.inaxes==self.imagepage1_axes2:
                if self.imagepage1_panelselect==1:
                    imgdatashow3 = self.imagdata1[int(event.xdata),:,:,v_number]
                    imgdatashow4 = self.imagdata1[:,:,int(event.ydata),v_number]
                elif self.imagepage1_panelselect==2:
                    imgdatashow3 = self.imagdata2[int(event.xdata),:,:,v_number]
                    imgdatashow4 = self.imagdata2[:,:,int(event.ydata),v_number]
                self.imagepage1_rightimage.set_data(imgdatashow3.T)
                self.imagepage1_leftimage.set_data(imgdatashow4.T)
                self.spinBox_x_3panel.setValue(int(event.xdata))
                self.spinBox_z_3panel.setValue(int(event.ydata))
                
        elif event.inaxes==self.imagepage1_axes3:
                if self.imagepage1_panelselect==1:
                    imgdatashow5 = self.imagdata1[:,int(event.xdata),:,v_number]
                    imgdatashow6 = self.imagdata1[:,:,int(event.ydata),v_number]
                elif self.imagepage1_panelselect==2:
                    imgdatashow5 = self.imagdata2[:,int(event.xdata),:,v_number]
                    imgdatashow6 = self.imagdata2[:,:,int(event.ydata),v_number]
                self.imagepage1_mid_image.set_data(imgdatashow5.T)
                self.imagepage1_leftimage.set_data(imgdatashow6.T)
                self.spinBox_y_3panel.setValue(int(event.xdata))
                self.spinBox_z_3panel.setValue(int(event.ydata))
        else:
            pass        
        
        self.pressed = 1
        self.imagepage1_canvas.draw() 
        #QApplication.setOverrideCursor(QCursor(QtCore.Qt.CrossCursor))

    def on_page1_mousemove(self, event):
        v_number=self.spinBox_v_3panel.value()   
        if self.pressed == 1:
            if event.inaxes==self.imagepage1_axes1:
                if self.imagepage1_panelselect==1:
                    imgdatashow1 = self.imagdata1[int(event.xdata),:,:,v_number]
                    
                    imgdatashow2 = self.imagdata1[:,int(event.ydata),:,v_number]
                elif self.imagepage1_panelselect==2:
                    imgdatashow1 = self.imagdata2[int(event.xdata),:,:,v_number]
                    imgdatashow2 = self.imagdata2[:,int(event.ydata),:,v_number]
                self.imagepage1_rightimage.set_data(imgdatashow1.T)
                self.imagepage1_mid_image.set_data(imgdatashow2.T)
                self.spinBox_x_3panel.setValue(int(event.xdata))
                self.spinBox_y_3panel.setValue(int(event.ydata))

                
            elif event.inaxes==self.imagepage1_axes2:
                if self.imagepage1_panelselect==1:
                    imgdatashow3 = self.imagdata1[int(event.xdata),:,:,v_number]
                    imgdatashow4 = self.imagdata1[:,:,int(event.ydata),v_number]
                elif self.imagepage1_panelselect==2:
                    imgdatashow3 = self.imagdata2[int(event.xdata),:,:,v_number]
                    imgdatashow4 = self.imagdata2[:,:,int(event.ydata),v_number]
                self.imagepage1_rightimage.set_data(imgdatashow3.T)
                self.imagepage1_leftimage.set_data(imgdatashow4.T)
                self.spinBox_x_3panel.setValue(int(event.xdata))
                self.spinBox_z_3panel.setValue(int(event.ydata))

            elif event.inaxes==self.imagepage1_axes3:
                if self.imagepage1_panelselect==1:
                    imgdatashow5 = self.imagdata1[:,int(event.xdata),:,v_number]
                    imgdatashow6 = self.imagdata1[:,:,int(event.ydata),v_number]
                elif self.imagepage1_panelselect==2:
                    imgdatashow5 = self.imagdata2[:,int(event.xdata),:,v_number]
                    imgdatashow6 = self.imagdata2[:,:,int(event.ydata),v_number]
                self.imagepage1_mid_image.set_data(imgdatashow5.T)
                self.imagepage1_leftimage.set_data(imgdatashow6.T)
                
                self.spinBox_y_3panel.setValue(int(event.xdata))
                self.spinBox_z_3panel.setValue(int(event.ydata))
            else:
                pass
        self.imagepage1_canvas.draw()


    def on_page1_mouserelease(self, event):
        self.pressed = 0
        #QApplication.restoreOverrideCursor()
        
    def on_mousescroll(self, event):
#        #spinBox_page1_ScrollSteps
        delta=self.spinBox_p1_Scroll.value()
        if event.inaxes==self.imagepage1_axes1:
            if event.button == 'up':
                currentvalue=self.spinBox_z_3panel.value()
                self.statusBar().showMessage(str(currentvalue))
                self.spinBox_z_3panel.setValue(currentvalue+delta)
            else:
                currentvalue=self.spinBox_z_3panel.value()
                self.statusBar().showMessage(str(currentvalue))
                self.spinBox_z_3panel.setValue(currentvalue-delta)
                
        elif event.inaxes==self.imagepage1_axes2:
            if event.button == 'up':
                currentvalue=self.spinBox_y_3panel.value()
                self.statusBar().showMessage(str(currentvalue))
                self.spinBox_y_3panel.setValue(currentvalue+delta)
            else:
                currentvalue=self.spinBox_y_3panel.value()
                self.statusBar().showMessage(str(currentvalue))
                self.spinBox_y_3panel.setValue(currentvalue-delta)
        else:
            if event.button == 'up':
                currentvalue=self.spinBox_x_3panel.value()
                self.statusBar().showMessage(str(currentvalue))
                self.spinBox_x_3panel.setValue(currentvalue+delta)
            else:
                currentvalue=self.spinBox_x_3panel.value()
                self.statusBar().showMessage(str(currentvalue))
                self.spinBox_x_3panel.setValue(currentvalue-delta)
                
    def on_ROI_view(self):
        self.statusBar().showMessage('ROI' )
        self.ControlPanel.show()
        self.stackedWidget_image.setCurrentIndex(2)
        self.stackedWidget_control.setCurrentIndex(1)
        
    def imagepage1_init(self):
        
        z_number=self.spinBox_z_3panel.value() 
        y_number=self.spinBox_y_3panel.value() 
        x_number=self.spinBox_x_3panel.value()   
        v_number=self.spinBox_v_3panel.value()   
        
        if self.imagepage1_panelselect == 1:
            imgdatashow1 = self.imagdata1[:,:,z_number,v_number]
            imgdatashow2 = self.imagdata1[:,y_number,:,v_number]
            imgdatashow3 = self.imagdata1[x_number,:,:,v_number]
        elif self.imagepage1_panelselect == 2:
            imgdatashow1 = self.imagdata2[:,:,z_number,v_number]
            imgdatashow2 = self.imagdata2[:,y_number,:,v_number]
            imgdatashow3 = self.imagdata2[x_number,:,:,v_number]

        self.pg1dataShape1=imgdatashow1.shape
        self.imagepage1_leftimage.set_data(imgdatashow1.T)
        self.pg1dataShape2=imgdatashow2.shape
        self.imagepage1_mid_image.set_data(imgdatashow2.T)
        self.pg1dataShape3=imgdatashow3.shape
        self.imagepage1_rightimage.set_data(imgdatashow3.T)
        self.imagepage1_leftimage.set_clim(np.amin(imgdatashow1),np.amax(imgdatashow1))
        self.imagepage1_mid_image.set_clim(np.amin(imgdatashow2),np.amax(imgdatashow2))
        self.imagepage1_rightimage.set_clim(np.amin(imgdatashow3),np.amax(imgdatashow3))
        self.axis_extent_adj('page1')
        
        ####### in revision
        self.imagepage1_canvas.draw()
         
    def onSpinBoxZ_Changed(self,value):
        v_number=self.spinBox_v_3panel.value()  
        if self.imagepage1_panelselect==1:
            imgdatashow1 = self.imagdata1[:,:,int(value[2:]),v_number]
        elif self.imagepage1_panelselect==2:
            imgdatashow1 = self.imagdata2[:,:,int(value[2:]),v_number]
        self.imagepage1_leftimage.set_data(imgdatashow1.T)

        self.imagepage1_canvas.draw()
    def onSpinBoxY_Changed(self,value):
        v_number=self.spinBox_v_3panel.value()
        if self.imagepage1_panelselect==1:
            imgdatashow1 = self.imagdata1[:,int(value[2:]),:,v_number]
        elif self.imagepage1_panelselect==2:
            imgdatashow1 = self.imagdata2[:,int(value[2:]),:,v_number]
        self.imagepage1_mid_image.set_data(imgdatashow1.T)
        self.imagepage1_canvas.draw()
    
    def onSpinBoxX_Changed(self,value):
        v_number=self.spinBox_v_3panel.value()
        if self.imagepage1_panelselect==1:
            imgdatashow1 = self.imagdata1[int(value[2:]),:,:,v_number]
        elif self.imagepage1_panelselect==2:
            imgdatashow1 = self.imagdata2[int(value[2:]),:,:,v_number]
        self.imagepage1_rightimage.set_data(imgdatashow1.T)

        self.imagepage1_canvas.draw()

    def onSpinBoxV_Changed(self,value):
        z_number=self.spinBox_z_3panel.value()
        x_number=self.spinBox_x_3panel.value()
        y_number=self.spinBox_y_3panel.value()
        if self.imagepage1_panelselect==1:
            imgdatashow1 = self.imagdata1[:,:,z_number,int(value[2:])]
            self.imagepage1_leftimage.set_data(imgdatashow1.T)
            
            imgdatashow1 = self.imagdata1[x_number,:,:,int(value[2:])]
            self.imagepage1_rightimage.set_data(imgdatashow1.T)
            
            imgdatashow1 = self.imagdata1[:,y_number,:,int(value[2:])]
            self.imagepage1_mid_image.set_data(imgdatashow1.T)
            
            
        elif self.imagepage1_panelselect==2:
#            imgdatashow1 = self.imagdata2[:,:,z_number,int(value[2:])]
            
            imgdatashow1 = self.imagdata2[:,:,z_number,int(value[2:])]
            self.imagepage1_leftimage.set_data(imgdatashow1.T)
            
            imgdatashow1 = self.imagdata2[x_number,:,:,int(value[2:])]
            self.imagepage1_rightimage.set_data(imgdatashow1.T)
            
            imgdatashow1 = self.imagdata2[:,y_number,:,int(value[2:])]
            self.imagepage1_mid_image.set_data(imgdatashow1.T)

        self.imagepage1_canvas.draw()
   
    def on_imagepage2_init(self):
        v_number=self.spinBox_v_3panel.value()   

        if self.actionAxial.isChecked():
            z_number=self.spinBox_z_3panel.value()     
            imagdata1show = self.imagdata1[:,:,z_number,v_number]
            self.leftimage.set_data(imagdata1show.T)
            self.p2img1shape=imagdata1show.shape
            imagdata2show = self.imagdata2[:,:,z_number,v_number]
            self.p2img2shape=imagdata1show.shape
            self.rightimage.set_data(imagdata2show.T)
        if self.actionSaggital.isChecked():   
            x_number=self.spinBox_x_3panel.value()     
            imagdata1show = self.imagdata1[x_number,:,:,v_number]
            self.leftimage.set_data(imagdata1show.T)
            self.p2img1shape=imagdata1show.shape
            imagdata2show = self.imagdata2[x_number,:,:,v_number]
            self.p2img2shape=imagdata1show.shape
            self.rightimage.set_data(imagdata2show.T)

        if self.actionCoronal.isChecked():   
            y_number=self.spinBox_y_3panel.value()     
            imagdata1show = self.imagdata1[:,y_number,:,v_number]
            self.leftimage.set_data(imagdata1show.T)
            self.p2img1shape=imagdata1show.shape
            imagdata2show = self.imagdata2[:,y_number,:,v_number]
            self.p2img2shape=imagdata1show.shape
            self.rightimage.set_data(imagdata2show.T)
            
        self.axis_extent_adj('page2')
        self.leftimage.set_clim(np.amin(imagdata1show),np.amax(imagdata1show))
        self.rightimage.set_clim(np.amin(imagdata2show),np.amax(imagdata2show))
        self.canvas.draw()
    
    def axis_extent_adj(self,control):
        if control=='page1':
            if self.imagepage1_leftimage.origin == 'upper':
                self.imagepage1_leftimage.set_extent((-0.5,self.pg1dataShape1[0]-.5,self.pg1dataShape1[1]-.5,-.5))
                self.imagepage1_axes1.set_xlim((-0.5,self.pg1dataShape1[0]-.5))
                self.imagepage1_axes1.set_ylim((self.pg1dataShape1[1]-.5,-.5))
            else:
                self.imagepage1_leftimage.set_extent((-0.5,self.pg1dataShape1[0]-.5,-.5,self.pg1dataShape1[1]-.5))
                self.imagepage1_axes1.set_xlim((-0.5,self.pg1dataShape1[0]-.5))
                self.imagepage1_axes1.set_ylim((-.5,self.pg1dataShape1[1]-.5))
            aspectratio=self.voxelsize[1]/float(self.voxelsize[0])
            self.imagepage1_axes2.set_aspect(aspectratio)
            self.imagepage1_axes1.set_aspect('equal')
    
            if self.imagepage1_mid_image.origin == 'upper':
                self.imagepage1_mid_image.set_extent((-0.5,self.pg1dataShape2[0]-.5,self.pg1dataShape2[1]-.5,-.5))
                self.imagepage1_axes2.set_xlim((-0.5,self.pg1dataShape2[0]-.5))
                self.imagepage1_axes2.set_ylim((self.pg1dataShape2[1]-.5,-.5))
            else:
                self.imagepage1_mid_image.set_extent((-0.5,self.pg1dataShape2[0]-.5,-.5,self.pg1dataShape2[1]-.5))
                self.imagepage1_axes2.set_xlim((-0.5,self.pg1dataShape2[0]-.5))
                self.imagepage1_axes2.set_ylim((-.5,self.pg1dataShape2[1]-.5))
            #self.imagepage1_axes2.set_aspect('equal')
            aspectratio=self.voxelsize[2]/float(self.voxelsize[1])
            self.imagepage1_axes2.set_aspect(aspectratio)
    
            if self.imagepage1_rightimage.origin == 'upper':
                self.imagepage1_rightimage.set_extent((-0.5,self.pg1dataShape3[0]-.5,self.pg1dataShape3[1]-.5,-.5))
                self.imagepage1_axes3.set_xlim((-0.5,self.pg1dataShape3[0]-.5))
                self.imagepage1_axes3.set_ylim((self.pg1dataShape3[1]-.5,-.5))
            else:
                self.imagepage1_rightimage.set_extent((-0.5,self.pg1dataShape3[0]-.5,-.5,self.pg1dataShape3[1]-.5))
                self.imagepage1_axes3.set_xlim((-0.5,self.pg1dataShape3[0]-.5))
                self.imagepage1_axes3.set_ylim((-.5,self.pg1dataShape3[1]-.5))
            aspectratio=self.voxelsize[2]/float(self.voxelsize[0])
            self.imagepage1_axes3.set_aspect(aspectratio)
        
        elif control=='page2':
            if self.leftimage.origin == 'upper':
                self.leftimage.set_extent((-0.5,self.p2img1shape[0]-.5,self.p2img1shape[1]-.5,-.5))
                self.imagepage2_axes1.set_xlim((-0.5,self.p2img1shape[0]-.5))
                self.imagepage2_axes1.set_ylim((self.p2img1shape[1]-.5,-.5))
            else:
                self.leftimage.set_extent((-0.5,self.p2img1shape[0]-.5,-.5,self.p2img1shape[1]-.5))
                self.imagepage2_axes1.set_xlim((-0.5,self.p2img1shape[0]-.5))
                self.imagepage2_axes1.set_ylim((-.5,self.p2img1shape[1]-.5))
    
            if self.rightimage.origin == 'upper':
                self.rightimage.set_extent((-0.5,self.p2img2shape[0]-.5,self.p2img2shape[1]-.5,-.5))
                self.page2_axes2.set_xlim((-0.5,self.p2img2shape[0]-.5))
                self.page2_axes2.set_ylim((self.p2img2shape[1]-.5,-.5))
            else:
                self.rightimage.set_extent((-0.5,self.p2img2shape[0]-.5,-.5,self.p2img2shape[1]-.5))
                self.page2_axes2.set_xlim((-0.5,self.p2img2shape[0]-.5))
                self.page2_axes2.set_ylim((-.5,self.p2img2shape[1]-.5))
        if self.actionAxial.isChecked():
            self.imagepage2_axes1.set_aspect('equal')
            self.page2_axes2.set_aspect('equal')
        if self.actionSaggital.isChecked():   
            aspectratio=self.voxelsize[2]/float(self.voxelsize[1])
            self.imagepage2_axes1.set_aspect(aspectratio)
            self.page2_axes2.set_aspect(aspectratio)
        else:
            pass
        
    def VScrollBarChangeValue(self, value):
        v_number=self.spinBox_v_3panel.value()   
        if self.actionAxial.isChecked():        
            imgdatashow1=self.imagdata1[:,:,value,v_number]
            self.leftimage.set_data(imgdatashow1.T)
            imgdatashow2=self.imagdata2[:,:,value,v_number]
            self.rightimage.set_data(imgdatashow2.T)
            self.lineEdit_z_2panel.setText(str(value))
            self.spinBox_z_3panel.setValue(value)
                
        if self.actionSaggital.isChecked():        
            imgdatashow1=self.imagdata1[value,:,:,v_number]
            self.leftimage.set_data(imgdatashow1.T)
            imgdatashow2=self.imagdata2[value,:,:,v_number]
            self.rightimage.set_data(imgdatashow2.T)
            self.lineEdit_z_2panel.setText(str(value))
            self.spinBox_x_3panel.setValue(value)
            
        if self.actionCoronal.isChecked():        
            imgdatashow1=self.imagdata1[:,value,:,v_number]
            self.leftimage.set_data(imgdatashow1.T)
            imgdatashow2=self.imagdata2[:,value,:,v_number]
            self.rightimage.set_data(imgdatashow2.T)
            self.lineEdit_z_2panel.setText(str(value))
            self.spinBox_y_3panel.setValue(value)

        self.canvas.draw()


    def onspinBox_v_2panel_Changed(self, value):
        zvalue=self.spinBox_z_3panel.value()   
        if self.actionAxial.isChecked(): 
#            zvalue=20
            imgdatashow1=self.imagdata1[:,:,zvalue,int(value[2:])]
            self.leftimage.set_data(imgdatashow1.T)
            imgdatashow2=self.imagdata2[:,:,zvalue,int(value[2:])]
            self.rightimage.set_data(imgdatashow2.T)
            self.lineEdit_z_2panel.setText(str(value))
        if self.actionSaggital.isChecked():        
            imgdatashow1=self.imagdata1[value,:,:,int(value[2:])]
            self.leftimage.set_data(imgdatashow1.T)
            imgdatashow2=self.imagdata2[value,:,:,int(value[2:])]
            self.rightimage.set_data(imgdatashow2.T)
            self.lineEdit_z_2panel.setText(str(value))
            self.spinBox_v_3panel.setValue(value)
            
        if self.actionCoronal.isChecked():        
            imgdatashow1=self.imagdata1[:,value,:,int(value[2:])]
            self.leftimage.set_data(imgdatashow1.T)
            imgdatashow2=self.imagdata2[:,value,:,int(value[2:])]
            self.rightimage.set_data(imgdatashow2.T)
            self.lineEdit_z_2panel.setText(str(value))
            self.spinBox_v_3panel.setValue(value)

        self.canvas.draw()
  
    def about(self):
        QtGui.QMessageBox.about(self, "About",
        """Susceptibility Imageing Toolkit (SITK). Copyright 2015 Wei Li.
This program is developed for medical image processing.
It can be used and modified for non-conmerical use."""
)

        
def main(): 
    a = QtGui.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(a.exec_())
    
if __name__ == "__main__":
    main()  
 


