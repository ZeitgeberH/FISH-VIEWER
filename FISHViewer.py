# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 15:43:55 2021

@author: ZeitgeberH
"""
from pathlib import Path
from PyQt5 import QtGui, QtCore, QtSvg
from pyqtgraph.Qt import QtWidgets
from PyQt5.QtWidgets import QMessageBox, QTableWidgetItem
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from pyqtgraph.parametertree import Parameter, ParameterTree
import gr_pars as AllMyPars
from pyqtgraph import GraphicsLayoutWidget
import sys
import os
import numpy as np
import pandas as pd
import itertools
import time as sysTime
from copy import deepcopy
import glob
from gr_viewClass import (FileView,ParameterView,PlotView,
    TabView,TreeView, TableView,showdialog,getfiles,ImageView_gv)
from gr_allenFISH import GeneRenderAPI
from gr_allenFISH_utils import (
    check_gene_cached,
    check_gene_cached_images,
    load_cached_gene,
    download_and_cache,
    download_and_cache_image,
    queryGene_ncbi,rgb2intensity
)
from gr_sys_utils import base_dir
import pdb
from PIL import Image
Image.MAX_IMAGE_PIXELS = 1000000000  ## silence decompression bomb warning
import tqdm
import webbrowser
from skimage import filters as SKIfilters
from skimage import measure as SKImeasure
from matplotlib.cm import jet as jetcm

class MainWindow(QtWidgets.QMainWindow):
    """
    Main frame.
    """
    def __init__(self, app,  parent = None):
        super(MainWindow, self).__init__(parent)
        self.app = app
        self.create_mainWindow()
        self.setWindowTitle("FISHViewer")
        self.setWindowIcon(pg.QtGui.QIcon('data\\icons\\Fish.png'))
        self.setWindowState(QtCore.Qt.WindowMaximized)

    def create_mainWindow(self):
        # Configure Qt GUI:
        self.make_layout()  
        self.add_menubar()
        self.add_toolbar()
        self.statusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.statusBar)
        self.MouseMode =  pg.ViewBox.RectMode
        self.mainFrame = pg.QtGui.QWidget()
        self.mainFrame.setLayout(self.MainLayout)
        self.setCentralWidget(self.mainFrame)
        self.otherInitStuff()
        self.setRightPanelStrechFactor()
        
    def make_layout(self):
        ## layout
        self.MainLayout = pg.QtGui.QGridLayout()
        # horizontal splitter
        self.frame_splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        # splitter for file browswer (top) and  pul view (bottom)
        self.right_panels = QtWidgets.QSplitter(QtCore.Qt.Vertical, parent = self.frame_splitter)
       
        ### left pannels for main plotting and data view
        self.visulization_view = TabView(self.frame_splitter)
        self.plot_splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical, parent = self.visulization_view)

        self.ISH_view = pg.ImageView(view=pg.PlotItem()) #, SynAtlasFunc = self.SynAtlas)
        self.ISH_view_vLine = pg.InfiniteLine(angle=90, movable=False, name='cursorV1')
        self.ISH_view_hLine = pg.InfiniteLine(angle=0, movable=False, name='cursorH1')
        self.ISH_view.addItem(self.ISH_view_vLine, ignoreBounds=True)
        self.ISH_view.addItem(self.ISH_view_hLine, ignoreBounds=True)        
        self.ISH_view.view.hoverEvent = self.imageHoverEvent
        
        self.ISH_view.view.scene().sigMouseClicked.connect(self.ISH_view_mouseClicked)
        # self.ISH_view.ui.histogram.hide()

        self.plot_splitter.addWidget(self.ISH_view)
        # self.updateFishImage(self.ISH_view, 'atlas_example.jpg')
        self.ISH_view.setColorMap(pg.colormap.get('CET-L1'))
        self.ISH_view.ui.roiBtn.hide()
        self.ISH_view.ui.menuBtn.hide()
        self.ISH_view.roi.hide()
        # self.ISH_view.ui.roiPlot.close()
        # self.ISH_view.ui.roiPlot.hide()
        self.scale1 = pg.ScaleBar(size=50, suffix='pixels')
        self.scale1.setParentItem(self.ISH_view.view.getViewBox())
        self.scale1.anchor((1, 1), (1, 1), offset=(-20, -20))  
        
        self.Expression_view = pg.ImageView(view=pg.PlotItem())
        self.plot_splitter.addWidget(self.Expression_view)
        self.Expression_view_vLine = pg.InfiniteLine(angle=90, movable=False, name='cursorV2')
        self.Expression_view_hLine = pg.InfiniteLine(angle=0, movable=False, name='cursorH2')
        self.Expression_view.addItem(self.Expression_view_vLine, ignoreBounds=True)
        self.Expression_view.addItem(self.Expression_view_hLine, ignoreBounds=True)   
        self.Expression_view.view.hoverEvent = self.imageHoverEvent
        # self.updateFishImage(self.Expression_view, 'atlas_example.jpg')
        self.Expression_view.setColorMap(pg.colormap.get("CET_L16", source='colorcet')) #pg.colormap.get('CET-L1'))
        self.Expression_view.ui.roiBtn.hide()
        self.Expression_view.ui.menuBtn.hide()
        self.Expression_view.roi.hide()
        self.Expression_view.ui.roiPlot.getPlotItem().hideAxis('bottom')
        self.scale2 = pg.ScaleBar(size=50, suffix='pixels')
        self.scale2.setParentItem(self.Expression_view.view.getViewBox())
        self.scale2.anchor((1, 1), (1, 1), offset=(-20, -20)) 
        # self.Expression_view.ui.roiPlot.close()
        # self.Expression_view.ui.roiPlot.hide()
        
        self.Expression_view.view.setXLink(self.ISH_view.view)
        self.Expression_view.view.setYLink(self.ISH_view.view)
        
        self.visulization_view.addTab(self.plot_splitter, 'ABI Views')

        self.frame_splitter.addWidget(self.visulization_view)
        
        ## right panels for parameters and input
        self.frame_splitter.addWidget(self.right_panels)
        self.files_view = TabView(self.right_panels)
                
        ## parameters for ABI FISH images
        self.globalPars = Parameter.create(name='params', type='group', children=AllMyPars.download_options, readonly=False)
        self.globalPars.sigTreeStateChanged.connect(self.event_parameters_stateChange) ## call back
        self.parTree = ParameterTree()
        self.parTree.setHeaderLabels(["Parameter                                  ", "Value"])
        self.parTree.setParameters(self.globalPars, showTop=False)
        self.parTree.setWindowTitle('Global parameter')
        self.files_view.addTab(self.parTree, 'ABA parameters')
                
        self.right_panels.addWidget(self.files_view)
        
        ### Middle panel: recording paramters
        self.experiments_view = TabView(self.right_panels) 
        self.right_panels.addWidget(self.experiments_view)
               
        self.experimentsInfo_Tab = TableView(self.experiments_view,editable=False, sortable=False, mdcFunc=self.loadSectionImageInfo)
        self.experiments_view.addTab(self.experimentsInfo_Tab, 'Experiments')
        
        self.experimentImage_Tab = TableView(self.experiments_view,editable=False, sortable=False) #, mdcFunc=self.render_sectionImage,kpFunc=self.render_sectionImage)
        self.experimentImage_Tab.currentCellChanged.connect(self.render_sectionImage)
        self.experiments_view.addTab(self.experimentImage_Tab, 'Section images')        
        
        ## bottom panel: extra informations
        self.trees_view = TabView(self.right_panels)
        self.geneInfo_Tab = TableView(self.trees_view,mdcFunc=self.querySelectedGeneInfo) ## Table to show queried gene information
        self.trees_view.addTab(self.geneInfo_Tab,'Gene information' )
        
        self.abAtlas_view = TableView(self.trees_view)
        self.trees_view.addTab(self.abAtlas_view,'Allen brain atlas ontology' ) ## brain regions names
        
        self.right_panels.addWidget(self.trees_view)  

        self.setRightPanelStrechFactor()
        self.frame_splitter.setStretchFactor(0, 8)
        self.frame_splitter.setStretchFactor(1, 1)
        self.frame_splitter.setCollapsible(0, False)
        self.frame_splitter.setCollapsible(1, True)
        
        # add frame_splitter to main layout
        self.MainLayout.addWidget(self.frame_splitter)
          
    def SynAtlas(self):
        print('syn in action')
    
    def setRightPanelStrechFactor(self):
        self.right_panels.setStretchFactor(0, 4)
        self.right_panels.setStretchFactor(1, 5)    
        self.right_panels.setStretchFactor(2,0)
        self.right_panels.setCollapsible(0, True)
        self.right_panels.setCollapsible(1, True)
        self.right_panels.setCollapsible(2, True)
        self.right_panels.setSizes([360, 360, 200])
        
    def add_menubar(self):
        """
        set up manu bar for the main window
        """
        self.mbar = pg.QtGui.QMenuBar()
        self.MainLayout.setMenuBar(self.mbar)       
        self.OptionMenu = self.mbar.addMenu('&Tools')
        self.OptionAction1 = pg.QtGui.QAction("&Pickle current images")
        self.OptionAction1.triggered.connect(self.pickeImage_clicked)
        self.OptionMenu.addAction(self.OptionAction1)
        self.OptionAction2 = pg.QtGui.QAction("&Batch downloading")
        self.OptionAction2.triggered.connect(self.batchDownloading_clicked)
        self.OptionMenu.addAction(self.OptionAction2)
        
        self.HelpMenu = self.mbar.addMenu('&Help')
        self.HelpAction = pg.QtGui.QAction("&LICENSE")
        self.HelpAction.setStatusTip('BSD-3')       
        self.HelpAction.triggered.connect(self.License_clicked)
        self.HelpMenu.addAction(self.HelpAction)
        
    def imageHoverEvent(self, event):
        """Show the position, pixel, and value under the mouse cursor.
        """
        if event.isExit():
            self.Expression_view.view.setTitle("")
            self.ISH_view.view.setTitle("")
            return

        # data = self.ISH_view.image
        if type(self.Expression_view.image)!=type(None):
            if self.ISH_view.view.sceneBoundingRect().contains(event.pos().x(), event.pos().y()):
                ppos = self.ISH_view.view.mapToView(event.pos())           
                i, j = int(ppos.x()), int(ppos.y())
                if i >=0 and j >=0 and i < self.Expression_view.image.shape[1] and j <self.Expression_view.image.shape[0]:               
                    fishVal = self.ISH_view.image[j,i, :]
                    self.ISH_view.view.setTitle(f"pos:({i},{j}), rgb: {fishVal}")
                    if len(self.Expression_view.image.shape)==3:
                        expressionVal = self.Expression_view.image[j, i, :]
                        pixelIntensity = rgb2intensity(np.reshape(expressionVal,(1,1,3)), jetcm)[0][0]
                        rgb_txt = f'RGB: {expressionVal}.'
                        
                    else:
                        rgb_txt = ''
                        pixelIntensity = self.Expression_view.image[j, i]
                        
                    if self.atlas_ImageID!=None:
                        expTitle = f"pos:({i},{j}),"+rgb_txt +f" Expresssion intensity:{pixelIntensity:.2f}. AtlasImage ID: {self.atlas_ImageID}"
                    else:
                        expTitle =  f"pos:({i},{j}),"+rgb_txt +f" Expresssion intensity:{pixelIntensity:.2f}. AtlasImage ID: {self.atlas_ImageID}"
                        
                        
                    self.Expression_view.view.setTitle(expTitle)
                    
                    self.ISH_view_vLine.setPos(ppos.x())
                    self.ISH_view_hLine.setPos(ppos.y())
                    self.ISH_view_vLine.setZValue(1000)
                    self.ISH_view_hLine.setZValue(1000)
                    
                    self.Expression_view_vLine.setPos(ppos.x())
                    self.Expression_view_hLine.setPos(ppos.y())
                    self.Expression_view_vLine.setZValue(1000)
                    self.Expression_view_hLine.setZValue(1000)
                    
                else:
                    self.Expression_view.view.setTitle("")
                    self.ISH_view.view.setTitle("")
        
    def License_clicked(self):
        with open("Data/LICENSE.txt") as f:
            BSD_3 = f.readlines()
        BSD_3b = ""
        for l in BSD_3:
            BSD_3b = BSD_3b+l
        showdialog("This program is under BSD-3 license.\nCopyright (c) 2020-2021, ZeitgeberH@github. All rights reserved.", BSD_3b,False)

    def batchDownloading_clicked(self):
        WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType('Data/uis/batchDownloadOptionsDialog.ui')
        self.bd_dialog  = TemplateBaseClass()
        self.bd_dialog_form = WindowTemplate()
        self.bd_dialog_form.setupUi(self.bd_dialog)
        self.bd_dialog_form.fileButton.clicked.connect(self.batchGenesFile)
        self.bd_dialog_form.buttonBox.clicked.connect(self.buttonBoxResponse)
        self.bd_dialog.exec_()

    def getBatchDownloadInfo(self,showDialog=False):
        fullResolution =self.bd_dialog.children()[2].children()[0].children()[1].isChecked()
        Downsampled =self.bd_dialog.children()[2].children()[0].children()[2].isChecked()
        downSampledFactor = self.bd_dialog.children()[2].children()[0].children()[0].children()[4].value()
        QFactor = self.bd_dialog.children()[2].children()[0].children()[0].children()[5].value()
        
        for j in range(3):
            if self.bd_dialog.children()[2].children()[1].children()[j].isChecked():
                planeOption = j
                break
        userDefinedGenesFiles = self.bd_dialog_form.filePath.text()
        with open(userDefinedGenesFiles,'r') as file:
            geneNames = file.read().splitlines()
        geneNames = [g for g in geneNames if len(g)>0]
        self.BatchDownloading_FileInfo(geneNames, planeOption, fullResolution, downSampledFactor, QFactor,showDialog)
                
    def buttonBoxResponse(self,evt):
        '''
        if 'OK' button pressed, parse values from the UI

        Parameters
        ----------
        evt : TYPE: buttonBox responses
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        if evt.text() == 'OK':
            if self.bd_dialog_form.filePath.text()!='':
                self.getBatchDownloadInfo(True)
            else:
                showdialog("Press 'Locate file' button to specify the file with names of genes you would like to know!")
            
    def BatchDownloading_FileInfo(self, geneNames, planeOption, fullResolution, downSampledFactor, QFactor, showDialog=False):
        p = ['all','coronal','sagittal']
        if fullResolution:
            ds = f'\nDownloading {p[planeOption]} planes section images for following genes at full resolution:'
        else: 
            ds = f'\nDownloading {p[planeOption]} section images at downsampling factor:{downSampledFactor}, quality:{QFactor} for following genes at reduced of {planeOption}: '
        print(ds)
        print(f'{geneNames}')

        geneDict = {}
        inforStr = 'Gene    Experiments\n'
        for g in geneNames:
            geneExpID_all, geneExpData_all = self.geapi.get_gene_experiments2(g)
            geneExpID = []
            geneExpData = []
            if geneExpID_all!=None:
                if planeOption!=0:
                    for gid, ged in zip(geneExpID_all,geneExpData_all):
                        if ged['plane_of_section_id']== planeOption: ## filter specified image planes!
                            geneExpID.append(gid)
                            geneExpData.append(ged)
                else:
                    geneExpID = geneExpID_all
                    geneExpData = geneExpData_all
                    
                inforStr = inforStr+ f'{g:8}  {len(geneExpID)}\n'
                geneDict[g]={'geneExpID': geneExpID, 'geneExpData':geneExpData}
        
        returnValue = showdialog(ds, details= inforStr)
            
        if returnValue == QtWidgets.QMessageBox.Ok:
            self.BatchDownload(geneDict, fullResolution, downSampledFactor, QFactor)
        
    def BatchDownload(self,geneDict, imageQuality, downSampleFactor, qualityFactor):
        Expression=False
        with pg.ProgressDialog("Downloading FISH images for gene..." , maximum=len(geneDict),\
                               busyCursor=True, nested=True) as dlg_gene:              
            for gene in geneDict:
                expIDs = geneDict[gene]['geneExpID']
                with pg.ProgressDialog("Experiments..." , maximum=len(expIDs),\
                                       busyCursor=True, nested=True) as dlg_exp:                           
                    for expID in expIDs:
                        cache = check_gene_cached_images(self.geapi.fish_images_cache, gene, expID,Expression, imageQuality, downSampleFactor, qualityFactor)
                        #### check cache and download if not availabe
                        if not cache:
                            imageIDs, sectionNumbers, expMetaData = self.geapi.get_gene_experiments_imageList(expID) ## query imageIDs
                            with pg.ProgressDialog("... ISH & Expression" , maximum=2,\
                                                   busyCursor=True, nested=True) as dlg_eh:
                
                                for Expression in [True, False]: ## expression images
                                    with pg.ProgressDialog("...... Images..." , maximum=len(imageIDs),\
                                                           busyCursor=True, nested=True) as dlg1:
                                        for j, sn in zip(imageIDs, sectionNumbers):
                                            self.geapi.download_gene_section_imageData(gene, expID, j, sn, Expression, imageQuality, downSampleFactor,qualityFactor)
                                            dlg1 += 1
                                            if dlg1.wasCanceled():
                                                print("Canceled stage")
                                                break
                                    dlg_eh += 1
                                    if dlg_eh.wasCanceled():
                                        print("Canceled images stage")
                                        break
                        dlg_exp += 1
                        if dlg_exp.wasCanceled():
                            print("Canceled experiment stage")
                            break
                                        
                dlg_gene += 1
                if dlg_gene.wasCanceled():
                    print("Canceled gene stage")
                    break
                    
    def batchGenesFile(self):        
        fileName = getfiles()
        self.bd_dialog_form.filePath.setText(fileName)
        # self.getBatchDownloadInfo(showDialog=True)


    def add_toolbar(self):
        self.toolbar = QtWidgets.QToolBar("Main toolbar")
        self.addToolBar(2, self.toolbar) # https://doc.qt.io/qt-5/qt.html#ToolBarArea-enum
        self.toolbar.addSeparator()
        self.toolbar.addSeparator()
        # self.tooglesettingsAction = pg.QtGui.QAction(pg.QtGui.QIcon("Data/icons/lsm.png"), "Open LSM images") 
        # self.tooglesettingsAction.triggered.connect(self.openLSM_clicked)     
        # self.toolbar.addAction(self.tooglesettingsAction) 
        
    def sycImageToRegion(self,expID, regionName):
        ##TODO
        try:
            imgData_meta = self.geapi.getExpTargetXY(expID, regionName)
            print('current region:  '+ regionName)

            imgPath = os.path.join(self.currentImageCacheDir, str(imgData_meta['section_number']) + '_'+str(imgData_meta['section_image_id']))
            imgPath_E = imgPath+'_E.jpg'
            imgPath_H = imgPath+'_H.jpg'            
            self.updateFishImage(self.ISH_view, imgPath_H)
            self.updateFishImage(self.Expression_view, imgPath_E, True)
            row = self.currentSearchISH_imageLst.index(imgPath_H)
            self.experimentImage_Tab.setCurrentCell(row, 0)               
            self.experimentImage_Tab.setFocus()
            return True
        
        except:
            showdialog('Syncronization failed!')
            return False


                    
    def event_parameters_stateChange(self, params, changes):
        '''
        Responding to user's input in the parameter tree

        Parameters
        ----------
        params : TYPE
            DESCRIPTION.
        changes : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        pv = self.globalPars.getValues()
        for param, change, data in changes:
            childName = param.name()
            if childName == 'Gene':
                self.experiments_view.setCurrentIndex(0)
                
            if childName=='region':
                if self.currentSearchDict!=None:
                    self.currentSearchDict['region'] = data
            if (childName == 'Sync') and (self.currentExpID!=None):
                regionName = pv['Gene & Brain region & species'][1]['Region'][0]
                if regionName!='':
                    print(f'Sync to region {regionName}')
                else:
                    showdialog("Specifiy a rgion to syn!")
            if (childName == 'Syn to atlas after double clicking'):
                if not data:
                    self.clearSVGItems() ## turn off atlas overlay
            if childName=='Expression':
                if data=='color mask':
                    self.updateFishImage(self.Expression_view,self.Expression_mask_file, True, False)
                else:
                    self.updateExpressionIntensity()
            if childName == 'atlas section offset':
                # pdb.set_trace()
                if self.atlasID!=None:
                    idx = self.svgItem_atlasE_index
                    atlas_ofs_idx = pv['image-to-Atlas'][1]['atlas section offset'][0] + idx
                    if atlas_ofs_idx>-1 and atlas_ofs_idx<len(self.atlas_sectionID_list):
                        atlas_imageID = self.atlas_sectionID_list[atlas_ofs_idx]                        
                        atlas_fn  = self.geapi.getAtlasBoundaryImage(atlas_imageID, self.atlasID, GraphicGroupLabel_id=28)
                        xP, yP = self.svgItem_atlasE.pos().x(),self.svgItem_atlasE.pos().y()
                        self.overlayAtlasSVG(xP+self.atlas_x,yP+self.atlas_y,self.atlas_x,self.atlas_y, atlas_fn)
                        self.atlas_ImageID = atlas_imageID
                else:
                    print('Double clicking to syn an atlas image first!')

            if childName == 'x offset':
                # pdb.set_trace()
                x, y = self.svgItem_atlasE_pos
                y_offset = pv['image-to-Atlas'][1]['y offset'][0]
                self.svgItem_atlasE.setPos(x+data, y+y_offset)
                # self.svgItem_atlasH.setPos(x+data, y+y_offset)
                
            if childName == 'y offset':
                # pdb.set_trace()
                x_offset = pv['image-to-Atlas'][1]['x offset'][0]
                x, y = self.svgItem_atlasE_pos
                self.svgItem_atlasE.setPos(x+x_offset, y+data)
                # self.svgItem_atlasH.setPos(x+x_offset, y+data)
            
            if childName == 'Expression mask opacity':
                opc = pv['image-to-Atlas'][1][childName][0]
                self.Expression_view.view.setOpacity(opc)                
          
        if (len(changes) == 1) and (childName == 'Sync') and (self.currentExpID!=None):
            if self.geapi!=[]:
                regionName = pv['Gene & Brain region & species'][1]['Region'][0]
                self.sycImageToRegion(self.currentExpID, regionName.upper())
            else:
                showdialog("Load an experiment first!")
            return
                
        #
        search_pars = {}
        search_pars['gene'] = pv['Gene & Brain region & species'][1]['Gene'][0]
        if search_pars['gene']=='':
            return
        search_pars['region'] = pv['Gene & Brain region & species'][1]['Region'][0]
        search_pars['species'] = pv['Gene & Brain region & species'][1]['Species'][0]
        search_pars['plane'] = pv['View options'][1]['Plane'][0]
        search_pars['Syn to atlas after double clicking'] = pv['image-to-Atlas'][1]['Syn to atlas after double clicking'][0]
        search_pars['maximal resolution'] = pv['View options'][1]['Maximal resolution'][0]
        search_pars['downsample'] = pv['View options'][1]['Image quality at lower resolution'][1]['Downsample'][0]
        search_pars['quality'] = pv['View options'][1]['Image quality at lower resolution'][1]['Downsample'][0]
        self.updateGeneOfInterest(search_pars)
            
    def resetMergeImage(self):
        for item in self.lsm_bottomRight_view.view.items:
            if hasattr(item, 'Name'):
                if item.Name[:3]!='br_':                 
                    self.lsm_bottomRight_view.view.removeItem(item)
        self.lsm_bottomRight_view.setImage(self.mergeImage, axes={'x':0,'y':1,'c':2,'t':None},autoRange=False,autoLevels=False)
              
    def overlay_segmentation(self, bds, imageViewHandle, noMask=False):
        '''
        Overlay segmented ROIs to given ImageView instance

        Parameters
        ----------
        bds : TYPE: ndarray
            contour boundaries of segemented results
        imageViewHandle : pyqtgraph.imageView instance
            imageview window to draw

        Returns
        -------
        None.
        '''
        for item in imageViewHandle.view.items:
            if hasattr(item, 'Name'):
                if item.Name[:3]=='sgb':                 
                    imageViewHandle.view.removeItem(item)
        if noMask: ## no mask
            return
        if len(bds.shape)==2:
            imItem = pg.ImageItem(bds)
            imItem.Name = 'sgb'
            imageViewHandle.view.addItem(imItem)
            imItem.setZValue(100)
        else:
            allChans = set([0,1,2])
            for nMask in range(bds.shape[2]):
                mask_ =  np.empty_like(bds)
                mask_[:,:,nMask] = bds[:,:,nMask]
                otherChans = list(allChans.difference([nMask]))
                for j in otherChans:
                    mask0 =  np.empty_like(bds[:,:,0])
                    mask0[bds[:,:,nMask]>0] = 80 ## set other channel's pixels intensity. so we can get a color image
                    mask_[:,:,j] = mask0
                imItem = pg.ImageItem(mask_)
                imageViewHandle.view.addItem(imItem)
                imItem.Name = 'sgb'+str(nMask)
                imItem.setZValue(100)
            
    def initGeneAPI(self):
        try:
            self.geapi = GeneRenderAPI()
            self.geapi.DOWNSAMPLE_FACTOR = 4 ## downsampling factor if not Full image, range(0 to 10), higher the lower quality
            self.geapi.QUALITY_FACTOR = 50 ## quality factor if not Full image, range(0 to 100), higher the better quality
            return True
        except:
            self.geapi = []
            print('Error to initiate gene query API. Try again')
            return False

    def loadSectionImageInfo(self,row, col):
        '''
        load section image metainformation       

        Parameters
        ----------
        row : TYPE
            DESCRIPTION.
        col : TYPE
            DESCRIPTION.

        Returns
        -------

        '''
        expID = self.experimentsInfo_Tab.item(row, 0).value
        expPlane = self.experimentsInfo_Tab.item(row, 1).value
        expGene = self.experimentsInfo_Tab.item(row, 2).value


        if self.geapi==[]:
            r=self.initGeneAPI()
            if not r:
                return
        ## query and return section image list from this experiment
        ##TODO
        ## CHECK cached image file first.If no, then request from ABA
        imageIDs, secNums, imageMeta = self.geapi.get_gene_experiments_imageList(expID)
        expDF = pd.DataFrame()
        expDF.loc[:,'Section'] = secNums
        expDF.loc[:,'imageID'] = imageIDs
        expDF['Experiment'] = expID
        expDF['Gene'] = expGene
        expDF['Plane'] = expPlane
        expDF['Height'] = [e['image_height'] for e in imageMeta]
        expDF['Width'] = [e['image_width'] for e in imageMeta]
        expDF['Resolution'] = [e['resolution'] for e in imageMeta]
        expDF['Failed'] = [e['failed'] for e in imageMeta]
        expDF = expDF.set_index('Section')
        if expPlane == 'coronal':
            expDF.sort_index(ascending=False, inplace=True)
        data = np.array(expDF.to_records(),
                              dtype=[('Section',object), ('imageID', object), ('Experiment', object), ('Gene', object),\
                                     ('Plane', object),('Height', object),('Width', object),('Resolution', object),\
                                        ('Failed', object) ])
        self.experimentImage_Tab.setData(data)
        self.experiments_view.setCurrentIndex(1)
        self.currentExpID = expID
        # self.currentSearchISH_imageLst = imageFilesList_H 
        # self.currentSearchExpression_imageLst = imageFilesList_E
        self.loadExp_image(expGene, expID, expPlane,0)

    def render_sectionImage(self,row, col, previousRow, previousColumn):
        '''
        render section images in imageItems panels

        Parameters
        ----------
        row : TYPE
            DESCRIPTION.
        col : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        if self.currentSearchISH_imageLst !=None:
            if self.currentISH_image != self.currentSearchISH_imageLst[row]:
                self.clearSVGItems()
                self.updateFishImage(self.ISH_view, self.currentSearchISH_imageLst[row])
                self.Expression_mask_file = self.currentSearchExpression_imageLst[row]
                self.Expression_intensityVal = None
                self.updateFishImage(self.Expression_view, self.currentSearchExpression_imageLst[row], True)
                self.Expression_view.view.setXLink(self.ISH_view.view)
                self.Expression_view.view.setYLink(self.ISH_view.view)
                self.currentISH_image = self.currentSearchISH_imageLst[row]
                self.statusBar.showMessage(self.currentISH_image , 5000)
        
    def loadExp_image(self, gene, expID,plane, idx=0):
        '''
        Load FISH image list for current gene and experiment ID with image index at index into table
        If not in cache, download it from ABI.

        Parameters
        ----------
        gene : TYPE, str
            Gene symbol.
        expID : TYPE, Int
            Experiment ID
        idx : TYPE, optional
            Section index to load for current experiment. The default is 0.
        fullImage : TYPE, bool
            DESCRIPTION. whether to download the full resolution images. Would take times if set true (and possible connection failure)

        Returns
        -------
        None.

        '''
        imageQuality = self.currentSearchDict['maximal resolution']
        downSampleFactor =self.currentSearchDict['downsample']
        qualityFactor = self.currentSearchDict['quality']
        Expression=False
        ##TODO
        ## check by image IDs. (for now, it just check if there is ANY images)
        cache = check_gene_cached_images(self.geapi.fish_images_cache, gene, expID,Expression, imageQuality, downSampleFactor, qualityFactor)
        #### check cache and download if not availabe
        if not cache:
            self.statusBar.showMessage(f'Downloading data for gene {gene}. Please wait...', 3000)
            imageIDs, sectionNumbers, expMetaData = self.geapi.get_gene_experiments_imageList(expID) ## query imageIDs
            with pg.ProgressDialog("Downloading FISH images..." , maximum=2,\
                                   busyCursor=True, nested=True) as dlg0:

                for Expression in [True, False]: ## expression images
                    with pg.ProgressDialog("Downloading ..." , maximum=len(imageIDs),\
                                           busyCursor=True, nested=True) as dlg1:
                        for j, sn in tqdm.tqdm(zip(imageIDs, sectionNumbers), total=len(imageIDs)):
                            self.geapi.download_gene_section_imageData(gene, expID, j, sn, Expression, imageQuality, downSampleFactor,qualityFactor)
                            dlg1 += 1
                            if dlg1.wasCanceled():
                                print("Canceled stage")
                                break
                        dlg0 += 1
                        if dlg0.wasCanceled():
                            print("Canceled stage")
                            break
                else:        
                    cache=True
            
        if cache: ## load cached data
            imageFilesList_H = self.geapi.get_gene_experiments_imageData(gene, expID,True, False, imageQuality,downSampleFactor, qualityFactor)## ISH IMAGES
            # self.updateFishImage(self.ISH_view, imageFilesList_H[idx])
            imageFilesList_E = self.geapi.get_gene_experiments_imageData(gene, expID,True, True,imageQuality,downSampleFactor, qualityFactor)## EXPRESSION MASKS
            # self.updateFishImage(self.Expression_view, imageFilesList_E[idx])
            if plane == 'coronal':
                imageFilesList_H  = imageFilesList_H [::-1]
                imageFilesList_E = imageFilesList_E[::-1]
            self.currentSearchISH_imageLst = imageFilesList_H 
            self.currentSearchExpression_imageLst = imageFilesList_E
            self.currentImageCacheDir = str(Path(imageFilesList_E[0]).parent)
            region = self.globalPars.getValues()['Gene & Brain region & species'][1]['Region'][0] 
            self.sycImageToRegion(self.currentExpID, region.upper())
        else:
            print('Glitches in cached FISH images...')

    def clearSVGItems(self):
        self.atlas_ImageID = None
        for item in self.Expression_view.view.items:
            if type(item) == QtSvg.QGraphicsSvgItem:
                self.Expression_view.view.removeItem(item)
                item.setParent(None)
                del item
        for item in self.ISH_view.view.items:
            if type(item) == QtSvg.QGraphicsSvgItem:
                self.ISH_view.view.removeItem(item)
                item.setParent(None)
                del item
                                        
    def otherInitStuff(self):
        self.initGeneAPI()
        self.initABAtlasTable()
        self.currentSearchDict = None
        self.currentSearchISH_imageLst = None
        self.currentSearchExpression_imageLst = None
        self.currentExpID = None
        self.currentISH_image = None
        self.defaultBrowser = None
        self.atlas_sectionID_list = None
        self.atlas_x = None
        self.atlas_y = None
        self.atlasID = None
        self.atlas_ImageID = None
        self.Expression_mask = None
        self.Expression_intensity = None
        ## LSM file related
        # self.currentLSMfile = None
        # self.currentLSMimage = None
        # self.segs_merge = {'R':None,'G':None,'B':None} ## store segmented value right before postprocessing steps
        # self.merge_r = True
        # self.merge_g = True
        # self.merge_b = True
        # self.DNNModel = None ## deal with tf.function retracing?
        
    def initABAtlasTable(self):
        ## load Allen Brain isntitute atlas acronym
        df = self.geapi.atlas_df[['id','acronym','name']]
        self.abAtlas_view.setData(df.to_records())        

    def getSynedAtlas(self,atlasID, imageID, x,y, groupID=28):
        imgSynDict = self.geapi.synImageToAtlas(atlasID, imageID, (x, y)) ## get  dictionary of syned atlas Image information
        #imgSynDict['section_image_id'],imgSynDict['x'],imgSynDict['y']
        if len(imgSynDict)>0:
            atlas_imageID, xa,ya = imgSynDict['section_image_id'], imgSynDict['x'],imgSynDict['y']
            print(f'Syn info: section imgeID {imageID}, pos:({x},{y}). Atlas imageID {atlas_imageID}, pos: ({xa},{ya})')
            fileName  = self.geapi.getAtlasBoundaryImage(atlas_imageID, atlasID, GraphicGroupLabel_id=groupID)
            return fileName, atlas_imageID,xa, ya
        return None
            
    def overlayAtlasSVG(self, x, y, xa,ya, atlas_fileName):
        '''
        Parameters
        ----------
        x : TYPE, float
            current section image point's x
        y : TYPE, float
            current section image point's x
        xa : TYPE, float
            atlas image synced X
        ya atlas_fileName : TYPE
            atlas image synced Y
        Returns
        -------
        None.

        '''
     
        self.clearSVGItems()
        self.svgItem_atlasE = QtSvg.QGraphicsSvgItem(atlas_fileName)

        self.svgItem_atlasE.setOpacity(1);#0.08)
        self.svgItem_atlasE.scale(1, 1)   
        self.svgItem_atlasE.setPos(x-xa, y-ya)
        self.svgItem_atlasE_pos = (x-xa, y-ya) ## for reference
        self.Expression_view.view.addItem(self.svgItem_atlasE)
        self.svgItem_atlasE.setZValue(-100)
                                            
    def ISH_view_mouseClicked(self, event):
        if type(self.ISH_view.image)!=type(None):
            if self.currentSearchDict['Syn to atlas after double clicking']: ## only syn when True
                if event.button()==1 and event.double(): ## if left double clicking~                 
                    if self.ISH_view.view.sceneBoundingRect().contains(event.pos().x(), event.pos().y()):
                            if self.currentSearchDict['species']=='Mouse':
                                row = self.experimentImage_Tab.currentRow()
                                imageID = self.experimentImage_Tab.item(row, 1).value ## get current section image ID
                                planeName = self.experimentImage_Tab.item(row, 4).value
                                if planeName == 'coronal':
                                    atlasID = 1
                                else:
                                    atlasID = 2
                                x = self.ISH_view_vLine.pos().x() ## position value from event is not consistent with hover event!!!
                                y = self.ISH_view_hLine.pos().y()
                                
                                atlas_fn,atlas_imageID, xa,ya = self.getSynedAtlas(atlasID, imageID, x, y)
                                self.overlayAtlasSVG(x,y,xa,ya, atlas_fn)
                                ##TODO
                                if self.atlas_sectionID_list==None:
                                    self.atlas_sectionID_list, _,_ = self.geapi.getAtlasImagesAnnoationList(atlasID)
                                self.svgItem_atlasE_index  = self.atlas_sectionID_list.index(atlas_imageID) ## store it for manul correction
                                self.atlas_x = xa ## same above
                                self.atlas_y = ya
                                self.atlasID = atlasID
                                self.atlas_ImageID = atlas_imageID
                            else:
                                print('Not supported yet for human data')
        
    def updateExperimentsTable(self, expDF):    
        data = np.array(expDF.to_records(),
                              dtype=[('Experiment', object),('Plane', object), ('Gene symbol', object),\
                                     ])
        self.experimentsInfo_Tab.setData(data)
             
    def updateGeneOfInterest(self, parDict):
        if self.geapi==[]:
            r=self.initGeneAPI()
            if not r:
                return
        if parDict['species']=='Mouse':
            parDict['gene'] = parDict['gene'].capitalize()
        else:
            parDict['gene'] = parDict['gene'].upper()
        self.updateGenessTable(parDict['gene'],  parDict['species'])
        try:
            geneExpID, geneExpData = self.geapi.get_gene_experiments2(parDict['gene'])
        except:
            print('No respondes from Server!')
            return
        if geneExpID != None:
            expDF = pd.DataFrame()
            expDF.loc[:,'Experiment'] = geneExpID
            expDF.loc[:,'Plane'] = [self.geapi.sectionPlaneDict[d['plane_of_section_id']] for d in geneExpData]
            expDF['Gene Symbol'] = parDict['gene']
            # expDF.loc[:,'Expression Summary'] = [] ## TODO
            expDF = expDF.set_index('Experiment')
            if parDict['plane']!='All':
                expDF=expDF[expDF['Plane']==parDict['plane'].lower()]
            self.updateExperimentsTable(expDF)
            self.statusBar.showMessage(f'Data for gene {parDict["gene"]} is listed!', 10000)
            self.currentSearchDict = parDict ## updating state variables!
            self.geapi.DOWNSAMPLE_FACTOR = parDict['downsample']
            self.geapi.QUALITY_FACTOR = parDict['quality']
        else:
            self.statusBar.showMessage(f'Data for gene {parDict["gene"]} is not found', 10000)
 
    def updateExpressionIntensity(self):
        # self.Expression_mask
        with pg.BusyCursor():
            if self.Expression_intensityVal is None:
                if self.Expression_view.image is None:
                    img = np.array(Image.open(self.Expression_mask_file))                            
                else:
                    img = self.Expression_view.image
                self.Expression_intensityVal = rgb2intensity(img)
                self.Expression_view.setImage(self.Expression_intensityVal, pos = [0,0], axes={'x':1,'y':0,'c':None,'t':None},autoRange=False)
            else:
               self.Expression_view.setImage(self.Expression_intensityVal, pos = [0,0], axes={'x':1,'y':0,'c':None,'t':None},autoRange=False) 

                
    def updateFishImage(self, imageItemHandle, imgFileName, Expression=False, autoRange=True):
        '''set image for imageView
        '''
        # imgFileName = base_dir / 'FishImagesCache' / 'Stac-69887313' / 'DQ' /'65_69866246_H.jpg'
        # img = np.transpose(Image.open(imgFileName),axes=(1,0,2))
        img = np.array(Image.open(imgFileName))
        imageItemHandle.clear()
        # w,h,_ = img.shape
        # imageItemHandle.setImage(img, pos = [h/2,w/2])
        imageItemHandle.setImage(img, pos = [0,0], axes={'x':1,'y':0,'c':2,'t':None},autoRange=autoRange)
        if Expression:
            imageItemHandle.view.setOpacity(0.85)            
        # imageItemHandle.setOpts(update=True, opacity=0.2)
        
    def findGenesInTable(self,geneSymbol):
        nrows = self.geneInfo_Tab.rowCount()
        for r in range(nrows):
            if geneSymbol==self.geneInfo_Tab.item(r,0).value:
                return True
        else:
            return False
        
    def updateGenessTable(self, geneSymbol, species):
        ##TODO
        #### CHECK if gene is in cache
        if not self.findGenesInTable(geneSymbol):
            geneInfo = queryGene_ncbi(geneSymbol, species)
            if len(geneInfo) > 0:
                expDF = pd.DataFrame([geneInfo])
                expDF['symbol']= geneSymbol
                expDF['genecard url'] = 'genecards.org/cgi-bin/carddisp.pl?gene='+geneSymbol
                expDF = expDF[['symbol','description','type',\
                       'entrez_id','chromosomes', 'number of transcripts','genecard url']]
                data = np.array(expDF.to_records(index=False),
                                      dtype=[('symbol', object),('description', object), ('type', object),\
                                             ('entrez_id', object),('chromosomes', object),('number of transcripts', object),('genecard url', object)
                                             ])
                self.geneInfo_Tab.appendData(data)

    def querySelectedGeneInfo(self, row, col):
        if col== self.geneInfo_Tab.columnCount()-1:  ## last col
            url =  self.geneInfo_Tab.item(row, col).value  
            self.openGeneCardPage(url)
            
    def openGeneCardPage(self, url):
        if self.defaultBrowser==None: ##  if default browser not set, try get default          
            try:
                from gr_allenFISH_utils import getDefaultBrowerPath
                browserPath = getDefaultBrowerPath().replace('\\','/') +' %s'
                self.defaultBrowser = browserPath
                webbrowser.get(using = self.defaultBrowser).open_new(url)
            except:
                webbrowser.open(url,new=2) ## fallback to IE if failed
        elif self.defaultBrowser == -1: ## failed already
            webbrowser.open(url,new=2) ## fallback to IE
        else: ## properly set
            webbrowser.get(using = self.defaultBrowser).open_new(url)

    def pickeImage_clicked(self):
        if self.currentExpID is not None:
            import pickle
            with pg.BusyCursor():
                fname = self.currentSearchDict['gene']+'_'+self.currentISH_image+'_ExpID-'+str(self.currentExpID)
                fishVal = self.ISH_view.image
                expressionVal = self.Expression_view.image
                with open(fname+'_expression.pickle','wb') as handle:
                    pickle.dump(expressionVal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
                with open(fname+'_fish.pickle','wb') as handle:
                    pickle.dump(fishVal,handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            print('No image loaded yet!')
                
    def closeFISH_ViewItems(self):
        for view in [self.ISH_view,self.Expression_view]:
            view.clear()
            view.close()
            
    def closeEvent(self, event):       
        reply = QMessageBox.question(self, 'Closing Window', 'Are you sure you want to close the window?',\
                                     QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            try: ## release resources
                self.closeFISH_ViewItems()
            except:
                print('clean up failed!')
            event.accept() 
        else:
            event.ignore()     
def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow(app)
    main.show()
    sys.exit(app.exec_())
  
if __name__ == '__main__':
    main()