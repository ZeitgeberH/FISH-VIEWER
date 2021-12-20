# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:24:07 2021

@author: MHu
"""
from PyQt5 import QtGui,QtWidgets, QtCore
from PyQt5.QtWidgets import QMessageBox, QTableWidgetItem
import pyqtgraph as pg
from pyqtgraph import GraphicsLayoutWidget
from pyqtgraph.parametertree import Parameter, ParameterTree
from gr_utilities._buildParamTypes import makeAllParamTypes
import os

class FileModel(QtWidgets.QFileSystemModel):
    """
    Class for file system model
    """
    def __init__(self, root_path):
        super(FileModel, self).__init__()

        # hide system files
        self.setFilter(QtCore.QDir.AllDirs |
                       QtCore.QDir.NoDotAndDotDot |
                       QtCore.QDir.AllEntries)

        # filter out non dats and disable showing
        self.setNameFilters(['*.lsm'])
        self.setNameFilterDisables(False)

        # set root
        self.setRootPath(root_path)
        
class FileView(QtWidgets.QTreeView):
    """
    Class for view of file system model
    """
    def __init__(self, parent, root_path,excute_func=None):
        super(FileView, self).__init__()
        # set model
        self.model = FileModel(root_path)
        self.setModel(self.model)
        self.excute_func = excute_func ## main function to excute
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        # self.clicked.connect(self.clickedItem)       
        # self.myMenu = QtWidgets.QMenu('Menu', self)
        # self.addFileAction = pg.QtGui.QAction("Import .dat/.abf files for this slice")
        # self.addFileAction.triggered.connect(self.importAllDats_clicked)
        # self.myMenu.addAction(self.addFileAction)
        self.currentPath = root_path
        # set root
        self.setRootIndex(self.model.index(root_path))

        # hide unnecessary columns
        self.hideColumn(1)
        self.hideColumn(2)
        self.hideColumn(3)

        # bind event
        self.doubleClicked.connect(self.on_double_click)
        self.customContextMenuRequested.connect(self.openContextMenu)
                                   
        
    def openContextMenu(self, index):
        self.myMenu.exec_(QtGui.QCursor.pos())
        
    @QtCore.pyqtSlot(QtCore.QModelIndex)
    def on_double_click(self, index):
        """
        Event for changing .lsm file
        :param index:
        :return:
        """
        if self.excute_func!=None:
            model_index = self.model.index(index.row(), 0, index.parent())
    
            # get file from index
            file_path = os.path.abspath(self.model.filePath(model_index))
            _, f = os.path.split(file_path)
    #        os.chdir(self.frame.root)
            # check extension
            _, ext = os.path.splitext(file_path)
            if ext in ['.lsm']:
                self.excute_func(file_path) ## excute function with file_path as paraemter
    
class TabView(QtWidgets.QTabWidget):
    """
    Class for tab view.
    """
    def __init__(self, parent, enableContextMenu = False):
        super().__init__()
        self.frame = parent
        self.enableContextMenu = enableContextMenu
        if enableContextMenu:
            self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
            self.customContextMenuRequested.connect(self.openContextMenu)
            self.actionItems = []
            
    def openContextMenu(self):
        try:
            myMenu = QtWidgets.QMenu('Menu', self)
            myMenu.addAction(self.actionItems[idx])
            myMenu.exec_(QtGui.QCursor.pos())
        except IndexError:
            print('action item not set for this tab')

    def insertTabWithAction(self, tabObj, actions, label='' ):
        ## inherit parent addTab function
        if isinstance(actions, list):
            for a in actions:
                assert isinstance(a, pg.QtGui.QAction), "second argument should be list of instance of QAction"
        else:
            assert isinstance(actions, pg.QtGui.QAction), "second argument should be instance of QAction"
        if self.enableContextMenu:
            super().addTab(tabObj, label)
            self.actionItems.append(actions)
        else:
            raise AttributeError("Set enableContextMenu as True before using this method")
            
class PlotView(GraphicsLayoutWidget):
    """
    Class for plot view.
    """
    def __init__(self, parent, title = '', background = 'k'):
        super(PlotView, self).__init__()
        self.frame = parent
        self.setWindowTitle(title)
#        self.setBackground(background)
#        self.layout = pg.GraphicsLayoutWidget(title = title)
#        self.layout.ci.setBorder()
#        self.layout.ci.setBorder(pg.mkPen(None))
#        self.layout.ci.setContentsMargins(10, 10, 10, 20)
#        self.setCentralWidget(self.layout.ci)
        self.setBackground(background)
#        self.layout = pg.GraphicsLayoutWidget(title = title)
        self.ci.setBorder()
        self.ci.setBorder(pg.mkPen(None))
        self.ci.setContentsMargins(10, 10, 10, 20)
        self.setCentralItem(self.ci)
        self.show()
        
    def refresh(self):
        self.ci.setContentsMargins(10, 10, 10, 20)
            
class ParameterView():
    """
    Class for enhance pyqtgraph.parametertreeParameter class
    """
    def __init__(self, parent, name='params', type='group', children=None, readonly=False):
        self.parent = parent
        self.params = children  ## actual paratmers
        self.state = [] ## save and restore
        self.p = Parameter.create(name= name, type= type, children= children, readonly=readonly)
        #self.p.sigTreeStateChanged.connect(self.change)
#        for child in self.p.children():
#            child.sigValueChanging.connect(self.valueChanging)
#            for ch2 in child.children():
#                ch2.sigValueChanging.connect(self.valueChanging)
        #self.p.param('Save/Restore parameters', 'Save State').sigActivated.connect(self.save)
        #self.p.param('Save/Restore parameters', 'Restore State').sigActivated.connect(self.restore)

    def setChildValue(self, childPath, childValue):
        '''
        Set value of a child
        ----------
        childPath : list of string
            the full path for this child
        childValue : 
        Returns
        -------
        None.
        '''
        vals = self.p.getValues()
        v = None
        depths = 0
        expression = 'vals'
        while v == None:
            child = vals[childPath[depths]]
            v = child[0]
            if v== None:
                index = '[1]'
            else:
                index = '[0]'
            expression = ''.join('['+childPath[depths]+']'+index)
            depths += 1
        eval(expression)
        self.p.setValue(vals)
        
    def change(self, param, changes):
        for param, change, data in changes:
            path = self.p.childPath(param)
            if path is not None:
                childName = '.'.join(path)
            else:
                childName = param.name()
            print('  Parameter: %s'% childName)
            print('  Value:      %s'% str(data))
            print('  ----------')
        self.parent.rePlotSpikePanels()
                
    def valueChanging(param, value):
        print("Value changing (not finalized): %s %s" % (param, value))
      
    def save(self):
        self.state = self.p.saveState()
        
    def restore(self):
        add = self.p['Save/Restore parameters', 'Restore State', 'Add missing items']
        rem = self.p['Save/Restore parameters', 'Restore State', 'Remove extra items']
        self.p.restoreState(self.state, addChildren=add, removeChildren=rem)
        
class TreeView(pg.QtGui.QTreeWidget):
    """
    Class for viewing tree of pul file
    """
    def __init__(self, parent):
        super(TreeView, self).__init__()
        self.frame = parent
        self.setHeaderLabels(['Node', 'Label'])
        self.setColumnWidth(0, 200)
        # allow multi selection
        self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
                
            
    def update_tree_recursive(self, root_item, index, dat_index = None):
        """Recursively read nested ABA Json file
        and add items into the GUI tree to allow browsing.
        """        
        root = self.bundle.pul
        self.filetype = '.dat'
        node = root
        if node == None:
            return
        for i in index:
            node = node[i]
        node_type = node.__class__.__name__
        if node_type.endswith('Record'):
            node_type = node_type[:-6]
        try:
            if node_type[:2] == 'V9':
                node_type += str(getattr(node, node_type[3:] + 'Count'))
            else:
                node_type += str(getattr(node, node_type + 'Count'))
        except AttributeError:
            pass
        
        try:
            node_label = node.Label
        except AttributeError:
            if node_type[:5]=='Pulse':
                node_label = self.bundle.header.Version
            else:
                node_label = ''
        if node_type[:2] == 'V9':
            item = pg.QtGui.QTreeWidgetItem([node_type[3:], node_label])
        else:
            item = pg.QtGui.QTreeWidgetItem([node_type, node_label])
        if len(index) == 2:
            self.seriesNode.append(node)
            
        if dat_index == None: ## update all
            root_item.addChild(item)
            item.node = node
            item.index = index

                #print(self.seriesNode)
            # if len(index) == 3:
            #     self.sweepCount +=1
            # elif len(index) ==4:
            #     self.traceCount +=1
            if len(index) < 2:
                item.setExpanded(True)
            for i in range(len(node.children)):
                self.update_tree_recursive(item, index + [i])
        else:
            if len(index)==2:  ## if this series in the selection list!
                self.seriesNode.append(node)
                if index == dat_index:  ## only add this series
        
                    root_item.addChild(item)
                    item.node = node
                    item.index = index
                    for i in range(len(node.children)):
                                self.update_tree_recursive(item, index + [i])       
            else:
                root_item.addChild(item)
                root_item.setExpanded(True)
                item.node = node
                item.index = index
                if len(index) < 2:
                    item.setExpanded(True)
                for i in range(len(node.children)):
                    self.update_tree_recursive(item, index + [i], dat_index)
               
    def get_plot_params(self):
        """
        Gets information for plotting.
        :return:
        """
        if self.indices is not None:

            ret = []

            for trace in self.indices:
                data = self.bundle.data[trace]

                pul = self.pul[trace[0]][trace[1]][trace[2]][trace[3]]
                y_label = pul.Label
                y_units = pul.YUnit
                x_interval = pul.XInterval

                ret.append((data, x_interval, y_label, y_units))

            return ret

        else:
            return [(None, None, None, None)]

    def on_selection_changed(self):
        """
        Event for browsing through wave traces
        :return:
        """
        selected = self.selectedItems()
        indices = []

        for j, item in enumerate(selected):    ### for muliple selection
            index = item.index
            # print(index)
            if len(index) == 2: ## grab series level of data
                indices.append(index)
        self.indices = indices
        self.frame.currentPulseTree = self
        self.frame.update_trace_plot()
        
class TableView(pg.TableWidget):
    """
    Class for a table
    """
    def __init__(self, parent, title = '', editable=False, sortable=True, mdcFunc=None, icFunc=None, kpFunc=None):
        super(TableView, self).__init__(editable=editable,sortable=sortable)
        self.frame = parent
        self.mouseDoubleClickFunc = mdcFunc
        self.installEventFilter(self)
        self.keyPressEventFunc = kpFunc
        self.name = title
        self.setWindowTitle(title)
        self.contextMenu.addAction(QtCore.QCoreApplication.translate("TableWidget", 'Load data')).triggered.connect(self.LoadData)
        self.contextMenu.addAction(QtCore.QCoreApplication.translate("TableWidget", 'Clear all')).triggered.connect(self.clear)
          
    def LoadData(self):
        file_name = pg.QtGui.QFileDialog.getOpenFileName(self, 'Load saved file', '',  "TSV Files (*.tsv)")
        #print(file_name)
        if file_name[0] == '':
            return
        else:
             x = pd.read_table(file_name[0])
             x = x.set_index('ChanID')
             data = x.to_records()
             self.setData(np.array(data,
                                   dtype=[('ChanID', object), ('X', object), ('Y', object)]))
             self.frame.roidata = list(data)
             self.frame.redrawAllROIs()
                    
    def mouseDoubleClickEvent(self, event):
        if self.mouseDoubleClickFunc!=None: ## if a function is attached, excute this function with this row & column's parameters
            self.mouseDoubleClickFunc(self.currentRow(), self.currentColumn())
        
    def keyPressEvent(self, event):
         key = event.key()
         if key == QtCore.Qt.Key_Return or key == QtCore.Qt.Key_Enter:
             # Process current item here            
             self.keyPressEventFunc(self.currentRow(), self.currentColumn())
         else:
             super(TableView, self).keyPressEvent(event)
             
def showdialog(msgText, details = '',withAbort=True):
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Information)
    msg.setText(msgText)
    msg.setWindowTitle("Help information")
    msg.setDetailedText(details)
    if withAbort:
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Abort)
    else:
        msg.setStandardButtons(QMessageBox.Ok)
             
    retval = msg.exec_()
    return retval

class ImageView_gv(pg.ImageView):
    """
    ## TODO
    ## NOT WORKING~~
    Wrapper around the ImageView to add more contex menus
    """

    def __init__(self, *args, SynAtlasFunc,**kwargs):
        super().__init__(*args,**kwargs)
        self.SynAtlas = SynAtlasFunc
        
    @classmethod
    def buildMenu(self):
        self.menu = super().buildMenu()
        self.SynAction = pg.QtGui.QAction(QtCore.QCoreApplication.translate('ImageView','Syn to Atlas'), self.menu)
        self.SynAction.triggered.connect(self.SynAtlas)
        self.menu.addAction(self.SynAction)           

        

def getfiles():
    dlg = QtWidgets.QFileDialog()
    dlg.setFileMode(QtWidgets.QFileDialog.AnyFile) 
    dlg.setNameFilter("Text files (*.txt *.csv)")		
    if dlg.exec_():
        filenames = dlg.selectedFiles()
        return filenames[0]
    else:
        return ''
    
      