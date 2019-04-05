#!/usr/bin/env python
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@physik.tu-freiberg.de>
# Created at: Sa 28. Jul 23:05:52 CEST 2018
# Computer: yopad 
# System: Linux 4.4.0-67-generic on x86_64
#
# Copyright (c) 2018 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------
import sys
import os
import random
import matplotlib
# Make sure that we are using QT5
matplotlib.use('Qt5Agg')

from PyQt5 import QtCore, QtGui
import PyQt5.QtWidgets as QW

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

progname = os.path.basename(sys.argv[0])

from pyxrr import structure,  Model
from pyxrr.measurement import Measurement
import numpy as np

protection = structure.Layer("Si", name="protection", thickness=100.4, roughness=4.514, density=2.33)


Multilayer = structure.Group([], periods=499, roughness=2.651)
Multilayer.append(structure.Layer("Mo", name="Moly", thickness=6.4727, roughness=2.4884, density=10.118))
Multilayer.append(structure.Layer("B4C", name="Boron Carbide", thickness=8.7981, roughness=4.461, density=2.456))

substrate = structure.Layer("Si", roughness=2., name="Silicon")

mystack = structure.Stack([protection, Multilayer], substrate=substrate)

theta = np.linspace(0, 2, 1001)
R = np.exp(-theta)
m1 = Measurement(theta, R)
m2 = Measurement(theta, R**1.5)
m3 = Measurement(theta, R**2)

ref = Model(mystack, (m1, m2, m3))
#for p in ref.params:
#    ref.params[p].vary = False
#mystack.update()



class PyXrrCanvas(FigureCanvas):
    def __init__(self, parent, model, width=5, height=5, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi, facecolor="0.95")
        self.model = model
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QW.QSizePolicy.Expanding,
                                   QW.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.redraw()
        self.setParent(parent)

    def redraw(self):
        #print("Redrawing Figure")
        self.fig.clf()
        measurements = self.model.measurements
        num = len(measurements)
        for idm in sorted(measurements):
            m = measurements[idm]
            ax = self.fig.add_subplot(num, 1, idm+1)
            ax.semilogy(m.get_theta(),
                        m.reflectivity,
                        "+-k",
                        label="measured")
            ax.semilogy(m.get_theta(),
                        self.model.reflectivity(idm=idm),
                        "-r",
                        label="fit")
            ax.semilogy(m.get_theta(),
                        self.model.reflectivity(idm=idm),
                        "-",
                        c="0.3",
                        alpha=0.6,
                        label="current")
            ax.set_ylabel("Reflectivity")
        ax.set_xlabel("glancing angle (deg)")
        self.fig.tight_layout()
        self.draw()


    def _update(self, which=None):
        _all = "measured", "fit", "current"
        if which is None:
            which = _all
        if isinstance(which, str):
            which = (which,)
        #print("Redrawing Figure")
        axes = self.fig.axes
        measurements = self.model.measurements
        for i, ax in enumerate(axes):
            idm = list(sorted(measurements))[i]
            m = measurements[idm]
            for label in which:
                il = _all.index(label)
                if label=="measured":
                    x = m.get_theta()
                    y = m.reflectivity
                elif label=="fit":
                    x = m.get_theta()
                    y = self.model.reflectivity(idm=idm)
                elif label=="current":
                    x = m.get_theta()
                    y = self.model.reflectivity(idm=idm)
                line = ax.lines[il]
                line.set_data(x,y)
        self.draw()







class PyXrrPlot(QW.QWidget):
    def __init__(self, parent, model, *args, **kwargs):
        super(PyXrrPlot, self).__init__(parent, *args, **kwargs)
        self.setLayout(QW.QVBoxLayout())
        self.canvas = PyXrrCanvas(self, model, width=6, height=7)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout().addWidget(self.toolbar)
        self.layout().addWidget(self.canvas)





class PyXrrStackDelegate(QW.QItemDelegate):
    def createEditor(self, parent, option, index):
        if index.column() != 0:
            return super(PyXrrStackDelegate, self).createEditor(parent, option, index)
        return None

class PyXrrParamDelegate(QW.QItemDelegate):
    def createEditor(self, parent, option, index):
        col = index.column()
        if col > 1:# and col != 2:
            return super(PyXrrParamDelegate, self).createEditor(parent, option, index)
        return None




class PyXrrStackItem(QW.QTreeWidgetItem):
    def __init__(self, pyxrrObject):
        self.pyxrrObject = pyxrrObject
        if isinstance(pyxrrObject, structure.Group):
            super(PyXrrStackItem, self).__init__([
                    "%i"%pyxrrObject.id,
                    pyxrrObject.name,
                    "%i"%pyxrrObject.periods
                ])
        elif isinstance(pyxrrObject, structure.Layer):
            super(PyXrrStackItem, self).__init__([
                    "%i"%pyxrrObject.id,
                    pyxrrObject.name,
                    pyxrrObject.composition
                ])
        else:
            raise ValueError("Unknown Treeview item: %s"
                             %str(pyxrrObject))
    def setData(self, column, role, value): 
        if column == 1:
            self.pyxrrObject.name = value
        elif column == 2:
            if isinstance(self.pyxrrObject, structure.Group):
                try:
                    self.pyxrrObject.periods = value
                except Exception as emsg:
                    print(emsg)
                    return
            elif isinstance(self.pyxrrObject, structure.Layer):
                try:
                    self.pyxrrObject.composition = str(value)
                except Exception as emsg:
                    print(emsg)
                    return

        return super(PyXrrStackItem, self).setData(column, role, value)
        #print(column, role, value)



class PyXrrParamItem(QW.QTreeWidgetItem):
    def __init__(self, param, name):
        self.param = param
        self.name = name
        super(PyXrrParamItem, self).__init__([
                                name,
                                param.name,
                                "%g"%param.value,
                                None,
                                "%g"%param.min,
                                "%g"%param.max
                            ])
    def setData(self, column, role, value): 
        if column == 3:
            if "inf" in self.text(2):
                self.param.vary = False
                return
            self.param.vary = bool(value)
        try:
            value = float(value)
        except Exception as emsg:
            print(emsg)
            return
        if column == 2:
            self.param.value = float(value)
        elif column == 4:
            self.param.min = float(value)
        elif column == 5:
            self.param.max = float(value)
        return super(PyXrrParamItem, self).setData(column, role, value)





class PyXrrStackTree(QW.QTreeWidget):
    def __init__(self, *args, **kwargs):
        super(PyXrrStackTree, self).__init__(*args, **kwargs)
        self.setColumnCount(3)
        self.setHeaderLabels(["#id", "Name", "Composition / Periods"])
        self.setItemDelegate(PyXrrStackDelegate(self))
        #self.resize(200, 200)

    def process_stack(self, stack):
        self.clear()
        for group in stack:
            if isinstance(group, structure.Group):
                item = PyXrrStackItem(group)#, parent=self)
                item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable)
                self.addTopLevelItem(item)
                for layer in group:
                    _item = PyXrrStackItem(layer)#, parent=self)
                    _item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable)
                    item.addChild(_item)
            elif isinstance(group, structure.Layer):
                item = PyXrrStackItem(group)#, parent=self)
                item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable)
                self.addTopLevelItem(item)




class PyXrrStackWidget(QW.QWidget):
    def __init__(self, stack, *args, **kwargs):
        super(PyXrrStackWidget, self).__init__(*args, **kwargs)
        self.stack = stack
        self.layout = QW.QHBoxLayout(self)
        self.stacktree = tree = PyXrrStackTree(self)
        tree.process_stack(stack) # TODO: make local variable
        self.layout.addWidget(tree)

        self.btns = dict()

        self.form = form = QW.QFrame(self)
        form.setFrameShape(QW.QFrame.StyledPanel)
        form.layout = QW.QFormLayout(form)
        hbox = QW.QVBoxLayout()
        for k in ("Move &Up", "Move &Down", "New Layer", "New Group", "&Remove", "&Apply\nChanges"):
            name = k.replace("&","").replace("\n", " ")
            btn = QW.QPushButton(k, self)
            btn.resize(btn.minimumSizeHint())
            btn.setStatusTip("Macht irgendwas")
            btn.clicked.connect(getattr(self, name.lower().replace(" ", "_")))
            self.btns[name] = btn
            hbox.addWidget(btn)
        form.layout.addRow(hbox)
        self.layout.addWidget(form)

    def move_up(self):
        return self._move(-1)

    def move_down(self):
        return self._move(1)

    def _move(self, step):
        getSelected = self.stacktree.selectedItems()
        if not getSelected:
            return
        item = getSelected[0]
        pyxrrObject = item.pyxrrObject
        parent = item.parent()
        stacktree = self.stacktree

        if parent is None:
            idx = self.stack.index(pyxrrObject)
            if step == "remove":
                self.stack.pop(idx)
                stacktree.takeTopLevelItem(idx)
            elif step == "newgroup":
                newgroup = structure.Group((), 1, 0.)
                groupitem = PyXrrStackItem(newgroup)#, parent=stacktree)
                groupitem.setFlags(groupitem.flags() | QtCore.Qt.ItemIsEditable)
                self.stack.insert(idx, newgroup)
                stacktree.insertTopLevelItem(idx, groupitem)

                newlayer = structure.Layer("Al", density=0., roughness=0., thickness=0., name="New")
                layeritem = PyXrrStackItem(newlayer)#, parent=groupitem)
                layeritem.setFlags(layeritem.flags() | QtCore.Qt.ItemIsEditable)
                newgroup.append(newlayer)
                groupitem.addChild(layeritem)
            elif step == "newlayer":
                newlayer = structure.Layer("Al", density=0., roughness=0., thickness=0., name="New")
                newitem = PyXrrStackItem(newlayer)#, parent=stacktree)
                newitem.setFlags(newitem.flags() | QtCore.Qt.ItemIsEditable)
                self.stack.insert(idx, newlayer)
                stacktree.insertTopLevelItem(idx, newitem)
            else:
                new_idx = idx+step
                if new_idx<0 or new_idx>=(len(self.stack)):
                    return
                self.stack.insert(new_idx, self.stack.pop(idx))
                stacktree.insertTopLevelItem(new_idx, stacktree.takeTopLevelItem(idx))
        else:
            group = parent.pyxrrObject
            idx = group.index(pyxrrObject)
            if step == "remove":
                group.pop(idx)
                parent.takeChild(idx)
                return
            elif step == "newlayer":
                newlayer = structure.Layer("Al", density=0., roughness=0., thickness=0., name="New")
                newitem = PyXrrStackItem(newlayer, parent=parent)
                newitem.setFlags(newitem.flags() | QtCore.Qt.ItemIsEditable)
                self.stack.insert(idx, newlayer)
                parent.insertChild(idx, newitem)
            elif step == "newgroup":
                return
            else:
                new_idx = idx+step
                if new_idx<0 or new_idx>=(len(group)):
                    return
                group.insert(new_idx, group.pop(idx))
                parent.insertChild(new_idx, parent.takeChild(idx))

    def remove(self):
        self._move("remove")

    def new_layer(self):
        self._move("newlayer")

    def new_group(self):
        self._move("newgroup")

    def apply_changes(self):
        self.stack.update()
        self.stacktree.process_stack(self.stack)
        mainwidget = self.parent().parent().parent()
        mainwidget.paramwidget.paramtree.process_model()
        mainwidget.plot.canvas.redraw()




### PARAMETERS LIST ###

class PyXrrParamTree(QW.QTreeWidget):
    def __init__(self, *args, **kwargs):
        super(PyXrrParamTree, self).__init__(*args, **kwargs)
        self.setColumnCount(6)
        self.setItemDelegate(PyXrrParamDelegate(self))
        self.setHeaderLabels(["Object", "Quantity", "Value", "Vary?", "Min", "Max"])

    def process_model(self, model=None):
        if model is None:
            model = self.parent().model
        self.clear()
        model.update()
        last = None
        for pname in (model.params):
            param = model.params[pname]
            name = param.parent.name
            item = PyXrrParamItem(param, name)
            item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(3, QtCore.Qt.Unchecked)
            self.addTopLevelItem(item)
            #if param.parent is not last:
            #    item = QW.QTreeWidgetItem()
            #    item.setFlags(QtCore.Qt.NoItemFlags)
            #    #item.setSizeHint(0, QtCore.QSize(1,1))
            #    self.addTopLevelItem(item)
            #last = param.parent



class PyXrrParamWidget(QW.QWidget):
    def __init__(self, model, *args, **kwargs):
        super(PyXrrParamWidget, self).__init__(*args, **kwargs)
        self.model = model
        self.layout = QW.QVBoxLayout(self)
        self.paramtree = tree = PyXrrParamTree(self)
        self.layout.addWidget(QW.QLabel("Parameters"))
        self.layout.addWidget(tree)
        tree.process_model(model)
        for column in range(tree.columnCount()):
            tree.resizeColumnToContents(column)

        self.btns = dict()

        self.form = form = QW.QFrame(self)
        form.setFrameShape(QW.QFrame.StyledPanel)
        form.layout = QW.QFormLayout(form)
        hbox = QW.QHBoxLayout()
        for k in ("<<", "<", ">", ">>"):
            name = k.replace("&","").replace("\n", " ")
            btn = QW.QPushButton(k, self)
            btn.resize(btn.minimumSizeHint())
            btn.setStatusTip("Macht irgendwas")
            self.btns[name] = btn
            hbox.addWidget(btn)
        self.btns["<<"].clicked.connect(lambda: self.tune(0.5))
        self.btns["<"].clicked.connect(lambda: self.tune(0.5**0.1))
        self.btns[">"].clicked.connect(lambda: self.tune(2**0.1))
        self.btns[">>"].clicked.connect(lambda: self.tune(2))
        form.layout.addRow(hbox)
        self.layout.addWidget(form)
    
    def tune(self, direction):
        getSelected = self.paramtree.selectedItems()
        if not getSelected:
            return
        item = getSelected[0]
        param = item.param
        if param.value != 0:
            newval = param.value*direction
        else:
            newval = 0.01 * ((direction>1)-0.5)
        #param.value *= self.factors[direction]
        item.setText(2, "%g"%newval)
        self.parent().parent().parent().plot.canvas._update("current")



class PyXrrMeasurementItem(QW.QTreeWidgetItem):
    def __init__(self, measurement):
        self.pyxrrObject = pyxrrObject
        if not isinstance(measurement, Measurement):
            raise ValueError("Invalid Treeview item: %s"
                             %str(measurement))
        key = "%i"%measurement.id
        if measurement.name:
            key += " (%s)"%measurement.name
        super(PyXrrMeasurementItem, self).__init__([
                key,
                "%g...%g"%tuple(measurement.x[[0,-1]]),
                measurement.xaxis,
                measurement.is_regular,
                str(measurement.fit_range),
                "%g"%measurement.energy,
                "%g"%measurement.scale,
                "%g"%measurement.offset,
                "%g"%measurement.resolution,
                "%g"%measurement.polarization,
                "%g"%measurement.background
                ])
        self.setFlags(self.flags() | QtCore.Qt.ItemIsUserCheckable)


class PyXrrMeasurementTree(QW.QTreeWidget):
    def __init__(self, *args, **kwargs):
        super(PyXrrParamTree, self).__init__(measurements,*args, **kwargs)
        self.measurements = measurements
        self.setColumnCount(9)
        self.setHeaderLabels(["#id",
                              "range",
                              "x-axis",
                              "regular",
                              "fit range",
                              "energy",
                              "scale",
                              "offset",
                              "resolution",
                              "polarization",
                              "background"])
        self.process()


    def process(self, measurements=None):
        if measurements is None:
            measurements = self.measurements
        self.clear()
        for idm in measurements:
            item = PyXrrMeasurementItem(param, name)
            self.addTopLevelItem(item)


class MeasurementsDialog(QW.QDialog): # Base Class
    Input = dict()
    def __init__(self, parent=None, measurements=dict()):
        super(MeasurementsDialog, self).__init__(parent)
        self.resize(300,200)
        self.measurements = measurements
        self.home()

    def home(self):
        font = QtGui.QFont()
        font.setPointSize(9)
        self.setFont(font)

        layout = QW.QHBoxLayout(self)
        tree = PyXrrMeasurementTree(self.measurements)
        layout.addWidget()
        self.setLayout(layout)


    def showEvent(self, event):
        geom = self.frameGeometry()
        geom.moveCenter(QtGui.QCursor.pos())
        self.setGeometry(geom)
        super(MeasurementsDialog, self).showEvent(event)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Enter:
            pass
        elif event.key() == QtCore.Qt.Key_Escape:
            self.hide()
            event.accept()
        else:
            super(MeasurementsDialog, self).keyPressEvent(event)


### MAIN WINDOW ###
class ApplicationWindow(QW.QMainWindow):
    def __init__(self, model):
        QW.QMainWindow.__init__(self)
        self.model = model
        self.stack = model.stack
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("PyXRR -- GUI")

        self.file_menu = QW.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.edit_menu = QW.QMenu('&Edit', self)
        self.edit_menu.addAction('&Measurements', self.measurements,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_M)

        self.help_menu = QW.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.edit_menu)
        self.menuBar().addMenu(self.help_menu)

        self.help_menu.addAction('&About', self.about)

        #self.main_widget = QW.QWidget(self)
        #self.main_widget.setFocus()
        #self.setCentralWidget(self.main_widget)

        #self.statusBar().showMessage("All hail matplotlib!", 2000)
        self.setGeometry(50, 50, 1200, 800)

        self.home()


    def home(self):
        self.splith = splitter = QW.QSplitter(QtCore.Qt.Horizontal)
        self.setCentralWidget(splitter)

        self.splitleft = QW.QSplitter(QtCore.Qt.Vertical, parent=self) 
        self.stackwidget = PyXrrStackWidget(self.stack, parent=self.splitleft)
        self.paramwidget = PyXrrParamWidget(self.model, parent=self.splitleft)
        self.splitleft.addWidget(self.stackwidget)
        self.splitleft.addWidget(self.paramwidget)
        self.splitleft.setStretchFactor(1, 1)
        #self.splitleft.setSizes((800, 1000))

        self.plot = PyXrrPlot(splitter, self.model)

        splitter.addWidget(self.splitleft)
        splitter.addWidget(self.plot)

    def fileQuit(self):
        self.close()

    def measurements(self):
        dialog = MeasurementsDialog(self)
        dialog.exec_()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        QMessageBox.about(self, "About", "PyXRR")


qApp = QW.QApplication(sys.argv)

aw = ApplicationWindow(ref)
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()
