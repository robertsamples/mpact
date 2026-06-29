"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
"""
import concurrent.futures
import concurrent
import numpy as np
import pandas as pd
import pickle

from csvcache import cached_read_csv, invalidate as invalidate_csv_cache
import ordination
import clusterpurity

import matplotlib
#matplotlib.style.use('ggplot')
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT)
from matplotlib.figure import Figure
import matplotlib.ticker as ticker
import seaborn as sns
import squarify # pip install squarify
import upsetplot
import matplotlib.colors as mc
import colorsys
#plt.rc('xtick',labelsize=16)
#plt.rc('ytick',labelsize=16)

import platform
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect, QSize, QTime, QUrl, Qt, QEvent)
from PyQt5.QtGui import (QBrush, QColor, QIcon, QPalette, QPainter, QPixmap)
from PyQt5.QtWidgets import *
from pathlib import Path

import scipy.cluster.hierarchy as shc
from sklearn.preprocessing import normalize
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Ellipse
from filter import listfilter
import time
import math

from pvclust import PvClust

class NavigationToolbar(NavigationToolbar2QT):
        """
        This is a subclass of the NavigationToolbar2QT that provides additional functionality. It is responsible for creating a custom button in the plot toolbar. The custom button is used to configure the plot.

        Constructor arguments:
        - canvas: a FigureCanvasQTAgg instance representing the plot canvas
        - parent: the parent widget
        
        Attributes:
        - clearButtons: a list of QToolButton instances representing the plot toolbar buttons
        - picker: a QAction instance representing the custom button
        
        Methods:
        - __init__(self, canvas, parent): constructs a NavigationToolbar instance
        """
        def __init__(self, canvas, parent):
            NavigationToolbar2QT.__init__(self,canvas,parent)
            self.clearButtons=[]
            next=None
            for c in self.findChildren(QToolButton):
                if next is None:
                    next=c
                if str(c.text()) in ('Pan','Zoom','Subplots'):
                    self.clearButtons.append(c)
                    next=None 
            # create custom button
            icon = QtGui.QIcon("cog.ico")
            picker=QAction("Pick",self)
            picker.setIcon(icon)
            picker.setToolTip("Configure")
            self.picker = picker
            button=QToolButton(self)
            button.setDefaultAction(self.picker) 
            button.clicked.connect(lambda: parent.dialog.show())
            self.addWidget(button)
        
def _is_duplicate_pick(parent, event):
    """True if this pick event's underlying mouse click already triggered a
    pick on another overlapping artist.

    Matplotlib fires one ``pick_event`` per artist whose hit-test matches a
    click, not one per click. When the same feature is plotted in more than
    one groupset/colour layer (separate ``ax.scatter()`` calls at the same
    coordinates), clicking it fires our callback once per overlapping layer
    for a single physical click. Each call reselects the same feature, which
    -- because selecting an already-selected feature toggles the highlight
    off -- made the net result of clicking an overplotted feature an
    immediate deselect. Deduplicating by the shared ``mouseevent`` (every pick
    event from one click carries the same ``mouseevent`` instance) makes one
    click produce exactly one selection regardless of how many layers it hits.
    """
    mouseevent = getattr(event, 'mouseevent', None)
    if mouseevent is None:
        return False
    if getattr(parent, '_last_pick_mouseevent', None) is mouseevent:
        return True
    parent._last_pick_mouseevent = mouseevent
    return False


#General plot method#
class ui_plot:
    """
    This class is used to define the plot layout and properties for the plot canvas. It is responsible for initializing the plot canvas and the plot axes, setting plot parameters, and resetting the plot.
                                                                        
    Constructor arguments:
    - parent: the parent widget
    - currplt: an integer representing the current plot number
    - frame: a widget used as a container for the plot
    
    Attributes:
    - fcsfont: a dictionary containing the font properties for text labels
    - fhfont: a dictionary containing the font properties for title labels
    - plotbackground: a tuple containing the plot background color
    - event: a string representing the type of event that triggered the plot update
    
    Methods:
    - __init__(self, parent, currplt, frame): constructs a ui_plot instance
    - onpick(self, event, parent, iondict, plotcols): a method that highlights a picked feature on the plot and updates the label information
    - reset(self, file, filtereddfs, groupsets): a method that clears the plot axes and updates the plot with new data.
    """
    
    def __init__(self, parent, currplt, frame):
        parent.fig[currplt] = Figure()
        parent.pltlayout[currplt] = QtWidgets.QVBoxLayout()
        parent.canvas[currplt] = FigureCanvas(parent.fig[currplt])
        parent.pltlayout[currplt].addWidget(parent.canvas[currplt])
        parent.toolbar[currplt] = NavigationToolbar(parent.canvas[currplt], parent)
        parent.toolbar[currplt].setStyleSheet("background-color:rgba(225,225,225,0);")
        parent.pltlayout[currplt].addWidget(parent.toolbar[currplt])
        frame.setLayout(parent.pltlayout[currplt])
        parent.ax[currplt] = parent.canvas[currplt].figure.subplots()
        parent.ax[currplt].set_axisbelow(True)
        self.event = ''

        self.fcsfont = {'fontname': 'Bahnschrift', 'weight': 'bold', 'size': 20}
        self.fhfont = {'fontname': 'Bahnschrift', 'weight': 'bold', 'size': 25}
        self.plotbackground = (.89, .89, .89, 0)
        parent.canvas[currplt].figure.set_facecolor(self.plotbackground)
        parent.ax[currplt].set_facecolor(self.plotbackground)

    def onpick(self, event, parent, iondict, plotcols):
        if _is_duplicate_pick(parent, event):
            return
        ind = event.ind
        coord = event.artist.get_offsets()[ind, :]
        pickedfeature = iondict.loc[iondict[plotcols[0]] == coord[0, 0], :].loc[iondict[plotcols[1]] == coord[0, 1], :]
        if pickedfeature.empty:
            return
        parent.ui.lbl_featurename.setText('Compound: ' + str(pickedfeature.iloc[0, 0]))
        parent.ui.lbl_featurert.setText('Retention time: ' + str(round(pickedfeature.iloc[0, 2], 4)))
        parent.ui.lbl_featuremz.setText('m/z: ' + str(round(pickedfeature.iloc[0, 1], 4)))
        pickedfeature = str(pickedfeature.iloc[0, 0])
        parent.highlight_feature(pickedfeature)

    def reset(self, file, filtereddfs, groupsets):
        self.parent.ax[self.currplt].clear()
        #self.parent.canvas[self.currplt].draw()
        self.parent.canvas[self.currplt].mpl_disconnect(self.event)
        self.plot(self.parent, file, filtereddfs, groupsets)


    
#subclassed plot methods
class plot_abund():
    """
    A class for generating an abundance plot made of a barchart and strip/point plots.
    
    Attributes:
    parent (QMainWindow): The parent QMainWindow object that holds the plot widget.
    currplt (int): The index of the current plot widget.
    fcsfont (dict): The font dictionary used for the figure caption.
    fhfont (dict): The font dictionary used for the figure title.
    plotbackground (tuple): The background color of the plot.
    
    Methods:
    __init__(self, parent, currplt):
        Initializes the plot widget and sets the attributes.
        
    plot(self, parent):
        Generates the abundance plot for the given feature using Seaborn bar, strip, and point plots.
        
    reset(self):
        Clears the plot and generates a new abundance plot.
    """
    def __init__(self, parent, currplt):
        self.parent = parent
        self.currplt = currplt
        parent.fig[currplt] = Figure()
        parent.pltlayout[currplt] = QtWidgets.QVBoxLayout()
        parent.canvas[currplt] = FigureCanvas(parent.fig[currplt])
        parent.pltlayout[currplt].addWidget(parent.canvas[currplt])
        parent.toolbar[currplt] = NavigationToolbar(parent.canvas[currplt], parent)
        parent.toolbar[currplt].setStyleSheet("background-color:rgba(225,225,225,0);")
        parent.pltlayout[currplt].addWidget(parent.toolbar[currplt])
        parent.ftrdialog.ui.frame_abundance.setLayout(parent.pltlayout[currplt])
        parent.ax[currplt] = parent.canvas[currplt].figure.subplots(ncols=2)  

        self.plotbackground = (.89, .89, .89, 0)
        parent.canvas[currplt].figure.set_facecolor(self.plotbackground)
        self.fcsfont = {'fontname':'Bahnschrift', 'weight': 'bold', 'size': 12}
        self.fhfont = {'fontname':'Bahnschrift', 'weight': 'bold', 'size': 25}

    
    def plot(self, parent):
        # Get header info
        currplt = self.currplt
        msdata = cached_read_csv(parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_filtered.csv'), sep=',', header=[0, 1, 2], index_col=[0]).iloc[:, 2:].loc[parent.pickedfeature]
        msdata = msdata.reset_index()
        msdata.columns = ['biolgroup','sample','injection','abundance']
        msdata = msdata.drop(columns=['injection'])

        # Get stats for the given ion
        summary = cached_read_csv(parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_summarydata.csv'), sep=',', header=[0, 1], index_col=[0]).iloc[:, 2:].loc[parent.pickedfeature]
        combasd = summary.loc[['combASD']].to_frame()
        combasd.index = combasd.index.droplevel(level=0)
        neff = summary.loc[['neff']].to_frame()
        neff.index = neff.index.droplevel(level=0)
        ionsummary = summary.loc[['average']]
        ionsummary.index = ionsummary.index.droplevel(level=0)
    
        ionsummary = ionsummary.to_frame()
        ionsummary.columns = ['average']
        ionsummary['combASD'] = combasd
        ionsummary['neff'] = neff
        if parent.analysis_paramsgui.blnkfltr:
            ionsummary = ionsummary.drop([parent.analysis_paramsgui.blnkgrp], axis=0)
        # Force the index to be unnamed before resetting it, so the new column
        # is always literally named "index" (used below by sns.barplot) --
        # newer pandas versions can carry a level name through droplevel()
        # where older ones silently dropped it, which otherwise renames this
        # column out from under the hardcoded x="index" reference.
        ionsummary.index.name = None
        ionsummary = ionsummary.reset_index()

        # Use one explicit, shared palette keyed by biolgroup so the bar
        # chart's per-group colors (left) match the stripplot's hue colors
        # (right) instead of each picking its own default palette.
        groups = sorted(msdata['biolgroup'].unique())
        palette = dict(zip(groups, sns.color_palette(n_colors=len(groups))))

        # Define the plotting function for stripplot
        def plot_stripplot():
            sns.stripplot(ax=parent.ax[currplt][1], x="sample", y="abundance", hue="biolgroup",
                          data=msdata, palette=palette, size=7, dodge=False, alpha=.5, zorder=2)
    
        # Define the plotting function for pointplot
        def plot_pointplot():
            sns.pointplot(ax=parent.ax[currplt][1], x="sample", y="abundance", color="#999999",
                          data=msdata, dodge=False, linestyle='none', markers="d", markersize=6,
                          zorder=1, capsize=.1, err_kws={'linewidth': 1}, errorbar=('ci', 95))
    
        # Define a list of plotting functions to be run in parallel
        plot_functions = [plot_stripplot, plot_pointplot]
    
        # Make plots for the features using multithreading
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = [executor.submit(func) for func in plot_functions]
    
        # Wait for all threads to finish
        for result in results:
            result.result()
    
        # ionsummary is already one pre-aggregated row per category, not raw
        # per-sample data -- use matplotlib's ax.bar (not sns.barplot) since
        # newer seaborn no longer supports passing a manual yerr array through
        # to barplot (it routes kwargs through an internal single-point
        # "scout" call that only accepts a scalar/1-matching yerr).
        barcolors = [palette[g] for g in ionsummary["index"]]
        parent.ax[currplt][0].bar(ionsummary["index"], ionsummary["average"], yerr=ionsummary["combASD"],
                                  color=barcolors, ecolor=".2", edgecolor=".2", zorder=1, capsize=3)
        parent.ax[currplt][1].set_xticklabels(parent.ax[currplt][1].get_xticklabels(), rotation=90, horizontalalignment='center')
    
        ylims = (0, parent.ax[currplt][1].get_ylim()[1] * 1.05)
        parent.ax[currplt][0].set_ylim(ylims)
        parent.ax[currplt][1].set_ylim(ylims)
        parent.ax[currplt][1].legend_.remove()
        plt.tight_layout
    
    
    def reset(self):
        self.parent.ax[self.currplt][0].clear()
        self.parent.ax[self.currplt][1].clear()
        self.plot(self.parent)
        self.parent.canvas[self.currplt].draw()

class show_spectrum(ui_plot):
    """
    A class for displaying a MS2 spectrum viewer plot widget.

    Attributes:
    parent (QMainWindow): The parent QMainWindow object that holds the plot widget.
    currplt (int): The index of the current plot widget.
    fcsfont (dict): The font dictionary used for the figure caption.
    fhfont (dict): The font dictionary used for the figure title.
    
    Methods:
    __init__(self, parent, currplt):
        Initializes the plot widget and sets the attributes.
        
    plot(self, frags):
        Generates the MS2 spectrum plot for the given fragments.
        
    reset(self, frags):
        Clears the plot and generates a new MS2 spectrum plot.
    """
    def __init__(self, parent, currplt):
        super().__init__(parent, currplt, parent.ftrdialog.ui.frame_spec)
        self.parent = parent
        self.currplt = currplt

    def plot(self, frags):
        ax = self.parent.ax[self.currplt]
        ax.vlines(frags[:, 0], 0, frags[:, 1], colors='k', linestyles='solid')
        ax.set_xlabel('m/z', **self.fcsfont)
        ax.set_ylabel('Abundance', **self.fcsfont)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        self.parent.canvas[self.currplt].draw()

    def reset(self, frags):
        self.parent.ax[self.currplt].clear()
        self.plot(frags)
        
class show_featureplt(ui_plot):
    """
    A class for displaying a filtered feature plot widget.

    Attributes:
    parent (QMainWindow): The parent QMainWindow object that holds the plot widget.
    currplt (int): The index of the current plot widget.
    fcsfont (dict): The font dictionary used for the figure caption.
    fhfont (dict): The font dictionary used for the figure title.
    
    Methods:
    __init__(self, parent, currplt, frame, iondictloc, filtereddfs, groupsets):
        Initializes the plot widget and sets the attributes.
        
    plot(self, parent, iondictloc, filtereddfs, groupsets):
        generates plot of m/z value and retention for each feature.
    """
    def __init__(self, parent, currplt, frame, iondictloc, filtereddfs, groupsets):
        ui_plot.__init__(self, parent, currplt, frame)
        self.parent = parent
        self.currplt = currplt
        self.plot(parent, iondictloc, filtereddfs, groupsets)
        
    def plot(self, parent, iondictloc, filtereddfs, groupsets):
        iondict = pd.read_csv(iondictloc, sep=',', header=[0], index_col=None)
        iondict['colour'] = '#000000' # colour map for features based on which filter they are caught by
        
        if parent.analysis_paramsgui.decon:
            iondict.loc[iondict['pass_insource'] == False, 'colour'] = '#00aa00'
        if parent.analysis_paramsgui.CVfil:
            iondict.loc[iondict[parent.analysis_paramsgui.cvparam] >= parent.analysis_paramsgui.cvthresh, 'colour'] = '#ff0000' 
        if parent.analysis_paramsgui.blnkfltr:
            iondict.loc[iondict['pass_blnkfil'] == False, 'colour'] = '#00aaaa'
        if parent.analysis_paramsgui.relfil:
            iondict.loc[iondict['pass_relfil'] == False, 'colour'] = '#0000ff'
    
        colours = ['#0000ff', '#00aaaa', '#ff0000', '#00aa00', '#000000']
        labels = ["Mispicked", "Blank features", "Nonreproducable", "In-source ions", "High Quality"]
        for i in range(len(colours)):
            parent.ax[self.currplt].scatter(iondict.loc[iondict['colour'] == colours[i], :]['Retention time (min)'], iondict.loc[iondict['colour'] == colours[i], :]['m/z'], color=colours[i], label=labels[i], picker=True, alpha=.5)
        
        parent.highlight[self.currplt], = parent.ax[self.currplt].plot([], [], 'o', markersize=12, color='yellow')
        parent.ax[self.currplt].set_xlabel('Retention time (min)', **self.fcsfont)
        parent.ax[self.currplt].set_ylabel('m/z', **self.fcsfont)
        parent.ax[self.currplt].set_xlim([0, 11])
        parent.ax[self.currplt].set_ylim([0, 2000])
        parent.ax[self.currplt].grid()
        parent.ax[self.currplt].legend()
        self.event = parent.canvas[self.currplt].figure.canvas.mpl_connect('pick_event', lambda event: self.onpick(event, parent, iondict, ('Retention time (min)', 'm/z')))
        parent.fig[self.currplt].subplots_adjust(left=.1, right=.95, bottom=0.15, top=0.9, hspace=0.2, wspace=0.2)
        parent.canvas[self.currplt].draw()


class plot_heatmap():
    """
    This class generates a heatmap and adds it to the GUI. The heatmap is created using seaborn's clustermap function
    and is displayed using matplotlib. The heatmap is interactive, allowing the user to select individual features and
    see additional information about them.

    Attributes:
        plotbackground (tuple): The background color of the heatmap.
        fcsfont (dict): A dictionary containing font information for labels on the heatmap.
        fhfont (dict): A dictionary containing font information for the heatmap title.
        event: The event that allows the user to select features on the heatmap.

    Args:
        parent (object): The parent object that the heatmap will be added to.
        currplt (str): The name of the plot.
        frame (object): The object that the heatmap will be added to.
        file (str): The file that contains the data for the heatmap.
    """
    def __init__(self, parent, currplt, frame, file):
        msdata = cached_read_csv(parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_filtered.csv'), sep = ',', header = [2], index_col = [0]).iloc[:,2:]
        cm = sns.clustermap(msdata, standard_scale=0, metric="euclidean", method="ward", cmap = parent.analysis_paramsgui.colorscheme) #viridis
        parent.cmind = cm.dendrogram_row.reordered_ind #saves reordered index so that we can increment selection up and down.
        parent.fig[currplt] = cm.fig
        parent.pltlayout[currplt] = QtWidgets.QVBoxLayout()
        parent.canvas[currplt] = FigureCanvas(parent.fig[currplt])
        parent.pltlayout[currplt].addWidget(parent.canvas[currplt])
        parent.toolbar[currplt] = NavigationToolbar(parent.canvas[currplt], parent)
        parent.toolbar[currplt].setStyleSheet("background-color:rgba(225,225,225,255);")
        parent.pltlayout[currplt].addWidget(parent.toolbar[currplt])
        parent.ui.frame_heatmap.setLayout(parent.pltlayout[currplt])
        
        self.plotbackground = (.89, .89, .89, 0)
        parent.canvas[currplt].figure.set_facecolor(self.plotbackground)
        self.fcsfont = {'fontname':'Bahnschrift',
                'weight' : 'bold',
                'size'   : 12}
        self.fhfont = {'fontname':'Bahnschrift',
                'weight' : 'bold',
                'size'   : 25}
        parent.highlight[currplt], = parent.canvas[currplt].figure.axes[2].plot([], [], color='yellow', linestyle='-', linewidth=1)
    
        def onpick8(event): #rename
                # this gets the position of the selected item, finds the item in the reodered index,
                # and highlights the feature by name on other plots
            if _is_duplicate_pick(parent, event):
                return
            iondict = cached_read_csv(parent.analysis_paramsgui.outputdir / 'iondict.csv', sep = ',', header = [0], index_col = [0])
            coord = [event.mouseevent.xdata,event.mouseevent.ydata]
            parent.heatind = int(np.floor(coord[1]))
            name = msdata.index.tolist()[parent.cmind[parent.heatind]]    
    
            parent.ui.lbl_featurename.setText('Compound: ' + name)
            parent.ui.lbl_featurert.setText('Retention time: ' + str(iondict.loc[name,'Retention time (min)']))
            parent.ui.lbl_featuremz.setText('m/z: ' + str(iondict.loc[name,'m/z']))
            parent.highlight_feature(name)
            
        def on_xlims_change(event_ax): #changes the limits of the x dendrogram when the heatmap is zoomed
            heatxlim = event_ax.get_xlim()
            yaxxlim = (10 * heatxlim[0], 10 * heatxlim[1])
            parent.canvas['heatmap'].figure.axes[1].set_xlim(yaxxlim)
        
        def on_ylims_change(event_ax): #changes the limits of the y dendrogram when the heatmap is zoomed
            heatylim = event_ax.get_ylim()
            xaxylim = (10 * heatylim[0], 10 * heatylim[1])
            parent.canvas['heatmap'].figure.axes[0].set_ylim(xaxylim)
        
        parent.canvas[currplt].figure.axes[2].callbacks.connect('xlim_changed', on_xlims_change)
        parent.canvas[currplt].figure.axes[2].callbacks.connect('ylim_changed', on_ylims_change)
        parent.canvas[currplt].figure.axes[2].set_picker(True)
        self.event = parent.canvas[currplt].figure.canvas.mpl_connect('pick_event', onpick8)
        parent.fig[currplt].subplots_adjust(left=0.015,right=0.986,
                    bottom=0.271,top=0.975,
                    hspace=0.005,wspace=0.005)
        #parent.canvas[currplt].figure.axes[2].set_xlabel('Sample',  **self.fcsfont)
        parent.canvas[currplt].figure.axes[2].set_ylabel('Feature',  **self.fcsfont)
        parent.canvas[currplt].figure.axes[3].set_aspect(2.5)
        parent.canvas[currplt].figure.axes[2].get_yaxis().set_ticks([])
        #colind = cm.dendrogram_col.reordered_ind
        #parent.canvas[currplt].figure.axes[2].set_xticks(range(len(colind)))
        #parent.canvas[currplt].figure.axes[2].set_xticklabels(colind, rotation=0)
        parent.canvas[currplt].draw()
        
    def reset(self, parent, currplt, frame, file):
        #makes new figure with updated heatmap and saves
        msdata = cached_read_csv(parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_filtered.csv'), sep = ',', header = [2], index_col = [0]).iloc[:,2:]
        cm2 = sns.clustermap(msdata, standard_scale=0, metric="euclidean", method="ward", cmap = parent.analysis_paramsgui.colorscheme) #viridis
        parent.cmind = cm2.dendrogram_row.reordered_ind
        updatedfig = cm2.fig
        parent.canvas[currplt].figure.clf()
            
        figaxes = {} #removes original heatmap components sequentially
        for x in reversed(range(0,4)):
            figaxes[x] = updatedfig.axes[x]
            updatedfig.axes[x].remove()
        for x in range(0,4): # readds new heatmap components sequentially
            figaxes[x].figure=parent.canvas[currplt].figure
            parent.canvas[currplt].figure.axes.append(figaxes[x])
            parent.canvas[currplt].figure.add_subplot(figaxes[x])
        parent.fig[currplt].subplots_adjust(left=0.014,right=0.967,
                    bottom=0.328,top=0.806,
                    hspace=0.005,wspace=0.005)
        parent.highlight[currplt], = parent.canvas[currplt].figure.axes[2].plot([], [], color='yellow', linestyle='-', linewidth=1)
        

 
        def on_xlims_change(event_ax):#make it soi can use these methods in the first call befre regeneration
            heatxlim = event_ax.get_xlim()
            yaxxlim = (10 * heatxlim[0], 10 * heatxlim[1])
            parent.canvas['heatmap'].figure.axes[1].set_xlim(yaxxlim)
    
        def on_ylims_change(event_ax):
            heatylim = event_ax.get_ylim()
            xaxylim = (10 * heatylim[0], 10 * heatylim[1])
            parent.canvas['heatmap'].figure.axes[0].set_ylim(xaxylim)
        """
        def onresize(event): # i have no idea why this is getting automatically sized back down
            newsize = parent.canvas[currplt].figure.get_size_inches()
            resizeratio = ([newsize[0]/self.figsize[0], newsize[1]/self.figsize[1]])
            self.figsize = parent.canvas[currplt].figure.get_size_inches()
            parent.fig[currplt].subplots_adjust(left=0.014,right=0.967*resizeratio[0],
                                                bottom=0.328,top=0.806*resizeratio[1],
                                                hspace=0.005,wspace=0.005)
            parent.canvas[currplt].draw()
        """
        #parent.canvas[currplt].figure.canvas.mpl_connect('resize_event', onresize)
        parent.canvas[currplt].figure.axes[2].callbacks.connect('xlim_changed', on_xlims_change)
        parent.canvas[currplt].figure.axes[2].callbacks.connect('xlim_changed', on_xlims_change)
        parent.canvas[currplt].figure.axes[2].callbacks.connect('ylim_changed', on_ylims_change)
        parent.canvas[currplt].figure.axes[2].set_picker(True)
        
        parent.canvas[currplt].figure.set_facecolor(self.plotbackground)
        #parent.canvas[currplt].figure.axes[2].set_xlabel('Sample',  **self.fcsfont)
        parent.canvas[currplt].figure.axes[2].set_ylabel('Feature',  **self.fcsfont)
        parent.canvas[currplt].figure.axes[3].set_aspect(2.5)
        parent.canvas[currplt].figure.axes[2].get_yaxis().set_ticks([])
        self.figsize = parent.canvas[currplt].figure.get_size_inches()
        parent.canvas[currplt].draw()

class plot_mzrt(ui_plot):
    """A subclass of ui_plot that plots a feature mass vs retention time plot.
    
    Args:
        parent (object): The parent object of the plot.
        currplt (int): The index of the current plot.
        frame (object): The frame of the plot.
        file (str): The path to the file that contains the ion dictionary.
        filtereddfs (dict): A dictionary containing filtered dataframes for each group.
        groupsets (dict): A dictionary containing group settings for each group.
    
    Attributes:
        parent (object): The parent object of the plot.
        currplt (int): The index of the current plot.
    
    Methods:
        __init__(self, parent, currplt, frame, file, filtereddfs, groupsets): Initializes the plot and calls the plot method.
        plot(self, parent, file, filtereddfs, groupsets): Plots the feature mass vs retention time plot.
    """
    def __init__(self, parent, currplt, frame, file, filtereddfs, groupsets):
        ui_plot.__init__(self, parent, currplt, frame)
        self.parent = parent
        self.currplt = currplt
        self.plot(parent, file, filtereddfs, groupsets)
        
    def plot(self, parent, file, filtereddfs, groupsets): # abundance tied opacity used here currently
        iondict = cached_read_csv(parent.analysis_paramsgui.outputdir / 'iondict.csv', sep=',', header=[0], index_col=None)
        
        for elem in filtereddfs:
            if parent.analysis_paramsgui.blnkfltr:
                filtereddfs[elem] = filtereddfs[elem][filtereddfs[elem]['pass_blnkfil']]
            
            plotcol = groupsets[elem].plotcol.lstrip('#')
            rgbcol = tuple(float(int(plotcol[i:i+2], 16)/255) for i in (0, 2, 4))
            rgbacol = np.asarray([(rgbcol[0], rgbcol[1], rgbcol[2], a) for a in filtereddfs[elem]['logmax']/parent.analysis_paramsgui.maxval])
            rgbacol = np.clip(rgbacol, 0, 1)
            
            parent.ax[self.currplt].scatter(filtereddfs[elem]['Retention time (min)'].to_list(), 
                                             filtereddfs[elem]['m/z'].to_list(), 
                                             c=rgbacol, 
                                             label=str(groupsets[elem].legendname), 
                                             picker=True) 
            #parent.ax[self.currplt].scatter(filtereddfs[elem]['Retention time (min)'].to_list(), filtereddfs[elem]['m/z'].to_list(), color = str(groupsets[elem].plotcol), label = str(groupsets[elem].legendname), picker=True, alpha=.5)

        parent.highlight[self.currplt], = parent.ax[self.currplt].plot([], [], 'o', markersize=12, color='yellow')
        parent.ax[self.currplt].set_xlabel("Retention time (min)", **self.fcsfont)
        parent.ax[self.currplt].set_ylabel('m/z', **self.fcsfont)
        parent.ax[self.currplt].set_xlim(-.5, 11.5)
        parent.ax[self.currplt].set_ylim(0, 1850)
        parent.ax[self.currplt].legend()
        
        self.event = parent.canvas[self.currplt].figure.canvas.mpl_connect('pick_event', lambda event: self.onpick(event, parent, iondict, ('Retention time (min)', 'm/z')))
        parent.fig[self.currplt].subplots_adjust(left=.1,right=.95,
                            bottom=0.15,top=0.9,
                            hspace=0.2,wspace=0.2)
        parent.canvas[self.currplt].draw() 

class plot_samplecorr(ui_plot):
    """
    The plot_samplecorr class generates a heatmap plot of the Spearman or Pearson correlation between samples.
    
    Parameters:
    
    parent: the parent widget for the plot
    currplt: the index of the current plot within the parent widget
    frame: the parent frame for the plot
    file: a path to the file containing the ion dictionary
    filtereddfs: a dictionary containing filtered dataframes for each group in the plot
    groupsets: a dictionary containing GroupSet objects for each group in the plot
    Methods:
    
    __init__(self, parent, currplt, frame, file, filtereddfs, groupsets): initializes the plot by calling the plot() method with the given parameters
    plot(self, parent, file, filtereddfs, groupsets): generates the plot with the given data. Reads the ion dictionary from a csv file and reads the filtered data from a csv file generated by the program. Calculates the Spearman correlation matrix and generates a heatmap plot using the Seaborn library. Adjusts the layout of the plot and draws it on the parent canvas.
    """
    def __init__(self, parent, currplt, frame, file, filtereddfs, groupsets):
        ui_plot.__init__(self, parent, currplt, frame)
        self.parent = parent
        self.currplt = currplt
        self.plot(parent, file, filtereddfs, groupsets)
        
    def plot(self, parent, file, filtereddfs, groupsets):
        iondict = cached_read_csv(self.parent.analysis_paramsgui.outputdir / 'iondict.csv', sep=',', header=[0], index_col=None)
        msdata = cached_read_csv(self.parent.analysis_paramsgui.outputdir / (self.parent.analysis_paramsgui.filename.stem + '_filtered.csv'), sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
        try:
            msdata = msdata.stack([0, 1, 2], future_stack=True).groupby(level=[0, 1, 2, 3, 4]).mean().droplevel(level=3, axis=0).unstack()
        except TypeError:
            msdata = msdata.stack([0, 1, 2]).groupby(level=[0, 1, 2, 3, 4]).mean().droplevel(level=3, axis=0).unstack()
        msdata.index = msdata.index.droplevel([1, 2])
        pmatrix = msdata.corr(method='spearman')
        fig = self.parent.fig[self.currplt]
        ax = self.parent.ax[self.currplt]
        # Remove any axes left over from a previous run (notably the colorbar
        # that sns.heatmap appends). Without this a new colour-legend bar is
        # stacked onto the figure every time the plot is regenerated.
        for extra_ax in list(fig.axes):
            if extra_ax is not ax:
                extra_ax.remove()
        ax.clear()
        sns.heatmap(pmatrix, ax=ax, cmap=self.parent.analysis_paramsgui.colorscheme, vmin=0, vmax=1)
        ax.tick_params(axis='both', which='both', labelsize=10)
        ax.set_xticks(range(len(pmatrix.columns)))
        ax.set_xticklabels(pmatrix.columns, rotation=90)
        ax.set_yticks(range(len(pmatrix.index)))
        ax.set_yticklabels(pmatrix.index, rotation=0)
        ax.axes.get_xaxis().get_label().set_visible(False)
        ax.axes.get_yaxis().get_label().set_visible(False)
        self.parent.fig[self.currplt].subplots_adjust(left=.1, right=.95, bottom=0.15, top=0.9, hspace=0.2, wspace=0.2)
        self.parent.canvas[self.currplt].draw()
       
    
class kendrick(ui_plot): 
    """
    The purpose of this class is to plot the mass defect versus the nominal mass of compounds based on the input files and parameters provided.

    The __init__ method initializes the object with the following parameters:
    
    parent: The parent object of the current class instance.
    currplt: The current plot index.
    frame: The frame object that contains the plot.
    file: The input file to be used for plotting.
    filtereddfs: A dictionary containing filtered dataframes.
    groupsets: A dictionary containing groupsets.
    """
    def __init__(self, parent, currplt, frame, file, filtereddfs, groupsets):
        super().__init__(parent, currplt, frame)
        self.parent = parent
        self.currplt = currplt
        self.plot(parent, file, filtereddfs, groupsets)

    def plot(self, parent, file, filtereddfs, groupsets):
        """
        The method first reads the input file as a pandas dataframe iondict. If the blnkfltr parameter in parent.analysis_paramsgui is True, then the dataframes in filtereddfs are filtered to exclude blank data.
        For each element in filtereddfs, the nominal mass (m/z) and mass defect (kmd) data is plotted as a scatter plot. The color of each point is determined by the plotcol parameter in groupsets. The legendname parameter in groupsets is used to label the points in the plot.
        If the mdguide parameter in parent.analysis_paramsgui is True, then reference lines are plotted.
        less
        """
        iondict = pd.read_csv(file, sep=',', header=0, index_col=None)

        if parent.analysis_paramsgui.blnkfltr:
            for elem in filtereddfs:
                filtereddfs[elem] = filtereddfs[elem][filtereddfs[elem]['pass_blnkfil']]

        for elem in filtereddfs:
            plotcol = groupsets[elem].plotcol.lstrip('#')
            rgbcol = tuple(int(plotcol[i:i+2], 16) / 255 for i in (0, 2, 4))
            rgbacol = [(rgbcol[0], rgbcol[1], rgbcol[2], a) for a in filtereddfs[elem]['logmax'] / parent.analysis_paramsgui.maxval]
            rgbacol = np.clip(rgbacol, 0, 1)
            parent.ax[self.currplt].scatter(filtereddfs[elem]['m/z'], filtereddfs[elem]['kmd'], c=rgbacol, label=str(groupsets[elem].legendname), picker=True)

        parent.highlight[self.currplt], = parent.ax[self.currplt].plot([], [], 'o', markersize=12, color='yellow')

        parent.ax[self.currplt].set_xlabel('m/z', **self.fcsfont)
        parent.ax[self.currplt].set_ylabel('Mass Defect', **self.fcsfont)
        parent.ax[self.currplt].set_ylim(0, 1)
        parent.ax[self.currplt].grid()
        parent.ax[self.currplt].legend()
        self.event = parent.canvas[self.currplt].figure.canvas.mpl_connect('pick_event', lambda event: self.onpick(event, parent, iondict, ('m/z', 'kmd')))
        parent.fig[self.currplt].subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.9, hspace=0.2, wspace=0.2)
        parent.canvas[self.currplt].draw()

        if parent.analysis_paramsgui.mdguide:
            parent.ax[self.currplt].plot([0, 733.314], [0.00112, 0.8408], color='dimgrey', linestyle='-', linewidth=1)
            parent.ax[self.currplt].plot([63.6623, 733.314], [1, 0.8408], color='dimgrey', linestyle='-', linewidth=1)
            parent.ax[self.currplt].plot([0, 465.456], [0.9884, 0.5408], color='dimgrey', linestyle='-', linewidth=1)

class plot_volcano(ui_plot): 
    """
    volcano plot of -log10p vs log2 abundance value
    """
    def __init__(self, parent, currplt, frame, file, filtereddfs, groupsets):
        super().__init__(parent, currplt, frame)
        self.parent = parent
        self.currplt = currplt
        self.plot(parent, file, filtereddfs, groupsets)

    def plot(self, parent, file, filtereddfs, groupsets):
        pqvar = '-logq' if parent.analysis_paramsgui.FDR else '-logp'
        parent.ui.lbl_volcanowarn.setText('' if parent.analysis_paramsgui.FDR else 'False discovery rate correction off')

        querylist = parent.analysis_paramsgui.querylist
        iondict = pd.read_csv(file, sep = ',', header = [0], index_col = None)
        

        for elem in querylist:
            if parent.analysis_paramsgui.blnkfltr: #this uses ionfilter object which may be deprecated long term, see comments on this class in stats module
                filtereddfs[elem] = filtereddfs[elem][filtereddfs[elem]['pass_blnkfil']]#check why this gets called in conditional? should be true if blankfil is off, does this even need to be here if these are already filtered?
            filtereddfs[elem] = listfilter(iondict, groupsets[elem].ionlist, True)
            filtereddfs[elem]['logfc'] = np.log2(filtereddfs[elem]['fc']) #not sure if best to move this bit to stats?
        iondict['logfc'] = np.log2(iondict['fc'])
        iondict.to_csv(parent.analysis_paramsgui.outputdir / 'iondict.csv', header = True, index = False)
        # iondict.csv just changed on disk -- drop every cached read of it so
        # other code (e.g. fillfttree(), _refresh_highlight()) doesn't keep
        # serving a stale pre-volcano-plot copy under some other cache key.
        invalidate_csv_cache()

        maxpval = 0
        for elem in filtereddfs:
            if parent.analysis_paramsgui.blnkfltr:
                filtereddfs[elem] = filtereddfs[elem][filtereddfs[elem]['pass_blnkfil']]
            max_val = filtereddfs[elem][pqvar].max()
            if max_val > maxpval:
                maxpval = max_val
                
            sig = filtereddfs[elem][filtereddfs[elem][pqvar] >= -np.log10(parent.analysis_paramsgui.pqthresh)]
            sig = sig[sig['logfc'].abs() >= np.log2(parent.analysis_paramsgui.fcthresh)]
            nonsig = filtereddfs[elem][~filtereddfs[elem]['Compound'].isin(sig['Compound'].to_list())]
            
            parent.ax[self.currplt].scatter(sig[sig['logfc'] > 0]['logfc'], sig[sig['logfc'] > 0][pqvar], color='red', picker=True, alpha=0.5)
            parent.ax[self.currplt].scatter(sig[sig['logfc'] <= 0]['logfc'], sig[sig['logfc'] <= 0][pqvar], color='blue', picker=True, alpha=0.5)
            parent.ax[self.currplt].scatter(nonsig['logfc'], nonsig[pqvar], color='black', picker=True, alpha=0.5)
        
        parent.highlight[self.currplt], = parent.ax[self.currplt].plot([], [], 'o', markersize=12, color='yellow')
        
        parent.ax[self.currplt].plot([-np.log2(parent.analysis_paramsgui.fcthresh), -np.log2(parent.analysis_paramsgui.fcthresh)], [0, maxpval*1.2], color='dimgrey', linestyle='-', linewidth=1)
        parent.ax[self.currplt].plot([np.log2(parent.analysis_paramsgui.fcthresh), np.log2(parent.analysis_paramsgui.fcthresh)], [0, maxpval*1.2], color='dimgrey', linestyle='-', linewidth=1)
        parent.ax[self.currplt].plot([-6.75, 6.75], [-np.log10(parent.analysis_paramsgui.pqthresh), -np.log10(parent.analysis_paramsgui.pqthresh)], color='dimgrey', linestyle='-', linewidth=1)
        
        parent.ax[self.currplt].set_xlabel("log2 fold change", **self.fcsfont)
        parent.ax[self.currplt].set_ylabel('-log10 ' + pqvar[-1] + '-value', **self.fcsfont)
        parent.ax[self.currplt].set_xlim(-6.75, 6.75)
        
        try:
            parent.ax[self.currplt].set_ylim(0, maxpval*1.1)
        except Exception:
            parent.error('Warning: No volcano plot features generated, check test groups')
            pass
        
        parent.ax[self.currplt].grid()
        
        def on_pick_event(event):
            self.onpick(event, parent, iondict, ('logfc', '-logq'))
        
        self.event = parent.canvas[self.currplt].figure.canvas.mpl_connect('pick_event', on_pick_event)
        
        parent.fig[self.currplt].subplots_adjust(left=.1, right=.95, bottom=0.15, top=0.9, hspace=0.2, wspace=0.2)
        parent.canvas[self.currplt].draw()
 
class plot_fc3d(ui_plot):
    """
    3D mass, retention time, fold change plot, not currently highlightable
    """

    def __init__(self, parent, currplt, frame, file, filtereddfs, groupsets):
        ui_plot.__init__(self, parent, currplt, frame)
        self.parent = parent
        self.currplt = currplt
        self.plot(parent, file, filtereddfs, groupsets)

    def plot(self, parent, file, filtereddfs, groupsets):
        iondict = pd.read_csv(file, sep=',', header=[0], index_col=None)
        for elem in parent.analysis_paramsgui.querylist:
            if parent.analysis_paramsgui.blnkfltr:
                filtereddfs[elem] = filtereddfs[elem][filtereddfs[elem]['pass_blnkfil']]
            filtereddfs[elem] = listfilter(iondict, groupsets[elem].ionlist, True)
            filtereddfs[elem]['logfc'] = np.log2(filtereddfs[elem]['fc'])

        parent.ax[self.currplt].remove()
        parent.ax[self.currplt] = parent.canvas[self.currplt].figure.add_subplot(111, projection='3d')
        parent.ax[self.currplt].set_axisbelow(True)
        for elem in filtereddfs:
            x = filtereddfs[elem]['logfc']
            z = filtereddfs[elem]['m/z']
            y = filtereddfs[elem]['Retention time (min)']
            parent.ax[self.currplt].scatter(x, y, z, color=str(groupsets[elem].plotcol), marker='o', label=str(groupsets[elem].legendname))
        parent.canvas[self.currplt].figure.set_facecolor(self.plotbackground)
        parent.ax[self.currplt].set_facecolor(self.plotbackground)
        parent.ax[self.currplt].set_xlabel('log2 fold change', **self.fcsfont, labelpad=20)
        parent.ax[self.currplt].set_ylabel('Retention time (min)', **self.fcsfont, labelpad=20)
        parent.ax[self.currplt].set_zlabel('m/z', **self.fcsfont, labelpad=25)
        parent.ax[self.currplt].tick_params(axis='z', pad=10)
        parent.ax[self.currplt].set_ylim(-0.5, 11.5)
        parent.ax[self.currplt].grid()
        parent.ax[self.currplt].legend()
        parent.fig[self.currplt].subplots_adjust(
            left=0.1, right=0.95, bottom=0.15, top=0.9, hspace=0.2, wspace=0.2)
        parent.canvas[self.currplt].draw()

class plot_dendrogram(ui_plot):
    """
    Dendrogram generation, with a combo-box switcher (same pattern as
    plot_ordination's method/view bar) between two purity-colored views:

    - "Technical Replicates": every injection is its own leaf, colored
      green wherever an entire Sample's injections cluster together before
      merging with anything else -- a quick visual QC for whether technical
      replicates are tight.
    - "Biological Replicates": injections are first averaged per Sample
      (same collapsing logic as the ordination tab's "Collapse Technical
      Replicates" checkbox, via ordination.load_ordination_matrix), then
      leaves are colored green wherever an entire Biolgroup's samples
      cluster together -- a quick visual QC for whether biological groups
      are separable at all, independent of technical noise.

    Either view can be regular or bootstrapped (PvClust), depending on
    parent.analysis_paramsgui.bootstrap, same as before this rework. The
    purity-coloring math lives in the Qt-free clusterpurity.py.
    """

    VIEWS = ('Technical Replicates', 'Biological Replicates')

    def __init__(self, parent, currplt, frame, file, filtereddfs, groupsets):
        ui_plot.__init__(self, parent, currplt, frame)
        self.parent = parent
        self.currplt = currplt
        # Default matches the plot's previous (injection-level) behaviour
        # exactly, so existing sessions see no change until they explicitly
        # switch to the biological-replicate view.
        self.view = 'Technical Replicates'
        self._build_switcher_bar(parent, currplt)
        self.plot(parent, file, filtereddfs, groupsets)

    def _build_switcher_bar(self, parent, currplt):
        bar = QtWidgets.QWidget()
        bar.setStyleSheet(_SWITCHER_BAR_STYLE)
        bar.setMaximumHeight(_SWITCHER_BAR_HEIGHT)
        layout = QtWidgets.QHBoxLayout(bar)
        layout.setContentsMargins(4, 2, 4, 2)

        layout.addWidget(QtWidgets.QLabel('View:'))
        view_combo = QtWidgets.QComboBox()
        view_combo.addItems(self.VIEWS)
        view_combo.setCurrentText(self.view)
        view_combo.currentTextChanged.connect(self._on_view_changed)
        layout.addWidget(view_combo)
        layout.addStretch()

        self.view_combo = view_combo
        parent.pltlayout[currplt].insertWidget(0, bar)

    def _on_view_changed(self, view):
        self.view = view
        self.reset(self._last_file, self._last_filtereddfs, self._last_groupsets)

    def plot(self, parent, file, filtereddfs, groupsets):
        self._last_file = file
        self._last_filtereddfs = filtereddfs
        self._last_groupsets = groupsets

        # PvClust (bootstrap path) expects "variables x objects" -- it
        # bootstraps over the rows (features) and transposes internally
        # before clustering the columns (the objects/leaves). shc.linkage
        # (regular path) expects the opposite, "objects x variables" --
        # build both orientations from the same scaled data below.
        if self.view == 'Biological Replicates':
            # Collapse technical replicates first -- leaves are Samples,
            # purity is judged against Biolgroup.
            raw_header = cached_read_csv(
                parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_filtered.csv'),
                sep=',', header=None, index_col=[0, 1, 2]).iloc[:3, :].transpose()
            x, biolgroup = ordination.load_ordination_matrix(file, raw_header.copy(), collapse_replicates=True)
            data_scaled = normalize(x.values, axis=1)  # normalize each sample's profile
            data_scaled = pd.DataFrame(data_scaled, columns=x.columns, index=x.index)  # samples x features
            textlabels = data_scaled.index.tolist()
            leaf_labels = [biolgroup[sample] for sample in textlabels]
            data_for_linkage = data_scaled
            data_for_pvclust = data_scaled.transpose()
            purity_noun = 'biological groups separable'
        else:
            heirarch = pd.read_csv(file, sep=',', header=[2], index_col=[0]).drop(['m/z', 'Retention time (min)'], axis=1)
            data_scaled = normalize(heirarch, axis=0)  # normalize features
            data_scaled = pd.DataFrame(data_scaled, columns=heirarch.columns, index=heirarch.index)  # features x injections
            textlabels = data_scaled.columns.tolist()
            raw_header = cached_read_csv(
                parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_filtered.csv'),
                sep=',', header=None, index_col=[0, 1, 2]).iloc[:3, :].transpose()
            raw_header.columns = ['Biolgroup', 'Sample', 'Injection']
            sample_of_injection = raw_header.set_index('Injection')['Sample'].to_dict()
            leaf_labels = [sample_of_injection[name] for name in textlabels]
            data_for_linkage = data_scaled.transpose()
            data_for_pvclust = data_scaled
            purity_noun = "samples' replicates clustered together"

        if parent.analysis_paramsgui.bootstrap:
            # bootstrap dendrogram
            pv = PvClust(data_for_pvclust, method="ward", metric="euclidean", nboot=1000, parallel=True)
            link_color_func = clusterpurity.purity_link_color_func(pv.linkage_matrix, leaf_labels)
            dend = pv.plot(parent.ax[self.currplt], labels=textlabels, link_color_func=link_color_func)
            Z = pv.linkage_matrix
        else:
            # regular dendrogram
            Z = shc.linkage(data_for_linkage, method='ward')
            link_color_func = clusterpurity.purity_link_color_func(Z, leaf_labels)
            dend = shc.dendrogram(Z, ax=parent.ax[self.currplt], leaf_rotation=90, above_threshold_color='black', link_color_func=link_color_func, labels=textlabels)  # default leaf label size 16

        n_pure, n_total = clusterpurity.purity_summary(Z, leaf_labels)
        parent.ax[self.currplt].set_title(f'{n_pure}/{n_total} {purity_noun}', fontsize=10)

        parent.fig[self.currplt].subplots_adjust(
            left=0.1, right=0.95, bottom=0.35, top=0.9, hspace=0.2, wspace=0.2)
        parent.canvas[self.currplt].draw()

# Shared by plot_dendrogram's and plot_ordination's combo-box switcher bars
# -- page_dend and page_pca both have the same light background
# (rgba(225,225,225,255), see ui_main.py), unlike searchtree.py's filter bar
# (a dark-themed tab) -- dark text on a light/white combo box, not
# searchtree's light-on-dark scheme.
_SWITCHER_BAR_HEIGHT = 32

_SWITCHER_BAR_STYLE = """
QWidget {
    background: transparent;
}
QComboBox {
    background-color: rgb(255,255,255);
    color: rgb(30,30,30);
    border: 1px solid rgb(150,150,150);
    border-radius: 2px;
    padding: 2px;
}
QLabel {
    color: rgb(30,30,30);
    background: transparent;
}
"""


class plot_ordination(ui_plot):
    """Multivariate ordination plot: PCA, NMDS, or PLS-DA, with a
    scores-vs-loadings view toggle.

    A combo-box switcher bar (built once in ``__init__``, inserted above the
    canvas the same way ``SearchTreePanel``'s filter bar is substituted into
    a Designer placeholder -- see searchtree.py) lets the user pick the
    ordination method and the scores/loadings view; both redraw onto the
    same axes via ``self.plot(...)`` rather than rebuilding the canvas.

    The actual math lives in the Qt-free ``ordination.py`` (PCA/NMDS/PLS-DA,
    technical-replicate collapsing, top-N loadings selection); this class is
    just the Qt plumbing and rendering on top of it.
    """

    METHODS = ('NMDS', 'PCA', 'PLS-DA')
    VIEWS = ('Scores', 'Loadings')

    def __init__(self, parent, currplt, frame, file, filtereddfs, groupsets):
        ui_plot.__init__(self, parent, currplt, frame)
        self.parent = parent
        self.currplt = currplt
        # Defaults match the plot's previous (NMDS-only, scores-only)
        # behaviour exactly, so existing sessions see no change until they
        # explicitly switch the new controls.
        self.method = 'NMDS'
        self.view = 'Scores'
        self.loadings_df = None
        self._build_switcher_bar(parent, currplt)
        self.plot(parent, file, filtereddfs, groupsets)

    def _build_switcher_bar(self, parent, currplt):
        bar = QtWidgets.QWidget()
        bar.setStyleSheet(_SWITCHER_BAR_STYLE)
        bar.setMaximumHeight(_SWITCHER_BAR_HEIGHT)
        layout = QtWidgets.QHBoxLayout(bar)
        layout.setContentsMargins(4, 2, 4, 2)

        layout.addWidget(QtWidgets.QLabel('Method:'))
        method_combo = QtWidgets.QComboBox()
        method_combo.addItems(self.METHODS)
        method_combo.setCurrentText(self.method)
        method_combo.currentTextChanged.connect(self._on_method_changed)
        layout.addWidget(method_combo)

        layout.addWidget(QtWidgets.QLabel('View:'))
        view_combo = QtWidgets.QComboBox()
        view_combo.addItems(self.VIEWS)
        view_combo.setCurrentText(self.view)
        view_combo.currentTextChanged.connect(self._on_view_changed)
        layout.addWidget(view_combo)
        layout.addStretch()

        self.method_combo = method_combo
        self.view_combo = view_combo
        parent.pltlayout[currplt].insertWidget(0, bar)

    def _on_method_changed(self, method):
        self.method = method
        self.reset(self._last_file, self._last_filtereddfs, self._last_groupsets)

    def _on_view_changed(self, view):
        self.view = view
        self.reset(self._last_file, self._last_filtereddfs, self._last_groupsets)

    def plot(self, parent, file, filtereddfs, groupsets):
        """(Re)draw the ordination plot for the current method/view.

        Args:
            parent (QWidget): Parent widget (MainWindow).
            file (str): Path to the ``_filtered.csv`` peak table.
            filtereddfs, groupsets: unused here (kept for the shared
                ``_create_or_reset``/``reset`` call signature every plot
                class follows).
        """
        parent = self.parent
        self._last_file = file
        self._last_filtereddfs = filtereddfs
        self._last_groupsets = groupsets

        collapse_replicates = parent.dialog.ui.checkBox_collapsereps.isChecked()
        raw_header = cached_read_csv(
            parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_filtered.csv'),
            sep=',', header=None, index_col=[0, 1, 2]).iloc[:3, :].transpose()
        x, biolgroup = ordination.load_ordination_matrix(file, raw_header.copy(), collapse_replicates)

        n_components = max(2, min(len(x) - 1, 10))

        colors = ['red', 'blue', 'black', 'grey', 'purple', 'orange', 'green', 'yellow', 'lime', 'plum', 'teal', 'olivedrab', 'sienna', 'maroon', 'navy', 'lightcoral', 'darkgoldenrod', 'seagreen', 'lightseagreen', 'aqua', 'lightsteelblue', 'slateblue', 'blueviolet', 'plum', 'burlywood', 'salmon', 'aquamarine', 'magenta', 'tan']
        colorpos, biolgroupmap = 0, {}
        for elem in biolgroup:
            if elem not in biolgroupmap and elem != parent.analysis_paramsgui.blnkgrp: ###### delete blank clause OR CHANGE TO THE BLNKFILTER OPTION
                biolgroupmap[elem] = colors[colorpos]
                colorpos += 1

        plot_title = None
        if self.method == 'PCA':
            scores, loadings, expvar = ordination.run_pca(x, n_components)
            axis_labels = [f'PC{i + 1} ({100 * expvar[i]:.1f}%)' for i in range(2)]
        elif self.method == 'PLS-DA':
            scores, loadings, expvar = ordination.run_plsda(x, biolgroup, n_components)
            axis_labels = [f'PLS{i + 1} ({100 * expvar[i]:.1f}%)' for i in range(2)]
        else:
            scores, expvar, stress = ordination.run_nmds(x, n_components)
            loadings = ordination.nmds_loading_proxy(x, scores)
            # NMDS doesn't canonically report percent-variance-explained the
            # way PCA/PLS-DA do (it's a rank-based embedding, not a linear
            # decomposition of the feature space) -- stress is the
            # conventional thing to report for NMDS instead.
            axis_labels = ['NMDS1', 'NMDS2']
            plot_title = f'Stress: {stress:.4f}'

        self.loadings_df = loadings
        principalDf = scores.copy()
        principalDf['Biolgroup'] = biolgroup

        if self.view == 'Loadings':
            self._plot_loadings(parent, loadings, axis_labels)
        else:
            self._plot_scores(parent, principalDf, biolgroupmap, axis_labels)

        if plot_title:
            parent.ax[self.currplt].set_title(plot_title, fontsize=10)

        parent.fig[self.currplt].subplots_adjust(left=.1, right=.9, bottom=0.1, top=0.9, hspace=0.2, wspace=0.2)
        parent.canvas[self.currplt].draw()

    def _plot_scores(self, parent, principalDf, biolgroupmap, axis_labels):
        for elem in biolgroupmap:
            scatterframe = principalDf[principalDf['Biolgroup'] == elem]
            points = scatterframe.iloc[:, [0, 1]].to_numpy()
            if np.shape(points)[0] > 2:
                self.plot_point_cov(points, nstd=2, ax=parent.ax[self.currplt], alpha=0.5, color=self.lighten_color(biolgroupmap[elem], 0.3))
            parent.ax[self.currplt].scatter(scatterframe.iloc[:, 0], scatterframe.iloc[:, 1], color=biolgroupmap[elem], marker='o', s=30, label=str(elem), picker=True)

        parent.highlight[self.currplt], = parent.ax[self.currplt].plot([], [], 'o', markersize=12, color='yellow')
        parent.ax[self.currplt].set_xlabel(axis_labels[0], **self.fcsfont)
        parent.ax[self.currplt].set_ylabel(axis_labels[1], **self.fcsfont)

        self.highlightcol = (0, 0, 0, 0)
        parent.pickedsample = pd.DataFrame(0, index=['empty'], columns=['empty'])

        def picksample(event):
            if _is_duplicate_pick(parent, event):
                return
            ind = event.ind
            coord = event.artist.get_offsets()[ind, :]
            newsample = principalDf.loc[principalDf.iloc[:, 0] == coord[0, 0], :].loc[principalDf.iloc[:, 1] == coord[0, 1], :].reset_index()
            if newsample.empty:
                return

            if newsample.iloc[0, 0] == parent.pickedsample.iloc[0, 0] and self.highlightcol != (0, 0, 0, 0):
                self.highlightcol = (0, 0, 0, 0)
            else:
                self.highlightcol = 'yellow'

            parent.pickedsample = newsample
            parent.ui.lbl_injname.setText('Injection/Sample: ' + str(parent.pickedsample.iloc[0, 0]))
            parent.highlight[self.currplt].set_data(coord[0, 0], coord[0, 1])
            parent.highlight[self.currplt].set_color(self.highlightcol)
            parent.canvas[self.currplt].draw_idle()

        self.event = parent.canvas[self.currplt].figure.canvas.mpl_connect('pick_event', picksample)
        parent.ax[self.currplt].legend()

    def _plot_loadings(self, parent, loadings, axis_labels):
        """Loadings (biplot-style) view: origin-anchored arrows for the
        top-N features by vector magnitude, plus -- regardless of
        magnitude -- whichever feature is currently highlighted elsewhere
        in the app (``parent.pickedfeature``), so a feature too small to
        make the default cut is still visible on demand.
        """
        always_include = [parent.pickedfeature] if getattr(parent, 'pickedfeature', '') else []
        # Rank by magnitude within the 2 displayed components only, not the
        # full (up to 10-component) loadings -- a feature could rank in the
        # overall top-25 purely from a large contribution to some other,
        # unplotted component while barely showing up here, displacing a
        # feature that's actually prominent in this 2D view.
        subset = ordination.top_loadings(loadings.iloc[:, :2], n=25, always_include=always_include)

        # ax.annotate()'s arrows don't reliably drive matplotlib's autoscale
        # the way ax.scatter()/ax.plot() do (confirmed empirically: points
        # can end up outside the auto-picked view limits), so the axis
        # range is set explicitly here instead of relying on autoscale.
        # Symmetric around 0 since loadings/correlations are naturally
        # origin-centered (a biplot convention).
        limit = subset.iloc[:, :2].abs().values.max() * 1.2 if len(subset) else 1.0
        parent.ax[self.currplt].set_xlim(-limit, limit)
        parent.ax[self.currplt].set_ylim(-limit, limit)

        for feature, row in subset.iterrows():
            xcoord, ycoord = row.iloc[0], row.iloc[1]
            parent.ax[self.currplt].annotate(
                '', xy=(xcoord, ycoord), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='steelblue', lw=1))
            parent.ax[self.currplt].annotate(
                str(feature), xy=(xcoord, ycoord), fontsize=8, color='black')

        # Pre-created empty artist for the highlighted-loading marker,
        # following the same convention as the scores view's
        # parent.highlight[currplt] -- updated on demand by
        # MainWindow._refresh_highlight() via self.highlight_loading(),
        # even when the highlighted feature isn't in the default top-25.
        self.loadings_highlight, = parent.ax[self.currplt].plot([], [], 'o', markersize=12, color='yellow', zorder=5)
        self.highlight_loading(getattr(parent, 'pickedfeature', ''), getattr(parent, 'highlightcol', (0, 0, 0, 0)))

        parent.ax[self.currplt].axhline(0, color='grey', lw=0.5)
        parent.ax[self.currplt].axvline(0, color='grey', lw=0.5)
        parent.ax[self.currplt].set_xlabel(axis_labels[0], **self.fcsfont)
        parent.ax[self.currplt].set_ylabel(axis_labels[1], **self.fcsfont)

    def highlight_loading(self, feature, colour):
        """Update the loadings-view highlight marker for ``feature`` (a
        no-op outside the loadings view or before it's been drawn once).

        Called from ``MainWindow._refresh_highlight()`` -- the same
        pre-create-empty-artist/update-via-set_data convention every other
        plot's highlight already follows, just driven by this plot's own
        last-computed loadings instead of ``iondict``.
        """
        if self.view != 'Loadings' or self.loadings_df is None or not hasattr(self, 'loadings_highlight'):
            return
        if not feature or feature not in self.loadings_df.index:
            self.loadings_highlight.set_data([], [])
        else:
            row = self.loadings_df.loc[feature]
            self.loadings_highlight.set_data([row.iloc[0]], [row.iloc[1]])
            self.loadings_highlight.set_color(colour)
        self.parent.canvas[self.currplt].draw_idle()

    def plot_point_cov(self, points, nstd=2, ax=None, **kwargs):
        """Generate an ellipse for the confidence interval.

        Args:
            points (numpy.ndarray): Array of points to calculate the ellipse for.
            nstd (int): Number of standard deviations to calculate the ellipse for. Default is 2.
            ax (matplotlib.axes.Axes): The Axes instance on which to plot the ellipse. Default is None.
            **kwargs: Additional keyword arguments to be passed to the matplotlib Ellipse constructor.
        
        Returns:
            matplotlib.patches.Ellipse: Ellipse for the confidence interval.
        """
        pos = points.mean(axis=0)
        cov = np.cov(points, rowvar=False)
        return self.plot_cov_ellipse(cov, pos, nstd, ax, **kwargs)

    def lighten_color(self, color, amount=0.5):
        """Lighten a given color by a given amount.

        Args:
            color (str or tuple): Color to be lightened, either as a string or a tuple of RGB values.
            amount (float): Amount to lighten the color by. Default is 0.5.

        Returns:
            tuple: Tuple of RGB values for the lightened color.
        """
        try:
            c = mc.cnames[color]
        except:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

    def plot_cov_ellipse(self, cov, pos, nstd=2, ax=None, **kwargs):
        """Generate an optimized ellipse for the confidence interval.

        Args:
            cov (numpy.ndarray): Array containing the covariance matrix.
            pos (numpy.ndarray): Array containing the mean values.
            nstd (int): Number of standard deviations to calculate the ellipse for. Default is 2.
            ax (matplotlib.axes.Axes): The Axes instance on which to plot the ellipse. Default is None.
            **kwargs: Additional keyword arguments to be passed to the matplotlib Ellipse constructor.

        Returns:
            matplotlib.patches.Ellipse: Ellipse for the confidence interval.
        """
        def eigsorted(cov):
            vals, vecs = np.linalg.eigh(cov)
            order = vals.argsort()[::-1]
            return vals[order], vecs[:,order]

        if ax is None:
            ax = plt.gca()

        vals, vecs = eigsorted(cov)
        theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
        width, height = 2 * nstd * np.sqrt(vals)
        ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)
        ax.add_artist(ellip)
        return ellip
 
class prev_cv(ui_plot):
    """
    Plot the CV rarefaction plot with mean and median CV.
    
    Args:
    parent: A parent object representing the main window.
    currplt: A string representing the current plot.
    frame: A widget frame containing the plot.
    file: A string representing the file name.
    filtereddfs: A dictionary containing filtered data.
    groupsets: A dictionary containing groups of sets.
    
    Returns:
    None
    """
    def __init__(self, parent, currplt, frame, file, filtereddfs, groupsets):
        super().__init__(parent, currplt, frame)
        self.parent = parent
        self.currplt = currplt
        self.plot(parent, file, filtereddfs, groupsets)

    def plot(self, parent, file, filtereddfs, groupsets):
        # Load and filter ion data
        iondict = cached_read_csv(parent.analysis_paramsgui.outputdir / 'iondict.csv', header=0, index_col=0)
        iondict = iondict[~np.isnan(iondict['average CV'])]

        # Calculate mean and median CV, and scale data
        iondictmean = iondict.sort_values(['average CV']).reset_index()
        iondictmed = iondict.sort_values(['median CV']).reset_index()
        iondictmean = iondictmean.reset_index()
        iondictmed = iondictmed.reset_index()
        iondictmean.iloc[:,0] = 100 * iondictmean.iloc[:,0]/len(iondictmean['average CV'])
        iondictmed.iloc[:,0] = 100 * iondictmed.iloc[:,0]/len(iondictmed['median CV'])

        # Calculate maximum theoretical CV based on neff
        msdata_header = cached_read_csv(parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_filtered.csv'), sep=',', header=None, index_col=[0,1,2]).iloc[:3,:].transpose()
        msdata_header.columns = ['Biolgroup', 'Sample', 'Injection']
        average_n = msdata_header['Injection'].nunique() / msdata_header['Sample'].nunique()
        modelstdevlist = [1] + [0] * (int(average_n) - 1)
        modelstdev = pd.Series(modelstdevlist).std() / pd.Series(modelstdevlist).mean()
        cv50 = iondictmean.iloc[(iondictmean.iloc[:,0] - 50).abs().argsort()[:1]]['average CV']
        sortedcv = iondictmean.iloc[(iondictmean.iloc[:,0]).argsort()]['average CV']
        prevav = 0
        aucav = 0
        prevmed = 0
        aucmed = 0
        for pos in range(0,len(iondictmean.iloc[:,0])):
            dist = iondictmean.iloc[pos,:]['average CV'] - prevav
            aucav += dist*iondictmean.iloc[pos,0]
            prevav = iondictmean.iloc[pos,:]['average CV']
            
            dist = iondictmed.iloc[pos,:]['median CV'] - prevmed
            aucmed += dist*iondictmed.iloc[pos,0]
            prevmed = iondictmed.iloc[pos,:]['median CV']
            
        meanav = 0
        meanmed = 0
        sumskew = 0
        if math.isnan(modelstdev):
            modelstdev = 1.7
        for val in range(1, int((modelstdev*100))):
            pos = val/100
            meanav = iondictmean[abs(iondictmean['average CV'] - pos-modelstdev/200) < modelstdev/200].iloc[:,0].mean()
            meanmed = iondictmed[abs(iondictmed['average CV'] - pos-modelstdev/200) < modelstdev/200].iloc[:,0].mean()
            skew = abs(meanmed-meanav)
            if not np.isnan(skew):
                sumskew += skew * modelstdev/100


        sumskew = sumskew/ ((aucmed+aucav)/2)
        rep = ((aucmed+aucav)/2)/(modelstdev*100)
        qualscore = (1-sumskew)*rep*100
        
        #qualscore = round(100 * (1 - cv50 / modelstdev), 1)

        # Update UI
        parent.ui.lbl_spllist_3.setText('Reproducibility:\n' + str(round(100*rep,1))  + '%\n' +
                                        'Skewnewss:\n' + str(round(100*sumskew,1))  + '%\n\n' +
                                        'Overall:\n' + str(round(qualscore,1))  + '%')

        # Plot data
        currplt = 'cvplt' #instead take this from input
        parent.ax[currplt].plot(iondictmed['median CV'], iondictmed.iloc[:,0], color = '#0000ff', label="Median CV")
        parent.ax[currplt].plot(iondictmean['average CV'], iondictmed.iloc[:,0], color = 'red', label="Mean CV")
        parent.ax[currplt].set_xlabel("CV",  **self.fcsfont)
        parent.ax[currplt].set_ylabel('Percentage of Features',  **self.fcsfont)
        parent.ax[currplt].legend()
        parent.ax[currplt].set_xlim([0, modelstdev])
        parent.ax[currplt].set_ylim([0, 100])
        parent.ax[currplt].plot([parent.analysis_paramsgui.cvthresh, parent.analysis_paramsgui.cvthresh], [0, 100], color='dimgrey', linestyle='-', linewidth=1)
        threshpercent = iondictmed.iloc[(iondictmed[parent.analysis_paramsgui.cvparam]-parent.analysis_paramsgui.cvthresh).abs().argsort()[:1],0]
        parent.ax[currplt].plot([0, modelstdev], [threshpercent, threshpercent] , color='dimgrey', linestyle='-', linewidth=1)
            
        
        parent.fig[currplt].subplots_adjust(left=.15,right=0.9,
                            bottom=.15,top=.85,
                            hspace=10,wspace=10)
        parent.canvas[currplt].draw()

    
def gen_upsetplt(parent):    #need to do something to handle groups with names that are substrings of other group names
    """
    Generate an upset plot to visualize sets of compounds in groups. This function also handles groups with names that are substrings of other group names.

    Parameters:
    parent (object): The parent object that the generated plot will be a child of.

    Returns:
    None
    """
    iondict = cached_read_csv(parent.analysis_paramsgui.outputdir / 'iondict.csv', sep=',', header=0, index_col=None)
                            
    # Apply filters if required
    if parent.analysis_paramsgui.relfil:
        iondict = iondict[iondict['pass_relfil']]
    if parent.analysis_paramsgui.decon:
        iondict = iondict[iondict['pass_insource']]
    if parent.analysis_paramsgui.blnkfltr:
        iondict = iondict[iondict['pass_blnkfil']]
    if parent.analysis_paramsgui.CVfil:
        iondict = iondict[iondict['pass_cvfil']]
    
    # Prepare data for upset plot
    iongroups = iondict['groups'].tolist()
    freq = {}
    biolgroups = []
    for item in iongroups:
        if item not in freq:
            freq[item] = 0
        freq[item] += 1
    
    header = cached_read_csv(parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_filtered.csv'), sep=',', header=None, index_col=[0, 1, 2]).iloc[0, :]
    for elem in header:
        if elem not in biolgroups:
            biolgroups.append(elem)
    
    sets = [' ' + elem for elem in list(freq.keys())]
    size = list(freq.values())
    setdf = pd.DataFrame({'groups': sets})
    for elem in biolgroups:  #have to do this if one group is a substring of another, add space
        setdf[elem] = setdf['groups'].str.contains(' ' + elem)
    setdf['size'] = size
    setdf = setdf.iloc[:, 1:]
    setdf = setdf.set_index(biolgroups)['size']
    
    # Plot and display the upset plot
    with plt.rc_context({"font.size": 8}):
        upsetplt = upsetplot.plot(setdf, show_counts='%d', show_percentages=True, sort_categories_by=None)
    
    figup = upsetplt['matrix'].figure
    figup.set_size_inches(5, 4)
    figup.set_facecolor((0, 0, 0, 0))
    upsetplt['intersections'].set_facecolor((1, 1, 1, .25))
    figup.savefig('test_upsetplt.png', dpi=150, bbox_inches='tight')
    pixmap = QPixmap('test_upsetplt.png')
    parent.ui.label_upset.setPixmap(pixmap)                                                                     

def gen_treemap(parent):
    #generate treemap for visualization of filtering levels
    #needed to refilter data and see how df row lengths change to avoid issues with one feature being in multiple filter lists
    """
    The gen_treemap function generates a treemap for visualizing filtering levels. The function reads a CSV file containing the
    filtered data and another CSV file containing information about the ions. The function then filters the ion data based on 
    various filter options and calculates the number of ions filtered by each filter. Finally, the function generates a treemap 
    to display the number of ions that passed each filter and saves it as a PNG file. The treemap is then displayed in a QLabel in the GUI.

    Args:
    
    parent: the parent widget where the treemap will be displayed

    """
    plt.clf()
    msdata_filtered = cached_read_csv(parent.analysis_paramsgui.outputdir / (parent.analysis_paramsgui.filename.stem + '_filtered.csv'), sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
    fltrcnt, color = {}, []
    iondict = cached_read_csv(parent.analysis_paramsgui.outputdir / 'iondict.csv', sep=',', header=[0], index_col=[0])
    total = len(iondict.index)
    current = total
    
    if parent.analysis_paramsgui.relfil:
        filteredsetsize = len(iondict[iondict['pass_relfil']].index)
        fltrcnt['Mispicked'] = current - filteredsetsize
        current = filteredsetsize
        color.append('#0000ff')
    
    if parent.analysis_paramsgui.blnkfltr:
        filteredsetsize = len(iondict[iondict['pass_blnkfil']].index)
        fltrcnt['Blank'] = current - filteredsetsize
        current = filteredsetsize
        color.append('#00aaaa')
    
    if parent.analysis_paramsgui.CVfil:
        fltrcnt['Nonreproducible'] = len(parent.ionfilters['cv'].ions)
        current = current - fltrcnt['Nonreproducible']
        color.append('#ff0000')
    
    if parent.analysis_paramsgui.decon:
        fltrcnt['Insource'] = len(parent.ionfilters['insource'].ions)
        color.append('#00aa00')
    
    fltrcnt['High Quality'] = len(msdata_filtered.index)
    color.append('#000000')

    sizes = list(fltrcnt.values())
    total_size = sum(fltrcnt.values())
    labels = [f"{label}\n{size}\n{round(100*size/total_size,1)}%" for label, size in fltrcnt.items()]

    squarify.plot(sizes=sizes, label=labels, color=color, alpha=0.3, text_kwargs={'fontsize': 10})
    plt.axis('off')
    plt.savefig('treemap.png', dpi=150, bbox_inches='tight')
    pixmap = QPixmap('treemap.png')
    parent.ui.label_treemap.setPixmap(pixmap)