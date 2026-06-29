"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
"""

import re

import sys
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import time
import string

import platform
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QSizeGrip, QGraphicsDropShadowEffect, QFileDialog, QListWidgetItem, QColorDialog
from PyQt5.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect, QSize, QTime, QUrl, Qt, QEvent)
from PyQt5.QtGui import QBrush, QColor, QIcon, QPalette, QPainter, QPixmap
from pathlib import Path

# Install/verify non-stock dependencies (epam.indigo, UpSetPlot, squarify)
# before the modules that import them are loaded below.
from importdependencies import ensure_dependencies
ensure_dependencies()

# GUI FILE
from ui_main import Ui_MainWindow
import ui_plotparam
import ui_featureinfo

# IMPORT FUNCTIONS AND RESOURCES
from ui_functions import *
import files

from MSFaST import run_MSFaST, analysis_parameters
from groupsets import GroupSet, GroupSetModel, build_query_dict
from plotslots import PlotSlotRegistry
from paramfields import save_checkbox_fields
from csvcache import cached_read_csv, invalidate as invalidate_csv_cache
from biogroups import compute_biological_groups
from dbsearch import search_npatlas
from searchtree import SearchTreePanel
from plotting import plot_abund, show_spectrum, show_featureplt, plot_heatmap, plot_mzrt, plot_samplecorr, kendrick, plot_volcano, plot_fc3d, plot_dendrogram, plot_PCA, prev_cv, gen_upsetplt, gen_treemap
import getfragdb

from indigo import Indigo
from indigo.renderer import IndigoRenderer
indigo = Indigo()
renderer = IndigoRenderer(indigo)
import os

import pickle


import pathlib
if sys.platform == 'win32':
    temp = pathlib.PosixPath
    pathlib.PosixPath = pathlib.WindowsPath
else:
    temp = pathlib.WindowsPath
    pathlib.WindowsPath = pathlib.PosixPath

low_memory=False

'''
#need to make sure rel vs absolute is stored, and that merge checkbox connected, add to load filey
Check if low_memory=False increases ram usage for average grps?

-add bypass for plots based on checkmark. possibly use if check: ... else: button.hide() then pass

- distribution of CVs on bottom of cvplt?
- add pca option and allow visualization of key features on multivar plt?

#TODO#
#Easy/high priority
- in source spectra
- do overall data quality score, AUC
- standardize method and class names
- database management, options
- fix up analysisinfo file output

#Medium priority
- mzmine msp file import
- add other ordination options
- add custom keyword arguments for each plot to make calling them easier
- add more conditionals so if one plot fails it doesn't kill everything else
- add runcheck before searching when switching to search tab
- Figure out way to have only active plot be updated and then to update others
    when plot is switched
- change treewidget in search tab to treeview for better search, add filter options
- make it so groups can be reordered
- make it so iondict and msdata are saved as a parent object so they dont need to
    be reopened each time a heatmap feature is scrolled. This may actually not
    be a good idea depending on ram demands
- make goto buttons just one class and lambda an index for the stacked widgets
    when connecting!

#low priority/long term
- Switch to MVC format for groupsets
- maybe have a comparison mode for many different strains with and without elicitor
- specificity/sensitivity plot
- other statistical models
- refactor the filtereddfs and groupset code. everything should probably
    be filtered from the main file in the analysis by referencing the groups
    column in iondict
    ~best way to do this may be to simplify the query/groupset flow to only have
    one object, and then to generate a list of ions that pass this filter in the
    object. That or add a color column to the iondict file
    ~refactor to run filtering based of a list of filters instead of several 
    conditionals each time OR do it so that the appropriate columns are still 
    generated each time if filtering is off but they are always true
'''

class query:
    """Legacy groupset shape, kept only so old ``.mpct`` files (which pickle
    this class by qualified name ``main.query``) can still be unpickled.

    Live code uses ``groupsets.GroupSet``/``GroupSetModel`` instead; on load,
    any ``query`` objects found in a save file are converted via
    ``GroupSet.from_legacy`` (see ``MainWindow.read_save``). Do not use this
    class for new code.
    """
    def __init__(self):
        self.name = ''
        self.src = ''  # or groups, feature can be in
        self.incl = ''  # and groups, feature must be in
        self.excl = ''  # not groups, feature must not be in
        self.colour = '#000000'


class AnalysisWorker(QObject):
    """Runs the heavy, Qt-free computation (``run_MSFaST``) on a worker thread.

    Only the pure pandas/numpy/file-I/O step is moved off the GUI thread;
    matplotlib/Qt plot generation must stay on the main thread and runs in
    ``MainWindow._finish_analysis`` once ``finished`` is emitted.
    """
    finished = QtCore.pyqtSignal()
    failed = QtCore.pyqtSignal(str)

    def __init__(self, window):
        QObject.__init__(self)
        self.window = window

    def run(self):
        try:
            result = run_MSFaST(self.window.analysis_paramsgui)
            self.window.ionfilters = result.ionfilters
            self.window.groupionlists = result.groupionlists
            self.window.groupsets = result.groupsets
            self.window.filtereddfs = result.filtereddfs
        except Exception:
            import traceback
            traceback.print_exc()
            self.failed.emit('Processing failed (see console)')
            return
        self.finished.emit()


fullruntime = 0


def start_functime():
    """
    Used to calculate runtime
    """
    global initfunctime
    initfunctime = time.time()


def stop_functime(text):
    """
    Prints runtime for step and increments global runtime
    """
    final = time.time()
    interval = final - initfunctime
    global fullruntime
    fullruntime = fullruntime + interval
    print(text)
    print(interval)
    print('')
    start_functime()


def reset_runtime():
    """
    Resets for new analysis
    """
    global fullruntime
    fullruntime = 0
        
    
#UI WINDOWS#
class ftrdialog(QMainWindow): #dialog window for feature level data (spec, abund, hits)
    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = ui_featureinfo.Ui_MainWindow()
        self.ui.setupUi(self)
        
        self.sizegrip = QSizeGrip(self.ui.frame_grip_corner)
        self.ui.frame_grip_corner.setStyleSheet("background: transparent; background-image: url(:/resources/icons/24x24/resize.png); background-position: center; background-repeat: no-repeat;" )
        self.sizegrip.setToolTip("Resize Window")
        
        # MOVE WINDOW
        def moveWindow(event):
            # IF LEFT CLICK MOVE WINDOW
            if event.buttons() == Qt.LeftButton:
                self.move(self.pos() + event.globalPos() - self.dragPos2)
                self.dragPos2 = event.globalPos()
                event.accept()
        
        self.ui.label_title_bar_top.mouseMoveEvent = moveWindow 
        
    def mousePressEvent(self, event):
        self.dragPos2 = event.globalPos()
        
class dialog(QMainWindow): #plot configuration dialog
    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = ui_plotparam.Ui_MainWindow()
        self.ui.setupUi(self)
        
        def moveWindow(event):
            if event.buttons() == Qt.LeftButton:
                self.move(self.pos() + event.globalPos() - self.dragPos2)
                self.dragPos2 = event.globalPos()
                event.accept()
        
        self.ui.label_title_bar_top.mouseMoveEvent = moveWindow
        
    def mousePressEvent(self, event):
        self.dragPos2 = event.globalPos()
        
class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = Ui_MainWindow()
        
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.WindowStaysOnTopHint)  # set always on top flag, makes window disappear
        self.show() # makes window reappear, but it's ALWAYS on top
        self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowStaysOnTopHint) # clear always on top flag, makes window disappear
        #self.show() # makes window reappear, acts like normal window now (on top now but can be underneath if you raise another window)
        
        self.ui.setupUi(self)
        self.ui.label_credits.setText('v1.00 r24.11.17')
                
        #initialize other dialog windows
        self.dialog = dialog()
        self.ftrdialog = ftrdialog()
        self.ftrdialog.ui.btn_masst.hide()

        #set defaults
        self.recentdir = '..'
        self.savedata = {}
        self.pickedfeature = ''
        self.highlightcol = 'yellow'
        self.dbsearchdone = False
        self.groupsetmodel = GroupSetModel()  # the groupset MVC model
        self.groups = [] #biological groups in the dataset
        self.filename = Path.cwd()
        self.extractmetadatafilename = Path.cwd()
        self.samplelistfilename = Path.cwd()
        self.ui.label_status.setText('') #running/idle status label
        self.analysisrun = False
        self.ui.combo_maxisoshift.setCurrentIndex(3)
        self.outputdir = ''
        self.fragfilename = ''
        self.colourdic = {} #dictionary of available colours
        #self.atlas = pd.read_csv('npatlas.csv', sep = ',', header = [0], index_col = [1])
        #self.atlas = pd.read_csv('pyranodb.csv', sep = ',', header = [0], index_col = [1])

        

        
        # Per-plot widget/render state, one PlotSlot per plot name internally;
        # these six attributes are dict-like views over that single registry
        # (see plotslots.py) so existing self.fig['heatmap']-style call sites
        # throughout plotting.py and main.py work unchanged.
        self._plotslots = PlotSlotRegistry()
        self.fig = self._plotslots.fig
        self.canvas = self._plotslots.canvas
        self.pltlayout = self._plotslots.pltlayout
        self.toolbar = self._plotslots.toolbar
        self.ax = self._plotslots.ax
        self.highlight = self._plotslots.highlight
                
        self.ui.btn_run.clicked.connect(self.run_analysis)
        # Replaces the Designer-created QTreeWidget with a real QTreeView +
        # per-column filter bar (searchtree.py) -- doesn't edit ui_main.py,
        # see SearchTreePanel's docstring for how the runtime swap works.
        self.searchtree = SearchTreePanel(self.ui.treeWidget)
        self.searchtree.view.selectionModel().selectionChanged.connect(self.on_tree_item_selection_changed)


        def moveWindow(event):
            if UIFunctions.returnStatus() == 1:
                UIFunctions.maximize_restore(self)

            if event.buttons() == Qt.LeftButton:
                self.move(self.pos() + event.globalPos() - self.dragPos)
                self.dragPos = event.globalPos()
                event.accept()

        self.ui.label_title_bar_top.mouseMoveEvent = moveWindow
        UIFunctions.uiDefinitions(self)
        self.show()


    #---Methods---

        
    
    def exportgnps(self):
        """Export filtered GNPS peak table and optionally filter MGF file with new aligned IDs."""
        
        # Load the GNPS/MZmine file
        df_gnps = pd.read_csv(self.gnpsfilename, sep=',')
        
        # Find the RT and m/z columns
        rt_col = None
        mz_col = None
        
        for col in df_gnps.columns:
            if col == 'row m/z':
                mz_col = col
            elif col == 'row retention time':
                rt_col = col
        
        # Fallback to position-based if exact names not found
        if rt_col is None or mz_col is None:
            cols = df_gnps.columns.tolist()
            if 'row ID' in cols[0] or cols[0].lower() == 'row id':
                mz_col = cols[1] if mz_col is None else mz_col
                rt_col = cols[2] if rt_col is None else rt_col
        
        print(f"GNPS columns - RT: '{rt_col}', m/z: '{mz_col}'")
        
        # Load the filtered file 
        filtered_file = self.analysis_paramsgui.outputdir / (self.analysis_paramsgui.filename.stem + '_filtered.csv')
        df_filtered = pd.read_csv(filtered_file, header=[0,1,2])
        
        # Extract compound identifiers from filtered file
        compounds_to_keep = set()
        filtered_rt_mz = []
        
        first_col_values = df_filtered.iloc[:, 0]
        for val in first_col_values:
            if pd.notna(val) and '_' in str(val):
                compounds_to_keep.add(str(val))
                # Parse RT and m/z
                parts = str(val).split('_')
                if len(parts) == 2:
                    try:
                        rt = float(parts[0])
                        mz = float(parts[1])
                        filtered_rt_mz.append((rt, mz))
                    except:
                        pass
        
        print(f"Found {len(compounds_to_keep)} compounds in filtered file")
        
        # Tolerance-based matching for GNPS data
        rt_tolerance = 0.01
        mz_tolerance = 0.01
        
        matched_indices = []
        matched_compounds_info = []  # Store (index, rt, mz) for matched compounds
        
        for idx, row in df_gnps.iterrows():
            try:
                rt_gnps = float(row[rt_col])
                mz_gnps = float(row[mz_col])
                
                # Check if this matches any filtered compound
                for rt_filt, mz_filt in filtered_rt_mz:
                    if abs(rt_gnps - rt_filt) <= rt_tolerance and abs(mz_gnps - mz_filt) <= mz_tolerance:
                        matched_indices.append(idx)
                        matched_compounds_info.append((idx, rt_gnps, mz_gnps))
                        break
            except:
                continue
        
        # Create filtered dataframe with new sequential row IDs
        df_gnps_filtered = df_gnps.loc[matched_indices].copy()
        
        # Assign new sequential row IDs starting from 1
        if 'row ID' in df_gnps_filtered.columns:
            df_gnps_filtered['row ID'] = range(1, len(df_gnps_filtered) + 1)
        
        print(f"Filtered {len(df_gnps_filtered)} rows from {len(df_gnps)} total rows")
        
        # Save the filtered GNPS file
        output_file = self.analysis_paramsgui.outputdir / (self.gnpsfilename.stem + '_gnps_filtered.csv')
        df_gnps_filtered.to_csv(output_file, sep=',', index=False)
        print(f"Saved filtered GNPS file to: {output_file}")
        
        # ===== PROCESS MGF FILE IF IT EXISTS =====
        # Look for MGF file in the same directory as the GNPS file
        mgf_file = None
        gnps_dir = Path(self.gnpsfilename).parent
        
        # Search for .mgf file
        mgf_files = list(gnps_dir.glob('*.mgf'))
        if mgf_files:
            mgf_file = mgf_files[0]  # Take the first MGF file found
            print(f"\nFound MGF file: {mgf_file.name}")
        else:
            print("\nNo MGF file found in GNPS directory")
        
        if mgf_file:
            print("Processing MGF file...")
            
            # Create lookup for matched compounds
            rt_mz_to_new_id = {}
            for new_id, (orig_idx, rt, mz) in enumerate(matched_compounds_info, start=1):
                key = (round(rt, 5), round(mz, 5))
                rt_mz_to_new_id[key] = new_id
            
            # Parse and filter MGF
            matched_mgf_entries = []
            total_mgf_entries = 0
            
            def parse_mgf_entry(entry_lines):
                """Parse MGF entry to extract RT and m/z."""
                entry_data = {'lines': entry_lines}
                for line in entry_lines:
                    if line.startswith('FEATURE_ID='):
                        entry_data['feature_id'] = int(line.split('=')[1])
                    elif line.startswith('PEPMASS='):
                        entry_data['mz'] = float(line.split('=')[1])
                    elif line.startswith('RTINSECONDS='):
                        # Convert RT from seconds to minutes
                        entry_data['rt'] = float(line.split('=')[1]) / 60.0
                    elif line.startswith('RTINMINUTES='):
                        entry_data['rt'] = float(line.split('=')[1])
                return entry_data
            
            def update_mgf_feature_id(entry_lines, new_id):
                """Update FEATURE_ID in MGF entry."""
                updated_lines = []
                for line in entry_lines:
                    if line.startswith('FEATURE_ID='):
                        updated_lines.append(f'FEATURE_ID={new_id}\n')
                    else:
                        updated_lines.append(line)
                return updated_lines
            
            # Read and process MGF file
            with open(mgf_file, 'r') as f:
                current_entry = []
                
                for line in f:
                    if line.strip() == 'BEGIN IONS':
                        if current_entry:
                            total_mgf_entries += 1
                            entry_data = parse_mgf_entry(current_entry)
                            
                            if 'rt' in entry_data and 'mz' in entry_data:
                                rt = entry_data['rt']
                                mz = entry_data['mz']
                                
                                # Debug first few entries
                                if total_mgf_entries <= 3:
                                    print(f"  MGF entry {total_mgf_entries}: RT={rt:.4f} min, m/z={mz:.4f}")
                                
                                # Check if matches any filtered compound
                                matched = False
                                key = (round(rt, 5), round(mz, 5))
                                if key in rt_mz_to_new_id:
                                    new_id = rt_mz_to_new_id[key]
                                    entry_data['new_id'] = new_id
                                    matched_mgf_entries.append(entry_data)
                                    matched = True
                                else:
                                    # Check with tolerance
                                    for (comp_rt, comp_mz), new_id in rt_mz_to_new_id.items():
                                        if abs(comp_rt - rt) <= rt_tolerance and abs(comp_mz - mz) <= mz_tolerance:
                                            entry_data['new_id'] = new_id
                                            matched_mgf_entries.append(entry_data)
                                            matched = True
                                            if len(matched_mgf_entries) <= 3:
                                                print(f"    Matched to compound RT={comp_rt:.4f}, m/z={comp_mz:.4f}")
                                            break
                            
                            current_entry = []
                        current_entry.append(line)
                    elif line.strip() == 'END IONS':
                        current_entry.append(line)
                    elif current_entry:
                        current_entry.append(line)
                
                # Process last entry
                if current_entry:
                    total_mgf_entries += 1
                    entry_data = parse_mgf_entry(current_entry)
                    if 'rt' in entry_data and 'mz' in entry_data:
                        rt = entry_data['rt']
                        mz = entry_data['mz']
                        key = (round(rt, 5), round(mz, 5))
                        if key in rt_mz_to_new_id:
                            new_id = rt_mz_to_new_id[key]
                            entry_data['new_id'] = new_id
                            matched_mgf_entries.append(entry_data)
            
            # Sort by new ID and write filtered MGF
            matched_mgf_entries.sort(key=lambda x: x['new_id'])
            
            mgf_output = self.analysis_paramsgui.outputdir / (mgf_file.stem + '_filtered.mgf')
            with open(mgf_output, 'w') as f:
                for entry_data in matched_mgf_entries:
                    new_id = entry_data['new_id']
                    updated_lines = update_mgf_feature_id(entry_data['lines'], new_id)
                    for line in updated_lines:
                        f.write(line.rstrip('\n') + '\n')
                    f.write('\n')
            
            print(f"Filtered MGF: {len(matched_mgf_entries)} entries from {total_mgf_entries} total")
            print(f"Saved filtered MGF to: {mgf_output}")
        
        # Update status
        status_text = f'Filtered GNPS table ({len(df_gnps_filtered)} compounds)'
        if mgf_file:
            status_text += f' and MGF ({len(matched_mgf_entries)} entries)'
        self.ui.label_status.setText(status_text)
        
        print("\n✓ Export completed successfully!")
                
    def error(self, message):
        self.ui.label_status.setText(message)
        self.ui.label_status.setStyleSheet('color: rgb(150,0,0);')
    
    def getgroups(self):
        """
        Get biological groups on input of all input files, fills comboboxes with these.
        """
        try:
            groups, unresolved = compute_biological_groups(
                self.extractmetadatafilename, self.samplelistfilename, self.filename)
        except ValueError as exc:
            self.error(str(exc))
            return

        if unresolved:
            # Report once, not once per row -- a dataset with many injections
            # missing from the metadata join would otherwise spam the console
            # with no visible indication anything is wrong, and self.groups
            # would silently end up incomplete with no signal to the combo
            # boxes/groupset editor that consume it.
            unique_unresolved = sorted(set(map(str, unresolved)))
            print(f"getgroups: {len(unique_unresolved)} injection(s) had no "
                  f"matching Biological_Group and were skipped: {unique_unresolved}")
            self.error(f"{len(unique_unresolved)} injection(s) missing from "
                       "metadata -- biological group list may be incomplete.")

        self.groups = groups
    
        # Set experimental and control grp defaults in ui
        self.ui.combo_blankfil_name.clear()
        self.ui.combo_blankfil_name.addItems(self.groups)
        self.ui.combo_expgrp.clear()
        self.ui.combo_expgrp.addItems(self.groups)
        self.ui.combo_ctrgrp.clear()
        self.ui.combo_ctrgrp.addItems(self.groups)
        
        # Set experimental and control group defaults in ui
        if hasattr(self, 'analysis_paramsgui') and self.analysis_paramsgui.statstgrps:
            self.ui.combo_expgrp.setCurrentText(str(self.analysis_paramsgui.statstgrps[0]))
            self.ui.combo_ctrgrp.setCurrentText(str(self.analysis_paramsgui.statstgrps[1]))
        else:
            pos = 0
            expset = False
            for group in self.groups:
                if 'blank' in group.lower():
                    self.ui.combo_blankfil_name.setCurrentIndex(pos)
                else:
                    if expset:
                        self.ui.combo_ctrgrp.setCurrentIndex(pos)
                    else:
                        self.ui.combo_expgrp.setCurrentIndex(pos)
                        expset = True
                pos += 1
    
        
    def fulldbsearch(self):
        """
        Run a full compound database search. Filter the database matches within a specified mass window.
        Concatenate the hits and sort them by parts-per-million (ppm).
        """
        self.hitdb, _ = search_npatlas(
            self.analysis_paramsgui.outputdir,
            self.analysis_paramsgui.filename.stem,
            self.atlas,
            self.analysis_paramsgui.ppmthresh,
        )

        


    def fillfttree(self):
        # Fill feature tree with database hits
        iondict = cached_read_csv(
            self.analysis_paramsgui.outputdir / 'iondict.csv',
            sep=',',
            header=[0],
            index_col=None
        )
        iondict = iondict[iondict['hits'] >= 0]

        # Rows in searchtree.COLUMNS order (Compound, m/z, TR, Max, Sets,
        # Groups, FC, Hits) -- raw values, not pre-formatted strings; the
        # model formats numeric columns for display and keeps the raw value
        # available for proper numeric sort/filter (see IonTableModel).
        rows = [
            (row['Compound'], row['m/z'], row['Retention time (min)'], row['max'],
             row['groups'], row['groups'], row['fc'], row['hits'])
            for _, row in iondict.iterrows()
        ]
        self.searchtree.set_rows(rows)



    def on_tree_item_selection_changed(self):
        name = self.searchtree.selected_compound()
        if name:
            # Call the method to highlight the feature
            self.highlight_feature(name)

            # Call the onItemClicked logic
            #self.onItemClicked(item)
                                     
    def runsearch(self, mass): # refactor if I can save a database of hits
        #possibly make a third common method used in both dbsearch and this method
        """
        Runs a search for a mass and displays the database hits in the feature tree dialog.
        
        Args:
            mass (float): The mass to search for.
        
        Returns:
            None.
        """
        ppmwindow = self.analysis_paramsgui.ppmthresh
        hits_h = self.atlas[abs(1000000*(self.atlas['compound_m_plus_h'] - mass)/self.atlas['compound_m_plus_h']) < ppmwindow]
        hits_h['ppm'] = abs(1000000*(hits_h['compound_m_plus_h'] - mass)/hits_h['compound_m_plus_h'])
        hits_na = self.atlas[abs(1000000*(self.atlas['compound_m_plus_na'] - mass)/self.atlas['compound_m_plus_na']) < ppmwindow]
        hits_na['ppm'] = abs(1000000*(hits_na['compound_m_plus_na'] - mass)/hits_na['compound_m_plus_na'])
        hits = pd.concat([hits_h, hits_na])
        hits = hits.sort_values(by=['ppm'])
        
        x=0 #rename this variable
        itemdict = {}
        smilesdict = {}
        self.ftrdialog.ui.treeWidget.clear()
        for index, row in hits.iterrows():
            itemdict[x] = QtWidgets.QTreeWidgetItem( [row['compound_name'], row['compound_molecular_formula'], str(row['compound_accurate_mass']), str(row['ppm']), (str(row['genus'] + ' ' + row['origin_species'])), row['compound_smiles']])
            self.ftrdialog.ui.treeWidget.addTopLevelItem(itemdict[x])
            x = x + 1
        
        def onItemClicked(): #show structure of selected match
        
            def clean_ascii(text):
                printable = set(string.printable)
                ascii_chars = filter(lambda x: x in printable, text)
                return ''.join(ascii_chars)
        
                
            if len(self.ftrdialog.ui.treeWidget.selectedItems()) > 0:
                item = self.ftrdialog.ui.treeWidget.selectedItems()[0]
                cmpd =  clean_ascii(item.text(0))
                cmpd = re.sub(r'\/<>.*{}\|', ' ', cmpd)

                if os.path.isfile('compoundimages/' + (cmpd + '.png')):
                    pixmap = QPixmap('compoundimages/' + (cmpd + '.png'))
                else:
                    m = indigo.loadMolecule(item.text(5))
                    indigo.setOption('render-image-size', '400,400')
                    renderer.renderToFile(m, 'compoundimages/' + (cmpd + '.png'))
                    pixmap = QPixmap('compoundimages/' + (cmpd + '.png'))
                
                pixmap2 = pixmap.scaled(400, 400, QtCore.Qt.KeepAspectRatio)
                self.ftrdialog.ui.label.setPixmap(pixmap2)
        
        self.ftrdialog.ui.treeWidget.itemSelectionChanged.connect(onItemClicked)
        if len(hits.index) == 0:
            pixmap = QPixmap('blank.png')
            pixmap2 = pixmap.scaled(400, 400, QtCore.Qt.KeepAspectRatio)
            self.ftrdialog.ui.label.setPixmap(pixmap2)
        else:
            self.ftrdialog.ui.treeWidget.setCurrentItem(itemdict[0])
                    
    def highlight_feature(self, newfeature):
        """
        Select (or, if the same feature is clicked twice, deselect) a feature
        and highlight it across all plots.

        This is the entry point for actual selection events (clicking a plot
        point, a tree row, or arrow-key heatmap navigation) -- it is where the
        "click the same feature twice to clear it" toggle belongs. To simply
        refresh the displays for the *currently* selected feature (e.g. when
        switching tabs in the feature-info dialog) call ``_refresh_highlight``
        instead; calling this method with ``newfeature == self.pickedfeature``
        would incorrectly toggle the highlight off.

        Args:
            newfeature (str): The name of the new feature to highlight.

        Returns:
            None.
        """
        # Deselect the highlighted feature if clicked twice
        if newfeature == self.pickedfeature and self.highlightcol != (0, 0, 0, 0):
            self.highlightcol = (0, 0, 0, 0)
        else:
            self.highlightcol = 'yellow'

        self.pickedfeature = newfeature
        self._refresh_highlight()

    def _refresh_highlight(self):
        """Redraw every plot/dialog for the current ``self.pickedfeature`` and
        ``self.highlightcol`` without changing selection state.

        Safe to call repeatedly (e.g. on every feature-info tab switch) since,
        unlike ``highlight_feature``, it never toggles ``highlightcol``. Also
        guards against the feature-info plot objects (``abundplt``/``spec``)
        not existing yet, and against the current feature not being present in
        the loaded fragmentation database.
        """
        for plot in self.highlight:
            self.highlight[plot].set_color(self.highlightcol)

        if not self.pickedfeature:
            return

        # Read iondict file to get ion data
        iondict = cached_read_csv(self.analysis_paramsgui.outputdir / 'iondict.csv',
                              sep=',', header=[0], index_col=[0])

        # Update volcano plot with the selected feature
        if self.analysis_paramsgui.Volcanoplt:
            # 'logfc' is never persisted to iondict.csv -- plot_volcano
            # computes it in-memory from 'fc' (np.log2(iondict['fc'])) each
            # time it draws, rather than writing it back to disk. Derive it
            # the same way here instead of indexing a column that doesn't
            # exist in the file.
            self.highlight['volcano'].set_data(
                [np.log2(iondict.loc[self.pickedfeature, 'fc'])],
                [iondict.loc[self.pickedfeature, '-logq']]
            )
            self.canvas['volcano'].draw_idle()

            # Update MZRT plot with the selected feature
            self.highlight['mzrt'].set_data(
                [iondict.loc[self.pickedfeature, 'Retention time (min)']],
                [iondict.loc[self.pickedfeature, 'm/z']]
            )
            self.canvas['mzrt'].draw_idle()

        # Update KMD plot with the selected feature
        if self.analysis_paramsgui.KMD:
            self.highlight['kmd'].set_data(
                [iondict.loc[self.pickedfeature, 'm/z']],
                [iondict.loc[self.pickedfeature, 'kmd']]
            )
            self.canvas['kmd'].draw_idle()

        # Update feature plot with the selected feature
        self.highlight['featureplt'].set_data(
            [iondict.loc[self.pickedfeature, 'Retention time (min)']],
            [iondict.loc[self.pickedfeature, 'm/z']]
        )
        self.canvas['featureplt'].draw_idle()

        msdata = pd.read_csv(
            self.analysis_paramsgui.outputdir / (self.analysis_paramsgui.filename.stem + '_filtered.csv'),
            sep=',', header=[2], index_col=[0]
        )

        # Set heatmap highlight based on view
        self.heatind = self.cmind.index(msdata.index.to_list().index(self.pickedfeature))
        xlim = int(self.canvas['heatmap'].figure.axes[2].get_xlim()[1])
        self.highlight['heatmap'].set_data(
            [0, xlim, xlim, 0, 0],
            [self.heatind, self.heatind, self.heatind+1, self.heatind+1, self.heatind]
        )
        self.canvas['heatmap'].draw_idle()

        # Run search when feature is selected
        if self.ftrdialog.ui.stackedWidget.currentIndex() == 0 and not self.ftrdialog.isHidden():
            self.runsearch(iondict.loc[self.pickedfeature, 'm/z'])

        # Reset spec if a fragmentation database with this feature is loaded.
        # spec may not exist yet (e.g. a prior dataset had no frag file) and
        # the feature may not be present in a newly-loaded fragdb, so both are
        # guarded rather than assumed.
        fragdb = getattr(self, 'fragdb', 'None')
        if self.fragfilename != '' and getattr(self, 'spec', None) is not None \
                and fragdb != 'None' and self.pickedfeature in getattr(fragdb, 'ions', {}):
            self.safe_generate('spectrum plot', self.spec.reset, fragdb.ions[self.pickedfeature].pattern)

        if self.ftrdialog.ui.stackedWidget.currentIndex() == 2 and not self.ftrdialog.isHidden() \
                and getattr(self, 'abundplt', None) is not None:
            self.safe_generate('abundance plot', self.abundplt.reset)


     # Move heatmap selection up or down with W/S key press
    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_W:
            self.mv_heatmap(-1)
        elif event.key() == QtCore.Qt.Key_S:
            self.mv_heatmap(1)
    
    def mv_heatmap(self, shift):
        """
        Moves heatmap selection up or down by a given shift value.
        
        Args:
            shift (int): The value to shift the heatmap selection by.
        
        Returns:
            None.
        """
        if self.ui.stackedWidget_plot.currentIndex() != 8 or self.ui.stackedWidget.currentIndex() != 3:
            return

        cmind = getattr(self, 'cmind', None)
        heatind = getattr(self, 'heatind', None)
        if not cmind or heatind is None:
            return

        # Clamp rather than wrap/raise: Python's negative indexing would
        # silently jump to the opposite end of the list on underflow, and an
        # unclamped index overflows with IndexError on the last row -- both
        # triggerable just by holding W/S past either end of the heatmap.
        new_heatind = max(0, min(len(cmind) - 1, heatind + shift))
        if new_heatind == heatind:
            return

        iondict = cached_read_csv(self.analysis_paramsgui.outputdir / 'iondict.csv', sep=',', header=[0], index_col=[0])
        msdata = cached_read_csv(self.analysis_paramsgui.outputdir / (self.analysis_paramsgui.filename.stem + '_filtered.csv'), sep=',', header=[2], index_col=[0]).iloc[:, 2:]
        index = cmind[new_heatind]
        name = msdata.index.tolist()[index]
    
        self.ui.lbl_featurename.setText('Compound: ' + name)
        self.ui.lbl_featurert.setText('Retention time: ' + str(iondict.loc[name, 'Retention time (min)']))
        self.ui.lbl_featuremz.setText('m/z: ' + str(iondict.loc[name, 'm/z']))
        self.highlight_feature(name)
     
    
    def read_save(self, savefile):
        """Reads a saved analysis session from a pickle file.

        Args:
            savefile (Path): The path to the saved session file.
        """
        # Loading a session may point at an entirely different dataset's
        # output directory -- drop any cached reads from whatever was
        # previously loaded/analysed.
        invalidate_csv_cache()

        print("read1")
        with open(savefile, 'rb') as read_pi:
            print("read3")
            self.savedata = pickle.load(read_pi)
            print("read2")

        # Assign the analysis parameters
        self.analysis_paramsgui = self.savedata['parameters']
            
        # Set output dir to save file directory, make rawdata folder
        self.analysis_paramsgui.outputdir = Path(savefile).parent
        Path(self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'rawdata').mkdir(parents=True, exist_ok=True)
            
        # Write raw files from save file
        self.savedata['peaktable'].to_csv(self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'rawdata' / self.analysis_paramsgui.filename.name, header=False, index=False)
        self.savedata['samplelist'].to_csv(self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'rawdata' / self.analysis_paramsgui.samplelistfilename.name, header=False, index=False)
        self.savedata['extractmetadata'].to_csv(self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'rawdata' / self.analysis_paramsgui.extractmetadatafilename.name, header=False, index=False)
        
        # Write frag file if it exists
        if self.analysis_paramsgui.fragfilename:
            fragmsp = open(self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'rawdata' / self.analysis_paramsgui.fragfilename.name, 'w')     
            fragmsp.write(self.savedata['fragdb'])
            fragmsp.close()
            
        # Assign file names and output directory
        self.outputdir = self.analysis_paramsgui.outputdir
        self.samplelistfilename = self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'rawdata' / self.analysis_paramsgui.samplelistfilename.name
        self.extractmetadatafilename = self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'rawdata' / self.analysis_paramsgui.extractmetadatafilename.name
        
        # Assign frag file name if it exists
        if self.analysis_paramsgui.fragfilename:
            self.fragfilename = self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'rawdata' / Path(self.analysis_paramsgui.fragfilename).name
        
        self.filename = self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'rawdata' / self.analysis_paramsgui.filename.name
    
        # Get groups and load groupsets (accepts both old pickled `query`
        # objects and current `GroupSet` objects -- see GroupSet.from_legacy).
        self.getgroups()
        self.groupsetmodel = GroupSetModel.from_legacy_list(self.analysis_paramsgui.queries)
        UIFunctions.updatesets(self)
        
        
    def write_save(self):
        """Write the current analysis session to a ``.mpct`` (pickle) file.

        Each component is captured defensively, so that one unreadable input
        (e.g. a file held open by another program) does not abort the whole
        save. The pickle is written to a temporary file that atomically
        replaces the target, so an interrupted save cannot corrupt an existing
        ``.mpct``.
        """
        params = self.analysis_paramsgui
        params.queries = self.groupsetmodel.to_legacy_list()  # plain GroupSet list

        def _capture(label, reader, default=None):
            try:
                return reader()
            except Exception:
                import traceback
                print('Warning: could not capture ' + label + ' for save:')
                traceback.print_exc()
                return default

        savedata = dict(self.savedata)
        savedata['parameters'] = params
        savedata['peaktable'] = _capture(
            'peak table',
            lambda: pd.read_csv(params.outputdir / params.filename.name, sep=',', header=None, index_col=None, low_memory=False))
        savedata['samplelist'] = _capture(
            'sample list',
            lambda: pd.read_csv(self.samplelistfilename, sep=',', header=None, index_col=None))
        savedata['extractmetadata'] = _capture(
            'extract metadata',
            lambda: pd.read_csv(self.extractmetadatafilename, sep=',', header=None, index_col=None))
        if params.fragfilename != '':
            savedata['fragdb'] = _capture(
                'fragment database',
                lambda: Path(params.fragfilename).read_text(), default='None')
        else:
            savedata['fragdb'] = 'None'
        self.savedata = savedata

        target = params.outputdir / (params.filename.stem + '.mpct')
        tmp = target.with_name(target.name + '.tmp')
        try:
            with open(tmp, 'wb') as file_pi:
                pickle.dump(savedata, file_pi)
            os.replace(tmp, target)  # atomic replace on the same filesystem
        except Exception:
            import traceback
            print('Error saving .mpct file:')
            traceback.print_exc()
            self.error('Analysis complete but save failed (see console)')
            try:
                if tmp.exists():
                    tmp.unlink()
            except Exception:
                pass
    
    def export_filtered_outputs(self):
        """Write post-filter exports next to the analysis output (best-effort).

        - The filtered peak table in its source layout (a row subset of the
          input file; matched by compound id, falling back to m/z + RT).
        - If a fragmentation file is loaded, an MSP/MGF re-indexed so its entry
          IDs match the filtered peak-table row order -- the form GNPS2 expects.

        Failures are logged and never abort the analysis.
        """
        import translators
        params = self.analysis_paramsgui
        filtered_internal = params.outputdir / (params.filename.stem + '_filtered.csv')

        try:
            source = Path(params.filename)
            out = params.outputdir / (params.filename.stem + '_filtered_source' + source.suffix)
            count = translators.filter_source_peaktable(source, filtered_internal, out)
            print('Wrote filtered source peak table (' + str(count) + ' features): ' + str(out))
        except Exception:
            import traceback
            print('Could not export filtered source peak table:')
            traceback.print_exc()

        fragfile = getattr(params, 'fragfilename', '')
        if fragfile:
            try:
                fragpath = Path(fragfile)
                out = params.outputdir / (fragpath.stem + '_reindexed' + fragpath.suffix)
                count = translators.reindex_fragments(filtered_internal, fragpath, out)
                print('Wrote GNPS2 re-indexed fragments (' + str(count) + ' entries): ' + str(out))
            except Exception:
                import traceback
                print('Could not re-index fragmentation file:')
                traceback.print_exc()

    def safe_generate(self, label, func, *args, **kwargs):
        """Run a plot/generation step so a single failure is non-fatal.

        Logs the traceback and notes the failure in the status bar, then
        returns None so the remaining plots are still generated.
        """
        try:
            return func(*args, **kwargs)
        except Exception:
            import traceback
            print('Error generating ' + label + ':')
            traceback.print_exc()
            self.ui.label_status.setText('Warning: ' + label + ' failed (see console)')
            return None

    def _create_or_reset(self, attr, label, create_call, reset_call=None):
        """Create a plot object the first time it's needed, or reset() it if
        it already exists.

        Plot (re)generation used to be gated entirely on ``self.analysisrun``
        (a whole-session "has any analysis completed" flag): the very first
        analysis created every enabled plot, and every later run/regenerate
        only ever called ``.reset()`` on them. That breaks as soon as an
        optional output toggles between runs in the same session -- e.g. the
        first dataset has no fragmentation file (so ``self.spec`` is never
        created) and a later dataset does: ``self.spec.reset(...)`` would then
        raise ``AttributeError`` because the object was never created. Keying
        the create-vs-reset decision off whether the attribute already exists
        (rather than off analysisrun) fixes this for every optional plot.
        """
        if getattr(self, attr, None) is None:
            result = self.safe_generate(label, create_call)
            if result is not None:
                setattr(self, attr, result)
        elif reset_call is not None:
            self.safe_generate(label, reset_call)

    def _generate_plots(self):
        """Create or refresh every plot for the current analysis parameters.

        Shared by the post-compute step (``_finish_analysis``, after a fresh
        ``run_MSFaST``) and the dialog's "Apply" button (``regenerateplts``),
        so both paths get the same create-if-missing safety net.
        """
        params = self.analysis_paramsgui
        pltfile = params.outputdir / (params.filename.stem + '_filtered.csv')
        iondictfile = params.outputdir / 'iondict.csv'
        dfs = self.filtereddfs
        grpsts = self.groupsets

        if params.CVfil:
            self._create_or_reset('prevcv', 'CV plot',
                lambda: prev_cv(self, 'cvplt', self.ui.frame_cvplt, 'none', 'none', 'none'),
                lambda: self.prevcv.reset(self, '', ''))
            stop_functime('cvplt complete')

        self._create_or_reset('ftplt', 'feature plot',
            lambda: show_featureplt(self, 'featureplt', self.ui.frame_featureplt, iondictfile, '', ''),
            lambda: self.ftplt.reset(iondictfile, '', ''))
        stop_functime('ftplt complete')

        self._create_or_reset('dend', 'dendrogram',
            lambda: plot_dendrogram(self, 'dend', self.ui.frame_dend, pltfile, '', ''),
            lambda: self.dend.reset(pltfile, '', ''))
        stop_functime('dendrogram complete')

        if params.PCA:
            self._create_or_reset('pca', 'PCA/NMDS plot',
                lambda: plot_PCA(self, 'pca', self.ui.frame_pca, pltfile, '', ''),
                lambda: self.pca.reset(pltfile, '', ''))
            stop_functime('nmds complete')

        if params.FC3Dplt:
            self._create_or_reset('fc3d', '3D fold-change plot',
                lambda: plot_fc3d(self, 'fc3d', self.ui.frame_fc3d, iondictfile, dfs, grpsts),
                lambda: self.fc3d.reset(iondictfile, dfs, grpsts))
            stop_functime('fc3d complete')

        if params.KMD:
            self._create_or_reset('kmd', 'Kendrick mass defect plot',
                lambda: kendrick(self, 'kmd', self.ui.frame_kmd, iondictfile, dfs, grpsts),
                lambda: self.kmd.reset(iondictfile, dfs, grpsts))
            stop_functime('kmd complete')

        if params.MZRTplt:
            self._create_or_reset('mzrt', 'm/z vs RT plot',
                lambda: plot_mzrt(self, 'mzrt', self.ui.frame_mzrt, iondictfile, dfs, grpsts),
                lambda: self.mzrt.reset(iondictfile, dfs, grpsts))
            stop_functime('mzrt complete')

        if params.Volcanoplt:
            self._create_or_reset('volcano', 'volcano plot',
                lambda: plot_volcano(self, 'volcano', self.ui.frame_volcano, iondictfile, dfs, grpsts),
                lambda: self.volcano.reset(iondictfile, dfs, grpsts))
            stop_functime('volcano complete')

        self._create_or_reset('heatmap', 'heatmap',
            lambda: plot_heatmap(self, 'heatmap', self.ui.frame_heatmap, pltfile),
            lambda: self.heatmap.reset(self, 'heatmap', self.ui.frame_heatmap, pltfile))
        stop_functime('heatmap complete')

        # abundplt/spec are only ever (re)drawn for the currently picked
        # feature via _refresh_highlight, not here -- just make sure they
        # exist when needed.
        self._create_or_reset('abundplt', 'abundance plot', lambda: plot_abund(self, 'abund'))
        if self.fragfilename != '':
            self._create_or_reset('spec', 'spectrum plot', lambda: show_spectrum(self, 'spec'))

        self._create_or_reset('samplecorr', 'sample correlation plot',
            lambda: plot_samplecorr(self, 'samplecorr', self.ui.frame_samplecorr, iondictfile, dfs, grpsts),
            lambda: self.samplecorr.reset(iondictfile, dfs, grpsts))
        stop_functime('samplecorr complete')

    def run_analysis(self):
        # Ignore re-clicks while an analysis is already running on the worker thread.
        if getattr(self, '_analysis_thread', None) is not None and self._analysis_thread.isRunning():
            return
        self.dbsearchdone = False
        # iondict.csv/_filtered.csv/_summarydata.csv from any previous run
        # are about to be superseded -- drop cached reads of them so feature
        # selection etc. doesn't serve stale data from before this run.
        invalidate_csv_cache()
        start_functime()
        self.enumerate_inputs()
        print('')
        stop_functime('inputs obtained')
        
        # Heavy, Qt-free computation runs on a worker thread so the GUI stays
        # responsive; plotting resumes on the main thread in _finish_analysis
        # once the worker emits 'finished'.
        self.ui.btn_run.setEnabled(False)
        self.ui.label_status.setText('Processing data...')
        self.ui.label_status.setStyleSheet('color: rgb(150,150,150);')

        self._analysis_thread = QtCore.QThread()
        self._analysis_worker = AnalysisWorker(self)
        self._analysis_worker.moveToThread(self._analysis_thread)
        self._analysis_thread.started.connect(self._analysis_worker.run)
        self._analysis_worker.finished.connect(self._on_compute_finished)
        self._analysis_worker.failed.connect(self._on_compute_failed)
        self._analysis_worker.finished.connect(self._analysis_thread.quit)
        self._analysis_worker.failed.connect(self._analysis_thread.quit)
        # The thread/worker are kept as attributes (replaced on the next run)
        # rather than deleteLater'd, so the re-entrancy guard can safely call
        # isRunning() afterwards without touching a deleted C++ object.
        self._analysis_thread.start()

    def _on_compute_failed(self, msg):
        stop_functime('calculations failed')
        self.error(msg)
        self.ui.btn_run.setEnabled(True)

    def _on_compute_finished(self):
        print('')
        stop_functime('calculations complete')
        try:
            self._finish_analysis()
        except Exception:
            import traceback
            traceback.print_exc()
            self.error('Error during plot generation (see console)')
        finally:
            self.ui.btn_run.setEnabled(True)

    def _finish_analysis(self):
        try:
            gen_treemap(self)  # move back to end
        except Exception:
            print("not generating tremap due to an error")
        stop_functime('treemap complete')
        
        # Used for point opacity based on abundance colouring
        iondict = cached_read_csv(self.analysis_paramsgui.outputdir / 'iondict.csv', sep=',', header=[0], index_col=None)
        self.analysis_paramsgui.maxval = iondict['logmax'].max()

        self._generate_plots()
        self.analysisrun = True
        
        #text = open(self.analysis_paramsgui.outputdir / 'analysisinfo.txt').read() #writes output text to report tab
        #self.ui.textBrowser_info.setPlainText(text)
    
    

        # Writes filtering statistics in data review summary
        msdata_unformatted = pd.read_csv(self.analysis_paramsgui.filename, sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
        msdata_formatted = pd.read_csv(self.analysis_paramsgui.outputdir / (self.analysis_paramsgui.filename.stem + '_formatted.csv'), sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
        msdata_filtered = cached_read_csv(self.analysis_paramsgui.outputdir / (self.analysis_paramsgui.filename.stem + '_filtered.csv'), sep=',', header=[0, 1, 2], index_col=[0, 1, 2])
        #test = open('testionfilters.pkl', 'wb') # this exports ionfilters as a pickle for debuging
        #pickle.dump(self.ionfilters, test)
        
        text = ''
        if self.analysis_paramsgui.relfil:
            text += 'Features failing peak correction filtering: ' + str(len(self.ionfilters['relfil'].ions)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(self.ionfilters['relfil'].ions) / len(msdata_unformatted.index), 2)) + '%\n'
        if self.analysis_paramsgui.blnkfltr:
            text += 'Features failing blank filtering: ' + str(len(self.groupionlists[self.analysis_paramsgui.blnkgrp])) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(self.groupionlists[self.analysis_paramsgui.blnkgrp]) / len(msdata_unformatted.index), 2)) + '%\n'
        if self.analysis_paramsgui.decon:
            text += 'Features in-source ion filtering: ' + str(len(self.ionfilters['insource'].ions)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(self.ionfilters['insource'].ions) / len(msdata_unformatted.index), 2)) + '%\n'
        if self.analysis_paramsgui.CVfil:
            text += 'Features failing CV filtering: ' + str(len(self.ionfilters['cv'].ions)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(self.ionfilters['cv'].ions) / len(msdata_unformatted.index), 2)) + '%\n'
        text += 'Features failing any filters: ' + str(len(msdata_unformatted.index) - len(msdata_filtered.index)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * (len(msdata_unformatted.index) - len(msdata_filtered.index)) / len(msdata_unformatted.index), 2)) + '%\n'
        text += 'Features passing all filters: ' + str(len(msdata_filtered.index)) + '/' + str(len(msdata_unformatted.index)) + ' ' + str(round(100 * len(msdata_filtered.index) / len(msdata_unformatted.index), 2)) + '%\n'
        
        self.ui.textBrowser_mp_prev.setPlainText(text)
        
        # Imports msp fragment database
        if self.fragfilename != '':
            self.fragdb = getfragdb.importfrag(self.analysis_paramsgui.fragfilename)
            self.ftrdialog.ui.btn_spectrum.show()
        else:
            self.fragdb = 'None'
            self.ftrdialog.ui.btn_spectrum.hide()
        
        # Runs DB search only if tab is active
        if self.ui.stackedWidget.currentIndex() == 5:
            self.fulldbsearch()
            self.fillfttree()
            self.dbsearchdone = True
        
        try:
            gen_upsetplt(self)
        except Exception:
            print("not generating upset plot due to an error")
        stop_functime('upsetplt complete')
        self.ui.label_status.setText('Analysis Complete')
        stop_functime('analysis complete')
        print('')
        print('full runtime:' + str(fullruntime))
        print('')
        print('')
        reset_runtime()


        self.write_save()
        self.export_filtered_outputs()
        
    def enumerate_inputs(self):
        self.analysis_paramsgui = analysis_parameters()
        if self.ui.checkBox_blankfilter.isChecked():
            self.analysis_paramsgui.blnkgrp = str(self.ui.combo_blankfil_name.currentText())
        else:
            self.analysis_paramsgui.blnkgrp = ''
        self.analysis_paramsgui.blnkfltr = self.ui.checkBox_blankfilter.isChecked()

        if len(self.groupsetmodel) == 0:
            if self.analysis_paramsgui.blnkfltr:
                UIFunctions.addgroup(self, 'Features not in blanks')
                src = list(self.groupsetmodel.selected.excl)
                src.remove(self.analysis_paramsgui.blnkgrp)
                self.groupsetmodel.update(self.groupsetmodel.selected_index,
                    src=src, excl=[self.analysis_paramsgui.blnkgrp])
            else:
                UIFunctions.addgroup(self, 'All features')
                self.groupsetmodel.update(self.groupsetmodel.selected_index,
                    src=self.groupsetmodel.selected.excl, excl=[])
            UIFunctions.updategroups(self)
        else:
            self.groupsetmodel.select(self.ui.listWidget_pltgrps.currentRow())
            UIFunctions.writegroups(self)
            UIFunctions.updategroups(self)
        
        self.analysis_paramsgui.kingdom = self.ui.combo_kingdom.currentText()
        self.analysis_paramsgui.genus = str(self.ui.lineEdit_genus.text())
        self.atlas = pd.read_csv('npatlas.tsv', sep='\t', header=0, index_col=1)
        if len(self.analysis_paramsgui.kingdom) > 3:
            self.atlas = self.atlas[self.atlas['origin_type'] == self.analysis_paramsgui.kingdom]
        if len(self.analysis_paramsgui.genus) > 3:
            self.atlas = self.atlas[self.atlas['genus'] == self.analysis_paramsgui.genus]

        
        # A list of active filter-name tokens, not a hand-built string --
        # consumers that need the old space-joined string shape (e.g. an
        # older .mpct save) go through groupsets.normalize_graphfilters().
        self.analysis_paramsgui.graphfilters = []
        if self.ui.checkBox_cv.isChecked():
            self.analysis_paramsgui.graphfilters.append('cv')
        if self.ui.checkBox_mp.isChecked():
            self.analysis_paramsgui.graphfilters.append('rel')
        if self.ui.checkBox_decon.isChecked():
            self.analysis_paramsgui.graphfilters.append('insource')
            self.analysis_paramsgui.deconthresh = float(self.ui.lineEdit_insourcethresh.text())
        
        

        
        # Build the {descriptive_name: GroupSet} mapping MSFaST/plotting
        # consume, merging in the active graph-level filters (cv/rel/insource).
        querydict = build_query_dict(self.groupsetmodel, self.analysis_paramsgui.graphfilters)
        self.analysis_paramsgui.querydict = querydict
        self.analysis_paramsgui.querylist = list(querydict.keys())
        
        # Set analysis parameters
        self.analysis_paramsgui.filename = self.filename
        self.analysis_paramsgui.samplelistfilename = self.samplelistfilename
        self.analysis_paramsgui.extractmetadatafilename = self.extractmetadatafilename 
        if self.outputdir == '':
            self.outputdir = (Path(self.analysis_paramsgui.filename).parent)
        self.analysis_paramsgui.outputdir = self.outputdir
        self.analysis_paramsgui.fragfilename = self.fragfilename
        Path(self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem).mkdir(parents=True, exist_ok=True)
        Path(self.analysis_paramsgui.outputdir / self.analysis_paramsgui.filename.stem / 'plots').mkdir(parents=True, exist_ok=True)
        self.analysis_paramsgui.outputdir /= self.analysis_paramsgui.filename.stem
    
        self.analysis_paramsgui.statstgrps = [str(self.ui.combo_expgrp.currentText()), str(self.ui.combo_ctrgrp.currentText())]
        self.analysis_paramsgui.cvthresh = float(self.ui.lineEdit_cvthresh.text())
        if self.ui.radioButton_meancv.isChecked():
            self.analysis_paramsgui.cvparam = 'average CV'
        else:
            self.analysis_paramsgui.cvparam = 'median CV'
            
        
        # Plain 1:1 checkbox fields (no branching/derived values) -- see
        # paramfields.py for the shared save/restore schema; the matching
        # restore side is loadsession()'s restore_checkbox_fields() call.
        save_checkbox_fields(self, self.analysis_paramsgui)

        self.analysis_paramsgui.pqthresh = float(self.dialog.ui.lineEdit_pqthresh.text())
        self.analysis_paramsgui.fcthresh = float(self.dialog.ui.lineEdit_fcthresh.text())
        self.analysis_paramsgui.colorscheme = self.dialog.ui.combo_colorscheme.currentText()
        self.analysis_paramsgui.relfil = self.ui.checkBox_mp.isChecked()
        self.analysis_paramsgui.CVfil = self.ui.checkBox_cv.isChecked()
        self.analysis_paramsgui.decon = self.ui.checkBox_decon.isChecked()
        self.analysis_paramsgui.merge = self.ui.checkBox_merge.isChecked()
        self.analysis_paramsgui.RTwin = float(self.ui.lineEdit_rtwin.text())
        self.analysis_paramsgui.ringingwin = float(self.ui.lineEdit_ringwin.text())
        self.analysis_paramsgui.isopeakwin = float(self.ui.lineEdit_isowin.text())
        self.analysis_paramsgui.dimerpeakwin = float(self.ui.lineEdit_isowin.text())
        self.analysis_paramsgui.maxisowin = float(self.ui.combo_maxisoshift.currentText())
        self.analysis_paramsgui.grpave = True
        self.analysis_paramsgui.prperr = True
        self.analysis_paramsgui.blnkfltr = self.ui.checkBox_blankfilter.isChecked()
        
        if self.ui.radioButton_blankfilter_abs.isChecked():
            self.analysis_paramsgui.blankfilparam = 'absolute'
            self.analysis_paramsgui.blankfilthresh = float(self.ui.lineEdit_blankfilter_absthresh.text())
        else:
            self.analysis_paramsgui.blankfilparam = 'relative'
            self.analysis_paramsgui.blankfilthresh = float(self.ui.lineEdit_blankfilter_relthresh.text())
        if self.ui.checkBox_blankfilter.isChecked() and self.analysis_paramsgui.blnkgrp == '':
            self.analysis_paramsgui.blnkgrp = 'Blanks'
        self.analysis_paramsgui.ppmthresh = float(self.ui.lineEdit_ppmthresh.text())
    
    
    def regenerateplts(self):
        """Re-render plots in place (e.g. from the dialog's "Apply" button)
        without re-running the analysis pipeline.

        Delegates to the same create-or-reset logic used after a fresh
        analysis, so a plot that doesn't exist yet (e.g. because the prior
        analysis had it disabled) is created rather than crashing on
        ``.reset()``.
        """
        self._generate_plots()
        self.ui.label_status.setText('Update Complete')

        

    ## APP EVENTS
    ########################################################################
    def mousePressEvent(self, event):
        self.dragPos = event.globalPos()
        
    
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    if sys.platform != 'win32':
        app.setStyle('Fusion')
    app.setStyleSheet("QFrame { border: 0px; }") #QToolTip { color: #999999; background-color: rgb(0, 255, 0); border: 1px solid grey; }")
    window = MainWindow()
    sys.exit(app.exec_())
