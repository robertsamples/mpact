"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
"""

## ==> GUI FILE
from main import MainWindow, start_functime, stop_functime, reset_runtime, ftrdialog, dialog
#import masstdriver #from old version of masst search push
import webbrowser #may not be needed now
from mzmineimport import format_check
from paramfields import restore_checkbox_fields

import sys
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import time

import platform
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QSizeGrip, QGraphicsDropShadowEffect, QFileDialog, QListWidgetItem, QColorDialog
from PyQt5.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint, QRect, QSize, QTime, QUrl, Qt, QEvent)
from PyQt5.QtGui import QBrush, QColor, QIcon, QPalette, QPainter, QPixmap
from pathlib import Path


# IMPORT FUNCTIONS AND RESOURCES
import files

from MSFaST import analysis_parameters

import os




## ==> GLOBALS

GLOBAL_STATE = 0

class UIFunctions(MainWindow):

    ## ==> MAXIMIZE RESTORE FUNCTION
    def maximize_restore(self):
        global GLOBAL_STATE
        status = GLOBAL_STATE

        # IF NOT MAXIMIZED
        if status == 0:
            self.showMaximized()

            # SET GLOBAL TO 1
            GLOBAL_STATE = 1

            # IF MAXIMIZED REMOVE MARGINS AND BORDER RADIUS
            #self.ui.drop_shadow_layout.setContentsMargins(0, 0, 0, 0)
            #self.ui.drop_shadow_frame.setStyleSheet("background-color: qlineargradient(spread:pad, x1:0, y1:1, x2:0, y2:0, stop:0 rgba(34, 34, 34, 255), stop:1 rgba(52, 52, 52, 255));")
            self.ui.btn_maximize.setToolTip("Restore")
        else:
            GLOBAL_STATE = 0
            self.showNormal()
            self.resize(self.width()+1, self.height()+1)
            #self.ui.drop_shadow_layout.setContentsMargins(10, 10, 10, 10) #not sure about this one throws error when run, may cut
            self.ui.drop_shadow_frame.setStyleSheet("background-color: qlineargradient(spread:pad, x1:0, y1:1, x2:0, y2:0, stop:0 rgba(34, 34, 34, 255), stop:1 rgba(52, 52, 52, 255));")
            self.ui.btn_maximize.setToolTip("Maximize")
            

    ## ==> UI DEFINITIONS
    def uiDefinitions(self):

        # REMOVE TITLE BAR
        self.setWindowFlag(QtCore.Qt.FramelessWindowHint)

        # SET DROPSHADOW WINDOW
        self.shadow = QGraphicsDropShadowEffect(self)
        self.shadow.setBlurRadius(20)
        self.shadow.setXOffset(0)
        self.shadow.setYOffset(0)
        self.shadow.setColor(QColor(0, 0, 0, 100))

        # APPLY DROPSHADOW TO FRAME
        #self.ui.drop_shadow_frame.setGraphicsEffect(self.shadow)

        self.ui.frame_plts.hide()
        self.ui.checkBox_fc.hide()
        self.ui.checkBox_ttest.hide()
        self.dialog.ui.checkBox_bootstrap.setChecked(True)

        # Top bar functions
        self.ui.btn_maximize.clicked.connect(lambda: UIFunctions.maximize_restore(self))
        self.ui.btn_minimize.clicked.connect(lambda: self.showMinimized())
        self.ui.btn_close.clicked.connect(lambda: self.close())

        #mainbar functions
        self.ui.stackedWidget.setCurrentIndex(6)
        self.ui.btn_import.clicked.connect(lambda: UIFunctions.goto_import(self))
        self.ui.btn_filter.clicked.connect(lambda: UIFunctions.goto_filter(self))
        self.ui.btn_parameters.clicked.connect(lambda: UIFunctions.goto_params(self))
        self.ui.btn_plots.clicked.connect(lambda: UIFunctions.goto_plot(self))
        self.ui.btn_info.clicked.connect(lambda: UIFunctions.goto_info(self))
        self.ui.btn_search.clicked.connect(lambda: UIFunctions.goto_search(self))

        #plotbar functions
        self.ui.stackedWidget_plot.setCurrentIndex(0)
        self.ui.btn_review.clicked.connect(lambda: UIFunctions.goto_review(self))
        self.ui.btn_upset.clicked.connect(lambda: UIFunctions.goto_upset(self))
        self.ui.btn_dend.clicked.connect(lambda: UIFunctions.goto_dend(self))
        self.ui.btn_pca.clicked.connect(lambda: UIFunctions.goto_pca(self))
        self.ui.btn_mzrt.clicked.connect(lambda: UIFunctions.goto_mzrt(self))
        self.ui.btn_kmd.clicked.connect(lambda: UIFunctions.goto_kmd(self))
        self.ui.btn_3dfc.clicked.connect(lambda: UIFunctions.goto_3dfc(self))
        self.ui.btn_volcano.clicked.connect(lambda: UIFunctions.goto_volcano(self))
        self.ui.btn_heatmap.clicked.connect(lambda: UIFunctions.goto_heatmap(self))
        
        self.ui.stackedWidget_review.setCurrentIndex(3)
        self.ui.btn_ftrplt.clicked.connect(lambda: self.ui.stackedWidget_review.setCurrentIndex(0))
        self.ui.btn_treemap.clicked.connect(lambda: self.ui.stackedWidget_review.setCurrentIndex(1))
        self.ui.btn_cvplt.clicked.connect(lambda: self.ui.stackedWidget_review.setCurrentIndex(2))
        self.ui.btn_datasummary.clicked.connect(lambda: self.ui.stackedWidget_review.setCurrentIndex(3))
        
        self.ui.btn_upsetplt.clicked.connect(lambda: self.ui.stackedWidget_grpanalysis.setCurrentIndex(0))
        self.ui.btn_samplecorr.clicked.connect(lambda: self.ui.stackedWidget_grpanalysis.setCurrentIndex(1))

        #feature info bar functions
        self.ftrdialog.ui.btn_close.clicked.connect(lambda: self.ftrdialog.hide())  
        self.ftrdialog.ui.btn_hits.clicked.connect(lambda: UIFunctions.goto_hits(self)) 
        self.ftrdialog.ui.btn_abund.clicked.connect(lambda: UIFunctions.goto_abund(self))
        self.ftrdialog.ui.btn_spectrum.clicked.connect(lambda: UIFunctions.goto_spectrum(self)) 
        self.ftrdialog.ui.btn_masst.clicked.connect(lambda: UIFunctions.masst(self)) 
        self.ui.btn_details.clicked.connect(lambda: UIFunctions.show_ftrdialog(self)) 
        self.ui.btn_details_2.clicked.connect(lambda: UIFunctions.show_ftrdialog(self)) 
        
        #input functions
        self.ui.btn_load_session.clicked.connect(lambda: UIFunctions.loadsession(self))
        self.ui.btn_import_pktbl.clicked.connect(lambda: UIFunctions.getfilename(self))
        self.ui.btn_import_spllist.clicked.connect(lambda: UIFunctions.getsamplelistfilename(self))
        self.ui.btn_import_splmdt.clicked.connect(lambda: UIFunctions.getextractmetadatafilename(self))
        self.ui.btn_import_outdir.clicked.connect(lambda: UIFunctions.getoutputdir(self))
        self.ui.btn_import_msp.clicked.connect(lambda: UIFunctions.getfragfilename(self))
        
        self.ui.btn_importgnps.clicked.connect(lambda: UIFunctions.getgnpsfilename(self))
        self.ui.btn_exportgnps.clicked.connect(lambda: UIFunctions.exportgnps(self))

        
        #colour picker functions
        self.ui.btn_col1.clicked.connect(lambda: UIFunctions.colour_picker1(self))

        self.ui.btn_addgroup.clicked.connect(lambda: UIFunctions.addgroup(self))
        self.ui.btn_removegroup.clicked.connect(lambda: UIFunctions.removegroup(self))
        self.ui.listWidget_pltgrps.currentItemChanged.connect(lambda: UIFunctions.writegroups(self))

        #dialog functions and parameters
        self.dialog.setWindowFlag(QtCore.Qt.FramelessWindowHint)
        self.dialog.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint)
        self.dialog.ui.btn_close.clicked.connect(lambda: self.dialog.hide())
        self.dialog.ui.btn_apply.clicked.connect(self.regenerateplts)
        
        self.ftrdialog.setWindowFlag(QtCore.Qt.FramelessWindowHint)
        self.ftrdialog.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint)
        self.ftrdialog.ui.treeWidget.hideColumn(5)
        self.searchtree.view.hideColumn(5) #hide sets column for now
        self.ftrdialog.hide()
        self.ftrdialog.ui.btn_related.hide()
        
        
        ## ==> CREATE SIZE GRIP TO RESIZE WINDOW
        self.sizegrip = QSizeGrip(self.ui.frame_grip_corner)
        self.ui.frame_grip_corner.setStyleSheet("background: transparent; background-image: url(:/resources/icons/24x24/resize.png); background-position: center; background-repeat: no-repeat;" )
        self.sizegrip.setToolTip("Resize Window")


    ## RETURN STATUS IF WINDOWS IS MAXIMIZE OR RESTAURED
    def returnStatus():
        return GLOBAL_STATE


    #ui buttons
    def goto_import(self):
        self.ui.stackedWidget.setCurrentIndex(0)
        UIFunctions.reset_mainbar(self)
        self.ui.btn_import.setStyleSheet(self.ui.mainbar_activebtn)
        
    def goto_filter(self):
        self.ui.stackedWidget.setCurrentIndex(1)
        UIFunctions.reset_mainbar(self)
        self.ui.btn_filter.setStyleSheet(self.ui.mainbar_activebtn)       
        
    def goto_params(self):
        self.ui.stackedWidget.setCurrentIndex(2)
        UIFunctions.reset_mainbar(self)
        self.ui.btn_parameters.setStyleSheet(self.ui.mainbar_activebtn)

    def goto_plot(self):
        self.ui.stackedWidget.setCurrentIndex(3)
        UIFunctions.reset_mainbar(self)
        self.ui.btn_plots.setStyleSheet(self.ui.mainbar_activebtn)

    def goto_info(self):
        self.ui.stackedWidget.setCurrentIndex(4)
        UIFunctions.reset_mainbar(self)
        self.ui.btn_info.setStyleSheet(self.ui.mainbar_activebtn)
        
    def goto_search(self):
        self.ui.stackedWidget.setCurrentIndex(5)
        UIFunctions.reset_mainbar(self)
        self.ui.btn_search.setStyleSheet(self.ui.mainbar_activebtn)
        #self.ftrdialog.show()
        if self.dbsearchdone == False and self.analysisrun:
            start_functime()
            self.fulldbsearch()
            self.fillfttree()
            self.dbsearchdone = True
            stop_functime('dbsearch complete')
            reset_runtime()
            
    #plotbar functions
    def goto_review(self):
        self.ui.stackedWidget_infobar.setCurrentIndex(0)
        self.ui.stackedWidget_plot.setCurrentIndex(1)
        UIFunctions.reset_plotbar(self)
        self.ui.btn_review.setStyleSheet(self.ui.plotbar_activebtn)
        
        self.dialog.ui.checkBox_applyfilter.hide()

    def goto_upset(self):
        self.ui.stackedWidget_infobar.setCurrentIndex(1)
        self.ui.stackedWidget_plot.setCurrentIndex(9)
        UIFunctions.reset_plotbar(self)
        self.ui.btn_upset.setStyleSheet(self.ui.plotbar_activebtn)
        
        self.dialog.ui.checkBox_applyfilter.hide()
        self.dialog.ui.frame_colorscheme.show()
        
    def goto_dend(self):
        self.ui.stackedWidget_infobar.setCurrentIndex(1)
        self.ui.stackedWidget_plot.setCurrentIndex(2)
        UIFunctions.reset_plotbar(self)
        self.ui.btn_dend.setStyleSheet(self.ui.plotbar_activebtn)
        
    def goto_pca(self):
        self.ui.stackedWidget_infobar.setCurrentIndex(1)
        self.ui.stackedWidget_plot.setCurrentIndex(3)
        UIFunctions.reset_plotbar(self)
        self.ui.btn_pca.setStyleSheet(self.ui.plotbar_activebtn)
        
        self.dialog.ui.frame_2.show()

    def goto_mzrt(self):
        self.ui.stackedWidget_infobar.setCurrentIndex(0)
        self.ui.stackedWidget_plot.setCurrentIndex(4)
        UIFunctions.reset_plotbar(self)
        self.ui.btn_mzrt.setStyleSheet(self.ui.plotbar_activebtn)
        
    def goto_kmd(self):
        self.ui.stackedWidget_infobar.setCurrentIndex(0)
        self.ui.stackedWidget_plot.setCurrentIndex(5)
        UIFunctions.reset_plotbar(self)
        self.ui.btn_kmd.setStyleSheet(self.ui.plotbar_activebtn)
        
        self.dialog.ui.frame_mdguide.show()
        
    def goto_3dfc(self):
        self.ui.stackedWidget_infobar.setCurrentIndex(0)
        self.ui.stackedWidget_plot.setCurrentIndex(6)
        UIFunctions.reset_plotbar(self)
        self.ui.btn_3dfc.setStyleSheet(self.ui.plotbar_activebtn)
        
    def goto_volcano(self):
        self.ui.stackedWidget_infobar.setCurrentIndex(0)
        self.ui.stackedWidget_plot.setCurrentIndex(7)
        UIFunctions.reset_plotbar(self)
        self.ui.btn_volcano.setStyleSheet(self.ui.plotbar_activebtn)
        
        self.dialog.ui.frame_volcanoparams.show()

    def goto_heatmap(self):
        self.ui.stackedWidget_infobar.setCurrentIndex(0)
        self.ui.stackedWidget_plot.setCurrentIndex(8)
        UIFunctions.reset_plotbar(self)
        self.ui.btn_heatmap.setStyleSheet(self.ui.plotbar_activebtn)

        self.dialog.ui.frame_colorscheme.show()
        
    #ftinfobar functions
    def goto_abund(self):
        self.ftrdialog.ui.btn_masst.hide()
        self.ftrdialog.ui.stackedWidget.setCurrentIndex(2)
        # Refresh (not re-select/toggle) the current feature's display now
        # that the abundance tab is active.
        self._refresh_highlight()
        UIFunctions.reset_ftrdialogbar(self)
        self.ftrdialog.ui.btn_abund.setStyleSheet(self.ui.ftbar_activebtn)

    def goto_hits(self):
        self.ftrdialog.ui.btn_masst.hide()
        self.ftrdialog.ui.stackedWidget.setCurrentIndex(0)

        self._refresh_highlight()
        UIFunctions.reset_ftrdialogbar(self)
        self.ftrdialog.ui.btn_hits.setStyleSheet(self.ui.ftbar_activebtn)
        
    def goto_spectrum(self):
        self.ftrdialog.ui.btn_masst.show()
        self.ftrdialog.ui.stackedWidget.setCurrentIndex(1)
        UIFunctions.reset_ftrdialogbar(self)
        self.ftrdialog.ui.btn_spectrum.setStyleSheet(self.ui.ftbar_activebtn)

    def masst(self): 
        """
        Launches the web browser and searches for mass spectra of selected feature in Mass Spectrometry Interactive Virtual Environment (MassIVE).
        
        Args:
            None
        
        Returns:
            None
        """
        if self.fragdb != 'None' and self.fragdb.ions[self.pickedfeature].pattern.shape[0] > 0:
            description = str("MPACT_SUBMISSION:" + self.analysis_paramsgui.filename.stem + "_" + self.pickedfeature)
            precursor = str(self.fragdb.ions[self.pickedfeature].fragparams['PrecursorMZ']).strip()
            fragments = ''
            for row in self.fragdb.ions[self.pickedfeature].pattern:
                fragments = fragments + str(row[0]) + r'\t' + str(row[1]) + r'\n'
            url = r'''https://gnps.ucsd.edu/ProteoSAFe/index.jsp#{%22workflow%22:%22SEARCH_SINGLE_SPECTRUM%22,%22precursor_mz%22:"'''
            url = url + precursor + r'''%22,%22spectrum_string%22:"'''
            url = url + fragments + r'''%22,%22desc%22:"'''
            url = url + description + r'''%22}'''
            #print(url)
            webbrowser.open(url)

            #used for external version with nonlogin masst            
            #masstdriver.push(description, precursor, fragments)

    #colour buttons
    def colour_picker1(self):
        """
        Allows user to pick a color from color dialog and sets the selected
        groupset's colour attribute.
        """
        color = QColorDialog.getColor()
        self.groupsetmodel.update(self.groupsetmodel.selected_index, colour=color.name())
        self.ui.btn_col1.setStyleSheet("QPushButton {border: 2px solid lightgrey; background-color: " + str(self.groupsetmodel.selected.colour) +";}")

    def updatesets(self):
        """
        Rebuilds the groupset list widget from the model and restores selection.
        """
        list_widget = self.ui.listWidget_pltgrps

        # Block signals while rebuilding: clear()/addItem()/setCurrentRow()
        # below all fire currentItemChanged, which is wired to writegroups()
        # (ui_functions.py:148) -- a GUI-to-model sync that would otherwise
        # run mid-rebuild, before the drag-list boxes are refreshed for the
        # new selection, and clobber the model with stale leftover box
        # contents from whatever groupset was selected a moment earlier.
        list_widget.blockSignals(True)
        list_widget.clear()
        for groupset in self.groupsetmodel:
            item = QListWidgetItem(groupset.name)
            item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable)
            list_widget.addItem(item)

        # Trust the model's own selected_index (already correctly maintained
        # by add()/remove()/from_legacy_list()) rather than re-deriving a
        # guess from the widget's currentRow() before it's cleared above --
        # that previously read the stale pre-rebuild row and shifted it by
        # one, so live selection drifted toward row 0 on every remove/load.
        list_widget.setCurrentRow(self.groupsetmodel.selected_index)
        list_widget.blockSignals(False)

        UIFunctions.updategroups(self)


    def updategroups(self):
        """
        Syncs the three group-membership lists (and colour swatch) in the GUI
        to the currently-selected groupset in the model.
        """
        selgroup = self.groupsetmodel.select(self.ui.listWidget_pltgrps.currentRow())

        self.ui.listWidget_allgrps.clear()
        self.ui.listWidget_orgrps.clear()
        self.ui.listWidget_andgrps.clear()

        groupset = self.groupsetmodel.selected
        if groupset is None:
            return

        for group in groupset.excl:
            self.ui.listWidget_allgrps.addItem(QListWidgetItem(group))
        for group in groupset.src:
            self.ui.listWidget_orgrps.addItem(QListWidgetItem(group))
        for group in groupset.incl:
            self.ui.listWidget_andgrps.addItem(QListWidgetItem(group))

        self.ui.btn_col1.setStyleSheet(
            "QPushButton {border: 2px solid lightgrey; background-color: " +
            str(groupset.colour) +";}"
        )

    def writegroups(self):
        """
        Writes the GUI's groupset name/included/excluded/source-group widgets
        back into the model for the currently-selected groupset.
        """
        selset = self.groupsetmodel.selected_index
        if 0 <= selset < self.ui.listWidget_pltgrps.count():
            try:
                name = self.ui.listWidget_pltgrps.item(selset).text()
            except Exception:
                name = None

            self.groupsetmodel.update(
                selset,
                name=name,
                excl=[self.ui.listWidget_allgrps.item(x).text() for x in range(self.ui.listWidget_allgrps.count())],
                src=[self.ui.listWidget_orgrps.item(x).text() for x in range(self.ui.listWidget_orgrps.count())],
                incl=[self.ui.listWidget_andgrps.item(x).text() for x in range(self.ui.listWidget_andgrps.count())],
            )

        UIFunctions.updategroups(self)


    def addgroup(self, name = 'New Feature Set'):
        """
        Adds a new groupset (excluding every known biological group by
        default) to the model and the GUI list.
        """
        item = QListWidgetItem(name)
        item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable)

        list_widget = self.ui.listWidget_pltgrps

        # Block signals for the same reason as updatesets(): setCurrentItem()
        # below fires currentItemChanged -> writegroups(), which would write
        # whatever's currently in the (stale, belonging to the previously
        # selected groupset) drag-list boxes into this brand-new groupset's
        # slot instead of leaving its real (exclude-everything) default data
        # in place. Update the model first, then sync the GUI from it
        # ourselves via updategroups(), rather than relying on the signal.
        list_widget.blockSignals(True)
        list_widget.addItem(item)
        self.groupsetmodel.add(item.text(), all_groups=self.groups)
        list_widget.setCurrentItem(item)
        list_widget.blockSignals(False)

        UIFunctions.updategroups(self)


    def removegroup(self):
        """
        Removes the selected groupset from the model and the GUI list.
        """
        self.groupsetmodel.remove(self.ui.listWidget_pltgrps.currentRow())
        UIFunctions.updatesets(self)
        
        
        
    #import buttons
    def loadsession(self):
        """Open a saved ``.mpct`` session and restore every input/parameter widget.

        Opening the file is fatal if it fails, but each parameter is then
        restored independently so a single missing/renamed field (e.g. from an
        older save file) cannot silently abort restoration of the rest. The
        previous implementation wrapped the whole body in ``try/except: pass``,
        so the first failing widget (notably ``combo_maxisoshift``, which was
        being given a float index that raises ``TypeError`` in PyQt5) left every
        later parameter un-restored.
        """
        savefile, _ = QFileDialog.getOpenFileName(self, 'Open file', self.recentdir, "*.mpct")
        if not savefile:
            return
        self.savefile = Path(savefile)
        self.recentdir = str(self.savefile.parent)
        try:
            self.read_save(self.savefile)
        except Exception:
            import traceback
            traceback.print_exc()
            self.error('Could not open session file (see console)')
            return

        p = self.analysis_paramsgui

        def restore(what, action):
            try:
                action()
            except Exception as exc:
                print('Warning: could not restore ' + what + ': ' + str(exc))

        def set_maxiso():
            # maxisowin is stored as a float (e.g. 3.0); map it to the matching
            # combo item by text rather than passing a float as an index.
            combo = self.ui.combo_maxisoshift
            for cand in (str(int(p.maxisowin)), str(p.maxisowin)):
                idx = combo.findText(cand)
                if idx >= 0:
                    combo.setCurrentIndex(idx)
                    return

        def set_blankfil():
            if 'absolute' in p.blankfilparam:
                self.ui.lineEdit_blankfilter_absthresh.setText(str(p.blankfilthresh))
                self.ui.radioButton_blankfilter_abs.setChecked(True)
            else:
                self.ui.lineEdit_blankfilter_relthresh.setText(str(p.blankfilthresh))
                self.ui.radioButton_blankfilter_rel.setChecked(True)

        # ----- input file path labels -----
        restore('peak table label', lambda: self.ui.lbl_pktbl.setText(self.filename.name))
        restore('sample list label', lambda: self.ui.lbl_spllist.setText(self.samplelistfilename.name))
        restore('metadata label', lambda: self.ui.lbl_splmdt.setText(self.extractmetadatafilename.name))
        restore('fragment file label', lambda: self.ui.lbl_msp.setText(Path(self.fragfilename).name if self.fragfilename else ''))
        restore('output dir label', lambda: self.ui.lbl_outdir.setText(
            str(self.outputdir) if len(str(self.outputdir)) < 40 else '...' + str(self.outputdir)[-40:]))

        # ----- filtering parameters -----
        restore('CV filter', lambda: self.ui.checkBox_cv.setChecked('cv' in p.graphfilters))
        restore('mean CV radio', lambda: self.ui.radioButton_meancv.setChecked('average' in p.cvparam))
        restore('median CV radio', lambda: self.ui.radioButton_medcv.setChecked('median' in p.cvparam))
        restore('CV threshold', lambda: self.ui.lineEdit_cvthresh.setText(str(p.cvthresh)))
        restore('mispick filter', lambda: self.ui.checkBox_mp.setChecked('rel' in p.graphfilters))
        restore('merge ions', lambda: self.ui.checkBox_merge.setChecked(getattr(p, 'merge', False)))
        restore('RT window', lambda: self.ui.lineEdit_rtwin.setText(str(p.RTwin)))
        restore('ringing window', lambda: self.ui.lineEdit_ringwin.setText(str(getattr(p, 'ringingwin', ''))))
        restore('isotope window', lambda: self.ui.lineEdit_isowin.setText(str(p.isopeakwin)))
        restore('max isotope shift', set_maxiso)
        restore('deconvolution filter', lambda: self.ui.checkBox_decon.setChecked('insource' in p.graphfilters))
        restore('deconvolution threshold', lambda: self.ui.lineEdit_insourcethresh.setText(str(getattr(p, 'deconthresh', ''))))

        # ----- blank filter -----
        restore('blank filter', lambda: self.ui.checkBox_blankfilter.setChecked(p.blnkfltr))
        restore('blank group', lambda: self.ui.combo_blankfil_name.setCurrentText(str(p.blnkgrp)))
        restore('blank filter mode', set_blankfil)

        # ----- database search -----
        restore('kingdom', lambda: self.ui.combo_kingdom.setCurrentText(str(getattr(p, 'kingdom', ''))))
        restore('genus', lambda: self.ui.lineEdit_genus.setText(str(getattr(p, 'genus', ''))))
        restore('ppm threshold', lambda: self.ui.lineEdit_ppmthresh.setText(str(getattr(p, 'ppmthresh', ''))))

        # ----- plot toggles -----
        # Plain 1:1 checkbox fields (no branching/derived values) -- see
        # paramfields.py for the shared save/restore schema; the matching
        # save side is enumerate_inputs()'s save_checkbox_fields() call.
        restore_checkbox_fields(self, p, restore)

        # ----- statistics / plot params -----
        restore('experimental group', lambda: self.ui.combo_expgrp.setCurrentText(str(p.statstgrps[0])))
        restore('control group', lambda: self.ui.combo_ctrgrp.setCurrentText(str(p.statstgrps[1])))
        restore('colour scheme', lambda: self.dialog.ui.combo_colorscheme.setCurrentText(p.colorscheme))
        restore('p/q threshold', lambda: self.dialog.ui.lineEdit_pqthresh.setText(str(p.pqthresh)))
        restore('fold-change threshold', lambda: self.dialog.ui.lineEdit_fcthresh.setText(str(p.fcthresh)))
                
    def getfilename(self):
            self.filename, _  = QFileDialog.getOpenFileName(self, 'Open file', self.recentdir ,
                                                                           "*.csv *.txt")
            self.filename = Path(self.filename)
            self.recentdir = str(self.filename.parent)
            self.ui.lbl_pktbl.setText(self.filename.name)
            format_check(self)
            if self.filename.suffix == '.csv' and self.extractmetadatafilename.suffix == '.csv' and self.samplelistfilename.suffix == '.csv':
                self.getgroups() 

    def getsamplelistfilename(self):
            self.samplelistfilename, _ = QFileDialog.getOpenFileName(self, 'Open file', self.recentdir ,
                                                                           "*.csv")
            self.samplelistfilename = Path(self.samplelistfilename)
            self.recentdir = str(self.samplelistfilename)
            self.ui.lbl_spllist.setText(self.samplelistfilename.name) 
            if self.filename.suffix == '.csv' and self.extractmetadatafilename.suffix == '.csv' and self.samplelistfilename.suffix == '.csv':
                self.getgroups() 

    def getextractmetadatafilename(self):
            self.extractmetadatafilename, _ = QFileDialog.getOpenFileName(self, 'Open file', self.recentdir ,
                                                                           "*.csv")
            self.extractmetadatafilename = Path(self.extractmetadatafilename)
            self.recentdir = str(self.extractmetadatafilename)
            self.ui.lbl_splmdt.setText(self.extractmetadatafilename.name)
            if self.filename.suffix == '.csv' and self.extractmetadatafilename.suffix == '.csv' and self.samplelistfilename.suffix == '.csv':
                self.getgroups() 
                
    def getoutputdir(self):
            self.outputdir = Path(str(QFileDialog.getExistingDirectory(self, 'Select Directory', self.recentdir)))
            self.recentdir = str(self.outputdir)
            if len(str(self.outputdir)) <40:
                try:
                    self.ui.lbl_outdir.setText(self.outputdir)
                except Exception:
                    self.error('No directory selected')
                    pass
                    return()
            else:
                self.ui.lbl_outdir.setText('...' + str(self.outputdir)[-40:])

    def getfragfilename(self):
            self.fragfilename, _ = QFileDialog.getOpenFileName(self, 'Open file', self.recentdir,
                                                                           "*.msp *.mgf")
            self.recentdir = str(self.outputdir)
            self.fragfilename = Path(self.fragfilename)
            self.ui.lbl_msp.setText(self.fragfilename.name)

    def getgnpsfilename(self):
            self.gnpsfilename, _ = QFileDialog.getOpenFileName(self, 'Open file', self.recentdir,
                                                                           "*.csv *.txt")
            self.gnpsfilename = Path(self.gnpsfilename)
            self.ui.lbl_gnpstable.setText(self.gnpsfilename.name)



    #reset uibar buttons
    def reset_mainbar(self):
        self.ui.mainbar_activebtn = "QPushButton {	border: none; background-color: transparent;}QPushButton:hover {background-color: transparent}QPushButton:pressed {	background-color: rgb(75, 75, 75);}"
        self.ui.mainbar_inactivebtn ="QPushButton {	border: none; background-color: rgb(35, 35, 35);}QPushButton:hover {background-color: transparent}QPushButton:pressed {	background-color: rgb(75, 75, 75);}"
        self.ui.btn_import.setStyleSheet(self.ui.mainbar_inactivebtn)
        self.ui.btn_filter.setStyleSheet(self.ui.mainbar_inactivebtn)
        self.ui.btn_parameters.setStyleSheet(self.ui.mainbar_inactivebtn)
        self.ui.btn_plots.setStyleSheet(self.ui.mainbar_inactivebtn)
        self.ui.btn_info.setStyleSheet(self.ui.mainbar_inactivebtn)
        self.ui.btn_search.setStyleSheet(self.ui.mainbar_inactivebtn)
        self.ui.label_status.setText('') #to reset analysis status after tabs are switched. eventually use a second thread to dynamically update status
        self.ui.label_status.setStyleSheet('color: rgb(150,150,150);')
        
    def reset_plotbar(self):
        self.ui.plotbar_activebtn = "QPushButton {	border: none;background-color: rgba(225,225,225, 255);	color: rgba(75,75,75,255)} QPushButton:hover {background-color: rgba(225,225,225, 255);	color: rgba(75,75,75,255)} QPushButton:pressed {	 background-color: rgba(225,225,225, 255);	color: rgba(75,75,75,255)}"
        self.ui.plotbar_inactivebtn ="QPushButton {	border: none;background-color: rgba(0,0,0, 0)} QPushButton:hover {background-color: rgba(225,225,225, 50)} QPushButton:pressed {	 background-color: rgba(225,225,225, 255);	color: rgba(75,75,75,255)}"
        self.ui.btn_review.setStyleSheet(self.ui.plotbar_inactivebtn)
        self.ui.btn_dend.setStyleSheet(self.ui.plotbar_inactivebtn)
        self.ui.btn_upset.setStyleSheet(self.ui.plotbar_inactivebtn)
        self.ui.btn_pca.setStyleSheet(self.ui.plotbar_inactivebtn)
        self.ui.btn_mzrt.setStyleSheet(self.ui.plotbar_inactivebtn)
        self.ui.btn_kmd.setStyleSheet(self.ui.plotbar_inactivebtn)
        self.ui.btn_3dfc.setStyleSheet(self.ui.plotbar_inactivebtn)
        self.ui.btn_volcano.setStyleSheet(self.ui.plotbar_inactivebtn)
        self.ui.btn_heatmap.setStyleSheet(self.ui.plotbar_inactivebtn)
        
        self.dialog.ui.frame_2.hide()
        self.dialog.ui.frame.hide() #hide apply button, button doesnt work right now, use run button
        self.dialog.ui.frame_mdguide.hide()
        self.dialog.ui.frame_volcanoparams.hide()
        self.dialog.ui.frame_colorscheme.hide()
        self.dialog.ui.checkBox_applyfilter.show()
        self.ui.label_status.setText('')

    def reset_ftrdialogbar(self):
        self.ui.ftbar_activebtn = "QPushButton {border: none;background-color: rgba(225,225,225, 255);	color: rgba(75,75,75,255)} QPushButton:hover {background-color: rgba(225,225,225, 255);	color: rgba(75,75,75,255)} QPushButton:pressed {	 background-color: rgba(225,225,225, 255);	color: rgba(75,75,75,255)}"
        self.ui.ftbar_inactivebtn ="QPushButton {border: none;background-color: rgba(0,0,0,0)} QPushButton:hover {background-color: rgba(225,225,225, 50)} QPushButton:pressed {	 background-color: rgba(225,225,225, 255);	color: rgba(75,75,75,255)}"
        self.ftrdialog.ui.btn_hits.setStyleSheet(self.ui.ftbar_inactivebtn)
        self.ftrdialog.ui.btn_spectrum.setStyleSheet(self.ui.ftbar_inactivebtn)
        self.ftrdialog.ui.btn_abund.setStyleSheet(self.ui.ftbar_inactivebtn)
    
    def show_ftrdialog(self):
        self.ftrdialog.show()
        if self.pickedfeature != '':
            self._refresh_highlight()

