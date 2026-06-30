"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas
"""

## ==> GUI FILE
from main import MainWindow, start_functime, stop_functime, reset_runtime
import webbrowser  # used by masst() to open the GNPS single-spectrum search
from mzmineimport import format_check
from paramfields import restore_checkbox_fields

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

from PyQt5 import QtCore
from PyQt5.QtWidgets import QSizeGrip, QGraphicsDropShadowEffect, QFileDialog, QListWidgetItem, QColorDialog
from PyQt5.QtGui import QColor
from pathlib import Path


# Imported for its side effect only (registers the Qt resource blobs used by
# the stylesheets/icons); do not remove even though it looks unused.
import files




## ==> GLOBALS

GLOBAL_STATE = 0

# Tab-switch tables. Each nav button maps to the stacked-widget page(s) it
# selects (plus, for plot tabs, any dialog-frame show/hide deltas to apply
# after the bar reset). goto_maintab()/goto_plottab() are the single
# data-driven dispatchers that replace the ~15 near-identical goto_* methods;
# the buttons are wired to them by key in uiDefinitions(). Module-level (not
# self/class attributes) because the dispatchers run with a MainWindow `self`,
# not a UIFunctions instance.
#
# _MAINBAR_TABS:  key -> (stackedWidget index, active-button attr on self.ui)
_MAINBAR_TABS = {
    'import': (0, 'btn_import'),
    'filter': (1, 'btn_filter'),
    'params': (2, 'btn_parameters'),
    'plot':   (3, 'btn_plots'),
    'info':   (4, 'btn_info'),
    'search': (5, 'btn_search'),
}
# _PLOTBAR_TABS: key -> (infobar index, plot index, active-button attr,
#                        [(action, dialog.ui widget attr), ...] applied after reset)
_PLOTBAR_TABS = {
    'review':  (0, 1, 'btn_review',  []),
    'upset':   (1, 9, 'btn_upset',   [('show', 'frame_colorscheme')]),
    'dend':    (1, 2, 'btn_dend',    []),
    'pca':     (1, 3, 'btn_pca',     [('show', 'frame_2')]),
    'mzrt':    (0, 4, 'btn_mzrt',    []),
    'kmd':     (0, 5, 'btn_kmd',     [('show', 'frame_mdguide')]),
    '3dfc':    (0, 6, 'btn_3dfc',    []),
    'volcano': (0, 7, 'btn_volcano', [('show', 'frame_volcanoparams')]),
    'heatmap': (0, 8, 'btn_heatmap', [('show', 'frame_colorscheme')]),
}

# Nav-button stylesheets. Defined once here instead of being rebuilt on
# every reset_*bar() call (and read back off self.ui, where they were
# only ever used internally). The plot tab bar and the feature-info dialog
# bar shared a visually identical light style (the two former copies
# differed only by insignificant CSS whitespace), so they now share one
# pair. The left main-nav bar has its own (dark/transparent) theme.
_MAINBAR_ACTIVE_STYLE = 'QPushButton {\tborder: none; background-color: transparent;}QPushButton:hover {background-color: transparent}QPushButton:pressed {\tbackground-color: rgb(75, 75, 75);}'
_MAINBAR_INACTIVE_STYLE = 'QPushButton {\tborder: none; background-color: rgb(35, 35, 35);}QPushButton:hover {background-color: transparent}QPushButton:pressed {\tbackground-color: rgb(75, 75, 75);}'
_LIGHT_ACTIVE_STYLE = 'QPushButton {\tborder: none;background-color: rgba(225,225,225, 255);\tcolor: rgba(75,75,75,255)} QPushButton:hover {background-color: rgba(225,225,225, 255);\tcolor: rgba(75,75,75,255)} QPushButton:pressed {\t background-color: rgba(225,225,225, 255);\tcolor: rgba(75,75,75,255)}'
_LIGHT_INACTIVE_STYLE = 'QPushButton {\tborder: none;background-color: rgba(0,0,0, 0)} QPushButton:hover {background-color: rgba(225,225,225, 50)} QPushButton:pressed {\t background-color: rgba(225,225,225, 255);\tcolor: rgba(75,75,75,255)}'

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
        # The per-feature-set "Use" checkbox (and its label) next to the colour
        # picker is inert -- nothing in the hand-written code ever connects or
        # reads checkBox_use1, so it does nothing. Hide both at runtime rather
        # than editing the generated ui_main.py.
        self.ui.checkBox_use1.hide()
        self.ui.lbl_use.hide()
        # "Bootstrap Analysis" and "Collapse Technical Replicates" moved
        # from this global plot-config dialog onto their one relevant
        # plot's own switcher bar (plot_dendrogram's "Bootstrap" checkbox,
        # plot_ordination's "Collapse Replicates" checkbox) -- hide the
        # now-orphaned dialog widgets rather than editing the generated
        # ui_plotparam.py.
        self.dialog.ui.frame_bootstrap.hide()
        self.dialog.ui.frame_2.hide()
        # "Apply Data Filtering" is inert (its checked state is never read and
        # it has no signal connection), so hide it permanently rather than
        # leaving a no-op control on the config dialog. Hidden once here; the
        # per-tab show() that used to reveal it in reset_plotbar() is removed.
        self.dialog.ui.checkBox_applyfilter.hide()

        # Top bar functions
        self.ui.btn_maximize.clicked.connect(lambda: UIFunctions.maximize_restore(self))
        self.ui.btn_minimize.clicked.connect(lambda: self.showMinimized())
        self.ui.btn_close.clicked.connect(lambda: self.close())

        #mainbar functions
        self.ui.stackedWidget.setCurrentIndex(6)
        # Simple tab buttons go through the table-driven dispatcher (keyed by
        # _MAINBAR_TABS); 'search' additionally kicks off the db search, so it
        # keeps its own handler. (*_a absorbs the clicked(bool) arg; k=key binds
        # per-iteration so the closure doesn't capture the final loop value.)
        for key in ('import', 'filter', 'params', 'plot', 'info'):
            getattr(self.ui, _MAINBAR_TABS[key][1]).clicked.connect(
                lambda *_a, k=key: UIFunctions.goto_maintab(self, k))
        self.ui.btn_search.clicked.connect(lambda: UIFunctions.goto_search(self))

        #plotbar functions
        self.ui.stackedWidget_plot.setCurrentIndex(0)
        for key in _PLOTBAR_TABS:
            getattr(self.ui, _PLOTBAR_TABS[key][2]).clicked.connect(
                lambda *_a, k=key: UIFunctions.goto_plottab(self, k))
        
        self.ui.stackedWidget_review.setCurrentIndex(3)
        self.ui.btn_ftrplt.clicked.connect(lambda: self.ui.stackedWidget_review.setCurrentIndex(0))
        self.ui.btn_treemap.clicked.connect(lambda: self.ui.stackedWidget_review.setCurrentIndex(1))
        self.ui.btn_cvplt.clicked.connect(lambda: self.ui.stackedWidget_review.setCurrentIndex(2))
        self.ui.btn_datasummary.clicked.connect(lambda: self.ui.stackedWidget_review.setCurrentIndex(3))
        
        self.ui.btn_upsetplt.clicked.connect(lambda: UIFunctions.switch_grpanalysis_tab(self, 0))
        self.ui.btn_samplecorr.clicked.connect(lambda: UIFunctions.switch_grpanalysis_tab(self, 1))

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
    def goto_maintab(self, key):
        """Switch the main left-nav stacked widget to the page for ``key`` and
        highlight its button. Data-driven via ``_MAINBAR_TABS``; replaces the
        former one-method-per-tab goto_import/filter/params/plot/info."""
        index, button_attr = _MAINBAR_TABS[key]
        self.ui.stackedWidget.setCurrentIndex(index)
        UIFunctions.reset_mainbar(self)
        getattr(self.ui, button_attr).setStyleSheet(_MAINBAR_ACTIVE_STYLE)

    def goto_search(self):
        """Main-nav 'Search' tab: switch like any other main tab, then run the
        one-time NPAtlas db search + feature-tree fill once an analysis exists."""
        UIFunctions.goto_maintab(self, 'search')
        #self.ftrdialog.show()
        if self.dbsearchdone == False and self.analysisrun:
            start_functime()
            self.fulldbsearch()
            self.fillfttree()
            self.dbsearchdone = True
            stop_functime('dbsearch complete')
            reset_runtime()
        elif not self.analysisrun:
            self.error('Run an analysis before searching.')

    def switch_grpanalysis_tab(self, idx):
        """Switch the Group Analysis sub-tab (UpSet Plot=0, Sample
        Correlations=1) and grey out plot_samplecorr's Method/View/Use-Names
        controls -- shared with btn_upsetplt/btn_samplecorr in frame_12 --
        whenever the UpSet Plot tab is active, since they don't apply there."""
        self.ui.stackedWidget_grpanalysis.setCurrentIndex(idx)
        if getattr(self, 'samplecorr', None) is not None:
            self.samplecorr.set_controls_enabled(idx == 1)

    #plotbar functions
    def goto_plottab(self, key):
        """Switch the plots area (infobar + plot stacked widgets) to the page
        for ``key``, highlight its button, and apply that tab's dialog-frame
        show/hide deltas. Data-driven via ``_PLOTBAR_TABS``; replaces the former
        one-method-per-tab goto_review/upset/dend/pca/mzrt/kmd/3dfc/volcano/
        heatmap. Order matches the originals: set pages -> reset bar (which
        hides every per-tab frame) -> mark the active button -> apply this
        tab's frame deltas on top."""
        infobar_index, plot_index, button_attr, frame_ops = _PLOTBAR_TABS[key]
        self.ui.stackedWidget_infobar.setCurrentIndex(infobar_index)
        self.ui.stackedWidget_plot.setCurrentIndex(plot_index)
        UIFunctions.reset_plotbar(self)
        getattr(self.ui, button_attr).setStyleSheet(_LIGHT_ACTIVE_STYLE)
        for action, widget_attr in frame_ops:
            getattr(getattr(self.dialog.ui, widget_attr), action)()

    #ftinfobar functions
    def goto_abund(self):
        self.ftrdialog.ui.btn_masst.hide()
        self.ftrdialog.ui.stackedWidget.setCurrentIndex(2)
        # Refresh (not re-select/toggle) the current feature's display now
        # that the abundance tab is active.
        self._refresh_highlight()
        UIFunctions.reset_ftrdialogbar(self)
        self.ftrdialog.ui.btn_abund.setStyleSheet(_LIGHT_ACTIVE_STYLE)

    def goto_hits(self):
        self.ftrdialog.ui.btn_masst.hide()
        self.ftrdialog.ui.stackedWidget.setCurrentIndex(0)

        self._refresh_highlight()
        UIFunctions.reset_ftrdialogbar(self)
        self.ftrdialog.ui.btn_hits.setStyleSheet(_LIGHT_ACTIVE_STYLE)
        
    def goto_spectrum(self):
        self.ftrdialog.ui.btn_masst.show()
        self.ftrdialog.ui.stackedWidget.setCurrentIndex(1)
        UIFunctions.reset_ftrdialogbar(self)
        self.ftrdialog.ui.btn_spectrum.setStyleSheet(_LIGHT_ACTIVE_STYLE)

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
        # select() also sets the model's selected_index (side effect we rely on);
        # its return value is unused here.
        self.groupsetmodel.select(self.ui.listWidget_pltgrps.currentRow())

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
            selected = QFileDialog.getExistingDirectory(self, 'Select Directory', self.recentdir)
            if not selected:
                # Dialog cancelled -- leave the previously-selected output dir
                # (and recentdir) untouched rather than clobbering them.
                self.error('No directory selected')
                return
            self.outputdir = Path(selected)
            self.recentdir = str(self.outputdir)
            # setText needs a str (passing a Path raises); show the full path
            # when short, otherwise truncate from the left so the informative
            # tail stays visible.
            shown = str(self.outputdir)
            self.ui.lbl_outdir.setText(shown if len(shown) < 40 else '...' + shown[-40:])

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
        self.ui.btn_import.setStyleSheet(_MAINBAR_INACTIVE_STYLE)
        self.ui.btn_filter.setStyleSheet(_MAINBAR_INACTIVE_STYLE)
        self.ui.btn_parameters.setStyleSheet(_MAINBAR_INACTIVE_STYLE)
        self.ui.btn_plots.setStyleSheet(_MAINBAR_INACTIVE_STYLE)
        self.ui.btn_info.setStyleSheet(_MAINBAR_INACTIVE_STYLE)
        self.ui.btn_search.setStyleSheet(_MAINBAR_INACTIVE_STYLE)
        self.ui.label_status.setText('') #to reset analysis status after tabs are switched. eventually use a second thread to dynamically update status
        self.ui.label_status.setStyleSheet('color: rgb(150,150,150);')
        
    def reset_plotbar(self):
        self.ui.btn_review.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ui.btn_dend.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ui.btn_upset.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ui.btn_pca.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ui.btn_mzrt.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ui.btn_kmd.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ui.btn_3dfc.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ui.btn_volcano.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ui.btn_heatmap.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        
        self.dialog.ui.frame_2.hide()
        self.dialog.ui.frame.hide() #hide apply button, button doesnt work right now, use run button
        self.dialog.ui.frame_mdguide.hide()
        self.dialog.ui.frame_volcanoparams.hide()
        self.dialog.ui.frame_colorscheme.hide()
        self.ui.label_status.setText('')

    def reset_ftrdialogbar(self):
        self.ftrdialog.ui.btn_hits.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ftrdialog.ui.btn_spectrum.setStyleSheet(_LIGHT_INACTIVE_STYLE)
        self.ftrdialog.ui.btn_abund.setStyleSheet(_LIGHT_INACTIVE_STYLE)
    
    def show_ftrdialog(self):
        self.ftrdialog.show()
        if self.pickedfeature != '':
            self._refresh_highlight()

