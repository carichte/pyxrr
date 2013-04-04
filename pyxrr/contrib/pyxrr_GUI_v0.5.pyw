#!/usr/bin/env python

# Copyright 2013 by Hartmut Stoecker
# Contact: hartmut.stoecker@physik.tu-freiberg.de
#
# PyXRR_GUI provides a graphical user interface to PyXRR by Carsten Richter.
#
# The present version is 0.5.
# It is only working with pyxrr 0.9.07.

import os, wx, sys
import pyxrr
import matplotlib
matplotlib.use('WXAgg')

from pylab import arange, array, log10, savetxt, setp, sqrt, vstack
from copy import deepcopy
from wx.lib.agw import ultimatelistctrl as ULC

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2WxAgg


class MainFrame(wx.Frame):
    """ The main frame of the application. """
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, 'Anpassen von Reflektometrie-Daten')
        
        # Attributes of class
        self.filename = ""
        self.path = ""
        self.tempfile = "Current_Model.param"
        self.angle = arange(0, 5, 0.01)
        self.int = array([])
        self.start = array([])
        self.fit = array([])
        self.measparams_changed = 0
        self.ranges = {}
        self.Ns = 20
        
        # Dictionaries to correlate numbers
        self.group_layer = {}
        self.layer_sigma = {}
        self.group_sigma = {}
        
        # Build GUI
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        
        # Load last model or default values
        try:
            self.open_model(self.tempfile)
        except:
            self.save_model(self.tempfile, ['Ambience: name=air, code=N0.78O.21, rho=0.00125', 'Group: name=multilayer, sigma=1, periods=1, grad_d=0', 'Layer: name=Anatase, code=TiO2, rho=4.26, d=250', 'Substrate: name=Glass, code=SiO2, rho=2.2, sigma=1', 'Measurement: x_axis=twotheta, fit_range=0.2->100, energy=8.04116, resolution=0.02, offset=0.0, scale=1.0, background=-7.0, pol=0.5'])
            self.open_model(self.tempfile)

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_load = menu_file.Append(-1, "&Daten/Modell laden...\tCtrl-O", "Daten aus Datei laden")
        self.Bind(wx.EVT_MENU, self.on_open_file, m_load)
        menu_file.AppendSeparator()
        m_runfit = menu_file.Append(-1, "&Fit starten!\tCtrl-F", "Fit-Prozedur starten")
        self.Bind(wx.EVT_MENU, self.on_run_fit, m_runfit)     
        m_redrawmodel = menu_file.Append(-1, "&Grafik neu zeichnen\tCtrl-R", "Grafik erneut zeichnen")
        self.Bind(wx.EVT_MENU, self.on_draw_model, m_redrawmodel)
        menu_file.AppendSeparator()
        m_saveplot = menu_file.Append(-1, "&Grafik speichern...\tCtrl-G", "Grafik als Datei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_saveplot)
        m_savemodel = menu_file.Append(-1, "&Modell speichern...\tCtrl-M", "Modell als Datei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_model, m_savemodel)
        m_savetext = menu_file.Append(-1, "&Ergebnis speichern...\tCtrl-S", "Ergebnis als Datei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_text, m_savetext)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "&Beenden\tCtrl-X", "Programm verlassen")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&Hilfe\tF1", "Hilfe zum Programm")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.menubar.Append(menu_file, "&Datei")
        self.menubar.Append(menu_help, "&Hilfe")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it. """
        self.panel = wx.Panel(self)
        
        self.b_load = wx.Button(self.panel, -1, "Daten/Modell laden...", size=(170,25))
        self.b_load.Bind(wx.EVT_BUTTON, self.on_open_file)
        
        self.t_anglelabel = wx.StaticText(self.panel, -1, "x-Achse:", size=(60,-1))
        self.cb_angle = wx.ComboBox(self.panel, -1, choices=['omega', '2theta', 'q_z (A)', 'q_z (nm)'], style=wx.CB_READONLY, size=(102,-1))
        self.cb_angle.SetValue('omega')
        self.cb_angle.Bind(wx.EVT_COMBOBOX, self.on_change_measparams)
        
        self.t_pollabel = wx.StaticText(self.panel, -1, "Polarisation:", size=(60,-1))
        self.cb_pol = wx.ComboBox(self.panel, -1, choices=['unpolarisiert', 'parallel', 'senkrecht'], style=wx.CB_READONLY, size=(102,-1))
        self.cb_pol.SetValue('unpolarisiert')
        self.cb_pol.Bind(wx.EVT_COMBOBOX, self.on_change_measparams)
             
        self.t_weight = wx.StaticText(self.panel, -1, "Wichtung:", size=(60,-1))
        self.cb_weight = wx.ComboBox(self.panel, -1, choices=['keine Wichtung', 'dritte Datenspalte', 'statistisch'], style=wx.CB_READONLY, size=(102,-1))
        self.cb_weight.SetValue('keine Wichtung')
        self.cb_weight.Bind(wx.EVT_COMBOBOX, self.on_change_measparams)
        
        self.t_algorithm = wx.StaticText(self.panel, -1, "Algorithmus:", size=(60,-1))
        self.cb_algorithm = wx.ComboBox(self.panel, -1, choices=['Least Squares (Std.)', 'Brute Force', 'Simulated Annealing', 'Simplex', 'fmin_bfgs', 'fmin_powell', 'fmin_cg'], style=wx.CB_READONLY, size=(102,-1))
        self.cb_algorithm.SetValue('Least Squares (Std.)')
        
        self.t_startlabel = wx.StaticText(self.panel, -1, "Start:")
        self.t_start = wx.TextCtrl(self.panel, -1, '0.0', style=wx.TE_PROCESS_ENTER, size=(45,20))    
        self.t_start.Bind(wx.EVT_TEXT, self.apply_fit_range)
        self.t_start.Bind(wx.EVT_TEXT_ENTER, self.on_draw_model)
        
        self.t_endlabel = wx.StaticText(self.panel, -1, " Ende:")
        self.t_end = wx.TextCtrl(self.panel, -1, '100.0', style=wx.TE_PROCESS_ENTER, size=(45,20))    
        self.t_end.Bind(wx.EVT_TEXT, self.apply_fit_range)
        self.t_end.Bind(wx.EVT_TEXT_ENTER, self.on_draw_model)
        
        self.lb_model = ULC.UltimateListCtrl(self.panel, -1, agwStyle=ULC.ULC_REPORT | ULC.ULC_HAS_VARIABLE_ROW_HEIGHT | ULC.ULC_SINGLE_SEL, size=(178,180))
        self.lb_model.InsertColumn(0, "Name", width=70)
        self.lb_model.InsertColumn(1, "Material", width=87)
        
        self.bm_new = wx.Button(self.panel, -1, "Neu", size=(45,25))
        self.bm_new.Bind(wx.EVT_BUTTON, self.on_model_new)
        
        self.bm_del = wx.Button(self.panel, -1, "Entf.", size=(45,25))
        self.bm_del.Bind(wx.EVT_BUTTON, self.on_model_del)
    
        self.bm_up = wx.Button(self.panel, -1, "Hoch", size=(45,25))
        self.bm_up.Bind(wx.EVT_BUTTON, self.on_model_up)
        
        self.bm_down = wx.Button(self.panel, -1, "Runter", size=(45,25))
        self.bm_down.Bind(wx.EVT_BUTTON, self.on_model_down)
        
        self.bm_density = wx.Button(self.panel, -1, "Grafik", size=(45,25))
        self.bm_density.Bind(wx.EVT_BUTTON, self.on_show_density)
        
        self.lb_table = ULC.UltimateListCtrl(self.panel, -1, agwStyle=ULC.ULC_REPORT | ULC.ULC_HAS_VARIABLE_ROW_HEIGHT | ULC.ULC_NO_HIGHLIGHT, size=(430,150))
        self.lb_table.InsertColumn(0, "Parameter (Fit?)", width=160)
        self.lb_table.InsertColumn(1, "Startwert", width=68)
        self.lb_table.InsertColumn(2, "", width=22)
        self.lb_table.InsertColumn(3, "", width=22)
        self.lb_table.InsertColumn(4, "Fitergebnis", width=77)
        self.lb_table.InsertColumn(5, "Fehler", width=62)
       
        self.b_runfit = wx.Button(self.panel, -1, "Fit starten!", size=(120,35))
        self.b_runfit.Bind(wx.EVT_BUTTON, self.on_run_fit)
        
        self.b_redrawfigure = wx.Button(self.panel, -1, "Grafik neu\nzeichnen", size=(70,35))
        self.b_redrawfigure.Bind(wx.EVT_BUTTON, self.on_draw_model)
        
        self.b_saveplot = wx.Button(self.panel, -1, "Grafik\nspeichern", size=(70,35))
        self.b_saveplot.Bind(wx.EVT_BUTTON, self.on_save_plot)
        
        self.b_savemodel = wx.Button(self.panel, -1, "Modell\nspeichern", size=(70,35))
        self.b_savemodel.Bind(wx.EVT_BUTTON, self.on_save_model)
        
        self.b_savetext = wx.Button(self.panel, -1, "Ergebnis\nspeichern", size=(70,35))
        self.b_savetext.Bind(wx.EVT_BUTTON, self.on_save_text)
        
        self.dpi = 100
        self.fig = Figure((7, 5), dpi=self.dpi)
        self.canvas = FigureCanvasWxAgg(self.panel, -1, self.fig)
        self.canvas.mpl_connect('motion_notify_event', self.on_UpdateCursor)
        self.canvas.mpl_connect('resize_event', self.on_Resize)
        
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        
        self.bm_log = wx.Button(self.panel, -1, "Log verbergen", size=(88,25))
        self.bm_log.Bind(wx.EVT_BUTTON, self.on_hide_log)
        
        self.log = wx.TextCtrl(self.panel, -1, size=(300,125), style = wx.TE_MULTILINE|wx.TE_READONLY)
        self.log.SetFont(wx.Font(8, wx.MODERN, wx.NORMAL, wx.NORMAL, 0, "Sans"))
        sys.stdout = RedirectText(self.log, 'black')
        sys.stderr = RedirectText(self.log, 'red')
       
        # Layout with box sizers ----------------------------------------------
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        
        self.vbox11 = wx.BoxSizer(wx.VERTICAL)
        self.vbox11.Add(self.b_load, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL)
        
        self.hbox111 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox111.Add(self.t_anglelabel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox111.Add(self.cb_angle, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.hbox112 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox112.Add(self.t_pollabel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox112.Add(self.cb_pol, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.hbox113 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox113.Add(self.t_weight, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox113.Add(self.cb_weight, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.hbox114 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox114.Add(self.t_algorithm, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox114.Add(self.cb_algorithm, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.hbox115 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox115.Add(self.t_startlabel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox115.Add(self.t_start, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox115.Add(self.t_endlabel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox115.Add(self.t_end, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.vbox11.AddSpacer(4)
        self.vbox11.Add(self.hbox111, 0, border=2, flag = wx.ALIGN_LEFT | wx.ALL)        
        self.vbox11.Add(self.hbox112, 0, border=2, flag = wx.ALIGN_LEFT | wx.ALL)
        self.vbox11.Add(self.hbox113, 0, border=2, flag = wx.ALIGN_LEFT | wx.ALL)
        self.vbox11.Add(self.hbox114, 0, border=2, flag = wx.ALIGN_LEFT | wx.ALL)
        self.vbox11.AddSpacer(4)
        self.vbox11.Add(self.hbox115, 0, border=2, flag = wx.ALIGN_LEFT | wx.ALL)
        
        self.vbox12 = wx.BoxSizer(wx.VERTICAL)
        self.vbox12.AddSpacer(25)
        self.vbox12.Add(self.bm_new, 0, border=3, flag = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
        self.vbox12.Add(self.bm_del, 0, border=3, flag = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
        self.vbox12.Add(self.bm_up, 0, border=3, flag = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
        self.vbox12.Add(self.bm_down, 0, border=3, flag = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
        self.vbox12.AddSpacer(13)
        self.vbox12.Add(self.bm_density, 0, border=3, flag = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
             
        self.hbox1.Add(self.vbox11, 0, border=5, flag = wx.ALIGN_LEFT | wx.ALL)
        self.hbox1.AddSpacer(5)
        self.hbox1.Add(self.lb_model, 0, border=5, flag = wx.ALIGN_LEFT | wx.ALL | wx.EXPAND)
        self.hbox1.Add(self.vbox12, 0, border=0, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.b_runfit, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(self.b_redrawfigure, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(self.b_saveplot, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(self.b_savemodel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(self.b_savetext, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.vbox1.Add(self.hbox1, 0, border=5, flag = wx.ALIGN_LEFT | wx.ALL | wx.EXPAND)
        self.vbox1.Add(self.lb_table, 1, border=5, flag = wx.ALIGN_LEFT | wx.ALL | wx.EXPAND)
        self.vbox1.Add(self.hbox2, 0, border=5, flag = wx.ALIGN_LEFT | wx.ALL | wx.EXPAND)
        
        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.vbox2.Add(self.canvas, 1, flag = wx.ALIGN_LEFT  | wx.TOP | wx.EXPAND)
        
        self.hbox21 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox21.Add(self.toolbar, 1, flag = wx.ALIGN_LEFT | wx.TOP | wx.EXPAND)
        self.hbox21.Add(self.bm_log, 0, flag = wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        self.hbox21.AddSpacer(5)
        
        self.vbox2.Add(self.hbox21, 0, flag = wx.EXPAND)
        self.vbox2.Add(self.log, 0, flag = wx.EXPAND)
             
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.vbox1, 0, flag = wx.ALIGN_LEFT | wx.TOP | wx.EXPAND)
        self.hbox.Add(self.vbox2, 1, flag = wx.ALIGN_LEFT | wx.TOP | wx.EXPAND)
        
        self.panel.SetSizer(self.hbox)
        self.hbox.Fit(self)
        
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetFieldsCount(3)
        self.statusbar.SetStatusWidths([-1,100,110])
    
    def update_model(self):
        for i in range(self.lb_model.GetItemCount()):
            self.lb_model.DeleteItem(0)
            
        self.new_layer(0, 'Umgebung', self.model[0])
        
        row = 1
        layer = 1
        for i in range(len(self.params['LayerCount']) - 2):
            self.new_layer(row, 'Gruppe ' + str(i+1), 'Gruppe')
            row += 1
            for j in range(self.params['LayerCount'][i+1]):
                self.new_layer(row, '> Schicht ' + str(layer), self.model[layer])
                row += 1
                layer += 1
            
        self.new_layer(row, 'Substrat', self.model[layer])
        
        self.lb_model.Update()
        self.update_table()
        
    def new_layer(self, row, label, name):
        row = self.lb_model.GetItemCount()
        self.lb_model.InsertStringItem(row, label)
        self.lb_model.SetItemData(row, label)
        if name == "Gruppe":
            name = str(self.params['N'][int(label.split()[-1])])
            self.lb_model.SetStringItem(row, 1, 'Perioden:')
            text = wx.TextCtrl(self.lb_model, value=name, id=1000+row, style=wx.TE_PROCESS_ENTER, size=(30,20))
            text.Bind(wx.EVT_TEXT, self.on_enter_periods)
            text.Bind(wx.EVT_TEXT_ENTER, self.on_draw_model)
            self.lb_model.SetItemWindow(row, 1, wnd=text)
        else:
            text = wx.TextCtrl(self.lb_model, value=name, id=1000+row, style=wx.TE_PROCESS_ENTER, size=(85,20))
            text.Bind(wx.EVT_TEXT, self.on_enter_material)
            text.Bind(wx.EVT_TEXT_ENTER, self.on_draw_model)
            self.lb_model.SetItemWindow(row, 1, wnd=text)

    def update_table(self):
        self.group_layer = {}
        self.layer_sigma = {}
        self.group_sigma = {}
        self.fit_keys = []

        for i in range(self.lb_table.GetItemCount()):
            if self.lb_table.GetItem(0).IsChecked():
                self.fit_keys.append(self.lb_table.GetItemData(0))
            self.lb_table.DeleteItem(0)
          
        self.new_param('Energie (keV)', 'energy0')
        button = wx.Button(self.lb_table, label="...", size=(23,20))
        button.Bind(wx.EVT_BUTTON, self.on_select_energy)
        self.lb_table.SetItemWindow(0, 0, wnd=button)
        
        self.new_param('Winkelgenauigkeit (deg)', 'resolution0')
        self.new_param('Winkelverschiebung (deg)', 'offset0')
        self.new_param('Skalierung', 'scale0')
        self.new_param('Konstanter Untergrund', 'background0')
        
        for i in range(len(self.params['LayerCount']) - 2):
            if self.params['N'][i+1] > 1:
                self.new_param('Dickengradient in Gruppe ' + str(i+1), 'grad_d_' + str(i+1))

        for i in range(sum(self.params['LayerCount'])):
            if i == 0:
                self.new_param('Dichte Umgebung (g/cm^3)', 'rho_' + str(i))
            elif i == sum(self.params['LayerCount']) - 1:
                self.new_param('Dichte Substrat (g/cm^3)', 'rho_' + str(i))
            else:
                self.new_param('Dichte Schicht ' + str(i) + ' (g/cm^3)', 'rho_' + str(i))
                
        for i in range(sum(self.params['LayerCount'])):
            if i > 0 and i < sum(self.params['LayerCount']) - 1:
                self.new_param('Dicke Schicht ' + str(i) + ' (A)', 'd_' + str(i))
        
        sigma_pos = 0
        layer_pos = 1
        for i in range(len(self.params['LayerCount']) - 2):
            self.group_layer[i+1] = []
            if self.params['N'][i+1] > 1:
                self.new_param('Rauigkeit Gruppe ' + str(i+1) + ' (A)', 'sigma_' + str(sigma_pos))
                self.group_sigma[i+1] = sigma_pos
                sigma_pos += 1
            multilayer_sigma = 0
            for j in range(self.params['LayerCount'][i+1]):
                self.group_layer[i+1].append(layer_pos)
                if j == 0:
                    if self.params['N'][i+1] == 1:
                        self.new_param('Rauigkeit Schicht ' + str(layer_pos) + ' (A)', 'sigma_' + str(sigma_pos))
                        self.layer_sigma[layer_pos] = sigma_pos
                    else:
                        special_sigma = sigma_pos + self.params['LayerCount'][i+1] - 1
                        self.new_param('Rauigkeit Schicht ' + str(layer_pos) + ' (A)', 'sigma_' + str(special_sigma))
                        self.layer_sigma[layer_pos] = special_sigma
                        sigma_pos -= 1
                        multilayer_sigma = 1
                else:
                    self.new_param('Rauigkeit Schicht ' + str(layer_pos) + ' (A)', 'sigma_' + str(sigma_pos))
                    self.layer_sigma[layer_pos] = sigma_pos
                layer_pos += 1
                sigma_pos += 1
            sigma_pos += multilayer_sigma
        self.new_param('Rauigkeit Substrat (A)', 'sigma_' + str(sigma_pos))
             
        self.lb_table.Update()
      
    def new_param(self, name, key):
        row = self.lb_table.GetItemCount()
        self.lb_table.InsertStringItem(row, ' ' + name, it_kind=1)
        self.lb_table.SetItemData(row, key)
        
        self.lb_table.SetStringItem(row, 4, '%.6g' %self.params_new[key])
        self.lb_table.SetStringItem(row, 5, '%.6g' %self.errors[key])
        text = wx.TextCtrl(self.lb_table, value=str(self.params[key]), id=2000+row, style=wx.TE_PROCESS_ENTER, size=(66,20))
        text.Bind(wx.EVT_TEXT, self.on_enter_value)
        text.Bind(wx.EVT_TEXT_ENTER, self.on_draw_model)
        self.lb_table.SetItemWindow(row, 1, wnd=text)
        
        button1 = wx.Button(self.lb_table, label="+", id=3000+row, size=(20,20))
        button1.Bind(wx.EVT_BUTTON, self.on_increase_value)
        self.lb_table.SetItemWindow(row, 2, wnd=button1)
        
        button2 = wx.Button(self.lb_table, label="-", id=4000+row, size=(20,20))
        button2.Bind(wx.EVT_BUTTON, self.on_decrease_value)
        self.lb_table.SetItemWindow(row, 3, wnd=button2)
        
        button3 = wx.Button(self.lb_table, label="<", id=5000+row, size=(20,20))
        button3.Bind(wx.EVT_BUTTON, self.on_load_new)
        self.lb_table.SetItemWindow(row, 4, wnd=button3)
        
    def on_select_energy(self, event):
        dlg = SelectEnergy(None, -1, 'Energie bestimmen')
        
        if dlg.ShowModal() == wx.ID_OK:
            value = dlg.GetValue()
            self.FindWindowById(2000).SetValue(value)

        dlg.Destroy()
   
    def draw_figure(self, redraw=1):
        """ Redraws the figure. """
        self.save_model(self.tempfile)
        
        if redraw:
            self.fig.clear()
            self.axes = self.fig.add_subplot(111)
            self.axes.grid(1)
            
            if self.int != array([]):
                self.plot_int = self.axes.semilogy(self.angle, self.int, 'b', label='Daten')
            else:
                self.plot_int = []
            if self.start != array([]):
                self.plot_start = self.axes.semilogy(self.angle, self.start, 'g', label='Startmodell')
            else:
                self.plot_start = []
            if self.fit != array([]):
                self.plot_fit = self.axes.semilogy(self.angle, self.fit, 'r', label='Fit')
            else:
                self.plot_fit = []
                
            prop = matplotlib.font_manager.FontProperties(size=10) 
            self.axes.legend(loc=0, prop=prop)
            self.axes.tick_params(axis='both', labelsize=10)
            self.axes.set_xlabel('Omega (deg)', fontsize=12)
            self.axes.set_ylabel('Intensitaet (normiert)', fontsize=12)
            if self.filename != '':
                self.titel = self.axes.set_title(self.filename, fontsize=12)
           
            self.on_Resize('')
            self.canvas.draw()
            
        else:
            if self.int != array([]):
                if self.plot_int != []:
                    setp(self.plot_int, ydata = self.int)
                else:
                    self.plot_int = self.axes.semilogy(self.angle, self.int, 'b', label='Daten')
            
            if self.start != array([]):
                if self.plot_start != []:
                    setp(self.plot_start, ydata = self.start)
                else:
                    self.plot_start = self.axes.semilogy(self.angle, self.start, 'g', label='Startmodell')
                    
            if self.fit != array([]):
                if self.plot_fit != []:
                    setp(self.plot_fit, ydata = self.fit)
                else:
                    self.plot_fit = self.axes.semilogy(self.angle, self.fit, 'r', label='Fit')
            
            prop = matplotlib.font_manager.FontProperties(size=10) 
            self.axes.legend(loc=0, prop=prop)
            self.canvas.draw()
 
    def on_open_file(self, event):
        file_choices = "Daten-Typen (*.asc, *.fio, *.njc, *.val, *.raw, *.x00, *.txt, *.dat, *.param)|*.asc;*.fio;*.njc;*.val;*.raw;*.x00;*.txt;*.dat;*.param|Modell-Datei (*.param)|*.param|ASC-Datei (*.asc)|*.asc|Hasylab (*.fio)|*.fio|Seifert (*.njc)|*.njc|Seifert (*.val)|*.val|Bruker (*.raw)|*.raw|Philips (*.x00)|*.x00|Text-Datei (*.txt)|*.txt|Daten-Datei (*.dat)|*.dat|Alle Dateien (*.*)|*.*"
        dlg = wx.FileDialog(self, "Daten laden", wildcard=file_choices, style=wx.OPEN)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.path = dlg.GetPath()
            
            try:
                if os.path.splitext(self.filename)[-1] == '.param':
                    self.open_model(self.path)
                else:
                    self.save_model(self.tempfile)
                    self.sample.__init__(self.tempfile)           
                    self.angle = self.sample.measured_data[0][:,0]
                    self.int = self.sample.measured_data[0][:,1]
                    self.start = self.sample.reflectogram(self.angle, 0)
                    self.fit = array([])
                    self.flash_status_message("%s geladen." % self.filename)
                    self.draw_figure()
                    
            except Exception as error:
                wx.MessageBox('Fehler beim Laden der Datei:\n' + self.path + '\n\nFehlermeldung:\n' + str(error), 'Fehler', style=wx.ICON_ERROR)

        dlg.Destroy()
        
    def open_model(self, openpath):        
        self.sample = pyxrr.multilayer(openpath)
        self.params = deepcopy(self.sample.parameters)
        self.params_new = deepcopy(self.params)
        self.reset_errors()
        self.model = self.sample.materials
        
        self.path = self.sample.paths[0]
        self.filename = os.path.split(self.path)[-1]
        
        value = self.sample.x_axes[0]
        if value == 'qz_a':
            self.cb_angle.SetValue('q_z (A)')
        elif value == 'qz_nm':
            self.cb_angle.SetValue('q_z (nm)')
        elif value == 'twotheta':
            self.cb_angle.SetValue('2theta')
        else:
            self.cb_angle.SetValue('omega')
    
        value = self.sample.pol[0]
        if value == 0:
            self.cb_pol.SetValue('senkrecht')
        elif value == 1:
            self.cb_pol.SetValue('parallel')
        else:
            self.cb_pol.SetValue('unpolarisiert')
        
        if self.filename != '':
            if self.sample.weightmethods.has_key(0):
                value = self.sample.weightmethods[0]
                if value == 'statistical':
                    self.cb_weight.SetValue('statistisch') 
                elif value == 'z':
                    self.cb_weight.SetValue('dritte Datenspalte')
            else:
                self.cb_weight.SetValue('keine Wichtung')

        self.t_start.SetValue(str(self.sample.fit_limits[0][0]))
        self.t_end.SetValue(str(self.sample.fit_limits[0][1]))
        
        if self.filename != '':
            self.angle = self.sample.measured_data[0][:,0]
            self.int = self.sample.measured_data[0][:,1]
        else:
            self.angle = arange(0, 5, 0.01)
            
        self.start = self.sample.reflectogram(self.angle, 0)
        self.fit = array([])
        self.measparams_changed = 0
                
        self.flash_status_message("%s geladen." % self.filename)
        self.update_model()
        self.draw_figure()
        
        if self.sample.number_of_measurements > 1:
            wx.MessageBox('Modell umfasst ' + str(self.sample.number_of_measurements) + ' Dateien mit Messdaten, aber nur die erste wurde geladen.', 'Hinweis', style=wx.ICON_INFORMATION)
            self.sample.number_of_measurements = 1
            
    def on_save_model(self, event):
        file_choices = "Parameter-Datei (*.param)|*.param"
        filename = os.path.splitext(self.filename)[0] + "_model.param"
        dlg = wx.FileDialog(self, message="Modell speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            try:
                self.save_model(path)
            except Exception as error:
                wx.MessageBox('Fehler beim Speichern der Datei:\n' + path + '\n\nFehlermeldung:\n' + str(error), 'Fehler', style=wx.ICON_ERROR)
            
        dlg.Destroy()
            
    def save_model(self, savepath, modell=''):
        if self.measparams_changed:
            value = self.cb_angle.GetValue()
            if value == 'q_z (A)':
                self.sample.x_axes[0] = 'qz_a'
            elif value == 'q_z (nm)':
                self.sample.x_axes[0] = 'qz_nm'
            elif value == '2theta':
                self.sample.x_axes[0] = 'twotheta'
            else:
                self.sample.x_axes[0] = 'theta'
        
            value = self.cb_pol.GetValue()
            if value == 'senkrecht':
                self.sample.pol[0] = 0
            elif value == 'parallel':
                self.sample.pol[0] = 1
            else:
                self.sample.pol[0] = 0.5
                
            value = self.cb_weight.GetValue()
            if value == 'statistisch':
                self.sample.weightmethods[0] = 'statistical'
            elif value == 'dritte Datenspalte':
                self.sample.weightmethods[0] = 'z'
            else:
                self.sample.weightmethods[0] = 'no'
                
            self.sample.fit_limits[0] = float(self.t_start.GetValue()), float(self.t_end.GetValue())
            
            self.sample.materials = deepcopy(self.model)
            
            self.measparams_changed = 0
        
        if modell == '':
            self.sample.paths[0] = self.path
            modell = self.sample.save_model().splitlines()
        
        f = open(savepath, 'w')
        for line in modell:
            f.write(line + '\n')
        f.close()
        
        self.flash_status_message("Modell gespeichert unter %s" % savepath)
        
    def make_model(self):
        string = self.sample.save_model()
        return filter(lambda x: x != '', string.splitlines())
        
    def apply_fit_range(self, event):
        control = event.GetEventObject()
        value = control.GetValue()
        try:
            float(value)
        except ValueError:
            pos = control.GetInsertionPoint()
            end = control.GetLastPosition()
            if pos < end:
                value = value[:pos-1] + value[pos:]
            else:
                value = value[:pos-1]
            try:
                float(value)
            except:
                value = '0'
            control.ChangeValue(value)
            control.SetInsertionPointEnd()
            
        self.sample.fit_limits[0] = float(self.t_start.GetValue()), float(self.t_end.GetValue())
        self.sample.process_fit_range()
        self.flash_status_message("Fit-Bereich aktualisiert.")
    
    def on_change_measparams(self, event):
        self.measparams_changed = 1
        self.on_draw_model(event)
    
    def update_fitvalues(self):  
        for i in range(self.lb_table.GetItemCount()):
            key = self.lb_table.GetItemData(i)
            self.lb_table.SetStringItem(i, 4, '%.6g' %self.params_new[key])
            self.lb_table.SetStringItem(i, 5, '%.6g' %self.errors[key])
    
    def reset_errors(self):
        self.errors = deepcopy(self.params_new)
        for key in self.errors:
            if key != 'N' and key != 'LayerCount':
                self.errors[key] = 0
    
    def on_run_fit(self, event):
        self.on_draw_model(event)
        self.reset_errors()
        
        self.fit_keys = []       
        for i in range(self.lb_table.GetItemCount()):
            if self.lb_table.GetItem(i, col=0).IsChecked():
                key = self.lb_table.GetItemData(i)
                self.fit_keys.append(key)
        
        if self.fit_keys == []:
            self.sample.var_names = []
        else:
            try:
                fit_ranges = []
                value = self.cb_algorithm.GetValue()
                if value == 'Least Squares (Std.)':
                    alg = 'leastsq'                            
                elif value == 'Simulated Annealing':
                    alg = 'anneal'
                elif value == 'Simplex':
                    alg = 'fmin'
                elif value == 'fmin_bfgs':
                    alg = 'fmin_bfgs'
                elif value == 'fmin_powell':
                    alg = 'fmin_powell'
                elif value == 'fmin_cg':
                    alg = 'fmin_cg'
                elif value == 'Brute Force':
                    alg = 'brute'
                    dlg = SelectRanges(None, -1, 'Fitbereiche bestimmen')
                    if dlg.ShowModal() == wx.ID_OK:
                        fit_ranges, self.Ns = dlg.GetValue()
                    dlg.Destroy()
                
                fitted_param = self.sample.fit(algorithm=alg, var_names=self.fit_keys, ranges=fit_ranges, Ns=self.Ns)
                self.errors = self.sample.fiterrors
            
            except Exception as error:
                wx.MessageBox('Problem beim Fitten:\n\n' + str(error), 'Fehler', style=wx.ICON_ERROR)
            
        self.params_new = deepcopy(self.sample.parameters)
        self.fit = self.sample.reflectogram(self.angle, 0)
        self.update_fitvalues()
        self.draw_figure(redraw=0)
        
        if self.fit_keys == []:
            self.flash_status_message("Zeichnen des Modells erfolgreich!")
        elif len(self.fit_keys) == 1:
            self.flash_status_message("Fitten eines Parameters erfolgreich!")
        else:
            self.flash_status_message("Fitten mit %.d Parametern erfolgreich!" %len(self.fit_keys))
            
    def on_draw_model(self, event):
        if self.measparams_changed:
            self.save_model(self.tempfile)
            self.open_model(self.tempfile)
            if self.filename == '':
                self.angle = arange(0, 5, 0.01)
            else:
                self.angle = self.sample.measured_data[0][:,0]
                self.int = self.sample.measured_data[0][:,1]
        else:
            self.sample.parameters = deepcopy(self.params)
        
        self.start = self.sample.reflectogram(self.angle, 0)
        
        try:
            x = event.GetId()
        except:
            x = -1
        
        if x > 0:
            self.draw_figure(redraw=0)
        else:
            self.draw_figure(redraw=1)
        
        self.flash_status_message("Zeichnen des Modells erfolgreich!")
        
    def on_model_new(self, event):
        index = self.lb_model.GetFirstSelected()
        
        if (index > 0 and index < self.lb_model.GetItemCount() - 1) or self.lb_model.GetItemCount() < 3:
            insert_group = 1
            if self.lb_model.GetItemCount() < 3:
                index = 1
            else:
                label = self.lb_model.GetItemData(index)
                if "Schicht" in label:
                    insert_group = 0
                
            modell = self.make_model()
            modell.insert(index, "Layer: name=Schicht_0, code=Si, rho=2.0, d=10.0, sigma=1.0\n")
            if "Gruppe_" in modell[index-1]:
                periods = int(modell[index-1].split("periods=")[-1].split(",")[0].strip())
                if periods == 1:
                    start, end = modell[index-1].split("sigma=")
                    numstr = end.split(",")[0]
                    sigma = numstr.strip()
                    modell[index-1] = start + "sigma=1.0" + end.lstrip(numstr)
                    modell[index+1] = modell[index+1].rstrip() + ", sigma=" + sigma + "\n"
            
            if insert_group:
                modell.insert(index, "\nGroup: name=Gruppe_0, sigma=1.0, periods=1, grad_d=0\n")
            
            self.save_model(self.tempfile, modell)
            self.open_model(self.tempfile)
            
        else:
            wx.MessageBox('Bitte Gruppe oder Schicht markieren!', 'Fehler', style=wx.ICON_ERROR)

    def on_model_del(self, event):
        index = self.lb_model.GetFirstSelected()
            
        if index > 0 and index < self.lb_model.GetItemCount() - 1:
            label = self.lb_model.GetItemData(index)
            modell = self.make_model()
            
            if "Gruppe" in label:
                del modell[index]
                while modell[index][:5] == 'Layer':
                    del modell[index]
                
            elif "Schicht" in label:
                del modell[index]
                if "Gruppe_" in modell[index-1]:
                    periods = int(modell[index-1].split("periods=")[-1].split(",")[0].strip())
                    if periods == 1 and "Schicht_" in modell[index]:
                        sigma = modell[index].split("sigma=")[-1].split(",")[0].strip()
                        start, end = modell[index-1].split("sigma=")
                        numstr = end.split(",")[0]
                        modell[index-1] = start + "sigma=" + sigma + end.lstrip(numstr)
                    elif "Gruppe_" in modell[index] or "Substrate" in modell[index]:
                        del modell[index-1]
                            
            self.save_model(self.tempfile, modell)
            self.open_model(self.tempfile)
            
        else:
            wx.MessageBox('Bitte Gruppe oder Schicht markieren!', 'Fehler', style=wx.ICON_ERROR)
       
    def on_model_up(self, event):
        index = self.lb_model.GetFirstSelected()
        
        if index > 0 and index < self.lb_model.GetItemCount() - 1:
            label = self.lb_model.GetItemData(index)
            num = int(label.split()[-1])
            modell = self.make_model()
            
            if "Gruppe" in label and num > 1:
                for i in range(index-1, 0, -1):
                    if "Gruppe_" in modell[i]:
                        lastindex = i
                        break
                for i in range(index+1, self.lb_model.GetItemCount()):
                    if "Gruppe_" in modell[i] or "Substrate" in modell[i]:
                        nextindex = i
                        break
                modell[lastindex:nextindex] = modell[index:nextindex] + modell[lastindex:index]
                self.save_model(self.tempfile, modell)
                self.open_model(self.tempfile)
                    
            elif "Schicht" in label and "Schicht_" in modell[index-1]:
                if "Gruppe_" in modell[index-2]:
                    periods = int(modell[index-2].split("periods=")[-1].split(",")[0].strip())
                    if periods == 1:
                        sigma_group = modell[index-2].split("sigma=")[-1].split(",")[0].strip()
                        sigma_layer = modell[index].split("sigma=")[-1].split(",")[0].strip()
                        start, end = modell[index-2].split("sigma=")
                        numstr = end.split(",")[0]
                        modell[index-2] = start + "sigma=" + sigma_layer + end.lstrip(numstr)
                        modell[index-1] = modell[index-1].rstrip() + ", sigma=" + sigma_group + "\n"
                modell[index-1], modell[index] = modell[index], modell[index-1]
                self.save_model(self.tempfile, modell)
                self.open_model(self.tempfile)
                    
            else:
                wx.MessageBox('Verschieben von ' + label.lstrip('> ') + ' nach oben nicht sinnvoll!', 'Fehler', style=wx.ICON_ERROR)
                
        else:
            wx.MessageBox('Bitte Gruppe oder Schicht markieren!', 'Fehler', style=wx.ICON_ERROR)
            
    def on_model_down(self, event):
        index = self.lb_model.GetFirstSelected()
        
        if index > 0 and index < self.lb_model.GetItemCount() - 1:
            label = self.lb_model.GetItemData(index)
            num = int(label.split()[-1])
            modell = self.make_model()
            
            if "Gruppe" in label and num < len(self.params['LayerCount'])-2:
                for i in range(index+1, self.lb_model.GetItemCount()):
                    if "Gruppe_" in modell[i]:
                        nextindex = i
                        break
                for i in range(nextindex+1, self.lb_model.GetItemCount()):
                    if "Gruppe_" in modell[i] or "Substrate" in modell[i]:
                        overnextindex = i
                        break
                modell[index:overnextindex] = modell[nextindex:overnextindex] + modell[index:nextindex]
                self.save_model(self.tempfile, modell)
                self.open_model(self.tempfile)
                    
            elif "Schicht" in label and "Schicht_" in modell[index+1]:
                if "Gruppe_" in modell[index-1]:
                    periods = int(modell[index-1].split("periods=")[-1].split(",")[0].strip())
                    if periods == 1:
                        sigma_group = modell[index-1].split("sigma=")[-1].split(",")[0].strip()
                        sigma_layer = modell[index+1].split("sigma=")[-1].split(",")[0].strip()
                        start, end = modell[index-1].split("sigma=")
                        numstr = end.split(",")[0]
                        modell[index-1] = start + "sigma=" + sigma_layer + end.lstrip(numstr)
                        modell[index] = modell[index].rstrip() + ", sigma=" + sigma_group + "\n"
                modell[index], modell[index+1] = modell[index+1], modell[index]
                self.save_model(self.tempfile, modell)
                self.open_model(self.tempfile)

            else:
                wx.MessageBox('Verschieben von ' + label.lstrip('> ') + ' nach unten nicht sinnvoll!', 'Fehler', style=wx.ICON_ERROR)
        
        else:
            wx.MessageBox('Bitte Gruppe oder Schicht markieren!', 'Fehler', style=wx.ICON_ERROR)
            
    def on_show_density(self, event):
        stack = self.sample.stack()
        z = arange(-20, stack[-1,1]+20)
        rho, delta, beta = self.sample.density(z)
        
        dlg = DensityPlot(None, -1, 'Darstellung des Modells')
        dlg.plot_density(z, rho)
        dlg.ShowModal()
            
    def on_increase_value(self, event):
        id = event.GetId()
        value = float(self.FindWindowById(id-1000).GetValue())
        key = self.lb_table.GetItemData(id-3000)

        absval = abs(value)
        if absval >= 1:
            add = 10**(int(log10(absval))-1)
        elif absval == 0:
            add = 0.1
        else:
            add = 10**(int(log10(absval))-2)
                
        self.FindWindowById(id-1000).SetValue(str(value+add))
        self.on_draw_model(event)
        
    def on_decrease_value(self, event):
        id = event.GetId()
        value = float(self.FindWindowById(id-2000).GetValue())
        key = self.lb_table.GetItemData(id-4000)
        
        absval = abs(value)
        if absval >= 1:
            add = 10**(int(log10(absval))-1)
        elif absval == 0:
            add = 0.1
        else:
            add = 10**(int(log10(absval))-2)
            
        self.FindWindowById(id-2000).SetValue(str(value-add))
        self.on_draw_model(event)
        
    def on_load_new(self, event):
        id = event.GetId()
        value = self.lb_table.GetItem(id-5000, 4).GetText()
        self.FindWindowById(id-3000).SetValue(value)
        self.on_draw_model(event)
        
    def on_enter_material(self, event):
        id = event.GetId()
        value = self.FindWindowById(id).GetValue()
        material = self.lb_model.GetItem(id-1000, 0).GetText()
        if material == 'Umgebung':
            self.model[0] = value
        elif material == 'Substrat':
            self.model[-1] = value
        else:
            num = int(material.split()[-1])
            self.model[num] = value
        self.measparams_changed = 1
        
    def on_enter_periods(self, event):
        id = event.GetId()
        index = id - 1000
        value = self.FindWindowById(id).GetValue()
        material = self.lb_model.GetItem(index, 0).GetText()
        num = int(material.split()[-1])

        try:
            new_value = int(float(value))
            old_value = self.params['N'][num]
            if new_value != float(value) or new_value < 1:
                self.FindWindowById(id).SetValue(str(old_value))
            else:
                if new_value == 1 and old_value != 1:
                    modell = self.make_model()
                    start, end = modell[index].split("periods=")
                    numstr = end.split(",")[0]
                    modell[index] = start + "periods=" + str(new_value) + end.lstrip(numstr)
                    self.save_model(self.tempfile, modell)
                    self.open_model(self.tempfile)
                elif new_value != 1 and old_value == 1:
                    modell = self.make_model()
                    start, end = modell[index].split("periods=")
                    numstr = end.split(",")[0]
                    modell[index] = start + "periods=" + str(new_value) + end.lstrip(numstr)
                    sigma = modell[index].split("sigma=")[-1].split(",")[0].strip()
                    modell[index+1] = modell[index+1].rstrip() + ", sigma=" + sigma + "\n"
                    self.save_model(self.tempfile, modell)
                    self.open_model(self.tempfile)
                else:
                    self.params['N'][num] = new_value
        except:
            self.FindWindowById(id).SetValue(str(self.params['N'][num]))
        
    def on_enter_value(self, event):
        id = event.GetId()
        value = self.FindWindowById(id).GetValue()
        key = self.lb_table.GetItemData(id-2000)
        try:
            self.params[key] = float(value)
        except:
            self.FindWindowById(id).SetValue(str(self.params[key]))
        
    def on_save_plot(self, event):
        file_choices = "Portable Network Graphics (*.png)|*.png|Encapsulated Postscript (*.eps)|*.eps|Portable Document Format (*.pdf)|*.pdf|Scalable Vector Graphics (*.svg)|*.svg"
        filename = os.path.splitext(self.filename)[0] + "_fit.png"
        dlg = wx.FileDialog(self, message="Grafik speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            
            try:
                self.canvas.print_figure(path, dpi=3*self.dpi)
                self.flash_status_message("Grafik gespeichert unter %s" % path)
            except Exception as error:
                wx.MessageBox('Fehler beim Speichern der Datei:\n' + path + '\n\nFehlermeldung:\n' + str(error), 'Fehler', style=wx.ICON_ERROR)
            
        dlg.Destroy()
        
    def on_save_text(self, event):
        if self.fit != array([]):
            file_choices = "TXT-Datei (*.txt)|*.txt|DAT-Datei (*.dat)|*.dat|Alle Dateien (*.*)|*.*"
            filename = os.path.splitext(self.filename)[0] + "_fit.txt"
            dlg = wx.FileDialog(self, message="Ergebnis speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                
                try:
                    f = open(path, 'w')
                    
                    f.write('--------------------------------------------------------------------------------\n')
                    f.write('Starting model:\n\n')
                    self.sample.parameters = deepcopy(self.params)
                    f.write(self.sample.print_parameter())
                    
                    if self.sample.var_names != []:
                        f.write('\n--------------------------------------------------------------------------------\n')
                        f.write('Model resulting from fit:\n\n')
                        self.sample.parameters = deepcopy(self.params_new)
                        f.write(self.sample.print_parameter())
                        self.sample.parameters = deepcopy(self.params)
                        
                        f.write('\n--------------------------------------------------------------------------------\n')
                        try: 
                            f.write("Total fitting error: " + str((self.sample.err**2).sum()))
                        except: 
                            f.write("Total fitting Error: " + str((self.sample.residuals({}, fitalg="leastsq")**2).sum()))
                        f.write("\n\n")
                        for key in self.sample.var_names:
                            f.write("Fitting of %s: old value: %.6g, new value: %.6g +- %.6g\n" %(key, self.params[key], self.params_new[key], self.errors[key]))
                    
                    f.write('\n--------------------------------------------------------------------------------\n\n\n')
                    
                    f.write("Omega\tData\tFit\n")
                    savetxt(f, vstack((self.angle, self.int, self.fit)).T, fmt='%11.6g', delimiter='\t')
                    f.close()
                    self.flash_status_message("Ergebnis gespeichert unter %s" % path)
                except Exception as error:
                    wx.MessageBox('Fehler beim Speichern der Datei:\n' + path + '\n\nFehlermeldung:\n' + str(error), 'Fehler', style=wx.ICON_ERROR)
                
            dlg.Destroy()
            
        else:
            wx.MessageBox('Kein Fit vorhanden - Speichern nicht erfolgt.', 'Fehler', style=wx.ICON_ERROR)
    
    def flash_status_message(self, msg, flash_len_ms=5000):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_flash_status_off, self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')
        
    def on_UpdateCursor(self, event):
        if event.inaxes:
            if abs(event.xdata) > 1e-3 and abs(event.xdata) < 1e5:
                text1 = 'x = %.5f' %event.xdata
            else:   
                text1 = 'x = %.4e' %event.xdata
            if abs(event.ydata) > 1e-1 and abs(event.ydata) < 1e1:
                text2 = 'y = %.5f' %event.ydata
            else:   
                text2 = 'y = %.4e' %event.ydata
            self.statusbar.SetStatusText(text1, 1)
            self.statusbar.SetStatusText(text2, 2)
        else:
            self.statusbar.SetStatusText('', 1)
            self.statusbar.SetStatusText('', 2)
            
    def on_Resize(self, event):
        try:
            x, y = self.fig.get_size_inches()
            self.fig.subplots_adjust(left=0.7/x, right=1-0.2/x, bottom=0.45/y, top=1-0.35/y)
            self.titel.set_position((0.5, 1+0.05/y))
        except Exception as error:
            return
            
    def on_hide_log(self, event):
        if self.bm_log.GetLabel() == 'Log verbergen':
            self.log.Hide()
            self.bm_log.SetLabel('Log anzeigen')
        else:
            self.log.Show()
            self.bm_log.SetLabel('Log verbergen')
        self.hbox.Layout()
            
    def on_about(self, event):
        msg = """Anpassen von Reflektometrie-Daten:
        
        - Laden einer Datei mit Messdaten
        - Modell aus param-Datei laden
        - Modell per GUI anpassen
        - Parameter per GUI anpassen
        - Modell zeichnen
        - Fit starten und wiederholen
        - Daten und Grafik speichern
                  
        (basiert auf wxPython und matplotlib)
        
        Version 0.5 - 04.04.2013
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def on_exit(self, event):
        self.Destroy()


class SelectEnergy(wx.Dialog):
    """ Dialog to select an energy value. """

    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(260, 240))
        
        panel = wx.Panel(self, -1)
        vbox = wx.BoxSizer(wx.VERTICAL)

        wx.StaticBox(panel, -1, 'Energie (keV)', (5, 5), (240, 155))
        self.rb1 = wx.RadioButton(panel, -1, 'Cu_K_alpha_1+2: 8.04116', (15, 30), style=wx.RB_GROUP)
        self.rb2 = wx.RadioButton(panel, -1, 'Cu-K_alpha_1: 8.04782', (15, 55))
        self.rb3 = wx.RadioButton(panel, -1, 'Cu-K_alpha_2: 8.02784', (15, 80))
        self.rb4 = wx.RadioButton(panel, -1, 'Cu-K_beta_1,3: 8.90541', (15, 105))
        self.rb5 = wx.RadioButton(panel, -1, 'Benutzerdefiniert:', (15, 130))
        self.tc = wx.TextCtrl(panel, -1, '', (125, 127))

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, wx.ID_OK, 'OK', size=(90, 25))
        self.SetAffirmativeId(wx.ID_OK)

        closeButton = wx.Button(self, wx.ID_CANCEL, 'Abbrechen', size=(90, 25))
        self.SetEscapeId(wx.ID_CANCEL)

        hbox.Add(okButton, 1)
        hbox.Add(closeButton, 1, wx.LEFT, 5)

        vbox.Add(panel)
        vbox.Add(hbox, 1, wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, 10)
        self.SetSizer(vbox)
        
    def GetValue(self):
        if self.rb1.GetValue():
            return str(8.04116)
        elif self.rb2.GetValue():
            return str(8.04782)
        elif self.rb3.GetValue():
            return str(8.02784)
        elif self.rb4.GetValue():
            return str(8.90541)
        elif self.rb5.GetValue():
            return self.tc.GetValue()


class SelectRanges(wx.Dialog):
    """ Dialog to select ranges for brute force fitting. """

    def __init__(self, parent, id, title):
        self.num = len(app.frame.fit_keys)
        height = 150 + 32 * self.num
    
        wx.Dialog.__init__(self, parent, id, title, size=(400, height))
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        box = wx.StaticBox(self, -1, 'Brute Force: Bereiche bestimmen', (5, 5), (380, height-85))
        bsizer = wx.StaticBoxSizer(box, wx.VERTICAL)

        for i in arange(self.num):
            if app.frame.ranges.has_key(app.frame.fit_keys[i]):
                name = wx.StaticText(self, -1, app.frame.sample.names[app.frame.fit_keys[i]] + ':', size=(170,-1))
                t_start = wx.TextCtrl(self, 100+i, app.frame.ranges[app.frame.fit_keys[i]][0], style=wx.TE_PROCESS_ENTER, size=(70,20))    
                zwischen = wx.StaticText(self, -1, " bis ", size=(20,-1))
                t_end = wx.TextCtrl(self, 200+i, app.frame.ranges[app.frame.fit_keys[i]][1], style=wx.TE_PROCESS_ENTER, size=(70,20))        
            else:
                name = wx.StaticText(self, -1, app.frame.sample.names[app.frame.fit_keys[i]] + ':', size=(170,-1))
                t_start = wx.TextCtrl(self, 100+i, str(app.frame.params[app.frame.fit_keys[i]] * 0.9), style=wx.TE_PROCESS_ENTER, size=(70,20))    
                zwischen = wx.StaticText(self, -1, " bis ", size=(20,-1))
                t_end = wx.TextCtrl(self, 200+i, str(app.frame.params[app.frame.fit_keys[i]] * 1.1), style=wx.TE_PROCESS_ENTER, size=(70,20))
            
            hbox1 = wx.BoxSizer(wx.HORIZONTAL)
            hbox1.Add(name, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            hbox1.Add(t_start, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)            
            hbox1.Add(zwischen, 0, border=3, flag = wx.ALIGN_CENTER | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            hbox1.Add(t_end, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            
            bsizer.Add(hbox1, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)       
        
        bsizer.AddSpacer(5)
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        label_t_anz = wx.StaticText(self, -1, "Anzahl Schritte: ", size=(170,-1))
        self.t_anz = wx.TextCtrl(self, -1, str(app.frame.Ns), style=wx.TE_PROCESS_ENTER, size=(70,20))
        hbox2.Add(label_t_anz, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        hbox2.Add(self.t_anz, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL) 
        bsizer.Add(hbox2, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
            
        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, wx.ID_OK, 'OK', size=(100, 25))
        self.SetAffirmativeId(wx.ID_OK)
        closeButton = wx.Button(self, wx.ID_CANCEL, 'Abbrechen', size=(100, 25))
        self.SetEscapeId(wx.ID_CANCEL)
        hbox3.Add(okButton, 1)
        hbox3.Add(closeButton, 1, wx.LEFT, 5)

        vbox.Add(bsizer, 0, border=7, flag = wx.ALIGN_CENTER | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        vbox.Add(hbox3, 1, border=7, flag = wx.ALIGN_CENTER | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.SetSizer(vbox)
        
    def GetValue(self):
        fit_ranges = []
        for i in arange(self.num):
            start = self.FindWindowById(100+i).GetValue()
            end = self.FindWindowById(200+i).GetValue()
            app.frame.ranges[app.frame.fit_keys[i]] = [start, end]
            fit_ranges.append([float(start), float(end)])
        Ns = int(self.t_anz.GetValue())
        return fit_ranges, Ns


class RedirectText(object):
    """ Redirects STDOUT and STDERR to log window. """
    def __init__(self, aWxTextCtrl, color):
        self.out = aWxTextCtrl
        self.color = color
 
    def write(self, string):
        self.out.SetDefaultStyle(wx.TextAttr(self.color))
        self.out.WriteText(string)
        
        
class DensityPlot(wx.Dialog):
    """ Window to show model. """

    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(600, 475))
        
        panel = wx.Panel(self, -1)
        
        self.fig = Figure((6, 4), dpi=100)
        self.canvas = FigureCanvasWxAgg(panel, -1, self.fig)
        
        okButton = wx.Button(self, wx.ID_OK, 'OK', size=(90, 25))
        self.SetAffirmativeId(wx.ID_OK)

        closeButton = wx.Button(self, wx.ID_CANCEL, 'Abbrechen', size=(90, 25))
        self.SetEscapeId(wx.ID_CANCEL)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(okButton, 1)
        hbox.AddSpacer(10)
        hbox.Add(closeButton, 1)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(panel)
        vbox.AddSpacer(10)
        vbox.Add(hbox, 1, wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM)
        
        self.SetSizer(vbox)
        
    def plot_density(self, z, rho):
        self.fig.clear()
        self.axes = self.fig.add_subplot(111)
        self.axes.grid(1)
        
        self.plot = self.axes.plot(z, rho, 'r')
        
        self.axes.tick_params(axis='both', labelsize=10)
        self.axes.set_xlabel('Tiefe (nm)', fontsize=12)
        self.axes.set_ylabel('rho (g/cm^3)', fontsize=12)
        self.titel = self.axes.set_title('Modell der Dichte', fontsize=12)
        
        self.canvas.draw()
           
        
if __name__ == '__main__':
    app = wx.PySimpleApp()
    app.frame = MainFrame()
    app.frame.Show()
    app.MainLoop()