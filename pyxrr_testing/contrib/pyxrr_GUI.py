#!/usr/bin/env python
import wx
import pyxrr

# The recommended way to use wx with mpl is with the WXAgg backend. 
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas, NavigationToolbar2WxAgg as NavigationToolbar
    
from pylab import *
from re import findall
from StringIO import StringIO
from wx.lib.agw import ultimatelistctrl as ULC


class MyFrame(wx.Frame):
    """ The main frame of the application. """
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, 'Anpassen von Reflektometrie-Daten')
        
        # Attributes of class
        self.filename = ""
        self.path = ""
        self.tempfile = "Current_Model.param"
        self.sample = []
        self.angle = array([0, 0])
        self.int = array([1, 1])
        self.fit = array([])
        self.lastfit = array([])
        
        # Default Parameters for initialization
        self.params = {'LayerCount': [1, 2, 1], 'N': [1, 2, 1], 'energy0': 8.04116, 'resolution0': 0.01, 'offset0': 0.0, 'scale0': 1.0, 'background0': -6.0, 'rho_0': 0.0012, 'rho_1': 1.0, 'rho_2': 2.0, 'rho_3': 3.0, 'd_0': 0.0, 'd_1': 1.0, 'd_2': 2.0, 'd_3': 3.0, 'sigma_0': 0.0, 'sigma_1': 1.0, 'sigma_2': 2.0, 'sigma_3': 3.0, 'grad_d_0': 0.0, 'grad_d_1': 0.0, 'grad_d_2': 0.0}
        self.params_last = self.params.copy()
        self.params_new = self.params.copy()
        self.model = ['N.78O.21', 'TiO2', 'SiO2', 'Si']
        
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
            self.update_model()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_load = menu_file.Append(-1, "&Datei laden...\tCtrl-O", "Daten aus Datei laden")
        self.Bind(wx.EVT_MENU, self.on_open_file, m_load)
        # m_loadmodel = menu_file.Append(-1, "&Modell laden...\tCtrl-M", "Modell aus Datei laden")
        # self.Bind(wx.EVT_MENU, self.on_load_model, m_loadmodel)
        m_savemodel = menu_file.Append(-1, "&Modell speichern...\tCtrl-N", "Modell als Datei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_model, m_savemodel)
        menu_file.AppendSeparator()
        m_runfit = menu_file.Append(-1, "&Fit starten!\tCtrl-F", "Fit-Prozedur starten")
        self.Bind(wx.EVT_MENU, self.on_run_fit, m_runfit)     
        m_savetext = menu_file.Append(-1, "&Ergebnis speichern...\tCtrl-S", "Ergebnis als Datei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_text, m_savetext)
        m_saveplot = menu_file.Append(-1, "&Grafik speichern...\tCtrl-G", "Grafik als Datei speichern")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_saveplot)
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
        
        self.dpi = 100
        self.fig = Figure((7, 5), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)      
        self.toolbar = NavigationToolbar(self.canvas)
    
        self.b_load = wx.Button(self.panel, -1, "Datei laden...", size=(130,25))
        self.b_load.Bind(wx.EVT_BUTTON, self.on_open_file)
        
        self.b_loadmodel = wx.Button(self.panel, -1, "Modell laden...", size=(130,25))
        self.b_loadmodel.Bind(wx.EVT_BUTTON, self.on_open_model)
        
        self.b_savemodel = wx.Button(self.panel, -1, "Modell speichern...", size=(130,25))
        self.b_savemodel.Bind(wx.EVT_BUTTON, self.on_save_model)
        
        self.b_runfit = wx.Button(self.panel, -1, "Fit starten!", size=(130,35))
        self.b_runfit.Bind(wx.EVT_BUTTON, self.on_run_fit)
        
        self.b_savetext = wx.Button(self.panel, -1, "Ergebnis speichern...", size=(130,25))
        self.b_savetext.Bind(wx.EVT_BUTTON, self.on_save_text)
        
        self.b_saveplot = wx.Button(self.panel, -1, "Grafik speichern...", size=(130,25))
        self.b_saveplot.Bind(wx.EVT_BUTTON, self.on_save_plot)
              
        self.t_pollabel = wx.StaticText(self.panel, -1, "Polarisation:", size=(60,-1))
        self.cb_pol = wx.ComboBox(self.panel, -1, choices=['unpolarisiert', 'parallel', 'senkrecht'], style=wx.CB_READONLY, size=(102,-1))
        self.cb_pol.SetValue('unpolarisiert')
        
        self.t_anglelabel = wx.StaticText(self.panel, -1, "Winkel:", size=(60,-1))
        self.cb_angle = wx.ComboBox(self.panel, -1, choices=['omega', '2theta'], style=wx.CB_READONLY, size=(102,-1))
        self.cb_angle.SetValue('2theta')
        
        self.t_startlabel = wx.StaticText(self.panel, -1, "Start:")
        self.t_start = wx.TextCtrl(self.panel, -1, '0.0', style=wx.TE_PROCESS_ENTER, size=(45,20))    
        self.t_start.Bind(wx.EVT_TEXT, self.check_float)
        
        self.t_endlabel = wx.StaticText(self.panel, -1, " Ende:")
        self.t_end = wx.TextCtrl(self.panel, -1, '100.0', style=wx.TE_PROCESS_ENTER, size=(45,20))    
        self.t_end.Bind(wx.EVT_TEXT, self.check_float)
        
        self.lb_model = ULC.UltimateListCtrl(self.panel, -1, agwStyle=ULC.ULC_REPORT | ULC.ULC_HAS_VARIABLE_ROW_HEIGHT | ULC.ULC_SINGLE_SEL, size=(175,180))
        self.lb_model.InsertColumn(0, "Name", width=70)
        self.lb_model.InsertColumn(1, "Material", width=85)
        
        self.bm_new = wx.Button(self.panel, -1, "Neu", size=(45,25))
        self.bm_new.Bind(wx.EVT_BUTTON, self.on_model_new)
        
        self.bm_del = wx.Button(self.panel, -1, "Entf.", size=(45,25))
        self.bm_del.Bind(wx.EVT_BUTTON, self.on_model_del)
    
        # self.bm_up = wx.Button(self.panel, -1, "Hoch", size=(45,25))
        # self.bm_up.Bind(wx.EVT_BUTTON, self.on_model_up)
        
        # self.bm_down = wx.Button(self.panel, -1, "Runter", size=(45,25))
        # self.bm_down.Bind(wx.EVT_BUTTON, self.on_model_down)
        
        self.lb_table = ULC.UltimateListCtrl(self.panel, -1, agwStyle=ULC.ULC_REPORT | ULC.ULC_HAS_VARIABLE_ROW_HEIGHT | ULC.ULC_NO_HIGHLIGHT, size=(430,150))
        self.lb_table.InsertColumn(0, "Parameter (Fit?)", width=160)
        self.lb_table.InsertColumn(1, "Startwert", width=90)
        self.lb_table.InsertColumn(2, "Letzter Fit")
        self.lb_table.InsertColumn(3, "Neuer Fit")

       
        # Layout with box sizers ----------------------------------------------
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        
        self.vbox11 = wx.BoxSizer(wx.VERTICAL)
        self.vbox11.Add(self.b_load, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL)
        
        self.hbox111 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox111.Add(self.t_pollabel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox111.Add(self.cb_pol, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.hbox112 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox112.Add(self.t_anglelabel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox112.Add(self.cb_angle, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.hbox113 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox113.Add(self.t_startlabel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox113.Add(self.t_start, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox113.Add(self.t_endlabel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox113.Add(self.t_end, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)

        self.vbox11.Add(self.hbox111, 0, border=2, flag = wx.ALIGN_LEFT | wx.ALL)        
        self.vbox11.Add(self.hbox112, 0, border=2, flag = wx.ALIGN_LEFT | wx.ALL)
        self.vbox11.Add(self.hbox113, 0, border=2, flag = wx.ALIGN_LEFT | wx.ALL)    
        self.vbox11.AddSpacer(2)
        self.vbox11.Add(self.b_loadmodel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL)
        self.vbox11.Add(self.b_savemodel, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL)
        
        self.vbox12 = wx.BoxSizer(wx.VERTICAL)        
        self.vbox12.Add(self.bm_new, 0, border=3, flag = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
        self.vbox12.Add(self.bm_del, 0, border=3, flag = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
        # self.vbox12.Add(self.bm_up, 0, border=3, flag = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
        # self.vbox12.Add(self.bm_down, 0, border=3, flag = wx.ALIGN_LEFT | wx.BOTTOM | wx.TOP | wx.RIGHT | wx.ALIGN_CENTER_VERTICAL)
             
        self.hbox1.Add(self.vbox11, 0, border=5, flag = wx.ALIGN_LEFT | wx.ALL)
        self.hbox1.AddSpacer(5)
        self.hbox1.Add(self.lb_model, 0, border=5, flag = wx.ALIGN_LEFT | wx.ALL | wx.EXPAND)
        self.hbox1.Add(self.vbox12, 0, border=0, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)            
        self.hbox2.Add(self.b_runfit, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(self.b_savetext, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(self.b_saveplot, 0, border=3, flag = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.vbox1.Add(self.hbox1, 0, border=5, flag = wx.ALIGN_LEFT | wx.ALL | wx.EXPAND)
        self.vbox1.Add(self.lb_table, 1, border=5, flag = wx.ALIGN_LEFT | wx.ALL | wx.EXPAND)
        self.vbox1.Add(self.hbox2, 0, border=5, flag = wx.ALIGN_LEFT | wx.ALL | wx.EXPAND)
        
        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.vbox2.Add(self.canvas, 1, flag = wx.LEFT | wx.TOP | wx.EXPAND)
        self.vbox2.Add(self.toolbar, 0, flag = wx.EXPAND)
             
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox.Add(self.vbox1, 0, flag = wx.ALIGN_LEFT | wx.TOP | wx.EXPAND)
        self.hbox.Add(self.vbox2, 1, flag = wx.ALIGN_LEFT | wx.TOP | wx.EXPAND)
        
        self.panel.SetSizer(self.hbox)
        self.hbox.Fit(self)
        
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()
    
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
            
        self.update_table()
        
    def new_layer(self, row, label, name):
        row = self.lb_model.GetItemCount()
        self.lb_model.InsertStringItem(row, label)
        self.lb_model.SetItemData(row, label)
        if name != "Gruppe":
            text = wx.TextCtrl(self.lb_model, value=name, id=1000+row, style=wx.TE_PROCESS_ENTER, size=(83,20))
            text.Bind(wx.EVT_TEXT, self.on_enter_material)
            self.lb_model.SetItemWindow(row, 1, wnd=text)
        self.lb_model.Update()

    def update_table(self):
        self.group_layer = {}
        self.layer_sigma = {}
        self.group_sigma = {}

        for i in range(self.lb_table.GetItemCount()):
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
            self.new_param('Perioden in Gruppe ' + str(i+1), 'N' + str(i+1))
        
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
       
    def new_param(self, name, key):
        row = self.lb_table.GetItemCount()
        self.lb_table.InsertStringItem(row, ' ' + name, it_kind=1)
        self.lb_table.SetItemData(row, key)
        if 'N' in key:
            num = int(key.strip('N'))
            self.lb_table.SetStringItem(row, 2, str(self.params_last['N'][num]))
            self.lb_table.SetStringItem(row, 3, str(self.params_new['N'][num]))
            text = wx.TextCtrl(self.lb_table, value=str(self.params['N'][num]), id=2000+row, style=wx.TE_PROCESS_ENTER, size=(88,20))     
        else:
            self.lb_table.SetStringItem(row, 2, str(self.params_last[key]))
            self.lb_table.SetStringItem(row, 3, str(self.params_new[key]))
            text = wx.TextCtrl(self.lb_table, value=str(self.params[key]), id=2000+row, style=wx.TE_PROCESS_ENTER, size=(88,20))
        text.Bind(wx.EVT_TEXT, self.on_enter_value)
        self.lb_table.SetItemWindow(row, 1, wnd=text)
        button1 = wx.Button(self.lb_table, label="<", id=3000+row, size=(20,20))
        button1.Bind(wx.EVT_BUTTON, self.on_load_last)
        self.lb_table.SetItemWindow(row, 2, wnd=button1)
        button2 = wx.Button(self.lb_table, label="<", id=4000+row, size=(20,20))
        button2.Bind(wx.EVT_BUTTON, self.on_load_new)
        self.lb_table.SetItemWindow(row, 3, wnd=button2) 
        self.lb_table.Update()
        
    def on_select_energy(self, event):
        dlg = SelectEnergy(None, -1, 'Energie bestimmen')
        
        if dlg.ShowModal() == wx.ID_OK:
            value = dlg.GetValue()
            self.FindWindowById(2000).SetValue(value)

        dlg.Destroy()
   
    def draw_figure(self):
        """ Redraws the figure. """
        self.fig.clear()
        self.axes = self.fig.add_subplot(111)

        self.axes.semilogy(self.angle, self.int, 'b', label='Daten')
        if self.lastfit != array([]):
            self.axes.semilogy(self.angle, self.lastfit, 'g', label='Letzter Fit')
        if self.fit != array([]):
            self.axes.semilogy(self.angle, self.fit, 'r', label='Neuer Fit')
            
        prop = matplotlib.font_manager.FontProperties(size=10) 
        self.axes.legend(loc=0, prop=prop)
        self.axes.tick_params(axis='both', labelsize=10)
        self.axes.set_xlabel('Omega (deg)', fontsize=12)
        self.axes.set_ylabel('Intensitaet (normiert)', fontsize=12)
        self.fig.suptitle(self.filename, fontsize=12)
        if self.filename != '':
            self.axes.set_title('Darstellung von Daten und Fit', fontsize=10)
        
        self.canvas.draw()
 
    def on_open_file(self, event):
        dlg = wx.FileDialog(self, "Daten laden", style=wx.OPEN)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.path = dlg.GetPath()
            
            try:
                self.save_model(self.tempfile)
                self.sample = pyxrr.multilayer(self.tempfile)           
                self.angle = self.sample.measured_data[0][:,0]
                self.int = self.sample.measured_data[0][:,1]
                self.fit = array([])
                self.lastfit = array([])
                                        
                self.flash_status_message("%s geladen." % self.filename)
                self.draw_figure()
            except Exception as error:
                wx.MessageBox('Fehler beim Laden der Datei:\n' + self.path + '\n(Pfad darf weder Umlaute noch Leerzeichen enthalten!)', str(error))

        dlg.Destroy()
        
    def on_open_model(self, event):
        file_choices = "PARAM (*.param)|*.param"
        dlg = wx.FileDialog(self, "Modell laden", style=wx.OPEN, wildcard=file_choices)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.filename = dlg.GetFilename()
            
            try:
                self.open_model(path)
            except Exception as error:
                wx.MessageBox('Fehler beim Laden der Datei:\n' + path + '\n(Pfad darf weder Umlaute noch Leerzeichen enthalten!)', str(error))

        dlg.Destroy()
        
    def open_model(self, openpath):
        changed_tt = 0
        changed_pol = 0
    
        f = open(openpath, 'r')
        for line in f:
            if "Measurement:" in line:
                for prop in line.replace("Measurement:", "").replace("\n","").replace(" ","").split(","):
                    if prop.find("=") > 0:
                        thisprop = prop.split("=")
                        if thisprop[0] == 'file':
                            self.path = thisprop[1]
                        elif thisprop[0] == 'fit_range':
                            start, end = thisprop[1].split("->")
                            self.t_start.SetValue(start)
                            self.t_end.SetValue(end)  
                        elif thisprop[0] == 'twotheta':
                            changed_tt = 1
                            if bool(int(thisprop[1])):
                                self.cb_angle.SetValue('2theta')
                            else:
                                self.cb_angle.SetValue('omega')
                        elif thisprop[0] == 'pol':
                            changed_pol = 1
                            if float(thisprop[1]) == 0.0:
                                self.cb_pol.SetValue('senkrecht')
                            elif float(thisprop[1]) == 1.0:
                                self.cb_pol.SetValue('parallel')
                            elif float(thisprop[1]) == 0.5:
                                self.cb_pol.SetValue('unpolarisiert')                
        f.close()
        
        if changed_tt == 0:
            self.cb_angle.SetValue('2theta')
        if changed_pol == 0:
            self.cb_pol.SetValue('unpolarisiert')
                    
        self.sample = pyxrr.multilayer(openpath)
        self.params = self.sample.parameters
        self.params_last = self.params.copy()
        self.params_new = self.params.copy()
        self.model = list(self.sample.materials)
        self.angle = self.sample.measured_data[0][:,0]
        self.int = self.sample.measured_data[0][:,1]
        self.fit = array([])
        self.lastfit = array([])
                                
        self.flash_status_message("%s geladen." % self.filename)
        self.update_model()
        self.draw_figure()

    def on_save_model(self, event):
        file_choices = "PARAM (*.param)|*.param"
        filename = self.filename.rstrip(".txt").rstrip(".param") + "_model.param"
        dlg = wx.FileDialog(self, message="Modell speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            
            try:
                self.save_model(path)
            except Exception as error:
                wx.MessageBox('Fehler beim Speichern der Datei:\n' + path + '\n(Pfad darf weder Umlaute noch Leerzeichen enthalten!)', str(error))
            
        dlg.Destroy()
            
    def save_model(self, savepath):
        f = open(savepath, 'w')      
        f.write("Ambience: name=Umgebung, code=" + self.model[0] + ", rho=" + str(self.params['rho_0']) + "\n")
        
        layer_pos = 0
        sigma_pos = 0
        for i in range(len(self.params['LayerCount']) - 2):
            f.write("Group: name=Gruppe_" + str(i+1) + ", sigma=" + str(self.params['sigma_' + str(sigma_pos)]) + ", periods=" + str(self.params['N'][i+1]) + ", grad_d=0\n")
            sigma_pos += 1
            multilayer_sigma = 0
            for j in range(self.params['LayerCount'][i+1]):
                layer_pos += 1
                if j == 0:
                    if self.params['N'][i+1] == 1:
                        f.write("Layer: name=Schicht_" + str(layer_pos) + ", code=" + self.model[layer_pos] + ", rho=" + str(self.params['rho_' + str(layer_pos)]) + ", d=" + str(self.params['d_' + str(layer_pos)]) + "\n")
                    else:
                        special_sigma = sigma_pos + self.params['LayerCount'][i+1] - 1
                        multilayer_sigma = 1
                        f.write("Layer: name=Schicht_" + str(layer_pos) + ", code=" + self.model[layer_pos] + ", rho=" + str(self.params['rho_' + str(layer_pos)]) + ", d=" + str(self.params['d_' + str(layer_pos)]) + ", sigma=" + str(self.params['sigma_' + str(special_sigma)]) + "\n")
                else:    
                    f.write("Layer: name=Schicht_" + str(layer_pos) + ", code=" + self.model[layer_pos] + ", rho=" + str(self.params['rho_' + str(layer_pos)]) + ", d=" + str(self.params['d_' + str(layer_pos)]) + ", sigma=" + str(self.params['sigma_' + str(sigma_pos)]) + "\n")
                    sigma_pos += 1
            sigma_pos += multilayer_sigma
        
        value = self.cb_pol.GetValue()
        if value == 'senkrecht':
            pol = '0'
        elif value == 'parallel':
            pol = '1'
        else:
            pol = '0.5'
        
        f.write("Substrate: name=Substrat, code=" + self.model[-1] + ", rho=" + str(self.params['rho_' + str(layer_pos+1)]) + ", sigma=" + str(self.params['sigma_' + str(sigma_pos)]) + "\n")          
        f.write("Measurement: file=" + self.path + ", twotheta=" + str(int(self.cb_angle.GetValue() == '2theta')) + ", fit_range=" + self.t_start.GetValue() + "->" + self.t_end.GetValue() + ", energy=" + str(self.params['energy0']) + ", resolution=" + str(self.params['resolution0']) + ", offset=" + str(self.params['offset0']) + ", scale=" + str(self.params['scale0']) + ", background=" + str(self.params['background0']) + ", pol=" + pol + "\n")
        f.close()
        self.flash_status_message("Modell gespeichert unter %s" % savepath)

    # def on_rb_omega(self, event):
        # self.draw_figure()
        
    # def on_rb_2theta(self, event):
        # self.draw_figure()
        
    # def on_cb_angle(self, event):
        # self.draw_figure()
        
    def check_float(self, event):
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
            
    def update_last_new(self):  
        for i in range(self.lb_table.GetItemCount()):
            key = self.lb_table.GetItemData(i)
            if 'N' in key:
                num = int(key.strip('N'))
                self.lb_table.SetStringItem(i, 2, str(self.params_last['N'][num]))
                self.lb_table.SetStringItem(i, 3, str(self.params_new['N'][num]))
            else:
                self.lb_table.SetStringItem(i, 2, str(self.params_last[key]))
                self.lb_table.SetStringItem(i, 3, str(self.params_new[key]))
        
    def on_run_fit(self, event):
        self.lastfit = self.fit
        self.params_last = self.params_new.copy()
    
        self.save_model(self.tempfile)
        sample = pyxrr.multilayer(self.tempfile)
        self.angle = sample.measured_data[0][:,0]
        self.int = sample.measured_data[0][:,1]
        
        fit_keys = []       
        for i in range(self.lb_table.GetItemCount()):
            if self.lb_table.GetItem(i, col=0).IsChecked():
                key = self.lb_table.GetItemData(i)
                if not 'N' in key:
                    fit_keys.append(key)
                              
        if fit_keys == []:
            self.params_new = self.params.copy()
        else:
            self.params_new = sample.fit_ref(fit_keys)
            
        self.fit = sample.reflectogram(self.angle, 0)
                      
        self.update_last_new()
        self.draw_figure()
        
        if fit_keys == []:
            self.flash_status_message("Zeichnen des Modells erfolgreich!")
        elif len(fit_keys) == 1:
            self.flash_status_message("Fitten eines Parameters erfolgreich!")
        else:
            self.flash_status_message("Fitten mit %.d Parametern erfolgreich!" %len(fit_keys))
        
    def on_model_new(self, event):
        index = self.lb_model.GetFirstSelected()
       
        if index < self.lb_model.GetItemCount() - 1:
            if index > 0:
                label = self.lb_model.GetItemData(index)
                num = int(label.split()[-1])
            else:
                label = ""
                
            if "Gruppe" in label or index == -1:
                self.insert_group()
            elif "Schicht" in label:
                self.insert_layer(num)
                
            self.params_last = self.params.copy()
            self.params_new = self.params.copy()
            self.update_model()
            
    def insert_layer(self, num):
        max = sum(self.params['LayerCount'])
        for i in reversed(range(num, max)):
            self.params['rho_' + str(i+1)] = self.params['rho_' + str(i)]
        self.params['rho_' + str(num)] = 1.0

        max = sum(self.params['LayerCount'])
        for i in reversed(range(num, max)):
            self.params['d_' + str(i+1)] = self.params['d_' + str(i)]
        self.params['d_' + str(num)] = 1.0

        for i in self.group_layer:
            for j in range(len(self.group_layer[i])):
                if num == self.group_layer[i][j]:
                    max = len(self.layer_sigma) + len(self.group_sigma) + 1
                    if j == 0 and self.params['N'][i] > 1:
                        save = self.params['sigma_' + str(self.layer_sigma[num])]
                        for k in reversed(range(self.group_sigma[i] + 1, max)):
                            self.params['sigma_' + str(k+1)] = self.params['sigma_' + str(k)]
                        self.params['sigma_' + str(self.group_sigma[i] + 1)] = save
                        self.params['sigma_' + str(self.layer_sigma[num] + 1)] = 0.0                        
                    else:
                        for k in reversed(range(self.layer_sigma[num], max)):
                            self.params['sigma_' + str(k+1)] = self.params['sigma_' + str(k)]
                        self.params['sigma_' + str(self.layer_sigma[num])] = 0.0
                    self.params['LayerCount'][i] += 1
                    break
       
        self.model.insert(num, '')
        
    def insert_group(self):
        num = sum(self.params['LayerCount']) - 1
        end = len(self.params['LayerCount']) - 1
    
        max = sum(self.params['LayerCount'])
        for i in reversed(range(num, max)):
            self.params['rho_' + str(i+1)] = self.params['rho_' + str(i)]
        self.params['rho_' + str(num)] = 1.0

        max = sum(self.params['LayerCount'])
        for i in reversed(range(num, max)):
            self.params['d_' + str(i+1)] = self.params['d_' + str(i)]
        self.params['d_' + str(num)] = 1.0
        
        max = len(self.layer_sigma) + len(self.group_sigma) + 1
        for i in range(max-1, max):
            self.params['sigma_' + str(i+1)] = self.params['sigma_' + str(i)]
        self.params['sigma_' + str(max-1)] = 0.0
       
        self.model.insert(num, '')
        self.params['LayerCount'].insert(end, 1)
        self.params['N'].insert(end, 1)

    def on_model_del(self, event):
        index = self.lb_model.GetFirstSelected()
            
        if index > 0 and index < self.lb_model.GetItemCount() - 1:
            label = self.lb_model.GetItemData(index)
            num = int(label.split()[-1])
            if "Gruppe" in label:
                self.remove_group(num)
            elif "Schicht" in label:
                for i in self.group_layer:
                    if num in self.group_layer[i] and self.params['LayerCount'][i] > 1:
                        self.remove_layer(num)
                
            self.params_last = self.params.copy()
            self.params_new = self.params.copy()
            self.update_model()

    def remove_layer(self, num):
        max = sum(self.params['LayerCount'])
        for i in range(num, max-1):
            self.params['rho_' + str(i)] = self.params['rho_' + str(i+1)]
        del self.params['rho_' + str(max-1)]

        max = sum(self.params['LayerCount'])
        for i in range(num, max-1):
            self.params['d_' + str(i)] = self.params['d_' + str(i+1)]
        del self.params['d_' + str(max-1)]   
        
        for i in self.group_layer:
            for j in range(len(self.group_layer[i])):
                if num == self.group_layer[i][j]:
                    max = len(self.layer_sigma) + len(self.group_sigma) + 1
                    if j == 0 and self.params['N'][i] > 1:
                        save = self.params['sigma_' + str(self.layer_sigma[num+1])]
                        for k in range(self.layer_sigma[num+1], max-1):
                            self.params['sigma_' + str(k)] = self.params['sigma_' + str(k+1)]
                        self.params['sigma_' + str(self.layer_sigma[num] - 1)] = save                    
                    else:
                        for k in range(self.layer_sigma[num], max-1):
                            self.params['sigma_' + str(k)] = self.params['sigma_' + str(k+1)]
                    del self.params['sigma_' + str(max-1)]
                    self.params['LayerCount'][i] -= 1
                    del self.group_layer[i][j]
                    break
        
        del self.model[num]
        del self.layer_sigma[num]
        
    def remove_group(self, groupnum):
        end_sigma = 0
        for num in reversed(self.group_layer[groupnum]):
            del self.model[num]
            del self.params['rho_' + str(num)]
            del self.params['d_' + str(num)]
            del self.params['sigma_' + str(self.layer_sigma[num])]
            end_sigma = max(end_sigma, self.layer_sigma[num])
            
        start = min(self.group_layer[groupnum])
        end = max(self.group_layer[groupnum])
        diff = end - start + 1
        
        maxi = sum(self.params['LayerCount'])
        for i in range(end+1, maxi):
            self.params['rho_' + str(i-diff)] = self.params['rho_' + str(i)]
            del self.params['rho_' + str(i)]
            
        maxi = sum(self.params['LayerCount'])
        for i in range(end+1, maxi):
            self.params['d_' + str(i-diff)] = self.params['d_' + str(i)]
            del self.params['d_' + str(i)]
            
        maxi = len(self.layer_sigma) + len(self.group_sigma) + 1
        for i in range(end_sigma+1, maxi):
            self.params['sigma_' + str(i-diff)] = self.params['sigma_' + str(i)]
            del self.params['sigma_' + str(i)]
            
        if self.params['N'][groupnum] > 1:
            maxi = len(self.layer_sigma) + len(self.group_sigma) + 1 - diff
            for i in range(self.group_sigma[groupnum], maxi-1):
                self.params['sigma_' + str(i)] = self.params['sigma_' + str(i+1)]
            del self.params['sigma_' + str(maxi-1)]
            del self.group_sigma[groupnum]   
            
        del self.params['LayerCount'][num]
        del self.params['N'][num]
        
    def on_model_up(self, event):
        index = self.lb_model.GetFirstSelected()
        if index >= 2:
            self.layers.insert(index-2, self.layers[index-1])
            self.layers_last.insert(index-2, self.layers_last[index-1])
            self.layers_new.insert(index-2, self.layers_new[index-1])
            del self.layers[index]
            del self.layers_last[index]
            del self.layers_new[index]
            self.update_lb()
            
    def on_model_down(self, event):
        index = self.lb_model.GetFirstSelected()
        if index >= 1 and index <= self.lb_model.GetItemCount()-2:
            self.layers.insert(index+1, self.layers[index-1])
            self.layers_last.insert(index+1, self.layers_last[index-1])
            self.layers_new.insert(index+1, self.layers_new[index-1])
            del self.layers[index-1]
            del self.layers_last[index-1]
            del self.layers_new[index-1]
            self.update_lb()
  
    def on_load_last(self, event):
        id = event.GetId()
        value = self.lb_table.GetItem(id-3000, 2).GetText()
        self.FindWindowById(id-1000).SetValue(value)
        
    def on_load_new(self, event):
        id = event.GetId()
        value = self.lb_table.GetItem(id-4000, 3).GetText()
        self.FindWindowById(id-2000).SetValue(value)      
        
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
        
    def on_enter_value(self, event):
        id = event.GetId()
        value = self.FindWindowById(id).GetValue()
        key = self.lb_table.GetItemData(id-2000)
        if 'N' in key:
            num = int(key.strip('N'))
            try:
                new_value = int(float(value))
                old_value = self.params['N'][num]
                self.params['N'][num] = new_value
                if new_value != float(value):
                    self.FindWindowById(id).SetValue(str(self.params['N'][num]))
                if new_value == 1 and old_value != 1:
                    self.remove_sigma(num)
                elif new_value != 1 and old_value == 1:
                    self.insert_sigma(num)
            except ValueError:
                self.FindWindowById(id).SetValue(str(self.params['N'][num]))
        else:
            try:
                self.params[key] = float(value)
            except:
                self.FindWindowById(id).SetValue(str(self.params[key]))
        
    def remove_sigma(self, groupnum):
        end_sigma = 0
        for num in self.group_layer[groupnum]:
            end_sigma = max(end_sigma, self.layer_sigma[num])
    
        maxi = len(self.layer_sigma) + len(self.group_sigma) + 1
        for i in range(end_sigma+1, maxi-1):
            self.params['sigma_' + str(i)] = self.params['sigma_' + str(i+1)]
            self.params_last['sigma_' + str(i)] = self.params_last['sigma_' + str(i+1)]
            self.params_new['sigma_' + str(i)] = self.params_new['sigma_' + str(i+1)]
        del self.params['sigma_' + str(maxi-1)]
        del self.params_last['sigma_' + str(maxi-1)]
        del self.params_new['sigma_' + str(maxi-1)]
    
        self.params['sigma_' + str(end_sigma)] = self.params['sigma_' + str(self.group_sigma[groupnum])]
        self.params_last['sigma_' + str(end_sigma)] = self.params_last['sigma_' + str(self.group_sigma[groupnum])]
        self.params_new['sigma_' + str(end_sigma)] = self.params_new['sigma_' + str(self.group_sigma[groupnum])]
        self.update_table()
    
    def insert_sigma(self, groupnum):
        end_sigma = 0
        for num in self.group_layer[groupnum]:
            end_sigma = max(end_sigma, self.layer_sigma[num])
            
        maxi = len(self.layer_sigma) + len(self.group_sigma) + 1
        for i in range(end_sigma+1, maxi):
            self.params['sigma_' + str(i+1)] = self.params['sigma_' + str(i)]
            self.params_last['sigma_' + str(i+1)] = self.params_last['sigma_' + str(i)]
            self.params_new['sigma_' + str(i+1)] = self.params_new['sigma_' + str(i)]
            
        self.params['sigma_' + str(end_sigma+1)] = 0.0
        self.params_last['sigma_' + str(end_sigma+1)] = 0.0
        self.params_new['sigma_' + str(end_sigma+1)] = 0.0
        self.update_table()
        
    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        filename = self.filename.rstrip(".txt").rstrip(".param") + "_fit.png"
        dlg = wx.FileDialog(self, message="Grafik speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            
            try:
                self.canvas.print_figure(path, dpi=3*self.dpi)
                self.flash_status_message("Grafik gespeichert unter %s" % path)
            except Exception as error:
                wx.MessageBox('Fehler beim Speichern der Datei:\n' + path + '\n(Pfad darf weder Umlaute noch Leerzeichen enthalten!)', str(error))
            
        dlg.Destroy()
        
    def on_save_text(self, event):
        if self.fit != array([]):
            file_choices = "TXT (*.txt)|*.txt"
            filename = self.filename.rstrip(".txt").rstrip(".param") + "_fit.txt"
            dlg = wx.FileDialog(self, message="Ergebnis speichern unter...", defaultFile=filename, wildcard=file_choices, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                
                try:
                    f = open(path, 'w')
                    
                    f.write('--------------------------------------------------------------------------------\n')
                    f.write('Starting model:\n\n')
                    f.write(self.sample.print_parameter())
                    
                    f.write('\n--------------------------------------------------------------------------------\n')
                    f.write('Model resulting from fit:\n\n')
                    self.sample.parameters = self.params_new
                    f.write(self.sample.print_parameter())
                    
                    # try: 
                        # f.write("\nError: " + str(sum(self.sample.err**2)))
                    # except: 
                        # f.write("\nError: " + str(sum(self.sample.residuals({})**2)))
                    f.write('--------------------------------------------------------------------------------\n\n\n')
                        
                    f.write("Omega\tData\tFit\n")
                    savetxt(f, vstack((self.angle, self.int, self.fit)).T, fmt='%.9f', delimiter='\t')
                    f.close()
                    self.flash_status_message("Ergebnis gespeichert unter %s" % path)
                except Exception as error:
                    wx.MessageBox('Fehler beim Speichern der Datei:\n' + path + '\n(Pfad darf weder Umlaute noch Leerzeichen enthalten!)', str(error))
                
            dlg.Destroy()
        
    def on_about(self, event):
        msg = """Anpassen von Reflektometrie-Daten:
        
        - Laden einer txt-Datei mit Daten
        - Modell per GUI anpassen
        - Parameter per GUI anpassen
        - Fit starten und wiederholen
        - Daten und Grafik speichern
                  
        (basiert auf wxPython und matplotlib)
        
        Version 0.2b2 - 23.08.2011
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
    
    def flash_status_message(self, msg, flash_len_ms=5000):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_flash_status_off, self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')

    def on_exit(self, event):
        self.Destroy()

class SelectEnergy(wx.Dialog):
    """ Dialog to select an energy value. """

    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title, size=(260, 240))
        
        self.value = 0.0

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
        
if __name__ == '__main__':
    app = wx.PySimpleApp()
    app.frame = MyFrame()
    app.frame.Show()
    app.MainLoop()