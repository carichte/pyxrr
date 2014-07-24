#!/usr/bin/env python

# Copyright 2013 by Carsten Richter
# Contact: carsten.richter@desy.de and
#          carsten.richter@physik.tu-freiberg.de
#
# This file is part of pyxrr.
#
# pyxrr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyxrr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyxrr inside the 'copying.txt' file.
# If not, see <http://www.gnu.org/licenses/>.


import os
os.environ["OPENBLAS_MAIN_FREE"] = '1'
import sys
import warnings
import pickle
import pylab as pl
import curses
import curses.textpad
import curses.ascii
import curses.panel
import signal


sys.path.insert(0, os.path.abspath(os.pardir)) # Path to pyxrr if it has not been installed as package
try:
    import pyxrr
    print "Using locally compiled version of xrr..."
except:
    print "Searching pyxrr in site packages..."
    sys.path.pop(0) # not there
    import pyxrr


warnings.filterwarnings("ignore")


"""
    This is a simple batch file for using pyxrr to fit measured data 
    interactively. It can be seen as an example
    
    see pyxrr.multilayer.__init__.__doc__ for an explanation.
"""

# database can be BrennanCowan, Chantler, CromerLiberman, EPDL97, Henke, Sasaki, Windt
# For more details see Databases on http://www.esrf.eu/computing/scientific/dabax
# and http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/






################################################################################



print("This is pyxrr version %s"%pyxrr.__version__)
print("Package location: %s"%os.path.abspath(pyxrr.__file__))

SEP = os.path.sep
lsep = os.linesep

if os.name == "posix": clcomm = 'clear'# Unix/Linux/MacOS/BSD/etc
elif os.name in ("nt", "dos", "ce"): clcomm='cls' # DOS/Windows
else: clcomm = ""

for directory in ["samples", "results", "tmp"]:
    if not os.path.isdir(directory):
        os.makedirs(directory)

print(lsep)
print("Choose sample file:")
for fname in sorted(os.listdir("samples")):
    if fname.endswith(".param"):
        print " %s"%fname

SAMPLENAME = raw_input(lsep + "Enter file name: samples" + SEP)
ROOT = SAMPLENAME.rpartition(".")[0]
if not ROOT:
    ROOT = SAMPLENAME

while True:
    try:
        sample=pyxrr.multilayer(os.path.join("samples",SAMPLENAME))
        break
    except pyxrr.pyxrrError as perr:
        if ("Density of material"   in str(perr.errmsg)) and \
           ("not found in database" in str(perr.errmsg)):
            print perr
            material = str(perr.errmsg).split("'")[1]
            rho = float(raw_input("Please enter density for %s: "%material))
            pyxrr.store_rho_to_db(material, rho, database = DATABASE_file)
        else: 
            raise(perr)

pl.ion()
fig, axes = pl.subplots(sample.number_of_measurements, num=SAMPLENAME)
axes = pl.array(axes, ndmin=1)
theta=[]
measured_plot=[]
initial_plot=[]
fitted_plot=[]
refl=[]
newy={}
xlabels = dict({'theta':'omega (deg)',
                'qz_nm':'q_z (nm)',
                'twotheta':'2theta (deg)',
                'qz_a':'q_z (A)'})

for i_M, ax in enumerate(axes):
    if len(sample.measured_data[i_M][:,0])>0:
        theta.append(sample.measured_data[i_M][:,0])
        measured_plot.append(ax.semilogy(theta[i_M], 
                                         sample.measured_data[i_M][:,1],
                                         c="k")[0])
    else:
        print("Only plotting supported for measurement %i"%i_M)
        while 1:
            try:
                thmin = float(raw_input("Start from %s = "%xlabels[sample.x_axes[i_M]]))
                thmax = float(raw_input("Stop at %s = "%xlabels[sample.x_axes[i_M]]))
                thnum = int(raw_input("Number of points = "))
                theta.append(pyxrr.np.linspace(thmin, thmax, thnum))
                break
            except KeyboardInterrupt:
                theta.append(sample.measured_data[i_M][:,0])
                break
            except Exception as errmsg:
                print("%sWrong input!%sMessage: %s" %(lsep,lsep,errmsg))
                continue
        measured_plot.append(
            ax.semilogy(theta[i_M],
                        pl.ones_like(theta[i_M]), 
                        visible=False)[0])
    
    refl = sample.reflectogram(theta[i_M], i_M)
    initial_plot.append(ax.semilogy(theta[i_M], refl, alpha=0.5, c="g")[0])
    fitted_plot.append( ax.semilogy(theta[i_M], refl, c="r")[0])
    ax.set_xlabel(xlabels[sample.x_axes[i_M]])
    ax.set_ylabel("Reflectivity")
    ax.grid(True)

axes[0].legend([measured_plot[0], fitted_plot[0], initial_plot[0]],
               ["measurement", "fit", "previous"])


fitoptions = dict(n = "leastsq refinement (Levenberg-Marquardt algorithm)", 
               nf = "fmin_bfgs (BFGS algorithm)",
               ns = "fmin (downhill simplex algorithm)",
               nb = "brute (brute force)",
               np = "fmin_powell (Powell's algorithm)",
               ncg = "fmin_cg (nonlinear conjugate gradient algorithm)")
options = dict(e = "edit parameter", 
               a = "add parameter", 
               q = "quit program",
               c = "set coupled parameters",
               t = "enter equation for thickness dependency",
               u = "undo last changes",
               s = "save parameters",
               l = "load parameters",
               w = "write ascii",
               d = "write density vs. depth profile",
               r = "refresh plot",
               p = "print parameters",
               h = "this message",
               setup = "set simulation options",
               help = "this message")
options["?"] = options["h"]
info = ["Fit options:"]
info += ["  %s - %s"%(k,fitoptions[k]) for k in sorted(fitoptions.keys())]
info += ["%s - %s"%(k,options[k]) for k in sorted(options.keys())]




def get_filename(directory="", default="", extension=""):
    try:
        fname = raw_input("  Filename (without extension)? [%s] "%default)
        if not fname:
            fname = default
        fpath = os.path.join(directory, fname + extension)
        if os.path.isfile(fpath):
            overwrite = raw_input("  File %s exists. Overwrite? (y/[n])"%fpath).lower() == "y"
        else:
            overwrite = True
    except:
        return
    if overwrite:
        return fpath


def reprint():
    pass

def validator(ch):
    if ch==curses.ascii.NL:
        return curses.ascii.BEL
    if not curses.ascii.isprint(ch):
        reprint()
    #print ch
    return ch


class Screen:
    def __init__(self, welcome=""):
        self.S = curses.initscr()
        curses.cbreak()
        curses.noecho()
        self.S.keypad(1)
        self.offsety = 2
        self._y, self._x = self.S.getmaxyx()
        if welcome:
            self.printlines(welcome)
            self.S.addstr(self._y-2, 2, "Hit any key...", curses.A_BOLD)
            self.S.border(0)
            self.S.getch()
    def __enter__(self):
        return self
    def __exit__(self,a,b,c):
        curses.nocbreak()
        self.S.keypad(0)
        curses.echo()
        curses.endwin()
    def printlines(self, lines, offsetx=2, offsety=None, fmt=0):
        if isinstance(lines, str):
            lines = [lines]
        if offsety==None:
            offsety = self.offsety
        line = 0
        width = self._x - offsetx
        for i, msg in enumerate(lines):
            #self.S.addnstr(offsety + line, offsetx - 1, "-", 1)
            if not ((offsety + line) < self._y) or not ((offsetx+1) < self._x):
                break
            while msg:
                try:
                    self.S.addnstr(offsety + line, offsetx, msg, width, fmt)
                except:
                    break
                msg = msg[width:]
                line+=1
        self.offsety = offsety + line
        self.S.refresh()
        return len(lines) - (i+1)
    def prompt(self, question, y=2, x=2, dy=3, dx=None, options=()):
        if dx==None:
            dx = self._x - x - 2
        def reprint():
            maxlen=0
            for option in options:
                self.printlines(option, x+maxlen+4, y+dy+3)
                maxlen += max(map(len, option))
            self.S.addstr(y, x, question)
            curses.textpad.rectangle(self.S, y+1, x, y+dy+2, x+dx)
        reprint()
        self._reprint = reprint
        self.S.refresh()
        self.win = curses.newwin(dy, dx-1, y+2, x+1)
        self.tb = curses.textpad.Textbox(self.win)
        text = self.tb.edit(validator)
        del self.tb
        del self._reprint
        #inp = self.S.getstr(y+dy, x+dx)
        return text
    
    def menu(self, items, position=0, info={}, showall=False, selected=None,
                   readonly = []):
        items.append('exit')
        if position in items:
            position = items.index(position)
        elif not isinstance(position, int):
            position = 0
        self.win = self.S.subwin(0,0)
        self.win.keypad(1)
        self.panel = curses.panel.new_panel(self.win)
        #self.panel.hide()
        curses.panel.update_panels()
        self.panel.top()
        self.panel.show()
        self.win.clear()
        values = list(items)
        items = map(str, items)
        maxlen = max(map(len, items))+1
        space = 5
        infolen = self._x - 1 - maxlen - space
        def navigate(n, pos):
            pos += n
            if pos < 0:
                pos = 0
            elif pos >= len(items):
                pos = len(items)-1
            return pos
        
        digits = [ord(str(i)) for i in xrange(min(len(items), 9))]

        while True:
            self.win.clear()
            self.win.refresh()
            curses.doupdate()
            if selected!=None:
                self.win.addstr(0, 0, 
                    "Use <space> to mark item/variable and <enter> to select",
                    curses.A_UNDERLINE)
            start = max(0, position - (self._y-7))
            #print start
            for index, item in enumerate(items):
                if index == position:
                    mode = curses.A_REVERSE
                else:
                    mode = curses.A_NORMAL 
                if index+1-start>=self._y-4:
                    break
                if index<start:
                    continue
                msg = '%d. %s' % (index, item)
                self.win.addstr(1+index-start, 1, msg, mode)
                if selected!=None and item in selected:
                    self.win.addstr(1+index-start, 0, "*", curses.A_BOLD)
                if item in info and showall:
                    iteminfo = info[item]
                    self.win.addstr(1+index-start, 1+space+maxlen, 
                                    iteminfo[:infolen], mode)
            if items[position] in info:
                iteminfo = info[items[position]]
                self.win.addstr(self._y-3, 0, iteminfo, curses.A_BOLD)

            key = self.win.getch()
            if key == curses.KEY_UP:
                position = navigate(-1, position)

            elif key == curses.KEY_DOWN:
                position = navigate( 1, position)
            elif key == curses.KEY_END:
                position = len(items)-1
            elif key == curses.KEY_HOME:
                position = 0
            elif key in digits:
                position = digits.index(key)
            elif values[position] in readonly:
                continue
            elif key in [curses.KEY_ENTER, ord('\n')]:
                break
            elif key in [ord(' ')] and selected!=None:
                if values[position]=="exit":
                    continue
                elif values[position] in selected:
                    selected.pop(selected.index(values[position]))
                else:
                    selected.append(values[position])
        
        self.win.clear()
        self.panel.hide()
        curses.panel.update_panels()
        curses.doupdate()
        return values[position]
    
    def draw_title(self, title=None):
        if title==None and hasattr(self, "title"):
            title = self.title
        elif title==None:
            return
        else:
            self.title = title
        self.S.addstr(0,0, title, curses.A_BOLD + curses.A_UNDERLINE)
    
    def resize(self, signum=None, frame=None):
        self.clear()
        curses.endwin()
        self.S = curses.initscr()
        self._y, self._x = self.S.getmaxyx()
        #curses.resizeterm(self._y, self._x)
        #self.S.resize(self._y, self._x)
        #self.S.addstr(9,5, str(resize))
        winsize = "rows: %i, cols: %i"%(self._y, self._x)
        self.S.addstr(self._y-1, self._x-len(winsize)-1, winsize)
        #self.S.addstr(0, self._x-len(winsize)-1, winsize)
        if hasattr(self, "_reprint"):
            self._reprint()
        if hasattr(self, "tb"):
            text = self.tb.gather()
            [self.tb.do_command(curses.ascii.BS) for ch in text]
            text = text.strip()
            [self.tb.do_command(ch) for ch in text]
        self.draw_title()
        self.S.refresh()
    def clear(self):
        self.S.clear()
        self.offsety = 2



# Main Loop:
print(lsep + lsep + "h for help")
old_param = sample.parameters.copy()
old_fiterr = sample.fiterrors.copy()
last_edited = None
variables = [] # varied during fit
while 1:
    #if sample.DB_Table=="User": print("a - add f1,f2 values to database")
    check = raw_input("$ ")
    if not check:
        continue
    os.system(clcomm)
    if check=="q":
        break
    elif check in ["?", "h", "help", "info"]:
        print(lsep.join(info))
    elif check=="u":
        sample.parameters.update(old_param)
        sample.fiterrors.update(old_fiterr)
        for i_M in xrange(sample.number_of_measurements):
            newy[i_M] = sample.reflectogram(theta[i_M], i_M)
            pl.setp(fitted_plot[i_M], xdata=theta[i_M], ydata=newy[i_M])
            pl.setp(measured_plot[i_M], xdata=theta[i_M])
    elif check=="s":
        fpath = get_filename("results", ROOT,".save")
        if fpath!=None:
            with open(fpath, "w") as f:
                pickle.dump((sample.parameters, sample.coupled_vars), f)
            print("Saved Parameters to %s"%fpath)
        else:
            print("nothing saved!")
        print("Saving param file:")
        fpath = get_filename("samples", ROOT + ".mod",".param")
        if fpath!=None:
            sample.save_model(fpath)
            print("Saved Model to %s"%fpath)
        else:
            print("nothing saved!")
    elif check=="l":
        flist = sorted(os.listdir("results"))
        flist = filter(lambda s: s.endswith(".save"), flist)
        print(lsep.join(flist))
        fpath = raw_input("Path to Parameters Filename? results%s"%SEP)
        fpath2 = os.path.join("results",fpath)
        if os.path.isfile(fpath2):
            fpath = fpath2
        try:
            with open(fpath, "r") as f:
                temp=pickle.load(f)
        except:
            print("Error while opening file at: %s"%fpath)
            continue
        if isinstance(temp, tuple): 
            params, coupled_vars = temp
        elif isinstance(temp, dict):
            params = temp
            coupled_vars = {}
        else:
            print("Corrupt file selected: %s"%fpath)
            continue
        if any([key not in params for key in sample.parameters]):
            print("Saved file does not fit to sample")
            continue
        else:
            sample.parameters.update(params)
            sample.coupled_vars.update(coupled_vars)
        print("Successfully loaded %s"%fpath)
    elif check=="E":
        print((sample.err**2).sum()/len(sample.err))
    elif check=="w":
        fpath = get_filename("results", ROOT,".txt")
        if fpath!=None:
            with open(fpath, "w") as f:
                f.write(sample.print_parameter())
            if sample.number_of_measurements>0:
                try:
                    f.write(lsep + "Error: " + str((sample.err**2).sum()/len(sample.err)))
                except:
                    pass
                print("Saved Parameters (ASCII) to %s"%fpath)
            for i_M in range(sample.number_of_measurements):
                if len(sample.measured_data[i_M][:,0]) == len(theta[i_M]):
                    mpath = fpath.replace(".txt", "_M%i.dat"%i_M)
                    data = pl.vstack((theta[i_M],
                                         sample.reflectogram(theta[i_M], i_M),
                                         sample.measured_data[i_M][:,1])).T
                    header="theta R_calc R_meas"
                else:
                    data = pl.vstack((theta[i_M],
                                         sample.reflectogram(theta[i_M], i_M))).T
                    header="theta R_calc"
                pl.savetxt(mpath, data, header=header)
                print("Saved Fit %i (ASCII) to %s"%(i_M, mpath))
        else:
            print("nothing saved!")
    elif check=="d":
        fpath = get_filename("results", ROOT, ".dat")
        if fpath!=None:
            from_z = raw_input("From (nm) ? ")
            to_z = raw_input("To (nm) ? ")
            try:
                z1 = float(from_z)
                z2 = float(to_z)
            except: 
                print("Enter Floats!")
                continue
            zpoints = pl.linspace(z1, z2, 2001)*10 #nm
            pl.savetxt(fpath,
                       pl.vstack((zpoints, sample.density(zpoints))).T, 
                       header="depth(A) density(g/cm3) delta beta")
            print("Saved density profile data to %s"%fpath)
        else:
            print("nothing saved!")
    elif check=="p":
        print(sample.print_parameter())
    elif check=="P":
        print sample.parameters.keys()
    elif check=="c":
        sample.set_coupled_parameters()
    elif check=="profile":
        sample.set_profile()
    elif check=="a":
        key = raw_input("Enter key of parameter to add: ")
        if key=="":
            continue
        elif key in sample.parameters:
            print("parameter already defined!")
            continue
        descr = raw_input("Enter description of parameter to add: ")
        if descr=="":
            descr = key
        sample.add_parameter(key, 0, descr)
    elif check=="r":
        for i_M in range(sample.number_of_measurements):
            newy[i_M]=sample.reflectogram(theta[i_M], i_M)
            pl.setp(fitted_plot[i_M], xdata=theta[i_M], ydata=newy[i_M])
            pl.setp(measured_plot[i_M], xdata=theta[i_M])
    elif check in ["nb", "n", "na", "ns", "nf", "np", "ncg"]:
        old_param = sample.parameters.copy()
        old_fiterr = sample.fiterrors.copy()
        if len(newy)==sample.number_of_measurements:
            for i_M in range(sample.number_of_measurements):
                pl.setp(initial_plot[i_M], xdata=theta[i_M], ydata=newy[i_M])
        if not variables:
            print("No variables defined. "
                  "Use edit menu `e' to select variables with <space>.")
            continue
        try:
            if check=="nb": 
                branges = []
                for var in variables:
                    while True:
                        print("Ctrl-c to abort...")
                        try:
                            ll = float(raw_input("Lower limit for %s (%.5g): "%(var, sample.parameters[var])))
                            ul = float(raw_input("Upper limit for %s (%.5g): "%(var, sample.parameters[var])))
                            branges.append((ll,ul))
                            break
                        except KeyboardInterrupt:
                            ul = None
                            break
                        except:
                            continue
                    if ul==None:
                        break
                if ul==None:
                    continue
                try:
                    Ns = int(raw_input("Number of points per parameter: "))
                except:
                    print("Invalid input. Need number.")
                    continue
                
                fitted_param = sample.fit("brute", variables, branges, Ns)
            elif check=="na":
                fitted_param = sample.fit("anneal", variables)
            elif check=="ns":
                fitted_param = sample.fit("fmin", variables)
            elif check=="nf":
                fitted_param = sample.fit("fmin_bfgs", variables)
            elif check=="np":
                fitted_param = sample.fit("fmin_powell", variables)
            elif check=="ncg":
                fitted_param = sample.fit("fmin_cg", variables)
            elif check=="n":
                fitted_param = sample.fit("leastsq", variables)
            else:
                raise AssertionError
        except KeyboardInterrupt:
            print(lsep + "Optimiziation process aborted...")
            answer = raw_input("Reset parameters? ([y]/n)  ")
            if answer=="n":
                pass
            else:
                sample.parameters.update(old_param)
                sample.fiterrors.update(old_fiterr.copy())
                print("Parameters reset.")
            continue
        except Exception as errmsg:
            print("An error occured!")
            print("Message: %s%sParameters reset." %(errmsg,lsep))
            sample.parameters.update(old_param)
            sample.fiterrors.update(old_fiterr)
            continue
        if check=="n":
            fpath = os.path.join("tmp","leastsq_temp.save")
        elif check=="nb":
            fpath = os.path.join("tmp","brute_temp.save")
        elif check=="na":
            fpath = os.path.join("tmp","anneal_temp.save")
        elif check in ["ns", "np", "ncg", "nf"]: 
            fpath = os.path.join("tmp","fmin_temp.save")
        else:
            raise AssertionError
        
        with open(fpath, "w") as f:
            pickle.dump((sample.parameters, sample.coupled_vars), f)
        for i_M in xrange(sample.number_of_measurements):
            newy[i_M]=sample.reflectogram(theta[i_M], i_M)
            pl.setp(fitted_plot[i_M],   xdata=theta[i_M], ydata=newy[i_M])
            pl.setp(measured_plot[i_M], xdata=theta[i_M])
        fig.canvas.draw()
        for key in sample.var_names:
            str1 = "old value for %s:%s%g" %(key, (12-len(key))*" ", old_param[key])
            if check=="n": 
                ind = sample.var_names.index(key)
                str2 = "new value: %g +-%g" %(fitted_param[key], sample.fiterrors[key])
            elif check!="n":
                str2 = "new value: %g" %fitted_param[key]
            else: raise AssertionError
            print("%s%s%s"%(str1, (27+10-len(str1))*" ", str2))
    elif check in ("e", "edit"):
        old_param = sample.parameters.copy()
        if len(newy)==sample.number_of_measurements:
            for i_M in range(sample.number_of_measurements):
                pl.setp(initial_plot[i_M], xdata=theta[i_M], ydata=newy[i_M])
        
        parlen = max(map(len, ["%s"%(str(sample.parameters[k])) for k in sample.parameters]))
        parfmt = "%%%i.5g - %%s"%parlen
        parinfo = dict([(k,parfmt%(sample.parameters[k], sample.names[k])) for k in sample.parameters])
        parreadonly = ["grad_d_%i"%i for i,N in enumerate(sample.parameters.N) if N==1]
        parreadonly += ["d_0", "N_0"]
        parreadonly.append("N_%i"%(len(sample.parameters.N)-1))
        parreadonly.append("d_%i"%(sample.total_layers-1))
        with Screen() as screen:
            signal.signal(signal.SIGWINCH, screen.resize)
            screen.clear()
            screen.draw_title("pyxrr - configuration")
            
            attr = screen.menu(sorted(parinfo), last_edited, info=parinfo, 
                               showall=True, selected=variables, 
                               readonly=parreadonly)
        if attr=="exit":
            continue
        currval = sample.parameters[attr]
        last_edited = attr
        atype = type(currval)
        value = raw_input("New %s [%s]: "\
                          %(attr, str(currval)))
        if value=="":
            continue
        try:
            value = atype(value)
        except:
            print("Wrong datatype for %s: %s"%(attr, value))
            continue
        sample.parameters[attr] = value
        for i_M in range(sample.number_of_measurements):
            newy[i_M]=sample.reflectogram(theta[i_M], i_M)
            pl.setp(fitted_plot[i_M],   xdata=theta[i_M], ydata=newy[i_M])
            pl.setp(measured_plot[i_M], xdata=theta[i_M])
        fig.canvas.draw()
    elif check.startswith("setup"):
        with Screen() as screen:
            signal.signal(signal.SIGWINCH, screen.resize)
            screen.clear()
            screen.draw_title("pyxrr - configuration")
            
            attr = screen.menu(sample.options.keys(), info=sample.info)
            if attr in sample.options and hasattr(sample.options[attr], "__iter__"):
                subopt = list(sample.options[attr])
                currind = subopt.index(getattr(sample, attr))
                value = screen.menu(subopt, currind, sample.info)
            else:
                value = None
        if value=="exit" or attr=="exit":
            continue
        if value==None:
            currval = getattr(sample, attr)
            atype = type(currval)
            value = raw_input("New %s [%s]: "\
                              %(attr, currval))
            if value=="":
                continue
            try:
                value = atype(value)
            except:
                print("Wrong datatype for %s: %s"%(attr, value))
                continue
        setattr(sample, attr, value)
        if attr=="fittype":
            sample.process_weights
        elif attr in ["DB", "DB_Table"]:
            if os.path.isfile(value):
                sample.fetch_optical_constants()
            else:
                sample.DB = currval
                print("File not found.")
            
