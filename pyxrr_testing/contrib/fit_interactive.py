#!/usr/bin/env python

# Copyright 2012 by Carsten Richter
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
import sys
import warnings
import pickle
import pylab


sys.path.insert(0, os.path.abspath("..")) # Path to pyxrr if it has not been installed as package
import pyxrr
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

FITTYPE = "log"
VERBOSE = 2
PENALTY = 1.
DATABASE_f1f2 = "Henke" # can be BrennanCowan, Chantler, CromerLiberman, EPDL97, Henke, Sasaki, Windt
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

print(lsep)
print("Choose sample file:")
for i in sorted(os.listdir("samples")):
    if ".param" in i: print " %s"%i

SAMPLENAME=raw_input(lsep + "Enter file name: samples" + SEP)
ROOT = SAMPLENAME.rpartition(".")[0]
if not ROOT: ROOT = SAMPLENAME


while True:
    try:
        sample=pyxrr.multilayer(os.path.join("samples",SAMPLENAME), DB_Table = DATABASE_f1f2, verbose=VERBOSE, penalty=PENALTY, fittype=FITTYPE)
        break
    except pyxrr.pyxrrError as perr:
        if ("Density of material" in str(perr.errmsg)) and ("not found in database" in str(perr.errmsg)):
            print perr
            material = str(perr.errmsg).split("'")[1]
            rho = float(raw_input("Please enter density for %s: "%material))
            pyxrr.store_rho_to_db(material, rho, database = DATABASE_file)
        else: 
            raise(perr)


pylab.ion()
fig=pylab.figure(1)
axis=[]
theta=[]
measured_plot=[]
initial_plot=[]
fitted_plot=[]
refl=[]
newy={}

for i_M in range(sample.number_of_measurements):
    axis.append(fig.add_subplot(sample.number_of_measurements, 1, i_M+1))
    if len(sample.measured_data[i_M][:,0])>0:
        theta.append(sample.measured_data[i_M][:,0])
        measured_plot.append(axis[i_M].semilogy(theta[i_M], sample.measured_data[i_M][:,1]))
    else:
        print("Only plotting supported for measurement %i"%i_M)
        while 1:
            try:
                thmin = float(raw_input("Start from theta = "))
                thmax = float(raw_input("Stop at theta = "))
                thnum = int(raw_input("Number of points = "))
                theta.append(pyxrr.np.linspace(thmin, thmax, thnum))
                break
            except KeyboardInterrupt:
                theta.append(sample.measured_data[i_M][:,0])
                break
            except Exception as errmsg:
                print("%sWrong input!%sMessage: %s" %(lsep,lsep,errmsg))
                continue
    refl = sample.reflectogram(theta[i_M], i_M)
    initial_plot.append(axis[i_M].semilogy(theta[i_M], refl))
    fitted_plot.append( axis[i_M].semilogy(theta[i_M], refl))
    axis[i_M].set_ylabel("Reflectivity")
    pylab.grid(True)

if len(axis)>0:
    axis[-1].set_xlabel("theta / degree")

print(lsep + lsep + "h for help")


# Main Loop:

while 1:
    #if sample.DB_Table=="User": print("a - add f1,f2 values to database")
    check = ""
    while not check: check=raw_input("$ ")
    os.system(clcomm)
    if check=="q": break
    elif check in ["?", "h", "help", "info"]:
        print """n - leastsq refinement (see scipy.optimize documentation)
    n   : leastsq (Levenberg-Marquardt algorithm)
    nf  : fmin_bfgs (BFGS algorithm)
    ns  : fmin (downhill simplex algorithm)
    nb  : brute (brute force)
    np  : fmin_powell (Powell's algorithm)
    ncg : fmin_cg (nonlinear conjugate gradient algorithm)
e - edit parameter
q - quit
c - set coupled parameters
u - undo last iteration
s - save parameters
l - load parameters
w - write ascii
d - write density vs. depth profile
r - refresh plot
p - print parameters
h/? - this message"""
    
    elif check=="u":
        sample.parameters=old_param.copy()
        sample.fiterrors=old_fiterr.copy()
        for i_M in range(sample.number_of_measurements):
            newy[i_M]=sample.reflectogram(theta[i_M], i_M)
            pylab.setp(fitted_plot[i_M], ydata=newy[i_M])
    elif check=="s":
        FILENAME=raw_input("Filename (without extension)? [" + ROOT + "] ")
        if not FILENAME: FILENAME = ROOT
        fpath = os.path.join("results",FILENAME + ".save")
        if os.path.isfile(fpath) and not (raw_input("File %s exists. Overwrite? (y/[n])"%fpath).lower() == "y"):
            pass
        else:
            with open(fpath, "w") as f:
                pickle.dump((sample.parameters, sample.coupled_vars), f)
            print "Saved Parameters to " + fpath
        fpath = os.path.join("samples",FILENAME + ".mod.param")
        if os.path.isfile(fpath) and not (raw_input("File %s exists. Overwrite? (y/[n])"%fpath) == "y").lower():
            continue
        else:
            sample.save_model(fpath)
            #print "Saved Model to " + fpath
    elif check=="l":
        for filename in sorted(os.listdir("results")):
            if ".save" in filename: print(filename)
        FILENAME=raw_input("Saved Parameters Filename ? results" + SEP)
        fpath = os.path.join("results",FILENAME)
        try:
            with open(fpath, "r") as f:
                temp=pickle.load(f)
        except:
            print "Bad File: %s" %fpath
            continue
        if isinstance(temp, tuple): params, coupled_vars = temp
        elif isinstance(temp, dict):
            params = temp
            coupled_vars = {}
        else:
            print "Corrupt file selected: %s"%fpath
            continue
        if [i in sample.parameters.keys() for i in params.keys()].count(False) > 0:
            print("Saved file does not fit to Sample")
            continue
        else:
            sample.parameters.update(params)
            sample.coupled_vars.update(coupled_vars)
        print("Loaded %s"%fpath)
    elif check=="E":
        print((sample.err**2).sum()/len(sample.err))
    elif check=="w":
        FILENAME=raw_input("Filename (without extension)? [" + ROOT + "] ")
        if not FILENAME: FILENAME = ROOT
        fpath = os.path.join("results",FILENAME + ".txt")
        if os.path.isfile(fpath) and not raw_input("File exists. Overwrite? (y/[n])").lower() == "y":
            continue
        with open(fpath, "w") as f:
            f.write(sample.print_parameter())
            if sample.number_of_measurements>0:
                try: f.write(lsep + "Error: " + str((sample.err**2).sum()/len(sample.err)))
                except: f.write(lsep + "Error: " + str((sample.residuals({}, fitalg="leastsq")**2).sum()/len(sample.err)))
        for i_M in range(sample.number_of_measurements):
            if len(sample.measured_data[i_M][:,0])==len(theta[i_M]):
                pylab.savetxt(os.path.join("results", FILENAME + "_M%i"%i_M + ".dat"), pylab.vstack((theta[i_M], sample.reflectogram(theta[i_M], i_M), sample.measured_data[i_M][:,1])).T)
            else:
                pylab.savetxt(os.path.join("results", FILENAME + "_M%i"%i_M + ".dat"), pylab.vstack((theta[i_M], sample.reflectogram(theta[i_M], i_M))).T)
        print "Saved Parameters (ASCII) to %s"%fpath
    elif check=="d":
        FILENAME=raw_input("Filename ? (without extension)? [" + ROOT + "] ")
        if not FILENAME: FILENAME = ROOT
        fpath = os.path.join("results", FILENAME + ".dat")
        if os.path.isfile(fpath) and not raw_input("File exists. Overwrite? (y/[n])").lower() == "y":
            continue
        from_z=raw_input("From (nm) ? ")
        to_z=raw_input("To (nm) ? ")
        try:
            z1=float(from_z)
            z2=float(to_z)
        except: 
            print "Enter Floats!"
            continue
        zpoints=pylab.linspace(z1, z2, 2001)*10 #nm
        pylab.savetxt(fpath, pylab.vstack((zpoints, sample.density(zpoints))).T)
        pyxrr.write_header(fpath, "depth(A) density(g/cm3) delta beta"+lsep)
        print "Saved density profile data to " + fpath
    elif check=="p": print(sample.print_parameter())
    elif check=="P": print sample.parameters.keys()
    elif check=="c": sample.set_coupled_parameters()
    elif check=="r":
        for i_M in range(sample.number_of_measurements):
            newy[i_M]=sample.reflectogram(theta[i_M], i_M)
            pylab.setp(fitted_plot[i_M], ydata=newy[i_M])
    elif check in ["nb", "n", "na", "ns", "nf", "np", "ncg"]:
        old_param = sample.parameters.copy()
        old_fiterr = sample.fiterrors.copy()
        if len(newy)==sample.number_of_measurements:
            for i_M in range(sample.number_of_measurements): pylab.setp(initial_plot[i_M], ydata=newy[i_M])
        try:
            if check=="nb": fitted_param = sample.fit("brute")
            elif check=="na": fitted_param = sample.fit("anneal")
            elif check=="ns": fitted_param = sample.fit("fmin")
            elif check=="nf": fitted_param = sample.fit("fmin_bfgs")
            elif check=="np": fitted_param = sample.fit("fmin_powell")
            elif check=="ncg": fitted_param = sample.fit("fmin_cg")
            elif check=="n": fitted_param = sample.fit("leastsq")
            else: raise AssertionError
        except KeyboardInterrupt:
            print(lsep + "Optimiziation process aborted...")
            answer = raw_input("Reset parameters? ([y]/n)  ")
            if answer=="n": pass
            else:
                sample.parameters=old_param.copy()
                sample.fiterrors = old_fiterr.copy()
                print("Parameters reset.")
            continue
        except Exception as errmsg:
            print("%sAn error occured!%sMessage: %s%sParameters reset." %(lsep,lsep,errmsg,lsep))
            sample.parameters=old_param.copy()
            sample.fiterrors = old_fiterr.copy()
            continue
        if check=="n": fpath = os.path.join("tmp","leastsq_temp.save")
        elif check=="nb": fpath = os.path.join("tmp","brute_temp.save")
        elif check=="na": fpath = os.path.join("tmp","anneal_temp.save")
        elif check in ["ns", "np", "ncg", "nf"]: fpath = os.path.join("tmp","fmin_temp.save")
        else: raise AssertionError
        with open(fpath, "w") as f:
            pickle.dump((sample.parameters, sample.coupled_vars), f)
        for i_M in range(sample.number_of_measurements):
            newy[i_M]=sample.reflectogram(theta[i_M], i_M)
            pylab.setp(fitted_plot[i_M], ydata=newy[i_M])
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
    elif check=="e":
        old_param=sample.parameters.copy()
        if len(newy)==sample.number_of_measurements:
            for i_M in range(sample.number_of_measurements): pylab.setp(initial_plot[i_M], ydata=newy[i_M])
        print(sample.print_parameter())
        edit_key=raw_input("Parameter Key ? ")
        if sample.parameters.has_key(edit_key):
            try:
                if edit_key=="N":
                    edit_value = []
                    ind=1
                    for num in sample.parameters["N"]: 
                        try: edit_value.append(int(raw_input("New number of periods for group %i (old:%i) :"%(ind,num))))
                        except: edit_value.append(num)
                        ind+=1
                else:
                    edit_value=float(raw_input("Value ? "))
                sample.parameters.update([(edit_key, edit_value)])
                sample.fiterrors.update([(edit_key, pyxrr.np.nan)])
                for i_M in range(sample.number_of_measurements):
                    newy[i_M]=sample.reflectogram(theta[i_M], i_M)
                    pylab.setp(fitted_plot[i_M], ydata=newy[i_M])
                fig.canvas.draw()
            except: print "Please enter numerical data"

