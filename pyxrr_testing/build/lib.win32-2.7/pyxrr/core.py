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

import pickle
import time
from copy import copy
from functions import *
from wrap4leastsq import *
from xrr import reflectivity
import scipy.optimize as sopt
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.ndimage.filters import gaussian_filter1d



class pyxrrError(Exception):
    def __init__(self, value, errmsg="", identifier=None):
        self.value = value
        self.errmsg = errmsg
    def __str__(self):
        return self.value + lsep + "Message: %s" %self.errmsg

class multilayer(object):
    """
        multilayer class.
    """
    def __init__(self, SampleFile, DB = DB_PATH, DB_Table = "Henke", penalty=1, fittype="log", verbose=2):
        """
            multilayer object as defined in the given 'SampleFile'.
            
            Optional inputs:
            
             - DB (string): sqlite database file where dispercion correction coefficients (f1, f2), 
                            atomic weights and densities are stored.
             
             - DB_Table (string): set of data to use for dispercion correction coefficients. Can be:
                => 'BrennanCowan' (see http://skuld.bmsc.washington.edu/scatter/)
                => 'Chantler' (see http://www.nist.gov/pml/data/ffast/index.cfm)
                => 'CromerLiberman'
                => 'EPDL97'
                => 'Henke' (see http://www-cxro.lbl.gov/optical_constants/)
                => 'Sasaki'
                => 'Windt' (also used by the Program IMD, see f.e. http://www.rxollc.com/idl/)
                For more details see Databases on http://www.esrf.eu/computing/scientific/dabax
                                 and http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/
             
             - penalty (float): value to increase (if >1) or decrease (if <1) the importance of simulation 
                                values larger than measured ones in comparison to those that are smaller.
             
             - fittype (string): way to calculate the residuals:
                => 'log':    err = (log10(y_meas) - log10(y_sim(x_meas))) * weights(x_meas)
                   'root':   err = ( sqrt(y_meas) -  sqrt(y_sim(x_meas))) * weights(x_meas)
                   'linear': err = (      y_meas  -       y_sim(x_meas) ) * weights(x_meas)
                others can be implemented.
             
             - verbose (int): sets verbosity during least squares iterations
                => 0: display nothing
                => 1: display sum of squares (sum(err**2))
                => 2: display sum of squares (sum(err**2)) AND afterwards time used by calculation
        """
        
        # initial parameters if none given:
        #self.parameters={"energy0":8., "scale0":1., "offset0": 0., "background0":-20., "resolution0": 0.}
        
        # initialize process variables
        #self.err=np.inf # starting value for brute force fit
        self.timeT=0. # starting values for timer
        self.timeC=0.
        self.fcalls = 0 # starting value for amount of function calls
        self.iterations = 0 # starting value for fit iterations
        self.verbose=verbose # verbose mode
        self.fittype=fittype # take residuals from logs of functions
        self.penalty=penalty # penalty for simulation values that are greater than measured ones
        
        if verbose: print "opening sample file " + SampleFile
        self.SampleFile = SampleFile
        self.coupled_vars=dict() # initialize dictionary of coupled parameters
        
        try:
            self.parameters, self.materials, self.dims, self.names, self.measured_data, self.weights,\
            self.pol, self.fit_limits, self.number_of_measurements, self.total_layers, self.x_axes, self.paths = parse_parameter_file(SampleFile)
        except Exception as errmsg:
            raise pyxrrError("An error occured while trying to parse the parameter file.", errmsg=errmsg)
        
        self.process_fit_range()
        
        
        self.weightmethods = self.weights.copy()
        self.process_weights()
        
        # One can add individual residual functions here
        if self.fittype=="log":
            self.ResidualFunction = lambda y_m, y_s, w: (np.log10(y_m) - np.log10(y_s))*w
        elif self.fittype=="root":
            self.ResidualFunction = lambda y_m, y_s, w: (np.sqrt(y_m) - np.sqrt(y_s))*w
        else:
            self.ResidualFunction = lambda y_m, y_s, w: (y_m - y_s)*w # linear
        
        
        if verbose:
            print("fetching optical constants from %s..." %DB)
            print(lsep)
        minE, maxE = 0, np.inf
        for i in range(self.total_layers):
            tc = get_components(self.materials[i])
            for element in tc[0]: 
                E, f1, f2 = get_f1f2_from_db(element, database = DB, table = DB_Table).T # Preload to save time later
                minE = max(minE, E.min())
                maxE = min(maxE, E.max())
        ind = (E>=minE) * (E<=maxE)
        newE = E[ind]
        delta, beta = [], [] # initialize some lists
        for i in range(self.total_layers):
            thisdelta, thisbeta = get_optical_constants(1., self.materials[i], newE, database = DB, table = DB_Table)
            delta.append(thisdelta)
            beta.append(thisbeta)
            
        self.optical_constants = interp1d(newE, (delta,beta), kind='linear', bounds_error=False) # Interpolated function for all layers
        
        
        self.fiterrors = self.parameters.copy()
        self.fiterrors.pop("LayerCount")
        self.fiterrors.pop("N")
        self.fiterrors.pop("d_0")
        self.fiterrors.pop("d_" + str(self.total_layers-1))
        for i in range(len(self.parameters["N"])):
            if self.parameters["N"][i] == 1:
                self.fiterrors.pop("grad_d_" + str(i))
        for key in self.fiterrors.keys():
            self.fiterrors[key] = np.nan
        
    
    
    def process_fit_range(self):
        if self.verbose: print("Processing ranges for fit...")
        self.fit_range = {}
        for key in self.fit_limits.iterkeys():
            try:
                a, b = self.fit_limits[key]
                a, b = float(a), float(b)
                a, b = min(a,b), max(a,b)
            except Exception as errmsg:
                raise pyxrrError("fit limits not understood (2-tuple of floats expected).", errmsg=errmsg)
            if self.verbose:
                print("  Measurement no %s: from %.2g to %.2g"%(key, a, b))
            data = self.measured_data[key]
            self.fit_range[key] = (data[:,0] > a)*\
                                  (data[:,0] < b)*\
                                  (data[:,1] > 0)
        if self.verbose: 
            print(lsep)
    
    def process_weights(self):
        if self.verbose: 
            print("Processing weights for fit...")
        for key in range(self.number_of_measurements):
            if self.weightmethods.has_key(key):
                if self.weightmethods[key] == "statistical":
                    if self.fittype == "log":
                        self.weights[key] = np.sqrt(self.measured_data[key][:,1])
                    elif self.fittype == "linear":
                        self.weights[key] = 1./np.sqrt(self.measured_data[key][:,1])
                    else:
                        self.weights[key] = 1.
                elif self.weightmethods[key] == "z":
                    try:
                        self.weights[key] = self.measured_data[key][:,2]
                    except IndexError:
                        raise pyxrrError("Error: no 3rd column found for weighting in measurement %i.%sDisable weighting or add data." %(key,lsep))
                else:
                    if self.verbose:
                        print("  weighting not understood for measurement %i."%key)
                    self.weights[key] = 1.
            else:
                if self.verbose:
                    print("  No weighting for measurement %i."%key)
                self.weights[key] = 1.
        if self.verbose: 
            print(lsep)
    
    def print_parameter(self):
        keylist=self.parameters.keys()
        keylist.sort()
        param_table="Sample: " + self.SampleFile + lsep
        for key in keylist:
            try: 
                if key in self.fiterrors.keys() and not np.isnan(self.fiterrors[key]):
                    errstr = "+-%.3g"%self.fiterrors[key]
                    errstr = (10-len(errstr))*" " + errstr
                    param_table += key + (11-len(key))*" " + "%12.5g"%self.parameters[key] + errstr + 2*" " + str(self.names[key]) + lsep
                else:
                    param_table += key + (11-len(key))*" " + "%12.5g"%self.parameters[key] + 12*" " + str(self.names[key]) + lsep
            except: param_table += key + (11-len(key))*" " + str(self.parameters[key]) + 12*" " + str(self.names[key]) + lsep
        return param_table
    
    def stack(self):
        """
            Returns the ideal stack with its parameters. Just try :-)
            (Does not take aperiodicities into account yet.)
        """
        N=self.parameters["N"]
        lc=self.parameters["LayerCount"]
        layerID=0
        sigmaID=-1
        z=0
        for N_i in range(len(N)):
            for i in range(N[N_i]):
                for j in range(lc[N_i]):
                    lastz=z
                    z+=self.parameters["d_" + str(layerID+j)]
                    if N_i==0 and i==0 and j==0: 
                        result=np.array([layerID, 
                            -np.inf, 
                            z,  
                            self.parameters["d_" + str(layerID)], 
                            self.parameters["rho_" + str(layerID)], 
                            0
                        ])
                    else: 
                        result=np.vstack((result, 
                             np.array([layerID+j, 
                                       lastz, z, 
                                       self.parameters["d_" + str(layerID+j)], 
                                       self.parameters["rho_" + str(layerID+j)], 
                                       self.parameters["sigma_" + str(sigmaID+j)]
                                     ])
                            ))
            print layerID, sigmaID
            layerID+=lc[N_i]
            sigmaID+=lc[N_i] + (lc[N_i]>1)*(N[N_i]>1)
        result[-1,2]=np.inf
        return result
    
    def concentration(self, depth):
        """
            calculates fraction that each layer contributes at given depth (float or array).
            (Does not take aperiodicities into account yet.)
        """
        from scipy.special import erf
        
        w2=np.sqrt(2)
        stack=self.stack()
        try: 
            z=float(depth)
            ind=(stack[:,1]<z) * (stack[:,2]>=z)
            result=np.zeros(self.total_layers)
            i=ind.argmax()
            if ind[0]==True:
                result[stack[i,0]]+=.5*(1-erf((z-stack[i,2])/(w2*stack[i+1,5])))
                result[stack[i+1,0]]+=.5*(1+erf((z-stack[i,2])/(w2*stack[i+1,5])))
            elif ind[-1]==True:
                result[stack[i-1,0]]+=.5*(1-erf((z-stack[i,1])/(w2*stack[i,5])))
                result[stack[i,0]]+=.5*(1+erf((z-stack[i,1])/(w2*stack[i,5])))
            else:
                result[stack[i-1,0]]+=.5*(1-erf((z-stack[i,1])/(w2*stack[i,5])))
                result[stack[i,0]]+=.5*(1+erf((z-stack[i,1])/(w2*stack[i,5])))-.5*(1+erf((z-stack[i,2])/(w2*stack[i+1,5])))
                result[stack[i+1,0]]+=.5*(1+erf((z-stack[i,2])/(w2*stack[i+1,5])))
            return result
        except:
            result=np.zeros((len(depth), self.total_layers+1))
            for j in range(len(depth)):
                ind=(stack[:,1]<depth[j]) * (stack[:,2]>=depth[j])
                result[j,0]=depth[j]
                i=ind.argmax()
                if ind[0]==True:
                    result[j,stack[i,0]+1]+=.5*(1-erf((depth[j]-stack[i,2])/(w2*stack[i+1,5])))
                    result[j,stack[i+1,0]+1]+=.5*(1+erf((depth[j]-stack[i,2])/(w2*stack[i+1,5])))
                elif ind[-1]==True:
                    result[j,stack[i-1,0]+1]+=.5*(1-erf((depth[j]-stack[i,1])/(w2*stack[i,5])))
                    result[j,stack[i,0]+1]+=.5*(1+erf((depth[j]-stack[i,1])/(w2*stack[i,5])))
                else:
                    result[j,stack[i-1,0]+1]+=.5*(1-erf((depth[j]-stack[i,1])/(w2*stack[i,5])))
                    result[j,stack[i,0]+1]+=.5*(1+erf((depth[j]-stack[i,1])/(w2*stack[i,5])))-.5*(1+erf((depth[j]-stack[i,2])/(w2*stack[i+1,5])))
                    result[j,stack[i+1,0]+1]+=.5*(1+erf((depth[j]-stack[i,2])/(w2*stack[i+1,5])))
            return result
    
    def density(self, depth):
        """
            calculates the mass density and optical constants at given depth (float or array).
            (Does not take aperiodicities into account yet.)
        """
        rho=np.array([self.parameters["rho_" + str(i)] for i in range(self.dims["rho"])])
        #delta, beta = get_optical_constants(rho, self.materials, self.parameters["energy0"]*1000.)
        delta,beta = self.optical_constants(self.parameters["energy0"]*1000)*rho
        try:
            z=float(depth)
            rho_out=np.array([self.concentration(depth)[i]*rho[i] for i in range(self.total_layers)]).sum()
            delta_out=np.array([self.concentration(depth)[i]*delta[i] for i in range(self.total_layers)]).sum()
            beta_out=np.array([self.concentration(depth)[i]*beta[i] for i in range(self.total_layers)]).sum()
        except:
            rho_out=(np.array([self.concentration(depth)[:,i+1]*rho[i] for i in range(self.total_layers)]).sum(0))
            delta_out=(np.array([self.concentration(depth)[:,i+1]*delta[i] for i in range(self.total_layers)]).sum(0))
            beta_out=(np.array([self.concentration(depth)[:,i+1]*beta[i] for i in range(self.total_layers)]).sum(0))
        return rho_out, delta_out, beta_out
    
    def set_coupled_parameters(self, equation=None, remove=None):
        """
            One can add equations to the dictionary of coupled parameters by
            giving an equation string here in the form: 'p1 = function(p2)+p3'.
            
            To delete coupled parameters, one has to give the name of the 
            dependent parameter as the 'remove' argument.
            
            If nothing is supplied, the function demands raw_input to work.
        """
        if equation==None and remove==None:
            print(self.print_parameter())
            while 1:
                command=raw_input("Enter equation: ")
                if command=="": break
                else:
                    command = command.replace(" ", "")
                    var, eq = command.split("=")
                    self.coupled_vars[var] = eq
            while 1:
                for var in self.coupled_vars.keys(): print("%s = %s" %(var,self.coupled_vars[var]))
                remove=raw_input("Decouple parameter? ")
                if remove=="": break
                elif self.coupled_vars.has_key(remove):
                    print("Removed equation: %s = %s" %(remove,self.coupled_vars.pop(remove)))
                else: continue
        elif remove==None: 
            try: var, eq = equation.split("=")
            except: raise pyxrrError("Given equation is not valid: %s" %str(equation))
            self.coupled_vars[var] = eq
        elif equation==None and self.coupled_vars.has_key(remove):
            if verbose: print("Removed equation: %s = %s" %(remove,self.coupled_vars.pop(remove)))
            else: self.coupled_vars.pop(remove)
    
    
    
    def couple(self):
        """
            update coupled parameters
        """
        loc_param =  self.parameters.copy()
        for var in self.coupled_vars.keys(): self.parameters[var] = eval(self.coupled_vars[var], loc_param)
    
    def reflectogram(self, theta_range=None, i_M=0):
        """
            Calculates a reflectivity curve for a given theta range.
            The theta array has to be equally spaced data and sorted.
            i_M specifies the measurement from which to take parameters like energy, offset, polarization etc.
        """
        if theta_range==None: theta_range = self.measured_data[i_M][:,0]
        self.couple() # call coupled parameters
        N=np.array(self.parameters["N"])
        LayerCount=np.array(self.parameters["LayerCount"])
        energy=self.parameters["energy" + str(i_M)]
        d=np.array([abs(self.parameters["d_" + str(i)]) for i in range(self.dims["d"])])
        grad_d=np.array([self.parameters["grad_d_" + str(i)] for i in range(self.dims["grad_d"])])
        sigma=np.array([abs(self.parameters["sigma_" + str(i)]) for i in range(self.dims["sigma"])])
        rho = np.array([(abs(self.parameters["rho_" + str(i)])) for i in range(self.dims["rho"])])
        delta,beta = self.optical_constants(energy*1000)*rho
        if len(theta_range)>1: blur_sigma=abs(self.parameters["resolution" + str(i_M)]/2.35482)/(theta_range[1]-theta_range[0])
        else: blur_sigma=0
        if self.verbose==2: timeC0=time.time()
        R = reflectivity(theta_range-self.parameters["offset" + str(i_M)], N, grad_d, LayerCount, d, delta, beta, sigma, 12.398/energy, self.pol[i_M])
        if self.verbose==2:
            self.timeC+=(time.time()-timeC0)
            self.fcalls+=1
        return gaussian_filter1d(abs(self.parameters["scale%i"%i_M])*R + 10**self.parameters["background%i"%i_M], blur_sigma)
    
    def residuals(self, new_parameters={}, fitalg=None):
        if fitalg is not None: self.fitalg=fitalg
        """
            Calculates the residuals (Array of deviations between measurements and simulations).
            Inputs:
                new_parameters - dictionary with names and values of
                                 parameters which shall be updated
                fitalg - string indicating the fit algorithm
                         Can be:
                             'brute', 'fmin', 'anneal', 'fmin_bfgs' etc.:
                                sum of squares of residuals is returned
                             
                             'leastsq': the array of deviations is returned
        """
        if self.verbose==2: timeT0=time.time()
        self.parameters.update(new_parameters)
        self.err=np.array([])
        for i_M in range(self.number_of_measurements):
            x_m = self.measured_data[i_M][self.fit_range[i_M],0]
            y_m = self.measured_data[i_M][self.fit_range[i_M],1]
            if hasattr(self.weights[i_M], "__iter__"):
                w = self.weights[i_M][self.fit_range[i_M]]
            else: w = self.weights[i_M]
            # get simulated reflectivity curve
            y_s = self.reflectogram(x_m, i_M)
            self.err = np.append(self.err, self.ResidualFunction(y_m, y_s, w))
        
        # discard simulation flaws
        ind = np.isnan(self.err)
        if ind.all():
            self.err[ind] = np.inf
        elif ind.any() and self.verbose==2:
            print("Warning: NaN values in residuals")
        else:
            # NaN is 10 times as bad as maximum deviation?
            self.err[ind] = 10*abs(self.err[~ind]).max() 
        if len(self.err)==0:
            raise pyxrrError("No measured data in fit range")
        if self.penalty != 1.:
            # better larger or smaller simulation values?
            self.err[self.err<0]*=self.penalty
        self.iterations+=1
        if self.fitalg =="leastsq":
            if self.verbose: 
                print("Iterations: %i%sValue: %.12f"%(self.iterations, 5*" ", (self.err**2).sum()/len(self.err)))
                if self.verbose==2: self.timeT+=(time.time()-timeT0)
            return self.err
        else:
            if self.verbose==2:
                if self.fitalg=="brute":
                    print("Iterations: %i/%i  (%.2f%%)" %(self.iterations, self.TotalCalls, (100.*self.iterations)/self.TotalCalls))
                else:
                    print("Iterations: %i%sValue: %.12f" %(self.iterations, 5*" ", (self.err**2).sum()/len(self.err)))
                self.timeT+=(time.time()-timeT0)
            return (self.err**2).sum()/len(self.err)
    
    def get_fit_parameters(self, limits = False, complete=False):
        """
            Asks for parameters that shall be varied during fit.
            This function should not belong to the class.
            It will be moved to the outer part in later versions.
        """
        ranges = ()
        varlist=list()
        keylist=self.fiterrors.keys()
        for p in self.coupled_vars:
            keylist.pop(keylist.index(p))
        keylist.sort()
        if complete: return keylist
        for i in keylist:
            answer=raw_input("Vary " + self.names[i] + "   " + str(self.parameters[i]) +  "    (y/[n])? ")
            if answer=="y":
                if limits:
                    while True:
                        try: 
                            minimum=float(raw_input("minimum? "))
                            maximum=float(raw_input("maximum? "))
                            break
                        except ValueError as errmsg:
                            print("No valid float entered: %s" %errmsg)
                            continue
                        except Exception as errmsg:
                            raise pyxrrError("Unexpected error.", errmsg = errmsg)
                    ranges += (sorted((minimum, maximum)),)
                varlist.append(i)
        if limits:
            while True:
                Ns = raw_input("Number of samples? ")
                try:
                    Ns = int(Ns)
                    break
                except KeyboardInterrupt:
                    break
                except:
                    pass
            return varlist, ranges, Ns
        else: return varlist
    
    def fit(self, algorithm="leastsq", var_names=[], ranges=[], Ns=20):
        """
            Method for fitting simulation to measured data by minimization of
            the residuals function.
            Inputs:
                
                algorithm : Optimization algorithm to be used. (see scipy.optimize)
                    - 'leastsq': Levenberg-Marquardt algorithm
                    - 'brute': brute force iteration through a parameter range
                    - 'anneal': ...
                    - 'fmin': simplex
                    - 'fmin_bfgs'
                    - 'fmin_powell'
                    - 'fmin_cg'
                
                var_names - list of names of parameters to vary (=variables)
                            see self.get_fit_parameters(complete=True)
                            for a complete list.
            
            Only if algorithm=='brute':
                ranges - a list of 2-tuples giving the ranges in which the 
                         variables will be sampled
                Ns - number of samples per variable
        """
        start_dict=self.parameters.copy()
        if var_names==[]:
            print "Ctrl+C to abort"
            try:
                if algorithm=="brute": 
                    var_names, ranges, Ns = self.get_fit_parameters(limits=True)
                elif algorithm in ["leastsq", "anneal", "fmin", "fmin_bfgs", "fmin_powell", "fmin_cg"]:
                    var_names = self.get_fit_parameters(limits=False)
            except KeyboardInterrupt:
                print("Aborted...")
                self.var_names=[]
                return start_dict, None
            except Exception as errmsg:
                self.var_names=[]
                raise pyxrrError("Failed to fetch variable parameters.", errmsg = errmsg)
            if var_names==[]:
                self.var_names=var_names
                return start_dict, None
        elif algorithm=="brute" and len(var_names)!=len(ranges):
            raise pyxrrError("For brute force fit length of variables list and ranges list has to be equal.")
        self.var_names=var_names
        func, start_val = wrap_for_fit_dict(self.residuals, start_dict, var_names)
        self.timeT=0. # starting values for timer
        self.timeC=0.
        self.fcalls = 0
        self.iterations = 0
        self.fitalg = algorithm
        
        # Do the fit...
        if algorithm=="brute":
            self.TotalCalls = Ns**len(ranges)
            param = sopt.brute(func, ranges, Ns=Ns, full_output=False, finish=None)
        elif algorithm=="leastsq":
            output = sopt.leastsq(func, start_val, full_output=True, ftol=2**-20, xtol=2**-20)
            param = output[0]
        elif algorithm=="anneal":
            output = sopt.anneal(func, start_val, schedule='boltzmann', full_output=True, learn_rate=0.5, dwell=50)
            param = output[0]
        elif algorithm=="fmin":
            output = sopt.fmin(func, start_val, full_output=True, maxfun=1000*len(start_val), maxiter=1000*len(start_val))
            param = output[0]
        elif algorithm=="fmin_bfgs":
            output = sopt.fmin_bfgs(func, start_val, full_output=True)
            param = output[0]
        elif algorithm=="fmin_powell":
            output = sopt.fmin_powell(func, start_val, full_output=True)
            param = output[0]
        elif algorithm=="fmin_cg":
            output = sopt.fmin_cg(func, start_val, full_output=True)
            param = output[0]
        if len(var_names)>1: fitresult = dict(zip(var_names, param))
        else: fitresult = dict({var_names[0]:param.item()})
        
        # some parameters are nonegative...
        for key in fitresult.keys():
            if key.startswith("d_") or key.startswith("rho_") or key.startswith("sigma_") or key.startswith("scale"):
                fitresult[key] = abs(fitresult[key])
        self.parameters.update(fitresult)
        if self.verbose==2:
            print("Python Time per function call: %f ms" %(1000.*(self.timeT-self.timeC)/self.fcalls))
            print("C Time per function call: %f ms" %(1000.*self.timeC/self.fcalls))
            print(str(self.fcalls) + " Function calls")
        
        for key in self.var_names:
            ind = self.var_names.index(key)
            if algorithm != "leastsq" or output[1]==None: self.fiterrors[key] = np.nan
            else: self.fiterrors[key] = np.sqrt(output[1][ind,ind])
        
        if algorithm=="leastsq" and output[1]==None: print("Covariance Matrix could not be estimated. (see scipy.optimize.leastsq)")
        return fitresult
    
    def save_model(self, savepath=None):
        """
            Saves the model into a .param file at a given path.
            
            Parameters
            ----------
            savepath : filename path or file handle
        """
        i_l = 0 # index for layers
        i_s = 0 # index for sigma
        names = copy(self.names["Names"])
        if isinstance(savepath, str):
            f = open(savepath, "w")
        elif hasattr(savepath, "readline"):
            f = savepath
        else:
            import StringIO
            f = StringIO.StringIO()
        try:
            f.write("Ambience: name=%s, code=%s, rho=%g%s" \
                  %(names.pop(0), self.materials[i_l], self.parameters["rho_%i"%i_l], lsep))
            for i in range(1, len(self.parameters["LayerCount"]) - 1): # iterate Groups but not Ambience nor Substrate
                f.write(lsep + "Group: name=%s, sigma=%g, periods=%i, grad_d=%g%s" \
                      %(names.pop(0), self.parameters["sigma_%i"%i_s], self.parameters["N"][i], self.parameters["grad_d_%i"%i], lsep))
                i_s += 1
                multilayer_sigma = 0
                for j in range(self.parameters["LayerCount"][i]):
                    i_l += 1
                    if j == 0:
                        if self.parameters["N"][i] == 1:
                            f.write("Layer: name=%s, code=%s, rho=%g, d=%g%s"\
                                  %(names.pop(0), self.materials[i_l], self.parameters["rho_%i"%i_l], self.parameters["d_%i"%i_l], lsep))
                        else:
                            special_sigma = i_s + self.parameters["LayerCount"][i] - 1
                            multilayer_sigma = 1
                            f.write("Layer: name=%s, code=%s, rho=%g, d=%g, sigma=%g%s"\
                                  %(names.pop(0), self.materials[i_l], self.parameters["rho_%i"%i_l], self.parameters["d_%i"%i_l], self.parameters["sigma_%i"%special_sigma], lsep))
                    else:    
                        f.write("Layer: name=%s, code=%s, rho=%g, d=%g, sigma=%g%s"\
                              %(names.pop(0), self.materials[i_l], self.parameters["rho_%i"%i_l], self.parameters["d_%i"%i_l], self.parameters["sigma_%i"%i_s], lsep))
                        i_s += 1
                i_s += multilayer_sigma
            
            i_l+=1
            f.write(lsep + "Substrate: name=%s, code=%s, rho=%g, sigma=%g%s"\
                  %(names.pop(0), self.materials[-1], self.parameters["rho_%i"%i_l], self.parameters["sigma_%i"%i_s], lsep))
            f.write(lsep + lsep)
            
            assert(len(names)==0)
            
            for i_M in range(self.number_of_measurements):
                newline =  "Measurement: "
                newline += "file=%s"%self.paths[i_M]
                newline += ", x_axis=%s"%self.x_axes[i_M]
                if self.weightmethods.has_key(i_M):
                    newline += ", weighting=%s"%self.weightmethods[i_M]
                newline += ", fit_range=%g->%g"%(self.fit_limits[i_M])
                newline += ", energy=%g"%self.parameters["energy%i"%i_M]
                newline += ", resolution=%g"%self.parameters["resolution%i"%i_M]
                newline += ", offset=%g"%self.parameters["offset%i"%i_M]
                newline += ", scale=%g"%self.parameters["scale%i"%i_M]
                newline += ", background=%g"%self.parameters["background%i"%i_M]
                newline += ", pol=%g"%self.pol[i_M]
                newline += lsep
                f.write(newline)
        
            if hasattr(f, "getvalue"):
                output = f.getvalue()
            else:
                output = None
        except Exception as errmsg:
            raise pyxrrError("An error occured while trying to save the model to %s."%f, errmsg=errmsg)
        finally:
            f.close()
        
        if self.verbose and output==None:
            print("Model saved in %s" %savepath)
        return output
