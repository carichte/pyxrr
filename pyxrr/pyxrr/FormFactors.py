#!/usr/bin/env python
#----------------------------------------------------------------------
# Description:
# Author: Carsten Richter <carsten.richter@desy.de>
# Created at: Di 20. Okt 15:45:24 CEST 2015
# Computer: haso227r 
# System: Linux 3.13.0-65-generic on x86_64
#
# Copyright (c) 2015 Carsten Richter  All rights reserved.
#----------------------------------------------------------------------
import os
import numpy as np
import scipy.interpolate
import scatfaccoef
import functions
import sympy


PATH_COPPENS = os.path.join(os.path.dirname(__file__), "coppens")

class FormFactors(object):
    _supported_tables = ["ITC", "Coppens"]
    def __init__(self, default_table = "ITC", mathmodule="math"):
        """
            default_table can be one of:
                - "ITC" (ITC Vol. C, ch. 6.1, Table 6.1.1.4, 
                         http://it.iucr.org/Cb/ch6o1v0001/#table6o1o1o4)
                - "Coppens" (http://harker.chem.buffalo.edu/group/groupindex.html)
        """
        assert default_table in self._supported_tables, "Invalid default table"
        self.f0func = dict()
        self.qmax = dict()
        self._currtable = dict()
        self.default_table = default_table
        self._mathmodule = np
        
        self.supported_species = dict(
            ITC = scatfaccoef.f0.keys(),
            Coppens = [s.strip(".npy") for s in os.listdir(PATH_COPPENS)]
            )
        
    def __call__(self, species, q):
        """
            Returns the form factor of an element with certain oxidation state 
            (species) and for a certain value q = 2 * sin(theta)/lambda
        """
        if species not in self.f0func:
            self.add_f0(species)
        qmax = self.qmax[species]
        if np.any(q<0)==True or np.any(q>qmax)==True:
            raise ValueError("The used table is limited to "\
                  "0 < 2 * sin(theta)/lambda < %.1f A^-1. "\
                  "The value entered (%g) exceeds this."%(qmax, np.max(q)))
        
        
        return self.f0func[species](q)
    
    def add_f0(self, species, table=None):
        if table not in self._supported_tables:
            found = False
            for table in [self.default_table] + self._supported_tables:
                if species not in self.supported_species[table]:
                    print("Warning: Species %s not found in %s, trying other " 
                          "table."%(species, table))
                else:
                    found = True
                    break
            if not found:
                raise ValueError("Species %s not found in databases"%species)
        
        r = getattr(self, "_add_f0_%s"%table)(species)
        self._currtable[species] = table
        return r
    
    def _add_f0_ITC(self, species):
        assert species in scatfaccoef.f0, \
            "%s not found in ITC table of form factors."%species
        
        coef = scatfaccoef.f0[species]
        a = sympy.Matrix(coef[0:4])
        b = sympy.Matrix(coef[4:8])
        c = coef[8]
        #print a, b, c
        #return  (a,(-b*(q/2)**2).applyfunc(sympy.exp)) # + c
        q = sympy.Symbol("q", real=True, nonnegative=True)
        f0 = sum(a.multiply_elementwise((-b*q**2/4).applyfunc(sympy.exp))) + c
        f0func = functions.makefunc(f0, self._mathmodule)
        # better agains sympy here? 
        # does it really make sense to have different q values for one reflection?
        
        self.qmax[species] = 4
        self.f0func[species] = f0func
    
    def _add_f0_Coppens(self, species):
        assert species in self.supported_species["Coppens"], \
            "%s not found in Coppens table of form factors."%species
        fname = "%s.npy"%species
        path = os.path.join(PATH_COPPENS, fname)
        q, f0 = np.load(path)
        q *= 2
        self.qmax[species] = q.max()
        self.f0func[species] = scipy.interpolate.interp1d(q, f0, kind="cubic")
    
    def clear(self):
        """
            Clear the cached form factors
        """
        self.f0func.clear()
    
    def __repr__(self):
        s = super(FormFactors, self).__repr__() + os.linesep
        s += "In cache:" + os.linesep
        for ff in self.f0func:
            s += "%6s   from %s database."%(ff, self._currtable[ff])
            s += "   2*sin(theta)/lambda < %.1f%s"%(self.qmax[ff], os.linesep)
        return s