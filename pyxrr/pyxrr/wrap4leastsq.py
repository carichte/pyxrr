#!/usr/bin/python
# vim:fileencoding=utf-8
# -*- coding: utf-8 -*-

# Copyright 2013 by Carsten Richter and Robert Mietrach
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





def wrap_for_fit(func, varlist):
    """
    Wraps func so that it only expects a tuple containing the values for
    the variables specified in varlist. The rest of the variables is
    supplied from argdict.

    func -- original function
    argdict -- complete dictionary of func's arguments with default values
    varlist -- variables that that returned function is still expecting


    returns a (wrapped_function, start_values) pair

    start_values -- list of values from argdict for the keys in varlist
    wrapped_function -- function which takes tuple of values corresponding
                        to the variables in varlist
    
    unpack : bool
        when True, unpacks the dictionary of parameters when passing it
        to the cost function. used by default 
    """

    def helper(t):
        if hasattr(t, "item") and t.size==1:
            t = (t.item(),)
        func_kwargs = dict(zip(varlist, t))
        return func(**func_kwargs)
    return helper
