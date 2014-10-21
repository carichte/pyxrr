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

import numpy as np


def load_njc(filename):
    with open(filename, 'rb') as fd:
        z = fd.read(2000)
    try:
        offset1 = int(z[1324])
    except:
        offset1 = int(z[1328])
    offset2 = ord(z[1348+offset1])
    offset = 1436 + offset1 + offset2
    
    fd = open(filename, 'rb')
    fd.read(offset)
    data32 = np.fromfile(file=fd, dtype=np.float32)
    fd.close()
    fd = open(filename, 'rb')
    fd.read(offset)
    data64 = np.fromfile(file=fd, dtype=np.float64)
    fd.close()
    
    start = round(data64[20], 6)
    th_range = round(data64[21], 6) - start
    step = round(data32[44], 6)
    number = int(round(th_range / step)) + 1
    time = round(data32[46], 6)
    
    angle = start + step * np.arange(number)
    intens = data32[-40-number:-40] * time
    return np.vstack((angle, intens)).T



def load_x00(filename):
    with open(filename, 'r') as f:
        header_length = 1
        for line in f:
            if "FirstAngle" in line:
                first = float(line.split()[1])
            if "StepWidth" in line:
                step = float(line.split()[1])
            if "TimePerStep" in line:
                time = float(line.split()[1])
            if "NrOfData" in line:
                number = int(line.split()[1])
            if "ScanData" in line:
                break
            else:
                header_length += 1
    angle = first + step * np.arange(number)
    intens = np.loadtxt(filename, skiprows=header_length) * time
    return np.vstack((angle, intens)).T



def load_fio(FILENAME):
    assert type(FILENAME)==str, AssertionError("Filename has to be of type string.")
    data=np.array([])
    with open(FILENAME, "r") as fiodata:
        while True:
            line = fiodata.readline()
            if not line: break
            if not line.find("!"): continue
            if "%c" in line:
                line = fiodata.readline()
                while not (line.find("!") == 0):
                    if not line: break
                    line = fiodata.readline()
            if "%p" in line:
                line = fiodata.readline()
                while not (line.find("!") == 0):
                    if not line: break
                    line = fiodata.readline()
            if "%d" in line:
                line = fiodata.readline()
                while not (line.find("!") == 0):
                    if not line: break
                    if "Col" in line: pass
                    elif not len(data): data=np.array(line.split())
                    else: data=np.vstack(( data, np.array(line.split()) ))
                    line = fiodata.readline()
    return data.astype(float)

def load_asc(FILENAME):
    with open(FILENAME, "r") as fn:
        zeilen=fn.readlines()
    y_array=np.array([])
    for zeile in zeilen:
        if not "data" in zeile:
            if "start angle" in zeile:
                start_angle=float(zeile.split()[0])
            if "step width" in zeile:
                step_width=float(zeile.split()[0])
            if "number of steps" in zeile:
                step_num=float(zeile.split()[0])
        else:
            y_array=np.append(y_array, float(zeile.split()[0]))
    x_array=np.arange(start_angle, step_num*step_width+start_angle, step_width)
    return np.vstack((x_array, y_array)).T


def load_val(FILENAME):
    with open(FILENAME, "r") as fn:
        seen_intens = False
        Intensity = np.array([])
        while True:
            zeile = fn.readline()
            if not zeile: break
            if "BerPar01" in zeile:
                zeile = fn.readline()
                om_start=float(zeile)
                zeile = fn.readline()
                twth_start=float(zeile)
                zeile = fn.readline()
                twth_stop=float(zeile)
                zeile = fn.readline()
                delta_twth=float(zeile)
                zeile = fn.readline()
                zeile = fn.readline()
                zeile = fn.readline()
                zeile = fn.readline()
                steps = float(zeile)
                continue
            if "Intens" in zeile:
                seen_intens = True
                continue
            if seen_intens:
                Intensity = np.append(Intensity, float(zeile))
    fn.close()
    twotheta=np.linspace(twth_start, twth_stop, steps)
    return np.vstack((twotheta, Intensity)).T

def load_raw(FILENAME):
    with open(FILENAME, 'rb') as fd:
        data_32 = np.fromfile(file=fd, dtype=np.float32)
    with open(FILENAME, 'rb') as fd:
        data_64 = np.fromfile(file=fd, dtype=np.float64)
    with open(FILENAME, 'rb') as fd:
        data_i = np.fromfile(file=fd, dtype=np.int32)
    with open(FILENAME, 'rb') as fd:
        data_s = np.fromfile(file=fd, dtype='a326')
    # Evaluate data
    start = data_64[91]
    step = data_64[111]
    length = data_i[179]
    index = 254 + data_i[242] / 4.
    twotheta = start + step * np.arange(length)
    intens = data_32[index:index+length]
    return np.vstack((twotheta, intens)).T


