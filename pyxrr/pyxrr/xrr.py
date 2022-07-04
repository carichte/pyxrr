import os
import glob
import ctypes as ct

import numpy as np
from numpy import uint32, int16, double, complex_, array


# _libxrr = np.ctypeslib.load_library('libxrr', os.path.dirname(__file__))
libfile = glob.glob(os.path.join(os.path.dirname(__file__), 'libxrr.*'))
_libxrr = ct.CDLL(libfile[0])


_libxrr.reflectivity.argtypes = [
                                 np.ctypeslib.ndpointer(dtype = complex_), # output
                                 np.ctypeslib.ndpointer(dtype = double),  # theta
                                 np.ctypeslib.ndpointer(dtype = double),  # thickness
                                 np.ctypeslib.ndpointer(dtype = double),  # roughness
                                 np.ctypeslib.ndpointer(dtype = complex_), # n
                                 ct.c_double,                             # lambda
                                 ct.c_double,                             # polarization
                                 np.ctypeslib.ndpointer(dtype = uint32),  # periods
                                 np.ctypeslib.ndpointer(dtype = uint32),  # groupsize
                                 np.ctypeslib.ndpointer(dtype = int16),  # periodic
                                 ct.c_int,  # ntheta
                                 ct.c_int,  # ngroups
                                 ct.c_int,  # nlayers
                                 ct.c_int,  # ninterfaces
                                 ct.c_int,  # numthreads
                            ]

def reflectivity(theta,
                 thickness,
                 roughness,
                 n,
                 wavelength,
                 polarization,
                 periods,
                 groupsize,
                 numthreads=0,
                 output=None,
                 check=True):
    """
        Wrapper for the libxrr C routine.
        Calculation of x-ray reflectivity.

        Inputs:
            theta : ndarray (1d)
                incidence angle of the x-rays

            roughness : ndarray (1d,double), RMS angstrom
                contains roughnesses of each unique interface (top to bottom)

            thickness : ndarray (1d,double), angstrom
                contains thicknesses of each unique layer

            n : ndarray (1d,complex)
                contains refactive index of each unique layer

            wavelength : double, angstrom
                x-ray wavelength

            polarization : double
                degree of polarization 
                (0 - perpendicular, 1 - parallel to the scattering plane)

            periods : sequence of int
                number of periods for each group

            groupsize : sequence of int
                number of unique layers per group

            numthreads : int
                number of thread to use for parallel computing (default all)

            output : ndarray (1d, theta.size)
                array to store result (perpendicular and/or  parallel)

            check : bool
                whether to convert the arguments (risk of segfault if not)

    """

    if check:
        theta = array(theta, dtype=double, ndmin=1)
        thickness = array(thickness, dtype=double, ndmin=1)
        roughness = array(roughness, dtype=double, ndmin=1)
        n = array(n, dtype=complex_, ndmin=1)
        periods = array(periods, dtype=uint32, ndmin=1)
        groupsize = array(groupsize, dtype=uint32, ndmin=1)
    assert len(thickness)==len(n)
    ntheta = uint32(len(theta))
    ngroups = uint32(len(periods))
    nlayers = uint32(len(thickness))
    ninterfaces = uint32(len(roughness))
    periodic = np.ones(ngroups, dtype=int16)

    if output is None:
        output = np.empty(ntheta, dtype=complex_)

    _libxrr.reflectivity(output,
                         theta,
                         thickness,
                         roughness,
                         n,
                         wavelength,
                         polarization,
                         periods,
                         groupsize,
                         periodic,
                         ntheta,
                         ngroups,
                         nlayers,
                         ninterfaces,
                         numthreads)

    if polarization>=0: # reflectivity
        return output.real

    return output # amplitude



#if __name__=="__main__":
#    th = np.linspace(0,1,10001)
#    R = reflectivity(th, [1,1,2], [.2,.2,.4,2], [1+1j, 1-.1j,1+.2j], 1.5, 0, [1,1,1], [1,1,1])


