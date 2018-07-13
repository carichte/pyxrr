import collections
import numpy as np
import lmfit
import itertools
from scipy.interpolate import interp1d
from structure import Parameter





class Measurement(object):
    _ids = itertools.count(0)
    _maxbins = 1e5
    def __init__(self, xvalues,
                       reflectivity,
                       sigma="poisson",
                       x_axis="theta",
                       energy=8.048,
                       scale=1., # corresponts to 1/I0
                       offset=0., # theta offset (deg)
                       resolution=0., # beam divergence
                       polarization=0, #0=perpendicular, 1=parallel, 0.5=unpolarized
                       background=-10, # log10 of background intensity
                       fitrange=(0,np.inf),
                       rebin=False, # interpolate or average
                       name=""
                ):
        self.id = _id = next(self._ids)
        self.name = ("Measurement #%i"%_id)
        if name:
            self.name += ": %s"%name

        self.x = np.array(xvalues, ndmin=1, dtype=float)
        self.reflectivity = np.array(reflectivity, ndmin=1, dtype=float)
        assert len(self.x)==len(self.reflectivity), \
            "lengths mismatch of xvalues and reflectivity arrays"

        if isinstance(sigma, np.ndarray):
            assert len(self.x)==len(sigma), "wrong length of uncertainties"
        elif sigma=="poisson":
            sigma = np.sqrt(reflectivity)
        elif sigma is None:
            sigma = np.ones_like(reflectivity)

        self._x_orig = self.x.copy()
        self._reflectivity_orig = self.reflectivity.copy()
        self._sigma_orig = self.sigma.copy()

        self.x_axis = x_axis
        # fittable parameters:
        self.energy = Parameter(name="energy_%i"%_id, value=energy, min=0)
        self.scale = Parameter(name="scale_%i"%_id, value=scale, min=0)
        self.offset = Parameter(name="offset_%i"%_id, value=offset)
        self.resolution = Parameter(name="resolution_%i"%_id, value=resolution, min=0)
        self.polarization = Parameter(name="polarization_%i"%_id, value=polarization, min=0)
        self.background = Parameter(name="background_%i"%_id, value=background)

        self.fitrange = fitrange
        self.rebin = rebin
        self.rebin_data()
        self.valid = (self.sigma > 0) * (self.x > fitrange[0]) * (self.x < fitrange[-1])


    def rebin_data(self):
        x = self._x_orig
        y = self._reflectivity_orig
        sigma = self._sigma_orig

        isort = x.argsort()
        x = x[isort]
        y = y[isort]
        sigma = sigma[isort]

        diff = np.diff(x)
        meddiff = np.median(diff)
        if np.allclose(diff, meddiff, rtol=0.01):
            self.x = x
            self.reflectivity = y
            self.sigma = sigma
            self.regular = True
            return

        elif self.rebin:
            if isinstance(self.rebin, bool):
                step = diff.max()
            else:
                step = self.rebin

            newx = np.arange(x[0], x[-1]+step/10., step)
            newborders = np.arange(x[0]-step/2., x[-1]+step/2.+step/10., step)

            num = np.histogram(x, bins=newborders, weights=None)[0]
            newy = np.histogram(x, bins=newborders, weights=y)[0]
            newsigma = np.histogram(x, bins=newborders, weights=sigma**2)[0]
            newsigma = np.sqrt(newsigma)
            newy /= num
            self.x = newx
            self.reflectivity = newy
            self.sigma = newsigma
            self.regular = True

        else:
            self.regular = False
            step = max(diff.min(), (x[-1]-x[0])/self._maxbins)
            self.x_regular = np.arange(x[0], x[-1]+step/10., step)


    def get_params(self):
        return (self.energy, 
                self.scale,
                self.offset,
                self.resolution,
                self.polarization,
                self.background)


