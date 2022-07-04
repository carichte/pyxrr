import collections
import numpy as np
import lmfit
import itertools
from .structure import Parameter

try:
    import pandas as pd
    PANDAS = True
except:
    PANDAS = False




class Measurement(object):
    _ids = itertools.count(0)
    _maxbins = 1e5
    const = 0.9866014922266592 # 12.398/(4*np.pi)
    xfunc = dict({
                  'theta':lambda x,E: x,
                  'twotheta':lambda x,E: x/2.,
                  'qz_nm':lambda x,E: np.degrees(np.arcsin(x/E*const/10.)),
                  'qz_a' :lambda x,E: np.degrees(np.arcsin(x/E*const))
            })

    sort_order =  ["scale", "offset", "resolution", "polarization", "background", "sample_length"]

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
                       sample_length=np.inf, # mm
                       beam_size=0.01, # mm
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

        self.sigma = sigma
        self._x_orig = self.x.copy()
        self._reflectivity_orig = self.reflectivity.copy()
        self._sigma_orig = self.sigma.copy()

        self.x_axis = x_axis
        # fittable parameters:
        self.energy = Parameter(parent=self,
                                name="energy_%i"%_id,
                                value=energy,
                                min=0)
        self.scale = Parameter(parent=self,
                               name="scale_%i"%_id,
                               value=scale,
                               min=0)
        self.offset = Parameter(parent=self,
                                name="offset_%i"%_id,
                                value=offset)
        self.resolution = Parameter(parent=self,
                                    name="resolution_%i"%_id,
                                    value=resolution,
                                    min=0)
        self.polarization = Parameter(parent=self,
                                      name="polarization_%i"%_id,
                                      value=polarization,
                                      min=0, max=1)
        self.background = Parameter(parent=self,
                                    name="background_%i"%_id,
                                    value=background)
        self.sample_length = Parameter(parent=self,
                                    name="sample_length_%i"%_id,
                                    value=sample_length)

        self.beam_size = beam_size
        self.fitrange = fitrange
        self.rebin = rebin
        self.rebin_data()
        self.valid = (self.sigma > 0) \
                   * (self.x >= fitrange[0]) \
                   * (self.x <= fitrange[-1])


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
            self.is_regular = True
            self.x_regular = self.x
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
            self.is_regular = True
            self.x_regular = self.x

        else:
            self.is_regular = False
            step = max(diff.min(), (x[-1]-x[0])/self._maxbins)
            self.x_regular = np.arange(x[0], x[-1]+step/10., step)

    def get_theta(self, x=None, regular=False): # base quantity
        if x is None:
            x = self.x if not regular else self.x_regular
        return self.xfunc[self.x_axis](x, self.energy.value)

    def get_params(self):
        return (self.energy, 
                self.scale,
                self.offset,
                self.resolution,
                self.polarization,
                self.background,
                self.sample_length)

    def _collect_data(self):
        data = collections.OrderedDict()
        data["name"] = self.name
        data["range"] = self.x.min(), self.x.max()
        data["x axis"] = self.x_axis
        for key in self.sort_order:
            data[key] = getattr(self, key).value
        data["rebining"] = self.rebin
        data["fit range"] = self.fitrange
        return data

    def to_DataFrame(self):
        data = self._collect_data()
        for k in data:
            data[k] = [data[k]]
        return pd.DataFrame(data, index=[self.id])

    def _repr_html_(self):
        if PANDAS:
            return self.to_DataFrame().to_html()
        else:
            return super(Stack, self).__repr__()

    def __repr__(self):
        if PANDAS:
            return self.to_DataFrame().to_string()
        else:
            return super(Stack, self).__repr__()


if __name__ == "__main__":
    theta = np.linspace(0, 2, 101)
    R = np.exp(-theta)
    m = Measurement(theta, R)

