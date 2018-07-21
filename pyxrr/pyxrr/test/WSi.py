import pylab as pl
from pyxrr import structure, measurement, Reflectivity


protection = structure.Layer("Si", name="protection", thickness=74.546, roughness=11.405, density=2.3513)


Multilayer = structure.Group(periods=250, roughness=2.3731)
#Multilayer.append(structure.Layer("W", name="Tungsten", thickness=6.9328, roughness=2.0, density=15.873))
#Multilayer.append(structure.Layer("Si", name="Silicon", thickness=13.306, roughness=2.1965, density=2.6676))
Multilayer.append(structure.Layer("W", name="Tungsten", thickness=6.9328, roughness=2.1965, density=15.873))
Multilayer.append(structure.Layer("Si", name="Silicon", thickness=13.306, roughness=2.3601, density=2.6676))

substrate = structure.Layer("Si", roughness=2., density=2.34, name="Wafer")

stack = structure.Stack([protection, Multilayer], substrate=substrate)
stack.update()


data = pl.loadtxt("WSi_M0.dat")
m = measurement.Measurement(data[:,0],
                            data[:,2],
                            energy=8.0,
                            scale=1., # corresponts to 1/I0
                            offset=-0.0039093, # theta offset (deg)
                            resolution=0.001, # beam divergence
                            polarization=0, #0=perpendicular, 1=parallel, 0.5=unpolarized
                            background=-6.6287, # log10 of background intensity
                            rebin=True, 
                            name="")


#m.polarization.value = 0.5

xrr = Reflectivity(stack, m)
pl.semilogy(m.get_theta(), xrr.reflectivity(), m._x_orig, data[:,1])

pl.show()
