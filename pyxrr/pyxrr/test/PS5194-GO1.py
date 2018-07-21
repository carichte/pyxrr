import pylab as pl
from pyxrr import structure, measurement, Reflectivity



protection = structure.Layer("Si", name="protection", thickness=100.4, roughness=4.514, density=2.33)


Multilayer = structure.Group([], periods=499, roughness=2.651)
Multilayer.append(structure.Layer("Mo", name="Moly", thickness=6.4727, roughness=2.4884, density=10.118))
Multilayer.append(structure.Layer("B4C", name="Boron Carbide", thickness=8.7981, roughness=4.461, density=2.456))

substrate = structure.Layer("Si", roughness=2., name="Silicon")

stack = structure.Stack([protection, Multilayer], substrate=substrate)
stack.update()


data = pl.loadtxt("PS5194-GO1_M0.dat")
m = measurement.Measurement(data[:,0],
                            data[:,2],
                            sigma="poisson",
                            x_axis="theta",
                            energy=8.0,
                            scale=1., # corresponts to 1/I0
                            offset=0., # theta offset (deg)
                            resolution=0.0025, # beam divergence
                            polarization=0, #0=perpendicular, 1=parallel, 0.5=unpolarized
                            background=-5.947, # log10 of background intensity
                            fitrange=(0,float("inf")),
                            rebin=False, 
                            name="")


#m.polarization.value = 0.5

xrr = Reflectivity(stack, m)
pl.semilogy(m.get_theta(), xrr.reflectivity(), m.get_theta(), data[:,1])

pl.show()
#Measurement: energy=8.0, file=measurements/PS5194-GO1_az_E2.dat, offset=0.00, resolution=0.0025, background=-6.044, scale=1.0, weighting=relative
#Measurement: energy=8.048, file=measurements/PS5194-GO1_az_E2.dat, offset=0.00, resolution=0.023, background=-6.744, scale=1.0, weighting=relative
#Measurement: energy=8.048, file=measurements/ps5194-R2_GO02697_qa.raw, offset=0.00, resolution=0.023, background=-6.744, scale=1.0, twotheta=1, weighting=relative
#Measurement: energy=8.905413, file=measurements/PS5194-GO1_az_hzg4.dat, offset=0.00, resolution=0.01, background=-4.0, scale=1.0, twotheta=1, weighting=relative

