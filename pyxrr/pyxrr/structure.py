import collections
import numpy as np
import lmfit
import itertools


class Parameter(lmfit.Parameter):
    _container = np.empty(1)
    @property
    def value(self):
        "The numerical value of the Parameter, with bounds applied"
        return self._getval()

    @value.setter
    def value(self, val):
        "Set the numerical Parameter value."
        self._val = val
        #print val
        self._container[:] = val
        if not hasattr(self, '_expr_eval'):  self._expr_eval = None
        if self._expr_eval is not None:
            self._expr_eval.symtable[self.name] = val



class Layer(object):
    _ids = itertools.count(0)
    periods=1
    def __init__(self, composition, density, roughness=0., thickness=np.inf, name=""):
        if not name:
            name = composition
        self.id = _id = next(self._ids)
        self.name = "Layer #%i: %s"%(_id, name)

        self.thickness = Parameter(name="Thickness_%i"%_id, 
                                   value=thickness, min=0)
        self.density = Parameter(name="Density_%i"%_id, 
                                 value=density, min=0)
        self.roughness = Parameter(name="Roughness_%i"%_id, 
                                   value=roughness, min=0)
        self.composition = composition
    
    def get_params(self):
        return self.thickness, self.density, self.roughness
    
    def __add__(self, nextLayer):
        if isinstance(nextLayer, Layer):
            return Group((self, nextLayer), 1, self.roughness)
        elif isinstance(nextLayer, Group):
            nextLayer.insert(0, self)
            return nextLayer
    


Air = Layer("N0.78O.21", 0.00125, name="Air")


class Group(list):
    _ids = itertools.count(0)
    def __init__(self, layers=(), periods=1, roughness=0., name=""):
        if any(not isinstance(layer, Layer) for layer in layers):
            raise ValueError("All group elements must be an instance "
                             "of %s"%str(Layer))
        super(Group, self).__init__(layers)
        self.id = _id = next(self._ids)
        self.roughness = Parameter(name="GroupRoughness_%i"%_id, 
                                   value=roughness, min=0)
        self.periods = periods
        self.name = name if name else "Group #%i"%self.id

    @property
    def periods(self):
        return self._periods
    @periods.setter
    def periods(self, val):
        # the roughness value is only significant if more than 1 period:
        self.roughness.vary = val>1
        self._periods = val

    def get_params(self):
        return self.roughness,

    def append(self, layer):
        if not isinstance(layer, Layer):
            raise ValueError("All group elements must be an instance "
                             "of %s"%str(Layer))
        super(Group, self).append(layer)
    def insert(self, index, layer):
        if not isinstance(layer, Layer):
            raise ValueError("All group elements must be an instance "
                             "of %s"%str(Layer))
        super(Group, self).insert(index, layer)



class Stack(list):
    def __init__(self, groups, substrate, ambience=Air, name=None):
        super(Stack, self).__init__(groups)
        self.insert(0, ambience)
        self.append(substrate)
        self.params = lmfit.Parameters()

    def iter_layers(self, interfaces=False):
        """
            Returns a generator that iterates over all layers
            in the stack
        """
        for group in self:
            if isinstance(group, Group):
                if interfaces and group.periods>1:
                    yield group
                for layer in group:
                    yield layer
            elif isinstance(group, Layer):
                yield group

    def numlayers(self, interfaces=False):
        i=0
        for group in self:
            if isinstance(group, Group):
                if interfaces and group.periods>1:
                    i+=1
                for layer in group:
                    i+=1
            elif isinstance(group, Layer):
                i+=1
        return i

    def update(self):
        self.params = lmfit.Parameters()
        nl = self.numlayers()
        ni = self.numlayers(interfaces=True)

        self._thicknesses = np.empty(nl)
        self._densities   = np.empty(nl)
        self._roughnesses = np.empty(ni)

        # some wild stuff
        for iL, layer in enumerate(self.iter_layers()):
            self._thicknesses[iL:iL+1] = layer.thickness.value
            self._densities[iL:iL+1] = layer.density.value

            layer.thickness._container = self._thicknesses[iL:iL+1]
            layer.density._container = self._densities[iL:iL+1]

        for iI, interf in enumerate(self.iter_layers(True)):
            self._roughnesses[iI:iI+1] = interf.roughness.value
            interf.roughness._container = self._roughnesses[iI:iI+1]

            print interf.name
            self.params.add_many(*interf.get_params())

        return self.params







Absorber = Layer("Mo", 5., 1., 5., "Moly")
Spacer = Layer("B4C",  2., 1., 10., "Boron Carbide")

Multilayer = Group((Absorber, Spacer), 100, 2.)

Sample = Stack([Multilayer], Layer("Si", 5, 1))



