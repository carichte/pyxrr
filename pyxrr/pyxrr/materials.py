import os
import appdirs
import shelve
import collections
from elements import elements
import re

datadir = appdirs.user_data_dir("pyxrr")
if not os.path.isdir(datadir):
    os.makedirs(datadir)

datafile = os.path.join(datadir, "materials.shelve")

Material = collections.namedtuple("Material", (
                                             "name",
                                             "composition",
                                             "density",
                                            ))



def get_components(compount, reduce_only=False):
    """
        Parse a chemical formulas of arbitrary complexity.
        
        Input:
            compount : str
                chemical formula
                e.g.: (SiO2)0.73(Na2O)0.05(K2O)0.17(CaO)0.03(Al2O3)0.02
            
            reduce_only : bool
                if True, returns sum formula only
        
        Returns:
            an elements, amount tuple where
                elements - list of element symbols
                amount   - list of quantity for each element
    """
    myfloat = lambda s: float(s) if s else 1.
    while True:
        #print compount
        parenthesis = compount.count("(")
        assert parenthesis==compount.count(")"), "Parenthesis error."
        if not parenthesis:
            break
        inner = re.findall(r'(\([^()]*\))(\d*\.*\d*)', compount)
        for (group, amount) in inner:
            fac = myfloat(amount)
            components = re.findall(r'([A-Z][a-z]?)(\d*\.*\d*)', group)
            newgroup = ["%s%g"%(k, myfloat(v)*fac) for (k,v) in components]
            compount = compount.replace(group+amount, "".join(newgroup))
    if reduce_only:
        return compount
    components = re.findall(r'([A-Z][a-z]?)(\d*\.*\d*)', compount)
    result = collections.defaultdict(float)
    for (k,v) in components:
        result[k] += myfloat(v)
    return result.keys(), result.values()


def check_compount(compount):
    components, number = get_components(compount)
    for component in components:
        if component in elements:
            compount = compount.replace(component, "")
    residue = filter(str.isalpha, compount)
    print residue
    return not bool(residue)



def add(composition, density, name=None):
    if name is None:
        name = composition
    d = shelve.open(datafile)
    material = (name, composition, density)
    d[name] = material
    d[composition] = material
    return Material(*material)


def get(key):
    if key in elements:
        element = elements[key]
        material = Material(element.symbol, element.symbol, element.density)
        return material
    d = shelve.open(datafile)
    material = Material(*d[key])
    return material

def keys():
    d = shelve.open(datafile)
    return d.keys()


if not os.path.isfile(datafile):
    add("AgBr", 6.473) # "from henke"
    add("AlAs", 3.81) # "from henke"
    add("AlP", 2.42) # "from henke"
    add("B4C", 2.52) # "from henke"
    add("BeO", 3.01) # "from henke"
    add("BN", 2.25) # "from henke"
    add("CdWO4", 7.9) # "from henke"
    add("CdS", 4.826) # "from henke"
    add("CoSi2", 5.3) # "from henke"
    add("Cr2O3", 5.21) # "from henke"
    add("CsI", 4.51) # "from henke"
    add("CuI", 5.63) # "from henke"
    add("InN", 6.88) # "from henke"
    add("In2O3", 7.179) # "from henke"
    add("InSb", 5.775) # "from henke"
    add("GaAs", 5.316) # "from henke"
    add("GaN", 6.1) # "from henke"
    add("GaP", 4.13) # "from henke"
    add("LiF", 2.635) # "from henke"
    add("LiH", 0.783) # "from henke"
    add("LiOH", 1.43) # "from henke"
    add("MgO", 3.58) # "from henke"
    add("Mg2Si", 1.94) # "from henke"
    add("MnO", 5.44) # "from henke"
    add("MnO2", 5.03) # "from henke"
    add("MoO2", 6.47) # "from henke"
    add("MoO3", 4.69) # "from henke"
    add("MoSi2", 6.31) # "from henke"
    add("NiO", 6.67) # "from henke"
    add("Ni2Si", 7.2) # "from henke"
    add("Ru2Si3", 6.96) # "from henke"
    add("SiC", 3.217) # "from henke"
    add("Si3N4", 3.44) # "from henke"
    add("TaN", 16.3) # "from henke"
    add("TiN", 5.22) # "from henke"
    add("Ta2Si", 14.0) # "from henke"
    add("UO2", 10.96) # "from henke"
    add("VN", 6.13) # "from henke"
    add("ZnO", 5.675) # "from henke"
    add("ZnS", 4.079) # "from henke"
    add("ZrN", 7.09) # "from henke"
    add("Mo", 10.28)
    add("Si", 2.33)
    add("LiNbO3", 4.64)
    add("SrTiO3", 5.13)
    add("SrO", 4.08)
    
    add("SiO2", 2.2, "Silica")
    add("SiO2", 2.65, "Qartz")
    add("TiO2", 4.26, "Rutile")
    add("Si.925Ti.075O2", 2.205, "ULE")
    add("H2O", 1.0, "Water")
    add("Si.56Al.5P.16Li.04Ti.02Zr.02Zn.03O2.46", 2.53, "Zerodur")
    add("ZrO2", 5.6, "Zirconia")
    add("Al2O3", 3.97, "Sapphire")
    add("C22H10N2O5", 1.43 , "Polyimide")
    add("C3H6", 0.9, "Polypropylene")
    add("C5H8O2", 1.19, "PMMA")
    add("C16H14O3", 1.2, "Polycarbonate")
    add("C16H14O3", 1.2, "Kimfol")
    add("C10H8O4", 1.4, "Mylar")
    add("C2F4", 2.2, "C2F4")
    add("C8H7Cl", 1.29, "Parylene-C")
    add("C8H8", 1.11, "Parylene-N")
    add("CaF2", 3.18, "Fluorite")
    add("KAl3Si3O12H2", 2.83, "Mica")
    add("NaCl", 2.165, "Salt")



