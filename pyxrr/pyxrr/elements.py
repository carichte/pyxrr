import collections


Element = collections.namedtuple("Element", (
                                             "Number",
                                             "Symbol",
                                             "density",
                                             "molar_mass"
                                            ))



class ElementDict(collections.OrderedDict):
    def __getitem__(self, key):
        if isinstance(key, int):
            for element in self.values():
                if element.Number == key:
                    return element
        return super(ElementDict, self).__getitem__(key)



Elements = ElementDict()

Elements["H"] = Element(1, "H", 9.0e-05, 1.008)
Elements["He"] = Element( 2, "He", 0.000179, 4.003)
Elements["Li"] = Element( 3, "Li", 0.534, 6.941)
Elements["Be"] = Element( 4, "Be", 1.85, 9.012)
Elements["B"]  = Element( 5, "B",  2.34, 10.811)
Elements["C"]  = Element( 6, "C",  2.2, 12.011)
Elements["N"]  = Element( 7, "N",  0.00125, 14.007)
Elements["O"]  = Element( 8, "O",  0.00143, 15.999)
Elements["F"]  = Element( 9, "F",  0.0017, 18.998)
Elements["Ne"] = Element(10, "Ne", 0.0009, 20.18)
Elements["Na"] = Element(11, "Na", 0.971, 22.99)
Elements["Mg"] = Element(12, "Mg", 1.74, 24.305)
Elements["Al"] = Element(13, "Al", 2.7, 26.982)
Elements["Si"] = Element(14, "Si", 2.33, 28.086)
Elements["P"]  = Element(15, "P",  2.2, 30.974)
Elements["S"]  = Element(16, "S",  2.05, 32.066)
Elements["Cl"] = Element(17, "Cl", 0.00321, 35.453)
Elements["Ar"] = Element(18, "Ar", 0.00178, 39.948)
Elements["K"]  = Element(19, "K",  0.862, 39.098)
Elements["Ca"] = Element(20, "Ca", 1.55, 40.078)
Elements["Sc"] = Element(21, "Sc", 2.99, 44.956)
Elements["Ti"] = Element(22, "Ti", 4.54, 47.867)
Elements["V"]  = Element(23, "V",  6.11, 50.942)
Elements["Cr"] = Element(24, "Cr", 7.19, 51.996)
Elements["Mn"] = Element(25, "Mn", 7.3, 54.938)
Elements["Fe"] = Element(26, "Fe", 7.87, 55.845)
Elements["Co"] = Element(27, "Co", 8.9, 58.933)
Elements["Ni"] = Element(28, "Ni", 8.9, 58.693)
Elements["Cu"] = Element(29, "Cu", 8.96, 63.546)
Elements["Zn"] = Element(30, "Zn", 7.13, 65.39)
Elements["Ga"] = Element(31, "Ga", 6.09, 69.723)
Elements["Ge"] = Element(32, "Ge", 5.32, 72.61)
Elements["As"] = Element(33, "As", 5.73, 74.922)
Elements["Se"] = Element(34, "Se", 4.5, 78.96)
Elements["Br"] = Element(35, "Br", 3.12, 79.904)
Elements["Kr"] = Element(36, "Kr", 0.00373, 83.8)
Elements["Rb"] = Element(37, "Rb", 1.53, 85.468)
Elements["Sr"] = Element(38, "Sr", 2.54, 87.62)
Elements["Y"]  = Element(39, "Y",  4.46, 88.906)
Elements["Zr"] = Element(40, "Zr", 6.51, 91.224)
Elements["Nb"] = Element(41, "Nb", 8.57, 92.906)
Elements["Mo"] = Element(42, "Mo", 10.2, 95.94)
Elements["Tc"] = Element(43, "Tc", 11.5, 98.0)
Elements["Ru"] = Element(44, "Ru", 12.4, 101.07)
Elements["Rh"] = Element(45, "Rh", 12.4, 102.906)
Elements["Pd"] = Element(46, "Pd", 12.0, 106.42)
Elements["Ag"] = Element(47, "Ag", 10.5, 107.868)
Elements["Cd"] = Element(48, "Cd", 8.65, 112.411)
Elements["In"] = Element(49, "In", 7.31, 114.818)
Elements["Sn"] = Element(50, "Sn", 7.3, 118.71)
Elements["Sb"] = Element(51, "Sb", 6.69, 121.76)
Elements["Te"] = Element(52, "Te", 6.24, 127.6)
Elements["I"]  = Element(53, "I",  4.93, 126.904)
Elements["Xe"] = Element(54, "Xe", 0.00589, 131.29)
Elements["Cs"] = Element(55, "Cs", 1.87, 132.905)
Elements["Ba"] = Element(56, "Ba", 3.5, 137.327)
Elements["La"] = Element(57, "La", 6.17, 138.906)
Elements["Ce"] = Element(58, "Ce", 6.77, 140.116)
Elements["Pr"] = Element(59, "Pr", 6.7, 140.908)
Elements["Nd"] = Element(60, "Nd", 6.9, 144.24)
Elements["Pm"] = Element(61, "Pm", 7.0, 145.0)
Elements["Sm"] = Element(62, "Sm", 7.5, 150.36)
Elements["Eu"] = Element(63, "Eu", 5.25, 151.964)
Elements["Gd"] = Element(64, "Gd", 7.9, 157.25)
Elements["Tb"] = Element(65, "Tb", 8.23, 158.925)
Elements["Dy"] = Element(66, "Dy", 8.54, 162.5)
Elements["Ho"] = Element(67, "Ho", 8.78, 164.93)
Elements["Er"] = Element(68, "Er", 9.05, 167.26)
Elements["Tm"] = Element(69, "Tm", 9.31, 168.934)
Elements["Yb"] = Element(70, "Yb", 6.7, 173.04)
Elements["Lu"] = Element(71, "Lu", 9.84, 174.967)
Elements["Hf"] = Element(72, "Hf", 13.3, 178.49)
Elements["Ta"] = Element(73, "Ta", 16.7, 180.948)
Elements["W"]  = Element(74, "W",  19.3, 183.84)
Elements["Re"] = Element(75, "Re", 21.0, 186.207)
Elements["Os"] = Element(76, "Os", 22.6, 190.23)
Elements["Ir"] = Element(77, "Ir", 22.4, 192.217)
Elements["Pt"] = Element(78, "Pt", 21.5, 195.078)
Elements["Au"] = Element(79, "Au", 19.3, 196.967)
Elements["Hg"] = Element(80, "Hg", 13.5, 200.59)
Elements["Tl"] = Element(81, "Tl", 11.9, 204.383)
Elements["Pb"] = Element(82, "Pb", 11.4, 207.2)
Elements["Bi"] = Element(83, "Bi", 9.75, 208.98)
Elements["Po"] = Element(84, "Po", 9.32, 209.0)
Elements["At"] = Element(85, "At", 0.0, 210.0)
Elements["Rn"] = Element(86, "Rn", 0.00973, 222.0)
Elements["Fr"] = Element(87, "Fr", 0.0, 223.0)
Elements["Ra"] = Element(88, "Ra", 5.0, 226.025)
Elements["Ac"] = Element(89, "Ac", 10.1, 227.028)
Elements["Th"] = Element(90, "Th", 11.7, 232.038)
Elements["Pa"] = Element(91, "Pa", 15.4, 231.036)
Elements["U"]  = Element(92, "U",  18.9, 238.029)