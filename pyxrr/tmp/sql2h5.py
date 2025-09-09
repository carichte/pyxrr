import pyxrr
import h5py

elements = pyxrr.elements.elements.keys()

tables = ['BrennanCowan', 'Chantler', 'CromerLiberman', 'EPDL97', 'Henke', 'Sasaki', 'Windt']

with h5py.File("f1f2c.h5") as fh:
    for table in tables:
        g = fh.require_group(table)
        for element in elements:
            try:
                E, f1, f2 = pyxrr.xi.get_f1f2_from_db(element, table=table)
            except:
                continue
            e = g.create_group(element)
            e.create_dataset("energy", data=E.astype('f4'), compression="gzip")
            e.create_dataset("f", data=(f1+1j*f2).astype('c8'), compression="gzip")
