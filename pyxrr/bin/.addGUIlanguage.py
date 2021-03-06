#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Erstellt die POT-Datei (Uebersetzungsvorlage).

Requirements
    - Python: http://www.python.org/
"""

import os
import sys
import subprocess

APPDIR = os.path.abspath(".")
LOCALEDIR = os.path.join(os.path.split(APPDIR)[0], "pyxrr", "locale")
LOCALEDOMAIN = "pyxrr_GUI"


#PYGETTEXT = os.path.join(sys.exec_prefix, "Tools/i18n/pygettext.py")
#EXECUTABLE = sys.executable
EXECUTABLE= "pygettext"

# Argumente
args = [
    EXECUTABLE,
    #PYGETTEXT,
    "--default-domain=%s" % LOCALEDOMAIN,
    "--add-location",
    "--output-dir=%s" % LOCALEDIR,
    "--output=%s.pot" % LOCALEDOMAIN,
    "--verbose",
    # Dateien
    os.path.join(APPDIR, "pyxrr_GUI_v0.7.pyw")
]

# Ausführen
print EXECUTABLE, args, APPDIR
subprocess.call(args = args, executable = EXECUTABLE, cwd = APPDIR)

print "Fertig"
