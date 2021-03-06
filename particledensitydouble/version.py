import os.path
import sys

title = "ParticleDensityDouble"
author = "Max Larsson"
version = "1.0.0"
date = ("June", "1", "2018")
email = "max.larsson@liu.se"
homepage = "www.hu.liu.se/forskning/larsson-max/software"
if hasattr(sys, 'frozen'):
    if '_MEIPASS2' in os.environ:
        path = os.environ['_MEIPASS2']
    else:
        path = sys.argv[0]
else:
    path = __file__
app_path = os.path.dirname(path)
icon = os.path.join(app_path, "pdd.ico")