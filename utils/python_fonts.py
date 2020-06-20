#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
font_paths = mpl.font_manager.findSystemFonts()
font_objects = mpl.font_manager.createFontList(font_paths)
font_names = [f.name for f in font_objects]

font_names = sorted(list(set(font_names)))

print('\n'.join(font_names))
