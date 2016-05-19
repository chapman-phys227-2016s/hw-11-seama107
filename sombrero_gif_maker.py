#! /usr/bin/env python

"""
File: sombrero_gif_maker.py
Copyright (c) 2016 Michael Seaman
License: MIT

Description: Draws upon the sombrero class to create an animated gif
via commandline
1st argument is the number of steps

"""
from sombrero import Sombrero
import sys
import subprocess
import numpy as np

num_frames = 100
step_size = 2 * np.pi
precision = 1000
h = step_size / precision
a = 0
b = step_size * 1000
n = int((b-a)/h)


s = Sombrero(.4, 0, 0)

print "Calculating points"
s.calc_n_rk4_steps(n, h)
print "Function approximated."

for theta in np.linspace(0, .99, num_frames):
    s.generate_poincare_graph(precision, theta = theta, xLim = (-1.5, 1.5), yLim = (-1.1, 1.1))
print "Frames Generated."
    
subprocess.call("convert -delay 5 -loop 0 *.png Sombrero.gif", shell=True)
subprocess.call("rm *.png", shell=True)
print "Process complete!"
