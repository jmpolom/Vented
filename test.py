#!/usr/bin/python

import matplotlib.pyplot as plt
import vented

# User specified design parameters
Fs = 26.2   # Driver resonance frequency, in Hertz
Qes = 0.39  # Driver electrical Q
Qms = 9.91  # Driver mechanical Q
Re = 3.5    # Driver DC resistance
Vas = 123.9 # Driver air compliance volume, in cubic feet
Lv = 10.0   # Enclosure vent length, in inches
D = 4.0     # Enclosure vent internal radius, in inches
Vb = 81.84  # Enclosure volume, in cubic feet
# Use default enclosure loss factor, Ql=7

vented.response_plot(Fs,Qes,Qms,Re,Vas,Vb,Lv,D,freq_max=500,plot=plt)
plt.show()
