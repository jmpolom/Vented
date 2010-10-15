#!/usr/bin/python

import matplotlib.pyplot as plt
import vented_box

# User specified design parameters
Fs = 32.0  # Driver resonance frequency, in Hertz
Qts = 0.26 # Total driver Q, dimensionless
Vas = 2.76 # Driver air compliance volume, in cubic feet
Lv = 4.0   # Enclosure vent length, in inches
R = 0.75   # Enclosure vent internal radius, in inches
Vb = 0.75  # Enclosure volume, in cubic feet
Ql = 7     # Enclosure loss factor

vented_box.response_plot(Fs,Qts,Vas,Lv,R,Vb,Ql,plot=plt)
plt.show()
