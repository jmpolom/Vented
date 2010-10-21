#!/usr/bin/python

import matplotlib.pyplot as plt
import vented

# Driver name
name = 'Dayton Audio QT305-4, 12 inch sub-woofer'

# User specified design parameters
Fs = 26.2   # Driver resonance frequency, in Hertz
Qes = 0.39  # Driver electrical Q
Qms = 9.91  # Driver mechanical Q
Re = 3.5    # Driver DC resistance
Vas = 123.9 # Driver air compliance volume, in liters
Lv = 10.0   # Enclosure vent length, in inches
D = 4.0     # Enclosure vent internal radius, in inches
Vb = 81.84  # Enclosure volume, in liters
# Use default enclosure loss factor, Ql=7

vented.response_plot(Fs,Qes,Qms,Re,Vas,D,Vb=Vb,Lv=Lv,name=name,freq_max=500,plot=plt)
vented.response_plot(Fs,Qes,Qms,Re,Vas,D,a=1.5109,h=1.0667,Vb=Vb,Lv=Lv,name=name,freq_max=500,plot=plt)
vented.impedance_plot(Fs,Qes,Qms,Re,Vas,D,Vb=Vb,Lv=Lv,name=name,freq_max=500,plot=plt)
vented.impedance_plot(Fs,Qes,Qms,Re,Vas,D,a=1.5109,h=1.0667,Vb=Vb,Lv=Lv,name=name,freq_max=500,plot=plt)
vented.displacement_plot(Fs,Qes,Qms,Re,Vas,D,Vb=Vb,Lv=Lv,name=name,freq_max=500,plot=plt)
vented.displacement_plot(Fs,Qes,Qms,Re,Vas,D,a=1.5109,h=1.0667,Vb=Vb,Lv=Lv,name=name,freq_max=500,plot=plt)
plt.show()
