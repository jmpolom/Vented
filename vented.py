#!/usr/bin/python

"""
vented.py provides functions that calculate the approximate (model predicted), theoretical frequency response of a vented box loudspeaker system. The script uses mathematical model analysis results from A. N. Theile and Richard Small, as published in 'Loudspeakers in Vented Boxes' and 'Vented Box Loudspeaker Systems', respectively.

Script Name:    vented_box.py
Author:         Jonathan Polom
Date Created:   8 October 2010
License:        GNU General Public License
"""

import datetime
import numpy as np

__author__ = 'Jonathan Polom <s0nic0nslaught@gmail.com>'
__date__ = datetime.date(2010, 10, 19)
__version__ = 0.8

def params_calc(Fs,Qes,Qms,Vas,D,Vb,Lv):
    """
    Calculates dependent system parameters: enclosure tuning frequency, system tuning ratio, system compliance ratio, and the inverted driver resonance angular frequency

    Returns
    -------
    Fb : float
        Enclosure tuning frequency
    Qt : float
        Total system Q
    Tb : float
        Inverted system tuning angual frequency
    Ts : float
        Inverted driver resonance angular frequency
    a : float
        System compliance ratio
    h : float
        System tuning ratio
    """
    Vb_ci = 1000/(2.54**3)*Vb # Convert enclosure volume from liters to in.^3
    R = D/2 # Convert vent diameter to radius

#    Fb = np.sqrt((1.463e7*R)/(Vb_ci * (Lv + 1.463*R))) # Bad correction factor
    Fb = np.sqrt((1.464e7*R)/(Vb_ci*Lv)) # Fb from port. Needs volume in in.^3
    Tb = 1/(2*np.pi*Fb) # Inverse of angular frequencies
    Ts = 1/(2*np.pi*Fs)

    h = Fb/Fs # Tuning ratio
    a = Vas/Vb # Compliance ratio

    Qt = Qes*Qms/(Qes + Qms) # Approximate total system Q (Qt) with total driver Q (Qts)

    return Fb,Qt,Tb,Ts,a,h

def params_align(Fs,Qes,Qms,Vas,D,a,h):
    """
    Calculates dependent system parameters: enclosure tuning frequency, system tuning ratio, system compliance ratio, and the inverted driver resonance angular frequency

    Returns
    -------
    Fb : float
        Enclosure tuning frequency
    Lv : float
        Vent length, in units of vent diameter
    Qt : float
        Total system Q
    Tb : float
        Inverted system tuning angual frequency
    Ts : float
        Inverted driver resonance angular frequency
    """

    R = D/2 # Convert vent diameter to radius

    Fb = h*Fs # Calculate enclosure tuning frequency from tuning ratio
    Tb = 1/(2*np.pi*Fb) # Inverse of angular frequencies
    Ts = 1/(2*np.pi*Fs)

    Vb = Vas/a # Compliance ratio

    Vb_ci = 1000/(2.54**3)*Vb # Convert enclosure volume from liters to in.^3
    Lv = (1.464e7*R)/(Vb_ci*Fb**2) # Vent length calculation, in units of R

    Qt = Qes*Qms/(Qes + Qms) # Approximate total system Q (Qt) with total driver Q (Qts)

    return Fb,Lv,Qt,Tb,Ts

def alignment_spec(a,h,Ql,Qt,Ts):
    """    
    Calculates alignment specification parameters T0, a1, a2, and a3

    Returns
    -------
    T0 : float
    a1 : float
    a2 : float
    a3 : float
    """
    T0 = Ts/np.sqrt(h)
    a1 = (Ql + h*Qt)/(np.sqrt(h)*Ql*Qt)
    a2 = (h + (a + 1 + h**2)*Ql*Qt)/(h*Ql*Qt)
    a3 = (h*Ql + Qt)/(np.sqrt(h)*Ql*Qt)

    return T0,a1,a2,a3

def response(f,T0,a1,a2,a3):
    """
    System frequency response function, equation (20) in 'Vented Box Loudspeaker Systems'

    Returns
    -------
    dB : float
        Relative system response gain, in decibels
    """
    s = complex(0,2*np.pi*f)

    num = s**4 * T0**4
    denom = s**4 * T0**4 + a1 * s**3 * T0**3 + a2 * s**2 * T0**2 + a3 * s * T0 + 1
    gain = np.sqrt((num/denom).real**2 + (num/denom).imag**2)
    dB = 10*np.log10(gain)

    return dB

def displacement(f,a,Ql,Qt,Tb,Ts):
    """
    Diaphragm displacement function, equation (14) in 'Vented Box Loudspeaker Systems'

    Returns
    -------
    dB : float
        Normalized displacement value
    """

    s = complex(0,2*np.pi*f)

    num = s**2*Tb**2 + s*Tb/Ql + 1
    denom = s**4*Tb**2*Ts**2 + s**3*(Tb**2*Ts/Qt + Tb*Ts**2/Ql) + s**2*((a + 1)*Tb**2 + (Tb*Ts)/(Ql*Qt) + Ts**2) + s*(Tb/Ql + Ts/Qt) + 1

    displacement = np.sqrt((num/denom).real**2 + (num/denom).imag**2)

    return displacement

def impedance(f,a,Ql,Qes,Qms,Re,Tb,Ts):
    """
    Voice coil impedance function, equation (16) in 'Vented Box Loudspeaker Systems'

    Returns
    -------
    impedance : float
        Voice coil impedance value, in ohms
    """

    s = complex(0,2*np.pi*f)

    num = s*(Ts/Qms)*(s**2*Tb**2 + s*Tb/Ql + 1)
    denom = s**4*Tb**2*Ts**2 + s**3*(Tb**2*Ts/Qms + Tb*Ts**2/Ql) + s**2*((a + 1)*Tb**2 + (Tb*Ts)/(Ql*Qms) + Ts**2) + s*(Tb/Ql + Ts/Qms) + 1

    Res = float(Re)*(Qms/Qes)

    impedance = float(Re) + Res*np.sqrt((num/denom).real**2 + (num/denom).imag**2)

    return impedance

def xtickmarks(xmin,xmax):
    """
    Return x-axis tickmark location and label arrays based on x-axis range

    Returns
    -------
    loc : array-like
        tickmark locations
    labels : array-like
        tickmark labels
    """
    # Locations and labels, ordered
    loc = [10,20,50,100,200,500,1000,2000,5000,10000,20000]
    labels = [10,20,50,100,200,500,'1k','2k','5k','10k','20k']

    # Default index min and max cover whole range
    imin = 0
    imax = len(loc) - 1

    # Find the largest location index not greater than given x-axis minimum
    for i in range(0,11):
        if xmin < loc[i]:
            if i > 0:
                imin = i - 1
                break

    # Find the first index location equal or greater than given x-axis maximum
    for i in range(0,11):
        if xmax <= loc[i]:
            imax = i
            break

    # Return location and label array slice
    return loc[imin:imax + 1],labels[imin:imax + 1]

def plotter(frequencies,values,ylabel,title,pyplot,fig,basex=10,label=None,suptitle=None,ymin=None,ymax=None,ystep=None):
    """
    Generalized plotting function. Standardizes x-axis tickmark values, range and a few other things.
    """
    # Figure for plot
    pyplot.figure(fig)
    if label is not None:
        pyplot.semilogx(frequencies,values,label=label,basex=10)
        pyplot.legend()
    else:
        pyplot.semilogx(frequencies,values,basex=10)

    # Figure titles
    if suptitle is not None:
        pyplot.suptitle(suptitle)
    pyplot.title(title)

    pyplot.grid(b=True,which='minor')
    pyplot.grid(b=True,which='major')

    # Y-axis properties
    if ymin is not None and ymax is not None:
        pyplot.ylim(ymin,ymax)
        if ystep is not None:
            pyplot.yticks(np.arange(ymin,ymax+ystep,ystep))
    pyplot.ylabel(str(ylabel))

    # X-axis properties
    pyplot.xlim(frequencies.min(),frequencies.max())
    loc,labels = xtickmarks(frequencies.min(),frequencies.max()) # Get tick locations and labels
    pyplot.xticks(loc,labels)
    pyplot.xlabel('Frequency (Hz.)')

def freq_vals(freq_min,freq_max,res,basex=10):
    """
    Returns logarithmically spaced frequencies, from freq_min to freq_min with res number of points in between
    
    Returns
    -------
    frequency_range : array-like
    """
    frequency_values = np.logspace(np.log10(freq_min),np.log10(freq_max),num=res,base=10)

    return frequency_values

def response_plot(Fs,Qes,Qms,Re,Vas,D,a=None,h=None,Lv=None,Vb=None,Ql=7,label=None,suptitle=None,freq_min=10,freq_max=20000,res=1000,pyplot=None,fig=1):
    """
    Calculates and returns a vented loudspeaker enclosure system's response (gain) values over the specified frequency range (default range is 10 Hz to 20 kHz).

    Returns
    -------
    frequencies : array-like
        Frequencies where system response calculated
    responses : array-like
        Relative system response gain values, calculated at frequencies in response_range
    impedance_values : array-like
    displacement_values : array-like
    """

    """
    Three cases:
        1) Fs, Qes, Qms, Vas, Vb, Lv, and D are user-specified -> create alignment accordingly starting with calculation of Fb from Lv and D
        2) Fs, Qes, Qms, Vas, and D are user-specified -> lookup relevant values of a and h according to specified alignment type (EX: SC4 and SBB4)
        3) Fs, Qes, Qms, Vas, and D are user-specified, along with a and h -> create alignment using a and h values specified
    """

    if a is None and h is None and Lv is None and Vb is None:
        Fb,Lv,Qt,Tb,Ts,a,h = params_lookup(Fs,Qes,Qms,Re,Vas,D)
    elif Lv is not None and Vb is not None:
        if a is None and h is None:
            # Generate dependent params via vent tuning equation
            Fb,Qt,Tb,Ts,a,h = params_calc(Fs,Qes,Qms,Vas,D,Vb,Lv)
        elif a is not None and h is not None:
            # Generate dependent params via specified a,h
            Fb,Lv,Qt,Tb,Ts = params_align(Fs,Qes,Qms,Vas,D,a,h)

    # Generate alignment specification parameters
    T0,a1,a2,a3 = alignment_spec(a,h,Ql,Qt,Ts)

    # Frequency response range
    frequencies = freq_vals(freq_min,freq_max,res)

    # Empty array for frequency response data calculated next
    responses = np.array([])

    # Generate frequency response values
    for f in frequencies:
        responses = np.append(responses,response(f,T0,a1,a2,a3))

    # Plot the response data if given matplotlib plotting object
    if pyplot != None:
        plotter(frequencies,responses,'Response Level (dB)','Frequency Response',pyplot,fig,label=label,suptitle=suptitle,ymin=-24,ymax=6,ystep=3)

    # Return system response data
    return frequencies,responses

def impedance_plot(Fs,Qes,Qms,Re,Vas,D,a=None,h=None,Lv=None,Vb=None,Ql=7,label=None,suptitle=None,freq_min=10,freq_max=20000,res=1000,pyplot=None,fig=2):
    """
    Calculates and returns a vented loudspeaker enclosure driver's impedance values over the specified frequency range (default range is 10 Hz to 20 kHz).

    Returns
    -------
    frequencies : array-like
        Frequencies where impedances calculated
    impedances : array-like
        Calculated vented enclosure driver voice coil impedances
    """

    """
    Three cases:
        1) Fs, Qes, Qms, Vas, Vb, Lv, and D are user-specified -> create alignment accordingly starting with calculation of Fb from Lv and D
        2) Fs, Qes, Qms, Vas, and D are user-specified -> lookup relevant values of a and h according to specified alignment type (EX: SC4 and SBB4)
        3) Fs, Qes, Qms, Vas, and D are user-specified, along with a and h -> create alignment using a and h values specified
    """

    if a is None and h is None and Lv is None and Vb is None:
        Fb,Lv,Qt,Tb,Ts,a,h = params_lookup(Fs,Qes,Qms,Re,Vas,D)
    elif Lv is not None and Vb is not None:
        if a is None and h is None:
            # Generate dependent params via vent tuning equation
            Fb,Qt,Tb,Ts,a,h = params_calc(Fs,Qes,Qms,Vas,D,Vb,Lv)
        elif a is not None and h is not None:
            # Generate dependent params via specified a,h
            Fb,Lv,Qt,Tb,Ts = params_align(Fs,Qes,Qms,Vas,D,a,h)

    # Frequency range
    frequencies = freq_vals(freq_min,freq_max,res)

    # Empty array for impedance data calculated next
    impedances = np.array([])

    # Generate frequency response values
    for f in frequencies:
        impedances = np.append(impedances,impedance(f,a,Ql,Qes,Qms,Re,Tb,Ts))

    # Plot the impedance data if given matplotlib plotting object
    if pyplot != None:
        plotter(frequencies,impedances,'Impedance (Ohms)','Voice Coil Impedance',pyplot,fig,label=label,suptitle=suptitle)

    return frequencies,impedances

def displacement_plot(Fs,Qes,Qms,Re,Vas,D,a=None,h=None,Lv=None,Vb=None,Ql=7,label=None,suptitle=None,freq_min=10,freq_max=20000,res=1000,pyplot=None,fig=3):
    """
    Calculates and returns a vented loudspeaker enclosure driver's impedance values over the specified frequency range (default range is 10 Hz to 20 kHz).

    Returns
    -------
    frequencies : array-like
        Frequencies where displacements calculated
    displacements : array-like
        Displacement function magnitude
    """

    """
    Three cases:
        1) Fs, Qes, Qms, Vas, Vb, Lv, and D are user-specified -> create alignment accordingly starting with calculation of Fb from Lv and D
        2) Fs, Qes, Qms, Vas, and D are user-specified -> lookup relevant values of a and h according to specified alignment type (EX: SC4 and SBB4)
        3) Fs, Qes, Qms, Vas, and D are user-specified, along with a and h -> create alignment using a and h values specified
    """

    if a is None and h is None and Lv is None and Vb is None:
        Fb,Lv,Qt,Tb,Ts,a,h = params_lookup(Fs,Qes,Qms,Re,Vas,D)
    elif Lv is not None and Vb is not None:
        if a is None and h is None:
            # Generate dependent params via vent tuning equation
            Fb,Qt,Tb,Ts,a,h = params_calc(Fs,Qes,Qms,Vas,D,Vb,Lv)
        elif a is not None and h is not None:
            # Generate dependent params via specified a,h
            Fb,Lv,Qt,Tb,Ts = params_align(Fs,Qes,Qms,Vas,D,a,h)

    # Frequency range
    frequencies = freq_vals(freq_min,freq_max,res)

    # Empty array for impedance data calculated next
    displacements = np.array([])

    # Generate frequency response values
    for f in frequencies:
        displacements = np.append(displacements,displacement(f,a,Ql,Qt,Tb,Ts))

    # Plot the impedance data if given matplotlib plotting object
    if pyplot != None:
        plotter(frequencies,displacements,'Displacement Function Magnitude','Diaphragm Displacement',pyplot,fig,label=label,suptitle=suptitle)

    return frequencies,displacements
