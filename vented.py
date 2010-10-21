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
__version__ = 0.6

def params_calc(Fs,Qes,Qms,Vas,D,Vb,Lv):
    """
    Calculates dependent system parameters: enclosure tuning frequency, system tuning ratio, system compliance ratio, and the inverted driver resonance angular frequency

    Returns
    -------
    a : float
        System compliance ratio
    Fb : float
        Enclosure tuning frequency
    h : float
        System tuning ratio
    Qt : float
        Total system Q
    Tb : float
        Inverted system tuning angual frequency
    Ts : float
        Inverted driver resonance angular frequency
    """
    Vb_ci = 1000/(2.54**3)*Vb # Convert enclosure volume from liters to in.^3
    R = D/2 # Convert vent diameter to radius

#    Fb = np.sqrt((1.463e7*R)/(Vb_ci * (Lv + 1.463*R))) # Needs volume in in.^3
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
#    dB = 10*np.log10(displacement)

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


def response_plot(Fs,Qes,Qms,Re,Vas,D,a=None,h=None,Lv=None,Vb=None,Ql=7,name='',freq_min=10,freq_max=20000,res=1000,plot=None):
    """
    Calculates and returns a vented loudspeaker enclosure system's response (gain) values over the specified frequency range (default range is 10 Hz to 20 kHz).

    Returns
    -------
    response_range : array-like
        Frequencies system response was calculated at
    response_values : array-like
        Relative system response gain values, calculated at frequencies in response_range
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
    response_range = np.logspace(np.log10(freq_min),np.log10(freq_max),num=res,base=10)

    # Empty array for frequency response data calculated next
    response_values = np.array([])
    """
    # Empty array for displacement data calculated next
    displacement_values = np.array([])

    # Empty array for impedance data calculated next
    impedance_values = np.array([])
    """
    # Generate frequency response values
    for f in response_range:
        response_values = np.append(response_values,response(f,T0,a1,a2,a3))
#        displacement_values = np.append(displacement_values,displacement(f,a,Ql,Qt,Tb,Ts))
#        impedance_values = np.append(impedance_values,impedance(f,a,Ql,Qes,Qms,Re,Tb,Ts))

    # Plot the response data if given matplotlib plotting object
    if plot != None:
        # Figure 1 for response plot
        plot.figure(1)
        plot.semilogx(response_range,response_values,basex=10)

        # Figure 1 title
        plot.suptitle(name)
        plot.title('Frequency Response')

        plot.grid(b=True,which='minor')
        plot.grid(b=True,which='major')

        # Frequency response Y-axis properties
        plot.ylim(-24,6)
        plot.yticks(range(-24,9,3))
        plot.ylabel('Response Level (dB)')

        # X-axis properties
        plot.xlim(freq_min,freq_max)
        loc,labels = xtickmarks(freq_min,freq_max) # Get tick locations and labels
        plot.xticks(loc,labels)
        plot.xlabel('Frequency (Hz.)')
    """
        # Create second Y-axis for driver displacement
#        plot2x = plot.twinx()
#        plot2x.semilogx(response_range,impedance_values,'k-',linewidth=0.5,basex=10)

        # Driver displacement Y-axis properties
#        plot2x.set_ylabel('Voice Coil Impedance (Ohms)')

        # X-axis labeling must come before twinx() called
#        plot.xlabel('Frequency (Hz.)')

        # Figure 2 for impedance plot
        plot.figure(2)
        plot.semilogx(response_range,impedance_values,basex=10)

        # Figure 1 title
        plot.suptitle(name)
        plot.title('Voice Coil Impedance')

        plot.grid(b=True,which='minor')
        plot.grid(b=True,which='major')

        # Frequency response Y-axis properties
        plot.ylabel('Impedance (Ohms)')

        # X-axis properties
        plot.xlim(freq_min,freq_max)
        loc,labels = xtickmarks(freq_min,freq_max) # Get tick locations and labels
        plot.xticks(loc,labels)
        plot.xlabel('Frequency (Hz.)')        

        # Figure 3 for diaphragm displacement plot
        plot.figure(3)
        plot.semilogx(response_range,displacement_values,basex=10)

        # Figure 3 title
        plot.suptitle(name)
        plot.title('Diaphragm Displacement')
        
        plot.grid(b=True,which='minor')
        plot.grid(b=True,which='major')
        
        # Y-axis properties
        plot.ylabel('Displacement Function Magnitude')
        
        # X-axis properties
        plot.xlim(freq_min,freq_max)
        loc,labels = xtickmarks(freq_min,freq_max) # Get tick locations and labels
        plot.xticks(loc,labels)
        plot.xlabel('Frequency (Hz.)')
    """
    # Return system response data
    return response_range,response_values

def impedance_plot(Fs,Qes,Qms,Re,Vas,D,a=None,h=None,Lv=None,Vb=None,Ql=7,name='',freq_min=10,freq_max=20000,res=1000,plot=None):
    """
    Calculates and returns a vented loudspeaker enclosure driver's impedance values over the specified frequency range (default range is 10 Hz to 20 kHz).

    Returns
    -------
    impedance_range : array-like
        Frequencies system response was calculated at
    impedance_values : array-like
        Relative system response gain values, calculated at frequencies in response_range
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
#    T0,a1,a2,a3 = alignment_spec(a,h,Ql,Qt,Ts)

    # Frequency range
    impedance_range = np.logspace(np.log10(freq_min),np.log10(freq_max),num=res,base=10)

    # Empty array for impedance data calculated next
    impedance_values = np.array([])

    # Generate frequency response values
    for f in impedance_range:
        impedance_values = np.append(impedance_values,impedance(f,a,Ql,Qes,Qms,Re,Tb,Ts))

    # Plot the impedance data if given matplotlib plotting object
    if plot != None:
        # Figure 2 for impedance plot
        plot.figure(2)
        plot.semilogx(impedance_range,impedance_values,basex=10)

        # Figure 2 title
        plot.suptitle(name)
        plot.title('Voice Coil Impedance')

        plot.grid(b=True,which='minor')
        plot.grid(b=True,which='major')

        # Frequency response Y-axis properties
        plot.ylabel('Impedance (Ohms)')

        # X-axis properties
        plot.xlim(freq_min,freq_max)
        loc,labels = xtickmarks(freq_min,freq_max) # Get tick locations and labels
        plot.xticks(loc,labels)
        plot.xlabel('Frequency (Hz.)')

    return impedance_range,impedance_values

def displacement_plot(Fs,Qes,Qms,Re,Vas,D,a=None,h=None,Lv=None,Vb=None,Ql=7,name='',freq_min=10,freq_max=20000,res=1000,plot=None):
    """
    Calculates and returns a vented loudspeaker enclosure driver's impedance values over the specified frequency range (default range is 10 Hz to 20 kHz).

    Returns
    -------
    displacement_range : array-like
        Frequencies system response was calculated at
    displacement_values : array-like
        Relative system response gain values, calculated at frequencies in response_range
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
#    T0,a1,a2,a3 = alignment_spec(a,h,Ql,Qt,Ts)

    # Frequency range
    displacement_range = np.logspace(np.log10(freq_min),np.log10(freq_max),num=res,base=10)

    # Empty array for impedance data calculated next
    displacement_values = np.array([])

    # Generate frequency response values
    for f in displacement_range:
        displacement_values = np.append(displacement_values,displacement(f,a,Ql,Qt,Tb,Ts))

    # Plot the impedance data if given matplotlib plotting object
    if plot != None:
        # Figure 3 for diaphragm displacement plot
        plot.figure(3)
        plot.semilogx(displacement_range,displacement_values,basex=10)

        # Figure 3 title
        plot.suptitle(name)
        plot.title('Diaphragm Displacement')
        
        plot.grid(b=True,which='minor')
        plot.grid(b=True,which='major')
        
        # Y-axis properties
        plot.ylabel('Displacement Function Magnitude')
        
        # X-axis properties
        plot.xlim(freq_min,freq_max)
        loc,labels = xtickmarks(freq_min,freq_max) # Get tick locations and labels
        plot.xticks(loc,labels)
        plot.xlabel('Frequency (Hz.)')

    return displacement_range,displacement_values
