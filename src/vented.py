#!/usr/bin/python

"""
vented loudspeaker enclosure response calculator
Copyright (C) 2010 Jonathan Polom <s0nic0nslaught@gmail.com>

vented.py provides functions that calculate model predicted frequency 
response of a vented box loudspeaker system. This script uses electro-
acoustic model analysis results from A. N. Theile and Richard Small, 
as published in 'Loudspeakers in Vented Boxes' and 'Vented Box 
Loudspeaker Systems', respectively.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import datetime
from numpy import append,array,float,log10,logspace,pi,sqrt

__author__ = "Jonathan Polom <s0nic0nslaught@gmail.com>"
__date__ = datetime.date(2011, 3, 9)
__version__ = '0.8'
__license__ = "GNU General Public License"

class Enclosure():
    def __init__(self, Vb, Lv, D, base='custom'):
        self.Vb = float(Vb)
        self.Lv = float(Lv)
        self.D = float(D)

        try:
            baseChoices = ['si','us', 'custom']
            i = baseChoices.index(base)
            self.base = baseChoices[i]
        except:
            raise ValueError

class Driver():
    def __init__(self, Fs, Qes, Qms, Vas, base='si'):
        self.Fs = float(Fs)
        self.Qes = float(Qes)
        self.Qms = float(Qms)
        self.Vas = float(Vas)
        
        try:
            baseChoices = ['si','us']
            i = baseChoices.index(base)
            self.base = baseChoices[i]
        except:
            raise ValueError

class System():
    def __init__(self, driver, enclosure, a=None, h=None, resp='calculated'):
        self.driver = driver
        self.enclosure = enclosure

        try:
            rChoices = ['qb3','sqb3','bb4','sbb4','c4','sc4','calculated']
            i = rChoices.index(resp)
            self.resp = rChoices[i]
        except:
            raise ValueError
        
        if a and h:
            self.a = float(a)
            self.h = float(h)

def filterCoeffs(a,h,Ql,Qt):
    """    
    Calculates filter coefficients a1, a2, and a3 from alignment parameters

    Returns
    -------
    a1 : float
    a2 : float
    a3 : float
    """
    a1 = (Ql + h*Qt)/(sqrt(h)*Ql*Qt)
    a2 = (h + (a + 1 + h**2)*Ql*Qt)/(h*Ql*Qt)
    a3 = (h*Ql + Qt)/(sqrt(h)*Ql*Qt)

    return a1,a2,a3

def timeConstant(h,Ts):
    """
    Calculates filter time constant T0 from tuning ratio and angular period
    of frequency driver free air resonance
    
    Returns
    -------
    T0 : float
    """
    T0 = Ts/sqrt(h)

    return T0

def paramsCalc(Fs,Qes,Qms,Vas,D,Vb,Lv):
    """
    Calculates dependent system parameters: enclosure tuning frequency, 
    system tuning ratio, system compliance ratio, and the inverted driver 
    resonance angular frequency

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

#    Fb = sqrt((1.463e7*R)/(Vb_ci * (Lv + 1.463*R))) # Bad correction factor
    Fb = sqrt((1.464e7*R)/(Vb_ci*Lv)) # Fb from port. Needs volume in in.^3
    Tb = 1/(2*pi*Fb) # Inverse of angular frequencies
    Ts = 1/(2*pi*Fs)

    h = Fb/Fs # Tuning ratio
    a = Vas/Vb # Compliance ratio

    Qt = Qes*Qms/(Qes + Qms) # Approximate total system Q (Qt) with total driver Q (Qts)

    return Fb,Qt,Tb,Ts,a,h

def paramsAlign(Fs,Qes,Qms,Vas,D,a,h):
    """
    Calculates dependent system parameters: enclosure tuning frequency, 
    system tuning ratio, system compliance ratio, and the inverted driver 
    resonance angular frequency

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
    Tb = 1/(2*pi*Fb) # Inverse of angular frequencies
    Ts = 1/(2*pi*Fs)

    Vb = Vas/a # Compliance ratio

    Vb_ci = 1000/(2.54**3)*Vb # Convert enclosure volume from liters to in.^3
    Lv = (1.464e7*R)/(Vb_ci*Fb**2) # Vent length calculation, in units of R

    Qt = Qes*Qms/(Qes + Qms) # Approximate total system Q (Qt) with total driver Q (Qts)

    return Fb,Lv,Qt,Tb,Ts

def paramsLookup(Fs,Qes,Qms,Vas,D,resp='qb3'):
    return NotImplementedError

def response(f,T0,a1,a2,a3):
    """
    System frequency response function, equation (20) in 'Vented Box 
    Loudspeaker Systems'

    Returns
    -------
    dB : float
        Relative system response gain, in decibels
    """
    s = complex(0,2*pi*f)

    num = s**4 * T0**4
    denom = s**4 * T0**4 + a1 * s**3 * T0**3 + a2 * s**2 * T0**2 + a3 * s * T0 + 1
    gain = sqrt((num/denom).real**2 + (num/denom).imag**2)
    dB = 10*log10(gain)

    return dB

def displacement(f,a,Ql,Qt,Tb,Ts):
    """
    Diaphragm displacement function, equation (14) in 'Vented Box 
    Loudspeaker Systems'

    Returns
    -------
    dB : float
        Normalized displacement value
    """

    s = complex(0,2*pi*f)

    num = s**2*Tb**2 + s*Tb/Ql + 1
    denom = s**4*Tb**2*Ts**2 + s**3*(Tb**2*Ts/Qt + Tb*Ts**2/Ql) + s**2*((a + 1)*Tb**2 + (Tb*Ts)/(Ql*Qt) + Ts**2) + s*(Tb/Ql + Ts/Qt) + 1

    displacement = sqrt((num/denom).real**2 + (num/denom).imag**2)

    return displacement

def impedance(f,a,Ql,Qes,Qms,Re,Tb,Ts):
    """
    Voice coil impedance function, equation (16) in 'Vented Box 
    Loudspeaker Systems'

    Returns
    -------
    impedance : float
        Voice coil impedance value, in ohms
    """

    s = complex(0,2*pi*f)

    num = s*(Ts/Qms)*(s**2*Tb**2 + s*Tb/Ql + 1)
    denom = s**4*Tb**2*Ts**2 + s**3*(Tb**2*Ts/Qms + Tb*Ts**2/Ql) + s**2*((a + 1)*Tb**2 + (Tb*Ts)/(Ql*Qms) + Ts**2) + s*(Tb/Ql + Ts/Qms) + 1

    Res = float(Re)*(Qms/Qes)

    impedance = float(Re) + Res*sqrt((num/denom).real**2 + (num/denom).imag**2)

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

def plotter(frequencies,values,ylabel,title,pyplot,fig,basex=10,
            label=None,suptitle=None,ymin=None,ymax=None,ystep=None):
    """
    Generalized plotting function. Standardizes x-axis tickmark values, 
    range and a few other things.
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
            pyplot.yticks(range(ymin,ymax,ystep))
    pyplot.ylabel(str(ylabel))

    # X-axis properties
    pyplot.xlim(frequencies.min(),frequencies.max())
    loc,labels = xtickmarks(frequencies.min(),frequencies.max()) # Get tick locations and labels
    pyplot.xticks(loc,labels)
    pyplot.xlabel('Frequency (Hz.)')

def freqVals(freq_min,freq_max,res,basex=10):
    """
    Returns logarithmically spaced frequencies, from freq_min to freq_min 
    with res number of points in between
    
    Returns
    -------
    frequency_range : array-like
    """
    frequency_values = logspace(log10(freq_min),log10(freq_max),num=res,base=10)

    return frequency_values

def responsePlot(Fs,Qes,Qms,Re,Vas,D,a=None,h=None,Lv=None,Vb=None,Ql=7,
                 resp='qb3',label=None,suptitle=None,freq_min=10,
                 freq_max=20000,res=1000,pyplot=None,fig=1):
    """
    Calculates and returns a vented loudspeaker enclosure system's 
    response (gain) values over the specified frequency range (default 
    range is 10 Hz to 20 kHz).

    Returns
    -------
    frequencies : array-like
        Frequencies where system response calculated
    responses : array-like
        Relative system response gain values, calculated at frequencies in 
        response_range
    impedance_values : array-like
    displacement_values : array-like
    """

    """
    Three cases:
        1) Fs, Qes, Qms, Vas, Vb, Lv, and D are specified -> create alignment accordingly starting with calculation of Fb from Lv and D
        2) Fs, Qes, Qms, Vas, and D are specified -> lookup relevant values of a and h according to specified alignment type (EX: SC4 and SBB4)
        3) Fs, Qes, Qms, Vas, and D are specified, along with a and h -> create alignment using a and h values specified
    """

    if a is None and h is None and Lv is None and Vb is None:
        Fb,Lv,Qt,Tb,Ts,a,h = params_lookup(Fs,Qes,Qms,Vas,D,resp=resp)
    elif Lv is not None and Vb is not None:
        if a is None and h is None:
            # Generate dependent params via vent tuning equation
            Fb,Qt,Tb,Ts,a,h = paramsCalc(Fs,Qes,Qms,Vas,D,Vb,Lv)
        elif a is not None and h is not None:
            # Generate dependent params via specified a,h
            Fb,Lv,Qt,Tb,Ts = paramsAlign(Fs,Qes,Qms,Vas,D,a,h)

    # Generate alignment specification parameters
    T0 = timeConstant(h,Ts)
    a1,a2,a3 = filterCoeffs(a,h,Ql,Qt)

    # Frequency response range
    frequencies = freqVals(freq_min,freq_max,res)

    # Empty array for frequency response data calculated next
    responses = array([])

    # Generate frequency response values
    for f in frequencies:
        responses = append(responses,response(f,T0,a1,a2,a3))

    # Plot the response data if given matplotlib plotting object
    if pyplot != None:
        plotter(frequencies,responses,'Response Level (dB)',
                'Frequency Response',pyplot,fig,label=label,suptitle=suptitle,
                ymin=-24,ymax=6,ystep=3)

    # Return system response data
    return frequencies,responses

def impedancePlot(Fs,Qes,Qms,Re,Vas,D,a=None,h=None,Lv=None,
                  Vb=None,Ql=7,label=None,suptitle=None,freq_min=10,
                  freq_max=20000,res=1000,pyplot=None,fig=2):
    """
    Calculates and returns a vented loudspeaker enclosure driver's 
    impedance values over the specified frequency range (default range is 
    10 Hz to 20 kHz).

    Returns
    -------
    frequencies : array-like
        Frequencies where impedances calculated
    impedances : array-like
        Calculated vented enclosure driver voice coil impedances
    """

    """
    Three cases:
        1) Fs, Qes, Qms, Vas, Vb, Lv, and D are specified -> create alignment accordingly starting with calculation of Fb from Lv and D
        2) Fs, Qes, Qms, Vas, and D are specified -> lookup relevant values of a and h according to specified alignment type (EX: SC4 and SBB4)
        3) Fs, Qes, Qms, Vas, and D are specified, along with a and h -> create alignment using a and h values specified
    """

    if a is None and h is None and Lv is None and Vb is None:
        Fb,Lv,Qt,Tb,Ts,a,h = params_lookup(Fs,Qes,Qms,Vas,D,resp=resp)
    elif Lv is not None and Vb is not None:
        if a is None and h is None:
            # Generate dependent params via vent tuning equation
            Fb,Qt,Tb,Ts,a,h = paramsCalc(Fs,Qes,Qms,Vas,D,Vb,Lv)
        elif a is not None and h is not None:
            # Generate dependent params via specified a,h
            Fb,Lv,Qt,Tb,Ts = paramsAlign(Fs,Qes,Qms,Vas,D,a,h)

    # Frequency range
    frequencies = freqVals(freq_min,freq_max,res)

    # Empty array for impedance data calculated next
    impedances = array([])

    # Generate frequency response values
    for f in frequencies:
        impedances = append(impedances,impedance(f,a,Ql,Qes,Qms,Re,Tb,Ts))

    # Plot the impedance data if given matplotlib plotting object
    if pyplot != None:
        plotter(frequencies,impedances,'Impedance (Ohms)',
                'Voice Coil Impedance',pyplot,fig,label=label,suptitle=suptitle)

    return frequencies,impedances

def displacementPlot(Fs,Qes,Qms,Re,Vas,D,a=None,h=None,
                     Lv=None,Vb=None,Ql=7,label=None,suptitle=None,
                     freq_min=10,freq_max=20000,res=1000,pyplot=None,fig=3):
    """
    Calculates and returns a vented loudspeaker enclosure driver's 
    impedance values over the specified frequency range (default range is 
    10 Hz to 20 kHz).

    Returns
    -------
    frequencies : array-like
        Frequencies where displacements calculated
    displacements : array-like
        Displacement function magnitude
    """

    """
    Three cases:
        1) Fs, Qes, Qms, Vas, Vb, Lv, and D are specified -> create alignment accordingly starting with calculation of Fb from Lv and D
        2) Fs, Qes, Qms, Vas, and D are specified -> lookup relevant values of a and h according to specified alignment type (EX: SC4 and SBB4)
        3) Fs, Qes, Qms, Vas, and D are specified, along with a and h -> create alignment using a and h values specified
    """

    if a is None and h is None and Lv is None and Vb is None:
        Fb,Lv,Qt,Tb,Ts,a,h = params_lookup(Fs,Qes,Qms,Vas,D,resp=resp)
    elif Lv is not None and Vb is not None:
        if a is None and h is None:
            # Generate dependent params via vent tuning equation
            Fb,Qt,Tb,Ts,a,h = paramsCalc(Fs,Qes,Qms,Vas,D,Vb,Lv)
        elif a is not None and h is not None:
            # Generate dependent params via specified a,h
            Fb,Lv,Qt,Tb,Ts = paramsAlign(Fs,Qes,Qms,Vas,D,a,h)

    # Frequency range
    frequencies = freqVals(freq_min,freq_max,res)

    # Empty array for impedance data calculated next
    displacements = array([])

    # Generate frequency response values
    for f in frequencies:
        displacements = append(displacements,displacement(f,a,Ql,Qt,Tb,Ts))

    # Plot the impedance data if given matplotlib plotting object
    if pyplot != None:
        plotter(frequencies,displacements,
                'Displacement Function Magnitude','Diaphragm Displacement',
                pyplot,fig,label=label,suptitle=suptitle)

    return frequencies,displacements
