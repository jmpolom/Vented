#!/usr/bin/python

"""
vented_box.py provides functions that calculate the approximate (model predicted), theoretical frequency response of a vented box loudspeaker system. The script uses mathematical model analysis results from A. N. Theile and Richard Small, as published in 'Loudspeakers in Vented Boxes' and 'Vented Box Loudspeaker Systems', respectively.

Script Name:    vented_box.py
Author:         Jonathan Polom
Date Created:   8 October 2010
License:        GNU General Public License
"""

import datetime
import numpy as np

__author__ = 'Jonathan Polom <s0nic0nslaught@gmail.com>'
__date__ = datetime.date(2010, 10, 13)
__version__ = 0.3

def dependent_params(Fs,Qts,Vas,Vb,Lv,R):
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
    Ts : float
        Inverted driver resonance angular frequency
    """
    Vb_ci = Vb*12**3 # Convert enclosure volume from cubic feet to cubic inches
    Fb = np.sqrt((1.463e7*R)/(Vb_ci * (Lv + 1.463*R))) # Needs volume in in.^3
    h = Fb/Fs # Tuning ratio
    a = Vas/Vb # Compliance ratio
    Ts = (2*np.pi*Fs)**-1

    # Approximate total system Q (Qt) with total driver Q (Qts)
    Qt = Qts

    return a,Fb,h,Qt,Ts

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


def response_plot(Fs,Qts,Vas,Lv,R,Vb,Ql,freq_min=10,freq_max=20000,res=1000,plot=None):
    """
    Calculates and returns the vented loudspeaker enclosure system's response (gain) values over the specified frequency range (default range is 10 Hz to 20 kHz).

    Returns
    -------
    response_range : array-like
        Frequencies system response was calculated at
    response_values : array-like
        Relative system response gain values, calculated at frequencies in response_range
    """
    # Generate dependent parameters
    a,Fb,h,Qt,Ts = dependent_params(Fs,Qts,Vas,Vb,Lv,R)

    # Generate alignment specification parameters
    T0,a1,a2,a3 = alignment_spec(a,h,Ql,Qt,Ts)

    # Frequency response range
    response_range = np.logspace(np.log10(freq_min),np.log10(freq_max),num=res,base=10)

    # Empty array for frequency response data calculated next
    response_values = np.array([])

    # Generate frequency response values
    for f in response_range:
        response_values = np.append(response_values,response(f,T0,a1,a2,a3))

    # Plot the response data if given matplotlib plotting object
    if plot != None:
        plot.semilogx(response_range,response_values,basex=10)
        #plot.grid(b=True,which='both') # Doesn't display major and minor gridlines
        plot.grid(b=True,which='minor') # Apparently the only way to get minor and major
        plot.grid(b=True,which='major') # grid lines on the plot...

        # Title
        plot.title('Vented Enclosure Frequency Response')

        # Y-axis properties
        plot.ylim(-24,6)
        plot.yticks(range(-24,9,3))
        plot.ylabel('Response Level (dB)')

        # X-axis properties
        plot.xlim(freq_min,freq_max)
        loc,labels = xtickmarks(freq_min,freq_max) # Get tick locations and labels
        plot.xticks(loc,labels)
        plot.xlabel('Frequency (Hz.)')

    # Return system response data
    return response_range,response_values
