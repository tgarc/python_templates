#!/usr/bin/python

'''
data availability analysis for python!

Thomas J Garcia

Note:
+ GPS system is assumed!
+ demodulatorStatus code 255 is ignored here (although it appears in the data) because it is undocumented as far as I know
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import dates, ticker
from datetime import datetime, timedelta

def mask2runs(mask):
  '''
  mask2runs

  Turns a boolean array into an array of indices indicating the start and end indices of each run of 1s
  '''
  runflip = np.nonzero(mask[1:] ^ mask[:-1])[0]+1

  # if the last run was a run of 1's, count the last epoch as the end of the run
  # similarly if the first bit started a run of 1's
  if mask[-1]: runflip = np.r_[runflip,len(mask)]
  if mask[0]: runflip = np.r_[0,runflip]

  return runflip[::2],runflip[1::2]


def mask2nums(mask, startTime, cadence, timeDenominator, xoffset=0):
    runstarts, runstops = mask2runs(mask)
    runstarts = [dates.date2num(timedelta(seconds=(rstart+xoffset)*cadence/timeDenominator)+startTime) for rstart in runstarts]
    runstops = [dates.date2num(timedelta(seconds=(rstop+xoffset)*cadence/timeDenominator)+startTime) for rstop in runstops]

    return runstarts, runstops


def plot_status(ax, data, times=None, values=None, colors=None, title=None
                , alpha=1.0, relsize=1.0, stack=True, barHeight=0.95):
    '''
    plot_status

    Plot a collection of 1d integer valued arrays as horizontal bars (using
    matplotlib's barh plot).  Plot style is defined by the keys, values, and
    colors.

    *datadict*
    Dictionary of dictionary items each containing a 'mask' and 'times' array.
    
    *keys*
    Text label descriptors for status values.

    *colors*
    Distinct colors applied to each status value.

    *values*
    Distinct integer values corresponding to each status.

    *stack*
    Two ways to set this parameter:    
    1)  Setting stack to True makes plot so that all unique values of an array
        are stacked on top of each other. When False, all values will be
        displayed at the same size and location (this works well if the status
        values are mutually exclusive).
    2)  Alternately, you can specify an array of booleans that determines which
        values are stacked and which aren't. This can be useful for example if you
        want to plot one value underneath a stack of different values.

    *relsize*
    Two ways to set this parameter:
    1)  A single value sets vertical size of _all_ status bars relative to 1.0
        (should be <= 1).
    2)  An array of values sets relative vertical size for each individual status
        value. Using this option together with the *stack* option gives some
        additional flexibility, allowing for example a dominant underlying
        status value with a stack of smaller sized status value bars fitting
        within it.
    
    *alpha*
    The alpha value to use for all status codes. Can also be specified as an
    array to set a unique alpha for each status code.

    *vspace*
    Amount of vertical space (relative to 1) to leave between horizontal bars.

    '''
    if colors is None:
      colors = tuple(mpl.cm.jet(i/float(len(values))) for i in range(len(values)))

    if not hasattr(stack,'__iter__'):   stack = [stack]*len(values)
    if not hasattr(alpha,'__iter__'):   alpha = [alpha]*len(values)
    if not hasattr(relsize,'__iter__'): relsize = [relsize]*len(values)

    # plot each status code
    j = 0
    for k,(s,v,a,scale) in enumerate(zip(stack,values,alpha,relsize)):
        codeheight = barHeight*scale
        codeoffset = 1 - codeheight/2.

        if s:
            codeheight /= sum(stack)
            codeoffset += codeheight*float(j)
            j += 1

        starts, stops = mask2runs(data & v)
        ax.barh([codeoffset]*len(starts), np.subtract(stops,starts), left=starts
                    , height=codeheight
                    , color=colors[k], linewidth=0, align='edge', alpha=a)
