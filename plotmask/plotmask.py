import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import dates
from datetime import datetime, timedelta

def mask2runs(mask):
  '''
  mask2runs

  Turns a boolean array into an array of indices indicating the start and end
  indices of each run of 1s
  '''
  runflip = np.nonzero(mask[1:] ^ mask[:-1])[0]+1

  # if the last run was a run of 1's, count the last epoch as the end of the run
  # similarly if the first bit started a run of 1's
  if mask[-1]: runflip = np.r_[runflip,len(mask)]
  if mask[0]: runflip = np.r_[0,runflip]

  return runflip[::2],runflip[1::2]


def plotmask(ax, data, values, xoffset=0, yoffset=1, delta=1, labels=None
             , colors=None , alpha=1, relsize=1.0, stack=True , height=0.9
             , linewidth=1, edgecolor='0.5'):
    '''
    plotmask

    Plot values of a 1d integer valued array as horizontal bars that are vertically
    stacked, overlayed, or some combination thereof, using matplotlib's barh
    plot.

    *data*
    1d array of integer values.

    *values*
    Values of array to plot.

    *xoffset*
    Left x limit of data.

    *yoffset*
    y value of center of bar plot.

    *delta*
    Fixed interval assumed between datapoints.
    
    *labels*
    Text label descriptors for *values*.

    *colors*
    Distinct colors applied to each of *values*.

    *alpha*
    The alpha transparency value. Can also be specified as an array to set an
    alpha for each value.

    *relsize*
    Specify scale of horizontal bars relative to *height* and their allocated
    area (should be <= 1). Can be specified as an array to set relative height
    for each value. If stacked, the bar height will be height*relsize/sum(stack)
    otherwise, bar height will be height*relsize.
    
    *stack*
    Specify whether to vertically stack horizontal bars or allow them to
    overlap. Two ways to set this parameter:
    1)  Setting stack to True plots all *values* of an array vertically stacked
        on top of each other. When False, all values will be displayed at the
        same size and location (this works well if the status values are
        mutually exclusive).
    2)  Alternately, you can specify an array of booleans that determines which
        values are stacked and which aren't. This is useful if you want to allow
        some *values* to overlap.

    '''
    if colors is None:
        colors = tuple(mpl.cm.jet(i/float(len(values))) for i in range(len(values)))
    if labels is None:
        labels = map(str,values)

    if not hasattr(stack,'__iter__'):   stack = [stack]*len(values)
    if not hasattr(alpha,'__iter__'):   alpha = [alpha]*len(values)
    if not hasattr(relsize,'__iter__'): relsize = [relsize]*len(values)

    if isinstance(xoffset,datetime) and isinstance(delta,timedelta):
        ax.xaxis_date()
        delta = dates.date2num(xoffset+delta) - dates.date2num(xoffset)
        xoffset = dates.date2num(xoffset)

    j = 0
    for s,v,a,scale,l,c in zip(stack,values,alpha,relsize,labels,colors):
        codeheight = height*scale
        codeoffset = yoffset+1 - codeheight/2.

        if s:
            codeheight /= sum(stack)
            codeoffset += codeheight*float(j)
            j += 1

        starts, stops = mask2runs((data & v) == v)
        starts = np.array(starts)*delta+xoffset
        stops = np.array(stops)*delta+xoffset

        ax.barh([codeoffset]*len(starts), np.subtract(stops,starts), left=starts
                , height=codeheight, linewidth=0, align='edge'
                , alpha=a, label=l, color=c)
    ax.barh(yoffset+1, len(data)*delta, left=xoffset
            , height=height*max(relsize)
            , align='center', color="None"
            , edgecolor=edgecolor, linewidth=linewidth)
