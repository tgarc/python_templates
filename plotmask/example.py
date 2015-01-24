#!/usr/bin/python

from matplotlib import dates, ticker
from datetime import timedelta,datetime
from codeplotter import plot_status


def finalizefigure(ax,fig,isDate=False):
    legdict = dict(zip(*reversed(ax.get_legend_handles_labels())))
    labels, handles = zip(*[(k,legdict[k]) for k in keys])
    leg = plt.legend(handles, labels
                     , bbox_to_anchor=(0, 1.015, 1., .05)
                     , mode='expand', ncol=len(legdict.keys())
                     , frameon=False, loc='upper left')
    ax.add_artist(leg)

    ax.yaxis.set_ticks(np.arange(len(datasets))+1)
    ax.yaxis.set_ticklabels(sorted(datasets))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())    

    if isDate:
      ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
      ax.xaxis.set_minor_formatter(dates.DateFormatter(':%M:'))
      plt.setp(ax.xaxis.get_minorticklabels(),rotation=90)
      fig.autofmt_xdate()
    
    ax.set_ylim(0.5,len(datasets)+0.5)
    fig.tight_layout()
    fig.subplots_adjust(top=0.925)

fig = plt.figure(figsize=(8.0,8.0))
ax = fig.add_subplot(111)

keys = ("Ephemeris", "Code Carrier Lock", "Unlocked", "Code Lock")
codes = tuple(1 << i for i in range(len(keys)))
colors = ('#2c7bb6','#1a9641','#d7191c','#fdae61')

# Generate some made up data
t0 = datetime.now()
datasets = {'PRN%02d'%i:np.random.choice(codes[1:],1000,p=(0.85,0.1,0.05)) for i in range(1,11)}

for i,prn in enumerate(sorted(datasets)):
    datasets[prn] |= codes[0]

    # stacked
    # plot_status(ax,datasets[prn],codes,labels=keys,colors=colors,stack=True,yoffset=i)

    # overlayed & no relsize
    # plot_status(ax,datasets[prn],codes,labels=keys,colors=colors,stack=False,yoffset=i)

    # stack + overlay 
    # plot_status(ax,datasets[prn],codes,labels=keys,colors=colors,yoffset=i
    #             ,stack=[False,True,True,True]
    #             ,alpha=[0.25,0.9,0.9,0.9])

    # stack + overlay + relsize
    # plot_status(ax, datasets[prn], codes , labels=keys, colors=colors,
    #             yoffset=i, stack=False , alpha=[0.65]+[0.9]*3,
    #             relsize=[1]+[0.9]*3, height=0.9)

    # stack + overlay + relsize + time
    # plot_status(ax, datasets[prn], codes, delta=timedelta(seconds=1), xoffset=t0
    #             , labels=keys, colors=colors, yoffset=i, stack=False
    #             , alpha=[0.65]+[0.9]*3, relsize=[1]+[0.9]*3, height=0.9)
finalizefigure(ax,fig,isDate=False)
