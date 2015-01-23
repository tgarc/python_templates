#!/usr/bin/python

from operator import itemgetter

def finalizefigure(ax,fig):
    legdict = dict(zip(*reversed(ax.get_legend_handles_labels())))
    labels, handles = zip(*[(k,legdict[k]) for k in keys])
    leg = plt.legend(handles, labels
                     , bbox_to_anchor=(0, 1.015, 1., .05)
                     , mode='expand', ncol=len(legdict.keys())
                     , frameon=False, loc='top left')
    ax.add_artist(leg)

    ax.yaxis.set_ticks(np.arange(len(datasets))+1)
    ax.yaxis.set_ticklabels(('set 1', 'set 2', 'set 3'))

    fig.tight_layout()
    fig.subplots_adjust(top=0.925)


fig = plt.figure(figsize=(8.0,8.0))
ax = fig.add_subplot(111)

keys = ('x <= 0.5','0 <= x < 0.2','0.4 <= x < 0.6','0.8 <= x < 1')
codes = tuple(1 << i for i in range(len(keys)))
colors = ('#1a9641','#2c7bb6','#d7191c','#fdae61')

datasets = [np.random.rand(100) for i in range(3)]

for i,data in enumerate(datasets):
    mask = np.zeros(data.shape,np.int32)

    mask[data <= 0.5]                   |= codes[0]
    mask[(0 <= data) & (data < 0.2)]    |= codes[1]
    mask[(0.4 <= data) & (data < 0.6)]  |= codes[2]
    mask[(0.8 <= data) & (data < 1)]    |= codes[3]

    # plot_status(ax,mask,codes,keys,colors=colors,stack=True,yoffset=i)
    # plot_status(ax,mask,codes,keys,colors=colors,stack=False,yoffset=i)
    # plot_status(ax,mask,codes,keys=keys,colors=colors,yoffset=i
    #             ,stack=[False,True,True,True]
    #             ,alpha=[0.25,0.9,0.9,0.9],linewidth=0.5)
    plot_status(ax,mask,codes,keys=keys,colors=colors,yoffset=i
                ,stack=False
                ,alpha=[0.25,1,1,1]
                ,relsize=[1,0.8,0.6,0.4],linewidth=0.5)

finalizefigure(ax,fig)
