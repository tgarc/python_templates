# python_templates
a set of python implementations for reference

### plotmask

Plot values of a 1d integer valued array as horizontal bars that are vertically
stacked, overlayed, or some combination thereof, using matplotlib's barh plot.

Example
```
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
      plt.setp(ax.xaxis.get_minorticklabels(),rotation=90,fontsize=9)
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
isDate=False
for i,prn in enumerate(sorted(datasets)):
    r = np.random.randint(len(datasets[prn]))
    datasets[prn] |= codes[0]
    datasets[prn][r:r+200+np.random.randint(-50,51)] ^= codes[0]

    # Example 3
    plotmask(ax, datasets[prn], codes, delta=timedelta(seconds=1), xoffset=t0
             , labels=keys, colors=colors, yoffset=i, stack=False
             , alpha=[0.65]+[0.8]*3, relsize=[1]+[0.5]*3, height=0.8)
finalizefigure(ax,fig,isDate=True)
```

Output

![Example](plotmask/ex3.png?raw=true "Example")

### tables

Classes for writing simple text and latex tables with automatic column width calculation.

Example
```
keys = ('x <= 0.5','0 <= x < 0.2','0.4 <= x < 0.6','0.8 <= x < 1')
codes = tuple(1 << i for i in range(len(keys)))

table = LatexTable(keys,caption='Array histograms')
sumrows = [0]*len(keys) # keep track of sums of rows
for i,data in enumerate(datasets):
    mask = np.zeros(data.shape,np.int32)
    mask[data <= 0.5]                   |= codes[0]
    mask[(0 <= data) & (data < 0.2)]    |= codes[1]
    mask[(0.4 <= data) & (data < 0.6)]  |= codes[2]
    mask[(0.8 <= data) & (data < 1)]    |= codes[3]

    for i in range(len(sumrows)):
        sumrows[i] += sum((mask & codes[i]) == codes[i])
    table.addRow(*map(lambda c: sum((mask & c) == c), codes))
table.addLineBreak()
table.addRow(*sumrows)
table.writeTable(out=sys.stdout)
```

Output
```
\begin{longtable}{|r|r|r|r|}
\caption{Array histograms}
\hline
\textbf{x <= 0.5} & \textbf{0 <= x < 0.2} & \textbf{0.4 <= x < 0.6} & \textbf{0.8 <= x < 1}\\
\hline
      55 &           28 &             20 &           22\\
      58 &           17 &             18 &           17\\
      47 &           24 &             18 &           17\\
\hline
     160 &           69 &             56 &           56\\
\hline
\end{longtable}
```
