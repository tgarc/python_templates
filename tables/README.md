# tables
Classes for writing simple text and latex tables with automatic column width
calculation.

## Text Example

```
x <= 0.5 | 0 <= x < 0.2 | 0.4 <= x < 0.6 | 0.8 <= x < 1
---------+--------------+----------------+-------------
      48 |           19 |             18 |           25
      40 |           16 |             22 |           19
      50 |           25 |             18 |           20
-------------------------------------------------------
     138 |           60 |             58 |           64
```

```
from tables import TextTable
import numpy as np

keys = ('x <= 0.5','0 <= x < 0.2','0.4 <= x < 0.6','0.8 <= x < 1')
codes = tuple(1 << i for i in range(len(keys)))
datasets = [np.random.rand(100) for i in range(3)]

table = TextTable(keys)
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

## Latex Example

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

```
from tables import LatexTable
import numpy as np

keys = ('x <= 0.5','0 <= x < 0.2','0.4 <= x < 0.6','0.8 <= x < 1')
codes = tuple(1 << i for i in range(len(keys)))
datasets = [np.random.rand(100) for i in range(3)]

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
