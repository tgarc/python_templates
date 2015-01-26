from tables import TextTable,LatexTable
import numpy as np
import sys

keys = ('x <= 0.5','0 <= x < 0.2','0.4 <= x < 0.6','0.8 <= x < 1')
codes = tuple(1 << i for i in range(len(keys)))
datasets = [np.random.rand(100) for i in range(3)]

# table = TextTable(keys)
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
