# python_templates
a set of python implementations for reference

### plotmask

Plot values of a 1d integer valued array as horizontal bars that are vertically
stacked, overlayed, or some combination thereof, using matplotlib's barh plot.

Example:
![Example](plotmask/ex3.png?raw=true "Example")

### tables

Classes for writing simple text and latex tables with automatic column width calculation.

Example:
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
