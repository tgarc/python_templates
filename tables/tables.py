import re, sys
import numpy as np
from collections import OrderedDict


class TextTable(object):
    def __init__(self, keys, fmt=None, delimiter='|', fill='-', padding=1
                 , linebreak=None):
        self._keys = tuple(keys)

        self.fill = fill
        self.padding = padding
        self.delimiter = delimiter.join([' '*padding]*2)
        self.padfill = fill.join([self.fill*self.padding]*2)
        self._linebreak = linebreak if linebreak is not None else self.fill

        self._content = []
        self._widths = np.zeros((1000,len(self._keys)),dtype=np.uint8)
        self.__index = 0
        self._updateWidths(keys)

        if isinstance(fmt,dict):
            self._fmt = [fmt[k] if k in fmt else "{!s:}" for k in self._keys]
        elif fmt is not None:
            self._fmt = list(fmt)
        else:
            self._fmt = []

        if len(self._fmt) < len(self._keys):
            self._fmt[len(self._fmt):] = ["{!s:}"]*(len(self._keys)-len(self._fmt))
        self._fmt = OrderedDict(zip(self._keys,self._fmt))

    def keys(self):
        return list(self._keys)

    def getRowDict(self,fill=None):
        return dict.fromkeys(self._keys, fill)

    def addLineBreak(self):
        self._content.append(self._linebreak)

    def _formatLine(self,*vals,**kwvals):
        line = [self._fmt[k].format(vals[i] if i < len(vals) else kwvals[k])
                if (i < len(vals) or k in kwvals) else "-"
                for i,k in enumerate(self._keys)]
        
        return line

    def _updateWidths(self,line):
        """ update the column width counter """
        if isinstance(line,str):
            return

        if self.__index >= len(self._widths):
            self._widths = np.r_[self._widths,np.zeros_like(self._widths)]
        self._widths[self.__index] = map(len,line[:len(self._keys)])
        self.__index += 1

    def addRow(self,*vals,**kwvals):
        line = self._formatLine(*vals,**kwvals)
        self._updateWidths(line)
        self._content.append(line)

    def addHeader(self,widths):
        header = []
        header.append(self.delimiter.join(k.center(widths[i])
                                          for i,k in enumerate(self._keys)))
        th_separator = '+'.join([self.fill*self.padding]*2)
        header.append(th_separator.join(self.fill*w for w in widths))

        self._content = header + self._content

    def addFooter(self,widths):
        return

    def _calculateWidths(self):
        return self._widths.max(0)

    def writeTable(self, out=sys.stdout):
        widths = self._calculateWidths()

        self.addHeader(widths)
        self.addFooter(widths)

        if self._linebreak == self.fill:
            linebreak = self.padfill.join(self.fill*w for w in widths)
        else:
            linebreak = self._linebreak

        for line in self._content:
            if line == self._linebreak:
                line = linebreak
            elif not isinstance(line, str):
                fline = self.delimiter.join([line[i].rjust(widths[i])
                                             for i in xrange(len(self._keys))])
                if len(line) > len(self._keys):
                    fline += str(line[len(self._keys)])
                line = fline
            # else:
            #     line = line
            print >> out, line


class LatexTable(TextTable):
    """
    A class which encapsulates crudely a latex table as a list of formatted
    strings.
    """
    _conv = { '&' : r'\&',
              '%' : r'\%',
              '$' : r'\$',
              '#' : r'\#',
              '_' : r'\_',
              '{' : r'\{',
              '}' : r'\}',
              '~' : r'\textasciitilde{}',
              '^' : r'\^{}',
              '\\': r'\textbackslash{}',
              '<' : r'\textless',
              '>' : r'\textgreater'}
    _regex = re.compile('|'.join(re.escape(unicode(key)) for key in _conv.keys()))
    _escape = lambda self, line: self._regex.sub(lambda m: self._conv[m.group()], line)
    _linebreak = r"\hline"
    _delimiter = "&"
    
    def __init__(self, keys, caption=None, label=None, **kwargs):
        kwargs['delimiter'] = kwargs.get('delimiter',self._delimiter)
        kwargs['linebreak'] = kwargs.get('linebreak',self._linebreak)
        super(self.__class__, self).__init__(keys, **kwargs)

        self.caption = caption
        self.label = label
        self.hdrFormat = ("|%c" % "r")*len(keys) +"|"

    def addHeader(self,*args,**kwargs):
        header = []
        header.append(r"\begin{longtable}{%s}" % self.hdrFormat)
        if self.caption:
            header.append(r"\caption{%s}" % self.caption)
        if self.label:
            header.append(self.label)
        header.append(r"\hline")

        rowheaders = self.delimiter.join(("\\textbf{{{}}}",)*len(self._keys)).format(*self._keys)
        header.append(rowheaders+r'\\')
        header.append(r"\hline")

        self._content = header + self._content

    def addFooter(self,*args,**kwargs):
        self._content += [r'\hline', r'\end{longtable}']

    def addRow(self, *vals,**kwvals):
        line = self._formatLine(*vals,**kwvals)
        line = map(self._escape, line) + [r'\\']
        self._updateWidths(line)
        self._content.append(line)
