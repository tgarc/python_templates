import mdh
import sys
import calendar
from datetime import datetime
import gpstk
import numpy as np
import dateutil as du
from collections import OrderedDict
from itertools import groupby
import itertools as it
import operator as op
import re
from abc import ABCMeta, abstractmethod

DEFAULT_TIME_FMT =                "%Y-%m-%d %H:%M:%S.%f"

GPS_EPOCH_IN_UNIXTIME=        315964800

RangingCode2h5t = {"CA": 1, "P": 2, "Y": 3, "CodelessY":4,
                   "CM":5, "CL":6, "I5":7, "Q5":8,
                   "CMCL":9, "CP":10, "CD":11, "CPCD":12,
                   "M":13, "Mprime":15}
CarrierCode2h5t = {'L1':1,'L2':2,'L5':5}
h5t2RangingCode = {v:k for k,v in RangingCode2h5t.iteritems()}
h5t2CarrierCode = {v:k for k,v in CarrierCode2h5t.iteritems()}

Code2h5t = lambda cc,rc: (CarrierCode2h5t[cc],RangingCode2h5t[rc])
VERBOSE = 0


# ============================================================
# MDH helper tools
# ============================================================
def mdhSpan(m,absolute=False,**attrs):
    dsstart, dsstop = zip(*[(ds.attrs['startTime']/float(ds.attrs['timeDenominator'] if absolute else 1)
                             ,getStopTime(ds)/float(ds.attrs['timeDenominator'] if absolute else 1))
                             for ds in m.matchingDatasets(**attrs)])
    return min(dsstart), max(dsstop)

def getStopTime(dataset):
    if "cadence" in dataset.attrs:
        return dataset.attrs['startTime']+len(dataset)*dataset.attrs['cadence']
    elif 'timeIndex' in dataset.attrs:
        return dataset.attrs['startTime']+dataset[-1]['timeIndex']
    else:
        return dataset.attrs['startTime']

def matchingSignals(m,SUBCLASSES=None,prnCodes=None,codeCarriers=None,carrierCodes=None,rangingCodes=None):
    sigs = []

    if codeCarriers: codeCarriers = map(tuple,codeCarriers)
    keys = ['SUBCLASS', 'prnCode','carrierCode','rangingCode','codeCarriers']
    optvalues = [SUBCLASSES, prnCodes, carrierCodes, rangingCodes, codeCarriers]
    attrs = {k:v[0] for k,v in zip(keys,optvalues) if v is not None and len(v) == 1 and k != 'codeCarriers'}
    opts = {k:v for k,v in zip(keys,optvalues) if v is not None and k not in attrs}

    for track in m.matchingDatasets(**attrs):
        prnkey = track.attrs['prnCode']
        cc = mdh.getEnumAttributeString(track,'carrierCode')
        rc = mdh.getEnumAttributeString(track,'rangingCode')
        trkattrs = track.attrs['SUBCLASS'], prnkey, cc, rc, (cc,rc)

        if any([attr not in opts[k] for attr,k in zip(trkattrs,keys) if k in opts]):
            continue

        k = (prnkey,cc,rc)
        if k in sigs: continue
        sigs.append(k)

    return tuple(sigs)


# ============================================================
# Report generation tools
# ============================================================
class TextTable(object):
    def __init__(self, keys, fmt=None, delimiter='|', fill='-', padding=1, linebreak=None):
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
    A class which encapsulates crudely a latex table as a list of formatted strings.
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


class Report(object):
    """
    Composes a report writing base class with a formatting class for quickly
    generating custom reports.
    """
    __metaclass__ = ABCMeta

    def __init__(self, formatter=TextTable, keys=None, **kwargs):
        if keys is None:
            try:
                keys = getattr(self.__class__,'_keys')
            except AttributeError as e:
                raise AttributeError("Missing required paramter: 'keys'")

        if 'fmt' not in kwargs and hasattr(self.__class__,'_fmt'):
            kwargs['fmt'] = self.__class__._fmt

        self._formatter = formatter(keys, **kwargs)

    # python goes here if an attribute can't be found
    # thus, this effectively delegates formatting specific attributes
    # to the formatting class
    def __getattr__(self,attr):
        return getattr(self._formatter, attr)

    @abstractmethod
    def generateContent(self, *args, **kwargs):
        """ Interface method that generates report content """
        return NotImplemented

    def writeReport(self, *args, **kwargs):
        self.writeTable(*args, **kwargs)


# ============================================================
# Misc. helpers
# ============================================================
def getCorrectInput(quest,default=None,opts=('y','n'),please="Please enter 'y' or 'n': "):
    quest += " ("
    quest += "/".join(opts)
    quest += (", default=%s" % default) if default is not None else ''
    quest += "): "

    rin = raw_input(quest)
    if default is not None: opts = list(opts)+['']
    while rin not in opts: rin = raw_input(please)

    return default if default is not None and rin == '' else rin

def parseDate(datestr,timefmt=DEFAULT_TIME_FMT):
    try:                return datetime.strptime(datestr,timefmt)
    except ValueError:  pass
    try:                return du.parser.parse(datestr)
    except:             return None

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

## TODO
# this should really be in another module since it's for ephemeris and not mdh files
# replace with SGLTk's EphReader class
class EphReader(gpstk.SP3EphemerisStore):
    def __init__(self, *ephFiles):
        super(EphReader,self).__init__()
        map(self.loadFile, ephFiles)
        self.cer = gpstk.CorrectedEphemerisRange()
        self.tm = gpstk.NBTropModel()
        self.dtype = [('azimuth',np.float32),('elevation',np.float32)
                      ,('pseudorange',np.float64),('dopplerFrequencyShift',np.float64)]

    def _setRxPos(self,rxPos):
        self.tm.setReceiverHeight(rxPos.getAltitude())
        self.tm.setReceiverLatitude(rxPos.getGeodeticLatitude())

    def _getExpObs(self,rxPos,satID,time):
        self.tm.setDayOfYear(gpstk.YDSTime(time).doy)

        trop = self.tm.correction(self.cer.elevation)

        sv = self.getXvt(satID,time)
        pseudorange = self.cer.ComputeAtReceiveTime(time,rxPos,satID,self)+trop
        azimuth = self.cer.azimuth
        elevation = self.cer.elevation
        dopplerL1 = sv.v.dot(self.cer.cosines)*gpstk.L1_FREQ_GPS/gpstk.C_MPS

        return azimuth,elevation,pseudorange,dopplerL1

    def getExpObs(self,rxPos,satID,time):
        """
        Compute expected obs.
        Input: rxPos - gpstk.Position
               satID - gpstk.SatID
               time  - gpstk.CommonTime, or array/list of CommonTimes
        Output: azimuth, elevation, pseudorange, doppler (L1)
        """
        self._setRxPos(rxPos)

        if hasattr(time,'__iter__'):
            if not hasattr(time,'__len__'): time = tuple(time)
            arr = np.ndarray(len(time),dtype=self.dtype)
            if VERBOSE > 1:
                print "Computing expected obs for",
                print "[%s,%s]..." % (commontime2datetime(time[0]).strftime(DEFAULT_TIME_FMT)
                                     ,commontime2datetime(time[-1]).strftime(DEFAULT_TIME_FMT))
            try:
                for i,t in enumerate(time): arr[i] = self._getExpObs(rxPos,satID,t)
            except KeyboardInterrupt:
                pass
            return arr

        return self._getExpObs(rxPos,satID,time)


# ============================================================
# Time conversion helpers
# ============================================================
def timegps(t):
    '''
    (analog to calendar.timegm())
    Generate a timestamp from datetime object assuming GPS time
    '''
    return calendar.timegm(t.utctimetuple()) + t.microsecond/1e6 - GPS_EPOCH_IN_UNIXTIME

def gpstime(t):
    '''
    (analog to datetime module's datetime.utcfromtimestamp())
    Generate a datetime object from a GPS timestamp
    '''
    return datetime.utcfromtimestamp(t+GPS_EPOCH_IN_UNIXTIME)

def _mdhtime2gpsweeksecond(time,timeDenominator):
    seconds, fraction = divmod(time,timeDenominator)
    week, sow         = divmod(seconds,gpstk.FULLWEEK)
    ms                = (fraction*1000)//timeDenominator
    r                 = (fraction - ms*1000)/float(timeDenominator)
    return week, sow, ms, r

def mdhtime2commontime(time,timeDenominator,timeSystem=gpstk.TimeSystem.GPS):
    if not isinstance(timeSystem, gpstk.TimeSystem):
        timeSystem = gpstk.TimeSystem(timeSystem)

    if hasattr(time,'__iter__'):
        def generateCommonTimes(time,timeDenominator,timeSystem):
            if hasattr(time,'__getitem__'): time = iter(time)
            t0 = time.next()
            
            startTime = mdhtime2commontime(t0,timeDenominator,timeSystem)

            t = t0
            while t:
                ct = gpstk.CommonTime(startTime)
                seconds, fraction = divmod(t-t0, timeDenominator)
                ms = (fraction*1000)//timeDenominator
                r = (fraction - ms*1000)/float(timeDenominator)

                ct.addSeconds(seconds+r)
                ct.addMilliseconds(int(ms))
                yield ct

                t = time.next()

        return generateCommonTimes(time,timeDenominator,timeSystem)

    week, sow, ms, r = _mdhtime2gpsweeksecond(time, timeDenominator)
    t = gpstk.GPSWeekSecond(int(week), np.float64(sow+r), timeSystem).toCommonTime()
    t.addMilliseconds(int(ms)) # avoid millisecond precision loss
    #t.addMicroseconds(us) # not implemented non-abstractly in the gpstk

    return t

def commontime2datetime(cmntime):
    return datetime.utcfromtimestamp(commontime2timestamp(cmntime) + GPS_EPOCH_IN_UNIXTIME)

def datetime2commontime(d,timesys=gpstk.TimeSystem.GPS):
    return gpstk.CivilTime(d.year,d.month,d.day,d.hour,d.minute,d.second+d.microsecond*1e-6,gpstk.TimeSystem(timesys)).toCommonTime()

def commontime2timestamp(cmntime):
    gws = gpstk.GPSWeekSecond(cmntime)
    return gws.getWeek()*86400L*7L + gws.getSOW()

def commontime2mdhtime(cmntime,timeDenominator):
    gws = gpstk.GPSWeekSecond(cmntime)
    return long((gws.getWeek()*86400L*7L + gws.getSOW())*timeDenominator)

# ============================================================
# Generic csv-like log reader/writer
# ============================================================
class LogFileWriter(object):
  def __init__(self,filename,keys,typemap,delimiter,start_date=None,end_date=None):
    self.name = filename
    self.file = open(self.name,'w')    

    self.delimiter = delimiter
    self.start_date = start_date
    self.end_date = end_date
    self._keys = keys
    self._buffer = ()
    self._index = 0
    self._bufsize = 1000

    if self.start_date: print >> self.file, "start_date=%s" % self.start_date.strftime(DEFAULT_TIME_FMT)
    if self.end_date:   print >> self.file, "end_date=%s" % self.end_date.strftime(DEFAULT_TIME_FMT)
    print >> self.file, "# " + self.delimiter.join(self._keys)

    if isinstance(typemap,dict):
      self.typecast = self.delimiter.join('{%s:%s}' % (str(i)+('' if k in typemap else "!s"), typemap[k] if k in typemap else "s") for i,k in enumerate(self._keys))
    else:
      self.typecast = self.delimiter.join(['{%s:%s}' % (str(i),f) for i,f in enumerate(self.typecast)])

  def write(self,*values,**kwvalues):
    if kwvalues: values = [kwvalues[k] for k in self._keys]
    print >> self.file, self.typecast.format(*values)

  def close(self):
    self.file.close()


class LogFileReader(object):

  """
  LogFile

  File interface that parses and interprets text log files acts like a file
  (opens,closes,write,read) but has keys that are consistent across each
  iteration (i.e., it acts like a list of dicts).

  Note:
  Requires a manually defined field->type relation (i.e. the typemap argument)
  """
  def __init__(self,filename,keys,typemap,delimiter,start_date=None,end_date=None):
    self.name = filename
    self.file = open(self.name,'r')    

    self.delimiter = delimiter
    self.start_date = start_date
    self.end_date = end_date
    self._keys = keys
    self._buffer = ()
    self._index = 0
    self._bufsize = 1000

    if isinstance(typemap,dict):
      self.typecast = tuple(typemap[k] for k in self._keys)
    else:
      self.typecast = tuple(typemap)

    self.firstline = 0
    self.__ParseHeader()

  def __ParseHeader(self):
    # skip the comments
    line = self.file.readline()
    while(line[0] == '#'):
      line = self.file.readline()

    lastread = self.file.tell()
    while(1):
        try:
            #re.match(r'^[A-Za-z]+.*=.*')
            k, v = line.strip().split('=')
        except ValueError:
            self.file.seek(lastread)
            break
        if k == 'start_date' and self.start_date is None:
            self.start_date = parseDate(v,DEFAULT_TIME_FMT)
        elif k == 'end_date' and self.end_date is None:
            self.end_date = parseDate(v,DEFAULT_TIME_FMT)
        lastread = self.file.tell()                
        line = self.file.readline()
    self.firstline = lastread

  def __SkipComments(self,n=1):
    line = self.file.readline()
    while(line):
      if line[0] != '#':
        break
      line = self.file.readline()
    return line.rstrip().split(self.delimiter)

  def __contains__(self, key):
    return key in self._keys    

  def getKeyIndex(self,key):
    return self._keys.index(key)

  def keys(self):
    return list(self._keys)

  def iterkeys(self):
    return (k for k in self._keys)

  def groupby(self,key):
    return groupby(self,key=op.itemgetter(*map(self.getKeyIndex,key)))

  def seek(self,offset,from_what=0):
    self.file.seek(self.firstline+offset,from_what)

  def __iter__(self):
    return self

  def readlines(self,n=1):
    if n > 0:
      return [tuple(cast(l) for cast,l in zip(self.typecast,line.rstrip().split(self.delimiter))) for i,line in zip(xrange(n),iter(self.file)) if line[0] != '#']
    if n < 0:
      return [tuple(cast(l) for cast,l in zip(self.typecast,line.rstrip().split(self.delimiter))) for line in self.file if line[0] != '#']

  def next(self):
    line = None
    while not line:
      line = self.readlines()        
    return line

  def close(self):
    self.file.close()
