import numpy as np
import calendar
import dateutil as du
import os,sys
from datetime import datetime,timedelta
from host import mdh,signals,icdgps812
from operator import attrgetter,itemgetter
import operator as op
import StringIO
from collections import OrderedDict,namedtuple
import itertools
from scipy.interpolate import LSQUnivariateSpline
import gpstk

MAX_PSEUDORANGE =             1
C_GPS_M =                     2.99792458e8
GPS_EPOCH_IN_UNIXTIME=        315964800
GPS_EPOCH_IN_MATPLOTLIB_DATE_FORMAT=722820.0
TIME_FMT_STR =                "%Y-%m-%d %H:%M:%S.%f"

VERBOSE =                     0

RangingCode2h5t = {'CA':icdgps812.rcCA,'CL':icdgps812.rcCL,'CM':icdgps812.rcCM
                   ,'I5':icdgps812.rcI5,'Q5':icdgps812.rcQ5}
h5t2RangingCode = {v:k for k,v in RangingCode2h5t.iteritems()}
CarrierCode2h5t = {'L1':1,'L2':2,'L5':5}
h5t2CarrierCode = {v:k for k,v in CarrierCode2h5t.iteritems()}

isvalid = lambda field,dic: bool(field in dic and dic[field] is not None)

LogDataType = {"event_id":int
              ,"event_date":du.parser.parse
              ,"sv_identifier":int
              ,"prn":int
              ,"carrier_code":str
              ,"ranging_code":str
              ,"detection_period":float
              ,"max_detection_value":float
              ,"threshold":float
              ,"file_path":str
              ,"start_date":du.parser.parse
              ,"end_date":du.parser.parse}

TrackKeys = ("sv_identifier","prn","carrier_code","ranging_code","start_date","end_date")
EventKeys = ("event_id","event_date","sv_identifier","prn","carrier_code","ranging_code","detection_period","max_detection_value","threshold","file_path")

EVENTSLOG_DEFAULT_FMT = ','.join(['%%(%s)s'%k for k in EventKeys])
TRACKSLOG_DEFAULT_FMT = ','.join(['%%(%s)s'%k for k in TrackKeys])

opstr2opfun = {'>=':op.ge,'<=':op.le,'<':op.lt, '>':op.gt, '==':op.eq
               , 'not in': lambda x,y: op.not_(op.contains(y,x))
               ,'in': lambda x,y: op.contains(y,x)
               ,'!=':op.ne}


def getStopTime(dataset):
    if "cadence" in dataset.attrs:
        return dataset.attrs['startTime']+len(dataset)*dataset.attrs['cadence']
    elif 'timeIndex' in dataset.attrs:
        return dataset.attrs['startTime']+dataset[-1]['timeIndex']
    else:
        return dataset.attrs['startTime']


def mdhSpan(m,absolute=False,**attrs):
    dsstart, dsstop = zip(*[(ds.attrs['startTime']/(ds.attrs['timeDenominator'] if absolute else 1)
                             ,getStopTime(ds)/(ds.attrs['timeDenominator'] if absolute else 1))
                             for ds in m.matchingDatasets(**attrs)])
    return min(dsstart), max(dsstop)
      

class TextTable(object):
    def __init__(self,keys,widths=None,fmt=None, delimiter='|', fill='-', padding=1, out=sys.stdout):
        self._keys = list(keys)
        self.out = out

        if isinstance(widths,dict):
            self.widths = [widths[k] if k in widths else 0 for k in self._keys]
        elif widths is not None:
            self.widths = list(widths)
        else:
            self.widths = []

        if len(self.widths) < len(self._keys):
            self.widths[len(self.widths):] = [0]*(len(self._keys)-len(self.widths))
        self.widths = OrderedDict(zip(self._keys,[max(w,len(k)) for k,w in zip(keys,self.widths)]))

        if isinstance(fmt,dict):
            self._fmt = [fmt[k] if k in fmt else "{!s:}" for k in self._keys]
        elif fmt is not None:
            self._fmt = list(fmt)
        else:
            self._fmt = []

        if len(self._fmt) < len(self._keys):
            self._fmt[len(self._fmt):] = ["{!s:}"]*(len(self._keys)-len(self._fmt))
        self._fmt = OrderedDict(zip(self._keys,self._fmt))

        self.fill = fill
        self.padding = padding
        self.delimiter = ' '*padding + delimiter + ' '*padding

        self.header = self.delimiter.join([str(k).center(w) for k,w in self.widths.items()])
        self.separator = (self.fill*padding + '+' + self.fill*padding)
        self.linebreak = (self.fill*padding + self.fill + self.fill*padding).join([self.fill*(w) for w in self.widths.values()])

    def printHeader(self,linebreak=True):
        print >> self.out, self.header
        if linebreak:
            print >> self.out, self.separator.join([self.fill*(w) for w in self.widths.values()])

    def printLineBreak(self):
        print >> self.out, self.linebreak

    def printLine(self,*vals,**kwvals):
        print >> self.out, self.delimiter.join([(self._fmt[k].format(vals[i] if i < len(vals) else kwvals[k]) if (i < len(vals) or k in kwvals) else "-").rjust(self.widths[k])
                                                for i,k in enumerate(self._keys)])


def parseDate(datestr,timefmt):
    try:
        return datetime.strptime(datestr,timefmt)
    except ValueError:
        pass
    try:
        return du.parser.parse(datestr)
    except:
        return None


def mask2runs(mask):
  '''
  mask2runs

  Turns a boolean array into an array of indices indicating the start and end indices of each run of 1s
  '''
  runflip = np.nonzero(mask[1:] ^ mask[:-1])[0]+1 # locations where a run of 1's starts

  # if the last run was a run of 1's, count the last epoch as the end of the run
  if mask[-1]: runflip = np.r_[runflip,len(mask)]
  if mask[0]: runflip = np.r_[0,runflip]

  return runflip[::2],runflip[1::2]


def printHeader(msg, fill='-'):
  sys.stdout.write(msg + '\n' + (fill*len(msg)) + '\n')
  sys.stdout.flush()


def CodeStrGen(prn=None,sv_identifier=None,carrier_code=None,ranging_code=None):
  if sv_identifier is not None: svkey =   ('SVN%02u' % int(sv_identifier))
  elif prn is not None: svkey = ('PRN%02u' % int(prn))
  else: svkey = ''

  codekey = ':'.join([('%s' % str(carrier_code)) if carrier_code is not None else '',('%s' % str(ranging_code)) if ranging_code is not None else ''])
  
  return ' '.join([svkey,codekey])


def CodeIDGen(prn=None,sv_identifier=None,carrier_code=None,ranging_code=None):
  if sv_identifier is not None: svkey =   ('SVN%02u' % int(sv_identifier))
  elif prn is not None: svkey = ('PRN%02u' % int(prn))
  else: svkey = ''

  codekey = ''.join([('%s' % str(carrier_code)) if carrier_code is not None else '',('%s' % str(ranging_code)) if ranging_code is not None else ''])
  
  return ''.join([svkey,codekey])


class LogFile(object):

  """
  LogFile

  File interface that parses and interprets text log files acts like a file
  (opens,closes,write,read) but has keys that are consistent across each
  iteration (i.e., it acts like a list of dicts).

  Note:
  Requires a manually defined field->type relation dictionary (LogDataType)
  """

  def __init__(self,filepath,mode='r',logfmt=None,delimiter=",",selectlines=False,start_date=None,end_date=None
               ,form={}):
    self.filepath = filepath
    self.delimiter = delimiter
    self.logfmt = logfmt
    self.selectlines = selectlines
    self.mode = mode
    self.comment= StringIO.StringIO() if mode == 'w' else ''
    self.start_date = start_date
    self.end_date = end_date
    self.firstline=None
    self.form = {k:v[1] for k,v in form.items()}
    self.cmpdict = {k:opstr2opfun[v[0]] for k,v in form.items()}

    if mode not in ('w','wb') and mode != 'r':
      raise Exception("LogFile Error: Only 'r','w' and 'wb' modes are currently supported")

    self.file = open(filepath,mode)
    if mode == 'r':
      # attempt to read in log format if one was not manually specified
      if self.logfmt is None and not self.__ParseHeader():
        self.file.close()          
        raise Exception("LogFile Error: Could not find key=value pair defining log format"
                        + " and read format not manually specified.")

      if self.start_date is None or self.end_date is None:
        print "LogFile Warning: Start or end date was not found or could not be parsed"

    elif mode in ('w','wb'):
      if not self.logfmt:
        self.file.close()
        raise Exception("LogFile Error: No log output formatting string specified.")

    self.__definekeys(*(self.__defineformatfun()))
    self.cast = tuple(LogDataType[k] for k in self._keys)
    self._thisline = OrderedDict(zip(self._keys,[None]*len(self._keys)))

  def __ParseHeader(self):
        # skip the comments
        line = self.file.readline()
        while(line[0] == '#'):
          line = self.file.readline()

        lastread = self.file.tell()
        while(1):
            try:
                k, v = line.strip().split('=')
            except ValueError:
                self.file.seek(lastread)
                break
            if k == 'start_date' and self.start_date is None:
                self.start_date = parseDate(v,TIME_FMT_STR)
            elif k == 'end_date' and self.end_date is None:
                self.end_date = parseDate(v,TIME_FMT_STR)
            elif k == 'logfmt' and self.logfmt is None:
                self.logfmt = v
            lastread = self.file.tell()                
            line = self.file.readline()
        self.firstline = lastread

        return self.logfmt is not None
    

  """
  definekeys
  
  Parses 'key' fields implied from the output log format string.  Requires
  parenthesis type used for formatting.
  """
  def __definekeys(self,ob,cb):
    self._keys = tuple(s[s.find(ob)+1:s.find(cb)] for s in self.logfmt.split(self.delimiter))


  """
  defineformatfun
  
  Infers an output format based on log formatting string.  Returns parenthesis
  type used for formatting.
  """  
  def __defineformatfun(self):
    # define a function for outputting formatted data
    if "%(" in self.logfmt and ")" in self.logfmt:
      self.strfline = lambda fmtstr,dic: fmtstr % dic            
      return '(',')'
    elif "{" in self.logfmt and "}" in self.logfmt:
      self.strfline = lambda fmtstr,dic: fmtstr.format(**dic)       
      return '{','}'
    else:
      self.file.close()      
      raise Exception("LogFile Error: string formatting only supports parenthesis or curly brace "
                       + "type format strings. e.g., %(event_id)%d or {event_id:d}.")


  def __contains__(self, key):
    return key in self._keys    


  def getKeyIndex(self,key):
    return self._keys.index(key)

  """
  addComment

  Add a comment string to top of log file
  """
  def addComment(self,comment):
    print >> self.comment, '\n'.join(["#%s" % line for line in comment.split('\n')])

  """
  writeDictList

  write an object that acts like a list of dictionaries out to file
  """
  def writeDictList(self,dictlist):
    self.file.write( '\n'.join([self.strfline(self.logfmt,row) for row in dictlist]) )


  def keys(self):
    return list(self._keys)

  def iterkeys(self):
    return (k for k in self._keys)

  def itervalues(self):
    self.seek(0)    
    while (1):
      try:
        yield self.next().values()
      except StopIteration:
        raise

  def values(self):
    return list(self.itervalues())

  def seek(self,offset,from_what=0):
    if self.mode == 'r': # seek from the first valid event line
      self.file.seek(self.firstline+offset,from_what)
    else:
      self.file.seek(offset,from_what)

  def __iter__(self):
    return self

  def __readNextValidLine(self):
    line = self.file.readline()
    while line and ((self.selectlines and (line[0] != '*')) or (line[0] == '#')):
#### comment above and use below if desired to plot questionable events (indicated by '?' in event log)
    # while line and ((self.selectlines and (line[0] != '*' or line[0] != '?')) or (line[0] == '#')):
####
      line = self.file.readline()
    if not line: raise StopIteration

    while(line[0] in ('*','?')): line = line[1:]
    self._thisline.update(zip( self._keys, [cast(v) for cast,v in zip(self.cast, line.strip().split(self.delimiter))] ))
    
  def next(self):
    self.__readNextValidLine()

    while self.form and not all(fcmp(self._thisline[k],self.form[k]) for k,fcmp in self.cmpdict.items()):
      self.__readNextValidLine()

    return OrderedDict(self._thisline)
    
  def close(self):
    self.file.close()

    # merge the comments with data file
    if self.mode == 'w':
      print >> self.comment, "start_date=%s" % self.start_date.strftime(TIME_FMT_STR)
      print >> self.comment, "end_date=%s" % self.end_date.strftime(TIME_FMT_STR)
      print >> self.comment, "logfmt=%s" % self.logfmt # output log format to file
      pushTextToFile(self.filepath,self.comment.getvalue()) # output comments to log
      self.comment.close()


"""
pushTextToFile

insert text into the top of a file
"""
def pushTextToFile(filepath, newtext):
  with open(filepath,'r+') as fh:
    tmp = fh.read()
    fh.seek(0)
    fh.write(newtext)
    fh.write(tmp)


def getCorrectInput(quest,default=None,opts=('y','n'),please="Please enter 'y' or 'n': "):
    quest += " ("
    quest += "/".join(opts)
    quest += (", default=%s" % default) if default is not None else ''
    quest += "): "

    rin = raw_input(quest)
    if default is not None: opts = list(opts)+['']
    while rin not in opts: rin = raw_input(please)

    return default if default is not None and rin == '' else rin

  
def mergelogs(lognames, **kwargs):
    times = []
    commonkeys = []
    logfiles = []
    timetuple = attrgetter('start_date','end_date')

    for name in lognames:
      logfiles += [LogFile(name, **kwargs)]
      commonkeys += [logfiles[-1].keys()]
      times += [timetuple(logfiles[-1])]

    allkeys = set(itertools.chain(*commonkeys))
    confkeys = allkeys.difference(set.intersection(*map(set,commonkeys)))
    if confkeys and VERBOSE:
        print "LogFile Error: Log files have conflicting fields: " + ', '.join(map(repr,confkeys))          
        if getCorrectInput("Drop conflicting fields? (y/n, default=y): ") not in ('y',''):
          raise Exception("LogFile: Could not merge. Log files have differing data.")

    data = []
    oldlen = 0
    commonkeys = [k for k in commonkeys[0] if k in set.intersection(*map(set,commonkeys))]
    for log in logfiles:
      if VERBOSE: print "Reading log file '%s'..." % os.path.basename(log.file.name) ;sys.stdout.flush()
        
      data += [[row[log.getKeyIndex(k)] for k in commonkeys] for row in log.itervalues()]

      if VERBOSE > 1: print "\t%d lines read." % (len(data)-oldlen)
      if VERBOSE > 2: print "\tContains data fields:", tuple(log.keys())
      oldlen = len(data)
      log.close(); del log

    return commonkeys, data, min(zip(*times)[0]), max(zip(*times)[1])
      

def mergeIntervals(intervals,startkey,stopkey,presorted=False,sortkey=None,eps=None):
    """
    Merge overlapping time intervals

    intervals -- list of elements containing a start and end time in each element

    presorted -- Defines whether intervals has already been sorted by start
    and end times

    startkey,stopkey -- accessor (index or key mapping) for items in each
    element which store start and stop times

    out -- list of disjoint time intervals represented by tuples of
    start and end time objects where the type is the same as the time
    type stored in the original list
    """
    key=sortkey if sortkey is not None else itemgetter(startkey,stopkey)
    itvs = map(key,intervals)
    if not presorted: itvs.sort()

    startitvs,stopitvs = np.array(zip(*itvs))
    if eps is None: eps = type(stopitvs[0]-startitvs[0])()

    isnewitv = (startitvs[1:] - stopitvs[:-1]) > eps
    enditv = list(stopitvs[:-1][isnewitv])
    newitv = list(startitvs[1:][isnewitv])

    if not enditv:                          enditv.append(stopitvs[-1])
    elif (startitvs[-1]-enditv[-1]) > eps:  enditv.append(stopitvs[-1])
    else:                                   enditv[-1] = max(stopitvs[-1], enditv[-1])

    if not newitv:                          newitv.insert(0,startitvs[0])
    elif (newitv[0]-stopitvs[0]) > eps:     newitv.insert(0,startitvs[0])
    else:                                   newitv[0] = min(startitvs[0], newitv[0])

    return zip(newitv,enditv)
    

def FindBestTrack(m,centerTime,prn=None,carrier_code=None,ranging_code=None,SUBCLASS='IQ_SUMS',associatedDataset=None):
  """
  FindBestTrack

  Return track in MDH file that best matches the specified
  time in the following sense:
  
  maximizes min(centerTime-trackStartTime, centerTime+trackStartTime)

  See http://sglwiki.arlut.utexas.edu/bin/view/SGL/SwRxMDHDataDictionary
  for valid file attributes.

  m -- MDH file handler
  centerTime -- requested time
  prn -- MDH prn attribute value to exclude from search
  ranging_code -- MDH ranging code attribute value
  SUBCLASS -- MDH subclass attribute value. Usually either IQ_SUMS or OBSERVATIONS.
  out -- dataset that best matches input time or None if none are found
  """
  max_coverage = 0
  bestTrack = None

  attrs = {}
  if prn is not None: attrs['prnCode'] = prn
  if ranging_code: attrs['rangingCode'] = ranging_code
  if carrier_code: attrs['carrierCode'] = carrier_code
  attrs['SUBCLASS']=SUBCLASS

  if associatedDataset is None:
    tracklist = m.matchingDatasets(**attrs)
  else:
    tracklist = m.associatedDatasets(associatedDataset,**attrs)

  for track in tracklist:
    if 'cadence' not in track.attrs: continue
    trackStartTime = tracktimefromslice(track,0)
    trackStopTime = tracktimefromslice(track,len(track))

    if SUBCLASS == 'OBSERVATIONS': thisTime = centerTime + track['pseudorange'][0]/C_GPS_M
    else: thisTime = centerTime

    if (thisTime >= trackStopTime or thisTime <= trackStartTime): continue

    coverage = min(thisTime-trackStartTime,trackStopTime-thisTime)
    if coverage > max_coverage:
      bestTrack = track
      max_coverage = coverage

  return bestTrack

def FindBestRefTracks(m, requestTime, detectionPeriod, prn, codeCarriers=None):
  """
  FindBestTrack

  Return track in MDH file corresponding to a different satellite than
  that specified that best matches the specified time in the following
  sense:
  
  maximizes min(requestTime-trackStartTime, requestTime+trackStartTime)
  maximizes min(snr)

  See http://sglwiki.arlut.utexas.edu/bin/view/SGL/SwRxMDHDataDictionary
  for valid file attributes.

  m -- MDH file handler
  requestTime -- requested time
  prn -- MDH prn attribute value to exclude from search
  detectionPeriod -- time in seconds of the detected events detection period
  out -- dataset that best matches input time or None if none are found
  """
  maxminsnr = 0
  obsTrack = None
  for track in m.matchingDatasets(SUBCLASS='OBSERVATIONS'):
    cc = mdh.getEnumAttributeString(track,'carrierCode')+':'+mdh.getEnumAttributeString(track,'rangingCode')
    if codeCarriers and cc not in codeCarriers: continue
    if track.attrs['prnCode'] == prn: continue

    trackStartTime = tracktimefromslice(track,0) - track['pseudorange'][0]/C_GPS_M
    trackStopTime = tracktimefromslice(track,len(track)) - track['pseudorange'][-1]/C_GPS_M

    if (requestTime >= trackStopTime or requestTime <= trackStartTime): continue

    minsnr = min(track['snr'])
    if min(requestTime-trackStartTime,trackStopTime-requestTime) > (1.5*detectionPeriod) and minsnr > maxminsnr:
      obsTrack = track
      maxminsnr = minsnr
  if maxminsnr == 0:
    return (None,None)

  iqTrack = FindBestTrack(m,requestTime,SUBCLASS='IQ_SUMS',associatedDataset=obsTrack)
      
  return iqTrack,obsTrack


def matchingSignals(m,startTime=None,stopTime=None,SUBCLASSES=None,prnCodes=None,codeCarriers=None,carrierCodes=None,rangingCodes=None):
    if startTime: startTime = timegps(startTime)
    if stopTime: stopTime = timegps(stopTime)

    sigs = []

    if any(isinstance(x,str) for x in (SUBCLASSES,prnCodes,codeCarriers,carrierCodes,rangingCodes)):
        if VERBOSE: print "Error: Expected a list of values; got a single value"
        return ()

    if codeCarriers: codeCarriers = map(tuple,codeCarriers)
    keys = ['SUBCLASS', 'prnCode','carrierCode','rangingCode','codeCarriers']
    optvalues = [SUBCLASSES, prnCodes, carrierCodes, rangingCodes, codeCarriers]
    attrs = {k:v[0] for k,v in zip(keys,optvalues) if v != None and len(v) == 1 and k != 'codeCarriers'}
    opts = {k:v for k,v in zip(keys,optvalues) if v != None and k not in attrs}

    for track in m.matchingDatasets(**attrs):
        if (startTime and (startTime*track.attrs['timeDenominator'] >= getStopTime(track))
            or stopTime and (stopTime*track.attrs['timeDenominator'] <= track.attrs['startTime'])):
            continue

        prnkey = track.attrs['prnCode']
        cc = mdh.getEnumAttributeString(track,'carrierCode')
        rc = mdh.getEnumAttributeString(track,'rangingCode')
        trkattrs = track.attrs['SUBCLASS'], prnkey, cc, rc, (cc,rc)

        if any([val not in opts[k] for val,k in zip(trkattrs,keys) if k in opts]):
            continue

        k = (prnkey,cc,rc)
        if k in sigs: continue
        sigs.append(k)

    return tuple(sigs)


"""
Get Times

Generates the array of timestamps where OBS and IQ are overlapping
(OBS may start slightly before and end after IQ) and the associated
track indices
"""
def GetTimes(iqTrack,obsTrack,tslice=None):
  iqSliceStart,iqSliceStop = tslice if tslice else (0,len(iqTrack))

  times = {}

  obsStartTime = tracktimefromslice(obsTrack,0) - (obsTrack['pseudorange'][0]/C_GPS_M)
  obsStopTime = tracktimefromslice(obsTrack,len(obsTrack)) - (obsTrack['pseudorange'][-1]/C_GPS_M)
  iqStartTime = tracktimefromslice(iqTrack,iqSliceStart)
  iqStopTime = tracktimefromslice(iqTrack,iqSliceStop)
  
  # Check if there is anything overlapping, crop to the overlap
  if (iqStartTime >= obsStopTime or iqStopTime <= obsStartTime):
    return {}

  # have to be careful matching times - the precision of timestamps are down to the microsecond
  # but stored track times are typically at millisecond precision
  # e.g., mapping an obs time with 20ms precision to an iq time with 1ms precision
  # doesn't necessarily give you an integer slice into iq time due to the variable
  # microsecond portion of the timestamp
  if obsStartTime > iqStartTime:
    iqSliceStart = long(slicefromtracktime(iqTrack,obsStartTime)+1) ## round up to the next slice
  if obsStopTime < iqStopTime:
    iqSliceStop = long(slicefromtracktime(iqTrack,obsStopTime)) ## round down to the previous slice
  times['iqSlice'] = slice(iqSliceStart,iqSliceStop)

  # Find the slice in the obstrack (in receiver time)
  # that corresponds to the iqSliceStart (where IQ is in sattelite time)
  obsSliceStart = long(FindIndexInOtherTrack(iqTrack,iqSliceStart,obsTrack))
  if obsSliceStart < 0 : obsSliceStart = 0
  obsTimeStop  = tracktimefromslice(iqTrack,iqSliceStop)+MAX_PSEUDORANGE # Assume pseudorange delay < MAX
  obsSliceStop = long(slicefromtracktime(obsTrack,obsTimeStop)+1)

  # Truncate iqSlice to the last sample of obs track if it is longer than obs
  if obsSliceStop > len(obsTrack):
    obsSliceStop = len(obsTrack)
    iqSliceStop -= long((trackrate(iqTrack)/trackrate(obsTrack))+1)
    times['iqSlice'] = slice(iqSliceStart,iqSliceStop)

  # Check for pseudorange exceeding the expected case
  if((tracktimefromslice(obsTrack,obsSliceStart)-(obsTrack['pseudorange'][obsSliceStart]/C_GPS_M))
     > tracktimefromslice(iqTrack,iqSliceStart)
     or
     (tracktimefromslice(obsTrack,obsSliceStop)
      -(obsTrack['pseudorange'][min(obsSliceStop,len(obsTrack)-1)]/C_GPS_M)) # min is needed to avoid out of bounds
     < tracktimefromslice(iqTrack,iqSliceStop)):
    if VERBOSE:
      print "Pseudorange offset exceeded maximum expected in track:\n%s. Ignoring OBSERVATIONS track..." % (
            obsTrack.name)
    return {}

  times['obsSlice'] = slice(obsSliceStart,obsSliceStop)

  times['iqTimes'] = ((iqTrack.attrs['startTime'] + np.arange(iqSliceStart,iqSliceStop)*iqTrack.attrs['cadence'])
                      /float(iqTrack.attrs['timeDenominator']))

  # Generate observation track timestamps subtracting out propagation delay to
  # get times in apparent sattelite time
  times['obsTimes'] = (((obsTrack.attrs['startTime']
                           + np.arange(obsSliceStart,obsSliceStop)*obsTrack.attrs['cadence'])
                        /float(obsTrack.attrs['timeDenominator']))                      
                        - obsTrack['pseudorange'][obsSliceStart:obsSliceStop]/C_GPS_M)

  return times

"""
GetDetrendedPhase

Get the detrended phase of a track over some time interval.
"""
def GetDetrendedPhase(iqTrack,obsTrack,tslice=None,cadence=None,obsoffset=None):
  iqSliceStart,iqSliceStop = tslice if tslice is not None else (0,len(iqTrack))

  result = {}

  times = GetTimes(iqTrack,obsTrack,(iqSliceStart,iqSliceStop))
  if times is None:
    return {}

  iqTimes = times['iqTimes']
  iqSlice = times['iqSlice']
  obsTimes = times['obsTimes']
  obsSlice = times['obsSlice']  

  obsADR = np.interp(iqTimes,obsTimes,obsTrack[obsSlice]['accumulatedDeltaRange'])
  obsCarrierPhaseAdj = obsADR*360.0

  iq = iqTrack[iqSlice]['i'] + 1j*iqTrack[iqSlice]['q']

  try:
    promptlagIdx = np.argmin(abs(iqTrack.attrs['lagOffsetsInSecs']))
  except KeyError:
    promptlagIdx = len(iqTrack.attrs['lagOffsetsInSecs'])/2 # assume N tap correlator and use center tap

  # Build up some SNR by integrating squared data, then unwrap phase
  squaredIQ = iq[:,promptlagIdx]**2
  squaredIQSmoothed = np.convolve(squaredIQ,np.ones(20)/20,'same') # Integrate out to 20 samples (milliseconds)
  squaredIQSmoothedAngle = np.unwrap(np.angle(squaredIQSmoothed))*0.5

  # Unwrap full rate data
  iqPhase = np.angle(squaredIQ)*0.5 # Squaring phase detector
  
  # Combine full rate data with cycle count from squaredIQSmoothedAngle
  #   such that we count full cycles with the higher SNR squaredIQSmoothedAngle,
  #   but get temporal resolution from unsmoothed iqPhase
  integerHalfCycleAdj = np.around((squaredIQSmoothedAngle-iqPhase)*2)*0.5
  iqPhase += integerHalfCycleAdj
  
  rawPhase=np.degrees(iqPhase)-obsCarrierPhaseAdj

  if cadence is None:
    result['detrendedPhase'] = rawPhase-np.polyval(np.polyfit(iqTimes-iqTimes[0],rawPhase,2),iqTimes-iqTimes[0])
  else:
    winsize = int(cadence*iqTrack.attrs['timeDenominator']/iqTrack.attrs['cadence'])
    rows = np.arange(len(rawPhase))
    knotRows = np.arange(winsize/2,len(rows)-winsize/2,winsize)

    result['detrendedPhase'] = rawPhase - LSQUnivariateSpline(rows,rawPhase,knotRows,k=3)(rows)
  
  result.update(times)

  return result


"""
computeSmoothedSNR

Compute moving average over interval_secs
"""
def computeSmoothedSNR(obsTrack,span=None,interval_secs=1.0):
  startSlice,stopSlice = span if span else (0,len(obsTrack))

  interval=int(np.around(interval_secs*trackrate(obsTrack)))
  snrCumSum = np.cumsum(obsTrack['snr'][startSlice:stopSlice])

  return (snrCumSum[interval:]-snrCumSum[:-interval])*1.0/interval


"""
GetWindow

Get a time window of data from an mdh file (meant to make plotting easier in ipython)
Note that the return type is a list since a requested time interval can span multiple
tracks.
return: a list of 3 tuples containing the track, start index in that track,
        and stop index in the track.
"""
def GetWindow(mdhfile,centerDate,prn=None,ranging_code=None,window=10,SUBCLASS='IQ_SUMS'):
  centerTime = timegps(centerDate)
  startTime = centerTime - window/2.
  stopTime = centerTime + window/2.
  attrs = {}
  attrs['SUBCLASS'] = SUBCLASS
  if prn is not None: attrs['prnCode'] = prn
  if ranging_code is not None: attrs['rangingCode'] = RangingCode2h5t[ranging_code]

  for track in mdhfile.matchingDatasets(**attrs):
    if (tracktimefromslice(track,0) > stopTime or
        tracktimefromslice(track,len(track)) < startTime):
      continue
    start = max(startTime*track.attrs['timeDenominator'],track.attrs['startTime'])
    stop = min(stopTime*track.attrs['timeDenominator'],getStopTime(track))
    if stop <= start: continue

    n = (stop-start)/track.attrs['cadence']
    sliceStart = (start-track.attrs['startTime'])/track.attrs['cadence']

    yield (track,sliceStart,sliceStart+n)


def getCarrierRate(iqTrack):
  return signals.findCarrierInfoBy812Code(int(iqTrack.attrs['carrierCode'])).nominalFreqHz


"""
phasefromdetectionmetric

Gives the phase jump in degree corresponding to a detection metric assuming a step jump
"""
def phasefromdetectionmetric(detectionMetric,tau,carrierRate):
  return detectionMetric*getXThreshold(tau)*carrierRate*360


"""
getYThreshold

Compute detection threshold for a given 'tau' (seconds) in terms of a two-sample (ie Allan variance)
See discussion in "A GPS Carrier Phase Stability Metric Draft Version 0.97"
by York et al dated 11 Dec 2010, particularly the limit line in Figure 12
"""
def getYThreshold(tau):
  return 1e-11*(tau**-1.5)+4e-11


"""
getXThreshold

Compute detection threshold for a given 'tau' (seconds) appropriate for
comparing to values of a second difference
"""
def getXThreshold(tau):
  return getYThreshold(tau)*tau*np.sqrt(2)


"""
computeSecondDifference

Compute second difference of array 'te' using an integration period of 'n' samples
(Computed via a third-difference approach, see C.A. Greenhall, "A shortcut for
estimating the modified Allan variance", Proc. 1992 IEEE Frequency Control Symp)
"""
def computeSecondDifference(te,n):
  w=np.zeros(len(te)+1,dtype=te.dtype)
  np.cumsum(te,out=w[1:])
  return (w[3*n:]-3*w[2*n:-1*n]+3*w[n:-2*n]-w[:-3*n])*(1./n)


"""
computeDetectionMetric

Compute the detection metric for a given 'observedPhaseInDegrees' array with a
'sampleRate' in Hertz using an integration period of 'tau' seconds.
"""
def computeDetectionMetric(observedPhaseInDegrees, tau, carrierFreqHz, sampleRate):
  timeError= observedPhaseInDegrees/360./carrierFreqHz
  two_sample_data=computeSecondDifference(timeError,int(tau*sampleRate))**2
  detectionMetric=np.sqrt(two_sample_data)/getXThreshold(tau)
  return detectionMetric


"""
FindIndexInOtherTrack

From an index from trackA, directly finds corresponding index in trackB
Note: returns a float expressing the exact fractional index
      so caller can decide whether to round up or down.
      Return value can also be negative or greater than the length of the track!

The equivalent expression is:
(indexA*cadenceA + tA,start - tB,start) / cadenceB
where cadenceX := X.attrs['cadence']/X.attrs['timeDenominator'] is the sample period
tA,start := X.attrs['startTime']/X.attrs['timeDenominator']
"""
def FindIndexInOtherTrack(trackA,indexA,trackB):
  return ( (indexA*trackA.attrs['cadence']+trackA.attrs['startTime'])
           *trackB.attrs['timeDenominator']
           -trackB.attrs['startTime']*trackA.attrs['timeDenominator']
           )/float(trackA.attrs['timeDenominator']*trackB.attrs['cadence'])


"""
getTimeDeltaBetweenTracks

"""
def getTimeDeltaBetweenTracks(trackA,trackB,indexA=0,indexB=0):
    return  ((trackA.attrs['startTime']+trackA.attrs['cadence']*indexA)
             *trackB.attrs['timeDenominator']          
             -
             (trackB.attrs['startTime']+trackB.attrs['cadence']*indexB)
             *trackA.attrs['timeDenominator']
             )/float(trackA.attrs['timeDenominator']*trackB.attrs['timeDenominator'])


def numfromslice(track,tslice):
  return ((trkSlice*track.attrs['cadence']+track.attrs['startTime'])
          /(track.attrs['timeDenominator']*86400.)) + GPS_EPOCH_IN_MATPLOTLIB_DATE_FORMAT


def tracktimefromslice(track,trkSlice):
  return ((trkSlice*track.attrs['cadence']+track.attrs['startTime'])
          /float(track.attrs['timeDenominator']))


def slicefromtracktime(track,time):
  return ((time*track.attrs['timeDenominator']-track.attrs['startTime'])
          /float(track.attrs['cadence']))


"""
strfgpstimestamp

Generate formatted date string from GPS timestamp
"""
def strfgpstimestamp(t,fmt=TIME_FMT_STR):
  return (gpstime(t)).strftime(fmt)


"""
strpgpstimestamp

Generate timestamp from formatted date string assuming GPS time
"""
def strpgpstimestamp(t,fmt=TIME_FMT_STR):
  return timegps(datetime.strptime(t,fmt))


"""
timegps (analogue to calendar.timegm())

Generate a timestamp from datetime object assuming GPS time
"""
def timegps(t):
  return calendar.timegm(t.utctimetuple()) + t.microsecond/1e6 - GPS_EPOCH_IN_UNIXTIME


"""
gpstime (analog function to datetime module's datetime.utcfromtimestamp())

Generate a datetime object from a GPS timestamp
"""
def gpstime(t):
  return datetime.utcfromtimestamp(t+GPS_EPOCH_IN_UNIXTIME)


"""
trackrate

Helper function that returns track sample rate in Hz (samples per second)
"""
def trackrate(track):
  return track.attrs['timeDenominator'] / float(track.attrs['cadence'])


def mdhtime2commontime(time,timeDenominator,timeSystem=gpstk.TimeSystem.GPS):
    seconds,fraction = divmod(time,timeDenominator)
    week,sow = divmod(seconds,gpstk.FULLWEEK)
    return gpstk.GPSWeekSecond(int(week),sow+np.float64(fraction)/timeDenominator,
                               gpstk.TimeSystem(timeSystem)).toCommonTime()

def commontime2mdhtime(cmntime,timeDenominator):
    gws = gpstk.GPSWeekSecond(cmntime)
    return long((gws.getWeek()*86400L*7L + gws.getSOW())*timeDenominator)

  
def commontime2datetime(cmntime):
    return datetime.utcfromtimestamp(commontime2timestamp(cmntime) + GPS_EPOCH_IN_UNIXTIME)


def datetime2commontime(d,timesys=gpstk.TimeSystem.GPS):
    return gpstk.CivilTime(d.year,d.month,d.day,d.hour,d.minute,d.second+d.microsecond*1e-6,gpstk.TimeSystem(timesys)).toCommonTime()


def commontime2timestamp(cmntime):
    gws = gpstk.GPSWeekSecond(cmntime)
    return gws.getWeek()*86400L*7L + gws.getSOW()
