#!/usr/bin/python

'''
data availability analysis for python!

Thomas J Garcia

Note:
+ GPS system is assumed!
+ demodulatorStatus code 255 is ignored here (although it appears in the data) because it is undocumented as far as I know
'''

import sys, os
import numpy as np
import mdh
import gpstk
import calendar, time
import psycopg2 as ps
import psycopg2.extras
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import dates, ticker
from datetime import datetime, timedelta
import operator as op
from collections import OrderedDict
import mdh_tools as mtools
from mdh_tools import timegps,gpstime,commontime2datetime,datetime2commontime
from mdh_tools import LatexTable, TextTable, Report
import h5py
import glob



VERBOSE = 0
DEBUG = 0

UNAVAILABLE = 2**6
lockCodes = OrderedDict(zip(("Code Carrier Lock", "Unlocked", "Code Lock", "No Code", "Invalid"),(2,0,1,3,4)))

statusColors = ('#2c7bb6','#1a9641','#d7191c','#fdae61','#984ea3','k')

statusCodes = ["Ephemeris"] + lockCodes.keys()
statusValues = [5,] + lockCodes.values()
statusCodes = OrderedDict(zip(statusCodes,map(lambda x: 2**x, statusValues)))


def genericPlot(plotter,args,title,start_date,stop_date,savedir=None,filetag=None,dpi=100,**kwargs):
  if VERBOSE: sys.stdout.write("Generating plot \"%s\"..." % title); sys.stdout.flush()

  axtitle = "%s to %s GPS" % (start_date,stop_date)
  callargs = [start_date,stop_date]+list(args)
  fig = plotter(*callargs,title=title + ' ' + axtitle,**kwargs)

  if not fig:
    if VERBOSE: sys.stdout.write("Failed.\n")
    return

  s = start_date.strftime("%Y%m%d-%H%M%S")
  e = stop_date.strftime("%Y%m%d-%H%M%S")

  if not filetag:
      filetag = '_'.join( title.lower().split(' ')) + ('_%s_%s' % (s,e))

  if savedir is not None:
      fig.savefig(os.path.join(savedir,filetag),dpi=dpi)
  else:
      plt.show(fig.number)
  plt.close(fig)

  if VERBOSE: print


class AvailabilityReport(Report):
    _keys    = ["PRN", "CC", "Type", "Elv", "PassStart", "PassStop", "Epochs"
                 , "Unlocked", "Missed"]
    keys = _keys # for public access

    def generateContent(self, svdict, summaryOnly=True):
        numpct = lambda num,denom: "%d (%3.2f%%)" % (num, 100.*num/denom)
        toffset = lambda t,ds: gpstime((startTime + t*ds['cadence'])/ds['timeDenominator'])
        counter = dict(zip(('EphEpochs','ObsEpochs','Missed','Unlocked'),[0]*4))

        for prn in sorted(svdict):
            if not svdict[prn]:
                self.addRow(PRN=prn)
                self.addLineBreak()
                continue
            if VERBOSE > 1: print "Processing PRN%02d..." % int(prn)
            svstatus = svdict[prn][0]
            startTime = timegps(svstatus['startTime'])*svstatus['timeDenominator']
            
            sigrow = self.getRowDict('-')
            sigrow['PRN'] = prn
            sigrow['Type']="Eph"

            ephmask = (svstatus['mask'] & statusCodes['Ephemeris']) != 0
            nEphEpochs = sum(ephmask)
            if not summaryOnly:
                for start, stop in zip(*mtools.mask2runs(ephmask)):
                    sigrow['Epochs']    = stop-start
                    sigrow['PassStart'] = toffset(start,svstatus)
                    sigrow['PassStop']  = toffset(stop,svstatus)
                    self.addRow(**sigrow)
            sigrow['Type']="EphTotal"
            sigrow['Epochs']=nEphEpochs
            sigrow['PassStart'] = toffset(0,svstatus)
            sigrow['PassStop']  = toffset(len(ephmask),svstatus)
            counter['EphEpochs'] += nEphEpochs
            self.addRow(**sigrow)

            if svstatus['type'] != "OBSERVATIONS":
                self.addLineBreak()
                continue

            for ds in svdict[prn]:
                sigrow['Type'] = "Obs"
                sigrow['CC'] = ds['carrierCode']+ds['rangingCode']
                sigrow['Epochs'] = sigrow['Unlocked'] = 0
                sigrow['Missed'] = '-'

                obslockstatus = (ds['mask'] & reduce(op.or_,lockCodes.values()))
                obsmask = obslockstatus != 0
                if not summaryOnly:
                    for start, stop in zip(*mtools.mask2runs(obsmask)):
                        sigrow['Epochs']    = stop-start
                        sigrow['Unlocked']  = sum(obslockstatus[start:stop][ephmask[start:stop]]
                                                  != statusCodes['Code Carrier Lock'])

                        sigrow['PassStart'] = toffset(start,ds)
                        sigrow['PassStop']  = toffset(stop,ds)
                        self.addRow(**sigrow)

                nLocked = sum(obslockstatus == statusCodes['Code Carrier Lock'])
                nMissed = nEphEpochs - sum(obslockstatus[ephmask] == statusCodes['Code Carrier Lock'])
                nObsEpochs = sum(obsmask)
                sigrow['Type']      = "ObsTotal"
                sigrow['Epochs']    = nObsEpochs
                sigrow['Unlocked']  = nObsEpochs-nLocked
                sigrow['PassStart'] = ds['obsStartTime']
                sigrow['PassStop']  = ds['obsStopTime']
                sigrow['Missed']    = numpct(nMissed,nEphEpochs)
                self.addRow(**sigrow)

                counter['ObsEpochs'] += sigrow['Epochs']
                counter['Unlocked']  += sigrow['Unlocked']
                counter['Missed']    += nMissed
            self.addLineBreak()

        self.addRow(Type="EphTotal",Epochs=counter['EphEpochs'])
        self.addRow(Type="ObsTotal",Epochs=counter['ObsEpochs'],Unlocked=counter['Unlocked']
                    ,Missed=numpct(counter['Missed'],counter['EphEpochs']))
        return self


def connectToDB(dbHost="verdi", dbUsr='hrtr_monitor' ,dbPwd ='hrtr_monitor',dbName="events",dbPort=5432):
  """Function that allows a connection to the events, or test_events, database.
     keyword arguments:
       dbHost="borodin"
       dbUsr='hrtr_monitor'
       dbPwd ='{insert pwd}'
       dbName="events"
       dbPort=5432
     returns:
        psycopg2 connection object
  """
  # Try to connect to the database. Note that the named_pipe argument is for windows only
  try:
    conn = ps.connect(host=dbHost, user=dbUsr, password = dbPwd, dbname=dbName, port=dbPort)
  except ps.Warning, w:
    if VERBOSE: print 'Warning: %s'%(w)
  except ps.Error, e:
    if VERBOSE: print "Error %s" % (e.args[:])
    conn = None
  return conn
	

def executeQuery(conn,qstr):
  """Function that allows a connection to the events, or test_events, database.
     Arguments:
        conn - a psycopg2 connection object
        qstr - a valid query
  """
  if VERBOSE > 1: sys.stdout.write("Querying..."); sys.stdout.flush()
  if VERBOSE > 2: sys.stdout.write("QUERY:\n%s\n" % qstr); sys.stdout.flush()

  try:
    queryCursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    t1 = time.time()
    queryCursor.execute(qstr)
    t2 = time.time()
  except ps.Error, e:
    if VERBOSE: print "Query error:", e
    queryCursor.close()
    queryCursor = None
  else:
    if VERBOSE > 1: sys.stdout.write("%d rows returned in %g s.\n" % (queryCursor.rowcount,t2-t1))
  
  return queryCursor


def mask2nums(mask, startTime, cadence, timeDenominator, xoffset=0):
    runstarts, runstops = mtools.mask2runs(mask)
    runstarts = [dates.date2num(timedelta(seconds=(rstart+xoffset)*cadence/timeDenominator)+startTime) for rstart in runstarts]
    runstops = [dates.date2num(timedelta(seconds=(rstop+xoffset)*cadence/timeDenominator)+startTime) for rstop in runstops]

    return runstarts, runstops


def plot_status(startTime, stopTime, datadict, cadence, timeDenominator, keys=(1,0), values=(1,0)
                , colors=None, title=None, alpha=1.0, relsize=1.0, stack=True, vspace=0.05):
    '''
    plot_status

    Plot a collection of 1d integer valued arrays as horizontal bars (using matplotlib's barh plot).
    Plot style is defined by the keys, values, and colors.

    *datadict*
    Dictionary of dictionary items each containing a 'mask' and 'times' array.
    
    *cadence, timeDenominator*
    Uniform cadence/timeDenominator assumed for the data sets.
    
    *keys*
    Text label descriptors for status values.

    *colors*
    Distinct colors applied to each status value.

    *values*
    Distinct integer values corresponding to each status.

    *stack*
    Two ways to set this parameter:    
    1) Setting stack to True makes plot so that all unique values of an array are stacked on top of
    each other. When False, all values will be displayed at the same size and location (this works
    well if the status values are mutually exclusive).
    2) Alternately, you can specify an array of booleans that determines which values are stacked
    and which aren't. This can be useful for example if you want to plot one value underneath a
    stack of different values.

    *relsize*
    Two ways to set this parameter:
    1) A single value sets vertical size of _all_ status bars relative to 1.0 (should be <= 1).
    2) An array of values sets relative vertical size for each individual status value. Using this
    option together with the *stack* option gives some additional flexibility, allowing for example
    a dominant underlying status value with a stack of smaller sized status value bars fitting
    within it.

    *alpha*
    The alpha value to use for all status codes. Can also be specified as an array to set a unique
    alpha for each status code.

    *vspace*
    Amount of vertical space (relative to 1) to leave between horizontal bars.

    '''
    if colors is None:  colors = tuple(mpl.cm.jet(i/float(len(keys))) for i in range(len(keys)))

    if not hasattr(stack,'__iter__'):   stack = [stack]*len(keys)
    if not hasattr(alpha,'__iter__'):   alpha = [alpha]*len(keys)
    if not hasattr(relsize,'__iter__'): relsize = [relsize]*len(keys)

    fig = plt.figure()
    fig.set_size_inches(10.24,10.24)
    ax = fig.add_subplot(111)

    barHeight = 1.0-vspace # relative size alotted for each signal
    for i, (key,data) in enumerate( sorted(datadict.iteritems(),key=op.itemgetter(0)) ): # for each carrier
        keyoffset = i+1
        timeoffset = (data['startTime']-startTime).total_seconds()*timeDenominator/cadence

        j = 0
        for k,(s,v,a,scale) in enumerate(zip(stack,values,alpha,relsize)):        # plot each status code
            codeheight = barHeight*scale
            codeoffset = keyoffset - codeheight/2.

            if s:
                codeheight /= sum(stack)
                codeoffset += codeheight*float(j)
                j += 1

            starts, stops = mask2nums((data['mask'] & v) == v, startTime, cadence
                                      , timeDenominator, timeoffset)
            ax.barh([codeoffset]*len(starts), np.subtract(stops,starts), left=starts
                        , height=codeheight
                        , color=colors[k], linewidth=0, align='edge', alpha=a)

    # create a custom legend for status flags
    leg = plt.legend([plt.Line2D([],[], linewidth=5, color=c) for c in colors]
                     , keys, prop={'size':'medium'}, frameon=False
                     , bbox_to_anchor=(0, 1.015, 1., .05), mode='expand',ncol=len(keys))
    ax.add_artist(leg)

    ax.yaxis.set_ticks(np.arange(len(datadict))+1)
    ax.yaxis.set_ticklabels(sorted(datadict.keys()),fontsize=9,family='monospace')
    ax.set_ylim(0,len(datadict)+1)
    ax.set_xlabel('Date')
    # ax.set_ylabel(ylabel)
    ax.grid()

    if title is not None: fig.suptitle(title,fontsize=14,fontweight='bold')
    fig.autofmt_xdate()
    ax.set_xlim(startTime,stopTime)
    plt.setp(ax.get_xticklabels(),fontsize=9, rotation=15)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    fig.subplots_adjust(bottom=0.08,left=0.1,right=0.95)

    return fig


def GetSVSignals(conn,prn,codes):
    q = ("select * from gps_satellite_signals"
         + " where prn=%d"  % (prn) 
         + " and ("
         + " or ".join("(carrier_code='%s' and ranging_code='%s')" % (cc,rc) for cc,rc in codes)
         + " )")
    try:
        curs = executeQuery(conn,q)
        num = curs.fetchall()
        curs.close()
    except:
        num = []

    return num


def plot_signalavailability(startTime, stopTime, ephPasses, samplePeriod, title=None, obsPasses=None
                            , codes=[('L1','CA')]):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    clr1 = '#2c7fb8'
    clr2 = '#7fcdbb'
    clr3 = '#d95f0e'
    clr4 = '#fec44f'

    conn = connectToDB()
    numEpochs = int((stopTime-startTime).total_seconds()/samplePeriod)
    ephNumSigs = np.zeros(numEpochs, dtype=np.uint8)
    ephNumSats = np.zeros(numEpochs, dtype=np.uint8)
    obsNumSigs = np.zeros(numEpochs, dtype=np.uint8)
    obsNumSats = np.zeros(numEpochs, dtype=np.uint8)

    # Calculate number of satellites/signals available over given time interval
    for prn, ephdict in ephPasses.iteritems():
        ephOffset = (ephdict['startTime']-startTime).total_seconds()/samplePeriod
        nsvsigs = len(GetSVSignals(conn,prn,codes=codes)) if conn else 1

        ephNumSats[ephOffset:ephOffset+len(ephdict['mask'])] += ephdict['mask']
        ephNumSigs[ephOffset:ephOffset+len(ephdict['mask'])] += ephdict['mask']*nsvsigs

        if obsPasses:
            allmask = np.zeros(numEpochs, dtype=np.uint8)
            for codeTrack in obsPasses[prn]:
                obsOffset = (codeTrack['startTime']-startTime).total_seconds()/samplePeriod
                demodSlice = slice(obsOffset, obsOffset+len(codeTrack['demodStatus']))
                allmask[demodSlice] |= (codeTrack['demodStatus'] == CODECARRIERLOCK)
                obsNumSigs[demodSlice] += (codeTrack['demodStatus'] == CODECARRIERLOCK)
            obsNumSats += allmask
    if conn: conn.close()

    if obsPasses:
        # plot obs signal availability
        obsstarts, obsstops = mtools.mask2runs(obsNumSigs > 0)
        for start, stop in zip(obsstarts,obsstops):
            times = dates.drange(timedelta(seconds=start*float(samplePeriod))+startTime
                                 , timedelta(seconds=stop*float(samplePeriod))+startTime
                                 , timedelta(seconds=samplePeriod))
            ax.plot(times, obsNumSigs[start:stop],color=clr3, drawstyle='steps-post')
            ax.plot(times, obsNumSats[start:stop],color=clr4, drawstyle='steps-post')            

    # plot ephemeris signal availability            
    ephstarts, ephstops = mtools.mask2runs(ephNumSigs > 0)
    for start, stop in zip(ephstarts,ephstops):
        times = dates.drange(timedelta(seconds=start*float(samplePeriod))+startTime
                             , timedelta(seconds=stop*float(samplePeriod))+startTime
                             , timedelta(seconds=samplePeriod))
        ax.plot(times, ephNumSigs[start:stop],color=clr1, drawstyle='steps-post')
        ax.plot(times, ephNumSats[start:stop],color=clr2, drawstyle='steps-post')        

    # create a custom legend
    clrs = [clr1,clr2] + ([clr3,clr4] if obsPasses else [])
    types = ['Ephemeris'] + (['Obs'] if obsPasses else [])
    labels = ['Ephemeris Signals Available', 'Ephemeris Satellites Available'] \
                + (['Obs Signals Available', 'Obs Satellites Available'] if obsPasses else [])
    legmarkers = tuple(plt.Line2D([],[], linewidth=5, color=c) for c in clrs)
    ax.legend(legmarkers, labels, loc='best', prop={'size':'small'})

    ax.set_ylim(0,max(ephNumSigs)+1)
    ax.set_xlabel('Date')
    ax.set_ylabel('Number of signals available')
    ax.grid()

    fig.autofmt_xdate()
    ax.set_xlim(startTime,stopTime)
    plt.setp(ax.get_xticklabels(),fontsize=9, rotation=15)      
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())    
    if title is not None: ax.set_title(title,fontsize=14,fontweight='bold')
    fig.tight_layout()

    return fig


# TODO
# + write C++ wrapper for SVAvailability to speed it up
def SVAvailability(eph, rxPos, sv, startTime=None, stopTime=None, elmask=10., samplePeriod=None):
    epoch = eph.getClockTimeStep(sv)
    if startTime is None or startTime < (eph.getInitialTime(sv)+epoch*50):
        currentTime = eph.getInitialTime(sv)+epoch*50
    else:
        currentTime = gpstk.CommonTime(startTime)
    if stopTime is None or stopTime > (eph.getFinalTime(sv)-epoch*50):
        stopTime = eph.getFinalTime(sv)-epoch*50
    if not samplePeriod: samplePeriod=epoch

    # skip to where the satellite is above elevation angle
    while(1):
        svXVT = eph.getXvt(sv,currentTime)
        el = rxPos.elvAngle(svXVT.x)

        # base case:
        # if we've skipped past the stop time while searching for next
        # time sv is above elevation angle then previous interval
        # was the last valid interval
        if currentTime >= stopTime: return []

        if el >= elmask: break
        currentTime.addSeconds(samplePeriod)

    # iterate through until satellite is below elevation angle
    start = gpstk.CommonTime(currentTime)
    while(currentTime < stopTime):
        svXVT = eph.getXvt(sv,currentTime)
        if rxPos.elvAngle(svXVT.x) < elmask: break
        currentTime.addSeconds(samplePeriod)
    else:
        # base case:
        # if we've skipped past the stop time this is the last valid interval
        return [(start,gpstk.CommonTime(stopTime))]
    stop = gpstk.CommonTime(currentTime)

    return [(start,stop)] + SVAvailability(eph, rxPos, sv, startTime=stop, stopTime=stopTime, elmask=elmask, samplePeriod=samplePeriod)


def ProcessEphemerisAvailability(ephemeris, sv, cadence, timeDenominator, rxPos
                                 , elmask=10., startTime=None, stopTime=None):
    # Collect data from all passes from this SV during the specified time in the ephemeris
    ephPasses = SVAvailability(ephemeris, rxPos, sv, startTime=startTime
                               , stopTime=stopTime, elmask=elmask)
    if not ephPasses:
        return None, None, np.array([]), None

    # Allocate mask for all passes (including in between passes)
    ephGroupStart_ct, ephGroupStop_ct = min(zip(*ephPasses)[0]), max(zip(*ephPasses)[1])
    ephAboveMask = np.zeros(long((ephGroupStop_ct - ephGroupStart_ct)*timeDenominator)/cadence,dtype=np.bool)

    # Map each pass into the 'global' ephemeris mask
    for ephPassStart_ct, ephPassStop_ct in ephPasses:
        startidx = long((ephPassStart_ct - ephGroupStart_ct)*timeDenominator)/cadence
        stopidx = long((ephPassStop_ct - ephGroupStart_ct)*timeDenominator)/cadence
        ephAboveMask[startidx:stopidx] = True

        # get elv at two points in time
        midTime = ephPassStart_ct + (ephPassStop_ct-ephPassStart_ct)/2.
        xvt = ephemeris.getXvt(sv,midTime)
        el0 = rxPos.elvAngle(xvt.x)
        xvt = ephemeris.getXvt(sv,midTime+2*ephemeris.getClockTimeStep(sv))
        el1 = rxPos.elvAngle(xvt.x)
            
    return ephGroupStart_ct, ephGroupStop_ct, ephAboveMask, (el0,el1)


def ProcessDataAvailibility(ephemeris, rxPos, cadence, timeDenominator, startTime=None, stopTime=None
                            , mdh_file=None, elmask=10., codeCarriers=None, exsats=None):

    minStartTime_ct = mtools.datetime2commontime(startTime) if startTime else ephemeris.getInitialTime()
    maxStopTime_ct = mtools.datetime2commontime(stopTime) if stopTime else ephemeris.getFinalTime()

    svStatus = {}
    for sv in ephemeris.getSatList():
        if exsats and sv.id in exsats: continue

        start, stop, ephmask, elv = ProcessEphemerisAvailability(ephemeris, sv, cadence, timeDenominator
                                                                 , rxPos, elmask, startTime=minStartTime_ct
                                                                 , stopTime=maxStopTime_ct)
        if not ephmask.size: continue
        
        statmask = np.zeros(ephmask.shape,dtype=np.uint32)
        statmask[ephmask] = statusCodes['Ephemeris']
        statmask[~ephmask] = UNAVAILABLE
        del ephmask
        
        ephGroupStart = mtools.commontime2mdhtime(start,timeDenominator)
        ephGroupStop = mtools.commontime2mdhtime(stop,timeDenominator)

        svStatus[sv.id] = []
        if not mdh_file: # only add this if there is no obs, otherwise it is redundant
            svStatus[sv.id] = [{'startTime': mtools.commontime2datetime(start)
                                ,'cadence': cadence,'timeDenominator': timeDenominator
                                ,'mask':statmask,'elevation': elv,'type':"EPHEMERIS"}]
            continue

        for (prn,cc,rc) in mtools.matchingSignals(mdh_file, SUBCLASSES=["OBSERVATIONS"]
                                                  , prnCodes=[sv.id], codeCarriers=codeCarriers):
            if VERBOSE > 1: print "Processing PRN%02d %2s%2s..." % (prn,cc,rc)
            obsGroupStart = obsGroupStop = None
            svmask = statmask.copy()
            for obsTrack in mdh_file.matchingDatasets(SUBCLASS="OBSERVATIONS", prnCode=prn, carrierCode=cc, rangingCode=rc):
                obsPassStart = max(obsTrack.attrs['startTime'], ephGroupStart)
                obsPassStop = min(mtools.getStopTime(obsTrack), ephGroupStop)
                if obsPassStop <= obsPassStart: continue

                if not obsGroupStart or obsPassStart < obsGroupStart: obsGroupStart = obsPassStart
                if not obsGroupStop or obsPassStop > obsGroupStop: obsGroupStop = obsPassStop

                nEpochs = (obsPassStop-obsPassStart)/obsTrack.attrs['cadence']
                ephoffset = (obsPassStart - ephGroupStart)/cadence
                trkoffset = (obsPassStart - obsTrack.attrs['startTime'])/obsTrack.attrs['cadence']

                # convert code to a unique bit position
                demodStatus = obsTrack['demodulatorStatus'][trkoffset:trkoffset+nEpochs]
                for key,lockcode in lockCodes.iteritems():
                   svmask[ephoffset:ephoffset+nEpochs] |= (demodStatus == lockcode)*statusCodes[key]

            # save the avaialability data for this code carrier
            svStatus[prn].append({'mask': svmask
                                  , 'startTime': gpstime(ephGroupStart/timeDenominator)
                                  , 'obsStartTime': gpstime(obsGroupStart/timeDenominator)
                                  , 'obsStopTime': gpstime(obsGroupStop/timeDenominator)
                                  , 'rangingCode': rc, 'carrierCode': cc
                                  , 'cadence': cadence, 'timeDenominator': timeDenominator
                                  , 'elevation': elv, 'type': 'OBSERVATIONS'})
    return {k:sorted(stat, key=op.itemgetter('startTime')) for k,stat in svStatus.iteritems()}


if __name__ == '__main__':
    # ============================================================
    # Define CLI
    # ============================================================
    def position(value):
        return gpstk.Position(*map(float, value.split(",")))
    def codecarrierTuples(value):
        return tuple(value.split(':'))

    from argparse import ArgumentParser
    parser = ArgumentParser(description="Data availability analysis for ephemeris"
                            +" and MDH obs files. If obs is passed in, results give a comparison"
                            +" between obs availability and availability predicted from provided ephemeris file.")

    parser.add_argument("-e", "--eph",dest='ephemeris',required=True,default=None
                        ,action='append'
                        ,help=("Path to ephemeris file. Accepts multiple"
                               +" ephemeris files when --eph is specified multiple times."))

    parser.add_argument("-p","--position", dest='rxpos', type=position, required=True
                        , help="Receiver X,Y,Z position.")

    parser.add_argument("--obs", dest='obs', default=None,action='append'
                        , help="Path to MDH observations file. (%(default)s)")

    parser.add_argument("--start", dest="start", default=None
                        , help="Consider only data after this time. (max(ObsStart,EphStart))")

    parser.add_argument("--stop", dest="stop", default=None
                        , help="Consider only data before this time. (min(ObsEnd,EphEnd)")

    parser.add_argument("--elmask", dest='elmask', default=10., type=float
                        , help="Minimum elevation angle in degrees. (%(default)s)")

    # parser.add_argument("--decimate", dest="decimate", default=1, type=int
    #                   , help="Decimate obs data by a factor. (1)")

    parser.add_argument("--exSats", dest="exsats", default=None, nargs='*', type=int
                        , help="Specify SVs to exclude from analysis. (%(default)s)")

    parser.add_argument("--savepath", dest='savepath', default=None
                        ,help="Save directory for plots. (%(default)s/Interactive plotting)")

    parser.add_argument("--log", dest="log", default='/dev/stdout'
                      ,help="Where to output analysis. (%(default)s)")

    parser.add_argument("--latex",dest='latex',action='store_true',default=False
                      ,help="Output report in latex. (%(default)s)")

    parser.add_argument("--detailed",dest='detailedReport',action='store_true',default=False
                        ,help="Give a detailed report - include results for each pass. (%(default)s)")

    parser.add_argument("--plot-signals-count",dest='plt_signalscount'
                        ,default=False, nargs='?', const=True
                        ,help=("Generate plot displaying number of available signals over time."
                               +" Optionally, pass in a custom filename (relative to savepath)."
                               +" (%(default)s) CURRENTLY UNAVAILABLE"))

    parser.add_argument("--plot-availability",dest='plt_availability'
                      ,default=False, const=True, nargs='?'
                      ,help=("Generate plot displaying satellite status over time."
                             +" Optionally, pass in a custom filename (relative to savepath). (%(default)s)"))

    parser.add_argument("-q","--quiet",dest='quiet',action='store_true',default=False
                      ,help="Disable verbose output. (%(default)s)")

    parser.add_argument("-d", "--debug", default=0, dest="debug", action="count"
                      ,help="Increase the debug output level. (%(default)s)")

    parser.add_argument("--cc", "--code-carriers", dest="codecarriers"
                        , type=codecarrierTuples, nargs='*', default=None
                        , help="Restrict data to specific carrier:code combination(s). "
                        + "Pass multiple signals in as a comma delimited string. "
                        + "(ALL available signals)")

    parser.add_argument("--timefmt", dest="timefmt", default="%Y-%m-%d %H:%M:%S",
                      help="Set time format string. (%(default)s)")

    parser.add_argument("-v", dest="verbose", action="count",default=1
                      ,help=("Print status information to command line. "
                             +"Enter multiple 'v's for more verbosity. (v)"))

    opts = parser.parse_args()
    
    # ==========================================================
    # Preprocessing, configuration
    # ==========================================================
    VERBOSE = 0 if opts.quiet else opts.verbose
    DEBUG = opts.debug
    mtools.VERBOSE = VERBOSE

    # Parse input
    if opts.start: opts.start = mtools.parseDate(opts.start,opts.timefmt)
    if opts.stop: opts.stop = mtools.parseDate(opts.stop,opts.timefmt)

    # Read in all ephemeris
    eph = mtools.EphReader(*opts.ephemeris)

    # Set analysis time and assumed cadence
    if opts.obs:
        filenames = filter(h5py.is_hdf5, [fn for path in opts.obs for fn in glob.glob(path)])
        mdh_fh = mdh.MDHFileList(*filenames,mode='r')
        minStartTime, maxStopTime = map(gpstime,mtools.mdhSpan(mdh_fh,absolute=True,SUBCLASS="OBSERVATIONS"))

        attrs = []
        obs = mdh_fh.matchingDatasets(SUBCLASS="OBSERVATIONS")
        while 'cadence' not in attrs: attrs = dict(obs.next().attrs)
        if not attrs: sys.exit("Error: No 'cadence' attribute found in datasets")
        timeDenominator = attrs['timeDenominator']
        cadence = attrs['cadence']
    else:
        minStartTime = commontime2datetime(eph.getInitialTime())
        maxStopTime = commontime2datetime(eph.getFinalTime())
        cadence = timeDenominator = 1

    if opts.start: minStartTime = max(minStartTime,opts.start)
    if opts.stop:  maxStopTime = min(maxStopTime,opts.stop)

    # force start/stop times to second boundaries
    roundSec = lambda t,r=int: int(r(t.second + 1e-6*t.microsecond))
    minStartTime = minStartTime.replace(microsecond=0, second=roundSec(minStartTime,np.ceil))
    maxStopTime = maxStopTime.replace(microsecond=0, second=roundSec(maxStopTime))
    
    if VERBOSE: print "Analysis span: %s through %s" % (minStartTime, maxStopTime)        

    # ==========================================================
    # Run tasks
    # ==========================================================    
    # Generate availability report
    if VERBOSE: print "Extracting availability data..."
    svStatus = ProcessDataAvailibility(eph, opts.rxpos, cadence, timeDenominator
                                       , startTime=minStartTime, stopTime=maxStopTime
                                       , mdh_file=mdh_fh if opts.obs else None
                                       , elmask=opts.elmask
                                       , codeCarriers=opts.codecarriers
                                       , exsats=opts.exsats)

    if VERBOSE: print "Generating report..."
    with open(opts.log,'w') as logger:
        fmt = LatexTable if opts.latex else TextTable
        keys = AvailabilityReport.keys
        keys.remove('Elv')
        
        report = AvailabilityReport(formatter=fmt, keys=keys)
        report.generateContent(svStatus, summaryOnly=not opts.detailedReport)
        report.writeReport(out=logger)

    # Generate plots
    if opts.plt_availability:
        #flatten the list of statuses to plot
        flatstat = {}
        for prn in svStatus:
            svkey = "PRN%02d" % prn
            for ds in svStatus[prn]:
                if ds['type'] == 'OBSERVATIONS':
                    key = svkey+(" %(carrierCode)2s%(rangingCode)2s" % ds)
                else:
                    key = svkey
                flatstat[key.ljust(10)] = ds

        genericPlot(plot_status, (flatstat, cadence, timeDenominator), "Data Availability"
                    , minStartTime, maxStopTime, opts.savepath
                    , filetag = opts.plt_availability if isinstance(opts.plt_availability,str) else None
                    , colors = statusColors if opts.obs else [statusColors[0],]
                    , keys = statusCodes.keys() if opts.obs else [statusCodes.keys()[0],]
                    , values = statusCodes.values() if opts.obs else [statusCodes.values()[0],]
                    , stack = False
                    , relsize = [1.0]+[0.8]*(len(statusCodes)-1) if opts.obs else 1.0
                    , alpha = [0.5]+[1.0]*(len(statusCodes)-1) if opts.obs else 0.5
                    , vspace = 0.1)

    # if opts.plt_signalscount:
    #     genericPlot(plot_signalavailability, (ccStatus, samplePeriod)
    #                 , "Signal Availability"
    #                 , minStartTime, maxStopTime, opts.savepath
    #                 , filetag = opts.plt_signalscount if isinstance(opts.plt_signalscount,str) else None
    #                 , codes=opts.codecarriers
    #                 , obsPasses=ccStatus)
