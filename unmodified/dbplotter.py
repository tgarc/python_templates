#! /usr/bin/python

'''
dbplotter

Command line interface for visualizing data from the HRTR events database

Tomas Garcia
'''

import os
import sys
import time
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib import dates
import mdh_tools as mtools
from host import mdh
from optparse import OptionParser
import dbtools
import mdh_plotevents as mplot

VERBOSE = 0


def genericQuery(queryfun,conn,form,*args,**kwargs):
  qlist,qstr,wc = queryfun(conn,form,*args,**kwargs)

  if not qlist:
    if VERBOSE: sys.stdout.write("\nFailed.\n")
    return None,qstr,wc

  if VERBOSE: print "filter: \n%s\n" % (wc)

  return qlist,qstr,wc


def genericLog(dictlist,title,form,savedir,fieldfmt,qstr='',wc='',fmt="%Y-%m-%d %H:%M:%S"):
  if VERBOSE: mtools.printHeader(title)
  if VERBOSE: sys.stdout.write("Generating log '%s'..." % title); sys.stdout.flush()

  if savedir:
    s = datetime.strptime(form['startdate'],fmt)
    e = datetime.strptime(form['enddate'],fmt)
    filetag = '_'.join( title.lower().split() )
    filetag += '_'.join(map(str, (s,e)))

    log = mtools.LogFile(os.path.join(savedir,filetag).replace(' ','_'),'w',logfmt=fieldfmt
                         ,start_date=s, end_date=e)
    log.addComment("Query 'where' clause" +'\n' + wc)

    log.writeDictList(dictlist)
    log.close()
  else:
    print "\n".join(", ".join(map(str,r.values())) for r in dictlist)
  if VERBOSE: print

parser = OptionParser(usage="./dbplotter.py [options] --savepath=/save/path")

parser.add_option("--savepath", dest='savepath',type='string', default=None,
                  help="Save plot(s) / log file(s) to the directory path specified.")

parser.add_option("--eventslog",dest='eventslog',action='store_true',default=False
                  ,help="Save an events log file.")

parser.add_option("--trackslog",dest='trackslog',action='store_true',default=False
                  ,help="Save total track time log file.")

parser.add_option("--plotall",dest='plotall',action='store_true',default=False
                  ,help="For each event, generate plots for all recorded tracks on the event's corresponding MDH file.")

parser.add_option("-r", "--plotreference",dest='plotreference',action='store_true',default=False
                  ,help="For each plot, overlay a plot of the highest SNR sattelite at event time.")

parser.add_option("-q","--quiet",dest='quiet',action='store_true',default=False
                ,help="Disable verbose output.")

parser.add_option("-v", dest="verbose", action='count',default=0
                  , help="Print status information to command line. Enter multiple 'v's for more verbosity.")

parser.add_option("--unclustered", dest="unclustered", action='store_true'
                  , help="Don't cluster events over multiple detection periods.")

parser.add_option("--start", dest="start", default=None, help="Consider events after this date time")

parser.add_option("--stop", dest="stop", default=None, help="Consider events before this date time")

parser.add_option("--svn", dest="svn", default='', help="Filter on SVN.")

parser.add_option("--sv_type", dest="svType", default='', help="Filter on satellite block type (e.g., IIF).")

parser.add_option("--cc", "--code-carriers", dest="codeCarriers", default=''
                  , help="Filter on carrier:ranging code combination. Enter multiple signals as a comma delimited list (e.g., --code-carriers L1:CA,L2:CM)")

parser.add_option("--event_status", dest="event_status", default=None, help="Filter on event status.")

parser.add_option("--min_signal_strength", dest="min_signal_strength", default=None, type='float'
                  ,help="Filter out events with signal strength below a set value.")

parser.add_option("--event_id", dest="event_id", default='' 
                  ,help="Choose specific events.")

parser.add_option("--max_detection_value", dest="max_detection_value", default=None, type='float'
                  ,help="Filter out events with phase event detection metric above a set value.")

parser.add_option("--max_threshold", dest="max_threshold", default=None, type='float'
                  ,help="Filter out events with threshold above a set value.")

parser.add_option("--timefmt", dest="timefmt", default="%Y-%m-%d %H:%M:%S"
                  ,help="Set input time format string. (DEFAULT=\"%s\")" % "%Y-%m-%d %H:%M:%S")

parser.add_option("--eventslogfmt", dest="evtlogfmt",
                  default=mtools.EVENTSLOG_DEFAULT_FMT
                  ,help="Set output format of events log.  Note that the bracket"
                  +"{} style format strings allow for more flexibility in the"
                  +"output.For example, it allows one to specify the output format"
                  +"of dates like so \"{startdate:%%Y-%%m-%%d"
                  +"%%H:%%M:%%S}\". (DEFAULT = \"%s\")" % mtools.EVENTSLOG_DEFAULT_FMT)

parser.add_option("--trackslogfmt", dest="trklogfmt", default=mtools.TRACKSLOG_DEFAULT_FMT
                  ,help="Set output format of track log. (DEFAULT = \"%s\")" % mtools.TRACKSLOG_DEFAULT_FMT)


mplot.addPlotOptions(parser)
opts,args = parser.parse_args()

if opts.savepath and not (os.path.exists(opts.savepath) and os.path.isdir(opts.savepath)):
  sys.exit("Error: savepath does not exist or is not a directory.")

VERBOSE = 0 if opts.quiet else (opts.verbose if opts.verbose else 1)
mtools.VERBOSE = dbtools.VERBOSE = mplot.VERBOSE = VERBOSE

form = {}
form['event_type'] = ['phase'] # fixed for now

### Connect to DB
conn = dbtools.connectToDB()

if opts.event_status:
  form['event_status'] = opts.event_status.split(',')
if opts.start:
  form['startdate'] = str(mtools.parseDate(opts.start,opts.timefmt))
if opts.stop:
  form['enddate'] = str(mtools.parseDate(opts.stop,opts.timefmt))
if opts.svn:
  form['svn'] = [int(svn) for svn in opts.svn.split(',')]
if opts.min_signal_strength:
  form['min_signal_strength'] = opts.min_signal_strength
if opts.codeCarriers:
  form['carrier_code'],form['ranging_code'] = zip(*[s.split(':') for s in opts.codeCarriers.split(',')])
if opts.svType:
  form['sv_type'] = [svt for svt in opts.svType.split(',')]  
if opts.event_id:
  form['event_id'] = [int(v) for v in opts.event_id.split(',')]  
if opts.max_threshold:
  form['max_threshold'] = opts.max_threshold
if opts.max_detection_value:
  form['max_detection_value'] = opts.max_detection_value

start_date = mtools.parseDate(opts.start,opts.timefmt) if opts.start else None
end_date = mtools.parseDate(opts.stop,opts.timefmt) if opts.stop else None

eventslist = None
hblist = None

# make necessary queries
if any((opts.plothb, opts.plottimelineplushb, opts.trackslog, opts.plothbcum)):
  if VERBOSE: mtools.printHeader("Query for heartbeats data")
  hblist,hbqstr,hbwc = genericQuery(dbtools.executeHeartbeatsQuery,conn,form)
  
if any((opts.ploteventcount, opts.plottimeline, opts.plottimelineplushb
       ,opts.eventslog, opts.plotevents)):
  if VERBOSE: mtools.printHeader("Query for events data")  
  eventslist,eqstr,ewc = genericQuery(dbtools.executeEventsQuery,conn,form,clusterEvents=not opts.unclustered)

# make any requested aggregate plots
if opts.plothbcum:
  mplot.genericPlot(mplot.plot_hbcum,hblist,opts.savepath
                    ,"Cumulative heartbeat presence"
                    ,start_date=start_date,end_date=end_date)
if opts.plothb:
  mplot.genericPlot(mplot.plot_hb,hblist,opts.savepath,'Heartbeat presence'
                    ,start_date=start_date,end_date=end_date, timeseries=True)
if opts.ploteventcount:
  mplot.genericPlot(mplot.plot_eventcount,eventslist,opts.savepath,"Event count"
                    ,start_date=start_date,end_date=end_date)
if opts.plottimeline:
  mplot.genericPlot(mplot.plot_timeline,eventslist,opts.savepath,'Event timeline'
                    ,start_date=start_date,end_date=end_date, timeseries=True)
if opts.plottimelineplushb:
  mplot.genericPlot(mplot.plot_timeline,eventslist,opts.savepath
                    ,'Event timeline and heartbeat'
                    ,start_date=start_date,end_date=end_date
                    ,timeseries=True,**{'hblist':hblist})

# Write out log files containing query results
if opts.trackslog:
  genericLog(hblist,'tracks',form,opts.savepath,opts.trklogfmt,qstr=hbqstr
             ,wc=hbwc,fmt=opts.timefmt)
if opts.eventslog:
  genericLog(eventslist,'events',form,opts.savepath,opts.evtlogfmt,qstr=eqstr
             ,wc=ewc,fmt=opts.timefmt)

# Plot out each event returned by query
if opts.plotevents:
  if VERBOSE: mtools.printHeader("plotevents")

  for r in eventslist:
    if not r['file_path']:
      if VERBOSE: print "Could not find matching dataset for event ID %d." % r['event_id']
      continue

    if VERBOSE: sys.stdout.write("\rProcessing EventID: %(event_id)d..." % r); sys.stdout.flush()

    fig = mplot.plot_event(r, savedir=opts.savepath, plotall=opts.plotall
                           , style=opts.style, plotreference=opts.plotreference)[0]

    if not opts.savepath:
      if VERBOSE: sys.stdout.write("Close window to continue..."); sys.stdout.flush()
      plt.show(fig.number)

    # if opts.plotreference:
    #   t = "Reference SV for Event ID: %d" % r['event_id']
    #   ftag = mplot.FileTagGen(event_date=r['event_date']
    #                           ,prn=r['prn']
    #                           ,carrier_code=r['carrier_code']
    #                           ,ranging_code=r['ranging_code']
    #                           ,event_id=r['event_id']) + "_refsv"
    #   mplot.plot_event(r,savepath,plotall=opts.plotall,style=opts.style,plotreference=opts.plotreference
    #                    ,title=t,filetag=ftag)

  if VERBOSE: print

# Done!
conn.close()

