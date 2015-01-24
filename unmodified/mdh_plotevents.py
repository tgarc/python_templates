#! /usr/bin/python

"""
Helper tools for plotting anomalous phase events captured in mdh files

TODO
+ simplify plot_event
+ fix some hacky stuff with options required for format plot functions
Thomas J Garcia
"""
import matplotlib.dates as dates
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import scipy as sp
import os
import calendar, time
from datetime import datetime,timedelta
from host import signals,mdh
import mdh_tools as mtools
from mdh_tools import isvalid
import sys
from itertools import groupby
import operator as op
from operator import itemgetter, attrgetter

### Constants
cvals = {'L1:CA':'b','L2:CM':'r','L2:CL':'c','L5:I5':'g','L5:Q5':'goldenrod'}
VERBOSE = 0
PRINT_SNR_REPORT = False #if true will add SNR to report style plots 


def addPlotOptions(parser):
  parser.add_option("--style", dest='style', default='std',type='str'
                    ,help="Specify events plotting style. (default=\'std\') \n \
                          'std': Standard style plot. \n \
                          'std2': Like standard but larger detrended phase plot. \n \
                          'long': Long time scale style plot. \n \
                          'report': Like std but plots CA,CM, and I5 on the same plot")

  parser.add_option("--plot-events",dest='plotevents',action='store_true',default=False
                  ,help="Generate plots for events.")

  parser.add_option("--plot-event-count",dest='ploteventcount',action='store_true',default=False
                    ,help="Generate histogram plot of aggregate SV event counts.")

  parser.add_option("--plot-timeline",dest='plottimeline',action='store_true',default=False
                    ,help="Generate event timeline plot.")

  parser.add_option("--plot-timeline-with-heartbeat",dest='plottimelineplushb',action='store_true',default=False
                    ,help="Generate event timeline plot overlayed on heartbeat plot")

  parser.add_option("--plot-event-rates",dest='ploteventrates',action='store_true',default=False
                    ,help="Generate plot of daily event rates.")

  parser.add_option("--plot-heartbeat",dest='plothb',action='store_true',default=False
                    ,help="Plot heartbeat presence for a given date range.")

  parser.add_option("--plot-heartbeat-cum",dest='plothbcum',action='store_true',default=False
                    ,help="Generate cumulative heartbeat presence plot.")

  
def genericPlot(plotter,dictlist,savedir=None,title='',filetag='',start_date=None,end_date=None
                ,system='GPS',fmt="%Y%m%d",timeseries=False,**plotargs):
  if VERBOSE:
    sys.stdout.write("Generating plot '%s'..." % title); sys.stdout.flush()

  if timeseries:
    plotargs['start_date'], plotargs['end_date'] = start_date, end_date
  fig = plotter(dictlist,title=title,**plotargs)

  if not fig:
    if VERBOSE: sys.stdout.write("Failed.\n")
    return

  if not filetag:
    s = start_date.strftime(fmt) if start_date else ''
    e = end_date.strftime(fmt) if end_date else ''
    filetag = '_'.join( title.lower().split(' ') )
    filetag += '_' if start_date or end_date else ''
    if s and e: filetag += '_'.join((s,e))
    elif s or e: filetag += '_' + (s if start_date else e)

  if savedir:
    path = os.path.join(savedir,filetag)+'.png'
    if not os.path.exists(path):
      fig.savefig(path,dpi=100)
    elif VERBOSE == 0 or mtools.getCorrectInput("%s already exists. Replace? (y/n, default=y): " % path) in ('','y'):
      fig.savefig(path,dpi=100)
  else:
    if VERBOSE: sys.stdout.write("Close window to continue..."); sys.stdout.flush()
    plt.show(fig.number)
    if VERBOSE: print
  del(fig)


def FileTagGen(event_date=None,prn='',carrier_code='',ranging_code='',event_id=None,date_fmt="%Y%m%d_%H%M%S"):
  return '_'.join( [event_date.strftime(date_fmt) if event_date else ''
                    ,mtools.CodeIDGen(prn=prn,carrier_code=carrier_code,ranging_code=ranging_code)
                    ,("EID%u" % event_id if event_id else '')] )


def setup_stdPlot(fig,figsize=(12.96,7.2)):
  ## Create plot container
  fig.set_size_inches(*figsize)

  axes = {}
  axes['phase']=fig.add_subplot(511)
  # axes['ref_phase']=fig.add_subplot(612,sharex=axes['phase'])
  axes['detection_metric']=fig.add_subplot(512,sharex=axes['phase'])
  axes['correlation_counts']=fig.add_subplot(513,sharex=axes['phase'])
  axes['nco']=fig.add_subplot(514,sharex=axes['phase'])
  axes['snr']=fig.add_subplot(515,sharex=axes['phase'])
  fig.subplots_adjust(top=0.9,left=0.095,right=0.95,hspace=0.3)

  return axes

def setup_std2Plot(fig,figsize=(12.96,7.2)):
  ## Create plot container
  fig.set_size_inches(*figsize)

  axes = {}
  axes['phase']=fig.add_subplot(211)
  # axes['_ref_phase=fig.add_subplot(612,sharex=axis_phase)
  axes['detection_metric']=fig.add_subplot(714,sharex=axes['phase'])
  axes['correlation_counts']=fig.add_subplot(715,sharex=axes['phase'])
  axes['nco']=fig.add_subplot(716,sharex=axes['phase'])
  axes['snr']=fig.add_subplot(717,sharex=axes['phase'])  

  fig.subplots_adjust(top=0.9,left=0.095,right=0.95,hspace=0.3)

  return axes

def setup_reportPlot(fig):

  if not PRINT_SNR_REPORT:
    fig.set_size_inches(12.96,7.2)
  else:
    fig.set_size_inches(12.96,9.0)
  
  axes = {}
  if not PRINT_SNR_REPORT:
    axes['phase']=fig.add_subplot(411);
    axes['detection_metric']=fig.add_subplot(412,sharex=axes['phase'])
    axes['corrcount_subpane']=fig.add_subplot(413)
    axes['nco_subpane']=fig.add_subplot(414)
  else:
    axes['phase']=fig.add_subplot(511);
    axes['detection_metric']=fig.add_subplot(512,sharex=axes['phase'])
    axes['corrcount_subpane']=fig.add_subplot(513)
    axes['nco_subpane']=fig.add_subplot(514)
    axes['SNR_subpane']=fig.add_subplot(515)
  
  axes['corrcount_subpane'].set_frame_on(False);
  axes['nco_subpane'].set_frame_on(False)
  axes['nco_subpane'].get_xaxis().set_ticks([])
  axes['nco_subpane'].get_yaxis().set_ticks([])
  axes['corrcount_subpane'].get_xaxis().set_ticks([])
  axes['corrcount_subpane'].get_yaxis().set_ticks([])
  
  if not PRINT_SNR_REPORT:
    axes['correlation_counts']={};axes['nco']={}
    axes['correlation_counts']['L1:CA']=fig.add_subplot(12,1,7,sharex=axes['phase'])
    axes['correlation_counts']['L2:CM']=fig.add_subplot(12,1,8,sharex=axes['phase'])
    axes['correlation_counts']['L5:I5']=fig.add_subplot(12,1,9,sharex=axes['phase'])
    axes['nco']['L1:CA']=fig.add_subplot(12,1,10,sharex=axes['phase'])
    axes['nco']['L2:CM']=fig.add_subplot(12,1,11,sharex=axes['phase'])
    axes['nco']['L5:I5']=fig.add_subplot(12,1,12,sharex=axes['phase'])
  else:
    axes['correlation_counts']={};axes['nco']={}
    axes['correlation_counts']['L1:CA']=fig.add_subplot(15,1,7,sharex=axes['phase'])
    axes['correlation_counts']['L2:CM']=fig.add_subplot(15,1,8,sharex=axes['phase'])
    axes['correlation_counts']['L5:I5']=fig.add_subplot(15,1,9,sharex=axes['phase'])
    axes['nco']['L1:CA']=fig.add_subplot(15,1,10,sharex=axes['phase'])
    axes['nco']['L2:CM']=fig.add_subplot(15,1,11,sharex=axes['phase'])
    axes['nco']['L5:I5']=fig.add_subplot(15,1,12,sharex=axes['phase'])
 
    axes['SNR_subpane'].set_frame_on(False)
    axes['SNR_subpane'].get_xaxis().set_ticks([])
    axes['SNR_subpane'].get_yaxis().set_ticks([])
    axes['snr']={}
    axes['snr']['L1:CA']=fig.add_subplot(15,1,13,sharex=axes['phase'])
    axes['snr']['L2:CM']=fig.add_subplot(15,1,14,sharex=axes['phase'])
    axes['snr']['L5:I5']=fig.add_subplot(15,1,15,sharex=axes['phase'])

  fig.subplots_adjust(top=0.88,left=0.1,bottom=0.08,right=0.95,hspace=0.75)

  return axes

'''
format_stdPlot

Apply plot labels, legend, and timestamp to stdPlot

'''

def format_stdPlot(fig,axes,eventTime):
  axes['phase'].legend(bbox_to_anchor=(0., 1.0, 1., .05),loc='lower left'
                       ,ncol=len(axes['phase'].lines),mode='expand',borderpad=0.2)

  axes['phase'].set_ylabel("Detrended\nPhase (deg)",multialignment='left')
  axes['detection_metric'].set_ylabel('Approximate\nDetection\nMetric',multialignment='left')
  axes['correlation_counts'].set_ylabel('Correlation\nCounts',multialignment='left')
  axes['nco'].set_ylabel('NCO (Hz)',multialignment='center')
  axes['snr'].set_ylabel('C/N$_0$ (dB-Hz)',multialignment='center')

  for subplt in fig.axes: subplt.grid()

  fig.axes[-1].set_xlabel('Time (ms)\nRelative to %s GPS' % eventTime,multialignment='center')


def format_std2Plot(fig,axes,eventTime):
  axes['phase'].legend(bbox_to_anchor=(0., 1.0, 1., .05),loc='lower left'
                       ,ncol=len(axes['phase'].lines),mode='expand',borderpad=0.2)

  ## Label plots
  axes['phase'].set_ylabel("Detrended\nPhase (deg)",multialignment='left')
  axes['detection_metric'].set_ylabel('Approximate\nDetection\nMetric',multialignment='left',fontsize=10)
  axes['correlation_counts'].set_ylabel('Correlation\nCounts',multialignment='left',fontsize=10)
  axes['nco'].set_ylabel('NCO (Hz)',multialignment='center',fontsize=10)
  axes['snr'].set_ylabel('C/N$_0$ (dB-Hz)',multialignment='center')  

  for subplt in fig.axes:
    subplt.grid()
    if subplt != axes['phase']:
      plt.setp(subplt.get_xticklabels(),fontsize=10)
      plt.setp(subplt.get_yticklabels(),fontsize=10)
      subplt.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))

  axes['nco'].xaxis.set_major_locator(ticker.MaxNLocator(prune='both',nbins=14))

  axes['phase'].yaxis.set_minor_locator(ticker.AutoMinorLocator())
  axes['detection_metric'].yaxis.set_minor_locator(ticker.AutoMinorLocator())  
  axes['nco'].xaxis.set_minor_locator(ticker.AutoMinorLocator())

  fig.axes[-1].set_xlabel('Time (ms)\nRelative to %s GPS' % eventTime,multialignment='center')


def format_reportPlot(fig,axes,eventTime):
  axes['phase'].legend(bbox_to_anchor=(0., 1.0, 1., .05),loc='lower left'
                       ,ncol=len(axes['phase'].lines),mode='expand',borderpad=0.2)

  axes['phase'].set_ylabel("Detrended\nPhase (deg)",multialignment='left')
  axes['detection_metric'].set_ylabel('Approximate\nDetection\nMetric',multialignment='left')

  for carrierCode in axes['nco'].keys():
    plt.setp(axes['nco'][carrierCode].get_xticklabels(),fontsize=9)
    plt.setp(axes['nco'][carrierCode].get_yticklabels(),fontsize=9)
    plt.setp(axes['correlation_counts'][carrierCode].get_xticklabels(),fontsize=9)
    plt.setp(axes['correlation_counts'][carrierCode].get_yticklabels(),fontsize=9)

    axes['nco'][carrierCode].set_ylabel(carrierCode,multialignment='center',fontsize=9)
    axes['correlation_counts'][carrierCode].set_ylabel(carrierCode,multialignment='center',fontsize=9)

    if PRINT_SNR_REPORT:
      plt.setp(axes['snr'][carrierCode].get_yticklabels(),fontsize=9)
      plt.setp(axes['snr'][carrierCode].get_xticklabels(),fontsize=9)
      axes['snr'][carrierCode].set_ylabel(carrierCode,multialignment='center',fontsize=9)

  axes['phase'].xaxis.set_major_locator(ticker.MaxNLocator(prune='lower'))
  axes['phase'].xaxis.set_minor_locator(ticker.AutoMinorLocator())

  axes['corrcount_subpane'].set_ylabel('Correlation\nCounts',multialignment='left')
  axes['corrcount_subpane'].yaxis.set_label_coords(-0.07,0.5)
  axes['nco_subpane'].set_ylabel('NCO (Hz)',multialignment='left')
  axes['nco_subpane'].yaxis.set_label_coords(-0.07,0.5)
  if PRINT_SNR_REPORT:
    plt.setp(axes['nco'][carrierCode].get_xticklabels(),fontsize=9)
    axes['SNR_subpane'].set_ylabel('C/N$_0$ (dB-Hz)',multialignment='left')
    axes['SNR_subpane'].yaxis.set_label_coords(-0.07,0.5)

  for axis in fig.axes: axis.grid()

  if not PRINT_SNR_REPORT:
     axis_list = axes['correlation_counts'].values()+axes['nco'].values()
  else:
     axis_list = axes['correlation_counts'].values()+axes['nco'].values()+axes['snr'].values()

  for axis in axis_list: 
    axis.yaxis.set_label_coords(-0.05,0.5);
    axis.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    axis.yaxis.get_major_formatter().set_useOffset(False)
#    axis.yaxis.set_minor_locator(ticker.AutoMinorLocator())

  fig.axes[-1].set_xlabel('Time (ms) Relative to %s GPS'% eventTime,multialignment='center')


def plot_stdPlot(fig,axes,iqTrack,eventTime,obsTrack,iqSlice,obsSlice,detrendedPhase,iqTimes,obsTimes,detectionPeriod,label=None,smoothed=0,plotMovingAvg=False,color=None):
  carrierCode = mdh.getEnumAttributeString(iqTrack,'carrierCode')
  rangingCode = mdh.getEnumAttributeString(iqTrack,'rangingCode')
  codeCarrier = carrierCode+':'+rangingCode
  if color is None: color=cvals[codeCarrier]

  # Make some intermediate calculations for plots
  carrierInfo=signals.findCarrierInfoBy812Code(int(iqTrack.attrs['carrierCode']))
  detectMetric = mtools.computeDetectionMetric(detrendedPhase
                                               ,detectionPeriod
                                               ,carrierInfo.nominalFreqHz
                                               ,mtools.trackrate(iqTrack))
  NCO = - np.diff(obsTrack[obsSlice]['accumulatedDeltaRange'])/np.diff(obsTimes)

  if label is None:
    label= mtools.CodeStrGen(prn=iqTrack.attrs['prnCode']
                             ,carrier_code=carrierCode
                             ,ranging_code=rangingCode)
  # convert time axis to ms
  obsTimes *= 1000.
  iqTimes *= 1000.

  axes['phase'].plot(iqTimes
                     ,sp.signal.medfilt(detrendedPhase,smoothed) if smoothed else detrendedPhase
                     ,alpha=0.8,label=label,color=color)

  # Center detection metric by shifting by half the integration window (3*tau)
  tau_discrete = int(detectionPeriod*mtools.trackrate(iqTrack))
  lostSamples = 3*tau_discrete-1
  slicestart_dm = (3*tau_discrete)/2
  slicestop_dm = slicestart_dm+len(iqTimes)-lostSamples
  axes['detection_metric'].plot(iqTimes[slicestart_dm:slicestop_dm],detectMetric,color=color)

  if plotMovingAvg: # Plot the 'tau' moving average
    slicestart_ma = tau_discrete/2
    slicestop_ma = slicestart_ma+len(iqTimes)-tau_discrete+1
    mov_avg = np.convolve(detrendedPhase,[ 1.0/tau_discrete ]*tau_discrete,mode='valid')
    axes['phase'].plot(iqTimes[slicestart_ma:slicestop_ma],mov_avg
                       ,color='g',linewidth=2.25,alpha=0.8)

  try:              promptlagIdx = np.argmin(iqTrack.attrs['lagOffsetsInSec'])    
  except KeyError:  promptlagIdx = 1 # use center tap

  if isinstance(axes['correlation_counts'],dict):
      axes['correlation_counts'][codeCarrier].plot(iqTimes,iqTrack[iqSlice]['i'][:,promptlagIdx],color=color
                                      ,alpha=0.8,drawstyle='steps-post')
      axes['correlation_counts'][codeCarrier].plot(iqTimes,iqTrack[iqSlice]['q'][:,promptlagIdx],color='k'
                                      ,alpha=0.8,drawstyle='steps-post')
      axes['nco'][codeCarrier].plot(obsTimes[:-1],NCO,color=color)
      if PRINT_SNR_REPORT:
        axes['snr'][codeCarrier].plot(obsTimes,obsTrack[obsSlice]['snr']/100.,color=color)

  else:
      axes['correlation_counts'].plot(iqTimes,iqTrack[iqSlice]['i'][:,promptlagIdx],color=color
                                      ,alpha=0.8,drawstyle='steps-post')
      axes['correlation_counts'].plot(iqTimes,iqTrack[iqSlice]['q'][:,promptlagIdx],color='k'
                                      ,alpha=0.8,drawstyle='steps-post')
      axes['nco'].plot(obsTimes[:-1],NCO,color=color)

      axes['snr'].plot(obsTimes,obsTrack[obsSlice]['snr']/100.,color=color)

  axes['phase'].set_xlim(iqTimes[slicestart_dm],iqTimes[slicestop_dm])


def plot_hb(hblist,start_date=None,end_date=None,title='',axis=None,eps=None):
  '''
  plot_hb

  Plot heartbeat presence time as a bar chart
  '''
  if not axis:
    fig = plt.figure()
    fig.set_size_inches(10.24,7.68)  
    ax = fig.add_subplot(111)
  else:
    fig = axis.figure
    ax = axis

  hblist.sort(key=itemgetter('sv_identifier','carrier_code','ranging_code'))
  hbgroup = groupby(hblist, itemgetter('sv_identifier','carrier_code','ranging_code'))
  svkeystrings = []
  for i, (svkey, group) in enumerate(hbgroup):
    svkeystrings.append('SVN{:02d} {:s}:{:s}'.format(*svkey))

    itvlist = [(dates.date2num(start), dates.date2num(end)-dates.date2num(start))
               for start,end in mtools.mergeIntervals(list(group),'start_date','end_date',eps=eps)]
    startvals,tracklengths = zip(*itvlist)    

    y = [i+1]*len(itvlist)
    ax.barh(y,tracklengths,align='center',left=startvals,color=plt.cm.RdYlBu(float(i)/64),linewidth=0.25)

  plt.setp(ax.get_xticklabels(),fontsize=9)  

  ax.yaxis.set_ticks(np.arange(len(svkeystrings))+1)
  ax.yaxis.set_ticklabels(svkeystrings,fontsize=9)

  x0,x1 = ax.get_xlim()
  ax.set_xlim(start_date if start_date else dates.num2date(x0)
              , end_date if end_date else dates.num2date(x1))
  ax.set_ylim(0,len(svkeystrings)+1)

  ax.grid()

  ax.set_xlabel('Date')
  ax.set_ylabel('SVN Carrier Code:Ranging Code')
  ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
  fig.tight_layout()

  if title: ax.set_title(title,fontsize=14,fontweight='bold')

  return fig


def plot_timeline(eventslist,start_date=None,end_date=None,title='',hblist=None,eps=None):
  hbdict = {}
  if hblist is not None:
    key = ('sv_identifier','carrier_code','ranging_code')
    # hblist.sort(key=itemgetter(*map(hblist.getKeyIndex,key)))
    hblist.sort(key=itemgetter(*key))
    hbgroup = groupby(hblist, key=itemgetter(*key))
    
    for svkey,group in hbgroup:
      svn,cc,rc = svkey
      if rc != 'CA': continue # only use 'CA' for now

      itvlist = [(dates.date2num(start), dates.date2num(end)-dates.date2num(start))
                 for start,end in mtools.mergeIntervals(list(group),'start_date','end_date',eps=eps)]
      hbdict[svn] = itvlist

  # Build up a dict of SVs containing all individual event dates partitioned by day
  svdict = {}
  for r in eventslist:
    eventDay = datetime(r['event_date'].year,r['event_date'].month,r['event_date'].day)

    if r['sv_identifier'] in svdict:
      if eventDay in svdict[r['sv_identifier']]:
        svdict[r['sv_identifier']][eventDay].append(r['event_date'])
      else:
        svdict[r['sv_identifier']][eventDay] = [r['event_date']]
    else:
      svdict[r['sv_identifier']] = {eventDay: [r['event_date']]}

  fig = plt.figure()
  ax = fig.add_subplot(111)
     
  mrkr = 'o'
  for i,sv in enumerate(sorted(svdict)):
    clr = plt.cm.RdYlBu(float(i)/len(svdict))

    for day in svdict[sv]:
      x = svdict[sv][day][:] # get all events for this day
      y = (i+1)*np.ones(len(x)) # align with the ylabel
      ax.plot(x,y,mrkr,color=clr,markersize=4.5,markeredgewidth=0.25,markeredgecolor='0.8')

    if sv in hbdict:
      #hbclr = plt.cm.RdYlBu(float(i+0.5)/len(svdict))
      
      itvlist = hbdict[sv]
      startvals,tracklengths = zip(*itvlist)
      y = [i+1]*len(itvlist)
      ax.barh(y,tracklengths,align='center',left=startvals,linewidth=0
              ,color='0.5',edgecolor='0.5',alpha=0.6)
      
  ax.set_xlabel('Date')
  ax.set_ylabel('SVN')
  fig.subplots_adjust(left=0.08,right=0.95,top=0.925)

  ax.yaxis.set_ticks(np.arange(len(svdict))+1)
  ax.yaxis.set_ticklabels(sorted(svdict),fontsize=10)

  ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
  fig.autofmt_xdate()
  plt.setp(ax.get_xticklabels(),fontsize=10)
  
  ax.set_ylim([0,len(svdict)+1])
  ax.yaxis.grid()

  if title: ax.set_title(title,fontsize=14,fontweight='bold')
  x0,x1 = ax.get_xlim()
  ax.set_xlim(start_date if start_date else x0, end_date if end_date else x1)
  fig.tight_layout()

  return fig
  
def plot_daily_eventrate(eventslist,title='',hblist=None):
  #find earliest and latest dates in hb list and create hb dict
  if hblist is not None:
    hbdict={}
    hbstart = datetime.max
    hbend = datetime.min
    for r in hblist:
      if 'CA' not in r['ranging_code']: continue # only use 'CA' for now
      if r['start_date'] < hbstart: hbstart = r['start_date']
      if r['end_date'] > hbend: hbend = r['end_date']
      
      interval_span = (r['start_date'],r['end_date'])
      if r['sv_identifier'] in hbdict: hbdict[r['sv_identifier']].append(interval_span)
      else: hbdict[r['sv_identifier']] = [interval_span]

  #round to nearest full days
  hbstart = hbstart + timedelta(days=1)
  hbstart = datetime(hbstart.year,hbstart.month,hbstart.day)
  hbend = hbend - timedelta(days=1)
  hbend = datetime(hbend.year,hbend.month,hbend.day)
  
  ##testing
  #hbstart = datetime(2013,11,27)
  #hbend = datetime(2013,12,18)
  ##
  
  #calc number of days in period                                                                                   
  ndays=(hbend-hbstart).days+1
  
  #build up a dict of SVs containing number of events per day
  svdict = {}
  for r in eventslist:
    eventDay = datetime(r['event_date'].year,r['event_date'].month,r['event_date'].day)
    if not (eventDay >= hbstart and eventDay <= hbend):
       continue
    if r['sv_identifier'] in svdict:
      svdict[r['sv_identifier']][eventDay]+=1
    else:
      svdict[r['sv_identifier']]={}
      for d in range(ndays):  #create a dict with all days in period
        currd = hbstart+timedelta(days=d)
        for h in hbdict[r['sv_identifier']]: #that overlap with active heartbeats during period
           currhb0 = datetime(h[0].year,h[0].month,h[0].day)
           currhb1 = datetime(h[1].year,h[1].month,h[1].day)
           if currd==currhb0 or currd==currhb1: #simplistic check for hb overlap
             svdict[r['sv_identifier']][currd]=0 
             break
      svdict[r['sv_identifier']][eventDay] = 1

  fig = plt.figure()
  ax = fig.add_subplot(111)
  
  ##testing
  #svdict_temp = svdict
  #svdict={}
  #svdict[63]=svdict_temp[63]
  #for s in svdict:
  #   for d in sorted(svdict[s]):
  #      print d,svdict[s][d]
  ##
  
  mrkr = '.'
  for i,sv in enumerate(sorted(svdict)):
     clr = plt.cm.RdYlBu(float(i)/len(svdict))
     d = svdict[sv].keys()
     v = svdict[sv].values()
     ax.plot(d,v,'.',color=clr)
    
  ax.set_xlabel('Date')
  ax.set_ylabel('Daily event rate')
  
  svlist=[]
  for sv in sorted(svdict):
    svlist.append(str(sv))
  plt.legend(svlist,'upper right')
## would prefer below for legend, but it gets cutoff unless using something like the below savefig command  
  #leg = ax.legend(svlist,bbox_to_anchor=(0., 1.0, 1., .05),loc='lower left'
  #         ,ncol=4,mode='expand',borderpad=0.2)
  
  ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
  fig.autofmt_xdate()
  plt.setp(ax.get_xticklabels(),fontsize=10)

  fig.tight_layout()

  if title: ax.set_title(title,fontsize=14,fontweight='bold')

## below will generate plot correctly - that is, not cutting off legend when placed above plot
## would need to return legend object to genericPlot for this to work though
#  fig.savefig('/v/sgggid/SMWG/Documents/SVAnomalyReports/2014/Q1/Q1_events_filt_new_edit/test/test.png', bbox_extra_artists=(leg,), bbox_inches='tight')
  
  return fig


'''
plot_eventcount

Plots a histogram of event counts on each SV based on constraints in 'form'.

Note:
  One event is defined as the case where at least one ranging code on an SV triggered an anomaly detection at a given point in time (i.e. we do not consider which ranging codes triggered an event).
'''
def plot_eventcount(eventslist=None,title='',svkey='sv_identifier'):
  # Cluster overlapping events occurring on the same SV
  svEventDates = {}

  itvdict = {}
  for r in eventslist:
    espan = (r['event_date']-timedelta(seconds=r['detection_period'])
             ,r['event_date']+timedelta(seconds=r['detection_period']))
    if r[svkey] in itvdict:
      itvdict[r[svkey]].append(espan)
    else:
      itvdict[r[svkey]] = [espan]
  
  for sv,itvlist in itvdict.iteritems():
    svEventDates[sv] = mtools.mergeIntervals(itvlist,startkey=0,stopkey=1)
  
  fig = plt.figure()

  ax = fig.add_subplot(111)

  x = np.arange(len(svEventDates),dtype='int')+1 # Create range for SVs
  y = [len(v) for k,v in sorted(svEventDates.items())] # Grab event count for each SV in sorted order
  ax.bar(x,y,align='center')

  # Label the range with SVNs in sorted order
  ax.xaxis.set_ticks(x)
  ax.xaxis.set_ticklabels(sorted(svEventDates),fontsize=9,rotation='vertical')
  ax.yaxis.set_major_locator(ticker.MaxNLocator(9))
  ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

  ax.set_ylabel('Number of Detections')
  ax.set_xlabel('SVN' if svkey == 'sv_identifier' else svkey.upper())
  ax.set_xlim(0,len(svEventDates)+1)

  if title: ax.set_title(title,fontsize=14,fontweight='bold')
  fig.tight_layout()

  return fig


'''
plot_hbcum

Plot cumulative heartbeat presence time as a bar chart
'''
def plot_hbcum(hblist,title=''):
  svdict = {}
  for r in hblist:
    codecarrier = r['carrier_code']+r['ranging_code']
    svn = r['sv_identifier']
    if svn not in svdict: svdict[svn] = {}
    if codecarrier not in svdict[svn] : svdict[svn][codecarrier] = timedelta()
    svdict[svn][codecarrier] += r['end_date']-r['start_date']

  # keep only the max tracking time for each SV 
  for sv in svdict:
    svdict[sv] = max(svdict[sv][cc] for cc in svdict[sv])

  fig = plt.figure()

  ax = fig.add_subplot(111)

  x = np.arange(len(svdict))+1
  y = [svdict[sv].total_seconds()/3600. for sv in sorted(svdict)]
  ax.bar(x,y,align='center')
  
  ax.xaxis.set_ticks(x)
  ax.xaxis.set_ticklabels(sorted(svdict),fontsize=9)
  ax.yaxis.set_major_locator(ticker.MaxNLocator(9))
  ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

  ax.set_ylabel('Cumulative tracking time (Hours)')
  ax.set_xlabel('SVN')
  ax.set_xlim(0,len(svdict)+1)

  if title: ax.set_title(title,fontsize=14,fontweight='bold')
  fig.tight_layout()

  return fig


def GetCommonData(iqTrack,obsTrack,eventTime,detectionPeriod,startTime=None,stopTime=None):
  if startTime is None: startTime = mtools.tracktimefromslice(iqTrack,0)
  if stopTime is None: stopTime = mtools.tracktimefromslice(iqTrack,len(iqTrack))

  # Note: detection metric is shorter than original signal length by 3*detection_period - 1 samples
  tauDiscrete = int(detectionPeriod*mtools.trackrate(iqTrack))
  lostSamples = 3*tauDiscrete-1

  # Expand the selection of data to compensate for the samples lost by calculating the detection metric
  # round down
  iqTrackStart = long(mtools.slicefromtracktime(iqTrack,startTime) - (lostSamples+1)/2 ) 
  # round up (+2 to avoid using ceil)
  iqTrackStop = long(mtools.slicefromtracktime(iqTrack,stopTime) + (lostSamples+2)/2 )

  # truncate to ends of track
  if iqTrackStart < 0: iqTrackStart = 0 
  if iqTrackStop > len(iqTrack): iqTrackStop = len(iqTrack)

  # Get the overlapping data times
  trackData = mtools.GetDetrendedPhase(iqTrack,obsTrack,(iqTrackStart,iqTrackStop))
  
  # make sure there is enough data to plot detection metric
  if lostSamples >= (1.5*tauDiscrete + len(trackData['iqTimes'])):
    if VERBOSE: print "Error: Insufficient data for plotting detection metric."
    return None
  
  # set times relative to event time
  trackData['iqTimes'] -= eventTime
  trackData['obsTimes'] -= eventTime

  return trackData


# dict of function pointers for different plot styles
Plot_Style = {'report':(setup_reportPlot,plot_stdPlot,format_reportPlot)
              ,'std':(setup_stdPlot,plot_stdPlot,format_stdPlot)
              ,'std2':(setup_std2Plot,plot_stdPlot,format_std2Plot)}

reqfields = ['file_path','prn','carrier_code','ranging_code','event_date','detection_period']  

def plot_event(eventdict,savedir='',plotall=False,style='std',tlim=None,smoothed=0,alternate_codes=None,title=None,filetag='',plotreference=False):
  '''
  plot_event

  Process generic dictionary like object to generate appropriate inputs for
  event plotting.

  By default just plots the event parameterized by the 'eventdict' input Use the
  alternate_codes input to custom specify signals to plot (the event mdh file
  will be searched for any existing tracks).
  '''
  if not all([isvalid(field,eventdict) for field in reqfields]):
    if VERBOSE:
      missing_fields = ','.join([k for k in reqfields if not isvalid(k,eventdict)])
      print "Error: Missing required field(s): %s" % missing_fields
    return (None,None)

  m=mdh.MDHFile(eventdict['file_path'],'r')
  eventTime = mtools.timegps(eventdict['event_date'])

  # Attempt to load a valid IQ track
  if not isvalid('data_set_iq',eventdict) or eventdict['data_set_iq'] not in m:
    iqTrack = mtools.FindBestTrack(m,eventTime,eventdict['prn'],eventdict['carrier_code'],eventdict['ranging_code'],SUBCLASS="IQ_SUMS")
  else:
    iqTrack = m[eventdict['data_set_iq']]

  obsTrack = mtools.FindBestTrack(m,eventTime,eventdict['prn'],eventdict['carrier_code'],eventdict['ranging_code'],SUBCLASS='OBSERVATIONS')

  # Compile a list of the signals to be plotted along with the pertinent tracks
  signal_list = []
  track_list = []
  if plotall: # (this option is to be removed) plot all tracks in this dataset
    alternate_codes = mtools.matchingSignals(m,SUBCLASSES=['OBSERVATIONS'])
  elif style == 'report' and not alternate_codes:
    alternate_codes = ('L1','CA'),('L2','CM'),('L5','I5')
  elif not alternate_codes:
    alternate_codes = [itemgetter('prn','carrier_code','ranging_code')(eventdict)]

  for codecarrier in alternate_codes:
    cc,rc = codecarrier[-2:]
    prn = eventdict['prn'] if (len(codecarrier) == 2) else codecarrier[0]

    if prn == eventdict['prn'] and cc == eventdict['carrier_code'] and rc == eventdict['ranging_code']:
      iq = iqTrack
      obs = obsTrack
    else:
      iq = mtools.FindBestTrack(m,eventTime,prn=prn,carrier_code=cc,ranging_code=rc)
      obs = mtools.FindBestTrack(m,eventTime,associatedDataset=iq,SUBCLASS='OBSERVATIONS')
    track_list.append((iq,obs))
    signal_list.append((prn,cc,rc))

  if plotreference:
    # Limit reference signals to these for the report style plot since report
    # style currently doesn't support plotting anything outside these signals
    refcc = ('L1:CA','L2:CM','L5:I5') if style == 'report' else None
    ref_iq,ref_obs = mtools.FindBestRefTracks(m,eventTime,eventdict['detection_period']
                                              ,eventdict['prn'],codeCarriers=refcc)

    signal_list.append((ref_iq.attrs['prnCode']
                        ,mdh.getEnumAttributeString(ref_iq,'carrierCode')
                        ,mdh.getEnumAttributeString(ref_iq,'rangingCode')))
    track_list.append((ref_iq,ref_obs))

  # setup figure and axes
  setupHandler,plotHandler,formatHandler = Plot_Style[style]
  fig = plt.figure()
  axes = setupHandler(fig)

  # Generate plots for all signals that are relevant to this event
  minSNR,medSNR,maxSNR = [float('nan')]*3
  for sig,track in zip(signal_list,track_list):
    prn,cc,rc = sig
    iq,obs = track

    if not iq or not obs:
      if VERBOSE: "WARNING: No data for %s. Skipping..." % str((prn,cc,rc))
      continue

    label = mtools.CodeStrGen(prn=prn,carrier_code=cc,ranging_code=rc)
    color = cvals[cc+':'+rc]
    if prn != eventdict['prn']:
      label = '(Reference SV) ' + label
      color = 'm'
    elif rc == eventdict['ranging_code']:
      label = '*' + label

    startTime, stopTime = map(mtools.timegps,tlim) if tlim else (None,None)
    commonData = GetCommonData(iq,obs,eventTime,eventdict['detection_period']
                               ,startTime=startTime,stopTime=stopTime)
    if not commonData: continue

    if prn == eventdict['prn'] and rc == eventdict['ranging_code']:
      obsSNR = obs['snr'][commonData['obsSlice']]
      minSNR = min(obsSNR/100.)
      medSNR = np.median(obsSNR)/100.
      maxSNR = max(obsSNR)/100.

    plotHandler(fig=fig
                ,axes=axes
                ,iqTrack=iq
                ,obsTrack=obs
                ,eventTime=eventTime
                ,detectionPeriod=eventdict['detection_period']
                ,label=label
                ,smoothed=smoothed
                ,plotMovingAvg=bool(prn == int(eventdict['prn']) and rc == eventdict['ranging_code'])
                ,color=color
                ,**commonData)

  m.close() # done using the mdhfile

  # apply some final formatting
  formatHandler(fig,axes,eventdict['event_date'])
  if title is None: title=("Event ID: %d" % eventdict['event_id']) if 'event_id' in eventdict else ''
  fig.suptitle(title,fontsize=14,fontweight='bold')

  # Annotate with detection metrics if this was a detected event (i.e. not a reference)
  if (eventdict['prn'],eventdict['carrier_code'],eventdict['ranging_code']) in signal_list:
    dm_annot = ("Detection_Period = %s, Threshold: %s"
               %( ("%d ms" % int(eventdict['detection_period']*1000.))
                  ,("%.2f" % eventdict['threshold'])
                  if isvalid('threshold',eventdict) else '-'))
  else:
    dm_annot = "Detection_Period = - ms, Threshold: -"
        
  axes['detection_metric'].annotate(dm_annot,xy=(0.01,0.825)
                                    ,xycoords='axes fraction'
                                    ,fontsize=10,color='r')
  if style != 'report':
    axes['phase'].annotate("%d-point median filter applied" % (smoothed) if smoothed else "No filtering applied"
                         ,color='r',xy=(0.01,0.825),xycoords='axes fraction',fontsize=10)
  if not PRINT_SNR_REPORT:
    axes['phase'].annotate("SNR = (%.2f Min,%.2f Med, %.2f Max) dB-Hz" % (minSNR,medSNR,maxSNR)
                         ,xy=(0.69,0.825),xycoords='axes fraction',color='r',fontsize=10)

  if savedir:
    if not os.path.exists(savedir):
      if VERBOSE: print "Error: savepath does not exist."
      plt.close()
      del(fig)
      return (None,None)
    if not filetag:
      tags = tuple(k for k in ('event_date','prn','carrier_code','ranging_code','event_id') if k in eventdict)
      filetag= FileTagGen(**dict(zip(tags,map(eventdict.get,tags))))
    plt.savefig(os.path.join(savedir,filetag+'.png'), dpi=100)
    
  return fig,axes


if __name__ == "__main__":
# ----------------------------------------------------------
# Define the CLI
# ----------------------------------------------------------    
  from optparse import OptionParser

  parser =OptionParser(usage="./mdh_plotevents [options] [--eventslog /path/to/events.log] or [--mdhfile=/path/to/mdhfile --prn=PRN --time TIME]")

  parser.add_option("--eventslog", dest='eventslog', type='string', action='append', default=None,
                    help="Path to text file listing events. (None)")

  parser.add_option("--trackslog", dest='trackslog', type='string', action='append', default=None,
                    help="Path to text file listing heartbeats. (None)")

  parser.add_option("--mdhfile", dest='mdhfile', type='string', default=None,
                    help="Path to MDH file for plotting. (None)")

  parser.add_option("--savepath", dest='savepath',type='string', default=None,
                    help="Save plots directly to the file path specified.")

  parser.add_option("-p", "--prn", dest="prn", metavar='PRN', default=None,
                    help="Specify SVs to process. (All SVs)")

  parser.add_option("--cc", "--code-carriers", dest="codeCarriers", default=None,
                    help=("Plot data for carrier:code combination(s). "
                          +"Pass multiple signals in as a comma delimited string. "
                          +"(default is code carrier that event occurred on or L1:CA if not known)."))

  parser.add_option("--exSats", dest="exsats", default=None,
                    help="Exclude satellites from analysis. Pass multiple PRNs as a comma delimited string. (None)")

  parser.add_option("--plot-reference",dest='plotreference',action='store_true',default=False
                    ,help="For each plot, overlay a plot of the highest SNR sattelite at event time. (False)")

  parser.add_option("--selectevents", dest="selectevents", action='store_true',default=False,
                    help='Use only events that are *\'d. (False)')

  parser.add_option("-t", "--time", dest="time", default=None, metavar="TIME",
                    help="Center plot at TIME. (Middle of file)")

  parser.add_option("-w", "--window", dest="window", default=10., metavar="SECs", type='float',
                    help="Force plot time to a window of 'SECs' seconds around the event time.")

  parser.add_option("--tlim", dest='tlim', default=None, type='string',
                    help="Plot data on the interval [t1,t2] where t1,t2 are passed in as a tuple t1,t2. (None)")

  parser.add_option("--event_id",dest="eventID",default=None,
                    help="Specify a single event number to plot. (None)")

  parser.add_option("--folders","-f",dest="folders",action='store_true',default=False,
                    help="When saving multiple plots, creates directories to organize plots by PRN. (False)")

  parser.add_option("--tau", dest='tau', type='float', default=0.128,
                    help="Specify integration period (in seconds) for computing detection metric (0.128).")

  parser.add_option("-s", "--smooth", dest="smoothed", type='int', default=None,
                    help="Apply nth order median filter to carrier phase. (None)")

  parser.add_option("--timefmt", dest="timefmt", default="%Y-%m-%d %H:%M:%S",
                    help="Set time format string. (%Y-%m-%d %H:%M:%S)")

  parser.add_option("-v", dest="verbose", action='count',default=1,
                    help="Print status information to command line. Enter multiple 'v's for more verbosity. (v)")

  parser.add_option("-q","--quiet",dest='quiet',action='store_true',default=False
                  ,help="Disable verbose output. (False)")

  addPlotOptions(parser)

# ----------------------------------------------------------
# Parse multi argument options and set defaults
# ----------------------------------------------------------    
  opts,args = parser.parse_args()
  VERBOSE = 0 if opts.quiet else opts.verbose
  mtools.VERBOSE = VERBOSE
    
  if opts.savepath and not (os.path.exists(opts.savepath) and os.path.isdir(opts.savepath)):
    sys.exit("Error: savepath does not exist or is not a directory.")

  if not opts.codeCarriers and opts.style == 'report':
    opts.codeCarriers = "L1:CA,L2:CM,L5:I5"

  opts.window = timedelta(seconds=opts.window)
  if opts.tlim:
    t1,t2 = opts.tlim.split(',')
    opts.tlim = (datetime.strptime(t1,opts.timefmt),datetime.strptime(t2,opts.timefmt))
    opts.time = opts.tlim[0]+timedelta(seconds=(opts.tlim[1]-opts.tlim[0]).total_seconds()/2.)
  elif opts.time:
    opts.time = datetime.strptime(opts.time,opts.timefmt)    
    opts.tlim = (opts.time-opts.window/2,opts.time+opts.window/2)

  # parse multi-arg options
  if opts.exsats:   opts.exsats = map(int,opts.exsats.split(","))
  elif opts.prn:    opts.prn = map(int,opts.prn.split(","))
  if opts.eventID:  opts.eventID = map(int,opts.eventID.split(","))

  # form a filter for the log file (much like a where clause for a query)
  form = {}
  if opts.exsats: form['prn'] = ('not in',opts.exsats)
  elif opts.prn:  form['prn'] = ('in', opts.prn)

  if opts.codeCarriers:
    cc,rc = zip(*[s.split(':') for s in opts.codeCarriers.split(',')])
    form['carrier_code'] = ('in', cc)
    form['ranging_code'] = ('in', rc)

  events = opts.eventslog
  if opts.eventslog:
    if opts.eventID:
      eform = form.copy()
      eform['event_id'] = ('==',opts.eventID)
    else:
      eform = form

    keys, events, start_events, end_events = mtools.mergelogs(opts.eventslog
                                                              ,selectlines=opts.selectevents
                                                              ,form=eform)
    # events = mtools.DictList(keys,events)
    events = [dict(zip(keys,vals)) for vals in events]
  tracks = opts.trackslog
  if opts.trackslog:
    keys, tracks, start_hb, end_hb = mtools.mergelogs(opts.trackslog,form=form)
    # tracks = mtools.DictList(keys,tracks)
    tracks = [dict(zip(keys,vals)) for vals in tracks]

# ----------------------------------------------------------
# option: --plotevents
# Plots events from an event log file.
# ----------------------------------------------------------    
  if opts.plotevents:
    for row in events:
      if VERBOSE: print "\rProcessing EventID: %(event_id)d..." % row

      # Create a directory for this prn if necessary
      if opts.savepath and opts.folders:
          savedir = os.path.join(opts.savepath,mtools.CodeIDGen(prn=prn))
          if not os.path.exists(savedir): os.mkdir(savedir)
      else:
          savedir = opts.savepath

      # List the attributes of this event
      if VERBOSE > 1:
        for k,v in row.items(): print "\t",k,':',v

      if opts.window:
         opts.tlim = (row['event_date']-opts.window/2,row['event_date']+opts.window/2)

      fig,axes = plot_event(row, savedir=savedir, style=opts.style
                            , tlim=opts.tlim, smoothed=opts.smoothed
                            , plotreference=opts.plotreference)

      if fig is None and VERBOSE: print 'Plotting failed.'
      if fig is not None and not opts.savepath: plt.show(fig.number)

# ----------------------------------------------------------
# Plot any aggregate events or heartbeat data plots
# ----------------------------------------------------------
  # eps = timedelta(hours=5) # merge heartbeats that are less than 5 hours apart
  eps = None
  
  if opts.eventslog:
    if opts.ploteventcount:
      genericPlot(plot_eventcount,events,opts.savepath,"Event count"
                  ,start_date=start_events,end_date=end_events)
    if opts.plottimeline:
      genericPlot(plot_timeline,events,opts.savepath,'Event timeline'
                  ,timeseries=True
                  ,start_date=start_events,end_date=end_events,eps=eps)
  if opts.trackslog:
    if opts.plothb:
        genericPlot(plot_hb,tracks,opts.savepath,'Heartbeat presence'
                    ,timeseries=True
                    ,start_date=start_hb,end_date=end_hb,eps=eps)
    if opts.plothbcum:
        genericPlot(plot_hbcum,tracks,opts.savepath,"Cumulative heartbeat presence"
                    ,start_date=start_hb,end_date=end_hb)
  if opts.eventslog and opts.trackslog:
      if opts.plottimelineplushb:
          genericPlot(plot_timeline,events,opts.savepath,'Event timeline and heartbeat'
                      ,timeseries=True
                      ,start_date=min(start_hb,start_events)
                      ,end_date=max(end_hb,end_events)
                      ,hblist=tracks,eps=eps)
      if opts.ploteventrates:
          genericPlot(plot_daily_eventrate,events,opts.savepath,'Daily event rate'
                      ,start_date=start_events,end_date=end_events
                      ,hblist=tracks)

# ----------------------------------------------------------
# option: --mdhfile
# Plots an MDH file at the specified time
# ----------------------------------------------------------    
  elif opts.mdhfile:
    if not opts.time: sys.exit('Need event time')
    if opts.codeCarriers:
      ac = [s.split(':') for s in opts.codeCarriers.split(',')]
    else:
      ac = [('L1', 'CA')]
    cc,rc = ac[0]
    eventdict = {'file_path': os.path.abspath(opts.mdhfile)
                 , 'prn': opts.prn[0]
                 , 'carrier_code': cc
                 , 'ranging_code': rc
                 , 'event_date': opts.time
                 , 'detection_period': opts.tau}

    fig,axes = plot_event(eventdict,savedir=opts.savepath,style=opts.style
                          ,tlim=opts.tlim,smoothed=opts.smoothed
                          ,plotreference=opts.plotreference,alternate_codes=ac)
    
    if fig is None and VERBOSE: print '\tPlotting failed.\n'
    if fig is not None:
        if not opts.savepath:
            plt.show(fig.number)
        else:  
            tags = ('event_date','prn','carrier_code','ranging_code')
            filetag= FileTagGen(**dict(zip(tags,map(eventdict.get,tags))))
            plt.savefig(os.path.join(savedir,filetag+'.png'), dpi=100)

