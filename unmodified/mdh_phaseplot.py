#!/usr/bin/python
import sys, os

import mdh_tools as mtools
from mdh_tools import timegps,gpstime
import mdh
import h5py
import gpstk

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as dates
from scipy.interpolate import LSQUnivariateSpline
import numpy as np

VERBOSE = 0


def showPlot(fig,save=None):
    for j,ax in enumerate(fig.axes):
        ax.grid()
        plt.setp(ax.get_yticklabels(),fontsize=9)
        if (j+1) < len(fig.axes):
            plt.setp(ax.get_xticklabels(),visible=False)
        else:
            plt.setp(ax.get_xticklabels(),fontsize=9)
            plt.setp(ax.get_xticklabels(minor=True),fontsize=9,visible=True)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.15)
    if save:
        fig.savefig(save,dpi=100)
    else:
        plt.show(fig.number)

        
def mdhSpan(m,subclass):
    dsstart, dsstop = zip(*[(ds.attrs['startTime']/ds.attrs['timeDenominator']
                               ,mtools.getStopTime(ds)/ds.attrs['timeDenominator']) for ds in m.matchingDatasets(SUBCLASS=subclass)])
    return min(dsstart), max(dsstop)
        

def GetDetrendedPhase(mdhfile,startTime=None,stopTime=None,prnCodes=None,codeCarriers=None,cadence=None):
    if startTime: startTime = mtools.timegps(startTime)
    if stopTime: stopTime = mtools.timegps(stopTime)

    if prnCodes and not hasattr(prnCodes,'__iter__'): prnCodes = [prnCodes]

    if codeCarriers: codeCarriers = map(tuple,codeCarriers)
    for iq in mdhfile.matchingDatasets(SUBCLASS="IQ_SUMS"):
        if codeCarriers and (mdh.getEnumAttributeString(iq,'carrierCode'),mdh.getEnumAttributeString(iq,'rangingCode')) not in codeCarriers: continue
        if prnCodes and iq.attrs['prnCode'] not in prnCodes: continue

        iqStart = max(startTime*iq.attrs['timeDenominator'],iq.attrs['startTime']) if startTime else iq.attrs['startTime']
        iqStop = min(stopTime*iq.attrs['timeDenominator'], mtools.getStopTime(iq)) if stopTime else mtools.getStopTime(iq)
        if iqStop <= iqStart: continue

        for obs in mdhfile.matchingDatasets(SUBCLASS="OBSERVATIONS"
                                            ,prnCode = iq.attrs['prnCode']
                                            ,carrierCode = iq.attrs['carrierCode']
                                            ,rangingCode = iq.attrs['rangingCode']):
            obsStart = obs.attrs['startTime']/obs.attrs['timeDenominator']-obs['pseudorange'][0]/gpstk.C_MPS
            obsStop = mtools.getStopTime(obs)/obs.attrs['timeDenominator']-obs['pseudorange'][-1]/gpstk.C_MPS
            iqObsStart = max(iqStart, obsStart*iq.attrs['timeDenominator'])
            iqObsStop = min(iqStop, obsStop*iq.attrs['timeDenominator'])
            if iqObsStart < iqObsStop:
                sliceStart = (iqObsStart-iq.attrs['startTime'])/iq.attrs['cadence']
                n = (iqObsStop-iqObsStart)/iq.attrs['cadence']
                dphase = mtools.GetDetrendedPhase(iq,obs,(int(sliceStart),int(sliceStart+n)),cadence=cadence)
                dphase.update(dict(iq.attrs))
                yield dphase


from argparse import ArgumentParser

parser = ArgumentParser(description="Tool for viewing detrended phase from an MDP in HDF5 file")

parser.add_argument("mdh_files",nargs='+')

parser.add_argument('--cadence', dest="cadence", default=10, type=int
                    ,help="Number of seconds to fit over for carrier phase trend fitting. (10)")

parser.add_argument("--prn", dest="prn", metavar='PRN', default=None,
                  help="Specify SVs to process. (All SVs)")

parser.add_argument("--cc", "--code-carriers", dest="codeCarriers", default='L1:CA',
                  help=("Plot data for carrier:code combination(s). "
                        +"Pass multiple signals in as a comma delimited string. (%(default)s)"))

parser.add_argument("--start", dest="start", default=None
                  , help="Consider only data after this time. (min{All OBS dataset's start time})")

parser.add_argument("--stop", dest="stop", default=None
                  , help="Consider only data before this time. (max{All OBS dataset's end time})")

parser.add_argument("--timefmt", dest="timefmt", default="%Y-%m-%d %H:%M:%S.%f"
                    ,help="Set time format string. (%(default)s)")

parser.add_argument("--savepath", dest='savepath',default=None
                    ,help="Save plots directly to the specified directory.")

opts = parser.parse_args()

if opts.codeCarriers: opts.codeCarriers = [code.split(':') for code in opts.codeCarriers.split(',')]
# if opts.exsats:   opts.exsats = map(int,opts.exsats.split(","))
if opts.prn:    opts.prn = map(int,opts.prn.split(","))

if opts.start: opts.start = mtools.parseDate(opts.start,opts.timefmt)
if opts.stop: opts.stop = mtools.parseDate(opts.stop,opts.timefmt)

filenames = filter(h5py.is_hdf5, opts.mdh_files)
if not filenames: sys.exit(-1)
m = mdh.MDHFileList(*filenames,mode='r')

minStartTime, maxStopTime = map(mtools.gpstime,mdhSpan(m,"OBSERVATIONS"))
if opts.start: minStartTime = max(opts.start,minStartTime)
if opts.stop: maxStopTime = min(opts.stop,maxStopTime)
if opts.cadence and (maxStopTime-minStartTime).total_seconds() < 15: opts.cadence = None

sigTuples = sorted(mtools.matchingSignals(m,minStartTime,maxStopTime,SUBCLASSES=["OBSERVATIONS"],codeCarriers=opts.codeCarriers,prnCodes=opts.prn))
if not sigTuples: sys.exit(0)
    
print "# Analysis time span is [%s,%s)" % (minStartTime,maxStopTime)
print "# Found signals:", sigTuples
fig = plt.figure()
fig.set_size_inches(14.40,8.00)

ax = fig.add_subplot(len(sigTuples),1,1)
for i,(prn,cc,rc) in enumerate(sigTuples):
    if i > 0: ax = fig.add_subplot(len(sigTuples),1,i+1,sharex=ax)

    for dphase in GetDetrendedPhase(m,minStartTime,maxStopTime,cadence=opts.cadence,prnCodes=[prn],codeCarriers=[(cc,rc)]):
        if not dphase: continue
        label = 'PRN%02d %s:%s' % (prn,cc,rc)
        iqTimes = dphase['iqTimes']/86400. + mdh.GPS_EPOCH_IN_MATPLOTLIB_DATE_FORMAT
        ax.plot(iqTimes,dphase['detrendedPhase'],alpha=0.8)
        ax.set_ylabel(label,fontsize=9)
fig.autofmt_xdate()
ax.set_xlim(minStartTime,maxStopTime)
# ax.set_xlim(mtools.timegps(minStartTime),mtools.timegps(maxStopTime))
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.xaxis.set_major_formatter(dates.DateFormatter("%H:%M:%S.%f"))

tag = "{:%Y%m%d-%H%M%S}-{:%Y%m%d-%H%M%S}_dphase.png".format(minStartTime,maxStopTime)
showPlot(fig,os.path.join(opts.savepath,tag) if opts.savepath else None)
