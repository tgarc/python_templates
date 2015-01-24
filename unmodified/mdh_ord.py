#!/usr/bin/python

import os,sys
from sys import stderr, exit
import glob
from collections import OrderedDict
from operator import itemgetter
import traceback
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.dates as dates
import time

import gpstk
import mdh,h5py
import numpy as np
from scipy.interpolate import LSQUnivariateSpline
import mdh_tools as mtools
from mdh_tools import timegps,gpstime,TextTable,LatexTable,Report

DEBUG = 0
VERBOSE = 0


gps_to_tai = 19;  # Offset, in seconds from GPS to TAI

def york2dt(t):
    global gps_to_tai
    t_gps = (t - gps_to_tai)/86400.0
    return datetime(1980, 1, 6).toordinal() + t_gps

def york2pt(t):
    global gps_to_tai
    t_gps = t - gps_to_tai
    gps_t0 = datetime.datetime(1980, 1, 6)
    return datetime(1980, 1, 6).toordinal() + datetime.timedelta(seconds=t_gps)

def MAD(a, c=0.6745):
    """Estimate the variance of a sequence by computing the median
    absolute deviation"""
    d = np.median(a)
    return np.median(abs(a - d))/c


def writeExpObs(startTime, stopTime, expobs_fh, ephReader, obs1, obs2, rxpos):
    mdhTimes = xrange(startTime, stopTime, obs1.attrs['cadence'])
    times = mtools.mdhtime2commontime(mdhTimes, obs1.attrs['timeDenominator'])

    cc1 = mdh.getEnumAttributeString(obs1,'carrierCode')+mdh.getEnumAttributeString(obs1,'rangingCode')
    cc2 = mdh.getEnumAttributeString(obs2,'carrierCode')+mdh.getEnumAttributeString(obs2,'rangingCode')        
    sv = gpstk.SatID(int(obs1.attrs['prnCode']),gpstk.SatID.systemGPS)

    tags = [os.path.basename(obs1.name),cc1.replace(':',''),cc2.replace(':',''),"EXPOBS"]
    dsname = "_".join(tags)
    if dsname in expobs_fh: del expobs_fh[dsname]

    try:
        data = ephReader.getExpObs(rxpos,sv,times)
    except KeyboardInterrupt:
        del data
        return None

    try:
        dsetOut = expobs_fh.create_dataset(dsname,data=data,maxshape=(None,),chunks=True)

        for attr,val in obs1.attrs.iteritems():
            dt = h5py.h5a.open(obs1.id,attr).dtype
            dsetOut.attrs.create(attr,val,dtype=dt)
        dsetOut.attrs['SUBCLASS'] = "EXPOBS"
        dsetOut.attrs['startTime'] = startTime
    except IOError:
        if DEBUG: print traceback.print_exc()
        del data, dsetOut
        return None

    return dsetOut


class ORDReport(Report):
    _keys = ("bin","n","mad","std","mean","max","min")
    _fmt = {k:'{:.4f}' for k in _keys if k not in ("bin","n")}
    _fmt['bin'] = "{0[0]}-{0[1]}"

    def generateContent(self,bins,b_data):
        for b in bins:
            self.addRow(b,len(b_data[b]),MAD(b_data[b]),np.std(b_data[b])
                         ,np.mean(b_data[b]),np.max(b_data[b]),np.min(b_data[b]))
        return self


def AverageRanges(ords,cadence,timeDenominator):
    '''
    aligns and averages ords
    NOTE: assumes same cadence/timeDenominator across datasets
    '''
    ords.sort(key=itemgetter('startTime','stopTime'))
    minstart = ords[0]['startTime']
    maxstop = max(map(itemgetter('stopTime'),ords))

    # create a timestamp array that extends over all possible times
    # in order to time align the different ords data    
    sampsize = (maxstop - minstart)/cadence
    merged_times = minstart + np.arange(sampsize)*cadence
    samplesAtEpoch = np.zeros(sampsize)
    avgrange = np.zeros(sampsize)

    # fill in the arrays with the available ords data
    for row in ords:
        # find the overlapping indices including those times where the data is masked out
        startidx = (row['startTime'] - minstart)/cadence
        stopidx = startidx + (row['stopTime']-row['startTime'])/cadence

        # only write to the masked data times
        samplesAtEpoch[startidx:stopidx][row['ordsTimeMask']] += 1
        avgrange[startidx:stopidx][row['ordsTimeMask']] += row['ordsData']

    sampledEpochs = (samplesAtEpoch != 0)
    avgrange[sampledEpochs] /= samplesAtEpoch[sampledEpochs]
        
    return merged_times, avgrange, sampledEpochs


# ============================================================
# Add options for command line interface
# ============================================================
import argparse
def position(value):
    return gpstk.Position(*map(float, value.split(",")))
def codeCarrierCombos(value):
    return tuple(value.split(','))

parser = argparse.ArgumentParser(usage="%prog --position <x,y,z> --obs <obsfile> --eph <ephfile> [--expobs <expobsfile>]"
                               ,description="Plots/summarizes iono-free ORDs from mdh data.")

parser.add_argument('--obs', dest="obs", metavar='mdh-fn', required=True, action='append'
                    , help="MDH file with obs data. (required)")

parser.add_argument('--expobs', dest="expobs", default='/tmp/expobs_%d.mdh' % int(time.time())
                    , help=("Path to MDH file with (partial) expected obs data."
                            +" If file does not exist an expected obs file will"
                            +" be written to the specified path. (%(default)s)"))

parser.add_argument("--read-only",dest='read_only',action='store_true'
                    , help="Set read only mode for expected obs file. (%(default)s)")

parser.add_argument("-e", "--eph", dest='ephemeris', required=True, action='append'
                    , help=("Path to ephemeris file." 
                            +" Accepts multiple ephemeris files when --eph is"
                            +" specified multiple times."))

parser.add_argument("-p","--position", dest='rxpos', type=position, required=True
                    , help="Receiver X,Y,Z position.")

parser.add_argument("--plot-ords",dest='plot_ords',default=True, nargs='?', const=None
                    , help=("Generate ORDs plot"
                            +" Optionally, pass in a custom filename (relative to savepath)."
                            +" (%(default)s)"))

parser.add_argument("--latex",dest='latex',action='store_true',default=False
                    , help="Output report in latex. (%(default)s)")

parser.add_argument("--code-carrier-pairs", "--cc-pairs", dest="cc", metavar='cc'
                    , default=["L1:CA,L2:CM"], nargs='*', type=codeCarrierCombos
                    , help="List of dual frequency pairs to use. (%(default)s)")

parser.add_argument("--exSats", dest="exsats", default=None, nargs='*', type=int
                    , help="Specify SVs to exclude from analysis. (%(default)s)")

parser.add_argument("--snr", dest="snr", type=float, default=40
                    , help="Specify minimum SNR in decibels. (%(default)s)")

parser.add_argument("--include-unlocked", dest="unlocked", action="store_true"
                    , help="Don't ignore data without code carrier lock status. (%(default)s)")

parser.add_argument("--start", dest="start", default=""
                    , help="Set custom start time. (obs file start)")

parser.add_argument("--stop", dest="stop", default=""
                    , help="Set custom stop time. (obs file end)")

parser.add_argument("--timefmt", dest="timefmt", default="%Y-%m-%d %H:%M:%S"
                    , help="Set input time format string. (%(default)s)")

parser.add_argument('--clockfit', dest="clockfit", metavar='clk', default='spline'
                    , choices=["spline","average"]
                    , help="Method to apply when estimating receiver clock offset. (%(default)s)")

parser.add_argument('--percentile', dest="percentile", metavar='pct', default=99.0, type=float
                    , help="Percentage of the data to plot. (%(default)s)",)

parser.add_argument('--cadence', dest="cadence", metavar='cad', default=1000, type=int
                    , help="Number of seconds to fit over for clock offset estimation. (%(default)s)")

# parser.add_argument("--title", dest="title"
#                     , help="Specify title for graph. (Pseudorange residuals/Range-Phase)")

parser.add_argument("--log", dest="log", default='/dev/stdout'
                    , help="Where to output analysis. (%(default)s)")

parser.add_argument("--savepath", dest="savepath", metavar='savepath', default=None
                    , help="Instead of displaying graph on screen, save to indicated file. (%(default)s)")

parser.add_argument("--rmp", dest="rmp", default = False, action='store_true'
                    , help="Plot the range-phase as opposed to expected-observed range. (%(default)s)")

parser.add_argument("-v", dest="verbose", action='count',default=1
                    , help=("Print status information to command line."
                           +" Enter multiple 'v's for more verbosity. (v)"))

parser.add_argument("-d", "--debug", default=0, dest="debug", action="count"
                    , help="Increase the debug. (%(default)s)")

parser.add_argument("-q","--quiet",dest='quiet',action='store_true'
                    , help="Disable verbose output. (%(default)s)")


# ============================================================
# Parse user options
# ============================================================
opts = parser.parse_args()

VERBOSE = 0 if opts.quiet else opts.verbose
mtools.VERBOSE = VERBOSE
DEBUG = opts.debug

if DEBUG > 1: print opts
if opts.start: opts.start = timegps(mtools.parseDate(opts.start,opts.timefmt))
if opts.stop:  opts.stop = timegps(mtools.parseDate(opts.stop,opts.timefmt))

# Open all necessary files    
ephReader = mtools.EphReader(*[fn for path in opts.ephemeris for fn in glob.glob(path)])
ephReader.rejectBadClocks(True)
ephReader.rejectBadPositions(True)
ephReader.rejectPredClocks(False)
ephReader.rejectPredPositions(False)

m = mdh.MDHFileList(*[fn for path in opts.obs for fn in glob.glob(path) if h5py.is_hdf5(fn)],mode='r')

if VERBOSE > 1 and '/tmp/' in opts.expobs:
    print "Creating temporary expected obs file: %r" % opts.expobs

try:
    expobs_fh = mdh.MDHFile(opts.expobs,'r' if opts.read_only else 'a')
except IOError:
    if not os.path.exists(opts.expobs):
        raise
    if VERBOSE:
        print "Warning: Insufficient priveleges to write to %s." % opts.expobs
        print "Opening as read-only."
    expobs_fh = mdh.MDHFile(opts.expobs,'r')

# Determine time span for analysis
minStartTime, maxStopTime = mtools.mdhSpan(m,absolute=True,SUBCLASS="OBSERVATIONS")
if opts.start:  minStartTime = max(opts.start,minStartTime)
if opts.stop:   maxStopTime = min(opts.stop,maxStopTime)

# force start/stop times to second boundaries
minStartTime = np.ceil(minStartTime)
maxStopTime = long(maxStopTime)

if VERBOSE:
    print "Analysis time span is [%s,%s)" % (gpstime(minStartTime),gpstime(maxStopTime))

# define some constants needed throughout the program
c = gpstk.C_MPS
f = {"L1":gpstk.L1_FREQ_GPS,"L2":gpstk.L2_FREQ_GPS,"L5":gpstk.L5_FREQ_GPS}
bins = ((0,5), (5,10), (10,30), (30,60), (60,90), (10,90))

# ============================================================
# Main
# ============================================================
ordsdata=[]
cadence = None
timeDenominator = None

for cc1,cc2 in opts.cc:
    a1 = dict(zip(('carrierCode', 'rangingCode'),cc1.split(':')))
    a2 = dict(zip(('carrierCode', 'rangingCode'),cc2.split(':')))

    for obs1 in sorted(m.matchingDatasets(SUBCLASS="OBSERVATIONS", **a1), key=lambda x: x.attrs['startTime']):
        if opts.exsats and obs1.attrs['prnCode'] in opts.exsats: continue

        obsStartTime = max(long(minStartTime*obs1.attrs['timeDenominator']),obs1.attrs['startTime'])
        obsStopTime = min(long(maxStopTime*obs1.attrs['timeDenominator']),mtools.getStopTime(obs1))
        if obsStopTime <= obsStartTime: continue

        if DEBUG:
            print "{} ({}), PRN:{:02d} (cc:{} rc:{})".format(
                york2dt(obs1.attrs['startTime']), len(obs1['pseudorange']),
                obs1.attrs['prnCode'], obs1.attrs['carrierCode'], obs1.attrs['rangingCode']),
            
        # Now to find another obs from this SV for this time.
        a2['prnCode'] = obs1.attrs['prnCode']
        for obs2 in m.matchingDatasets(SUBCLASS="OBSERVATIONS"
                                       ,timeDenominator = obs1.attrs['timeDenominator']
                                       ,cadence = obs1.attrs['cadence']
                                       ,**a2):
            startTime = max(obsStartTime, obs2.attrs['startTime'])
            stopTime = min(obsStopTime, mtools.getStopTime(obs2))
            if startTime < stopTime: break
        else:
            if DEBUG: print "No second OBSERVATION found."
            continue

        c1 = mdh.getEnumAttributeString(obs1, 'carrierCode')
        c2 = mdh.getEnumAttributeString(obs2, 'carrierCode')
        prn = obs1.attrs['prnCode']
        if DEBUG:
            print "- (cc:{} rc:{})".format(obs2.attrs['carrierCode'], obs2.attrs['rangingCode'])
        if DEBUG>1:
            print "OBS1:", obs1.attrs.items()
            print "OBS2:", obs2.attrs.items()

        # Find the previously written expobs dataset
        for expobs in expobs_fh.matchingDatasets(SUBCLASS="EXPOBS",
                                                 prnCode = obs1.attrs['prnCode'],
                                                 timeDenominator = obs1.attrs['timeDenominator'],
                                                 cadence = obs1.attrs['cadence'] ):
            ephStartTime = max(startTime, expobs.attrs['startTime'])
            ephStopTime = min(stopTime, mtools.getStopTime(expobs))
            if ephStartTime < ephStopTime:
                startTime = ephStartTime
                stopTime = ephStopTime
                break
        else:
            if '+' in expobs_fh.mode:
                expobs = writeExpObs(startTime, stopTime, expobs_fh, ephReader
                                     , obs1, obs2, opts.rxpos)
                if expobs is None: continue
            else:
                if DEBUG: print "No EXPOBS found."
                continue
        if DEBUG > 1: print "EXPOBS:", len(expobs), expobs.attrs.items()
            
        # Now find the common chunk of data between the tracks
        cadence = obs1.attrs['cadence']
        timeDenominator = obs1.attrs['timeDenominator']
        n = (stopTime - startTime)/cadence
        i = (startTime - obs1.attrs['startTime'])/cadence
        obs1 = obs1[i:i+n]
        i = (startTime - obs2.attrs['startTime'])/cadence
        obs2 = obs2[i:i+n]
        i = (startTime - expobs.attrs['startTime'])/cadence
        expobs = expobs[i:i+n]

        mask = np.ones(n,dtype=np.bool)
        if not opts.unlocked: #removes all data without code & carrier lock
            mask &= (obs1['demodulatorStatus']==2) & (obs1['lockCount']>0) 
            mask &= (obs2['demodulatorStatus']==2) & (obs2['lockCount']>0)
        mask &= obs1['snr']>opts.snr
        mask &= obs2['snr']>opts.snr

        # Kill the bad data!
        obs1 = obs1[mask]
        obs2 = obs2[mask]
        expobs = expobs[mask]
        maskidx = np.arange(n)[mask]
        if len(maskidx) == 0: continue

        # Now we have the data cleaned up. Do math!
        eobs = expobs['pseudorange']
        elev = expobs['elevation']
        pr1 = obs1['pseudorange']
        adr1 = obs1['accumulatedDeltaRange']
        pr2 = obs2['pseudorange']
        adr2 = obs2['accumulatedDeltaRange']
        lamda1 = c / f[c1]
        lamda2 = c / f[c2]
        alpha = f[c1]**2/(f[c1]**2-f[c2]**2)
        beta = f[c2]**2/(f[c1]**2-f[c2]**2)
        icpr =  alpha*pr1 - beta*pr2 # iono free pseudorange
        icrmp = alpha*(pr1-adr1*lamda1) - beta*(pr2-adr2*lamda2)  # iono free range minus phase
        icrmp = icrmp - np.median(icrmp) # debias
        ord = icpr - eobs # remove expected pseudorange
        # ord_med = np.median(ord)
        # ord -= ord_med #don't want to do this on PR since it makes it harder to remove common clock offset
        if opts.rmp: data = icrmp
        else:        data = ord

        mask = mask[maskidx[0]:maskidx[-1]+1]   # crop the mask to the valid start and stop times
        ordsdata.append({'startTime':maskidx[0]*cadence+startTime
                         , 'stopTime':(maskidx[-1]+1)*cadence+startTime
                         , 'ordsData':data, 'ordsTimeMask': mask, 'prn':prn, 'elev':elev
                         , 'frequencyPair':(cc1,cc2)})

# ============================================================
# Plotting and text output
# ============================================================
if ordsdata:
    mergedTimes, avgRange, mergedMask = AverageRanges(ordsdata,cadence,timeDenominator)
    maskedFit = avgRange[mergedMask]
    mplTimes = mdh.mdhtime2num(mergedTimes[mergedMask],timeDenominator)

    # set up figure
    cmap = plt.cm.spectral
    legendprops = {'fancybox':True,'markerscale':10,'numpoints':1,'prop':{'size':'small'}}

    fig = plt.figure(figsize=(10.24,8.0))
    axes = {}
    axes['ords'] = fig.add_subplot(311)
    axes['clockEst'] = fig.add_subplot(312,sharex=axes['ords'])
    axes['clockRes'] = fig.add_subplot(313,sharex=axes['ords'])

    # generate and plot clock fit
    winsize = int(opts.cadence*timeDenominator/float(cadence))
    slicestart_ma = winsize/2
    slicestop_ma = slicestart_ma + len(maskedFit) - winsize + 1

    if opts.clockfit == "spline":
        rows = np.arange(len(maskedFit))
        knotRows = np.arange(winsize/2,len(maskedFit)-winsize/2,winsize)
        clockFit = LSQUnivariateSpline(rows,maskedFit,knotRows,k=3)
        axes['ords'].plot_date(mplTimes,clockFit(rows),'-',color='orange'
                               ,alpha=0.8,label="Spline Fit")
    elif opts.clockfit == "average":
        clockFit = np.convolve(maskedFit,[1.0/winsize]*winsize,mode='valid')

        axes['ords'].plot_date(mplTimes[slicestart_ma:slicestop_ma]
                               , clockFit,'-',color='orange',alpha=0.8,label="Running Avg")

    # plot clock estimate in its own axis
    if opts.clockfit == "spline":    
        axes['clockEst'].plot_date(mplTimes, maskedFit, '-', color='red', alpha=0.8
                                   , label="Est. Clock")
        axes['clockEst'].plot_date(mplTimes, clockFit(rows), '-', color='orange'
                                   , alpha=0.8, label="Spline Fit")
    elif opts.clockfit == "average":
        axes['clockEst'].plot_date(mplTimes[slicestart_ma:slicestop_ma]
                                   , maskedFit[slicestart_ma:slicestop_ma]
                                   ,'-',color='red',alpha=0.8, label="Est. Clock")
        axes['clockEst'].plot_date(mplTimes[slicestart_ma:slicestop_ma]
                                   , clockFit
                                   , '-', color='orange', alpha=0.8, label="Running Avg")
    axes['clockEst'].set_ylabel("est. clock and fit (m)")
    
    # plot ORDS and residuals = ORDS - fit to est. clock
    sumsqerr = 0
    sumlen = 0
    b_data= { b: np.array(0) for b in bins }
    
    for ords in ordsdata:
        startidx = (ords['startTime'] - mergedTimes[0])/cadence
        stopidx = startidx + (ords['stopTime']-ords['startTime'])/cadence
        maskedFit = avgRange[startidx:stopidx][ords['ordsTimeMask']]
        mplTimes = mdh.mdhtime2num(mergedTimes[startidx:stopidx][ords['ordsTimeMask']],timeDenominator)

        freqpair = '-'.join(map(lambda x: str.replace(x,':',''),ords['frequencyPair']))
        axes['ords'].plot_date(mplTimes, ords['ordsData'], '.', ms=1, color=cmap(ords['prn']/32.0)
                               , label="PRN%2d %s"%(ords['prn'],freqpair))

        if opts.clockfit == "spline":
            rows = np.arange(len(maskedFit))
            knotRows = np.arange(winsize/2,len(rows)-winsize/2,winsize)

            try:
                clockFit = LSQUnivariateSpline(rows,maskedFit,knotRows,k=3)
            except:
                if DEBUG: print "ERROR: Too few points to fit clock (%d)." % len(maskedFit)
                continue
            res = ords['ordsData'] - clockFit(rows)

            elev = ords['elev']
        elif opts.clockfit == "average":
            slicestop_ma = slicestart_ma + len(maskedFit) - winsize + 1

            avgfit = np.convolve(maskedFit,[1.0/winsize]*winsize,mode='valid')
            res = ords['ordsData'][slicestart_ma:slicestop_ma] - avgfit
            elev = ords['elev'][slicestart_ma:slicestop_ma]            
            mplTimes = mplTimes[slicestart_ma:slicestop_ma]

        if res.size:
            axes['clockRes'].plot_date(mplTimes, res, '.', ms=1, color=cmap(ords['prn']/32.0))

        # accumulate error for this fit
        sumlen += len(res)
        sumsqerr += sum(res**2)            

        # generate binned residual data
        for b in bins:
           b_data[b] = np.append(b_data[b], res[(elev > b[0]) & (elev <= b[1])])
    rmse = np.sqrt(sumsqerr/sumlen)

    # generate elevation binned ORDs report
    with open(opts.log,'w') as logger:
        report = ORDReport(LatexTable if opts.latex else TextTable)
        report.generateContent(bins,b_data).writeReport(out=logger)

    # 90% of the data means all data between 5% and 95%
    alldata = np.concatenate([ords['ordsData'] for ords in ordsdata])
    half_pct = 100.0 - (100.0 - opts.percentile)/2.0
    llim = np.percentile(alldata, 100-half_pct)
    ulim = np.percentile(alldata, half_pct)
    axes['ords'].set_ylim(llim, ulim)
    if opts.rmp:    axes['ords'].set_ylabel(" RMP (m)")
    else:           axes['ords'].set_ylabel(" ORDS (m)")
    axes['clockRes'].set_ylabel("ORDS-clock fit (m)")
    axes['clockRes'].annotate("RMSE=%.2f" % rmse,xy=(0.01,0.05)
                              ,xycoords='axes fraction',fontsize=10,color='r')

    for key,ax in axes.iteritems():
        if key != 'ords': ax.legend(**legendprops)
        ax.grid()
        ax.xaxis.set_major_formatter(dates.DateFormatter("%H:%M"))
        ax.xaxis.set_major_locator(dates.HourLocator(interval=2))

    txt = "{:3.2f}% data: [{:.3f}, {:.3f}]:".format(opts.percentile, llim, ulim)
    if opts.unlocked:       txt += ",  unlocked data included"
    else:                   txt += ",  unlocked data stripped"
    plt.figtext(0.01, 0.046, txt, size="small")
    b = bins[-1]
    plt.figtext(0.01, 0.018, "{:2d}-{:2d} stddev:{:6.3f} m".format(b[0], b[1], np.std(b_data[b]))
                , size="small")
    plt.figtext(0.3 , 0.018, "MAD:{:6.3f} m".format(MAD(b_data[b])), size="small")

    handles, labels = axes['ords'].get_legend_handles_labels()
    by_label = OrderedDict(sorted(zip(labels, handles)))
    axes['ords'].legend(by_label.values(), by_label.keys(),bbox_to_anchor=(0., 1.0, 1., .05)
                        ,ncol=min(len(by_label),4), borderpad=0.1, loc='lower left'
                        ,mode='expand',handletextpad=0.05,**legendprops)
    axes['ords'].set_xlim(gpstime(minStartTime),gpstime(maxStopTime))

    # if opts.title is None:
    #     if opts.rmp:
    #         opts.title="Range-Phase"
    #     else:
    #         opts.title="Pseudorange residuals"
    # fig.suptitle(opts.title)

    if opts.savepath is None:
        plt.show()
    else:
        if opts.plot_ords:
            filetag = opts.plot_ords
        else:
            filetag = gpstime(minStartTime).strftime("%Y%m%d-%H%M%S")
            filetag += '-' + '_'.join( opts.title.lower().split(' ') )
        filetag += '.png' if '.png' not in filetag else ''
        plt.savefig(os.path.join(opts.savepath,filetag), dpi=100)

expobs_fh.close()
if '/tmp/' in opts.expobs:
    if VERBOSE > 1: print "Removing temporary expobs file %r..." % opts.expobs
    os.remove(opts.expobs)
