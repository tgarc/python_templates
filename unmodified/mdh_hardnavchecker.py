#!/usr/bin/python
'''
hardnavchecker

Analysis for verifying correctness of HARD_BITS type mdh datasets
'''

import sys, os
import mdh
import mdh_tools as mtools
from mdh_tools import Report,LatexTable,TextTable
import h5py
import numpy as np
import mdh_daa as daa
import gpstk
from datetime import datetime, timedelta
from itertools import groupby, product
import matplotlib.pyplot as plt
import random
from numpy.lib.stride_tricks import as_strided

from collections import OrderedDict
import operator as op

VERBOSE = 0
DEBUG = 0

statusCodes = ('Ephemeris','HardNav','frameAligned','parityPassed','timeAligned')
statusCodes = OrderedDict(zip(statusCodes,2**np.arange(len(statusCodes))))


class MDPNavReader(object):
    def __init__(self,fn):
        keys = ['time','310', 'prn', 'carrier_code', 'range_code', 'nav_code']
        keys += ["word%d" % i for i in range(1,11)]
        typecasts = {0: lambda x: mtools.timegps(datetime.strptime(x,"%Y/%j/%H:%M:%S.%f"))
                     , 3: lambda x: mtools.h5t2CarrierCode.get(int(x),"UNKNOWN")
                     , 4: lambda x: mtools.h5t2RangingCode.get(int(x),"UNKNOWN")}

        types = [ float, np.uint16, np.uint8, 'S8', 'S8', np.uint32 ] + [ np.uint32 ]*10
        typecasts.update({i: (lambda x: int(x,16)) for i in range(6,len(keys))})

        self.values = np.loadtxt(fn,dtype=zip(keys,types),converters=typecasts,delimiter=", ")

        self.values['word1'][1:] = self.values['word1'][:-1]
        self.values = self.values[:-1]

    def toMDH(self, mdhIn, mdhOut, **attrs):
        attrs['SUBCLASS'] = "HARD_BITS"
        for navTrack in mdhIn.matchingDatasets(**attrs):
            prn = navTrack.attrs['prnCode']
            cc = mdh.getEnumAttributeString(navTrack,'carrierCode')
            rc = mdh.getEnumAttributeString(navTrack,'rangingCode')
            startTime= navTrack.attrs['startTime']/navTrack.attrs['timeDenominator']
            stopTime = mtools.getStopTime(navTrack)/navTrack.attrs['timeDenominator']

            if DEBUG:
                print "(MDH) {}, ({}) epochs, PRN:{:02d} (cc:{} rc:{}).".format(
                    mtools.gpstime(startTime),len(navTrack),prn,cc,rc),

            startTime, data = self.DecodeRawMDPNav(prn,cc,rc,startTime,stopTime)

            if not data.size:
                if DEBUG: print "No MDP Nav data found..."
                continue
            if DEBUG:
                print "(MDP+MDH) %s, (%d) epochs, PRN%02d(%2s:%2s)" % (mtools.gpstime(startTime),len(data),prn,cc,rc)

            dsetOut = mdhOut.create_dataset(os.path.basename(navTrack.name)+"_MDPNAV",maxshape=(None,)
                                            ,data=data,chunks=True)

            for attr,val in navTrack.attrs.iteritems():
                dt = h5py.h5a.open(navTrack.id,attr).dtype
                dsetOut.attrs.create(attr,val,dtype=dt)
            dsetOut.attrs['SUBCLASS'] = "MDPNAV"
            dsetOut.attrs['startTime'] = startTime*navTrack.attrs['timeDenominator']

    def DecodeRawMDPNav(self,prn,cc,rc,startTime=None,stopTime=None):
        # create a mask to capture the bits used for computing parity (D25-D30)
        # paritymask  = np.zeros((6,24),dtype=np.bool)
        # parityidx    = [ [0,1,2,4,5,9 ,10,11,12,13,16,17,19,22]
        #                 ,[1,2,3,5,6,10,11,12,13,14,17,18,20,23]
        #                 ,[0,2,3,4,6,7 ,11,12,13,14,15,18,19,21]
        #                 ,[1,3,4,5,7,8 ,12,13,14,15,16,19,20,22]
        #                 ,[0,2,4,5,6,8 ,9 ,13,14,15,16,17,20,21,23]
        #                 ,[2,4,5,7,8,9 ,10,12,14,18,21,22,23] ]
        # for i,row in enumerate(parityidx): paritymask[i,row] = True
        # paritybits = np.hstack([np.repeat(paritymask,6,axis=1) for word in framebits[i].reshape(10,30)])

        sigmask =   (self.values['prn'] == prn) & (self.values['carrier_code'] == cc) & (self.values['range_code'] == rc)
        if startTime:   sigmask &= (self.values['time'] >= startTime)
        if stopTime:    sigmask &= (self.values['time'] < stopTime)
        sigdata = self.values[sigmask]

        framebits = np.zeros(len(sigdata),dtype=[('data',np.uint8,(300,))])
        # framebits = np.zeros((len(sigdata),),dtype=[('data',np.uint8,(300,)),('parityPassed',np.uint8)])
        for i,subframe in enumerate(sigdata[['word%d' % i for i in range(1,11)]]):
            # Decode each word of the subframe based on section 20.3.5.2 of IS-GPS-200H

            # Convert 10 words into flat array of binaries
            framebits['data'][i] = np.array([0x1 & (w >> (29-k)) for w,k in product(subframe,xrange(30))],dtype=np.uint8)
            # framebits['data'][i] = (subframe.repeat(30) >> np.tile(np.arange(29,-1,-1,dtype=np.uint32),10)) & 0x1

            # upright the bits first
            if np.all(framebits['data'][i,:8] ^ np.unpackbits(np.array([0x8B],dtype=np.uint8))):
                # 0x1^ on integer arrays of zero's and one's is the equivalent of ~ on boolean masks
                framebits['data'][i] = 0x1^framebits['data'][i] 

            # grab the data bits for each word (as a memory view)
            datawords = as_strided(framebits['data'][i],shape=(10,24),strides=(30,1))

            # xor the last 9 words with the 30th bit of the previous word for the first 9 words of the subframe
            datawords[1:] ^= framebits['data'][i,29::30].reshape(10,1)[:-1]

            # perform parity check        
            # framebits['parityPassed'][i] = True
            # d29,d30 = framebits['data'][i,28:30]
            # for word in framebits[i]['data'].reshape(10,30)[1:]: # take a word level view of this subframe
            #     paritybits = word[24:].copy()
            #     databits = np.repeat([word[:24]],6,axis=0)
            #     databits[~paritymask] = 0
            #     word[24:] = np.bitwise_xor.reduce( np.hstack( (np.array([d29,d30,d29,d30,d30,d29],dtype=np.bool).reshape(-1,1), databits) ), axis=1)
            #     framebits[i]['parityPassed'] &= np.all(paritybits == word[24:])
            #     d29,d30 = word[28:30]

        return (min(sigdata['time']), framebits) if sigdata.size else (None, np.array([]))


class NavBitCompareReport(Report):
    _keys = ("PRN","CC","RC","rows","all_passed","overlapping_subframes"
             ,"overlapping_subframes&all_passed"
             # ,"tow_match"
             ,"matching_subframes")

    def generateContent(self, mdhfile, mdpfile, sigTuples=None):
        self.bits = self.bitCompare(mdhfile,mdpfile)

        if not sigTuples:
            sigTuples = [(sig[:2],sig[2:4],sig[4:6]) for sig in sorted(self.bits)]

        row=self.getRowDict()
        bitcmp = lambda x: 0x1^(x['mdhdata']^x['mdpdata']).any(1) | (x['mdhdata']^x['mdpdata']).all(1)
        for prn,cc,rc in sigTuples:
            svkey = "%02d%2s%2s" % (prn,cc,rc)
            cmpbits = self.bits.get(svkey,[])

            # row['tow_match'] = sum(subframes['tow_match'] for subframes in cmpbits[sig])
            row['PRN'],row['CC'],row['RC'] = prn,cc,rc
            row['rows'] = sum(subframes['rows'] for subframes in cmpbits)
            row['all_passed'] = sum(subframes['all_passed'] for subframes in cmpbits)
            row['overlapping_subframes'] = sum(len(subframes['mdhdata']) for subframes in cmpbits)
            row['overlapping_subframes&all_passed'] = sum(sf['overlapping_subframes&all_passed'] for sf in cmpbits)

            # accept both polarity matched and opposite polarity subframes as matches
            row['matching_subframes'] = sum(np.sum(bitcmp(subframes)) for subframes in cmpbits)
            self.addRow(**row)

        return self

    @staticmethod
    def bitPlotter(startTime,stopTime,bitmasks):
        fig = plt.figure()
        fig.set_size_inches(14.40,6.00)
        ax = fig.add_subplot(111)

        for i,key in enumerate(sorted(bitmasks)):
            ax.plot(bitmasks[key]['times'],bitmasks[key]['mask']+1.25*i+0.25,drawstyle='steps-post')

        ax.yaxis.set_ticks(np.arange(len(bitmasks))*1.25+0.75)
        ax.yaxis.set_ticklabels(sorted(bitmasks.keys()),fontsize=9,family='monospace')
        ax.set_ylim(0,len(bitmasks)*1.25+0.25)

        fig.autofmt_xdate()    
        ax.set_xlim(startTime,stopTime)    
        ax.set_xticks([startTime + timedelta(milliseconds=600*x) for x in range(0, 10)], minor=True)
        ax.xaxis.grid(which='minor')

        fig.subplots_adjust(bottom=0.08,left=0.09,right=0.99,top=0.95)
        plt.setp(ax.get_xticklabels(),fontsize=9, rotation=15)

        return fig

    def plot_bits(self,savepath=None,n=5,sigTuples=None):
        if not sigTuples:
            sigTuples = [(sig[:2],sig[2:4],sig[4:6]) for sig in sorted(self.bits)]
        
        for prn,cc,rc in sigTuples:
            cmpbits = self.bits[prn+cc+rc]
            cmpbits = cmpbits[random.randint(0,len(cmpbits)-1)]

            bitcmp = cmpbits['mdpdata'] ^ cmpbits['mdhdata']
            mdpdata = cmpbits['mdpdata']
            mdhdata = cmpbits['mdhdata']

            matches = 0x1^bitcmp.any(1) | bitcmp.all(1)
            p = min( len(bitcmp)-np.sum(matches) , n )
            mismatchidx = np.nonzero(~matches)[0]

            while p > 0:      # select (up to) n random mismatched subframes to plot
                idx = mismatchidx[random.randrange(len(mismatchidx))]

                plotStartTime = (mtools.timegps(cmpbits['startTime'])*cmpbits['timeDenominator']
                                 + idx*cmpbits['cadence'])
                times = mdh.mdhtime2num(plotStartTime+np.arange(300)*(cmpbits['cadence']/300.)
                                        ,cmpbits['timeDenominator'])
                fig = self.bitPlotter(mtools.gpstime(plotStartTime/cmpbits['timeDenominator'])
                                      ,mtools.gpstime((plotStartTime+1*cmpbits['cadence'])/cmpbits['timeDenominator'])
                                      ,{'MDP':{'mask':mdpdata[idx].flatten(),'times':times}
                                        ,'MDH':{'mask':mdhdata[idx].flatten(),'times':times}
                                        ,'MDP XOR MDH':{'mask':mdhdata[idx].flatten()^mdpdata[idx].flatten(),'times':times}})
                if savepath:
                    ftag = (mtools.gpstime(plotStartTime/cmpbits['timeDenominator']).strftime("%Y%m%d-%H%M%S"),prn,cc,rc)
                    fig.savefig(os.path.join(savepath,'%s_PRN%2s%2s%2s.png' % ftag),dpi=100)
                else:
                    plt.show(fig.number)
                p -= 1

    @staticmethod
    def bitCompare(mdhfile, mdpfile, defaultCadence=6000):
        cmpbits = {}
        for mdhnav in mdhfile.matchingDatasets(SUBCLASS="HARD_BITS"):
            if 'cadence' in mdhnav.attrs:
                cadence = mdhnav.attrs['cadence']
            else:
                cadence = defaultCadence
            timeDenominator = mdhnav.attrs['timeDenominator']

            if (mdhnav.attrs['startTime']%cadence) != 0: continue #skip non-time-aligned ds's

            if DEBUG:
                print "(MDH) {}, ({}) epochs, PRN:{:02d} (cc:{} rc:{}).".format(
                    mtools.gpstime(mdhnav.attrs['startTime']/mdhnav.attrs['timeDenominator']),len(mdhnav)
                    ,mdhnav.attrs['prnCode'], mdhnav.attrs['carrierCode'], mdhnav.attrs['rangingCode']),

            sig = {'prnCode':mdhnav.attrs['prnCode']
                   ,'carrierCode':mdhnav.attrs['carrierCode']
                   ,'rangingCode':mdhnav.attrs['rangingCode']
                   ,'cadence':cadence}
            for mdpnav in mdpfile.matchingDatasets(SUBCLASS="MDPNAV",**sig):
                startTime = max(mdpnav.attrs['startTime'], mdhnav.attrs['startTime'])
                stopTime = min(mtools.getStopTime(mdpnav), mtools.getStopTime(mdhnav))
                n = int((stopTime - startTime)/cadence)
                if n > 0: break
            else:
                if DEBUG: print "No MDP Nav data found..."
                continue

            if DEBUG:
                print "(MDP+MDH) {}, ({}) epochs".format(mtools.gpstime(startTime/timeDenominator),n)
            sig['rows'] = len(mdhnav)
            sig['all_passed'] = sum((mdhnav['frameAligned'] == 1) & (mdhnav['parityPassed'] == 1))

            prn = mdhnav.attrs['prnCode']
            cc = mdh.getEnumAttributeString(mdhnav,'carrierCode')
            rc = mdh.getEnumAttributeString(mdhnav,'rangingCode')

            # select the overlapping portion of data
            i = int((startTime - mdhnav.attrs['startTime'])/cadence)
            mdhnav = mdhnav[i:i+n]
            i = int((startTime - mdpnav.attrs['startTime'])/cadence)
            mdpnav = mdpnav[i:i+n]

            # mask out bad data
            mask = (mdhnav['frameAligned'] == 1) & (mdhnav['parityPassed'] == 1)
            mdpnav = mdpnav[mask]
            mdhnav = mdhnav[mask]
            maskidx = np.arange(n,dtype=np.uint32)[mask]
            if len(maskidx) == 0: continue
            stopTime = (maskidx[-1]+1)*cadence+startTime
            startTime += maskidx[0]*cadence

            # TOWcnt_ts = (((int(startTime)+maskidx*cadence)/timeDenominator) % gpstk.FULLWEEK) // 6 + 1
            # TOWcnt_nav = np.sum(mdpnav['data'][:,30:47] << np.arange(16,-1,-1,dtype=np.uint32),axis=1)            
            # sig['tow_match'] = np.sum(TOWcnt_nav != TOWcnt_ts)
            # if (TOWcnt_ts != TOWcnt_nav):
            #     print mtools.gpstime(startTime/timeDenominator)
            #     print TOWcnt_ts, TOWcnt_nav
            #     print np.binary_repr(TOWcnt_ts ^ TOWcnt_nav)

            sig['mdpdata'] = mdpnav['data']
            sig['mdhdata'] = mdhnav['data']
            sig['cadence'] = cadence
            sig['timeDenominator'] = timeDenominator
            sig['overlapping_subframes'] = n
            sig['overlapping_subframes&all_passed'] = sum(mask)
            sig['startTime'] = mtools.gpstime(startTime/timeDenominator)
            sig['stopTime'] = mtools.gpstime(stopTime/timeDenominator)

            svkey = "%02d%2s%2s" % (prn,cc,rc)
            if svkey not in cmpbits: cmpbits[svkey] = []    
            cmpbits[svkey].append(sig)

        return cmpbits


class NavStatusReport(Report):
    _keys = ["PRN","CC","RC","Datasets","Rows","TimeAligned","TOWmatched"
             ,"ParityPassed","FrameAligned","AllPassed"]

    def generateContent(self, mdhfile, sigTuples=None):
        """
        Prints out a summary report on the status of the nav bits dataset
        """        
        sigrow = self.getRowDict()
        counter = dict(zip(self._keys[3:],[0]*len(self._keys)))
        numpct = lambda num,denom: "%d (%3.0f%%)" % (num, 100.*num/denom)

        if not sigTuples: sigTuples = mtools.matchingSignals(self.mdhfile,SUBCLASSES=["HARD_BITS"])

        navStatus = {}
        for prn,cc,rc in sigTuples:
            sig = {'prnCode':prn,'carrierCode':cc,'rangingCode':rc}
            sigrow = self.getRowDict(fill=0)
            sigrow.update(zip(("PRN","CC","RC"),(prn,cc,rc)))

            for trk in mdhfile.matchingDatasets(SUBCLASS="HARD_BITS", **sig):
                cadence = 1 if 'cadence' not in trk.attrs else trk.attrs['cadence']
                timeAligned = (trk.attrs['startTime']%cadence) == 0

                TOWcnt_ts = (((mtools.getStopTime(trk)/trk.attrs['timeDenominator']) % gpstk.FULLWEEK) // 6) + 1
                TOWcnt_nav = np.sum(trk['data'][:,30:47] << np.arange(16,-1,-1,dtype=np.uint32),axis=1)
                sigrow['TOWmatched']    += np.sum(TOWcnt_nav != TOWcnt_ts)
                sigrow['Datasets']      += 1
                sigrow['TimeAligned']   += timeAligned
                sigrow['Rows']          += len(trk)
                sigrow['ParityPassed']  += np.sum(trk['parityPassed'])
                sigrow['FrameAligned']  += np.sum(trk['frameAligned'])

                if timeAligned:
                    sigrow['AllPassed'] += np.sum(trk['parityPassed'] & trk['frameAligned'])

            # increment the table counter
            for k in counter: counter[k]+=sigrow[k]

            sigrow['TOWmatched']    = numpct(sigrow['TOWmatched'],sigrow['Rows'])
            sigrow['ParityPassed']  = numpct(sigrow['ParityPassed'],sigrow['Rows'])
            sigrow['FrameAligned']  = numpct(sigrow['FrameAligned'],sigrow['Rows'])
            sigrow['AllPassed']     = numpct(sigrow['AllPassed'],sigrow['Rows'])
            self.addRow(**sigrow)

        # Print totals
        counter['TimeAligned'] = numpct(counter['TimeAligned'],counter['Datasets'])
        counter['ParityPassed'] = numpct(counter['ParityPassed'],counter['Rows'])
        counter['FrameAligned'] = numpct(counter['FrameAligned'],counter['Rows'])
        counter['AllPassed'] = numpct(counter['AllPassed'],counter['Rows'])

        self.addLineBreak()
        self.addRow(**counter)
        return self


class NavDetailReport(Report):
    _keys = ("Filename","Dataset","StartTime","TimeAligned","Rows"
              ,"ParityFail","FrameAlignFail")
    _fmt = {}
    _fmt['ParityFail'] = _fmt['FrameAlignFail'] = "[{0[0]:4d},{0[1]:4d})"

    def generateContent(self, mdhfile):
        dsrow = {}
        for trk in mdhfile.matchingDatasets(SUBCLASS="HARD_BITS"):
            dsrow.clear()
            cadence = 1 if 'cadence' not in trk.attrs else trk.attrs['cadence']

            dsrow['Filename']      = os.path.basename(trk.file.filename)
            dsrow['Dataset']       = trk.name
            dsrow['StartTime']     = mtools.gpstime(trk.attrs['startTime']/trk.attrs['timeDenominator'])
            dsrow['TimeAligned']   = (trk.attrs['startTime']%cadence) == 0
            dsrow['Rows']          = len(trk)

            if dsrow['TimeAligned'] and np.all((trk['parityPassed'] == 1) & (trk['frameAligned'] == 1)):
                continue

            pfail_runs = zip(*mtools.mask2runs(trk['parityPassed'] != 1))
            ffail_runs = zip(*mtools.mask2runs(trk['frameAligned'] != 1))

            # print the first run of failed epochs on the same line as other output
            if len(ffail_runs) > 0:
                dsrow['FrameAlignFail']    = ffail_runs[0]
                ffail_runs = ffail_runs[1:]
            if len(pfail_runs) > 0:
                dsrow['ParityFail']        = pfail_runs[0]
                pfail_runs = pfail_runs[1:]
            self.addRow(**dsrow)

            # print the remaining runs of failed epochs on their own line
            for pitv,fitv in zip(pfail_runs,ffail_runs):
                self.addRow(ParityFail=pitv, FrameAlignFail=fitv)

            remainder = len(pfail_runs)-len(ffail_runs)
            if remainder != 0:
                fail_runs   = pfail_runs      if remainder>0 else ffail_runs
                kw          = 'ParityFail'    if remainder>0 else 'FrameAlignFail'
                for itv in fail_runs[-abs(remainder):]: self.addRow(**{kw:itv})

        return self


def getStatus(mdhfile,prn,cc,rc,startTime=None,stopTime=None):
    sig = dict(zip(('prnCode','carrierCode','rangingCode'),(prn,cc,rc)))
    groupCadence = min(t.attrs['cadence'] for t in mdhfile.matchingDatasets(SUBCLASS="HARD_BITS", **sig) if 'cadence' in t.attrs)
    groupStart,groupStop = mtools.mdhSpan(mdhfile,SUBCLASS="HARD_BITS", **sig)

    if startTime:   groupStart = max(groupStart,startTime)
    if stopTime:    groupStop = min(groupStop,stopTime)
    if groupStop <= groupStart: return ()

    statusMask = np.zeros((groupStop-groupStart)/groupCadence,dtype=np.uint32)

    for trk in mdhfile.matchingDatasets(SUBCLASS="HARD_BITS", **sig):
        cadence = 1 if 'cadence' not in trk.attrs else trk.attrs['cadence']
        timeAligned = (trk.attrs['startTime']%cadence) == 0

        startTime = max(groupStart, trk.attrs['startTime'])
        stopTime = min(trk.attrs['startTime']+len(trk)*cadence, groupStop)
        nEpochs = (stopTime - startTime)/cadence
        if nEpochs <= 0: continue

        trkOffset = (startTime-trk.attrs['startTime'])/cadence
        grpOffset = (startTime - groupStart)/cadence

        statusMask[grpOffset:grpOffset+nEpochs] |= (trk['parityPassed'][trkOffset:trkOffset+nEpochs] == 1)*statusCodes['parityPassed']
        statusMask[grpOffset:grpOffset+nEpochs] |= (trk['frameAligned'][trkOffset:trkOffset+nEpochs] == 1)*statusCodes['frameAligned']
        statusMask[grpOffset:grpOffset+nEpochs] |= timeAligned*statusCodes['timeAligned']
        statusMask[grpOffset:grpOffset+nEpochs] |= statusCodes['HardNav']

    return groupStart, statusMask


if __name__ == '__main__':
    # ==========================================================
    # Define CLI
    # ==========================================================
    def position(value):
        return gpstk.Position(*map(float, value.split(",")))

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Tool for reviewing/analyzing HARD_BITS raw navigation"
                          +" message datasets. By default just lists a summary of the"
                          +" data sets and subframe status checks.")

    parser.add_argument("mdhfile", nargs='*'
                        , help="MDH files to process.")

    parser.add_argument("-p","--position", dest='rxpos', type=position
                        , help="Receiver X,Y,Z position. (required for ephemeris)")

    parser.add_argument("-e", "--eph",dest='ephemeris',default=None,action='append'
                        , help="Path to ephemeris file."
                        +" Accepts multiple ephemeris files when --eph is specified multiple times."
                        +" (required for plotting)")

    parser.add_argument("--log", dest="log", default='/dev/stdout'
                        , help="Where to output analysis. (%(default)s)")

    parser.add_argument("--latex",dest='latex',action='store_true',default=False
                        , help="Output report in latex. (%(default)s)")

    parser.add_argument("--mdp-compare",dest="mdp_compare", default=None
                        , help="Compare HARD_BITS to raw nav bits from an MDP file. (%(default)s)")

    parser.add_argument("--summary",dest="summary", default=False, action="store_true"
                        , help="List summary table of HARD_BITS status over all tracks in the mdh files. (%(default)s)")

    parser.add_argument("--list-bad-data",dest="listbaddata",action="store_true",default=False
                        , help="List detailed information about which files did not pass at least one of the checks. (%(default)s)")

    parser.add_argument("--plot-status",dest='plot_status',action='store_true',default=False
                        , help="Generate timeline plot displaying nav data status."
                        +" Include ephemeris file with --eph to include in plot. (%(default)s)")

    parser.add_argument("--plot-bits",dest='plot_bits'
                        , action='store_true',default=False
                        , help="Generate plots displaying comparison of mdp nav data"
                        +" to mdh hardnav data. (%(default)s)")
                        
    parser.add_argument("--savepath", dest='savepath', default=None
                        , help="Save directory for plot. (%(default)s/Interactive plotting)")
    
    parser.add_argument("-q","--quiet",dest="quiet",action="store_true",default=False
                        , help="Disable verbose output. (%(default)s)")

    parser.add_argument("-v", dest="verbose", action="count", default=1
                        , help="Print status information to command line."
                        + " Enter multiple 'v's for more verbosity. (%(default)s)")

    parser.add_argument("-d", "--debug", default=0, dest="debug", action="count"
                        , help="Increase the debug. (%(default)s)")
    opts = parser.parse_args()

    if not any((opts.plot_status,opts.listbaddata,opts.summary,opts.mdp_compare)):
        opts.summary = True

    VERBOSE = 0 if opts.quiet else opts.verbose
    mtools.VERBOSE = VERBOSE
    DEBUG = opts.debug

    # ==========================================================
    # Initialize
    # Open files, do some preprocessing, and define globals for analysis
    # ==========================================================
    logger = open(opts.log,'w')

    filenames = filter(h5py.is_hdf5, opts.mdhfile)
    fileshandler = mdh.MDHFileList(*filenames,mode="r")

    # a few helpers for accessing dataset attributes
    gethardbits = lambda : fileshandler.matchingDatasets(SUBCLASS="HARD_BITS")

    rowCadence = gethardbits().next().attrs['cadence']
    rowTimeDenominator = gethardbits().next().attrs['timeDenominator']

    minStartTime, maxStopTime = mtools.mdhSpan(fileshandler,absolute=True,SUBCLASS="HARD_BITS")
    sigTuples = sorted(mtools.matchingSignals(fileshandler,SUBCLASSES=["HARD_BITS"]))

    if opts.mdp_compare:
        if not h5py.is_hdf5(opts.mdp_compare): # For text MDP read, convert file to MDH first
            import re
            if VERBOSE > 1: print "Detected %r as a non-hdf5 file." % opts.mdp_compare
            if VERBOSE:     print "Converting MDP text data to HDF5 format. This may take several minutes..."

            mdptext = opts.mdp_compare
            opts.mdp_compare = re.sub(r'\..*','.mdh',mdptext) if re.search(r'\..*',mdptext) else mdptext+'.mdh'
            if os.path.exists(opts.mdp_compare):
                q = "ERROR: An mdh nav file already exists in the default save path (%r). Overwrite?"
                if mtools.getCorrectInput(q % opts.mdp_compare,default='n').lower() != 'y':
                    sys.exit("mdh nav savepath %r already exists" % opts.mdp_compare)
                
            with mdh.MDHFile(opts.mdp_compare,'w') as mdh_nav:
                if VERBOSE > 1: print "Reading MDP data..."
                mdpnav = MDPNavReader(mdptext)
                if DEBUG:       print "MDP file contains data for signals:\n"+str(np.unique(mdpnav.values[['prn','carrier_code','range_code']]))

                if VERBOSE > 1: print "Decoding raw MDP messages and writing to HDF5 format..."
                mdpnav.toMDH(fileshandler,mdh_nav)
                if VERBOSE: print "Conversion successful. File written to %r." % (os.path.abspath(opts.mdp_compare))
                del mdpnav                
        mdpnav_fh = mdh.MDHFile(opts.mdp_compare, 'r')
        if DEBUG:       print "MDP file contains signals:"
        if DEBUG:       print mtools.matchingSignals(mdpnav_fh,SUBCLASSES=["MDPNAV"])

    # ==========================================================
    # Run
    # Perform the selected tasks
    # ==========================================================
    print >> logger, "# Analysis time span is [%s,%s)" % (mtools.gpstime(minStartTime),mtools.gpstime(maxStopTime))
    print >> logger, "# Selecting nav data for signals:"
    print >> logger, "#\t"+"\n#\t".join(map(str,sigTuples))

    style = LatexTable if opts.latex else TextTable
    if opts.summary:
        print >> logger
        report = NavStatusReport(style)
        report.generateContent(fileshandler,sigTuples).writeReport(out=logger)

    if opts.listbaddata:
        print >> logger
        report = NavDetailReport(style)
        report.generateContent(fileshandler,sigTuples).writeReport(out=logger)

    if opts.mdp_compare:
        if VERBOSE: "Preparing Nav Bit comparison report..."
        print >> logger
        report = NavBitCompareReport(style)
        report.generateContent(fileshandler,mdpnav_fh,sigTuples).writeReport(out=logger)
        if opts.plot_bits: report.plot_bits(opts.savepath)
        
    if opts.plot_status:
        eph = mtools.EphReader(*opts.ephemeris)
        startDate = max(mtools.gpstime(minStartTime),mtools.commontime2datetime(eph.getInitialTime()))
        stopDate = min(mtools.gpstime(maxStopTime),mtools.commontime2datetime(eph.getFinalTime()))

        navStatus = {}
        for prn in set(prn for prn,cc,rc in sigTuples):
            start_ct = mtools.datetime2commontime(startDate)
            stop_ct = mtools.datetime2commontime(stopDate)
            start, stop, ephAboveMask, elv = daa.ProcessEphemerisAvailability(eph
                                                                              ,gpstk.SatID(int(prn))
                                                                              ,rowCadence/float(rowTimeDenominator)
                                                                              ,opts.rxpos
                                                                              ,startTime=start_ct
                                                                              ,stopTime=stop_ct)
            if not ephAboveMask.size:
                if DEBUG: print "Could not get ephemeris for prn",prn
                continue
            ephGroupStart = mtools.commontime2mdhtime(start,rowTimeDenominator)
            ephGroupStop = mtools.commontime2mdhtime(stop,rowTimeDenominator)

            for cc,rc in ((cc,rc) for sv,cc,rc in sigTuples if sv == prn):
                try:
                    startTime, statusMask = getStatus(fileshandler,prn,cc,rc
                                                      ,startTime=ephGroupStart
                                                      ,stopTime=ephGroupStop)
                except ValueError:
                    continue

                stopTime = startTime+len(statusMask)*rowCadence
                n = min((stopTime-startTime)/rowCadence,len(ephAboveMask)) # this is a hack - needs further digging
                ephoffset = (startTime-ephGroupStart)/rowCadence
                statusMask[:n] |= ephAboveMask[ephoffset:ephoffset+n]*statusCodes['Ephemeris']

                navStatus[("PRN%02d %2s%2s"%(prn,cc,rc)).ljust(10)] = {'mask': statusMask, 'startTime': mtools.gpstime(startTime/rowTimeDenominator)}

        daa.genericPlot(daa.plot_status, (navStatus, rowCadence/rowTimeDenominator)
                        , "Hard Bits Status"
                        , startDate, stopDate, opts.savepath
                        , keys=statusCodes.keys()
                        , values=statusCodes.values()
                        , colors=('#2c7bb6','#1b9e77','#d95f02','#7570b3','#e6ab02')
                        , relsize = [1.0]+[0.85]*(len(statusCodes)-1)
                        , stack=[False] + [True]*(len(statusCodes)-1)
                        , alpha=[0.4] + [1.0]*(len(statusCodes)-1))

    logger.close()
