#!/usr/bin/python

# Looks for tracks with detection metric above a given threshold and returns the time
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from optparse import OptionParser
from host import mdh
import sys,os
import numpy as np
import mdh_tools as mt

EVENT_PREALLOC_SIZE =         500
WINSIZE_SECONDS =             15.   # Size of window to detrend over
MINFILESIZE_SECONDS =         600  # don't analyze tracks shorter than 10 minutes

FILESTART_LENGTH_SECONDS =    10.
FILEEND_LENGTH_SECONDS =      0.1
DROPPEDFRAMES_MASK_SECONDS=   1.0 
LOWSNR_MASK_SECONDS =         1.0

DQ_FILESTART =               (1<<0) # Within FILESTART_LEN_SECONDS of beginning of file
DQ_FILEEND =                 (1<<1) # Within FILEEND_LENGTH_SECONDS of end of file
DQ_LOWSNR =                  (1<<2)
DQ_NEARLOWSNR =              (1<<3) # We are within LOWSNR_MASK_SECONDS of a low SNR period
DQ_DROPPEDFRAMES =           (1<<4)
DQ_NEARDROPPEDFRAMES=        (1<<5) # We are within DROPPEDFRAMES_MASK_SECONDS
                                    # of a region of dropped frames

integration_steps=1<<np.arange(15)
TAUS=0.0012288*integration_steps

def movingaverage(interval, window_size):
  window= np.ones(int(window_size))/float(window_size)
  return np.convolve(interval, window, 'same')

### TODO: make sure obs indices are offset correctly, right now they
### are relative to the span (if it is passed in)
def DisqualifyObs(obsTrack,snrthreshold=3000.,span=None):
  obsStart,obsStop =  span if span else (0,len(obsTrack))
            
  # allocate boolean array
  obsDisqualified = np.zeros(obsStop-obsStart,np.uint8)

  # use ~(statement) to include NaN values
  lowSNRobs = ~(mt.computeSmoothedSNR(obsTrack,(obsStart,obsStop),10.0) > snrthreshold)

  # Create a mask to disqualify times near LowSNR events
  nearLowSNRcumsum = np.cumsum(lowSNRobs)
  nearLowSNR = (nearLowSNRcumsum[ int(LOWSNR_MASK_SECONDS*mt.trackrate(obsTrack)) : ] 
                != nearLowSNRcumsum[ : -int(LOWSNR_MASK_SECONDS*mt.trackrate(obsTrack)) ])

  # Disqualify low snr obs times
  obsDisqualified[lowSNRobs] |= DQ_LOWSNR
  obsDisqualified[ : -int(LOWSNR_MASK_SECONDS*mt.trackrate(obsTrack)) ][nearLowSNR] |= DQ_NEARLOWSNR
  obsDisqualified[ int(LOWSNR_MASK_SECONDS*mt.trackrate(obsTrack)) : ][nearLowSNR] |= DQ_NEARLOWSNR

  # Find where codecarrier lock is lost in track
  droppedFrames = (obsTrack[obsStart:obsStop]['demodulatorStatus'] != 2)

  # Create a mask to disqualify times near dropped frames
  nearDroppedFramesCumsum = np.cumsum(droppedFrames)
  nearDroppedFrames = (nearDroppedFramesCumsum[int(DROPPEDFRAMES_MASK_SECONDS*mt.trackrate(obsTrack)):] 
           != nearDroppedFramesCumsum[:-int(DROPPEDFRAMES_MASK_SECONDS*mt.trackrate(obsTrack))])

  # Disqualify dropped frames times
  obsDisqualified[droppedFrames] |= DQ_DROPPEDFRAMES
  obsDisqualified[ : -int(DROPPEDFRAMES_MASK_SECONDS*mt.trackrate(obsTrack)) ][nearDroppedFrames] |= DQ_NEARDROPPEDFRAMES
  obsDisqualified[ int(DROPPEDFRAMES_MASK_SECONDS*mt.trackrate(obsTrack)) : ][nearDroppedFrames] |= DQ_NEARDROPPEDFRAMES

  return obsDisqualified

def FindDetects(iqTrack,tslice=None,snrthreshold=4000.,threshold=1.0,verbose=False,plot=False):
  # Offset needed for calculating detection metric
  metricOffsets = [ int( np.ceil( (3*int(tau*mt.trackrate(iqTrack))-1)/2 ) ) for tau in TAUS]
  fileStartOffset = long( np.ceil(FILESTART_LENGTH_SECONDS*mt.trackrate(iqTrack)) )
  fileEndOffset = long( FILEEND_LENGTH_SECONDS*mt.trackrate(iqTrack) )
  maxOffset = metricOffsets[-1]

  if tslice:
    iqTrackStart = max(tslice[0],fileStartOffset) + maxOffset
    iqTrackStop = min(tslice[1],len(iqTrack)-fileEndOffset) - maxOffset
  else:
    iqTrackStart = fileStartOffset + maxOffset
    iqTrackStop = len(iqTrack) - fileEndOffset - maxOffset
  
  # Disqualify start and end section of iq track
  # NOTE:may need to modify this if a trackslice is specified
  iqDisqualified = np.zeros(len(iqTrack),np.uint8)
  iqDisqualified[:fileStartOffset] |= DQ_FILESTART
  iqDisqualified[-fileEndOffset:]  |= DQ_FILEEND

  # preallocate a list for collecting events  
  eventDetects = [None]*EVENT_PREALLOC_SIZE 
  eventCnt = 0

  ## Search each observation track for the iqtrack
  for obsTrack in mdhfile.associatedDatasets(iqTrack,SUBCLASS='OBSERVATIONS'):
    times = mt.GetTimes(iqTrack,obsTrack,(fileStartOffset,len(iqTrack)-fileEndOffset))

    if not times: continue

    obsTrackStart,obsTrackStop = times['obsSlice'].start,times['obsSlice'].stop
    iqSliceStart,iqSliceStop = times['iqSlice'].start,times['iqSlice'].stop

    # Disqualify bad obs times and translate to iq times

    # get the disqualification mask for this obs track
    obsDisqualified = DisqualifyObs(obsTrack,snrthreshold=snrthreshold,span=(obsTrackStart,obsTrackStop))

    # get the indices of disqualified frames
    obsDisqualifiedTimes, = np.nonzero(obsDisqualified)

    upsampleF = int(mt.trackrate(iqTrack)/mt.trackrate(obsTrack))
    for obsIdx in obsDisqualifiedTimes:
      iqIdx = int(mt.FindIndexInOtherTrack(obsTrack,obsIdx,iqTrack))
      if iqIdx >= 0 and iqIdx <= (len(iqTrack)-upsampleF):
        iqDisqualified[iqIdx:iqIdx+upsampleF] = [obsDisqualified[obsIdx]]*upsampleF

    # Calculate window size and number of windows to plot
    lengthInSecs = (iqSliceStop-iqSliceStart)/mt.trackrate(iqTrack)
    winBlocksize= long(np.floor(WINSIZE_SECONDS*mt.trackrate(iqTrack)))
    numwindows = int(np.ceil((iqSliceStop-iqSliceStart)/float(winBlocksize)))
    winList = ( (iqSliceStart+i*winBlocksize-maxOffset,iqSliceStart+(i+1)*winBlocksize+maxOffset)
                if i < (numwindows-1)
                else (iqSliceStart+i*winBlocksize-maxOffset,iqSliceStop+maxOffset)
                for i in xrange(numwindows) )

    lastMaxDetect = None
    ## Search the data where iq overlaps with this obsTrack
    for winIdx,(winSliceStart,winSliceStop) in enumerate(winList):
      if verbose:
        sys.stdout.write("\rSearching window %d of %d..." % (winIdx+1,numwindows))
        sys.stdout.flush()

      FitData = mt.GetDetrendedPhase(iqTrack,obsTrack,(max(0,winSliceStart),min(winSliceStop,len(iqTrack))))
      detrendedPhase = FitData['detrendedPhase']
      carrierRate = mt.getCarrierRate(iqTrack)
      iqTimes = FitData['iqTimes']

      fig = plt.figure()
      # phaseAxes = fig.add_subplot(2,1,1); metricAxes = fig.add_subplot(2,1,2,sharex=phaseAxes)
      # fig.suptitle("PRN %02u, File: %s" % (iqTrack.attrs['prnCode'],os.path.basename(iqTrack.file.filename)),
      #              fontsize=14,fontweight='bold')

      ### Compute the detection metric over all taus on this window
      ## TODO: Add ability to detect multiple events in the same window
      maxDetection = (0,0,0,0,0)
      for idx,(metricOffset,tau) in enumerate(zip(metricOffsets,TAUS)):

        # find indices of this window relative to largest window
        subWinOffset = maxOffset-metricOffset
        detectMetric = mt.computeDetectionMetric(detrendedPhase[subWinOffset:-subWinOffset],tau,
                                              carrierRate,mt.trackrate(iqTrack))
        detectIndices = np.where(detectMetric >= threshold)[0]
        if detectIndices.size == 0: continue

        ## Pick the largest event in this (sub)window
        subWinArgmax = detectIndices[np.argmax(detectMetric[detectIndices])]
        trackDetectIdx = winSliceStart+metricOffset+subWinArgmax+maxOffset
        magnitude = detectMetric[subWinArgmax]

        # argmax = (subWinArgmax-tau*mt.trackrate(iqTrack) +
        #           np.argmax(detectMetric[subWinArgmax-tau*mt.trackrate(iqTrack):subWinArgmax+tau*mt.trackrate(iqTrack)]))
        # if argmax > subWinArgmax:
        #   argmin = (subWinArgmax-tau*mt.trackrate(iqTrack) +
        #             np.argmin(detectMetric[subWinArgmax-tau*mt.trackrate(iqTrack):subWinArgmax]))
        # else:
        #   argmin = (subWinArgmax + 
        #             np.argmin(detectMetric[subWinArgmax:subWinArgmax+tau*mt.trackrate(iqTrack)]))

        # eventZC = winSliceStart+metricOffset+argmax+maxOffset
        
        # print "\r",tau,magnitude,tracktimefromslice(iqTrack,trackDetectIdx),tracktimefromslice(iqTrack,eventZC)

        # slicestart = subWinOffset
        # slicestop = len(detectMetric)+slicestart
        # phaseAxes.plot(iqTimes[slicestart:slicestop],detrendedPhase[slicestart:slicestop])
        # metricAxes.plot(iqTimes[slicestart:slicestop]+1.5*tau,detectMetric,
        #                 label='tau=%.3f, estimate=%.2f'%(tau,phasefromdetectionmetric(magnitude,
        #                                                                               tau,carrierRate)))
        # metricAxes.plot(iqTimes[slicestart+argmin]+1.5*tau,detectMetric[argmin],'o')
        
        # fig.show()
        # s = raw_input("Press <enter> to find next event, or 'q' to quit: ")
        # if s == 'q':
        #   exit(0)

        # Check if disqualified or is the same event as one in last window
        if ( (iqDisqualified[trackDetectIdx]!=0) or
             (lastMaxDetect and (trackDetectIdx == lastMaxDetect))):
          continue

        if(magnitude > maxDetection[0]):
          maxDetection = (magnitude,tau,metricOffset,subWinArgmax,detectMetric)

      if maxDetection == (0,0,0,0,0): continue # no events found in this window

      # This was part of an attempt to get a more accurate measure of
      # phase jump magnitude but is purely heuristic
      # # # TODO: Make sure case where event is close to end of track is working
      # try:
      #   # Define the 1.5 tau (or less if this is near end of window) windows initially
      #   # adjacent to the max peak
      #   leftWinStart = max(argmax-metricOffset,0)
      #   rightWinStop = min(argmax+metricOffset,len(detectMetric))

      #   # Find the minimums of the 1.5tau windows to the left and right of the peak
      #   if leftWinStart != argmax:
      #     leftArgmin = leftWinStart + np.argmin(detectMetric[leftWinStart:argmax])
      #   rightArgmin = (argmax
      #                  + np.argmin(detectMetric[argmax:rightWinStop]))

      #   # Now redefine the 1.5 tau windows relative to the argmins
      #   leftWinStart = max(leftArgmin-metricOffset,0)
      #   rightWinStop = min(rightArgmin+metricOffset,len(detectMetric))

      #   # Find the max peak in the windows to the left and right
      #   # The 2nd max peak is the largest of these two
      #   leftArgmax = leftWinStart + np.argmax(detectMetric[leftWinStart:leftArgmin])
      #   rightArgmax = rightArgmin + np.argmax(detectMetric[rightArgmin:rightWinStop])
      #   if (detectMetric[rightArgmax] > detectMetric[leftArgmax]):
      #     (argmax2,eventCenter) = (rightArgmax,rightArgmin)
      #   else:
      #     (argmax2,eventCenter) = (leftArgmax,leftArgmin)
      # except : # sometimes get overindexing from events near end of window
      #   continue

      # # time of first occurring peak relative to this window
      # tPeak1 = min(argmax,argmax2)+subWinOffset+metricOffset
      # # time of 2nd occurring peak
      # tPeak2 = max(argmax,argmax2)+subWinOffset+metricOffset 
      # eventDuration = abs(tPeak1-tPeak2)/trackrate(iqTrack) # in seconds

      # # # Raw phase jump estimate
      # phaseDiff = detrendedPhase[tPeak2] - detrendedPhase[tPeak1]
      # medPhaseDiff = smoothedPhase[tPeak2] - smoothedPhase[tPeak1]

      magnitude,tau,metricOffset,subWinArgmax,detectMetric = maxDetection
      subWinOffset = maxOffset - metricOffset
                      # Start of window+maxOffset = actual time where window starts
                      # + metricOffset time aligns the detect metric with phase
                      # + subWinArgmax is the index relative to the detection metric
                      # = index in track where event occurs
      trackDetectIdx = winSliceStart+maxOffset+metricOffset+subWinArgmax
      lastMaxDetect = trackDetectIdx
      metricPhaseDiff = mt.phasefromdetectionmetric(detectMetric[subWinArgmax],tau,carrierRate)

      # print "\r***Metric hits max of %.3f, with phase estimate: %.3f degrees" % (detectMetric[subWinArgmax],
      #                                                                            metricPhaseDiff)
      
      # Add this event to the list
      eventDetects[eventCnt] = (trackDetectIdx,tau,metricPhaseDiff)
      eventCnt+=1
      if (eventCnt%EVENT_PREALLOC_SIZE) == 0: # expand the list
        eventDetects += [None]*EVENT_PREALLOC_SIZE 

      if verbose:
        detectTime = mt.tracktimefromslice(iqTrack,trackDetectIdx)
        print "\r%d,%s,%02u,%.6f,%.2f,%s" % (eventCnt,mt.strfgpstimestamp(detectTime),iqTrack.attrs['prnCode'],tau,
                                             metricPhaseDiff,args[0])
      if not plot: continue

      # Plot the phase and detection metric of event
      fig = plt.figure()
      phaseAxes = fig.add_subplot(2,1,1); metricAxes = fig.add_subplot(2,1,2,sharex=phaseAxes)
      fig.suptitle("PRN %02u, File: %s" % (iqTrack.attrs['prnCode'],os.path.basename(iqTrack.file.filename)),
                   fontsize=14,fontweight='bold')

      slicestart = subWinOffset      
      slicestop = len(detectMetric)+slicestart
      # times = FitData['iqTimes'] - (subWinArgmax+subWinOffset)/mt.trackrate(iqTrack)
      times = iqTimes - tracktimefromslice(iqTrack,trackDetectIdx)
      phaseAxes.plot(times[slicestart:slicestop],detrendedPhase[slicestart:slicestop],label='Detrended phase',alpha=0.9)
      phaseAxes.plot(times[slicestart:slicestop],movingaverage(detrendedPhase,tau*mt.trackrate(iqTrack)),
                     label='Tau (%.3fs) Moving average'%tau,linewidth=2,color="g",alpha=0.8)
      # phaseAxes.plot(times[slicestart:slicestop],detrendedPhase[slicestart:slicestop]
      #                ,label='Detrended phase',alpha=0.9)
      # phaseAxes.plot(times[slicestart:slicestop],movingaverage(detrendedPhase[slicestart:slicestop],metricOffset)
      #                ,label='Tau Moving average',linewidth=2,color="g",alpha=0.8)
      phaseAxes.axvspan(times[slicestart:slicestop][subWinArgmax-tau*mt.trackrate(iqTrack)],
                        times[slicestart:slicestop][subWinArgmax+tau*mt.trackrate(iqTrack)],
                        facecolor='r',alpha=0.15)
      # plot heuristically found start and end of phase jump
      # phaseAxes.plot(times[tPeak1],smoothedPhase[tPeak1],
      #                'mD',label='Event Start',ms=8,color="#9900FF")
      # phaseAxes.plot(times[tPeak2],smoothedPhase[tPeak2]
      #                ,'mD',label='Event Stop',ms=8,color="#9900FF")

      metricAxes.plot(times[slicestart:slicestop],detectMetric)

      phaseAxes.set_ylabel('Detrended Phase (deg)')
      metricAxes.set_ylabel('Approximate Detection Metric')
      # metricAxes.set_xlabel('Time (s)\nRelative to %s' % strfromgpstime(tracktimefromslice(iqTrack,trackDetectIdx)))
      metricAxes.set_xlabel('Time (s)\nRelative to %s' % strfromgpstimestamp(tracktimefromslice(iqTrack,winSliceStart)))
      #metricAxes.legend(numpoints=1)
      phaseAxes.legend(numpoints=1,bbox_to_anchor=(0., 1.0, 1., .102),ncol=8,loc='lower left',mode='expand')
      # phaseAxes.annotate("Metric based Estimate: %.2fdeg"%(metricPhaseDiff),xy=(10,-50),xycoords='axes points')
      metricAxes.annotate("Max Detection with Tau = %.3f seconds"%(tau),xy=(10,-50),xycoords='axes points')
                         
      fig.show()
      s = raw_input("Press <enter> to find next event, or 'q' to quit: ")
      if s == 'q':
        exit(0)
      plt.close()
#        if s == 's':
#          break
#      break

  if verbose:
    print "###--- Track Stats"
    print '%d Total frames in IQ track' % len(iqTrack)
    print '%d Total frames disqualified' % np.sum(iqDisqualified != 0)
    print '%d frames near ends of track disqualified' \
          % np.sum( (iqDisqualified & (DQ_FILEEND|DQ_FILESTART)) != 0 )
    print '%d Low SNR or near low SNR frames disqualified' \
          % np.sum( (iqDisqualified & (DQ_LOWSNR|DQ_NEARLOWSNR)) != 0 )
    print '%d Loss of code lock frames, or near loss of lock frames disqualified' \
          % np.sum( (iqDisqualified & (DQ_DROPPEDFRAMES|DQ_NEARDROPPEDFRAMES)) != 0 )
    print   "###---\n"

  return eventDetects[:eventCnt]


parser =OptionParser(usage="./mdh_findevents.py [options] MDHfilename")
parser.add_option("--satTime",dest="satTime", action='store_true',default=False,
                  help="Use new time format where IQ_SUMS is in apparent sattelite time.")
parser.add_option("-t", "--time",dest="time", default=None, metavar="TIME",
                  help="Look over a window centered at a given time")
parser.add_option("--snrthreshold", dest='snrthreshold', type='float',default=40,
                  help="Set the SNR threshold in dB below which events will be rejected. (default=40dB)")
parser.add_option("--plot", dest='plot', action='store_true', default=False, 
                  help='Plot the detection metric with detrended phase for each event')
parser.add_option("--log", dest='log', type='string', default=None, 
                  help='File to log events into')
parser.add_option("--threshold", dest='threshold', type='float', default=1.0, 
                  help='Detection metric threshold value (default=1)')
parser.add_option("-c","--code", dest="rangingCode", default='L1:CA', 
                  help="Plot data for carrier:code combination(s). Pass multiple signals in with ','.Ex: -c L1:CA,L2:CM. (default=L1:CA)")
parser.add_option("-p", "--prn", dest="prnCode", default=None, 
                  help="Plot data for PRN(s). Pass multiple PRNs in with ','.Ex: -p 25,04.")
parser.add_option("-q", dest="quiet", action='store_true',default=False, 
                  help="Suppress printing any status messages to output")
options,args = parser.parse_args()

if not args:            sys.exit("Need an MDH file")
options.snrthreshold *= 100 # convert to 0.01 dB units
options.rangingCode =   [x.split(':')[1] for x in options.rangingCode.split(',')]
if options.prnCode:     prnList = [int(x) for x in options.prnCode.split(",")]
if options.time:        trackTime = gpstimefromstr(options.time)

headerInfo = {} # Stores track info for logging
prnDict = {} # Stores the events for each PRN:ranging code combo

filepath = os.path.abspath(args[0])
mdhfile = mdh.MDHFile(filepath,'r')

for iqTrack in mdhfile.matchingDatasets(SUBCLASS='IQ_SUMS'):
  if options.prnCode and iqTrack.attrs['prnCode'] not in prnList:
    continue
  if options.rangingCode and mdh.getEnumAttributeString(iqTrack,'rangingCode') not in options.rangingCode:
    continue
  trackStart = mt.tracktimefromslice(iqTrack,0)
  trackStop = mt.tracktimefromslice(iqTrack,len(iqTrack))

  if (trackStop-trackStart < MINFILESIZE_SECONDS):
    if not options.quiet:
      print "Skipping Track \'%s\' (too short)..." % (iqTrack.name)
    continue

  if options.time:
    startTime = max(trackTime - 40,trackStart)
    stopTime = min(trackTime + 40,trackStop)

    if startTime > stopTime : continue
    if (trackStart > stopTime or trackStop < startTime): continue
  else:
    startTime,stopTime = trackStart,trackStop

  # Add this track to the list
  signalCode = "%02u:%s:%s" % (iqTrack.attrs['prnCode'],mdh.getEnumAttributeString(iqTrack,'carrierCode'),mdh.getEnumAttributeString(iqTrack,'rangingCode'))
  if signalCode not in headerInfo:
    headerInfo[signalCode] = [(trackStart,trackStop)]
  else:
    headerInfo[signalCode].append((trackStart,trackStop))

  if not options.quiet:
    print "Searching PRN%s track starting at %s..." % (signalCode,mt.strfgpstimestamp(trackStart))

  trackslice = ((slicefromtracktime(iqTrack,startTime),slicefromtracktime(iqTrack,stopTime))
                if options.time else None)

  eventDetects = FindDetects(iqTrack,tslice=trackslice,snrthreshold=options.snrthreshold,threshold=options.threshold,         verbose=not options.quiet,plot=options.plot)

  if eventDetects:
    if signalCode not in prnDict:
      prnDict[signalCode] = [(mt.tracktimefromslice(iqTrack,idx),tau,mag) for idx,tau,mag in eventDetects]
    else:
      prnDict[signalCode] += [(mt.tracktimefromslice(iqTrack,idx),tau,mag) for idx,tau,mag in eventDetects]

## Print the events log
if options.log:
  log = open(options.log,'w')
  print >> log, "# ---Search Settings--- #"
  print >> log, "# PRN(s): %s" % options.prnCode
  print >> log, "# Frequency Band:Carrier Code(s): %s" % options.rangingCode

  ## Print track info for each signal searched
  for signalCode in headerInfo:
    print >> log, "# --- PRN %s Ranging Code %s--- #" % (signalCode.split(':')[0],
                                                         signalCode.split(':')[1])
    print >> log, "# Start,Stop Times: "
    for (t1,t2) in headerInfo[signalCode]:
      print >> log, "# Starts at:%s, Ends at: %s" % (mt.strfgpstimestamp(t1),mt.strfgpstimestamp(t2))
    sumTime = sum([(t2-t1)/60. for (t1,t2) in headerInfo[signalCode]])
    print >> log, "# Tracks Total Sum Time in minutes: %.1f" % sumTime
    print >> log, "# Total Events Found in Tracks: %d" % (len(prnDict[signalCode]) if signalCode in prnDict else 0)
#        print >> log, "# Estimated Hourly Event Rate: %d" % len(prnDict[signalCode])/(sumTime*60.)
  print >> log, "# --- --- #"

  ## Finally, print out info for each event
  eventNum = 1 # arbitrarily start at 1
  for signalCode in prnDict:
    prn,cc,rc = signalCode.split(':')
    for (eventTime,tau,eventMag) in prnDict[signalCode]:

      # don't have access to SVN here so put a '-' for compatibility with other scripts
      print >> log, "%d,%s,%s,%s,%s,%s,%f,%.2f,%s" % (eventNum,mt.strfgpstimestamp(eventTime)
                                                      ,'-',prn,cc,rc,tau,eventMag,filepath)
      eventNum += 1
  log.close()
