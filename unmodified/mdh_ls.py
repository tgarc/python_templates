#!/usr/bin/env python

""" Simple program to examine MDH files specified on the command line and list the tracks contained within """

GPS_EPOCH_IN_UNIXTIME=315964800

import sys, os, os.path, h5py, time
from host import mdh
from argparse import ArgumentParser

parser = ArgumentParser(description="Prints attributes on groups/datasets in MDH files.")
parser.add_argument("-v","--verbose", dest="verbose",action='store_true',default=False,
                    help="Print any errors that may have occurred")
parser.add_argument("-R", dest="recursive", action='store_true', default=False, 
                    help="Look recursively inside any specified directories")
parser.add_argument("-i", dest="ignoreothers", action='store_true', default=False, 
                    help="Ignore any non-MDH files")
parser.add_argument("--timefmt", dest="timefmt", default="%Y-%m-%d %H:%M:%S",
                    help="Set time format string")
parser.add_argument("--fmt",dest="fmt",
                    default="%(startTime)s %(stopTime)s %(duration)5u %(prnCode)02u "
                    "%(carrierCode)s %(rangingCode)s %(SUBCLASS)s %(fileName)s", 
                    help='Set format string. Default is \"%(default)s\".')
parser.add_argument("--subclass",dest="subclasslist",default=[],action='append',
                    help="Only report subclass of type (may be specified multiply)")
parser.add_argument(dest="files", nargs="+", metavar="fn",
                    help="MDH files to scan. Directories will be be ignored unless -R is specified.")


args=parser.parse_args()

filesExamined=set()

for fn in args.files:

    if args.recursive and os.path.isdir(fn):
        fileNamesToExamine=[os.path.join(dirpath,filename)
                            for dirpath,dirnames, filenames in os.walk(fn, followlinks=True)
                            for filename in filenames]
    else:
        fileNamesToExamine=[fn]

    for fileName in fileNamesToExamine:
        # Skip if we have already looked at this file?
        absFileName=os.path.abspath(fileName)
        if absFileName in filesExamined: continue
        filesExamined.add(absFileName)

        # See if we have any reason to skip this file
        if not os.path.isfile(fileName): continue
        if args.ignoreothers and not h5py.is_hdf5(fileName): continue

        # Okay, now examine the file.
        try:
            f=mdh.MDHFile(fileName,'r')
        except:
            if args.verbose: print "# Could not open file:", fileName
            continue

        for track in f.matchingDatasets():
            if args.ignoreothers and 'startTime' not in track.attrs: continue
            if args.ignoreothers and 'cadence' not in track.attrs: continue
            if args.ignoreothers and 'timeDenominator' not in track.attrs: continue
            if args.ignoreothers and 'SUBCLASS' not in track.attrs: continue
            if args.subclasslist and track.attrs['SUBCLASS'] not in args.subclasslist: continue

            startTime=track.attrs['startTime']/float(track.attrs['timeDenominator'])
            stopTime=(track.attrs['startTime']+track.attrs['cadence']*len(track))/float(track.attrs['timeDenominator'])

            d={}
            for k,v in track.attrs.iteritems():
                  t=h5py.h5a.open(track.id,k).get_type()
                  if isinstance(t,h5py.h5t.TypeAtomicID):
                      d[k]=v
                  elif isinstance(t,h5py.h5t.TypeStringID):
                      d[k]=v
                  elif isinstance(t,h5py.h5t.TypeEnumID):
                      d[k]=mdh.getEnumAttributeString(track,k)

            d.update({
              'startTime': time.strftime(args.timefmt,time.gmtime(startTime+GPS_EPOCH_IN_UNIXTIME)),
              'stopTime':  time.strftime(args.timefmt,time.gmtime(stopTime+GPS_EPOCH_IN_UNIXTIME)),
              'duration': stopTime-startTime,
              'absFileName':  absFileName,
              'fileName':  fileName,
              'cadence': track.attrs['cadence']/float(track.attrs['timeDenominator'])
              })

            if track.attrs['SUBCLASS']=='IQ_SUMS' and len(track)>0:
              d['numLags']=len(track[0])
            print args.fmt % d

        f.close()
