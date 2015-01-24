"""
Helper tools for getting data from the HRTR events database. (mostly used by dbplotter.py)

TODO
+Make sure that all queries will give the correct results if there are multiple GNSS receiver locations (i.e. more than one uniq_receiver_id)

Note: All functions which are independent of the table structure should be placed at the top of the file.
All other functions should go below.  This will help separate the two classes.


Thomas J Garcia and David Munton
"""

import sys
import psycopg2 as ps
import psycopg2.extras
import dateutil as du
import time
import mdh_tools as mtools
from operator import itemgetter,attrgetter
from itertools import groupby

VERBOSE = 0

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
	
def executeQuery(conn,qstr,qDict=None):
  """Function that allows a connection to the events, or test_events, database.
     Arguments:
        conn - a psycopg2 connection object
        qstr - a valid query
        qDict - a dictionary of named variable values to go with the cursor.  This should
        always be used to avoid sql injection attacks
  """
  if VERBOSE: sys.stdout.write("Querying..."); sys.stdout.flush()
  if VERBOSE > 1: sys.stdout.write("QUERY:\n%s\n" % qstr); sys.stdout.flush()
  ## print type(qDict)

  try:
    queryCursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    t1 = time.time()
    if qDict == None:
      queryCursor.execute(qstr)
    elif isinstance(qDict, dict):
      queryCursor.execute(qstr, qDict)
    t2 = time.time()
  except ps.Error, e:
    if VERBOSE: print "Query error:", e
    queryCursor.close()
    queryCursor = None
  else:
    if VERBOSE: sys.stdout.write("%d rows returned in %g s.\n" % (queryCursor.rowcount,t2-t1))
  
  return queryCursor


###############################################################################################
#  Database dependent functions follow this line
##############################################################################################

#  First a bunch of standard queries.  We should look at alternate ways to do this.

# filter out events triggered by ephemeris cutout
ephem_cutout_filter = '''\
SELECT
  *
FROM
(
%(query)s
) AS t
WHERE
(
t.event_date <= '2013-10-09 00:00:00'
AND
    date_part('hour', t.event_date )::integer %% 2 != 0
AND NOT
    date_part('second', t.event_date) BETWEEN 18 and 22
)
OR
    t.event_date > '2013-10-09 00:00:00' '''


raw_events_query = '''\
SELECT
    e.event_id
    ,e.connection_id
    ,e.uniq_receiver_id
    ,e.signal_id
    ,e.demodulator_no
    ,e.event_type
    ,d.description
    ,CASE 
    WHEN dhb.detector_version_major='0' AND dhb.detector_version_minor<'2' THEN 
        e.event_date-(1.5*pe.detection_period || ' second')::interval 
    ELSE 
        e.event_date 
    END AS event_date 
    ,CASE 
    WHEN dhb.detector_version_major='0' AND dhb.detector_version_minor<'2' THEN 
        e.event_date-(2.5*pe.detection_period || ' second')::interval 
    ELSE 
        e.event_date-(pe.detection_period || ' second')::interval 
    END AS event_start 
    ,CASE 
    WHEN dhb.detector_version_major='0' AND dhb.detector_version_minor<'2' THEN 
        e.event_date-(0.5*pe.detection_period || ' second')::interval 
    ELSE 
        e.event_date+(pe.detection_period || ' second')::interval 
    END AS event_end 
    ,max_detection_value
    ,detection_period
    ,threshold
    ,prn
    ,sv_type
    ,cast(sv_identifier as int) AS sv_identifier
    ,trim(ranging_code) AS ranging_code
    ,carrier_code
    ,max_signal_strength
    ,min_signal_strength
FROM          
    events AS e
JOIN
    phase_events AS pe
USING(event_id)
JOIN
    event_status_description AS d
USING(event_status)
JOIN 
    gps_satellite_signals AS gss
ON
    e.signal_id=gss.signal_id
AND
    e.event_date >= gss.valid_start
AND
    (gss.valid_end IS NULL OR e.event_date <= gss.valid_end)
JOIN 
    detector_heartbeat AS dhb 
ON 
    e.event_type = dhb.event_type
AND
    e.connection_id = dhb.connection_id 
AND 
    e.uniq_receiver_id = dhb.uniq_receiver_id 
AND 
    e.signal_id = dhb.signal_id 
AND 
    e.event_date BETWEEN dhb.start_date and dhb.end_date 
%(whereclause)s'''


cluster_events_query = '''\
SELECT
    *
FROM
(
SELECT
    s3.*
    ,max(max_detection_value)
    OVER
    (PARTITION BY left_edge,s3.signal_id,s3.uniq_receiver_id,s3.connection_id)
    AS max_cluster_detection_value
FROM
(
SELECT
    s2.*
    ,max(new_start)
    OVER (PARTITION BY s2.signal_id,s2.uniq_receiver_id,s2.connection_id ORDER BY event_start,event_end)
    AS left_edge
FROM
(
SELECT
    s1.*
    ,CASE
    WHEN event_start < max(le)
    OVER (PARTITION BY s1.signal_id,s1.uniq_receiver_id,s1.connection_id ORDER BY event_start,event_end)
    THEN
	NULL
    ELSE
	event_start
    END
    AS new_start
FROM
(
SELECT
    e.*
    ,lag(event_end)
    OVER
    (PARTITION BY e.signal_id,e.uniq_receiver_id,e.connection_id ORDER BY event_start,event_end)
    AS le
FROM
(\n'''+raw_events_query+'''\n)
AS e
) AS s1
) AS s2
) AS s3
) AS s4
WHERE
  max_detection_value = max_cluster_detection_value'''


events_dataset_query = '''\
WITH e AS (\n%(query)s\n)
SELECT DISTINCT ON (event_id)
    e.*
    ,rds.file_path
    ,rds.data_set AS data_set_iq
FROM
    e
LEFT JOIN
    raw_data_sets AS rds
ON
    rds.signal_id=e.signal_id
AND
    rds.uniq_receiver_id=e.uniq_receiver_id
AND
    rds.connection_id=e.connection_id
AND
    e.event_date BETWEEN rds.start_date AND rds.end_date
WHERE
    date_part('second', LEAST(e.event_date-rds.start_date,rds.end_date-e.event_date)) > 1.5*e.detection_period
AND
    data_set_type = '403'
ORDER BY 
    event_id'''


# Query for plotting heartbeat message presence
heartbeat_query = '''\
SELECT 
    GREATEST(dhb.start_date,'%(startdate)s') as start_date
    ,LEAST(dhb.end_date,'%(enddate)s') as end_date
    ,cast(sv_identifier as int) 
    ,sv_type 
    ,prn 
    ,carrier_code 
    ,trim(ranging_code) as ranging_code 
FROM
    detector_heartbeat as dhb
JOIN
    gps_satellite_signals as gss
ON
    dhb.signal_id = gss.signal_id
AND
    start_date >= valid_start
AND
    (valid_end IS NULL OR end_date <= valid_end)
%(whereclause)s
ORDER BY
    prn,carrier_code,start_date,end_date'''


def executeEventsQuery(conn, form, clusterEvents=True):
    # don't cluster when looking for a particular event id
    if 'event_id' in form: clusterEvents = False 

    wc = ("WHERE\n    %s" % buildStdWhereClause(form)) if form else ''

    # only consider the events with the clustered flag set true for now
    wc += "\nAND\n    " if wc else "WHERE\n    "
    wc += "clustered = false"
    
    eventsform = {'whereclause':wc}
    if 'startdate' in form: eventsform['startdate'] = form['startdate']
    if 'enddate' in form: eventsform['enddate'] = form['enddate']
                  
    eventsquery = (cluster_events_query
                   if clusterEvents else raw_events_query) % eventsform
    qstr = events_dataset_query % {'query':eventsquery}

    if ('startdate' in form and du.parser.parse(form['startdate']) < du.parser.parse('2013-10-09 00:00:00')):
        qstr = ephem_cutout_filter % {'query':qstr}

    curs = executeQuery(conn,qstr)
    if not curs:
        return None,qstr,wc
    qlist = curs.fetchall()
    curs.close()

    return qlist,qstr,wc


def executeHeartbeatsQuery(conn,form):
    wc = ("WHERE\n    %s" % buildStdWhereClause(form,detector_heartbeat=True)) if form else ''

    hbform = {}
    hbform['enddate'] = form['enddate'] if 'enddate' in form else datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    hbform['startdate'] = form['startdate'] if 'startdate' in form else '2013-06-01'
    hbform['whereclause'] = wc

    qstr = heartbeat_query % hbform

    curs = executeQuery(conn,qstr)
    if not curs:
        return None,qstr,wc
    hblist = curs.fetchall()
    curs.close()

    # merge overlapping intervals
    svkey = 'sv_identifier', 'sv_type', 'prn','carrier_code','ranging_code'
    carrierCode = itemgetter(*svkey)
    hblist.sort(key=carrierCode)

    mergedHBlist = []
    cumulativerow = dict(zip(svkey,[None]*len(svkey)))
    for keys, rows in groupby(hblist, carrierCode):
      cumulativerow.update(zip(svkey,keys))
      for start,end in mtools.mergeIntervals(list(rows),'start_date','end_date'):
        cumulativerow.update( {'start_date': start, 'end_date': end} )
        mergedHBlist.append(cumulativerow.copy())

    return mergedHBlist,qstr,wc


def buildStdWhereClause(form,detector_heartbeat=False):
    """
    Builds a standard where clause given the cgi form that is
    sent.Returns a string where clause compatible with the database
    interface.
    """
    clauseList = []
    lb = '('
    rb = ')'
   
    # This is kind of a hack and currently not well defined
    # specifically, fields in the form should be associated with a particular table since
    # they may exist in multiple tables.
    if detector_heartbeat:
        if "startdate" in form and "enddate" in form:
            clauseList.append("(start_date,end_date) OVERLAPS ('%s'::timestamp,'%s'::timestamp)"
                              % (form['startdate'],form['enddate']))
        else:
            return ''
        if "event_type" in form: # needed here since event_type exists in both heartbeat and events table
            clauseList.append('event_type in '
                              + lb + ','.join(["'%s'" % p for p in form["event_type"]]) + rb)
    else:
        if "startdate" in form : 
            clauseList.append("event_date > '%s'" % form["startdate"])

        if "enddate" in form: 
            clauseList.append("event_date < '%s'" % form["enddate"])

        if "event_id" in form:
            clauseList.append('event_id in '
                              + lb + ','.join(["'%s'" % p for p in form["event_id"]]) + rb)

        if "max_detection_value" in form:
            clauseList.append("max_detection_value < %f" % form["max_detection_value"])

        if "max_threshold" in form:
            clauseList.append("threshold < %f" % form["max_threshold"])

        if "event_type" in form:
            clauseList.append('e.event_type in '
                              + lb + ','.join(["'%s'" % p for p in form["event_type"]]) + rb)
        if "event_status" in form:
            clauseList.append('description in '
                              + lb + ','.join(["'%s'" % p for p in form["event_status"]]) + rb)

    if "systems" in form:
        clauseList.append('system_code in '
                          + lb + ','.join(["'%s'" % p for p in form["systems"]]) + rb)

    if "svn" in form: 
        clauseList.append('sv_identifier in '
                          + lb + ','.join(["'%s'" % p for p in form["svn"]])  + rb)

    if "prn" in form: 
        clauseList.append('prn in '
                          + lb + ','.join(["'%s'" % p for p in form["prn"]])  + rb)

    if "ranging_code" and "carrier_code" in form:
        clauseList.append(lb + '\nOR\n    '.join("(ranging_code='%s' AND carrier_code='%s')" % (rc,cc)
                          for cc,rc in zip(form['carrier_code'],form['ranging_code'])) + rb)
            
    elif "ranging_code" in form: 
        clauseList.append('ranging_code in '
                          + lb + ','.join(["'%s'" % p for p in form["ranging_code"]])  + rb)

    elif "carrier_code" in form: 
        clauseList.append('carrier_code in '
                          + lb + ','.join(["'%s'" % p for p in form["carrier_code"]])  + rb)

    if "sv_type" in form:
        clauseList.append('sv_type in '
                          + lb + ','.join(["'%s'" % p for p in form["sv_type"]])  + rb)

    if "min_signal_strength" in form:
        clauseList.append("min_signal_strength > %f" % form["min_signal_strength"])

    if len(clauseList) > 0:
        whereclause = '\nAND\n    '.join(clauseList)
    else:
        whereclause =''

    return whereclause


def getEventRefDataset(conn,eid):
  curs = executeQuery(conn,'''\
SELECT
  *
FROM
  get_event_ref_data_set_iq(%d)
;''' % eid)

  if curs:
    eventds = curs.fetchone()
    curs.close()

  return eventds


def getEventDataset(conn,eid):
  curs = executeQuery(conn,'''\
SELECT
  *
FROM
  get_event_data_set_iq(%d)
;''' % eid)
  ds = curs.fetchone()
  curs.close()

  return ds


def getTableKeys(conn,dbtable,schema='hrtr_monitor'):
    qstr = '''\
    SELECT
        column_name
    FROM
        information_schema.columns
    WHERE
        table_schema='%(schema)s'
    AND
        table_name = '%(name)'
    ''' % {'schema':schema, 'name':dbtable}

    try:
        curs = executeQuery(conn,qstr)
        rows = curs.fetchall()
        curs.close()
    except Exception, e:
        print '%s'%e
        print 'An error in getTableKeys has occurred.'
        return []
    return rows


def getSystemList(myConn, systems='ALL'):
  """Function which returns a list of all GNSS systems in the database meeting the
     keyword restriction.
     Arguements:
          myConn - a psycopg2 connection created with connectToDB
     Keywords:
          systems = "{string}" default: "ALL"
  """
  try:
    # Get the query to execute
    query = getSystemListQuery(system=systems)
    myCur = executeQuery(myConn,query)
    rows = myCur.fetchall()
    myCur.close()
  except Exception, e:
    print '%s'%e
    print 'An error in getSystemList has occurred.'
    rows = []
  return rows


def getGPSPRNSVList(myConn):
  """Function which returns a list of all PRN codes in the database meeting.
     Arguments:
          myConn - a psycopg2 connection created with connectToDB
  """
  try:
    distinctSVPRN = getGPSSVNPRNListQuery()
    myCur = executeQuery(myConn,distinctSVPRN)
    rows = myCur.fetchall()
    myCur.close()
  except Exception, e:
    print 'An error in getGPSSystemList has occurred.'
    print '%s'%e
    rows = []
  return rows


def getGPSEventData(myConn, whereclause = '', groupclause ='', orderclause='event_date'):
  """
    Arguments:
      myConn - a valid psycopg2 connection
    Keywords:
      whereclause - A valid SQL query.  Default is to select all records.
      groupclause - A valid SQL group by clause.  No default.
      orderClause - A valid order by statement. The default is to order on event_date.
    Returns:
      List of rows returned,by default ordered by date.
  """    
  # Set up the database connection information
  colList = []
  try:
      # added the distinct on clause to eliminate events with duplicate instances
      query = "select * from hrtr_monitor.gps_phase_events "
      # print 'Query 1:', query
    
      if whereclause != "":
        query = "%s where %s" % (query,whereclause)
      if groupclause != '':
        query = "%s group by %s" % (query,groupclause)
      if orderclause != '':
        query = "%s order by %s" % (query,orderclause)
        
      query = query +';'
      # print 'Query 2:', query
      myCur = executeQuery(myConn, query)
      if myCur is not None:
        colList = [desc[0] for desc in myCur.description]
        data = myCur.fetchall()
        myCur.close()
      else:
        print 'No cursor returned.'
        print 'Query used is: ', query
        colList = []
        data = []
      
      
  except Exception, e:
      print 'Error in getGPSEventData'
      print 'Query used is: ', query
      print '%s'%e

  return (colList, data)


def getGPSSVNPRNListQuery():
  query = 'Select distinct sv_identifier, prn from all_gps_satellite_signals order by sv_identifier;'
  return query


# returns SQL query to get list of systems
def getSystemListQuery(system = 'ALL'):
  if system=='ALL':
    query="Select system_code, description from hrtr_monitor.systems;"
  else:
    query="Select system_code, description from hrtr_monitor.systems where system_code = '%s';"%(system)

  return query

  
def getEventStatusList(myConn):
  """ Function to get a list of available event status meanings
      Arguments:
	      myConn - a psyopg2 database connection
	  Returns:
	      A dictionary of event status definitions
  """
  statusDict = {}
  try:
      query = "select event_status, description from hrtr_monitor.event_status_description;"
      myCur = executeQuery(myConn, query)
      for row in myCur:
        statusDict[row['event_status']] = row['description']
      myCur.close()
      
  except Exception, e:
      print 'Error in getEventStatusList'
      print 'Query used is: ', query
      print '%s'%e

  return statusDict


def getEventStatus(myConn, eventIds):
  """ Function to return the status of one or more events
      Arguments:
	      myConn - a psyopg2 database connection
	      eventIds - either a list of one or more eventIds
	  Returns:
	      A dictionary of event status definitions
  """
  eventStatusDict = {}
  if len(eventIds) > 1:
    queryCondition = 'event_id in (' + ','.join(eventIds) + ')'
  elif len(eventIds) == 1:
    queryCondition = 'event_id=%i'%(eventIds[0])
    
  try:
      query = "select event_id, event_status from hrtr_monitor.events where %s ;"%(queryCondition)
      myCur = executeQuery(myConn, query)
      for row in myCur:
        eventStatusDict[row['event_id']] = row['event_status']
      myCur.close()
      
  except Exception, e:
      print 'Error in getEventStatus'
      print 'Query used is: ', query
      print '%s'%e

  return eventStatusDict


def updateEventStatus(myConn, eventStatusDict):
  """ Function to update the status of one event
      Arguments:
	      myConn - a psyopg2 database connection
	      eventStatusDict - a dictionary of the form {'event_status': val, 'event_id':int}
	  Returns:
	      0 on successful completion
  """
  # This is the naive approach to the update.  A more sophisticated approach would be to use the VALUES clause to create a mapping table.
  query = "update hrtr_monitor.events set event_status = %(event_status)s where event_id = %(event_id)s"
  #print query + '\n'
  
  try:
      myCur = executeQuery(myConn, query, qDict = eventStatusDict)
      # Commit the changes to the database
      myConn.commit()
      myCur.close()
      returnVal = 0
   
  except Exception, e:
      print 'Error in getEventStatus\n'
      print 'Query used is:%s\n'%(query)
      print '%s'%e
      returnVal = -1

  return returnVal


if __name__=='__main__':
  # Values need to be set for test_events or for events.
  # Values for events are the default: 'name':'Production Database','dbHost':"borodin",'dbUsr':'hrtr_monitor','dbPwd':'hrtr_monitor','dbName':"events",'dbPort':5432
  # Values for test_events are: 'name':'Test Database','dbHost':"borodin",'dbUsr':'hrtr_monitor','dbPwd':'hrtr_monitor','dbName':"test_events",'dbPort':5433
  dbLoginDict = {'name':'Test Database','dbHost':"borodin",'dbUsr':'hrtr_monitor','dbPwd':'hrtr_monitor','dbName':"test_events",'dbPort':5433}

  myConn = connectToDB(dbHost=dbLoginDict['dbHost'], dbUsr=dbLoginDict['dbUsr'] ,dbPwd =dbLoginDict['dbPwd'],dbName=dbLoginDict['dbName'], dbPort=dbLoginDict['dbPort'])
  wclause = "system_code in ('GPS') and description in ('Not reviewed') and sv_identifier in ('60') and event_date >'2013-02-13' and event_date < '2014-02-27'"
  myColumns, myData = getGPSEventData(myConn, whereclause = wclause, groupclause ='', orderclause='event_date')

  myConn.close()
