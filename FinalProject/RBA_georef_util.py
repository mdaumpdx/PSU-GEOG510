# **********************************************************************
# 
# NAME: Marilyn Daum
# DATE: 15 Mar 2015
# CLASS: GEOG510
# ASSIGNMENT: Final Project
# 
# DESCRIPTION:  Utility classes and functions for working with RBA survey
# georeferencing files and data adjustment object structures.  This file
# is for import by top-level scripts only.
#
# SOURCE(S): http://resources.arcgis.com/en/help/
#            https://docs.python.org/
#            "Python in a Nutshell" by Alex Martelli
#            http://chimera.labs.oreilly.com/books/1230000000393/ch06.html#_solution_94
#            http://stackoverflow.com/questions/9338507/converting-a-string-with-scientific-notation-to-an-int-in-python
#            http://stackoverflow.com/questions/1535327/python-how-to-print-a-class-or-objects-of-class-using-print
#            http://stackoverflow.com/questions/867115/python-instance-variables-as-optional-arguments
#            http://stackoverflow.com/questions/2937114/python-is-not-sequence
#            http://stackoverflow.com/questions/10564801/how-to-unpack-multiple-tuples-in-function-call
#            http://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
#
# **********************************************************************

# ********** IMPORT STATEMENTS **********
import sys
import os
import logging
import argparse
import csv
from collections import namedtuple
import arcpy


# ********** GLOBAL CONSTANTS **********

LLID = "LocationID"
STREAMNAME = "Stream_Name"
TRIB_TO = "Trib_To"
BEGIN_ = "Begin_"
END_ = "End_"
X_COORD = "X_coord"
Y_COORD = "Y_coord"
XY_NOTE = "XY_Note"
SURVEY_CUM_DISTANCE = "Survey_Cum_Dist"
STREAMLINE_CUM_DISTANCE = "Streamline_Cum_Dist"
SURVEY_COMMENT = "Comment"
ADJ_FACTOR = "Adj_Factor"

STREAMS_FC_NAME = "streams"

DEFAULT_ADJ_FACTOR = 1.0  # use when adjustment factor cannot be computed
DEFAULT_BEGIN_DIST = 0  # min cummulative distance for stream survey data
DEFAULT_END_DIST = 999999  # max cummulative distance for stream survey data

LAT_LONG_CRS = arcpy.SpatialReference(4326)

#LOG_LEVEL = logging.DEBUG
LOG_LEVEL = logging.INFO  # may be overwritten by importing module


# ********** CLASSES **********

class SyncPoint(object):
    """
    x,y and distance information for a synchronization point along a stream.
    SyncPoint can be used at the start of a stream, at an intermediate point
    where x and y coordinates are given, and at the end of a survey/stream.
    """

    def __init__(self,
                 x_coord=None,
                 y_coord=None,
                 xy_note=None,
                 survey_cum_dist=None,
                 streamline_cum_dist=None,
                 survey_comment=None):
        pass

    def __repr__(self):
        return "SyncPoint ({}, {}), survey dist {}, streamline dist {}".\
            format(self.x_coord, self.y_coord, self.survey_cum_dist,
                   self.streamline_cum_dist)

class StreamDistanceInfo(object):
    """
    Distance adjustment factors and related information for a single stream.
    """

    def __init__(self,
                 llid=None,
                 name="",
                 trib_to="",
                 adj_factors=[]):
        pass

    def __repr__(self):
        return "StreamDistanceInfo {}, {} trib to {}, {} adj factors ".\
            format(self.llid, self.name, self.trib_to, len(self.adj_factors))


# ********** FUNCTIONS **********


def valid_gdb(gdb_path):
    """
    Verifies gdb_path is a a path to an existing and valid geodatabase
    containing a polyline feature class.
    :param gdb_path: full path to geodatabase
    :return: verified gdb_path
    An argparse.ArgumentTypeError is raised if gdb_path is not valid.
    """
    try:
        desc = arcpy.Describe(gdb_path)
        assert desc.DataType == 'Workspace'
        assert desc.workspacetype == 'LocalDatabase'
    except IOError:
        raise argparse.ArgumentTypeError\
            ("Geodatabase {} does not exist.".format(gdb_path))
    except AssertionError:
        raise argparse.ArgumentTypeError\
            ("Geodatabase {} is not valid.".format(gdb_path))

    return gdb_path


def valid_file(filepathname):
    """
    Verifies filepathname is a path to an existing file.
    :param filepathname: full path to file
    :return: verified filepathname
    An argparse.ArgumentTypeError is raised if filepathname is not valid.
    """
    if not os.path.exists(filepathname):
        raise argparse.ArgumentTypeError\
            ("Cannot find file {}".format(filepathname))

    return filepathname


def valid_filedir(filepathname):
    """
    Verifies the directory in filepathname is a path to an existing directory.
    :param filepathname: full path to file
    :return: verified filepathname
    An argparse.ArgumentTypeError is raised if filepathname is not valid.
    """
    dirname = os.path.dirname(filepathname)
    if not os.path.exists(dirname):
        raise argparse.ArgumentTypeError\
            ("Cannot find directory {}".format(dirname))

    return filepathname


def valid_gdb_file(gdb_file_path):
    """
    Verifies gdb_file_path is a a path to an existing and valid geodatabase
    file.
    :param gdb_file_path: full path to geodatabase file
    :return: verified gdb_file_path
    An argparse.ArgumentTypeError is raised if gdb_file_path is not valid.
    """
    try:
        desc = arcpy.Describe(gdb_file_path)
        assert desc.DataType == 'FeatureClass'
    except IOError:
        raise argparse.ArgumentTypeError\
            ("Geodatabase file {} does not exist.".format(gdb_file_path))
    except AssertionError:
        raise argparse.ArgumentTypeError\
            ("Geodatabase file {} is not valid.".format(gdb_file_path))

    return gdb_file_path


def get_valid_polyline_pathname(gdb_path, line_fc_name):
    """
    Verifies a polyline feature class exists at location given by
    gdb_path + stream_fc.
    :param gdb_path: full path to a geodatabase
    :param line_fc_name: name of a polyline feature class in the geodatabase
    :return: full path to polyline feature class
    """
    polyline_pathname = os.path.join(gdb_path, line_fc_name)
    try:
        desc = arcpy.Describe(polyline_pathname)
        assert desc.DataType == 'FeatureClass'
        assert desc.featureType == 'Simple'
        assert desc.shapeType == 'Polyline'
    except IOError:
        raise argparse.ArgumentTypeError\
            ("Line feature class {} does not exist.".format(polyline_pathname))
    except AssertionError:
        raise argparse.ArgumentTypeError\
            ("Line feature class {} is not valid.".format(polyline_pathname))

    return polyline_pathname


def get_stream_geom(streams_pathname, llid):
    """
    Finds geometry object for stream matching input LLID.
    :param streams_pathname: feature class containing streams
    :param llid: Location ID for stream
    :return: stream geometry object
    """
    where_clause = """{} = '{}'""".\
        format(arcpy.AddFieldDelimiters(streams_pathname, LLID), llid)
    selected = [sr[0] for sr in arcpy.da.SearchCursor(streams_pathname,
                                                      ["SHAPE@"],
                                                      where_clause)]
    if len(selected) > 1:
        logging.warning((" Multiple matches for stream with {} {}. ".
                         format(LLID, llid) +
                         "Using first match."))

    return selected[0]


def new_sdi_object(llid, streamname, trib_to, adj_factors):
    """
    Create a new StreamDistanceInfo object with the attributes given
    :param llid: location ID
    :param streamname: stream name
    :param trib_to: name of stream into which this stream flows
    :param adj_factors: sequence of tuples with adjustment factors
        for this stream; each tuple contains
        (begin_sync_point, end_sync_point, adj_factor)
    :return: new StreamDistanceInfo object
    """
    sdi_obj = StreamDistanceInfo()
    sdi_obj.llid = llid
    sdi_obj.name = streamname
    sdi_obj.trib_to = trib_to
    sdi_obj.adj_factors = adj_factors
    return sdi_obj


def add_adj_factors_to_sdi(adj_factor_dictionary, llid, adj_factors):
    """
    Adds given adj_factors sequence to the appropriate StreamDistanceInfo
    object in the given adj_factor_dictionary.
    :param adj_factor_dictionary: dictionary of stream distance adjustment
        information, keyed on stream LLID.  Each value contains a
        StreamDistanceInfo object.
    :param llid: stream location id, key for adj_factor_dictionary
    :param adj_factors: sequence of adjustment factors to add to the given
        llid's StreamDistanceInfo object.
    :return: N/A, adj_factor_dictionary and StreamDistanceInfo object are
        modified by this function
    """
    sdi_obj = adj_factor_dictionary[llid]
    sdi_obj.adj_factors = adj_factors


def write_sdi_to_csv_file(sdi_dict, sdi_filepath):
    """
    Writes the contents of the given stream distance info dictionary to the
    given filepath.  One row is written for each adjustment factor, yielding
    possibly multiple rows for the same LLID.
    :param sdi_dict:  dictionary of stream distance adjustment
        information, keyed on stream LLID.  Each value contains a
        StreamDistanceInfo object, with a sequence of adjustment factors.
    :param sdi_filepath: full path to csv file where dictionary contents
        will be written
    :return: N/A, file at sdi_filepath is created and populated
    """
    # Column headers
    fieldnames = [LLID,
                  STREAMNAME,
                  TRIB_TO,
                  BEGIN_+SURVEY_CUM_DISTANCE,
                  BEGIN_+STREAMLINE_CUM_DISTANCE,
                  BEGIN_+X_COORD, BEGIN_+Y_COORD,
                  BEGIN_+XY_NOTE,
                  BEGIN_+SURVEY_COMMENT,
                  END_+SURVEY_CUM_DISTANCE,
                  END_+STREAMLINE_CUM_DISTANCE,
                  END_+X_COORD, END_+Y_COORD,
                  END_+XY_NOTE,
                  END_+SURVEY_COMMENT,
                  ADJ_FACTOR]

    with open(sdi_filepath, 'wb') as csvfile:
        adj_fact_writer = csv.writer(csvfile)
        adj_fact_writer.writerow(fieldnames)
        # Loop through each stream LLID in dictionary
        for stream_id in sorted(sdi_dict):
            sdi_obj = sdi_dict[stream_id]
            # Write a row for each adjustment factor, including
            # begin and end sync points
            for begin_sync_pt, end_sync_pt, adj_factor in sdi_obj.adj_factors:
                logging.debug(" {}".format
                      (chain_data_two_levels(stream_id,
                                             sdi_obj.name,
                                             sdi_obj.trib_to,
                                             expand_sync_pt(begin_sync_pt,
                                                         DEFAULT_BEGIN_DIST),
                                             expand_sync_pt(end_sync_pt,
                                                         DEFAULT_END_DIST),
                                             adj_factor)))
                adj_fact_writer.writerow\
                    (chain_data_two_levels
                     ("'{}'".format(stream_id),  # need to force LLID to be text
                                                 # (vs sci notation)
                      sdi_obj.name,
                      sdi_obj.trib_to,
                      expand_sync_pt(begin_sync_pt, DEFAULT_BEGIN_DIST),
                      expand_sync_pt(end_sync_pt, DEFAULT_END_DIST),
                      adj_factor))


def read_sdi_from_csvfile(sdi_filepath):
    """
    Reads the contents of the given filepath into a stream distance info
    dictionary.
    :param sdi_filepath: full path to csv file containing stream distance
        information.  There is one row for each adjustment factor, so there
         may be multiple rows for the same dictionary key LLID.
    :return: dictionary of stream distance adjustment information,
        keyed on stream LLID.  Each value contains a StreamDistanceInfo object,
        with a sequence of adjustment factors.  Each adjustment factor is a
        tuple: (begining_SycnPoint, ending_SyncPoint, adjustment_factor)
    """
    stream_distance_info_dict = {}
    adj_factors = []
    prev_llid = ""
    with open(sdi_filepath, 'rb') as sdi_file:
        # Read and process each row in csv file as namedtuple
        sdi_file_reader = csv.reader(sdi_file)
        headings = sdi_file_reader.next()
        Row = namedtuple('Row',headings)
        for r in sdi_file_reader:
            row = Row(*r)
            new_llid = str(row.LocationID).replace("'","")
            if new_llid != prev_llid:
                # New Stream
                if prev_llid != "":
                    # Store all adjustment factors for prev_llid
                    add_adj_factors_to_sdi(stream_distance_info_dict,
                                           prev_llid, adj_factors)
                # Reset stream variables and initialize dictionary entry
                # for new llid
                adj_factors = []
                # Gather basic data for new StreamDistanceInfo object
                streamname = str(row.Stream_Name)
                trib_to = str(row.Trib_To)
                # Add basic stream info to dictionary
                stream_distance_info_dict[new_llid] = \
                    new_sdi_object(new_llid, streamname, trib_to, adj_factors)
                prev_llid = new_llid

            # Gather data for one adjustment factor
            begin_sync_pt = create_syncpoint\
                (parse_float_or_NA(row.Begin_X_coord),
                 parse_float_or_NA(row.Begin_Y_coord),
                 str(row.Begin_XY_Note),
                 int(row.Begin_Survey_Cum_Dist),
                 float(row.Begin_Streamline_Cum_Dist),
                 str(row.Begin_Comment))
            end_sync_pt = create_syncpoint\
                (parse_float_or_NA(row.End_X_coord),
                 parse_float_or_NA(row.End_Y_coord),
                 str(row.End_XY_Note),
                 int(row.End_Survey_Cum_Dist),
                 float(row.End_Streamline_Cum_Dist),
                 str(row.End_Comment))
            adj_factor = float(row.Adj_Factor)
            #  Add this adjustment factor to the list for this SDI object
            adj_factors.append((begin_sync_pt, end_sync_pt, adj_factor))

        # End of file, no more rows
        # Tie up processing for last stream distance info
        add_adj_factors_to_sdi(stream_distance_info_dict, prev_llid,
                               adj_factors)

    sdi_file.close()
    return stream_distance_info_dict


def create_syncpoint(in_x_coord, in_y_coord, xy_note,
                    survey_cum_dist, streamline_cum_dist,
                    survey_comment):
    """
    Creates a new SyncPoint object with the given attributes
    :param in_x_coord: x coordinate for sync point
    :param in_y_coord: y coordinate for sync point
    :param xy_note: text field with notes about xy sync (may exist
        regardless of values in other fields)
    :param survey_cum_dist: survey-reported cumulative distance
    :param streamline_cum_dist: cumulative distance calculated along streamline
    :param survey_comment: text fields with notes about survey data/point
    :return: new SyncPoint object with all fields populated
    """
    syncpt = SyncPoint()
    syncpt.x_coord = in_x_coord
    syncpt.y_coord = in_y_coord
    syncpt.xy_note = xy_note
    syncpt.survey_cum_dist = survey_cum_dist
    syncpt.streamline_cum_dist = streamline_cum_dist
    syncpt.survey_comment = survey_comment
    logging.debug(" created syncpt {}".format(syncpt))
    return syncpt


def expand_sync_pt(sync_pt, default_distance):
    """
    Returns a tuple representing the values of the attribues for the given
    SyncPoint object.  If input SyncPoint is None, a minimal tuple of
    empty/default values is returned.
    :param sync_pt: Either a populated SyncPoint object or None.
    :param default_distance: Distance to use for values of
        survey_cum_dist and streamline_cum_dist, when input sync_pt is None.
    :return: tuple of (survey_cum_dist, streamline_cum_dist,
        x_coord, y_coord, xy_note, survey_comment) values for sync_pt
    """
    if sync_pt is None:
        return (default_distance, default_distance, "", "", "", "")
    else:
        return (sync_pt.survey_cum_dist, \
                sync_pt.streamline_cum_dist, \
                sync_pt.x_coord, \
                sync_pt.y_coord, \
                sync_pt.xy_note, \
                sync_pt.survey_comment)

def chain_data_two_levels(*elements):
    """
    Flattens given iterable into a list.  Function dives into
    list and tuple components, down to two levels.
    :param elements: iterable possibly containing nested lists and tuples
    :return: list of flattened components of input elements
    """
    result = []
    for elem in elements:
        if type(elem) in [list, tuple]:
            for el in elem:
                result.append(el)
        else:
            result.append(elem)
    return result

def parse_float_or_NA(read_value):
    """
    Utility to parse possibly-empty input value as a float.
    :param read_value: Input value, either empty or something that can
        be cast to a float.
    :return: Either float value or None.
    """
    if str(read_value) == "":
        return None
    else:
        return float(read_value)


# ********** MAIN **********

def main():
    logging.error(" Not intended for top-level use.")
    return 1


# ********** MAIN CHECK **********

if __name__ == '__main__':
    sys.exit(main())
