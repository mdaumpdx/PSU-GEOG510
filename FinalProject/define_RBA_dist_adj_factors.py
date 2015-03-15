# **********************************************************************
# 
# NAME: Marilyn Daum
# DATE: 15 Mar 2015
# CLASS: GEOG510
# ASSIGNMENT: Final Project
# 
# DESCRIPTION: This script creates a set of stream-segment adjustment factors
# for data collected during a Rapid Bio_Assessment (RBA) survey of the streams
# in a watershed.  Adjustment factors are based on X,Y coordinates added to
# the stream data and the reported cumulative distance at that point, compared
# to the streamline distance of the point.  Adjustment factors are written to
# a CSV file to facilitate review and direct update.
#
# INSTRUCTIONS:
#       Run the script at the command line. Use "-h" to view the input
#       arguments.
#
#       Input:
#          geodatabase: full path location of geodatabase containing a
#              "streams" feature class.  Each stream should be represented by a
#              single polyline feature, with linear distance originating at
#              the stream mouth.
#          survey_data_filepath: full path location of csv file containing
#              survey data, with x,y coordinates and optional notes added for
#              some of the pools.
#          sdi_filepath: name of csv file where stream distance information
#              will be written.
#          --sync_lat_long: indicates whether x,y coordinates in
#              survey_data_filepath are in Lat/Long (decimal degrees) or in
#              the coordinates of the stream layer (default)
#
#       Output:
#          Script returns 0 if it completes successfully, 1 if it does not.
#          It creates or overwrites the csv file at sdi_filepath.
#
#          Informational messages are logged to the console.  Debug-level
#          logging is available.
#
#       Exceptions:
#          Problem locating given files are handled and reported.
#          Other exceptions are not handled.
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
import csv
import argparse
import arcpy
from collections import namedtuple
import logging
import RBA_georef_util as rgutil


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
LOG_LEVEL = logging.INFO


# ********** CLASSES **********

# See RBA_georef_utl


# ********** FUNCTIONS **********

def parse_args(argv):
    """
    Defines and parses input arguments.
    :param argv: Input arguments, excluding the script name.
    :return: Argument values:
        geodatabase: path to geodatabase containing stream polylines
        csv_data_filepath: path to CSV file containing survey data with
            x,y coordinates for some pools
        sdi_filepath: path to possibly new CSV file where adjustment factors
            will be written
        sync_coords_in_lat_long: indicates whether x,y coordinates are in
            Lat/Long decimal degrees or in the coordinates of the stream
            layer (default)
    """
    parser = argparse.ArgumentParser\
        (description="Create a table of distance adjustment factors for survey data.")
    # positional arguments
    parser.add_argument("geodatabase", type=rgutil.valid_gdb,
                        help="full path location of geodatabase containing streams")
    parser.add_argument("survey_data_filepath", type=rgutil.valid_file,
                        help="full path location of csv file containing survey data " +
                              "with x,y coordinates for some of the pools")
    parser.add_argument("sdi_filepath", type=rgutil.valid_filedir,
                        help="stream distance information will be saved in this file")
    # optional arguments
    parser.add_argument("--sync_lat_long", dest="sync_coords_in_lat_long",
                        action='store_true')
    parser.set_defaults(sync_coords_in_lat_long=False)
    args = parser.parse_args(argv)
    return args.geodatabase, args.survey_data_filepath, args.sdi_filepath, \
           args.sync_coords_in_lat_long


def build_streamlength_adjustment_factor_dictionary(in_csv_filename,
                                                    streams_pathname,
                                                    sync_coords_in_lat_long):
    """
    Builds dictionary containing adjustment factors for stream segments,
    based on survey cumulative distance vs. stream polyline distance
    between synchronization points.
    :param in_csv_filename: CSV file containing RBA data plus XY sync point
        fields X, Y, and XY_Note.  File is assumed to be sorted by stream
        location ID (LLID) and cumulative distance.
    :param streams_pathname: stream polyline feature class, containing one
        polyline per stream, identified by location ID (LLID), with
        distance oriented from mouth to source.
    :param sync_coords_in_lat_long: True if XY data in in_csv_filename
        is in lat/long decimal degrees, False if XY data is in same reference
        system as streams_pathname
    :return: dictionary of stream distance adjustment information, keyed on
        stream LLID.  Each value contains a sequence of tuples:
        (begining_SycnPoint, ending_SyncPoint, adjustment_factor)
    """
    # We need to create a dictionary covering all reported cumulative distances.
    # However, the input data will contain x,y points for only some of these.
    # In order to ensure the dictionary begin/end points cover the full range,
    # we loop through the input csv file, one row at a time, keeping track of
    # which stream we're on and making sure this row is covered by a
    # (beginSyncPoint, endSyncPiont, adjustment_factor) tuple.
    stream_distance_info_dict = {}
    adj_factors = []
    prev_llid = ""
    need_adj_factor = False
    stream_geom = None
    end_sync_point = None
    with open(in_csv_filename, 'rb') as pts_file:
        # Read and process each row in csv file as namedtuple
        pts_file_reader = csv.reader(pts_file)
        headings = pts_file_reader.next()
        Row = namedtuple('Row',headings)
        for r in pts_file_reader:
            row = Row(*r)
            new_llid = str(row.LLID_num)
            if new_llid == "":  # if LLID is not given, log this and continue
                logging.warning(" No Location ID given for input data: " +
                                "stream {}, trib to {}, pool {}. Skipping entry".
                                format(row.STREAM, row.TRIB_TO, row.Pool_num))
            else:
                logging.debug(" read row = {}".format(row))
                pool_x_coord, pool_y_coord = get_pool_XY_coords(row)
                pool_cum_dist = int(row.CUM_DIST)
                xy_note = row.XY_Note
                survey_comment = row.COMMENT
                streamname = row.STREAM
                trib_to = row.TRIB_TO

                if new_llid != prev_llid:
                    # New Stream data
                    if prev_llid != "":
                        # Wrap up adjustment factor processing for previous stream
                        if need_adj_factor:
                            # if we're here, previous end_point was None, so we
                            # need an entry for the final adjustment factor
                            adj_factor = compute_adj_factor(begin_sync_point, None)
                            adj_factors.append((begin_sync_point, None, adj_factor))
                            need_adj_factor = False
                        # store all adj_factors for prev_llid
                        rgutil.add_adj_factors_to_sdi(stream_distance_info_dict,
                                                      prev_llid, adj_factors)
                    # Reset stream variables and initialize dictionary entry
                    # for new stream
                    adj_factors = []
                    stream_distance_info_dict[new_llid] = \
                        rgutil.new_sdi_object(new_llid, streamname, trib_to,
                                              adj_factors)

                    # Find geometry object for new stream
                    stream_geom = rgutil.get_stream_geom(streams_pathname,
                                                         new_llid)

                    # Find begin sync point for new stream
                    if (pool_x_coord is None) | (pool_y_coord is None):
                        begin_sync_point = new_syncpt_using_survey_dist\
                                (pool_cum_dist, xy_note, survey_comment)
                    else:
                        begin_sync_point = compute_xy_sync_point\
                            (stream_geom, pool_x_coord, pool_y_coord,
                             pool_cum_dist, xy_note, survey_comment,
                             sync_coords_in_lat_long)
                    prev_llid = new_llid
                    end_sync_point = begin_sync_point  # in case there is only
                                                       # one pool on this stream
                    need_adj_factor = True
                else:
                    # Data for another pool on the same stream
                    if (pool_x_coord is None) | (pool_y_coord is None):
                        # No xy coord given for this pool
                        end_sync_point = None
                        need_adj_factor = True
                    else:
                        end_sync_point = compute_xy_sync_point\
                            (stream_geom, pool_x_coord, pool_y_coord,
                             pool_cum_dist, xy_note, survey_comment,
                             sync_coords_in_lat_long)
                        # calculate adjustment factor
                        adj_factor = compute_adj_factor(begin_sync_point,
                                                        end_sync_point)
                        # add adjustment factor to list for this stream
                        adj_factors.append((begin_sync_point, end_sync_point,
                                            adj_factor))
                        # set begin sync point for next survey stream segment
                        begin_sync_point = end_sync_point
                        need_adj_factor = False  # not needed unless there is
                                                 # another row of data

        # End of file, no more rows
        # Tie up processing for last stream
        if need_adj_factor:
            adj_factor = compute_adj_factor(begin_sync_point, None)
            adj_factors.append((begin_sync_point, None, adj_factor))
        rgutil.add_adj_factors_to_sdi(stream_distance_info_dict, prev_llid,
                                      adj_factors)

    pts_file.close()
    return stream_distance_info_dict


def new_syncpt_using_survey_dist(pool_cum_dist, xy_note, survey_comment):
    """
    Create a new SyncPoint object with both survey and streamline
    cumulative distances set to input pool cumulative distance.
    :param pool_cum_dist: cumulative distance to use for both survey and
        streamline cum_dist in SyncPoint object.
    :param xy_note: text field with notes about xy sync (may exist
        regardless of values in other fields)
    :param survey_comment: text fields with notes about survey data/point
    :return: new SyncPoint object with x and y coordinates set to None and
        other fields populated based on input values
    """
    syncpt = rgutil.SyncPoint()
    syncpt.survey_cum_dist = pool_cum_dist
    syncpt.streamline_cum_dist = pool_cum_dist
    syncpt.x_coord = None
    syncpt.y_coord = None
    syncpt.xy_note = xy_note
    syncpt.survey_comment = survey_comment
    return syncpt


def compute_xy_sync_point(stream_geom, in_x_coord, in_y_coord, survey_cum_dist,
                          xy_note, survey_comment, sync_coords_in_lat_long):
    """
    Creates a new SyncPoint object with streamline distance based on
    input x and y coordinates.
    :param stream_geom: geometry object for stream
    :param in_x_coord: x coordinate for sync point
    :param in_y_coord: y coordinate for sync point
    :param survey_cum_dist: survey-reported cumulative distance
    :param xy_note: text field with notes about xy sync (may exist
        regardless of values in other fields)
    :param survey_comment: text fields with notes about survey data/point
    :param sync_coords_in_lat_long: True if X and Y coordinates are in lat/long
        (decimal degrees), False if XY data is in same reference system as
        stream_geom
    :return: new SyncPoint object with all fields populated, including
        streamline_cum_dist based on stream distance to x and y coordinates
    """
    logging.debug(" compute_xy_sync_point called " +
                  "with point ({}, {}) and survey_cum_dist {}".
                  format(in_x_coord, in_y_coord, survey_cum_dist))
    xy_point = arcpy.Point(in_x_coord, in_y_coord)
    if sync_coords_in_lat_long:
        xy_pt_geom = arcpy.PointGeometry(xy_point, LAT_LONG_CRS)
        xy_pt_geom = xy_pt_geom.projectAs(stream_geom.spatialReference)
    else:
        xy_pt_geom = arcpy.PointGeometry(xy_point)
        # assume same CRS as streams
    logging.debug(" xy_pt_geom = ({}, {})".format(xy_pt_geom.firstPoint.X,
                                                  xy_pt_geom.firstPoint.Y))

    # Make point for coordinates snapped to stream, and get
    # cumulative distance, etc.
    snapped_in_point, streamline_cum_dist, offset_dist, side = \
        stream_geom.queryPointAndDistance(xy_pt_geom)

    logging.debug(" calculated streamline_cum_dist = {}, offset_dist = {}".
                  format(streamline_cum_dist, offset_dist))
    syncpt = rgutil.SyncPoint()
    syncpt.survey_cum_dist = survey_cum_dist
    syncpt.streamline_cum_dist = streamline_cum_dist
    syncpt.x_coord = in_x_coord
    syncpt.y_coord = in_y_coord
    syncpt.xy_note = xy_note
    syncpt.survey_comment = survey_comment
    return syncpt


def compute_adj_factor(begin_sync_point, end_sync_point):
    """
    Computes an multiplicative adjustment factor to apply to reported survey
    cumulative distance measurements between two points, based on the
    relative difference between the reported survey and measured stream
    distances.
    :param begin_sync_point: SyncPoint object holding data for beginning of
        segment
    :param end_sync_point:  SyncPoint object holding data for end of segment
    :return: adjustment factor (double/float)
    """
    if end_sync_point:
        # compute adjustment factor from inputs:
        # linear geom diff betw start and end divided by
        # reported survey diff betw start and end
        adj_factor = (end_sync_point.streamline_cum_dist -
                      begin_sync_point.streamline_cum_dist) \
                    / (end_sync_point.survey_cum_dist -
                       begin_sync_point.survey_cum_dist)
    else:
        # estimate adjustment factor for end of stream
        adj_factor = DEFAULT_ADJ_FACTOR
    logging.debug(" adjustment factor {} calculated for {}, {}".
                  format(adj_factor, begin_sync_point, end_sync_point))
    return adj_factor

def get_pool_XY_coords(row):
    """
    Parses X and Y coordinates from input row.
    :param row: named tuple containing 1 row of data;
        empty fields contain empty strings.
    :return: tuple of (x_coord, y_coord), where coords may be None
    """
    if (row.X == ""):
        pool_x_coord = None
    else:
        pool_x_coord = float(row.X)
    if (row.Y == ""):
        pool_y_coord = None
    else:
        pool_y_coord = float(row.Y)

    return pool_x_coord, pool_y_coord


# ********** MAIN **********

def main(gdb_path, survey_data_filename, sdi_filepath,
         sync_coords_in_lat_long):

    # Initialize
    logging.basicConfig(level=LOG_LEVEL)
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = gdb_path

    # Get streams feature class path
    streams_pathname = rgutil.get_valid_polyline_pathname(gdb_path,
                                                          STREAMS_FC_NAME)

    # Build dictionary of stream distance information, including
    # adjustment factors for segments with x,y coordinates
    stream_distance_info = build_streamlength_adjustment_factor_dictionary\
        (survey_data_filename, streams_pathname, sync_coords_in_lat_long)

    # Write stream distance info to named csv file
    rgutil.write_sdi_to_csv_file(stream_distance_info,
                                sdi_filepath)
    logging.info(" developed adjustment factors for {} streams, saved to {}".
                 format(len(stream_distance_info.keys()),
                        sdi_filepath))
    return 0


# ********** MAIN CHECK **********

if __name__ == '__main__':
    sys.exit(main(*parse_args(sys.argv[1:])))
