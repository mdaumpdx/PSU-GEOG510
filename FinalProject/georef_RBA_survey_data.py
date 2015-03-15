# **********************************************************************
# 
# NAME: Marilyn Daum
# DATE: 15 March 2015
# CLASS: GEOG510
# ASSIGNMENT: Final Project
# 
# DESCRIPTION: This script creates a point feature class containing the
# data collected during a Rapid Bio_Assessment (RBA) survey of the streams
# in a watershed.  Point locations are determined as distance upstream,
# adjusted per a set of stream-segment adjustment factors.
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
#          sdi_filepath: name of csv file containing stream distance
#              information, including possibly modified adjustment factors.
#          survey_data_fc_name: name of feature class where survey data
#              will be written (within geodatabase)
#          survey_data_template: file with field definitions to use as
#              template for survey data feature class
#
#       Output:
#          Script returns 0 if it completes successfully, 1 if it does not.
#          It creates or overwrites the survey data feature class.
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
#            http://stackoverflow.com/questions/26502775/pycharm-is-telling-me-i-need-to-simplify-chained-comparison-how-can-this-comp
#
# **********************************************************************

# ********** IMPORT STATEMENTS **********
import sys
import os
import csv
import argparse
import arcpy
from arcpy import env
from collections import namedtuple
import logging
import RBA_georef_util as rgutil


# ********** GLOBAL CONSTANTS **********

STREAM_FC_LLID = "LocationID"

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
EXCLUDED_NEW_FIELD_NAMES = [u'FID', u'OBJECTID', u'Shape']

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
        csv_data_filepath: path to CSV file containing survey data with x,y coordinates for some pools
        sdi_filepath: path to possibly new CSV file where adjustment factors will be written
        survey_data_fc_name: name of feature class where survey data will be stored (in gdb)
    """
    parser = argparse.ArgumentParser\
        (description="Create a table of distance adjustment factors for survey data.")
    # positional arguments
    parser.add_argument("geodatabase", type=rgutil.valid_gdb,
                        help="full path location of geodatabase containing streams")
    parser.add_argument("survey_data_filepath", type=rgutil.valid_file,
                        help="full path location of csv file containing survey data")
    parser.add_argument("sdi_filepath", type=rgutil.valid_file,
                        help="full path location of csv file containing " +
                             "stream distance information")
    parser.add_argument("survey_data_fc_name", type=str,
                        help="name of feature class where survey data will " +
                             "be stored (within geodatabase)")
    parser.add_argument("survey_data_template", type=rgutil.valid_file,
                        help="file with field definitions to use as template for survey data feature class")
    args = parser.parse_args(argv)
    return args.geodatabase, args.survey_data_filepath, args.sdi_filepath, \
           args.survey_data_fc_name, args.survey_data_template


def georeference_survey_data(survey_data_filename, stream_dist_info_dict,
                             streams_pathname, survey_data_fc,
                             survey_data_template):
    """
    Creates points in survey_data_fc for rows in survey_data_filename,
    with points located at calculated distances on streams in streams_pathname.
    Input stream_dist_info_dict is used to adjust reported cumulative distance.
    :param survey_data_filename: CSV file containing RBA data plus XY sync point
        fields X, Y, and XY_Note.
    :param stream_dist_info_dict: Dictionary of stream distance information
        keyed on stream LLID.  Each value contains a list of tuples:
        (begining_SycnPoint, ending_SyncPoint, adjustment_factor)
    :param streams_pathname: stream polyline feature class, containing one
        polyline per stream, identified by location ID (LLID), with
        distance oriented from mouth to source.
    :param survey_data_fc: Feature class to which new survey data points are
        added.
    :param survey_data_template: Feature class containing schema for fields
        in survey_data_filename
    :return: N/A; survey_data_fc is update by this function.
    """
    stream_geom = None
    prev_llid = ""
    # List of fields added to survey_data_fc, based on survey_data_template
    insert_fields = [desc_field.name for
                     desc_field in arcpy.Describe(survey_data_template).fields
                     if desc_field.name not in EXCLUDED_NEW_FIELD_NAMES]
    # Create InsertCursor for adding new survey data points
    with arcpy.da.InsertCursor (survey_data_fc, ["SHAPE@"] + insert_fields) \
            as insertCursor:
        with open(survey_data_filename, 'rb') as pts_file:
            # Read and process each row in csv file as namedtuple
            pts_file_reader = csv.reader(pts_file)
            headings = pts_file_reader.next()
            Row = namedtuple('Row',headings)
            for r in pts_file_reader:
                row = Row(*r)
                logging.debug(" read row = {}".format(row))
                new_llid = str(row.LLID_num)
                streamname = str(row.STREAM)
                trib_to = str(row.TRIB_TO)
                pool_cum_dist = int(row.CUM_DIST)

                if new_llid == "":  # if LLID is not given, log this and continue
                    logging.warning(" No Location ID given for input data: " +
                                    "stream {}, trib to {}, pool {}. Skipping entry".
                                    format(streamname, trib_to, row.Pool_num))
                else:
                    if new_llid != prev_llid:
                        # New Stream
                        logging.info(" Georeferencing data for {} trib to {}".
                                     format(streamname, trib_to))
                        stream_adj_factors = \
                            stream_dist_info_dict[new_llid].adj_factors
                        # get stream geometry object
                        stream_geom = rgutil.get_stream_geom(streams_pathname,
                                                             new_llid)
                        prev_llid = new_llid

                    # Compute adjusted distance for this row
                    adjusted_distance = \
                        adjust_stream_distance(pool_cum_dist, stream_adj_factors)
                    # Georeference the survey data for this row
                    create_point_upstream(stream_geom, adjusted_distance,
                                          row, insertCursor)

def adjust_stream_distance(survey_dist, stream_adj_factors):
    """
    Computes an adjusted cumulative distance, based on survey_dist
    and stream adjustment factors.
    :param survey_dist: input survey distance, in feet
    :param stream_adj_factors: sequence of tuples with adjustment factors
        for this stream; each tuple contains
        (begin_sync_point, end_sync_point, adj_factor)
    :return: adjusted survey distance = input value adjusted per
        the adj_factor for the stream segment applicable to
        the input survey_dist
    """
    # return survey_dist # uncomment to apply 1 to all distances
    adj_factor = rgutil.DEFAULT_ADJ_FACTOR
    for factor in stream_adj_factors:
        begin_sync_pt_dist = factor[0].survey_cum_dist
        begin_streamline_pt_dist = factor[0].streamline_cum_dist
        end_sync_pt_dist = factor[1].survey_cum_dist
        # Find begin and end sync points that include survey_dist;
        # end sync point is inclusive to ensure final end points
        # get picked up (intermediate sync points can be applied
        # at either beginning or end of segment).
        if begin_sync_pt_dist <= survey_dist <= end_sync_pt_dist:
            adj_factor = factor[2]
            break
    new_dist = begin_streamline_pt_dist + \
               ((survey_dist - begin_sync_pt_dist) * adj_factor)
    logging.debug(" for survey_dist {}, ".format(survey_dist) +
                  "computed adjusted distance as " +
                  " {} + ({} - {}) * {}, yielding {}".
                  format(begin_streamline_pt_dist, survey_dist,
                         begin_sync_pt_dist, adj_factor, new_dist))
    return new_dist


def create_point_upstream(line_geom, distance, data_row, insertCursor):
    """
    Creates a new Point geometry object, located at the given distance
    upstream along line_geom. Inserts a row for this point, with fields
    from data_row.  Fields in data_row are assumed to contain RBA data
    plus XY sync point fields in the following order:
    ENTRY, YEAR, DATE, BASIN, TRIB_TO, STREAM, LLID_num, s_GUID,
    Channel_Type, Pool_num, s_Lineage, RAND, TYPE, LENGTH, WIDTH, VIS, COMP,
    DIST, CUM_DIST, COHO, Zero_plus, STHD, CUT, CHIN, RES_RB, CUL, KNOT, LONG,
    BEAVER_DAMS, Num_BEAVER_DAMS, GRAVEL_COUNT, COMMENT, XY_Text, X, Y, XY_Note
    Column numbers corresponding to these are hardcoded in this function.
    :param line_geom: Polyline Geometry object for stream
    :param distance: Distance from mouth of stream to locate new point
    :param data_row: namedtuple containing fields listed above
    :param insertCursor: cursor for inserting new Point Geometry
    :return: N/A, insertCursor is updated as a result of this function.
    """

    # Create new point at given distance along line_geom
    pt_geom = line_geom.positionAlongLine(distance)

    # Prepare data fields for use with insertCursor.

    # Remove sync_point xy fields, retaining original survey data
    data_row_survey_fields = list(data_row)[:-5]
    data_row_survey_fields.append(data_row.COMMENT)

    # Change blank fish counts to 0 (columms 19:24)
    for i in range(19, 25):
        if data_row_survey_fields[i] == '':
            data_row_survey_fields[i] = '0'

    # Limit comments (last column) to 254 chars.  Note this
    # limitation can be removed by using a fileGDB template
    # for the fields, instead of a shapefile.
    comment = data_row_survey_fields[-1]
    if len(comment) > 254:
        data_row_survey_fields[-1] = comment[:252] + '..'
    logging.debug(" data_row_survey_fields= {}".format(data_row_survey_fields))

    # Add row for new point geometry and fields to insertCursor
    insertCursor.insertRow([pt_geom] + list(data_row_survey_fields))


# ********** MAIN **********

def main(gdb_path, survey_data_filename, sdi_filepath,
         survey_data_fc_name, survey_data_template):

    # Initialize
    logging.basicConfig(level=LOG_LEVEL)
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = gdb_path

    # Get streams feature class
    streams_pathname = rgutil.get_valid_polyline_pathname\
        (gdb_path, STREAMS_FC_NAME)
    streams_spat_ref = arcpy.Describe(streams_pathname).SpatialReference

    # Create new feature class for survey data
    if arcpy.Exists(survey_data_fc_name):
        arcpy.Delete_management(survey_data_fc_name)
    survey_data_fc = arcpy.CreateFeatureclass_management\
        (gdb_path, survey_data_fc_name, "POINT",
         survey_data_template,
         spatial_reference=streams_spat_ref)

    # Populate dictionary of stream distance adjustment factors
    stream_dist_info_dict = rgutil.read_sdi_from_csvfile(sdi_filepath)
    logging.debug(" stream_dist_info_dict = {}".format(stream_dist_info_dict))

    # Create points for survey data
    georeference_survey_data(survey_data_filename, stream_dist_info_dict,
                             streams_pathname, survey_data_fc,
                             survey_data_template)

    return 0


# ********** MAIN CHECK **********

if __name__ == '__main__':
   sys.exit(main(*parse_args(sys.argv[1:])))