# -*- coding: utf-8 -*-

"""
RADIS includes [HAPI]_ for calculation of partition functions
at equilibrium, and comparaison with the RADIS LBL code in validation 
cases.  

References
----------

[HAPI]_


--------------------------------------------------------------------------------

This module provides an access to the HITRAN data.
Data is downloaded and cached.

This module serves as a simple database manager frontend.

API is aimed to be RESTful, which means that interaction
between local API and remote data-server will be held 
via sending RESTful queries (API->remote) and
receiving data preferrably in text format (remote->API).

Object are supposed to be implemented by structures/dicts
as they present in almost any programming language.

Trying to retain functional style for this API. 
"""

from __future__ import print_function, absolute_import, division, unicode_literals

# import httplib
# import urllib2
import json
import os
import os.path
import re
from os import listdir
from numpy import zeros, array, zeros, setdiff1d, ndarray, arange
from numpy import place, where, insert, real, polyval
from numpy import complex128, complex64, int64, int32, float64, float32
from numpy import sqrt, abs, exp, pi, log, sin, cos, tan
from numpy import convolve
from numpy import flipud
from numpy.fft import fft, fftshift
from numpy import linspace, floor
from numpy import any, minimum, maximum
from numpy import modf
from numpy import sort as npsort
from bisect import bisect

# from collections import OrderedDict
from warnings import warn, simplefilter
import pydoc
import six
from six.moves import range
from six.moves import zip

# Enable warning repetitions
simplefilter("always", UserWarning)

# Python 3 compatibility
try:
    import urllib.request as urllib2
except ImportError:
    import six.moves.urllib.request, six.moves.urllib.error, six.moves.urllib.parse

HAPI_VERSION = "1.1.0.6"
# CHANGES:
# FIXED GRID BUG (ver. 1.1.0.1)
# FIXED OUTPUT FORMAT FOR CROSS-SECTIONS (ver. 1.1.0.1)
# ADDED CPF BY SCHREIER (JQSRT_112_2011) (ver. 1.1.0.2)
# OPTIMIZED EXPRESSION EVALUATIONS FOR SELECT (ver. 1.1.0.3)
# ADDED SUPPORT FOR MIXTURES (ver. 1.1.0.4)
# ADDED SUPPORT FOR USER-DEFINED ENV DEPENDENCES (ver. 1.1.0.5)
# ADDED PROFILE SELECTION (ALPHA) (ver. 1.1.0.6)

# version header
# print('HAPI version: %s' % HAPI_VERSION)
# print('To get the most up-to-date version please check http://hitran.org/hapi')

# define precision
__ComplexType__ = complex128
__IntegerType__ = int64
__FloatType__ = float64

# define zero
cZero = __FloatType__(0.0)

# physical constants
cBolts = 1.380648813e-16  # erg/K, CGS
cc = 2.99792458e10  # cm/s, CGS
hh = 6.626196e-27  # erg*s, CGS

# computational constants
cSqrtLn2divSqrtPi = 0.469718639319144059835
cLn2 = 0.6931471805599
cSqrtLn2 = 0.8325546111577
cSqrt2Ln2 = 1.1774100225

# declare global variables

GLOBAL_DEBUG = False
if GLOBAL_DEBUG:
    warn("GLOBAL_DEBUG is set to True!")

GLOBAL_CURRENT_DIR = "."

GLOBAL_HITRAN_APIKEY = "e20e4bd3-e12c-4931-99e0-4c06e88536bd"

GLOBAL_USER = "user"
GLOBAL_REQUISITES = []

GLOBAL_CONNECTION = []
GLOBAL_DATABASE = "hitran"

LOCAL_HOST = "http://localhost"

# DEBUG switch
if GLOBAL_DEBUG:
    GLOBAL_HOST = LOCAL_HOST + ":8000"  # localhost
else:
    GLOBAL_HOST = "http://hitran.org"

# this is a backup url in the case GLOBAL_HOST does not work
GLOBAL_HOST_BACKUP = "http://hitranazure.cloudapp.net/"

# In this "robust" version of arange the grid doesn't suffer
# from the shift of the nodes due to error accumulation.
# This effect is pronounced only if the step is sufficiently small.


def arange_(lower, upper, step):
    npnt = int(floor((upper - lower) / step) + 1)
    upper_new = lower + step * (npnt - 1)
    if abs((upper - upper_new) - step) < 1e-10:
        upper_new += step
        npnt += 1
    return linspace(lower, upper_new, npnt)


# interface for checking of variable's existance
def empty(Instance):
    return True if Instance else False


# general interface for getattr


def getAttribute(Object, Attribute):
    return getattr(Object, Attribute)


# general interface for setattr


def setAttribute(Object, Attribute, Value):
    setattr(Object, Attribute, Value)
    return


# UNPARSED QUERY OBJECT
# uses formal language (SQL, noSQL, custom...)
GlobalQueryString = ""

# PARSED QUERY OBJECT
# = prototype for a Query instance
# there should be a getAttrbute/setSettribute functions defined
# For Django: Query=QuerySet (as  an example)
Query = {}

# prototype for cache storage
# there must be function for record/retrieve
# caching is performed by the value of Query
# cache parameters: (query+table_name)
# if there is already table with such query, copy it
# if there is already tble with such query AND table_name,
# return it as is => IT MAY DEPEND ON CERTAIN QUERY TYPE!!
TABLES = {}  # hash/dictionary

# ---------- CONNECTION MANAGEMENT-------------

# interface for establishing HTTP connection
# can return object/structure/handle


def setupConnection(Host=GLOBAL_HOST):
    Connection = httplib.HTTPConnection(Host)
    if not empty(Connection):
        return Connection
    else:
        raise Exception("can" "t setup connection")


# interface for HTTP-get method
# Connection must be established before use


def httpGet(URL, Connection=GLOBAL_CONNECTION):
    Method = "get"
    ServerResponse = Connection.request(Method, URL)
    return ServerResponse


# parse local data language to remote frontend


def parseToFrontend(Query, Host=GLOBAL_HOST):
    # convert Query object to server frontend's
    # query language
    pass


def prepareURL(Query, Connection=GLOBAL_CONNECTION):
    # make full URL from server name and it's parameters
    # considering server's frontend query language
    Host = getAttribute(Connection, "host")
    HostQuery = parseToFrontend(Query)
    URL = Host + HostQuery
    return URL


# stream raw data from the server
# the data is assumed to be very large that
# ordinary get is unefficient


def streamRawDataRemote(Query, Connection=GLOBAL_CONNECTION):
    pass


# collect raw data in whatever format server gives it


def getRawDataRemote(Query, Connection=GLOBAL_CONNECTION):
    URL = prepareURL(Query, Connection)
    ServerResponse = httpGet(URL, Connection)
    return ServerResponse


# parse raw data
# def parseRawData(RawData)
#    pass

# ---------- CONNECTION MANAGEMEND END --------

# Two types of interaction between API and DB:
# 1) via API library
# 2) via REST http protocol (torrent-like)

# ---------- NODE MANAGEMENT ------------------

# An interface for a node manager will follow soon.
# This is an implementation in Python
# Different implementations are language-specific.

# Default node with simple DB engine
# Prototype for a global nodelist for a given host

# Each node has it's unique ID, host name and
#   node name within it's host


NODE_NAME = "local"

GLOBAL_NODENAMES = {0: "hitran-main", 1: "local"}

GLOBAL_NODELIST = {
    0: {  # main HITRAN node
        "host": GLOBAL_HOST,
        "ACCESS_KEY": "9b6a7975-2a84-43d8-920e-f4dea9db6805",  # guest
    },
    1: {  # local node prototype
        "host": LOCAL_HOST,
        "ACCESS_KEY": "6cfd7040-24a6-4197-81f9-6e25e50005b2",  # admin
    },
}


def createNode(NodeID, NodeList=GLOBAL_NODELIST):
    # create a node, throw if exists
    node = NodeList.get(NodeID)
    if node:
        raise Exception("node %s already exists" % NodeName)
    NodeList[NodeID] = {}
    pass


def getNodeIDs(NodeList=GLOBAL_NODELIST):
    # return list of all available nodes
    return list(NodeList.keys())


def getNodeProperty(NodeID, PropName, NodeList=GLOBAL_NODELIST):
    # get a property for certain node
    # if not found throw exception
    node = NodeList.get(NodeName)
    if node:
        prop = node.get(PropName)
        if prop:
            return prop
        else:
            raise Exception("node %s doesn" "t have property %s" % (ModeName, Propname))
    else:
        raise Exception("no such node %s" % Nodename)


def setNodeProperty(NodeID, PropName, PropValue, NodeList=GLOBAL_NODELIST):
    # set a property for certain node
    # throw exception if node not found
    # if the property doesn't exist it will appear
    node = NodeList.get(NodeID)
    if not node:
        raise Exception("no such node %s " % NodeName)
    NodeList[PropName] = PropValue
    return


def resolveNodeID(NodeName, NodeNames=GLOBAL_NODENAMES):
    for NodeID in NodeNames.keys():
        if NodeNames[NodeID] == NodeName:
            return NodeID


def checkAccess(
    DBName,
    TableName,
    NodeName,
    UserName,
    Requisites,
    NodeList=GLOBAL_NODELIST,
    NodeNames=GLOBAL_NODENAMES,
):
    # simple node-level authentication (bridge to AUTH system)
    NodeID = resolveNodeID(NodeName, NodeNames)
    Node = NodeList[NodeID]
    if Requisites.key in Node["keys_allowed"]:
        return True
    else:
        return False


# ---------- NODE MANAGEMENT END --------------

# ---------- NODE AUTH SYSTEM -----------------

# AUTH SYSTEM is tightly connected to Node manager.

# Prototype for authentication  system.
# AUTH is responsible for giving an access privileges to all users.
# Each users has a key ACCESS_KEY which is stored in
#  a special database HOST:ACCESS_KEYS on a host.
# Every node has a separate privileges list connected with
#  each key.

# The current auth system is based on secret keys of access
# Default key is 'admin', it's created seamlessly for a local admin.


GLOBAL_PRIVILEGES = {
    "admin": {"ACCESS_KEY": "6cfd7040-24a6-4197-81f9-6e25e50005b2", "LEVEL": "ADMIN"},
    "guest": {"ACCESS_KEY": "9b6a7975-2a84-43d8-920e-f4dea9db6805", "LEVEL": "USER"},
}


def addUser():
    pass


def deleteUser():
    pass


def authenticate(UserName, Requisites, Privileges=GLOBAL_PRIVILEGES):
    # Authentication
    key_list = [Privileges[User]["ACCESS_KEY"] for User in Privileges.keys]
    return True if Requisites.AccessKey in key_list else False


def checkPrivileges(
    Path,
    UserName=GLOBAL_USER,
    Requisites=GLOBAL_REQUISITES,
    Privileges=GLOBAL_PRIVILEGES,
    NodeList=GLOBAL_NODELIST,
    Nodenames=GLOBAL_NODENAMES,
):
    # Privileges are checked before executing every query (needs optimization)
    # Path example: SOME_DB::SOME_TABLE::SOME_NODE
    if not authenticate(UserName, Requisites, Privileges):
        return False
    (DBName, TableName, NodeName) = Path.split("::")
    # loop on all nodes , use NODE_MANAGER's functions instead of
    #   working with GLOBAL_NODELIST directly
    if not checkAccess(
        DBName, TableName, NodeName, UserName, Requisites, NodeList, NodeNames
    ):
        return False
    return True


# ---------- NODE AUTH SYSTEM END -------------

# ---------- DATABASE FRONTEND ----------------

# Structure:
#   DB::TABLE::NODE
#     DB - distributed database
#     TABLE - table within the current database
#     NODE - instance of this API with fixed DB backend
#      !! parameter HOST is deprecated

#     HOST - computer at which the NODE/ENGINE is deployed

# TABLE should be considered as schema-free collection
#  (e.g. MongoDB-type)

# Two databases (DB) - GLOBAL (one) and LOCAL (many)

# Every DB has an ACCESS_KEY providing an access to it
# User can create a database and it will contain
#  a list of ACCESS_KEY's for authentication.

# GLOBAL AND LOCAL are distributed databases.
# A user can create his GLOBAL database and open an access to it.
# GLOBAL access implementation:

# GLOBAL is a distributed database

# The DB frontend contains interfaces to
#  the standard procedures of data creation and
#  retrieval of an "average" DBMS.
#   ("collection" = table)
#
#   Levels of access: (DB permissions implementation)
#   0:USER     read-only operations ("select")
#   1:MANAGER  manage single DB (create/delete docs)
#   2:ADMIN    manage multiple DB's (create/delete DB)
#
#   Every ACCESS_KEY has it's own access level.
#
#   Commands to implement:
#
#   ) create DATABASE
#   ) create ACCESS_KEY
#      (seamlessly for the local user)
#   ) select from LOCAL/GLOBAL doc (cached!)
#   ) access database
#      (seamlessly for the local user)
#   ) create/delete doc
#   ) copy/clone LOCAL doc
#   ) "create collection as select * from HOST:ENGINE:DB:COLLECTION"
#      (other types of table creations are forbidden)

# DB frontend is adapted to denormalized
#  schema-fixed tables or schema-independent documents.

# DB frontend is connected to multiple backends
#  which are largely language-specific.

# ATTENTION: since the system is distributed,
# the table/document caching is supposed to
# be in the frontend.
# Current higher-level implementation
# implies the query-based caching, i.e.
# cache lookup is performed by the value
# of Query structure/object.


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# LOCAL DATABASE MANAGEMENT SYSTEM
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# DATABASE BACKEND: simple text files, parsed into a python lists
# Use a directory as a database. Each table is stored in a
# separate text file. Parameters in text are position-fixed.

BACKEND_DATABASE_NAME_DEFAULT = "."

VARIABLES = {}
VARIABLES["BACKEND_DATABASE_NAME"] = BACKEND_DATABASE_NAME_DEFAULT

# For this node local DB is schema-dependent!
LOCAL_TABLE_CACHE = {
    "sampletab": {  # table
        "header": {  # header
            "order": ("column1", "column2", "column3"),
            "format": {"column1": "%10d", "column2": "%20f", "column3": "%30s"},
            "default": {"column1": 0, "column2": 0.0, "column3": ""},
            "number_of_rows": 3,
            "size_in_bytes": None,
            "table_name": "sampletab",
            "table_type": "strict",
        },  # /header
        "data": {
            "column1": [1, 2, 3],
            "column2": [10.5, 11.5, 12.5],
            "column3": ["one", "two", "three"],
        },  # /data
    }  # /table
}  # hash-map of tables

# FORMAT CONVERSION LAYER

# converts between TRANSPORT_FORMAT and OBJECT_FORMAT
HITRAN_FORMAT_160 = {
    "M": {"pos": 1, "len": 2, "format": "%2d"},
    "I": {"pos": 3, "len": 1, "format": "%1d"},
    "nu": {"pos": 4, "len": 12, "format": "%12f"},
    "S": {"pos": 16, "len": 10, "format": "%10f"},
    "R": {"pos": 26, "len": 0, "format": "%0f"},
    "A": {"pos": 26, "len": 10, "format": "%10f"},
    "gamma_air": {"pos": 36, "len": 5, "format": "%5f"},
    "gamma_self": {"pos": 41, "len": 5, "format": "%5f"},
    "E_": {"pos": 46, "len": 10, "format": "%10f"},
    "n_air": {"pos": 56, "len": 4, "format": "%4f"},
    "delta_air": {"pos": 60, "len": 8, "format": "%8f"},
    "V": {"pos": 68, "len": 15, "format": "%15s"},
    "V_": {"pos": 83, "len": 15, "format": "%15s"},
    "Q": {"pos": 98, "len": 15, "format": "%15s"},
    "Q_": {"pos": 113, "len": 15, "format": "%15s"},
    "Ierr": {"pos": 128, "len": 6, "format": "%6s"},
    "Iref": {"pos": 134, "len": 12, "format": "%12s"},
    "flag": {"pos": 146, "len": 1, "format": "%1s"},
    "g": {"pos": 147, "len": 7, "format": "%7f"},
    "g_": {"pos": 154, "len": 7, "format": "%7f"},
}

# This should be generating from the server's response
HITRAN_DEFAULT_HEADER = {
    "table_type": "column-fixed",
    "size_in_bytes": -1,
    "table_name": "###",
    "number_of_rows": -1,
    "order": [
        "molec_id",
        "local_iso_id",
        "nu",
        "sw",
        "a",
        "gamma_air",
        "gamma_self",
        "elower",
        "n_air",
        "delta_air",
        "global_upper_quanta",
        "global_lower_quanta",
        "local_upper_quanta",
        "local_lower_quanta",
        "ierr",
        "iref",
        "line_mixing_flag",
        "gp",
        "gpp",
    ],
    "format": {
        "a": "%10.3E",
        "gamma_air": "%5.4f",
        "gp": "%7.1f",
        "local_iso_id": "%1d",
        "molec_id": "%2d",
        "sw": "%10.3E",
        "local_lower_quanta": "%15s",
        "local_upper_quanta": "%15s",
        "gpp": "%7.1f",
        "elower": "%10.4f",
        "n_air": "%4.2f",
        "delta_air": "%8.6f",
        "global_upper_quanta": "%15s",
        "iref": "%12s",
        "line_mixing_flag": "%1s",
        "ierr": "%6s",
        "nu": "%12.6f",
        "gamma_self": "%5.3f",
        "global_lower_quanta": "%15s",
    },
    "default": {
        "a": 0.0,
        "gamma_air": 0.0,
        "gp": "FFF",
        "local_iso_id": 0,
        "molec_id": 0,
        "sw": 0.0,
        "local_lower_quanta": "000",
        "local_upper_quanta": "000",
        "gpp": "FFF",
        "elower": 0.0,
        "n_air": 0.0,
        "delta_air": 0.0,
        "global_upper_quanta": "000",
        "iref": "EEE",
        "line_mixing_flag": "EEE",
        "ierr": "EEE",
        "nu": 0.0,
        "gamma_self": 0.0,
        "global_lower_quanta": "000",
    },
    "description": {
        "a": "Einstein A-coefficient in s-1",
        "gamma_air": "Air-broadened Lorentzian half-width at half-maximum at p = 1 atm and T = 296 K",
        "gp": "Upper state degeneracy",
        "local_iso_id": "Integer ID of a particular Isotopologue, unique only to a given molecule, in order or abundance (1 = most abundant)",
        "molec_id": "The HITRAN integer ID for this molecule in all its isotopologue forms",
        "sw": "Line intensity, multiplied by isotopologue abundance, at T = 296 K",
        "local_lower_quanta": "Rotational, hyperfine and other quantum numbers and labels for the lower state of a transition",
        "local_upper_quanta": "Rotational, hyperfine and other quantum numbers and labels for the upper state of a transition",
        "gpp": "Lower state degeneracy",
        "elower": "Lower-state energy",
        "n_air": "Temperature exponent for the air-broadened HWHM",
        "delta_air": "Pressure shift induced by air, referred to p=1 atm",
        "global_upper_quanta": "Electronic and vibrational quantum numbers and labels for the upper state of a transition",
        "iref": "Ordered list of reference identifiers for transition parameters",
        "line_mixing_flag": "A flag indicating the presence of additional data and code relating to line-mixing",
        "ierr": "Ordered list of indices corresponding to uncertainty estimates of transition parameters",
        "nu": "Transition wavenumber",
        "gamma_self": "Self-broadened HWHM at 1 atm pressure and 296 K",
        "global_lower_quanta": "Electronic and vibrational quantum numbers and labels for the lower state of a transition",
    },
}

PARAMETER_META = {
    "global_iso_id": {
        "id": 1,
        "name": "global_iso_id",
        "name_html": "Global isotopologue ID",
        "table_name": "",
        "description": "Unique integer ID of a particular isotopologue: every global isotopologue ID is unique to a particular species, even between different molecules. The number itself is, however arbitrary.",
        "description_html": "Unique integer ID of a particular isotopologue: every global isotopologue ID is unique to a particular species, even between different molecules. The number itself is, however arbitrary.",
        "default_fmt": "%5d",
        "default_units": "",
        "data_type": "int",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "molec_id": {
        "id": 2,
        "name": "molec_id",
        "name_html": "Molecule ID",
        "table_name": "",
        "description": "The HITRAN integer ID for this molecule in all its isotopologue forms",
        "description_html": "The HITRAN integer ID for this molecule in all its isotopologue forms",
        "default_fmt": "%2d",
        "default_units": None,
        "data_type": "int",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "local_iso_id": {
        "id": 3,
        "name": "local_iso_id",
        "name_html": "Isotopologue ID",
        "table_name": "",
        "description": "Integer ID of a particular Isotopologue, unique only to a given molecule, in order or abundance (1 = most abundant)",
        "description_html": "Integer ID of a particular Isotopologue, unique only to a given molecule, in order or abundance (1 = most abundant)",
        "default_fmt": "%1d",
        "default_units": "",
        "data_type": "int",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "nu": {
        "id": 4,
        "name": "nu",
        "name_html": "<em>&nu;</em>",
        "table_name": "prm_nu",
        "description": "Transition wavenumber",
        "description_html": "Transition wavenumber",
        "default_fmt": "%12.6f",
        "default_units": "cm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "sw": {
        "id": 5,
        "name": "sw",
        "name_html": "<em>S</em>",
        "table_name": "prm_sw",
        "description": "Line intensity, multiplied by isotopologue abundance, at T = 296 K",
        "description_html": "Line intensity, multiplied by isotopologue abundance, at T&nbsp;=&nbsp;296&nbsp;K",
        "default_fmt": "%10.3e",
        "default_units": "cm-1/(molec.cm-2)",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "a": {
        "id": 6,
        "name": "a",
        "name_html": "<em>A</em>",
        "table_name": "prm_a",
        "description": "Einstein A-coefficient in s-1",
        "description_html": "Einstein <em>A</em>-coefficient",
        "default_fmt": "%10.3e",
        "default_units": "s-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "gamma_air": {
        "id": 7,
        "name": "gamma_air",
        "name_html": "<em>&gamma;</em><sub>air</sub>",
        "table_name": "prm_gamma_air",
        "description": "Air-broadened Lorentzian half-width at half-maximum at p = 1 atm and T = 296 K",
        "description_html": "Air-broadened Lorentzian half-width at half-maximum at p&nbsp;=&nbsp;1&nbsp;atm and T&nbsp;=&nbsp;296&nbsp;K",
        "default_fmt": "%6.4f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "gamma_self": {
        "id": 8,
        "name": "gamma_self",
        "name_html": "<em>&gamma;</em><sub>self</sub>",
        "table_name": "prm_gamma_self",
        "description": "Self-broadened HWHM at 1 atm pressure and 296 K",
        "description_html": "Self-broadened HWHM at 1&nbsp;atm pressure and 296&nbsp;K",
        "default_fmt": "%5.3f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "n_air": {
        "id": 9,
        "name": "n_air",
        "name_html": "<em>n</em><sub>air</sub>",
        "table_name": "prm_n_air",
        "description": "Temperature exponent for the air-broadened HWHM",
        "description_html": "Temperature exponent for the air-broadened HWHM",
        "default_fmt": "%7.4f",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "delta_air": {
        "id": 10,
        "name": "delta_air",
        "name_html": "<em>&delta;</em><sub>air</sub>",
        "table_name": "prm_delta_air",
        "description": "Pressure shift induced by air, referred to p=1 atm",
        "description_html": "Pressure shift induced by air, referred to <em>p</em>=1&nbsp;atm",
        "default_fmt": "%9.6f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "elower": {
        "id": 11,
        "name": "elower",
        "name_html": '<em>E"</em>',
        "table_name": "",
        "description": "Lower-state energy",
        "description_html": "Lower-state energy",
        "default_fmt": "%10.4f",
        "default_units": "cm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "gp": {
        "id": 12,
        "name": "gp",
        "name_html": "<em>g</em>'",
        "table_name": "",
        "description": "Upper state degeneracy",
        "description_html": "Upper state degeneracy",
        "default_fmt": "%5d",
        "default_units": "",
        "data_type": "int",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "gpp": {
        "id": 13,
        "name": "gpp",
        "name_html": '<em>g</em>"',
        "table_name": "",
        "description": "Lower state degeneracy",
        "description_html": "Lower state degeneracy",
        "default_fmt": "%5d",
        "default_units": "",
        "data_type": "int",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "global_upper_quanta": {
        "id": 14,
        "name": "global_upper_quanta",
        "name_html": "Global upper quanta",
        "table_name": "",
        "description": "Electronic and vibrational quantum numbers and labels for the upper state of a transition",
        "description_html": "Electronic and vibrational quantum numbers and labels for the upper state of a transition",
        "default_fmt": "%15s",
        "default_units": None,
        "data_type": "str",
        "selectable": 0,
        "has_reference": 0,
        "has_error": 0,
    },
    "global_lower_quanta": {
        "id": 15,
        "name": "global_lower_quanta",
        "name_html": "Global lower quanta",
        "table_name": "",
        "description": "Electronic and vibrational quantum numbers and labels for the lower state of a transition",
        "description_html": "Electronic and vibrational quantum numbers and labels for the lower state of a transition",
        "default_fmt": "%15s",
        "default_units": None,
        "data_type": "str",
        "selectable": 0,
        "has_reference": 0,
        "has_error": 0,
    },
    "local_upper_quanta": {
        "id": 16,
        "name": "local_upper_quanta",
        "name_html": "Local upper quanta",
        "table_name": "",
        "description": "Rotational, hyperfine and other quantum numbers and labels for the upper state of a transition",
        "description_html": "Rotational, hyperfine and other quantum numbers and labels for the upper state of a transition",
        "default_fmt": "%15s",
        "default_units": None,
        "data_type": "str",
        "selectable": 0,
        "has_reference": 0,
        "has_error": 0,
    },
    "local_lower_quanta": {
        "id": 17,
        "name": "local_lower_quanta",
        "name_html": "Local lower quanta",
        "table_name": "",
        "description": "Rotational, hyperfine and other quantum numbers and labels for the lower state of a transition",
        "description_html": "Rotational, hyperfine and other quantum numbers and labels for the lower state of a transition",
        "default_fmt": "%15s",
        "default_units": None,
        "data_type": "str",
        "selectable": 0,
        "has_reference": 0,
        "has_error": 0,
    },
    "line_mixing_flag": {
        "id": 18,
        "name": "line_mixing_flag",
        "name_html": "Line mixing flag",
        "table_name": "",
        "description": "A flag indicating the presence of additional data and code relating to line-mixing",
        "description_html": "A flag indicating the presence of additional data and code relating to line-mixing",
        "default_fmt": "%1s",
        "default_units": "",
        "data_type": "str",
        "selectable": 0,
        "has_reference": 0,
        "has_error": 0,
    },
    "ierr": {
        "id": 19,
        "name": "ierr",
        "name_html": "Error indices",
        "table_name": "",
        "description": "Ordered list of indices corresponding to uncertainty estimates of transition parameters",
        "description_html": "Ordered list of indices corresponding to uncertainty estimates of transition parameters",
        "default_fmt": "%s",
        "default_units": "",
        "data_type": "str",
        "selectable": 0,
        "has_reference": 0,
        "has_error": 0,
    },
    "iref": {
        "id": 20,
        "name": "iref",
        "name_html": "References",
        "table_name": "",
        "description": "Ordered list of reference identifiers for transition parameters",
        "description_html": "Ordered list of reference identifiers for transition parameters",
        "default_fmt": "%s",
        "default_units": None,
        "data_type": "str",
        "selectable": 0,
        "has_reference": 0,
        "has_error": 0,
    },
    "deltap_air": {
        "id": 21,
        "name": "deltap_air",
        "name_html": "<em>&delta;'</em><sub>air</sub>",
        "table_name": "prm_deltap_air",
        "description": "Linear temperature dependence coefficient for air-induced pressure shift",
        "description_html": "Linear temperature dependence coefficient for air-induced pressure shift",
        "default_fmt": "%10.3e",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "n_self": {
        "id": 22,
        "name": "n_self",
        "name_html": "<em>n</em><sub>self</sub>",
        "table_name": "prm_n_self",
        "description": "Temperature exponent for the self-broadened HWHM",
        "description_html": "Temperature exponent for the self-broadened HWHM",
        "default_fmt": "%7.4f",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "delta_self": {
        "id": 23,
        "name": "delta_self",
        "name_html": "<em>&delta;</em><sub>self</sub>",
        "table_name": "prm_delta_self",
        "description": "Self-induced pressure shift, referred to p=1 atm",
        "description_html": "Self-induced pressure shift, referred to <em>p</em>=1&nbsp;atm",
        "default_fmt": "%9.6f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "deltap_self": {
        "id": 24,
        "name": "deltap_self",
        "name_html": "<em>&delta;'</em><sub>self</sub>",
        "table_name": "prm_deltap_self",
        "description": "Linear temperature dependence coefficient for self-induced pressure shift",
        "description_html": "Linear temperature dependence coefficient for self-induced pressure shift",
        "default_fmt": "%10.3e",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "sd_air": {
        "id": 28,
        "name": "SD_air",
        "name_html": "SD</sub>air</sub>",
        "table_name": "prm_sd_air",
        "description": "Speed-dependence parameter, air-broadened lines",
        "description_html": "Speed-dependence parameter, air-broadened lines",
        "default_fmt": "%9.6f",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "sd_self": {
        "id": 29,
        "name": "SD_self",
        "name_html": "SD</sub>self</sub>",
        "table_name": "prm_sd_self",
        "description": "Speed-dependence parameter, self-broadened lines",
        "description_html": "Speed-dependence parameter, self-broadened lines",
        "default_fmt": "%9.6f",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "beta_g_air": {
        "id": 30,
        "name": "beta_g_air",
        "name_html": "<em>&beta;</em><sub>G, air</sub>",
        "table_name": "prm_beta_g_air",
        "description": "Dicke narrowing parameter for the air broadened Galatry line profile",
        "description_html": "Dicke narrowing parameter for the air broadened Galatry line profile",
        "default_fmt": "%9.6f",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "y_self": {
        "id": 31,
        "name": "y_self",
        "name_html": "<em>Y</em><sub>self</sub>",
        "table_name": "prm_y_self",
        "description": "First-order (Rosenkranz) line coupling coefficient; self-broadened environment",
        "description_html": "First-order (Rosenkranz) line coupling coefficient; self-broadened environment",
        "default_fmt": "%10.3e",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "y_air": {
        "id": 32,
        "name": "y_air",
        "name_html": "<em>Y</em><sub>air</sub>",
        "table_name": "prm_y_air",
        "description": "First-order (Rosenkranz) line coupling coefficient; air-broadened environment",
        "description_html": "First-order (Rosenkranz) line coupling coefficient; air-broadened environment",
        "default_fmt": "%10.3e",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "statep": {
        "id": 33,
        "name": "statep",
        "name_html": "qns'",
        "table_name": "",
        "description": "Upper state quantum numbers",
        "description_html": "Upper state quantum numbers",
        "default_fmt": "%256s",
        "default_units": "",
        "data_type": "str",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "statepp": {
        "id": 34,
        "name": "statepp",
        "name_html": 'qns"',
        "table_name": "",
        "description": "Lower state quantum numbers",
        "description_html": "Lower state quantum numbers",
        "default_fmt": "%256s",
        "default_units": "",
        "data_type": "str",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "beta_g_self": {
        "id": 35,
        "name": "beta_g_self",
        "name_html": "<em>&beta;</em><sub>G, self</sub>",
        "table_name": "prm_beta_g_self",
        "description": "Dicke narrowing parameter for the self-broadened Galatry line profile",
        "description_html": "Dicke narrowing parameter for the self-broadened Galatry line profile",
        "default_fmt": "%9.6f",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "trans_id": {
        "id": 36,
        "name": "trans_id",
        "name_html": "Transition ID",
        "table_name": "",
        "description": "Unique integer ID of a particular transition entry in the database. (The same physical transition may have different IDs if its parameters have been revised or updated).",
        "description_html": "Unique integer ID of a particular transition entry in the database. (The same physical transition may have different IDs if its parameters have been revised or updated).",
        "default_fmt": "%12d",
        "default_units": "",
        "data_type": "int",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "par_line": {
        "id": 37,
        "name": "par_line",
        "name_html": ".par line",
        "table_name": "",
        "description": "Native 160-character formatted HITRAN line",
        "description_html": "Native 160-character formatted HITRAN line",
        "default_fmt": "%160s",
        "default_units": "",
        "data_type": "str",
        "selectable": 1,
        "has_reference": 0,
        "has_error": 0,
    },
    "gamma_h2": {
        "id": 38,
        "name": "gamma_H2",
        "name_html": "<em>&gamma;</em><sub>H2</sub> ",
        "table_name": "prm_gamma_H2",
        "description": "Lorentzian lineshape HWHM due to pressure broadening by H2 at 1 atm pressure",
        "description_html": "Lorentzian lineshape HWHM due to pressure broadening by H<sub>2</sub> at 1&nbsp;atm pressure",
        "default_fmt": "%6.4f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "n_h2": {
        "id": 39,
        "name": "n_H2",
        "name_html": "<em>n</em><sub>H2</sub>",
        "table_name": "prm_n_H2",
        "description": "Temperature exponent for the H2-broadened HWHM",
        "description_html": "Temperature exponent for the H<sub>2</sub>-broadened HWHM",
        "default_fmt": "%7.4f",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "delta_h2": {
        "id": 40,
        "name": "delta_H2",
        "name_html": "<em>&delta;</em><sub>H2</sub>",
        "table_name": "prm_delta_H2",
        "description": "Pressure shift induced by H2, referred to p=1 atm",
        "description_html": "Pressure shift induced by H<sub>2</sub>, referred to <em>p</em>=1&nbsp;atm",
        "default_fmt": "%9.6f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "deltap_h2": {
        "id": 41,
        "name": "deltap_H2",
        "name_html": "<em>&delta;'</em><sub>H2</sub>",
        "table_name": "prm_deltap_H2",
        "description": "Linear temperature dependence coefficient for H2-induced pressure shift",
        "description_html": "Linear temperature dependence coefficient for H<sub>2</sub>-induced pressure shift",
        "default_fmt": "%10.3e",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "gamma_he": {
        "id": 42,
        "name": "gamma_He",
        "name_html": "<em>&gamma;</em><sub>He</sub> ",
        "table_name": "prm_gamma_He",
        "description": "Lorentzian lineshape HWHM due to pressure broadening by He at 1 atm pressure",
        "description_html": "Lorentzian lineshape HWHM due to pressure broadening by He at 1&nbsp;atm pressure",
        "default_fmt": "%6.4f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "n_he": {
        "id": 43,
        "name": "n_He",
        "name_html": "<em>n</em><sub>He</sub>",
        "table_name": "prm_n_He",
        "description": "Temperature exponent for the He-broadened HWHM",
        "description_html": "Temperature exponent for the He-broadened HWHM",
        "default_fmt": "%7.4f",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "delta_he": {
        "id": 44,
        "name": "delta_He",
        "name_html": "<em>&delta;</em><sub>He</sub>",
        "table_name": "prm_delta_He",
        "description": "Pressure shift induced by He, referred to p=1 atm",
        "description_html": "Pressure shift induced by He, referred to <em>p</em>=1&nbsp;atm",
        "default_fmt": "%9.6f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "gamma_co2": {
        "id": 45,
        "name": "gamma_CO2",
        "name_html": "<em>&gamma;</em><sub>CO2</sub> ",
        "table_name": "prm_gamma_CO2",
        "description": "Lorentzian lineshape HWHM due to pressure broadening by CO2 at 1 atm pressure",
        "description_html": "Lorentzian lineshape HWHM due to pressure broadening by CO<sub>2</sub> at 1&nbsp;atm pressure",
        "default_fmt": "%6.4f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "n_co2": {
        "id": 46,
        "name": "n_CO2",
        "name_html": "<em>n</em><sub>CO2</sub>",
        "table_name": "prm_n_CO2",
        "description": "Temperature exponent for the CO2-broadened HWHM",
        "description_html": "Temperature exponent for the CO<sub>2</sub>-broadened HWHM",
        "default_fmt": "%7.4f",
        "default_units": "",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "delta_co2": {
        "id": 47,
        "name": "delta_CO2",
        "name_html": "<em>&delta;</em><sub>CO2</sub>",
        "table_name": "prm_delta_CO2",
        "description": "Pressure shift induced by CO2, referred to p=1 atm",
        "description_html": "Pressure shift induced by CO<sub>2</sub>, referred to <em>p</em>=1&nbsp;atm",
        "default_fmt": "%9.6f",
        "default_units": "cm-1.atm-1",
        "data_type": "float",
        "selectable": 1,
        "has_reference": 1,
        "has_error": 1,
    },
    "gamma_HT_0_self_50": {"default_fmt": "%6.4f",},
    "n_HT_self_50": {"default_fmt": "%9.6f",},
    "gamma_HT_2_self_50": {"default_fmt": "%6.4f",},
    "delta_HT_0_self_50": {"default_fmt": "%9.6f",},
    "deltap_HT_self_50": {"default_fmt": "%9.6f",},
    "delta_HT_2_self_50": {"default_fmt": "%9.6f",},
    "gamma_HT_0_self_150": {"default_fmt": "%6.4f",},
    "n_HT_self_150": {"default_fmt": "%9.6f",},
    "gamma_HT_2_self_150": {"default_fmt": "%6.4f",},
    "delta_HT_0_self_150": {"default_fmt": "%9.6f",},
    "deltap_HT_self_150": {"default_fmt": "%9.6f",},
    "delta_HT_2_self_150": {"default_fmt": "%9.6f",},
    "gamma_HT_0_self_296": {"default_fmt": "%6.4f",},
    "n_HT_self_296": {"default_fmt": "%9.6f",},
    "gamma_HT_2_self_296": {"default_fmt": "%6.4f",},
    "delta_HT_0_self_296": {"default_fmt": "%9.6f",},
    "deltap_HT_self_296": {"default_fmt": "%9.6f",},
    "delta_HT_2_self_296": {"default_fmt": "%9.6f",},
    "gamma_HT_0_self_700": {"default_fmt": "%6.4f",},
    "n_HT_self_700": {"default_fmt": "%9.6f",},
    "gamma_HT_2_self_700": {"default_fmt": "%6.4f",},
    "delta_HT_0_self_700": {"default_fmt": "%9.6f",},
    "deltap_HT_self_700": {"default_fmt": "%9.6f",},
    "delta_HT_2_self_700": {"default_fmt": "%9.6f",},
    "nu_HT_self": {"default_fmt": "%6.4f",},
    "kappa_HT_self": {"default_fmt": "%9.6f",},
    "eta_HT_self": {"default_fmt": "%9.6f",},
}


def transport2object(TransportData):
    pass


def object2transport(ObjectData):
    pass


def getFullTableAndHeaderName(TableName):
    # print('TableName=',TableName)
    fullpath_data = VARIABLES["BACKEND_DATABASE_NAME"] + "/" + TableName + ".data"
    if not os.path.isfile(fullpath_data):
        fullpath_data = VARIABLES["BACKEND_DATABASE_NAME"] + "/" + TableName + ".par"
        if not os.path.isfile(fullpath_data) and TableName != "sampletab":
            raise Exception('Lonely header "%s"' % fullpath_data)
    fullpath_header = VARIABLES["BACKEND_DATABASE_NAME"] + "/" + TableName + ".header"
    return fullpath_data, fullpath_header


def getParameterFormat(ParameterName, TableName):
    return LOCAL_TABLE_CACHE[TableName]["header"]["format"]


def getTableHeader(TableName):
    return LOCAL_TABLE_CACHE[TableName]["header"]


# RowObject = list of tuples like (name,value,format)


def addRowObject(RowObject, TableName):
    # add RowObject to TableObject in CACHE
    # check consistency first
    if [p[0] for p in RowObject] != LOCAL_TABLE_CACHE[TableName]["header"]["order"]:
        raise Exception("The row is not consistent with the table")
    for par_name, par_value, par_format in RowObject:
        LOCAL_TABLE_CACHE[TableName]["data"][par_name] += par_value
    pass


def getRowObject(RowID, TableName):
    # return RowObject from TableObject in CACHE
    RowObject = []
    for par_name in LOCAL_TABLE_CACHE[TableName]["header"]["order"]:
        par_value = LOCAL_TABLE_CACHE[TableName]["data"][par_name][RowID]
        par_format = LOCAL_TABLE_CACHE[TableName]["header"]["format"][par_name]
        RowObject.append((par_name, par_value, par_format))
    return RowObject


# INCREASE ROW COUNT


def addRowObject(RowObject, TableName):
    # print 'addRowObject: '
    # print 'RowObject: '+str(RowObject)
    # print 'TableName:'+TableName
    for par_name, par_value, par_format in RowObject:
        # print 'par_name,par_value,par_format: '+str((par_name,par_value,par_format))
        # print '>>> '+ str(LOCAL_TABLE_CACHE[TableName]['data'][par_name])
        # LOCAL_TABLE_CACHE[TableName]['data'][par_name] += [par_value]
        LOCAL_TABLE_CACHE[TableName]["data"][par_name].append(par_value)


def setRowObject(RowID, RowObject, TableName):
    number_of_rows = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]
    if RowID >= 0 and RowID < number_of_rows:
        for par_name, par_value, par_format in RowObject:
            LOCAL_TABLE_CACHE[TableName]["data"][par_name][RowID] = par_value
    else:
        # !!! XXX ATTENTION: THIS IS A TEMPORARY INSERTION XXX !!!
        LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"] += 1
        addRowObject(RowObject, TableName)


def getDefaultRowObject(TableName):
    # get a default RowObject from a table
    RowObject = []
    for par_name in LOCAL_TABLE_CACHE[TableName]["header"]["order"]:
        par_value = LOCAL_TABLE_CACHE[TableName]["header"]["default"][par_name]
        par_format = LOCAL_TABLE_CACHE[TableName]["header"]["format"][par_name]
        RowObject.append((par_name, par_value, par_format))
    return RowObject


def subsetOfRowObject(ParameterNames, RowObject):
    # return a subset of RowObject according to
    # RowObjectNew = []
    # for par_name,par_value,par_format in RowObject:
    #     if par_name in ParameterNames:
    #        RowObjectNew.append((par_name,par_value,par_format))
    # return RowObjectNew
    dct = {}
    for par_name, par_value, par_format in RowObject:
        dct[par_name] = (par_name, par_value, par_format)
    RowObjectNew = []
    for par_name in ParameterNames:
        RowObjectNew.append(dct[par_name])
    return RowObjectNew


# FORMAT_PYTHON_REGEX = '^\%([0-9]*)\.?([0-9]*)([dfs])$'
FORMAT_PYTHON_REGEX = r"^\%(\d*)(\.(\d*))?([edfsEDFS])$"

# Fortran string formatting
#  based on a pythonic format string


def formatString(par_format, par_value, lang="FORTRAN"):
    # Fortran format rules:
    #  %M.NP
    #        M - total field length (optional)
    #             (minus sign included in M)
    #        . - decimal ceparator (optional)
    #        N - number of digits after . (optional)
    #        P - [dfs] int/float/string
    # PYTHON RULE: if N is abcent, default value is 6
    regex = FORMAT_PYTHON_REGEX
    (lng, trail, lngpnt, ty) = re.search(regex, par_format).groups()
    result = par_format % par_value
    if ty.lower() in set(["f", "e"]):
        lng = int(lng) if lng else 0
        lngpnt = int(lngpnt) if lngpnt else 0
        result = par_format % par_value
        res = result.strip()
        if lng == lngpnt + 1:
            if res[0:1] == "0":
                result = "%%%ds" % lng % res[1:]
        if par_value < 0:
            if res[1:2] == "0":
                result = "%%%ds" % lng % (res[0:1] + res[2:])
    return result


def formatGetLength(fmt, lang="FORTRAN"):
    regex = FORMAT_PYTHON_REGEX


def putRowObjectToString(RowObject):
    # serialize RowObject to string
    # TODO: support different languages (C,Fortran)
    output_string = ""
    for par_name, par_value, par_format in RowObject:
        # Python formatting
        # output_string += par_format % par_value
        # Fortran formatting
        # print 'par_name,par_value,par_format: '+str((par_name,par_value,par_format))
        output_string += formatString(par_format, par_value)
    return output_string


# Parameter nicknames are hardcoded.
PARAMETER_NICKNAMES = {
    "a": "A",
    "gamma_air": "gair",
    "gp": "g",
    "local_iso_id": "I",
    "molec_id": "M",
    "sw": "S",
    "local_lower_quanta": "Q_",
    "local_upper_quanta": "Q",
    "gpp": "g_",
    "elower": "E_",
    "n_air": "nair",
    "delta_air": "dair",
    "global_upper_quanta": "V",
    "iref": "Iref",
    "line_mixing_flag": "f",
    "ierr": "ierr",
    "nu": "nu",
    "gamma_self": "gsel",
    "global_lower_quanta": "V_",
}


def putTableHeaderToString(TableName):
    output_string = ""
    regex = FORMAT_PYTHON_REGEX
    for par_name in LOCAL_TABLE_CACHE[TableName]["header"]["order"]:
        par_format = LOCAL_TABLE_CACHE[TableName]["header"]["format"][par_name]
        (lng, trail, lngpnt, ty) = re.search(regex, par_format).groups()
        fmt = "%%%ss" % lng
        try:
            par_name_short = PARAMETER_NICKNAMES[par_name]
        except:
            par_name_short = par_name
        # output_string += fmt % par_name
        output_string += (fmt % par_name_short)[: int(lng)]
    return output_string


def getRowObjectFromString(input_string, TableName):
    # restore RowObject from string, get formats and names in TableName
    # print 'getRowObjectFromString:'
    pos = 0
    RowObject = []
    # print 'Header: '+str(LOCAL_TABLE_CACHE[TableName]['header'])
    for par_name in LOCAL_TABLE_CACHE[TableName]["header"]["order"]:
        # print 'ITERATION\npos: '+str(pos) #
        # print 'par_name: '+par_name #
        par_format = LOCAL_TABLE_CACHE[TableName]["header"]["format"][par_name]
        # print 'par_format: '+par_format #
        regex = r"^\%([0-9]+)\.?[0-9]*([dfs])$"
        regex = FORMAT_PYTHON_REGEX
        # print 'par_name: '+par_name #
        (lng, trail, lngpnt, ty) = re.search(regex, par_format).groups()
        lng = int(lng)
        # print 'lng,ty:'+str((lng,ty)) #
        par_value = input_string[pos : (pos + lng)]
        # print 'par_value: '+par_value #
        if ty == "d":  # integer value
            par_value = int(par_value)
        elif ty.lower() in set(["e", "f"]):  # float value
            par_value = float(par_value)
        elif ty == "s":  # string value
            # par_value = par_value.strip() # strip spaces and tabs
            pass  # don't strip string value
        else:
            print("err1")
            raise Exception('Format "%s" is unknown' % par_format)
        RowObject.append((par_name, par_value, par_format))
        pos += lng
    # Do the same but now for extra (comma-separated) parameters
    if "extra" in set(LOCAL_TABLE_CACHE[TableName]["header"]):
        csv_chunks = input_string.split(
            LOCAL_TABLE_CACHE[TableName]["header"].get("extra_separator", ",")
        )
        # Disregard the first "column-fixed" container if it presents:
        if LOCAL_TABLE_CACHE[TableName]["header"].get("order", []):
            pos = 1
        else:
            pos = 0
        for par_name in LOCAL_TABLE_CACHE[TableName]["header"]["extra"]:
            par_format = LOCAL_TABLE_CACHE[TableName]["header"]["extra_format"][
                par_name
            ]
            regex = r"^\%([0-9]+)\.?[0-9]*([dfs])$"
            regex = FORMAT_PYTHON_REGEX
            (lng, trail, lngpnt, ty) = re.search(regex, par_format).groups()
            lng = int(lng)
            par_value = csv_chunks[pos]
            if ty == "d":  # integer value
                try:
                    par_value = int(par_value)
                except:
                    par_value = 0
            elif ty.lower() in set(["e", "f"]):  # float value
                try:
                    par_value = float(par_value)
                except:
                    par_value = 0.0
            elif ty == "s":  # string value
                # par_value = par_value.strip() # strip spaces and tabs
                pass  # don't strip string value
            else:
                print("err")
                raise Exception('Format "%s" is unknown' % par_format)
            RowObject.append((par_name, par_value, par_format))
            pos += 1
    return RowObject
    # LOCAL_TABLE_CACHE[TableName]['data'][par_name] += par_value # or append()?


# Conversion between OBJECT_FORMAT and STORAGE_FORMAT
# This will substitute putTableToStorage and getTableFromStorage


def cache2storage(TableName):
    # print 'cache2storage:'
    try:
        os.mkdir(VARIABLES["BACKEND_DATABASE_NAME"])
    except:
        pass
    fullpath_data, fullpath_header = getFullTableAndHeaderName(TableName)
    # print 'fullpath_data:'+fullpath_data
    # print 'fullpath_header'+fullpath_header
    # check if file exists and throw an exception
    # if isfile(fullpath_data): raise Exception('Table \"%s\" already exists',NewTableName)
    # if isfile(fullpath_header): raise Exception('SCHEMA IS BROKEN')
    OutfileData = open(fullpath_data, "w")
    OutfileHeader = open(fullpath_header, "w")
    # write table data
    line_count = 1
    line_number = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]
    for RowID in range(0, LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]):
        # print '%d line from %d' % (line_count,line_number)
        line_count += 1
        RowObject = getRowObject(RowID, TableName)
        # print 'RowObject:'+str(RowObject)
        raw_string = putRowObjectToString(RowObject)
        # print 'RowObject_string:'+raw_string
        OutfileData.write(raw_string + "\n")
    # write table header
    TableHeader = getTableHeader(TableName)
    OutfileHeader.write(json.dumps(TableHeader, indent=2))


def storage2cache(TableName):
    # print 'storage2cache:'
    # print('TableName',TableName)
    fullpath_data, fullpath_header = getFullTableAndHeaderName(TableName)
    InfileData = open(fullpath_data, "r")
    InfileHeader = open(fullpath_header, "r")
    # try:
    header_text = InfileHeader.read()
    try:
        Header = json.loads(header_text)
    except:
        print("HEADER:")
        print(header_text)
        raise Exception("Invalid header")
    # print 'Header:'+str(Header)
    LOCAL_TABLE_CACHE[TableName] = {}
    LOCAL_TABLE_CACHE[TableName]["header"] = Header
    LOCAL_TABLE_CACHE[TableName]["data"] = {}
    # Check if Header['order'] and Header['extra'] contain
    #  parameters with same names, raise exception if true.
    # intersct = set(Header['order']).intersection(set(Header.get('extra',[])))
    intersct = set(Header.get("order", [])).intersection(set(Header.get("extra", [])))
    if intersct:
        raise Exception("Parameters with the same names: {}".format(intersct))
    # initialize empty data to avoid problems
    glob_order = []
    glob_format = {}
    glob_default = {}
    if "order" in list(LOCAL_TABLE_CACHE[TableName]["header"].keys()):
        glob_order += LOCAL_TABLE_CACHE[TableName]["header"]["order"]
        glob_format.update(LOCAL_TABLE_CACHE[TableName]["header"]["format"])
        glob_default.update(LOCAL_TABLE_CACHE[TableName]["header"]["default"])
        for par_name in LOCAL_TABLE_CACHE[TableName]["header"]["order"]:
            LOCAL_TABLE_CACHE[TableName]["data"][par_name] = []
    if "extra" in list(LOCAL_TABLE_CACHE[TableName]["header"].keys()):
        glob_order += LOCAL_TABLE_CACHE[TableName]["header"]["extra"]
        glob_format.update(LOCAL_TABLE_CACHE[TableName]["header"]["extra_format"])
        for par_name in LOCAL_TABLE_CACHE[TableName]["header"]["extra"]:
            glob_default[par_name] = PARAMETER_META[par_name]["default_fmt"]
            LOCAL_TABLE_CACHE[TableName]["data"][par_name] = []
    line_count = 0
    # line_number = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    for line in InfileData:
        # print '%d line from %d' % (line_count,line_number)
        # print 'line: '+line #
        try:
            RowObject = getRowObjectFromString(line, TableName)
            line_count += 1
        except:
            continue
        # print 'RowObject: '+str(RowObject)
        addRowObject(RowObject, TableName)
    # except:
    #    raise Exception('TABLE FETCHING ERROR')
    LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"] = line_count
    # Delete all character-separated values, treat them as column-fixed.
    try:
        del LOCAL_TABLE_CACHE[TableName]["header"]["extra"]
        del LOCAL_TABLE_CACHE[TableName]["header"]["extra_format"]
        del LOCAL_TABLE_CACHE[TableName]["header"]["extra_separator"]
    except:
        pass
    # Update header.order/format with header.extra/format if exist.
    LOCAL_TABLE_CACHE[TableName]["header"]["order"] = glob_order
    LOCAL_TABLE_CACHE[TableName]["header"]["format"] = glob_format
    LOCAL_TABLE_CACHE[TableName]["header"]["default"] = glob_default
    InfileData.close()
    InfileHeader.close()
    print("                     Lines parsed: %d" % line_count)
    pass


# / FORMAT CONVERSION LAYER


def getTableNamesFromStorage(StorageName):
    file_names = listdir(StorageName)
    table_names = []
    for file_name in file_names:
        # search all files with "header" extensions
        # matchObject = re.search('(\w+)\.header$',file_name)
        matchObject = re.search(r"(.+)\.header$", file_name)
        if matchObject:
            # print('matchObject.group(1)=',matchObject.group(1))
            table_names.append(matchObject.group(1))
    return table_names


# FIX POSSIBLE BUG: SIMILAR NAMES OF .PAR AND .DATA FILES
# BUG FIXED BY INTRODUCING A PRIORITY:
#   *.data files have more priority than *.par files
#   See getFullTableAndHeaderName function for explanation


def scanForNewParfiles(StorageName):
    file_names = listdir(StorageName)
    headers = {}  # without extensions!
    parfiles_without_header = []
    for file_name in file_names:
        # create dictionary of unique headers
        try:
            # fname,fext = re.search('(\w+)\.(\w+)',file_name).groups()
            fname, fext = re.search(r"(.+)\.(\w+)", file_name).groups()
        except:
            continue
        if fext == "header":
            headers[fname] = True
    for file_name in file_names:
        # check if extension is 'par' and the header is absent
        try:
            # fname,fext = re.search('(\w+)\.(\w+)',file_name).groups()
            fname, fext = re.search(r"(.+)\.(\w+)", file_name).groups()
        except:
            continue
        if fext == "par" and fname not in headers:
            parfiles_without_header.append(fname)
    return parfiles_without_header


def createHeader(TableName):
    fname = TableName + ".header"
    fp = open(VARIABLES["BACKEND_DATABASE_NAME"] + "/" + fname, "w")
    if os.path.isfile(TableName):
        raise Exception('File "%s" already exists!' % fname)
    fp.write(json.dumps(HITRAN_DEFAULT_HEADER, indent=2))
    fp.close()


def loadCache():
    # print 'loadCache:'
    print("Using " + VARIABLES["BACKEND_DATABASE_NAME"] + "\n")
    LOCAL_TABLE_CACHE = {}  # ?????
    table_names = getTableNamesFromStorage(VARIABLES["BACKEND_DATABASE_NAME"])
    # print('table_names=',table_names)
    parfiles_without_header = scanForNewParfiles(VARIABLES["BACKEND_DATABASE_NAME"])
    # create headers for new parfiles
    for tab_name in parfiles_without_header:
        # get name without 'par' extension
        createHeader(tab_name)
        table_names.append(tab_name)
    for TableName in table_names:
        print(TableName)
        storage2cache(TableName)


def saveCache():
    # print 'saveCache:'
    try:
        # delete query buffer
        del LOCAL_TABLE_CACHE[QUERY_BUFFER]
    except:
        pass
    for TableName in LOCAL_TABLE_CACHE:
        print(TableName)
        cache2storage(TableName)


# DB backend level, start transaction


def databaseBegin(db=None):
    if db:
        VARIABLES["BACKEND_DATABASE_NAME"] = db
    else:
        VARIABLES["BACKEND_DATABASE_NAME"] = BACKEND_DATABASE_NAME_DEFAULT
    # print 'databaseBegin:'
    # print(os.path.isdir("/home/el"))
    # print(os.path.exists("/home/el/myfile.txt"))
    if not os.path.exists(VARIABLES["BACKEND_DATABASE_NAME"]):
        os.mkdir(VARIABLES["BACKEND_DATABASE_NAME"])
    loadCache()


# DB backend level, end transaction


def databaseCommit():
    # print 'databaseCommit:'
    saveCache()


# def saveCache():
#    for TableName in LOCAL_TABLE_CACHE.keys():
#        putTableToStorage(TableName)

# ----------------------------------------------------
# ----------------------------------------------------
# CONDITIONS
# ----------------------------------------------------
# ----------------------------------------------------
# ----------------------------------------------------
# hierarchic query.condition language:
# Conditions: CONS = ('and', ('=','p1','p2'), ('<','p1',13))
# String literals are distinguished from variable names
#  by using the operation ('STRING','some_string')
# ----------------------------------------------------


# necessary conditions for hitranonline:
SAMPLE_CONDITIONS = (
    "AND",
    ("SET", "internal_iso_id", [1, 2, 3, 4, 5, 6]),
    (">=", "nu", 0),
    ("<=", "nu", 100),
)

# sample hitranonline protocol
# http://hitran.cloudapp.net/lbl/5?output_format_id=1&iso_ids_list=5&numin=0&numax=100&access=api&key=e20e4bd3-e12c-4931-99e0-4c06e88536bd

CONDITION_OPERATIONS = set(
    [
        "AND",
        "OR",
        "NOT",
        "RANGE",
        "IN",
        "<",
        ">",
        "<=",
        ">=",
        "==",
        "!=",
        "LIKE",
        "STR",
        "+",
        "-",
        "*",
        "/",
        "MATCH",
        "SEARCH",
        "FINDALL",
    ]
)

# Operations used in Condition verification
# Basic scheme: operationXXX(args),
# where args - list/array of arguments (>=1)


def operationAND(args):
    # any number if arguments
    for arg in args:
        if not arg:
            return False
    return True


def operationOR(args):
    # any number of arguments
    for arg in args:
        if arg:
            return True
    return False


def operationNOT(arg):
    # one argument
    return not arg


def operationRANGE(x, x_min, x_max):
    return x_min <= x <= x_max


def operationSUBSET(arg1, arg2):
    # True if arg1 is subset of arg2
    # arg1 is an element
    # arg2 is a set
    return arg1 in arg2


def operationLESS(args):
    # any number of args
    for i in range(1, len(args)):
        if args[i - 1] >= args[i]:
            return False
    return True


def operationMORE(args):
    # any number of args
    for i in range(1, len(args)):
        if args[i - 1] <= args[i]:
            return False
    return True


def operationLESSOREQUAL(args):
    # any number of args
    for i in range(1, len(args)):
        if args[i - 1] > args[i]:
            return False
    return True


def operationMOREOREQUAL(args):
    # any number of args
    for i in range(1, len(args)):
        if args[i - 1] < args[i]:
            return False
    return True


def operationEQUAL(args):
    # any number of args
    for i in range(1, len(args)):
        if args[i] != args[i - 1]:
            return False
    return True


def operationNOTEQUAL(arg1, arg2):
    return arg1 != arg2


def operationSUM(args):
    # any numbers of arguments
    if type(args[0]) in set([int, float]):
        result = 0
    elif type(args[0]) in set([str, six.text_type]):
        result = ""
    else:
        raise Exception("SUM error: unknown arg type")
    for arg in args:
        result += arg
    return result


def operationDIFF(arg1, arg2):
    return arg1 - arg2


def operationMUL(args):
    # any numbers of arguments
    if type(args[0]) in set([int, float]):
        result = 1
    else:
        raise Exception("MUL error: unknown arg type")
    for arg in args:
        result *= arg
    return result


def operationDIV(arg1, arg2):
    return arg1 / arg2


def operationSTR(arg):
    # transform arg to str
    if type(arg) != str:
        raise Exception("Type mismatch: STR")
    return arg


def operationSET(arg):
    # transform arg to list
    if type(arg) not in set([list, tuple, set]):
        raise Exception("Type mismatch: SET")
    return list(arg)


def operationMATCH(arg1, arg2):
    # Match regex (arg1) and string (arg2)
    # return bool(re.match(arg1,arg2)) # works wrong
    return bool(re.search(arg1, arg2))


def operationSEARCH(arg1, arg2):
    # Search regex (arg1) in string (arg2)
    # Output list of entries
    group = re.search(arg1, arg2).groups()
    result = []
    for item in group:
        result.append(("STR", item))
    return result


def operationFINDALL(arg1, arg2):
    # Search all groups of a regex
    # Output a list of groups of entries
    # XXX: If a group has more than 1 entry,
    #    there could be potential problems
    list_of_groups = re.findall(arg1, arg2)
    result = []
    for item in list_of_groups:
        result.append(("STR", item))
    return result


def operationLIST(args):
    # args is a list: do nothing (almost)
    return list(args)


# /operations

# def parse(Conditions):
#    pass

# GROUPING ----------------------------------------------


GROUP_INDEX = {}
# GROUP_INDEX has the following structure:
#  GROUP_INDEX[KEY] = VALUE
#    KEY = table line values
#    VALUE = {'FUNCTIONS':DICT,'FLAG':LOGICAL,'ROWID':INTEGER}
#      FUNCTIONS = {'FUNC_NAME':DICT}
#            FUNC_NAME = {'FLAG':LOGICAL,'NAME':STRING}

# name and default value
GROUP_FUNCTION_NAMES = {
    "COUNT": 0,
    "SUM": 0,
    "MUL": 1,
    "AVG": 0,
    "MIN": +1e100,
    "MAX": -1e100,
    "SSQ": 0,
}


def clearGroupIndex():
    # GROUP_INDEX = {}
    for key in GROUP_INDEX.keys():
        del GROUP_INDEX[key]


def getValueFromGroupIndex(GroupIndexKey, FunctionName):
    # If no such index_key, create it and return a value
    if FunctionName not in GROUP_FUNCTION_NAMES:
        raise Exception('No such function "%s"' % FunctionName)
    # In the case if NewRowObjectDefault is requested
    if not GroupIndexKey:
        return GROUP_FUNCTION_NAMES[FunctionName]
    if FunctionName not in GROUP_INDEX[GroupIndexKey]["FUNCTIONS"]:
        GROUP_INDEX[GroupIndexKey]["FUNCTIONS"][FunctionName] = {}
        GROUP_INDEX[GroupIndexKey]["FUNCTIONS"][FunctionName]["FLAG"] = True
        GROUP_INDEX[GroupIndexKey]["FUNCTIONS"][FunctionName][
            "VALUE"
        ] = GROUP_FUNCTION_NAMES[FunctionName]
    return GROUP_INDEX[GroupIndexKey]["FUNCTIONS"][FunctionName]["VALUE"]


def setValueToGroupIndex(GroupIndexKey, FunctionName, Value):
    GROUP_INDEX[GroupIndexKey]["FUNCTIONS"][FunctionName]["VALUE"] = Value


def initializeGroup(GroupIndexKey):
    if GroupIndexKey not in GROUP_INDEX:
        print("GROUP_DESC[COUNT]=" + str(GROUP_DESC["COUNT"]))
        GROUP_INDEX[GroupIndexKey] = {}
        GROUP_INDEX[GroupIndexKey]["FUNCTIONS"] = {}
        GROUP_INDEX[GroupIndexKey]["ROWID"] = len(GROUP_INDEX) - 1
    for FunctionName in GROUP_FUNCTION_NAMES:
        # initialize function flags (UpdateFlag)
        if FunctionName in GROUP_INDEX[GroupIndexKey]["FUNCTIONS"]:
            GROUP_INDEX[GroupIndexKey]["FUNCTIONS"][FunctionName]["FLAG"] = True
    print("initializeGroup: GROUP_INDEX=" + str(GROUP_INDEX))


def groupCOUNT(GroupIndexKey):
    FunctionName = "COUNT"
    Value = getValueFromGroupIndex(GroupIndexKey, FunctionName)
    if GroupIndexKey:
        if GROUP_INDEX[GroupIndexKey]["FUNCTIONS"][FunctionName]["FLAG"]:
            GROUP_INDEX[GroupIndexKey]["FUNCTIONS"][FunctionName]["FLAG"] = False
            Value = Value + 1
            setValueToGroupIndex(GroupIndexKey, FunctionName, Value)
    return Value


def groupSUM():
    pass


def groupMUL():
    pass


def groupAVG():
    pass


def groupMIN():
    pass


def groupMAX():
    pass


def groupSSQ():
    pass


OPERATORS = {  # List
    "LIST": lambda args: operationLIST(args),
    # And
    "&": lambda args: operationAND(args),
    "&&": lambda args: operationAND(args),
    "AND": lambda args: operationAND(args),
    # Or
    "|": lambda args: operationOR(args),
    "||": lambda args: operationOR(args),
    "OR": lambda args: operationOR(args),
    # Not
    "!": lambda args: operationNOT(args[0]),
    "NOT": lambda args: operationNOT(args[0]),
    # Between
    "RANGE": lambda args: operationRANGE(args[0], args[1], args[2]),
    "BETWEEN": lambda args: operationRANGE(args[0], args[1], args[2]),
    # Subset
    "IN": lambda args: operationSUBSET(args[0], args[1]),
    "SUBSET": lambda args: operationSUBSET(args[0], args[1]),
    # Less
    "<": lambda args: operationLESS(args),
    "LESS": lambda args: operationLESS(args),
    "LT": lambda args: operationLESS(args),
    # More
    ">": lambda args: operationMORE(args),
    "MORE": lambda args: operationMORE(args),
    "MT": lambda args: operationMORE(args),
    # Less or equal
    "<=": lambda args: operationLESSOREQUAL(args),
    "LESSOREQUAL": lambda args: operationLESSOREQUAL(args),
    "LTE": lambda args: operationLESSOREQUAL(args),
    # More or equal
    ">=": lambda args: operationMOREOREQUAL(args),
    "MOREOREQUAL": lambda args: operationMOREOREQUAL(args),
    "MTE": lambda args: operationMOREOREQUAL(args),
    # Equal
    "=": lambda args: operationEQUAL(args),
    "==": lambda args: operationEQUAL(args),
    "EQ": lambda args: operationEQUAL(args),
    "EQUAL": lambda args: operationEQUAL(args),
    "EQUALS": lambda args: operationEQUAL(args),
    # Not equal
    "!=": lambda args: operationNOTEQUAL(args[0], args[1]),
    "<>": lambda args: operationNOTEQUAL(args[0], args[1]),
    "~=": lambda args: operationNOTEQUAL(args[0], args[1]),
    "NE": lambda args: operationNOTEQUAL(args[0], args[1]),
    "NOTEQUAL": lambda args: operationNOTEQUAL(args[0], args[1]),
    # Plus
    "+": lambda args: operationSUM(args),
    "SUM": lambda args: operationSUM(args),
    # Minus
    "-": lambda args: operationDIFF(args[0], args[1]),
    "DIFF": lambda args: operationDIFF(args[0], args[1]),
    # Mul
    "*": lambda args: operationMUL(args),
    "MUL": lambda args: operationMUL(args),
    # Div
    "/": lambda args: operationDIV(args[0], args[1]),
    "DIV": lambda args: operationDIV(args[0], args[1]),
    # Regexp match
    "MATCH": lambda args: operationMATCH(args[0], args[1]),
    "LIKE": lambda args: operationMATCH(args[0], args[1]),
    # Regexp search
    "SEARCH": lambda args: operationSEARCH(args[0], args[1]),
    # Regexp findal
    "FINDALL": lambda args: operationFINDALL(args[0], args[1]),
    # Group count
    "COUNT": lambda args: groupCOUNT(GroupIndexKey),
}

# new evaluateExpression function,
#  accounting for groups
"""
def evaluateExpression(root,VarDictionary,GroupIndexKey=None):
    # input = local tree root
    # XXX: this could be very slow due to passing
    #      every time VarDictionary as a parameter
    # Two special cases: 1) root=varname
    #                    2) root=list/tuple
    # These cases must be processed in a separate way
    if type(root) in set([list,tuple]):
       # root is not a leaf
       head = root[0].upper()
       # string constants are treated specially
       if head in set(['STR','STRING']): # one arg
          return operationSTR(root[1])
       elif head in set(['SET']):
          return operationSET(root[1])
       tail = root[1:]
       args = []
       # evaluate arguments recursively
       for element in tail: # resolve tree by recursion
           args.append(evaluateExpression(element,VarDictionary,GroupIndexKey))
       # call functions with evaluated arguments
       if head in set(['LIST']): # list arg
          return operationLIST(args)
       elif head in set(['&','&&','AND']): # many args 
          return operationAND(args)
       elif head in set(['|','||','OR']): # many args
          return operationOR(args)
       elif head in set(['!','NOT']): # one args
          return operationNOT(args[0])
       elif head in set(['RANGE','BETWEEN']): # three args
          return operationRANGE(args[0],args[1],args[2])
       elif head in set(['IN','SUBSET']): # two args
          return operationSUBSET(args[0],args[1])
       elif head in set(['<','LESS','LT']): # many args
          return operationLESS(args)
       elif head in set(['>','MORE','MT']): # many args
          return operationMORE(args)
       elif head in set(['<=','LESSOREQUAL','LTE']): # many args
          return operationLESSOREQUAL(args)
       elif head in set(['>=','MOREOREQUAL','MTE']): # many args
          return operationMOREOREQUAL(args)
       elif head in set(['=','==','EQ','EQUAL','EQUALS']): # many args
          return operationEQUAL(args)
       elif head in set(['!=','<>','~=','NE','NOTEQUAL']): # two args
          return operationNOTEQUAL(args[0],args[1])
       elif head in set(['+','SUM']): # many args
          return operationSUM(args)
       elif head in set(['-','DIFF']): # two args
          return operationDIFF(args[0],args[1])
       elif head in set(['*','MUL']): # many args
          return operationMUL(args)
       elif head in set(['/','DIV']): # two args
          return operationDIV(args[0],args[1])
       elif head in set(['MATCH','LIKE']): # two args
          return operationMATCH(args[0],args[1])
       elif head in set(['SEARCH']): # two args
          return operationSEARCH(args[0],args[1])
       elif head in set(['FINDALL']): # two args
          return operationFINDALL(args[0],args[1])
       # --- GROUPING OPERATIONS ---
       elif head in set(['COUNT']):
          return groupCOUNT(GroupIndexKey)
       else:
          raise Exception('Unknown operator: %s' % root[0])
    elif type(root)==str:
       # root is a par_name
       return VarDictionary[root]
    else: 
       # root is a non-string constant
       return root
"""


def evaluateExpression(root, VarDictionary, GroupIndexKey=None):
    # input = local tree root
    # XXX: this could be very slow due to passing
    #      every time VarDictionary as a parameter
    # Two special cases: 1) root=varname
    #                    2) root=list/tuple
    # These cases must be processed in a separate way
    if type(root) in set([list, tuple]):
        # root is not a leaf
        head = root[0].upper()
        # string constants are treated specially
        if head in set(["STR", "STRING"]):  # one arg
            return operationSTR(root[1])
        elif head in set(["SET"]):
            return operationSET(root[1])
        tail = root[1:]
        args = []
        # evaluate arguments recursively
        for element in tail:  # resolve tree by recursion
            args.append(evaluateExpression(element, VarDictionary, GroupIndexKey))
        # call functions with evaluated arguments
        try:
            return OPERATORS[head](args)
        except KeyError:
            raise Exception("Unknown operator: %s" % head)
    elif type(root) == str:
        # root is a par_name
        return VarDictionary[root]
    else:
        # root is a non-string constant
        return root


def getVarDictionary(RowObject):
    # get VarDict from RowObject
    # VarDict: par_name => par_value
    VarDictionary = {}
    for par_name, par_value, par_format in RowObject:
        VarDictionary[par_name] = par_value
    return VarDictionary


def checkRowObject(RowObject, Conditions, VarDictionary):
    # VarDictionary = getVarDictionary(RowObject)
    if Conditions:
        Flag = evaluateExpression(Conditions, VarDictionary)
    else:
        Flag = True
    return Flag


# ----------------------------------------------------
# /CONDITIONS
# ----------------------------------------------------


# ----------------------------------------------------
# PARAMETER NAMES (includeing creation of new ones)
# ----------------------------------------------------

# Bind an expression to a new parameter
#   in a form: ('BIND','new_par',('some_exp',...))
def operationBIND(parname, Expression, VarDictionary):
    pass


# This section is for more detailed processing of parlists.

# Table creation must include not only subsets of
#   existing parameters, but also new parameters
#   derived from functions on a special prefix language
# For this reason subsetOfRowObject(..) must be substituted
#   by newRowObject(ParameterNames,RowObject)

# For parsing use the function evaluateExpression

# Get names from expression.
#  Must merge this one with evaluateExrpression.
# This is VERY LIMITED version of what will be
#  when make the language parser is implemented.
# For more ideas and info see LANGUAGE_REFERENCE

# more advansed version of expression evaluator


def evaluateExpressionPAR(ParameterNames, VarDictionary=None):
    # RETURN: 1) Upper-level Expression names
    #         2) Upper-level Expression values
    # Is it reasonable to pass a Context to every parse function?
    # For now the function does the following:
    #   1) iterates through all UPPER-LEVEL list elements
    #   2) if element is a parname: return parname
    #      if element is an BIND expression: return bind name
    #              (see operationBIND)
    #   3) if element is an anonymous expression: return #N(=1,2,3...)
    # N.B. Binds can be only on the 0-th level of Expression
    pass


def getContextFormat(RowObject):
    # Get context format from the whole RowObject
    ContextFormat = {}
    for par_name, par_value, par_format in RowObject:
        ContextFormat[par_name] = par_format
    return ContextFormat


def getDefaultFormat(Type):
    if Type is int:
        return "%10d"
    elif Type is float:
        return "%25.15E"
    elif Type is str:
        return "%20s"
    elif Type is bool:
        return "%2d"
    else:
        raise Exception("Unknown type")


def getDefaultValue(Type):
    if Type is int:
        return 0
    elif Type is float:
        return 0.0
    elif Type is str:
        return ""
    elif Type is bool:
        return False
    else:
        raise Exception("Unknown type")


# VarDictionary = Context (this name is more suitable)

# GroupIndexKey is a key to special structure/dictionary GROUP_INDEX.
# GROUP_INDEX contains information needed to calculate streamed group functions
#  such as COUNT, AVG, MIN, MAX etc...


def newRowObject(
    ParameterNames, RowObject, VarDictionary, ContextFormat, GroupIndexKey=None
):
    # Return a subset of RowObject according to
    # ParameterNames include either parnames
    #  or expressions containing parnames literals
    # ContextFormat contains format for ParNames
    anoncount = 0
    RowObjectNew = []
    for expr in ParameterNames:
        if type(expr) in set([list, tuple]):  # bind
            head = expr[0]
            if head in set(["let", "bind", "LET", "BIND"]):
                par_name = expr[1]
                par_expr = expr[2]
            else:
                par_name = "#%d" % anoncount
                anoncount += 1
                par_expr = expr
            par_value = evaluateExpression(par_expr, VarDictionary, GroupIndexKey)
            try:
                par_format = expr[3]
            except:
                par_format = getDefaultFormat(type(par_value))
        else:  # parname
            par_name = expr
            par_value = VarDictionary[par_name]
            par_format = ContextFormat[par_name]
        RowObjectNew.append((par_name, par_value, par_format))
    return RowObjectNew


# ----------------------------------------------------
# /PARAMETER NAMES
# ----------------------------------------------------


# ----------------------------------------------------
# OPERATIONS ON TABLES
# ----------------------------------------------------

QUERY_BUFFER = "__BUFFER__"


def getTableList():
    return list(LOCAL_TABLE_CACHE.keys())


def describeTable(TableName):
    """
    INPUT PARAMETERS: 
        TableName: name of the table to describe
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Print information about table, including 
        parameter names, formats and wavenumber range.
    ---
    EXAMPLE OF USAGE:
        describeTable('sampletab')
    ---
    """
    print("-----------------------------------------")
    print(TableName + " summary:")
    try:
        print("-----------------------------------------")
        print("Comment: \n" + LOCAL_TABLE_CACHE[TableName]["header"]["comment"])
    except:
        pass
    print(
        "Number of rows: "
        + str(LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"])
    )
    print("Table type: " + str(LOCAL_TABLE_CACHE[TableName]["header"]["table_type"]))
    print("-----------------------------------------")
    print("            PAR_NAME           PAR_FORMAT")
    print("")
    for par_name in LOCAL_TABLE_CACHE[TableName]["header"]["order"]:
        par_format = LOCAL_TABLE_CACHE[TableName]["header"]["format"][par_name]
        print("%20s %20s" % (par_name, par_format))
    print("-----------------------------------------")


# Write a table to File or STDOUT


def outputTable(TableName, Conditions=None, File=None, Header=True):
    # Display or record table with condition checking
    if File:
        Header = False
        OutputFile = open(File, "w")
    if Header:
        headstr = putTableHeaderToString(TableName)
        if File:
            OutputFile.write(headstr)
        else:
            print(headstr)
    for RowID in range(0, LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]):
        RowObject = getRowObject(RowID, TableName)
        VarDictionary = getVarDictionary(RowObject)
        VarDictionary["LineNumber"] = RowID
        if not checkRowObject(RowObject, Conditions, VarDictionary):
            continue
        raw_string = putRowObjectToString(RowObject)
        if File:
            OutputFile.write(raw_string + "\n")
        else:
            print(raw_string)


# Create table "prototype-based" way


def createTable(TableName, RowObjectDefault):
    # create a Table based on a RowObjectDefault
    LOCAL_TABLE_CACHE[TableName] = {}
    header_order = []
    header_format = {}
    header_default = {}
    data = {}
    for par_name, par_value, par_format in RowObjectDefault:
        header_order.append(par_name)
        header_format[par_name] = par_format
        header_default[par_name] = par_value
        data[par_name] = []
    # header_order = tuple(header_order) # XXX ?
    LOCAL_TABLE_CACHE[TableName]["header"] = {}
    LOCAL_TABLE_CACHE[TableName]["header"]["order"] = header_order
    LOCAL_TABLE_CACHE[TableName]["header"]["format"] = header_format
    LOCAL_TABLE_CACHE[TableName]["header"]["default"] = header_default
    LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"] = 0
    LOCAL_TABLE_CACHE[TableName]["header"]["size_in_bytes"] = 0
    LOCAL_TABLE_CACHE[TableName]["header"]["table_name"] = TableName
    LOCAL_TABLE_CACHE[TableName]["header"]["table_type"] = "column-fixed"
    LOCAL_TABLE_CACHE[TableName]["data"] = data


# simple "drop table" capability
def dropTable(TableName):
    """
    INPUT PARAMETERS: 
        TableName:  name of the table to delete
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Deletes a table from local database.
    ---
    EXAMPLE OF USAGE:
        dropTable('some_dummy_table')
    ---
    """
    # delete Table from both Cache and Storage
    try:
        # LOCAL_TABLE_CACHE[TableName] = {}
        del LOCAL_TABLE_CACHE[TableName]
    except:
        pass
    # delete from storage
    pass  # TODO


# Returns a column corresponding to parameter name


def getColumn(TableName, ParameterName):
    """
    INPUT PARAMETERS: 
        TableName:      source table name     (required)
        ParameterName:  name of column to get (required)
    OUTPUT PARAMETERS: 
        ColumnData:     list of values from specified column 
    ---
    DESCRIPTION:
        Returns a column with a name ParameterName from
        table TableName. Column is returned as a list of values.
    ---
    EXAMPLE OF USAGE:
        p1 = getColumn('sampletab','p1')
    ---
    """
    return LOCAL_TABLE_CACHE[TableName]["data"][ParameterName]


# Returns a list of columns corresponding to parameter names


def getColumns(TableName, ParameterNames):
    """
    INPUT PARAMETERS: 
        TableName:       source table name           (required)
        ParameterNames:  list of column names to get (required)
    OUTPUT PARAMETERS: 
        ListColumnData:   tuple of lists of values from specified column 
    ---
    DESCRIPTION:
        Returns columns with a names in ParameterNames from
        table TableName. Columns are returned as a tuple of lists.
    ---
    EXAMPLE OF USAGE:
        p1,p2,p3 = getColumns('sampletab',('p1','p2','p3'))
    ---
    """
    Columns = []
    for par_name in ParameterNames:
        Columns.append(LOCAL_TABLE_CACHE[TableName]["data"][par_name])
    return Columns


def addColumn(
    TableName,
    ParameterName,
    Before=None,
    Expression=None,
    Type=None,
    Default=None,
    Format=None,
):
    if ParameterName in LOCAL_TABLE_CACHE[TableName]["header"]["format"]:
        raise Exception('Column "%s" already exists' % ParameterName)
    if not Type:
        Type = float
    if not Default:
        Default = getDefaultValue(Type)
    if not Format:
        Format = getDefaultFormat(Type)
    number_of_rows = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]
    # Mess with data
    if not Expression:
        LOCAL_TABLE_CACHE[TableName]["data"][ParameterName] = [
            Default for i in range(0, number_of_rows)
        ]
    else:
        data = []
        for RowID in range(0, number_of_rows):
            RowObject = getRowObject(RowID, TableName)
            VarDictionary = getVarDictionary(RowObject)
            VarDictionary["LineNumber"] = RowID
            par_value = evaluateExpression(Expression, VarDictionary)
            data.append(par_value)
            LOCAL_TABLE_CACHE[TableName]["data"][ParameterName] = data
    # Mess with header
    header_order = LOCAL_TABLE_CACHE[TableName]["header"]["order"]
    if not Before:
        header_order.append(ParameterName)
    else:
        # i = 0
        # for par_name in header_order:
        #    if par_name == Before: break
        #    i += 1
        i = header_order.index(Before)
        header_order = header_order[:i] + [ParameterName,] + header_order[i:]
    LOCAL_TABLE_CACHE[TableName]["header"]["order"] = header_order
    LOCAL_TABLE_CACHE[TableName]["header"]["format"][ParameterName] = Format
    LOCAL_TABLE_CACHE[TableName]["header"]["default"][ParameterName] = Default


def deleteColumn(TableName, ParameterName):
    if ParameterName not in LOCAL_TABLE_CACHE[TableName]["header"]["format"]:
        raise Exception('No such column "%s"' % ParameterName)
    # Mess with data
    i = LOCAL_TABLE_CACHE[TableName]["header"]["order"].index(ParameterName)
    del LOCAL_TABLE_CACHE[TableName]["header"]["order"][i]
    del LOCAL_TABLE_CACHE[TableName]["header"]["format"][ParameterName]
    del LOCAL_TABLE_CACHE[TableName]["header"]["default"][ParameterName]
    if not LOCAL_TABLE_CACHE[TableName]["header"]["order"]:
        LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"] = 0
    # Mess with header
    del LOCAL_TABLE_CACHE[TableName]["data"][ParameterName]


def deleteColumns(TableName, ParameterNames):
    if type(ParameterNames) not in set([list, tuple, set]):
        ParameterNames = [ParameterNames]
    for ParameterName in ParameterNames:
        deleteColumn(TableName, ParameterName)


def renameColumn(TableName, OldParameterName, NewParameterName):
    pass


def insertRow():
    pass


def deleteRows(TableName, ParameterNames, Conditions):
    pass


# select from table to another table


def selectInto(DestinationTableName, TableName, ParameterNames, Conditions):
    # TableName must refer to an existing table in cache!!
    # Conditions = Restrictables in specific format
    # Sample conditions: cond = {'par1':{'range',[b_lo,b_hi]},'par2':b}
    # return structure similar to TableObject and put it to QUERY_BUFFER
    # if ParameterNames is '*' then all parameters are used
    # table_columns = LOCAL_TABLE_CACHE[TableName]['data'].keys()
    # table_length = len(TableObject['header']['number_of_rows'])
    # if ParameterNames=='*':
    #   ParameterNames = table_columns
    # check if Conditions contain elements which are not in the TableObject
    # condition_variables = getConditionVariables(Conditions)
    # strange_pars = set(condition_variables)-set(table_variables)
    # if strange_pars:
    #   raise Exception('The following parameters are not in the table \"%s\"' % (TableName,list(strange_pars)))
    # do full scan each time
    if DestinationTableName == TableName:
        raise Exception("Selecting into source table is forbidden")
    table_length = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]
    row_count = 0
    for RowID in range(0, table_length):
        RowObject = getRowObject(RowID, TableName)
        VarDictionary = getVarDictionary(RowObject)
        VarDictionary["LineNumber"] = RowID
        ContextFormat = getContextFormat(RowObject)
        RowObjectNew = newRowObject(
            ParameterNames, RowObject, VarDictionary, ContextFormat
        )
        if checkRowObject(RowObject, Conditions, VarDictionary):
            addRowObject(RowObjectNew, DestinationTableName)
            row_count += 1
    LOCAL_TABLE_CACHE[DestinationTableName]["header"]["number_of_rows"] += row_count


def length(TableName):
    tab_len = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]
    # print(str(tab_len)+' rows in '+TableName)
    return tab_len


# Select parameters from a table with certain conditions.
# Parameters can be the names or expressions.
# Conditions contain a list of expressions in a special language.
# Set Output to False to suppress output
# Set File=FileName to redirect output to a file.


def select(
    TableName,
    DestinationTableName=QUERY_BUFFER,
    ParameterNames=None,
    Conditions=None,
    Output=True,
    File=None,
):
    """
    INPUT PARAMETERS: 
        TableName:            name of source table              (required)
        DestinationTableName: name of resulting table           (optional)
        ParameterNames:       list of parameters or expressions (optional)
        Conditions:           list of logincal expressions      (optional)
        Output:   enable (True) or suppress (False) text output (optional)
        File:     enable (True) or suppress (False) file output (optional)
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Select or filter the data in some table 
        either to standard output or to file (if specified)
    ---
    EXAMPLE OF USAGE:
        select('sampletab',DestinationTableName='outtab',ParameterNames=(p1,p2),
                Conditions=(('and',('>=','p1',1),('<',('*','p1','p2'),20))))
        Conditions means (p1>=1 and p1*p2<20)
    ---
    """
    # TODO: Variables defined in ParameterNames ('LET') MUST BE VISIBLE IN Conditions !!
    # check if table exists
    if TableName not in list(LOCAL_TABLE_CACHE.keys()):
        raise Exception(
            "%s: no such table. Check tableList() for more info." % TableName
        )
    if not ParameterNames:
        ParameterNames = LOCAL_TABLE_CACHE[TableName]["header"]["order"]
    # clear QUERY_BUFFER for the new result
    LOCAL_TABLE_CACHE[DestinationTableName] = {}
    RowObjectDefault = getDefaultRowObject(TableName)
    VarDictionary = getVarDictionary(RowObjectDefault)
    ContextFormat = getContextFormat(RowObjectDefault)
    RowObjectDefaultNew = newRowObject(
        ParameterNames, RowObjectDefault, VarDictionary, ContextFormat
    )
    dropTable(DestinationTableName)  # redundant
    createTable(DestinationTableName, RowObjectDefaultNew)
    selectInto(DestinationTableName, TableName, ParameterNames, Conditions)
    if DestinationTableName != QUERY_BUFFER:
        if File:
            outputTable(DestinationTableName, File=File)
    elif Output:
        outputTable(DestinationTableName, File=File)


# SORTING ===========================================================


def arrangeTable(TableName, DestinationTableName=None, RowIDList=None):
    # print 'AT/'
    # print 'AT: RowIDList = '+str(RowIDList)
    # make a subset of table rows according to RowIDList
    if not DestinationTableName:
        DestinationTablename = TableName
    if DestinationTableName != TableName:
        dropTable(DestinationTableName)
        LOCAL_TABLE_CACHE[DestinationTableName]["header"] = LOCAL_TABLE_CACHE[
            TableName
        ]["header"]
        LOCAL_TABLE_CACHE[DestinationTableName]["data"] = {}
    LOCAL_TABLE_CACHE[DestinationTableName]["header"]["number_of_rows"] = len(RowIDList)
    # print 'AT: RowIDList = '+str(RowIDList)
    for par_name in LOCAL_TABLE_CACHE[DestinationTableName]["header"]["order"]:
        par_data = LOCAL_TABLE_CACHE[TableName]["data"][par_name]
        LOCAL_TABLE_CACHE[DestinationTableName]["data"][par_name] = [
            par_data[i] for i in RowIDList
        ]


def compareLESS(RowObject1, RowObject2, ParameterNames):
    # print 'CL/'
    # arg1 and arg2 are RowObjects
    # Compare them according to ParameterNames
    # Simple validity check:
    # if len(arg1) != len(arg2):
    #   raise Exception('Arguments have different lengths')
    # RowObject1Subset = subsetOfRowObject(ParameterNames,RowObject1)
    # RowObject2Subset = subsetOfRowObject(ParameterNames,RowObject2)
    # return RowObject1Subset < RowObject2Subset
    row1 = []
    row2 = []
    # n = len(RowObject1)
    # for i in range(0,n):
    #    par_name1 = RowObject1[i][0]
    #    if par_name1 in ParameterNames:
    #       par_value1 = RowObject1[i][1]
    #       par_value2 = RowObject2[i][1]
    #       row1 += [par_value1]
    #       row2 += [par_value2]
    VarDictionary1 = getVarDictionary(RowObject1)
    VarDictionary2 = getVarDictionary(RowObject2)
    for par_name in ParameterNames:
        par_value1 = VarDictionary1[par_name]
        par_value2 = VarDictionary2[par_name]
        row1 += [par_value1]
        row2 += [par_value2]
    Flag = row1 < row2
    # print 'CL: row1 = '+str(row1)
    # print 'CL: row2 = '+str(row2)
    # print 'CL: Flag = '+str(Flag)
    return Flag


def quickSort(index, TableName, ParameterNames, Accending=True):
    # print ''
    # print 'QS/'
    # print 'QS: index = '+str(index)
    # print index
    # ParameterNames: names of parameters which are
    #  taking part in the sorting
    if index == []:
        return []
    else:
        # pivot = lst[0]
        # lesser = quickSort([x for x in lst[1:] if x < pivot])
        # greater = quickSort([x for x in lst[1:] if x >= pivot])
        PivotID = index[0]
        Pivot = getRowObject(PivotID, TableName)
        lesser_index = []
        greater_index = []
        for RowID in index[1:]:
            RowObject = getRowObject(RowID, TableName)
            if compareLESS(RowObject, Pivot, ParameterNames):
                lesser_index += [RowID]
            else:
                greater_index += [RowID]
        # print 'QS: lesser_index = '+str(lesser_index)
        # print 'QS: greater_index = '+str(greater_index)
        lesser = quickSort(lesser_index, TableName, ParameterNames, Accending)
        greater = quickSort(greater_index, TableName, ParameterNames, Accending)
        # return lesser + [pivot_index] + greater
        if Accending:
            return lesser + [PivotID] + greater
        else:
            return greater + [PivotID] + lesser


# Sorting must work well on the table itself!


def sort(
    TableName,
    DestinationTableName=None,
    ParameterNames=None,
    Accending=True,
    Output=False,
    File=None,
):
    """
    INPUT PARAMETERS: 
        TableName:                name of source table          (required)
        DestinationTableName:     name of resulting table       (optional)
        ParameterNames:       list of parameters or expressions to sort by    (optional)
        Accending:       sort in ascending (True) or descending (False) order (optional)
        Output:   enable (True) or suppress (False) text output (optional)
        File:     enable (True) or suppress (False) file output (optional)
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Sort a table by a list of it's parameters or expressions.
        The sorted table is saved in DestinationTableName (if specified).
    ---
    EXAMPLE OF USAGE:
        sort('sampletab',ParameterNames=(p1,('+',p1,p2)))
    ---
    """
    number_of_rows = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]
    index = list(range(0, number_of_rows))
    # print 'num = '+str(number_of_rows)
    if not DestinationTableName:
        DestinationTableName = TableName
    # if names are not provided use all parameters in sorting
    if not ParameterNames:
        ParameterNames = LOCAL_TABLE_CACHE[TableName]["header"]["order"]
    elif type(ParameterNames) not in set([list, tuple]):
        # fix of stupid bug where ('p1',) != ('p1')
        ParameterNames = [ParameterNames]
    # print 'SRT: ParameterNames = '+str(ParameterNames)
    # print 'parnames: '+str(ParameterNames)
    index_sorted = quickSort(index, TableName, ParameterNames, Accending)
    arrangeTable(TableName, DestinationTableName, index_sorted)
    if Output:
        outputTable(DestinationTableName, File=File)


# /SORTING ==========================================================


# GROUPING ==========================================================

# GROUP_INDEX global auxillary structure is a Dictionary,
#   which has the following properties:
#      1) Each key is a composite variable:
#          [array of values of ParameterNames variable
#           STREAM_UPDATE_FLAG]
#      2) Each value is an index in LOCAL_TABLE_CACHE[TableName]['data'][...],
#          corresponding to this key
#   STREAM_UPDATE_FLAG = TRUE if value in GROUP_INDEX needs updating
#                      = FALSE otherwise
#   If no grouping variables are specified (GroupParameterNames==None)
#    than the following key is used: "__GLOBAL__"


def group(
    TableName,
    DestinationTableName=QUERY_BUFFER,
    ParameterNames=None,
    GroupParameterNames=None,
    Output=True,
):
    """
    INPUT PARAMETERS: 
        TableName:                name of source table          (required)
        DestinationTableName:     name of resulting table       (optional)
        ParameterNames:       list of parameters or expressions to take       (optional)
        GroupParameterNames:  list of parameters or expressions to group by   (optional)
        Accending:       sort in ascending (True) or descending (False) order (optional)
        Output:   enable (True) or suppress (False) text output (optional)
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        none
    ---
    EXAMPLE OF USAGE:
        group('sampletab',ParameterNames=('p1',('sum','p2')),GroupParameterNames=('p1'))
        ... makes grouping by p1,p2. For each group it calculates sum of p2 values.
    ---
    """
    # Implements such functions as:
    # count,sum,avg,min,max,ssq etc...
    # 1) ParameterNames can contain group functions
    # 2) GroupParameterNames can't contain group functions
    # 3) If ParameterNames contains parameters defined by LET directive,
    #    it IS visible in the sub-context of GroupParameterNames
    # 4) Parameters defined in GroupParameterNames are NOT visible in ParameterNames
    # 5) ParameterNames variable represents the structure of the resulting table/collection
    # 6) GroupParameterNames can contain either par_names or expressions with par_names
    # Clear old GROUP_INDEX value
    clearGroupIndex()
    # Consistency check
    if TableName == DestinationTableName:
        raise Exception("TableName and DestinationTableName must be different")
    # if not ParameterNames: ParameterNames=LOCAL_TABLE_CACHE[TableName]['header']['order']
    # Prepare the new DestinationTable
    RowObjectDefault = getDefaultRowObject(TableName)
    VarDictionary = getVarDictionary(RowObjectDefault)
    ContextFormat = getContextFormat(RowObjectDefault)
    RowObjectDefaultNew = newRowObject(
        ParameterNames, RowObjectDefault, VarDictionary, ContextFormat
    )
    dropTable(DestinationTableName)  # redundant
    createTable(DestinationTableName, RowObjectDefaultNew)
    # Loop through rows of source Table
    # On each iteration group functions update GROUP_INDEX (see description above)
    number_of_rows = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]
    # STAGE 1: CREATE GROUPS
    print("LOOP:")
    for RowID in range(0, number_of_rows):
        print("--------------------------------")
        print("RowID=" + str(RowID))
        # RowObject from source table
        RowObject = getRowObject(RowID, TableName)
        VarDictionary = getVarDictionary(RowObject)
        print("VarDictionary=" + str(VarDictionary))
        # This is a trick which makes evaluateExpression function
        #   not consider first expression as an operation
        GroupParameterNames_ = ["LIST"] + list(GroupParameterNames)
        GroupIndexKey = evaluateExpression(GroupParameterNames_, VarDictionary)
        # List is an unhashable type in Python!
        GroupIndexKey = tuple(GroupIndexKey)
        initializeGroup(GroupIndexKey)
        print("GROUP_INDEX=" + str(GROUP_INDEX))
        ContextFormat = getContextFormat(RowObject)
        RowObjectNew = newRowObject(
            ParameterNames, RowObject, VarDictionary, ContextFormat, GroupIndexKey
        )
        RowIDGroup = GROUP_INDEX[GroupIndexKey]["ROWID"]
        setRowObject(RowIDGroup, RowObjectNew, DestinationTableName)
    # Output result if required
    if Output and DestinationTableName == QUERY_BUFFER:
        outputTable(DestinationTableName, File=File)


# /GROUPING =========================================================

# EXTRACTING ========================================================


REGEX_INTEGER = r"[+-]?\d+"
REGEX_STRING = r"[^\s]+"
REGEX_FLOAT_F = r"[+-]?\d*\.?\d+"
REGEX_FLOAT_E = r"[+-]?\d*\.?\d+[eEfF]?[+-]?\d+"


def REGEX_INTEGER_FIXCOL(n):
    return r"\d{%d}" % n


def REGEX_STRING_FIXCOL(n):
    return r"[^\s]{%d}" % n


def REGEX_FLOAT_F_FIXCOL(n):
    return r"[\+\-\.\d]{%d}" % n


def REGEX_FLOAT_E_FIXCOL(n):
    return r"[\+\-\.\deEfF]{%d}" % n


# Extract sub-columns from string column


def extractColumns(
    TableName, SourceParameterName, ParameterFormats, ParameterNames=None, FixCol=False
):
    """
    INPUT PARAMETERS: 
        TableName:             name of source table              (required)
        SourceParameterName:   name of source column to process  (required)
        ParameterFormats:      c formats of unpacked parameters  (required)
        ParameterNames:        list of resulting parameter names (optional)
        FixCol:      column-fixed (True) format of source column (optional)
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Note, that this function is aimed to do some extra job on
        interpreting string parameters which is normally supposed
        to be done by the user.
    ---
    EXAMPLE OF USAGE:
        extractColumns('sampletab',SourceParameterName='p5',
                        ParameterFormats=('%d','%d','%d'),
                        ParameterNames=('p5_1','p5_2','p5_3'))
        This example extracts three integer parameters from
        a source column 'p5' and puts results in ('p5_1','p5_2','p5_3').
    ---
    """
    # ParameterNames = just the names without expressions
    # ParFormats contains python formats for par extraction
    # Example: ParameterNames=('v1','v2','v3')
    #          ParameterFormats=('%1s','%1s','%1s')
    # By default the format of parameters is column-fixed
    if type(
        LOCAL_TABLE_CACHE[TableName]["header"]["default"][SourceParameterName]
    ) not in set([str, six.text_type]):
        raise Exception("Source parameter must be a string")
    i = -1
    # bug when (a,) != (a)
    if ParameterNames and type(ParameterNames) not in set([list, tuple]):
        ParameterNames = [ParameterNames]
    if ParameterFormats and type(ParameterFormats) not in set([list, tuple]):
        ParameterFormats = [ParameterFormats]
    # if ParameterNames is empty, fill it with #1-2-3-...
    if not ParameterNames:
        ParameterNames = []
        # using naming convension #i, i=0,1,2,3...
        for par_format in ParameterFormats:
            while True:
                i += 1
                par_name = "#%d" % i
                fmt = LOCAL_TABLE_CACHE[TableName]["header"]["format"].get(
                    par_name, None
                )
                if not fmt:
                    break
            ParameterNames.append(par_name)
    # check if ParameterNames are valid
    Intersection = set(ParameterNames).intersection(
        LOCAL_TABLE_CACHE[TableName]["header"]["order"]
    )
    if Intersection:
        raise Exception("Parameters %s already exist" % str(list(Intersection)))
    # loop over ParameterNames to prepare LOCAL_TABLE_CACHE
    i = 0
    for par_name in ParameterNames:
        par_format = ParameterFormats[i]
        LOCAL_TABLE_CACHE[TableName]["header"]["format"][par_name] = par_format
        LOCAL_TABLE_CACHE[TableName]["data"][par_name] = []
        i += 1
    # append new parameters in order list
    LOCAL_TABLE_CACHE[TableName]["header"]["order"] += ParameterNames
    # cope with default values
    i = 0
    format_regex = []
    format_types = []
    # print 'ParameterNames='+str(ParameterNames)
    for par_format in ParameterFormats:
        par_name = ParameterNames[i]
        regex = FORMAT_PYTHON_REGEX
        # print 'par_name: '+par_name
        # print 'par_format: '+par_format
        (lng, trail, lngpnt, ty) = re.search(regex, par_format).groups()
        ty = ty.lower()
        if ty == "d":
            par_type = int
            if FixCol:
                format_regex_part = REGEX_INTEGER_FIXCOL(lng)
            else:
                format_regex_part = REGEX_INTEGER
        elif ty == "s":
            par_type = str
            if FixCol:
                format_regex_part = REGEX_STRING_FIXCOL(lng)
            else:
                format_regex_part = REGEX_STRING
        elif ty == "f":
            par_type = float
            if FixCol:
                format_regex_part = REGEX_FLOAT_F_FIXCOL(lng)
            else:
                format_regex_part = REGEX_FLOAT_F
        elif ty == "e":
            par_type = float
            if FixCol:
                format_regex_part = REGEX_FLOAT_E_FIXCOL(lng)
            else:
                format_regex_part = REGEX_FLOAT_E
        else:
            raise Exception("Unknown data type")
        format_regex.append("(" + format_regex_part + ")")
        format_types.append(par_type)
        def_val = getDefaultValue(par_type)
        LOCAL_TABLE_CACHE[TableName]["header"]["default"][par_name] = def_val
        i += 1
    format_regex = r"\s*".join(format_regex)
    # print 'format_regex='+str(format_regex)
    # return format_regex
    # loop through values of SourceParameter
    for SourceParameterString in LOCAL_TABLE_CACHE[TableName]["data"][
        SourceParameterName
    ]:
        try:
            ExtractedValues = list(
                re.search(format_regex, SourceParameterString).groups()
            )
        except:
            raise Exception('Error with line "%s"' % SourceParameterString)
        i = 0
        # loop through all parameters which are supposed to be extracted
        for par_name in ParameterNames:
            # print 'ExtractedValues[i]='+ExtractedValues[i]
            # print 'par_name='+par_name
            par_value = format_types[i](ExtractedValues[i])
            LOCAL_TABLE_CACHE[TableName]["data"][par_name].append(par_value)
            i += 1
    # explicitly check that number of rows are equal
    number_of_rows = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]
    number_of_rows2 = len(LOCAL_TABLE_CACHE[TableName]["data"][SourceParameterName])
    number_of_rows3 = len(LOCAL_TABLE_CACHE[TableName]["data"][ParameterNames[0]])
    if not (number_of_rows == number_of_rows2 == number_of_rows3):
        raise Exception("Error while extracting parameters: check your regexp")


# Split string columns into sub-columns with given names


def splitColumn(TableName, SourceParameterName, ParameterNames, Splitter):
    pass


# /EXTRACTING =======================================================

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# /LOCAL DATABASE MANAGEMENT SYSTEM
# ---------------------------------------------------------------
# ---------------------------------------------------------------


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# GLOBAL API FUNCTIONS
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


def mergeParlist(*arg):
    # Merge parlists and remove duplicates.
    # Argument contains a list of lists/tuples.
    container = []
    for a in arg:
        container += list(a)
    result = []
    index = set()
    for par_name in container:
        if par_name not in index:
            index.add(par_name)
            result.append(par_name)
    return result


# Define parameter groups to simplify the usage of fetch_
PARLIST_DOTPAR = [
    "par_line",
]
PARLIST_ID = [
    "trans_id",
]
PARLIST_STANDARD = [
    "molec_id",
    "local_iso_id",
    "nu",
    "sw",
    "a",
    "elower",
    "gamma_air",
    "delta_air",
    "gamma_self",
    "n_air",
    "n_self",
    "gp",
    "gpp",
]
PARLIST_LABELS = ["statep", "statepp"]
PARLIST_LINEMIXING = ["y_air", "y_self"]

PARLIST_VOIGT_AIR = ["gamma_air", "delta_air", "deltap_air", "n_air"]
PARLIST_VOIGT_SELF = ["gamma_self", "delta_self", "deltap_self", "n_self"]
PARLIST_VOIGT_H2 = ["gamma_H2", "delta_H2", "deltap_H2", "n_H2"]
PARLIST_VOIGT_CO2 = ["gamma_CO2", "delta_CO2", "n_CO2"]
PARLIST_VOIGT_HE = ["gamma_He", "delta_He", "n_He"]
PARLIST_VOIGT_ALL = mergeParlist(
    PARLIST_VOIGT_AIR,
    PARLIST_VOIGT_SELF,
    PARLIST_VOIGT_H2,
    PARLIST_VOIGT_CO2,
    PARLIST_VOIGT_HE,
)

PARLIST_SDVOIGT_AIR = ["gamma_air", "delta_air", "deltap_air", "n_air", "SD_air"]
PARLIST_SDVOIGT_SELF = ["gamma_self", "delta_self", "deltap_self", "n_self", "SD_self"]
PARLIST_SDVOIGT_H2 = []
PARLIST_SDVOIGT_CO2 = []
PARLIST_SDVOIGT_HE = []
PARLIST_SDVOIGT_ALL = mergeParlist(
    PARLIST_SDVOIGT_AIR,
    PARLIST_SDVOIGT_SELF,
    PARLIST_SDVOIGT_H2,
    PARLIST_SDVOIGT_CO2,
    PARLIST_SDVOIGT_HE,
)

PARLIST_GALATRY_AIR = ["gamma_air", "delta_air", "deltap_air", "n_air", "beta_g_air"]
PARLIST_GALATRY_SELF = [
    "gamma_self",
    "delta_self",
    "deltap_self",
    "n_self",
    "beta_g_self",
]
PARLIST_GALATRY_H2 = []
PARLIST_GALATRY_CO2 = []
PARLIST_GALATRY_HE = []
PARLIST_GALATRY_ALL = mergeParlist(
    PARLIST_GALATRY_AIR,
    PARLIST_GALATRY_SELF,
    PARLIST_GALATRY_H2,
    PARLIST_GALATRY_CO2,
    PARLIST_GALATRY_HE,
)

PARLIST_ALL = mergeParlist(
    PARLIST_ID,
    PARLIST_DOTPAR,
    PARLIST_STANDARD,
    PARLIST_LABELS,
    PARLIST_LINEMIXING,
    PARLIST_VOIGT_ALL,
    PARLIST_SDVOIGT_ALL,
    PARLIST_GALATRY_ALL,
)

PARAMETER_GROUPS = {
    "par_line": PARLIST_DOTPAR,
    "160-char": PARLIST_DOTPAR,
    ".par": PARLIST_DOTPAR,
    "id": PARLIST_ID,
    "standard": PARLIST_STANDARD,
    "labels": PARLIST_LABELS,
    "linemixing": PARLIST_LINEMIXING,
    "voigt_air": PARLIST_VOIGT_AIR,
    "voigt_self": PARLIST_VOIGT_SELF,
    "voigt_h2": PARLIST_VOIGT_H2,
    "voigt_co2": PARLIST_VOIGT_CO2,
    "voigt_he": PARLIST_VOIGT_HE,
    "voigt": PARLIST_VOIGT_ALL,
    "sdvoigt_air": PARLIST_SDVOIGT_AIR,
    "sdvoigt_self": PARLIST_SDVOIGT_SELF,
    "sdvoigt_h2": PARLIST_SDVOIGT_H2,
    "sdvoigt_co2": PARLIST_SDVOIGT_CO2,
    "sdvoigt_he": PARLIST_SDVOIGT_HE,
    "sdvoigt": PARLIST_SDVOIGT_ALL,
    "galatry_air": PARLIST_GALATRY_AIR,
    "galatry_self": PARLIST_GALATRY_SELF,
    "galatry_h2": PARLIST_GALATRY_H2,
    "galatry_co2": PARLIST_GALATRY_CO2,
    "galatry_he": PARLIST_GALATRY_HE,
    "galatry": PARLIST_GALATRY_ALL,
    "all": PARLIST_ALL,
}


def prepareParlist(pargroups=[], params=[], dotpar=True):
    # Apply defaults
    parlist_default = []
    if dotpar:
        parlist_default += ["par_line"]
    # parlist_default += PARAMETER_GROUPS['id']

    # Make a dictionary of "assumed" parameters.
    ASSUMED_PARAMS = {}
    if "par_line" in set(parlist_default):
        ASSUMED_PARAMS = HITRAN_DEFAULT_HEADER["format"]

    parlist = parlist_default

    # Iterate over parameter groups.
    for pargroup in pargroups:
        pargroup = pargroup.lower()
        parlist += PARAMETER_GROUPS[pargroup]

    # Iterate over single parameters.
    for param in params:
        param = param.lower()
        parlist.append(param)

    # Clean up parameter list.
    parlist = mergeParlist(parlist)
    result = []
    for param in parlist:
        if param not in ASSUMED_PARAMS:
            result.append(param)

    return result


def prepareHeader(parlist):
    HEADER = {
        "table_name": "",
        "number_of_rows": -1,
        "format": {},
        "default": {},
        "table_type": "column-fixed",
        "size_in_bytes": -1,
        "order": [],
        "description": {},
    }

    # Add column-fixed 160-character part, if specified in parlist.
    if "par_line" in set(parlist):
        HEADER["order"] = HITRAN_DEFAULT_HEADER["order"]
        HEADER["format"] = HITRAN_DEFAULT_HEADER["format"]
        HEADER["default"] = HITRAN_DEFAULT_HEADER["default"]
        HEADER["description"] = HITRAN_DEFAULT_HEADER["description"]

    # Insert all other parameters in the "extra" section of the header.
    # while 'par_line' in parlist: parlist.remove('par_line')
    plist = [v for v in parlist if v != "par_line"]
    HEADER["extra"] = []
    HEADER["extra_format"] = {}
    HEADER["extra_separator"] = ","
    for param in plist:
        param = param.lower()
        HEADER["extra"].append(param)
        HEADER["extra_format"][param] = PARAMETER_META[param]["default_fmt"]

    return HEADER


def queryHITRAN(
    TableName,
    iso_id_list,
    numin,
    numax,
    pargroups=[],
    params=[],
    dotpar=True,
    head=False,
):
    # import httplib
    # conn = httplib.HTTPConnection('hitranazure.cloudapp.com')
    # conn.Request('')
    # r = conn.getresponse()
    # print r.status, r.reason
    # data1 = data1.read
    # TableHeader = HITRAN_DEFAULT_HEADER # deprecated
    ParameterList = prepareParlist(pargroups=pargroups, params=params, dotpar=dotpar)
    TableHeader = prepareHeader(ParameterList)
    TableHeader["table_name"] = TableName
    DataFileName = VARIABLES["BACKEND_DATABASE_NAME"] + "/" + TableName + ".data"
    HeaderFileName = VARIABLES["BACKEND_DATABASE_NAME"] + "/" + TableName + ".header"
    # if TableName in LOCAL_TABLE_CACHE.keys():
    #   raise Exception('Table \"%s\" exists' % TableName)
    # if os.path.isfile(DataFileName):
    #   raise Exception('File \"%s\" exists' % DataFileName)
    # if os.path.isfile(HeaderFileName):
    #   raise Exception('!!File \"%s\" exists' % HeaderFileName)
    # create URL
    iso_id_list_str = [str(iso_id) for iso_id in iso_id_list]
    iso_id_list_str = ",".join(iso_id_list_str)
    # url = 'http://hitran.cloudapp.net' + '/lbl/5?' + \
    # url = 'http://hitranazure.cloudapp.net' + '/lbl/5?' + \
    # 'iso_ids_list=' + iso_id_list_str + '&' + \
    # 'numin=' + str(numin) + '&' + \
    # 'numax=' + str(numax) + '&' + \
    # 'access=api' + '&' + \
    #'key=' + GLOBAL_HITRAN_APIKEY
    if pargroups or params:  # custom par search
        url = (
            GLOBAL_HOST
            + "/lbl/api?"
            + "iso_ids_list="
            + iso_id_list_str
            + "&"
            + "numin="
            + str(numin)
            + "&"
            + "numax="
            + str(numax)
            + "&"
            + "head="
            + str(head)
            + "&"
            + "fixwidth=0&sep=[comma]&"
            + "request_params="
            + ",".join(ParameterList)
        )
    else:  # old-fashioned .par search
        url = (
            GLOBAL_HOST
            + "/lbl/api?"
            + "iso_ids_list="
            + iso_id_list_str
            + "&"
            + "numin="
            + str(numin)
            + "&"
            + "numax="
            + str(numax)
        )
    # raise Exception(url)
    # Download data by chunks.
    try:
        req = six.moves.urllib.request.urlopen(url)
    except six.moves.urllib.error.HTTPError:
        raise Exception("Failed to retrieve data for given parameters.")
    except six.moves.urllib.error.URLError:
        raise Exception(
            "Cannot connect to %s. Try again or edit GLOBAL_HOST variable."
            % GLOBAL_HOST
        )
    # CHUNK = 16 * 1024 # default value
    CHUNK = 64 * 1024
    print("BEGIN DOWNLOAD: " + TableName)
    with open(DataFileName, "w") as fp:
        while True:
            chunk = req.read(CHUNK)
            if not chunk:
                break
            fp.write(chunk.decode("utf-8"))
            print("  %d bytes written to %s" % (CHUNK, DataFileName))
    with open(HeaderFileName, "w") as fp:
        fp.write(json.dumps(TableHeader, indent=2))
        print("Header written to %s" % HeaderFileName)
    print("END DOWNLOAD")
    # Set comment
    # Get this table to LOCAL_TABLE_CACHE
    storage2cache(TableName)
    print("PROCESSED")


# NODE CODE
NODE_READY = False

# Node initialization


def nodeInit():
    # very unoptimal, since it loads all tables in memory!!
    # loadCache()
    databaseBegin()  # DB backend level, start transaction
    NODE_READY = True


# returns a table instance created from Query object


def globalSelectInto(NewTablePath, SourceTablePath, ParameterNames, Conditions):
    # creates table from parsed data
    # and store it in the database DB

    dbname, tablename, nodename = NewTablePath.split("::")
    dbname1, tablename1, nodename1 = SourceTablePath.split("::")

    if not NODE_READY:
        raise Exception('Node "%s" is not ready. Call nodeInit()' % NODE_NAME)

    # should get rid of selectLocal as planning to use network interface
    # ......    selectLocal OR selectRemote

    pass


# ---------------------------------------------------------------

# query_string - query written in the
# formal language of local database frontend


def makeQuery(query_string, Connection=GLOBAL_CONNECTION):
    # makes a query to remote server
    # using connection instance
    pass


# ---------- DATABASE FRONTEND END -------------

# simple implementation of getting a line list from a remote server


def getLinelist(local_name, query, api_key):
    return makeQuery(local_name)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# / GLOBABL API FUNCTIONS
# -------------------------------------------------------------------
# -------------------------------------------------------------------


# ---------------- FILTER ---------------------------------------------


def filter(TableName, Conditions):
    select(TableName=TableName, Conditions=Conditions, Output=False)


# ---------------------- ISO.PY ---------------------------------------


ISO_ID_INDEX = {"M": 0, "I": 1, "iso_name": 2, "abundance": 3, "mass": 4, "mol_name": 5}

#    id           M    I    iso_name                    abundance           mass        mol_name

ISO_ID = {
    1: [1, 1, "H2(16O)", 0.997317, 18.010565, "H2O"],
    2: [1, 2, "H2(18O)", 0.00199983, 20.014811, "H2O"],
    3: [1, 3, "H2(17O)", 0.000372, 19.01478, "H2O"],
    4: [1, 4, "HD(16O)", 0.00031069, 19.01674, "H2O"],
    5: [1, 5, "HD(18O)", 0.000000623, 21.020985, "H2O"],
    6: [1, 6, "HD(17O)", 0.000000116, 20.020956, "H2O"],
    7: [2, 1, "(12C)(16O)2", 0.9842, 43.98983, "CO2"],
    8: [2, 2, "(13C)(16O)2", 0.01106, 44.993185, "CO2"],
    9: [2, 3, "(16O)(12C)(18O)", 0.0039471, 45.994076, "CO2"],
    10: [2, 4, "(16O)(12C)(17O)", 0.000734, 44.994045, "CO2"],
    11: [2, 5, "(16O)(13C)(18O)", 0.00004434, 46.997431, "CO2"],
    12: [2, 6, "(16O)(13C)(17O)", 0.00000825, 45.9974, "CO2"],
    13: [2, 7, "(12C)(18O)2", 0.0000039573, 47.998322, "CO2"],
    14: [2, 8, "(17O)(12C)(18O)", 0.00000147, 46.998291, "CO2"],
    15: [2, 0, "(13C)(18O)2", 0.000000044967, 49.001675, "CO2"],
    120: [2, 11, "(18O)(13C)(17O)", 0.00000001654, 48.00165, "CO2"],
    121: [2, 9, "(12C)(17O)2", 0.0000001368, 45.998262, "CO2"],
    16: [3, 1, "(16O)3", 0.992901, 47.984745, "O3"],
    17: [3, 2, "(16O)(16O)(18O)", 0.00398194, 49.988991, "O3"],
    18: [3, 3, "(16O)(18O)(16O)", 0.00199097, 49.988991, "O3"],
    19: [3, 4, "(16O)(16O)(17O)", 0.00074, 48.98896, "O3"],
    20: [3, 5, "(16O)(17O)(16O)", 0.00037, 48.98896, "O3"],
    21: [4, 1, "(14N)2(16O)", 0.990333, 44.001062, "N2O"],
    22: [4, 2, "(14N)(15N)(16O)", 0.0036409, 44.998096, "N2O"],
    23: [4, 3, "(15N)(14N)(16O)", 0.0036409, 44.998096, "N2O"],
    24: [4, 4, "(14N)2(18O)", 0.00198582, 46.005308, "N2O"],
    25: [4, 5, "(14N)2(17O)", 0.000369, 45.005278, "N2O"],
    26: [5, 1, "(12C)(16O)", 0.98654, 27.994915, "CO"],
    27: [5, 2, "(13C)(16O)", 0.01108, 28.99827, "CO"],
    28: [5, 3, "(12C)(18O)", 0.0019782, 29.999161, "CO"],
    29: [5, 4, "(12C)(17O)", 0.000368, 28.99913, "CO"],
    30: [5, 5, "(13C)(18O)", 0.00002222, 31.002516, "CO"],
    31: [5, 6, "(13C)(17O)", 0.00000413, 30.002485, "CO"],
    32: [6, 1, "(12C)H4", 0.98827, 16.0313, "CH4"],
    33: [6, 2, "(13C)H4", 0.0111, 17.034655, "CH4"],
    34: [6, 3, "(12C)H3D", 0.00061575, 17.037475, "CH4"],
    35: [6, 4, "(13C)H3D", 0.0000049203, 18.04083, "CH4"],
    36: [7, 1, "(16O)2", 0.995262, 31.98983, "O2"],
    37: [7, 2, "(16O)(18O)", 0.00399141, 33.994076, "O2"],
    38: [7, 3, "(16O)(17O)", 0.000742, 32.994045, "O2"],
    39: [8, 1, "(14N)(16O)", 0.993974, 29.997989, "NO"],
    40: [8, 2, "(15N)(16O)", 0.0036543, 30.995023, "NO"],
    41: [8, 3, "(14N)(18O)", 0.00199312, 32.002234, "NO"],
    42: [9, 1, "(32S)(16O)2", 0.94568, 63.961901, "SO2"],
    43: [9, 2, "(34S)(16O)2", 0.04195, 65.957695, "SO2"],
    44: [10, 1, "(14N)(16O)2", 0.991616, 45.992904, "NO2"],
    45: [11, 1, "(14N)H3", 0.9958715, 17.026549, "NH3"],
    46: [11, 2, "(15N)H3", 0.0036613, 18.023583, "NH3"],
    47: [12, 1, "H(14N)(16O)3", 0.98911, 62.995644, "HNO3"],
    117: [12, 2, "H(15N)(16O)3", 0.003636, 63.99268, "HNO3"],
    48: [13, 1, "(16O)H", 0.997473, 17.00274, "OH"],
    49: [13, 2, "(18O)H", 0.00200014, 19.006986, "OH"],
    50: [13, 3, "(16O)D", 0.00015537, 18.008915, "OH"],
    51: [14, 1, "H(19F)", 0.99984425, 20.006229, "HF"],
    110: [14, 2, "D(19F)", 0.000115, 21.0125049978, "HF"],
    52: [15, 1, "H(35Cl)", 0.757587, 35.976678, "HCl"],
    53: [15, 2, "H(37Cl)", 0.242257, 37.973729, "HCl"],
    107: [15, 3, "D(35Cl)", 0.000118005, 36.9829544578, "HCl"],
    108: [15, 4, "D(37Cl)", 0.000037735, 38.9800043678, "HCl"],
    54: [16, 1, "H(79Br)", 0.50678, 79.92616, "HBr"],
    55: [16, 2, "H(81Br)", 0.49306, 81.924115, "HBr"],
    111: [16, 3, "D(79Br)", 0.0000582935, 80.9324388778, "HBr"],
    112: [16, 4, "D(81Br)", 0.0000567065, 82.9303923778, "HBr"],
    56: [17, 1, "H(127I)", 0.99984425, 127.912297, "HI"],
    113: [17, 2, "D(127I)", 0.000115, 128.918574778, "HI"],
    57: [18, 1, "(35Cl)(16O)", 0.75591, 50.963768, "ClO"],
    58: [18, 2, "(37Cl)(16O)", 0.24172, 52.960819, "ClO"],
    59: [19, 1, "(16O)(12C)(32S)", 0.93739, 59.966986, "OCS"],
    60: [19, 2, "(16O)(12C)(34S)", 0.04158, 61.96278, "OCS"],
    61: [19, 3, "(16O)(13C)(32S)", 0.01053, 60.970341, "OCS"],
    62: [19, 4, "(16O)(12C)(33S)", 0.01053, 60.966371, "OCS"],
    63: [19, 5, "(18O)(12C)(32S)", 0.00188, 61.971231, "OCS"],
    64: [20, 1, "H2(12C)(16O)", 0.98624, 30.010565, "H2CO"],
    65: [20, 2, "H2(13C)(16O)", 0.01108, 31.01392, "H2CO"],
    66: [20, 3, "H2(12C)(18O)", 0.0019776, 32.014811, "H2CO"],
    67: [21, 1, "H(16O)(35Cl)", 0.75579, 51.971593, "HOCl"],
    68: [21, 2, "H(16O)(37Cl)", 0.24168, 53.968644, "HOCl"],
    69: [22, 1, "(14N)2", 0.9926874, 28.006147, "N2"],
    118: [22, 2, "(14N)(15N)", 0.0072535, 29.997989, "N2"],
    70: [23, 1, "H(12C)(14N)", 0.98511, 27.010899, "HCN"],
    71: [23, 2, "H(13C)(14N)", 0.01107, 28.014254, "HCN"],
    72: [23, 3, "H(12C)(15N)", 0.0036217, 28.007933, "HCN"],
    73: [24, 1, "(12C)H3(35Cl)", 0.74894, 49.992328, "CH3Cl"],
    74: [24, 2, "(12C)H3(37Cl)", 0.23949, 51.989379, "CH3Cl"],
    75: [25, 1, "H2(16O)2", 0.994952, 34.00548, "H2O2"],
    76: [26, 1, "(12C)2H2", 0.9776, 26.01565, "C2H2"],
    77: [26, 2, "(12C)(13C)H2", 0.02197, 27.019005, "C2H2"],
    105: [26, 3, "(12C)2HD", 0.00030455, 27.021825, "C2H2"],
    78: [27, 1, "(12C)2H6", 0.97699, 30.04695, "C2H6"],
    106: [27, 2, "(12C)H3(13C)H3", 0.021952611, 31.050305, "C2H6"],
    79: [28, 1, "(31P)H3", 0.99953283, 33.997238, "PH3"],
    80: [29, 1, "(12C)(16O)(19F)2", 0.98654, 65.991722, "COF2"],
    119: [29, 2, "(13C)(16O)(19F)2", 0.0110834, 66.995083, "COF2"],
    81: [31, 1, "H2(32S)", 0.94988, 33.987721, "H2S"],
    82: [31, 2, "H2(34S)", 0.04214, 35.983515, "H2S"],
    83: [31, 3, "H2(33S)", 0.007498, 34.987105, "H2S"],
    84: [32, 1, "H(12C)(16O)(16O)H", 0.983898, 46.00548, "HCOOH"],
    85: [33, 1, "H(16O)2", 0.995107, 32.997655, "HO2"],
    86: [34, 1, "(16O)", 0.997628, 15.994915, "O"],
    87: [36, 1, "(14N)(16O)+", 0.993974, 29.997989, "NOp"],
    88: [37, 1, "H(16O)(79Br)", 0.5056, 95.921076, "HOBr"],
    89: [37, 2, "H(16O)(81Br)", 0.4919, 97.919027, "HOBr"],
    90: [38, 1, "(12C)2H4", 0.9773, 28.0313, "C2H4"],
    91: [38, 2, "(12C)H2(13C)H2", 0.02196, 29.034655, "C2H4"],
    92: [39, 1, "(12C)H3(16O)H", 0.98593, 32.026215, "CH3OH"],
    93: [40, 1, "(12C)H3(79Br)", 0.5013, 93.941811, "CH3Br"],
    94: [40, 2, "(12C)H3(81Br)", 0.48766, 95.939764, "CH3Br"],
    95: [41, 1, "(12C)H3(12C)(14N)", 0.97482, 41.026549, "CH3CN"],
    96: [42, 1, "(12C)(19F)4", 0.9893, 87.993616, "CF4"],
    116: [43, 1, "(12C)4H2", 0.955998, 50.01565, "C4H2"],
    109: [44, 1, "H(12C)3(14N)", 0.9646069, 51.01089903687, "HC3N"],
    103: [45, 1, "H2", 0.999688, 2.01565, "H2"],
    115: [45, 2, "HD", 0.00022997, 3.021825, "H2"],
    97: [46, 1, "(12C)(32S)", 0.939624, 43.971036, "CS"],
    98: [46, 2, "(12C)(34S)", 0.0416817, 45.966787, "CS"],
    99: [46, 3, "(13C)(32S)", 0.0105565, 44.974368, "CS"],
    100: [46, 4, "(12C)(33S)", 0.00741668, 44.970399, "CS"],
    114: [47, 1, "(32S)(16O)3", 0.9423964, 79.95682, "SO3"],
    101: [1001, 1, "H", None, None, "H"],
    102: [1002, 1, "He", None, None, "He"],
    104: [1018, 1, "Ar", None, None, "Ar"],
}

ISO_INDEX = {"id": 0, "iso_name": 1, "abundance": 2, "mass": 3, "mol_name": 4}

#        M    I             id    iso_name                    abundance           mass        mol_name

ISO = {
    (1, 1): [1, "H2(16O)", 0.997317, 18.010565, "H2O"],
    (1, 2): [2, "H2(18O)", 0.00199983, 20.014811, "H2O"],
    (1, 3): [3, "H2(17O)", 0.000372, 19.01478, "H2O"],
    (1, 4): [4, "HD(16O)", 0.00031069, 19.01674, "H2O"],
    (1, 5): [5, "HD(18O)", 0.000000623, 21.020985, "H2O"],
    (1, 6): [6, "HD(17O)", 0.000000116, 20.020956, "H2O"],
    (2, 1): [7, "(12C)(16O)2", 0.9842, 43.98983, "CO2"],
    (2, 2): [8, "(13C)(16O)2", 0.01106, 44.993185, "CO2"],
    (2, 3): [9, "(16O)(12C)(18O)", 0.0039471, 45.994076, "CO2"],
    (2, 4): [10, "(16O)(12C)(17O)", 0.000734, 44.994045, "CO2"],
    (2, 5): [11, "(16O)(13C)(18O)", 0.00004434, 46.997431, "CO2"],
    (2, 6): [12, "(16O)(13C)(17O)", 0.00000825, 45.9974, "CO2"],
    (2, 7): [13, "(12C)(18O)2", 0.0000039573, 47.998322, "CO2"],
    (2, 8): [14, "(17O)(12C)(18O)", 0.00000147, 46.998291, "CO2"],
    (2, 0): [15, "(13C)(18O)2", 0.000000044967, 49.001675, "CO2"],
    (2, 11): [120, "(18O)(13C)(17O)", 0.00000001654, 48.00165, "CO2"],
    (2, 9): [121, "(12C)(17O)2", 0.0000001368, 45.998262, "CO2"],
    (3, 1): [16, "(16O)3", 0.992901, 47.984745, "O3"],
    (3, 2): [17, "(16O)(16O)(18O)", 0.00398194, 49.988991, "O3"],
    (3, 3): [18, "(16O)(18O)(16O)", 0.00199097, 49.988991, "O3"],
    (3, 4): [19, "(16O)(16O)(17O)", 0.00074, 48.98896, "O3"],
    (3, 5): [20, "(16O)(17O)(16O)", 0.00037, 48.98896, "O3"],
    (4, 1): [21, "(14N)2(16O)", 0.990333, 44.001062, "N2O"],
    (4, 2): [22, "(14N)(15N)(16O)", 0.0036409, 44.998096, "N2O"],
    (4, 3): [23, "(15N)(14N)(16O)", 0.0036409, 44.998096, "N2O"],
    (4, 4): [24, "(14N)2(18O)", 0.00198582, 46.005308, "N2O"],
    (4, 5): [25, "(14N)2(17O)", 0.000369, 45.005278, "N2O"],
    (5, 1): [26, "(12C)(16O)", 0.98654, 27.994915, "CO"],
    (5, 2): [27, "(13C)(16O)", 0.01108, 28.99827, "CO"],
    (5, 3): [28, "(12C)(18O)", 0.0019782, 29.999161, "CO"],
    (5, 4): [29, "(12C)(17O)", 0.000368, 28.99913, "CO"],
    (5, 5): [30, "(13C)(18O)", 0.00002222, 31.002516, "CO"],
    (5, 6): [31, "(13C)(17O)", 0.00000413, 30.002485, "CO"],
    (6, 1): [32, "(12C)H4", 0.98827, 16.0313, "CH4"],
    (6, 2): [33, "(13C)H4", 0.0111, 17.034655, "CH4"],
    (6, 3): [34, "(12C)H3D", 0.00061575, 17.037475, "CH4"],
    (6, 4): [35, "(13C)H3D", 0.0000049203, 18.04083, "CH4"],
    (7, 1): [36, "(16O)2", 0.995262, 31.98983, "O2"],
    (7, 2): [37, "(16O)(18O)", 0.00399141, 33.994076, "O2"],
    (7, 3): [38, "(16O)(17O)", 0.000742, 32.994045, "O2"],
    (8, 1): [39, "(14N)(16O)", 0.993974, 29.997989, "NO"],
    (8, 2): [40, "(15N)(16O)", 0.0036543, 30.995023, "NO"],
    (8, 3): [41, "(14N)(18O)", 0.00199312, 32.002234, "NO"],
    (9, 1): [42, "(32S)(16O)2", 0.94568, 63.961901, "SO2"],
    (9, 2): [43, "(34S)(16O)2", 0.04195, 65.957695, "SO2"],
    (10, 1): [44, "(14N)(16O)2", 0.991616, 45.992904, "NO2"],
    (11, 1): [45, "(14N)H3", 0.9958715, 17.026549, "NH3"],
    (11, 2): [46, "(15N)H3", 0.0036613, 18.023583, "NH3"],
    (12, 1): [47, "H(14N)(16O)3", 0.98911, 62.995644, "HNO3"],
    (12, 2): [117, "H(15N)(16O)3", 0.003636, 63.99268, "HNO3"],
    (13, 1): [48, "(16O)H", 0.997473, 17.00274, "OH"],
    (13, 2): [49, "(18O)H", 0.00200014, 19.006986, "OH"],
    (13, 3): [50, "(16O)D", 0.00015537, 18.008915, "OH"],
    (14, 1): [51, "H(19F)", 0.99984425, 20.006229, "HF"],
    (14, 2): [110, "D(19F)", 0.000115, 21.0125049978, "HF"],
    (15, 1): [52, "H(35Cl)", 0.757587, 35.976678, "HCl"],
    (15, 2): [53, "H(37Cl)", 0.242257, 37.973729, "HCl"],
    (15, 3): [107, "D(35Cl)", 0.000118005, 36.9829544578, "HCl"],
    (15, 4): [108, "D(37Cl)", 0.000037735, 38.9800043678, "HCl"],
    (16, 1): [54, "H(79Br)", 0.50678, 79.92616, "HBr"],
    (16, 2): [55, "H(81Br)", 0.49306, 81.924115, "HBr"],
    (16, 3): [111, "D(79Br)", 0.0000582935, 80.9324388778, "HBr"],
    (16, 4): [112, "D(81Br)", 0.0000567065, 82.9303923778, "HBr"],
    (17, 1): [56, "H(127I)", 0.99984425, 127.912297, "HI"],
    (17, 2): [113, "D(127I)", 0.000115, 128.918574778, "HI"],
    (18, 1): [57, "(35Cl)(16O)", 0.75591, 50.963768, "ClO"],
    (18, 2): [58, "(37Cl)(16O)", 0.24172, 52.960819, "ClO"],
    (19, 1): [59, "(16O)(12C)(32S)", 0.93739, 59.966986, "OCS"],
    (19, 2): [60, "(16O)(12C)(34S)", 0.04158, 61.96278, "OCS"],
    (19, 3): [61, "(16O)(13C)(32S)", 0.01053, 60.970341, "OCS"],
    (19, 4): [62, "(16O)(12C)(33S)", 0.01053, 60.966371, "OCS"],
    (19, 5): [63, "(18O)(12C)(32S)", 0.00188, 61.971231, "OCS"],
    (20, 1): [64, "H2(12C)(16O)", 0.98624, 30.010565, "H2CO"],
    (20, 2): [65, "H2(13C)(16O)", 0.01108, 31.01392, "H2CO"],
    (20, 3): [66, "H2(12C)(18O)", 0.0019776, 32.014811, "H2CO"],
    (21, 1): [67, "H(16O)(35Cl)", 0.75579, 51.971593, "HOCl"],
    (21, 2): [68, "H(16O)(37Cl)", 0.24168, 53.968644, "HOCl"],
    (22, 1): [69, "(14N)2", 0.9926874, 28.006147, "N2"],
    (22, 2): [118, "(14N)(15N)", 0.0072535, 29.997989, "N2"],
    (23, 1): [70, "H(12C)(14N)", 0.98511, 27.010899, "HCN"],
    (23, 2): [71, "H(13C)(14N)", 0.01107, 28.014254, "HCN"],
    (23, 3): [72, "H(12C)(15N)", 0.0036217, 28.007933, "HCN"],
    (24, 1): [73, "(12C)H3(35Cl)", 0.74894, 49.992328, "CH3Cl"],
    (24, 2): [74, "(12C)H3(37Cl)", 0.23949, 51.989379, "CH3Cl"],
    (25, 1): [75, "H2(16O)2", 0.994952, 34.00548, "H2O2"],
    (26, 1): [76, "(12C)2H2", 0.9776, 26.01565, "C2H2"],
    (26, 2): [77, "(12C)(13C)H2", 0.02197, 27.019005, "C2H2"],
    (26, 3): [105, "(12C)2HD", 0.00030455, 27.021825, "C2H2"],
    (27, 1): [78, "(12C)2H6", 0.97699, 30.04695, "C2H6"],
    (27, 2): [106, "(12C)H3(13C)H3", 0.021952611, 31.050305, "C2H6"],
    (28, 1): [79, "(31P)H3", 0.99953283, 33.997238, "PH3"],
    (29, 1): [80, "(12C)(16O)(19F)2", 0.98654, 65.991722, "COF2"],
    (29, 2): [119, "(13C)(16O)(19F)2", 0.0110834, 66.995083, "COF2"],
    (31, 1): [81, "H2(32S)", 0.94988, 33.987721, "H2S"],
    (31, 2): [82, "H2(34S)", 0.04214, 35.983515, "H2S"],
    (31, 3): [83, "H2(33S)", 0.007498, 34.987105, "H2S"],
    (32, 1): [84, "H(12C)(16O)(16O)H", 0.983898, 46.00548, "HCOOH"],
    (33, 1): [85, "H(16O)2", 0.995107, 32.997655, "HO2"],
    (34, 1): [86, "(16O)", 0.997628, 15.994915, "O"],
    (36, 1): [87, "(14N)(16O)+", 0.993974, 29.997989, "NOp"],
    (37, 1): [88, "H(16O)(79Br)", 0.5056, 95.921076, "HOBr"],
    (37, 2): [89, "H(16O)(81Br)", 0.4919, 97.919027, "HOBr"],
    (38, 1): [90, "(12C)2H4", 0.9773, 28.0313, "C2H4"],
    (38, 2): [91, "(12C)H2(13C)H2", 0.02196, 29.034655, "C2H4"],
    (39, 1): [92, "(12C)H3(16O)H", 0.98593, 32.026215, "CH3OH"],
    (40, 1): [93, "(12C)H3(79Br)", 0.5013, 93.941811, "CH3Br"],
    (40, 2): [94, "(12C)H3(81Br)", 0.48766, 95.939764, "CH3Br"],
    (41, 1): [95, "(12C)H3(12C)(14N)", 0.97482, 41.026549, "CH3CN"],
    (42, 1): [96, "(12C)(19F)4", 0.9893, 87.993616, "CF4"],
    (43, 1): [116, "(12C)4H2", 0.955998, 50.01565, "C4H2"],
    (44, 1): [109, "H(12C)3(14N)", 0.9646069, 51.01089903687, "HC3N"],
    (45, 1): [103, "H2", 0.999688, 2.01565, "H2"],
    (45, 2): [115, "HD", 0.00022997, 3.021825, "H2"],
    (46, 1): [97, "(12C)(32S)", 0.939624, 43.971036, "CS"],
    (46, 2): [98, "(12C)(34S)", 0.0416817, 45.966787, "CS"],
    (46, 3): [99, "(13C)(32S)", 0.0105565, 44.974368, "CS"],
    (46, 4): [100, "(12C)(33S)", 0.00741668, 44.970399, "CS"],
    (47, 1): [114, "(32S)(16O)3", 0.9423964, 79.95682, "SO3"],
    (1001, 1): [101, "H", None, None, "H"],
    (1002, 1): [102, "He", None, None, "He"],
    (1018, 1): [104, "Ar", None, None, "Ar"],
}


def print_iso():
    print('The dictionary "ISO" contains information on isotopologues in HITRAN\n')
    print(
        "   M    I          id                  iso_name   abundance      mass        mol_name"
    )
    for i in ISO:
        ab = ISO[i][ISO_INDEX["abundance"]]
        ma = ISO[i][ISO_INDEX["mass"]]
        ab = ab if ab else -1
        ma = ma if ma else -1
        print(
            "%4i %4i     : %5i %25s %10f %10f %15s"
            % (
                i[0],
                i[1],
                ISO[i][ISO_INDEX["id"]],
                ISO[i][ISO_INDEX["iso_name"]],
                ab,
                ma,
                ISO[i][ISO_INDEX["mol_name"]],
            )
        )


def print_iso_id():
    print(
        'The dictionary "ISO_ID" contains information on "global" IDs of isotopologues in HITRAN\n'
    )
    print(
        "   id            M    I                    iso_name       abundance       mass        mol_name"
    )
    for i in ISO_ID:
        ab = ISO_ID[i][ISO_ID_INDEX["abundance"]]
        ma = ISO_ID[i][ISO_ID_INDEX["mass"]]
        ab = ab if ab else -1
        ma = ma if ma else -1
        print(
            "%5i     :   %4i %4i   %25s %15.10f %10f %15s"
            % (
                i,
                ISO_ID[i][ISO_ID_INDEX["M"]],
                ISO_ID[i][ISO_ID_INDEX["I"]],
                ISO_ID[i][ISO_ID_INDEX["iso_name"]],
                ab,
                ma,
                ISO_ID[i][ISO_ID_INDEX["mol_name"]],
            )
        )


profiles = "profiles"


def print_profiles():
    print("Profiles available:")
    print("  HT        : PROFILE_HT")
    print("  SDRautian : PROFILE_SDRAUTIAN")
    print("  Rautian   : PROFILE_RAUTIAN")
    print("  SDVoigt   : PROFILE_SDVOIGT")
    print("  Voigt     : PROFILE_VOIGT")
    print("  Lorentz   : PROFILE_LORENTZ")
    print("  Doppler   : PROFILE_DOPPLER")


slit_functions = "slit_functions"


def print_slit_functions():
    print("  RECTANGULAR : SLIT_RECTANGULAR")
    print("  TRIANGULAR  : SLIT_TRIANGULAR")
    print("  GAUSSIAN    : SLIT_GAUSSIAN")
    print("  DIFFRACTION : SLIT_DIFFRACTION")
    print("  MICHELSON   : SLIT_MICHELSON")
    print("  DISPERSION/LORENTZ : SLIT_DISPERSION")


tutorial = "tutorial"
units = "units"
index = "index"
data = "data"
spectra = "spectra"
plotting = "plotting"
python = "python"

python_tutorial_text = """
THIS TUTORIAL IS TAKEN FROM http://www.stavros.io/tutorials/python/
AUTHOR: Stavros Korokithakis


----- LEARN PYTHON IN 10 MINUTES -----


PRELIMINARY STUFF

So, you want to learn the Python programming language but can't find a concise 
and yet full-featured tutorial. This tutorial will attempt to teach you Python in 10 minutes. 
It's probably not so much a tutorial as it is a cross between a tutorial and a cheatsheet, 
so it will just show you some basic concepts to start you off. Obviously, if you want to 
really learn a language you need to program in it for a while. I will assume that you are 
already familiar with programming and will, therefore, skip most of the non-language-specific stuff. 
The important keywords will be highlighted so you can easily spot them. Also, pay attention because, 
due to the terseness of this tutorial, some things will be introduced directly in code and only 
briefly commented on.


PROPERTIES

Python is strongly typed (i.e. types are enforced), dynamically, implicitly typed (i.e. you don't 
have to declare variables), case sensitive (i.e. var and VAR are two different variables) and 
object-oriented (i.e. everything is an object). 


GETTING HELP

Help in Python is always available right in the interpreter. If you want to know how an object works, 
all you have to do is call help(<object>)! Also useful are dir(), which shows you all the object's methods, 
and <object>.__doc__, which shows you its documentation string: 

>>> help(5)
Help on int object:
(etc etc)

>>> dir(5)
['__abs__', '__add__', ...]

>>> abs.__doc__
'abs(number) -> number

Return the absolute value of the argument.'


SYNTAX

Python has no mandatory statement termination characters and blocks are specified by indentation. 
Indent to begin a block, dedent to end one. Statements that expect an indentation level end in a colon (:). 
Comments start with the pound (#) sign and are single-line, multi-line strings are used for multi-line comments. 
Values are assigned (in fact, objects are bound to names) with the _equals_ sign ("="), and equality testing is 
done using two _equals_ signs ("=="). You can increment/decrement values using the += and -= operators respectively 
by the right-hand amount. This works on many datatypes, strings included. You can also use multiple variables on one 
line. For example: 

>>> myvar = 3
>>> myvar += 2
>>> myvar
5

>>> myvar -= 1
>>> myvar
4

\"\"\"This is a multiline comment.
The following lines concatenate the two strings.\"\"\"

>>> mystring = "Hello"
>>> mystring += " world."
>>> print mystring
Hello world.

# This swaps the variables in one line(!).
# It doesn't violate strong typing because values aren't
# actually being assigned, but new objects are bound to
# the old names.
>>> myvar, mystring = mystring, myvar


DATA TYPES

The data structures available in python are lists, tuples and dictionaries. 
Sets are available in the sets library (but are built-in in Python 2.5 and later). 
Lists are like one-dimensional arrays (but you can also have lists of other lists), 
dictionaries are associative arrays (a.k.a. hash tables) and tuples are immutable 
one-dimensional arrays (Python "arrays" can be of any type, so you can mix e.g. integers, 
strings, etc in lists/dictionaries/tuples). The index of the first item in all array types is 0. 
Negative numbers count from the end towards the beginning, -1 is the last item. Variables 
can point to functions. The usage is as follows:

>>> sample = [1, ["another", "list"], ("a", "tuple")]
>>> mylist = ["List item 1", 2, 3.14]
>>> mylist[0] = "List item 1 again" # We're changing the item.
>>> mylist[-1] = 3.21 # Here, we refer to the last item.
>>> mydict = {"Key 1": "Value 1", 2: 3, "pi": 3.14}
>>> mydict["pi"] = 3.15 # This is how you change dictionary values.
>>> mytuple = (1, 2, 3)
>>> myfunction = len
>>> print myfunction(mylist)
3


You can access array ranges using a colon (:). Leaving the start index empty assumes the first item, 
leaving the end index assumes the last item. Negative indexes count from the last item backwards 
(thus -1 is the last item) like so:

>>> mylist = ["List item 1", 2, 3.14]
>>> print mylist[:]
['List item 1', 2, 3.1400000000000001]

>>> print mylist[0:2]
['List item 1', 2]

>>> print mylist[-3:-1]
['List item 1', 2]

>>> print mylist[1:]
[2, 3.14]

# Adding a third parameter, "step" will have Python step in
# N item increments, rather than 1.
# E.g., this will return the first item, then go to the third and
# return that (so, items 0 and 2 in 0-indexing).
>>> print mylist[::2]
['List item 1', 3.14]


STRINGS

Its strings can use either single or double quotation marks, and you can have quotation 
marks of one kind inside a string that uses the other kind (i.e. "He said 'hello'." is valid). 
Multiline strings are enclosed in _triple double (or single) quotes_ (\"\"\"). 
Python supports Unicode out of the box, using the syntax u"This is a unicode string". 
To fill a string with values, you use the % (modulo) operator and a tuple. 
Each %s gets replaced with an item from the tuple, left to right, and you can also use 
dictionary substitutions, like so:

>>>print "Name: %s\
Number: %s\
String: %s" % (myclass.name, 3, 3 * "-")

Name: Poromenos
Number: 3
String: ---

strString = \"\"\"This is
a multiline
string.\"\"\"

# WARNING: Watch out for the trailing s in "%(key)s".
>>> print "This %(verb)s a %(noun)s." % {"noun": "test", "verb": "is"}
This is a test.


FLOW CONTROL STATEMENTS

Flow control statements are if, for, and while. There is no select; instead, use if. 
Use for to enumerate through members of a list. To obtain a list of numbers, 
use range(<number>). These statements' syntax is thus:

rangelist = range(10)
>>> print rangelist
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

>>> for number in rangelist:
        # Check if number is one of
        # the numbers in the tuple.
        if number in (3, 4, 7, 9):
            # "Break" terminates a for without
            # executing the "else" clause.
            break
        else:
            # "Continue" starts the next iteration
            # of the loop. It's rather useless here,
            # as it's the last statement of the loop.
            continue
    else:
        # The "else" clause is optional and is
        # executed only if the loop didn't "break".
        pass # Do nothing

>>> if rangelist[1] == 2:
        print "The second item (lists are 0-based) is 2"
    elif rangelist[1] == 3:
        print "The second item (lists are 0-based) is 3"
    else:
        print "Dunno"

>>> while rangelist[1] == 1:
        pass


FUNCTIONS

Functions are declared with the "def" keyword. Optional arguments are set in 
the function declaration after the mandatory arguments by being assigned a default 
value. For named arguments, the name of the argument is assigned a value. 
Functions can return a tuple (and using tuple unpacking you can effectively return 
multiple values). Lambda functions are ad hoc functions that are comprised of 
a single statement. Parameters are passed by reference, but immutable types (tuples, 
ints, strings, etc) *cannot be changed*. This is because only the memory location of 
the item is passed, and binding another object to a variable discards the old one, 
so immutable types are replaced. For example:

# Same as def funcvar(x): return x + 1
>>> funcvar = lambda x: x + 1
>>> print funcvar(1)
2

# an_int and a_string are optional, they have default values
# if one is not passed (2 and "A default string", respectively).
>>> def passing_example(a_list, an_int=2, a_string="A default string"):
        a_list.append("A new item")
        an_int = 4
        return a_list, an_int, a_string

>>> my_list = [1, 2, 3]
>>> my_int = 10
>>> print passing_example(my_list, my_int)
([1, 2, 3, 'A new item'], 4, "A default string")

>>> my_list
[1, 2, 3, 'A new item']

>>> my_int
10


CLASSES

Python supports a limited form of multiple inheritance in classes. 
Private variables and methods can be declared (by convention, this is not enforced 
by the language) by adding at least two leading underscores and at most one trailing 
one (e.g. "__spam"). We can also bind arbitrary names to class instances. 
An example follows:

>>> class MyClass(object):
        common = 10
        def __init__(self):
            self.myvariable = 3
        def myfunction(self, arg1, arg2):
            return self.myvariable

# This is the class instantiation
>>> classinstance = MyClass()
>>> classinstance.myfunction(1, 2)
3

# This variable is shared by all classes.
>>> classinstance2 = MyClass()
>>> classinstance.common
10

>>> classinstance2.common
10

# Note how we use the class name
# instead of the instance.
>>> MyClass.common = 30
>>> classinstance.common
30

>>> classinstance2.common
30

# This will not update the variable on the class,
# instead it will bind a new object to the old
# variable name.
>>> classinstance.common = 10
>>> classinstance.common
10

>>> classinstance2.common
30

>>> MyClass.common = 50
# This has not changed, because "common" is
# now an instance variable.
>>> classinstance.common
10

>>> classinstance2.common
50

# This class inherits from MyClass. The example
# class above inherits from "object", which makes
# it what's called a "new-style class".
# Multiple inheritance is declared as:
# class OtherClass(MyClass1, MyClass2, MyClassN)
>>> class OtherClass(MyClass):
        # The "self" argument is passed automatically
        # and refers to the class instance, so you can set
        # instance variables as above, but from inside the class.
        def __init__(self, arg1):
            self.myvariable = 3
            print arg1

>>> classinstance = OtherClass("hello")
hello

>>> classinstance.myfunction(1, 2)
3

# This class doesn't have a .test member, but
# we can add one to the instance anyway. Note
# that this will only be a member of classinstance.
>>> classinstance.test = 10
>>> classinstance.test
10


EXCEPTIONS

Exceptions in Python are handled with try-except [exceptionname] blocks:

>>> def some_function():
        try:
            # Division by zero raises an exception
            10 / 0
        except ZeroDivisionError:
            print "Oops, invalid."
        else:
            # Exception didn't occur, we're good.
            pass
        finally:
            # This is executed after the code block is run
            # and all exceptions have been handled, even
            # if a new exception is raised while handling.
            print "We're done with that."

>>> some_function()
Oops, invalid.

We're done with that.


IMPORTING:

External libraries are used with the import [libname] keyword. 
You can also use from [libname] import [funcname] for individual functions. 
Here is an example:

>>> import random
>>> from time import clock

>>> randomint = random.randint(1, 100)
>>> print randomint
64


FILE I/O

Python has a wide array of libraries built in. As an example, here is how serializing 
(converting data structures to strings using the pickle library) with file I/O is used:

>>> import pickle
>>> mylist = ["This", "is", 4, 13327]
# Open the file C:\\binary.dat for writing. The letter r before the
# filename string is used to prevent backslash escaping.
>>> yfile = open(r"C:\\binary.dat", "w")
>>> pickle.dump(mylist, myfile)
>>> myfile.close()

>>> myfile = open(r"C:\\text.txt", "w")
>>> myfile.write("This is a sample string")
>>> myfile.close()

>>> myfile = open(r"C:\\text.txt")
>>> print myfile.read()
'This is a sample string'

>>> myfile.close()

# Open the file for reading.
>>> myfile = open(r"C:\\binary.dat")
>>> loadedlist = pickle.load(myfile)
>>> myfile.close()
>>> print loadedlist
['This', 'is', 4, 13327]


MISCELLANEOUS

    -> Conditions can be chained. 1 < a < 3 checks 
       that a is both less than 3 and greater than 1.
    -> You can use del to delete variables or items in arrays.
    -> List comprehensions provide a powerful way to create 
       and manipulate lists. They consist of an expression 
       followed by a for clause followed by zero or more 
       if or for clauses, like so:

>>> lst1 = [1, 2, 3]
>>> lst2 = [3, 4, 5]
>>> print [x * y for x in lst1 for y in lst2]
[3, 4, 5, 6, 8, 10, 9, 12, 15]

>>> print [x for x in lst1 if 4 > x > 1]
[2, 3]

# Check if a condition is true for any items.
# "any" returns true if any item in the list is true.
>>> any([i % 3 for i in [3, 3, 4, 4, 3]])
True

# This is because 4 % 3 = 1, and 1 is true, so any()
# returns True.

# Check for how many items a condition is true.
>>> sum(1 for i in [3, 3, 4, 4, 3] if i == 4)
2

>>> del lst1[0]
>>> print lst1
[2, 3]

>>> del lst1



    -> Global variables are declared outside of functions 
       and can be read without any special declarations, 
       but if you want to write to them you must declare them 
       at the beginning of the function with the "global" keyword, 
       otherwise Python will bind that object to a new local 
       variable (be careful of that, it's a small catch that can 
       get you if you don't know it). For example:

>>> number = 5

>>> def myfunc():
        # This will print 5.
        print number

>>> def anotherfunc():
        # This raises an exception because the variable has not
        # been bound before printing. Python knows that it an
        # object will be bound to it later and creates a new, local
        # object instead of accessing the global one.
        print number
        number = 3

>>> def yetanotherfunc():
        global number
        # This will correctly change the global.
        number = 3


EPILOGUE

This tutorial is not meant to be an exhaustive list of all (or even a subset) of Python. 
Python has a vast array of libraries and much much more functionality which you will 
have to discover through other means, such as the excellent book Dive into Python. 
I hope I have made your transition in Python easier. Please leave comments if you believe 
there is something that could be improved or added or if there is anything else 
you would like to see (classes, error handling, anything). 

"""


def print_python_tutorial():
    pydoc.pager(python_tutorial_text)


data_tutorial_text = """

ACCESS YOUR DATA!

Welcome to tutorial on retrieving and processing the data from HITRANonline.


  ///////////////
 /// PREFACE ///
///////////////

HITRANonline API is a set of routines in Python which is aimed to 
provide a remote access to functionality and data given by a new project 
HITRANonline (http://hitranazure.cloudapp.net).

At the present moment the API can download, filter and process data on 
molecular and atomic line-by-line spectra which is provided by HITRANonline portal.

One of the major purposes of introducing API is extending a functionality 
of the main site, particularly providing a possibility to calculate several 
types of high- and low-resolution spectra based on a flexible HT lineshape. 

Each feature of API is represented by a Python function with a set of parameters 
providing a flexible approach to the task.


  ///////////////////////
 /// FEATURE SUMMARY ///
///////////////////////

1) Downloading line-by-line data from the HITRANonline site to local database.
2) Filtering and processing the data in SQL-like fashion.
3) Conventional Python structures (lists, tuples, dictionaries) for representing 
   a spectroscopic data.
4) Possibility to use a large set of third-party Python libraries to work with a data
5) Python implementation of an HT (Hartmann-Tran [1]) lineshape which is used in spectra.
   simulations. This lineshape can also be reduced to a number of conventional 
   line profiles such as Gaussian (Doppler), Lorentzian, Voigt, Rautian, 
   Speed-dependent Voigt and Rautian.
6) Python implementation of total internal partition sums (TIPS-2011 [2]) 
   which is used in spectra simulations.
7) High-resolution spectra simulation accounting pressure, 
   temperature and optical path length. The following spectral functions 
   can be calculated:
      a) absorption coefficient
      b) absorption spectrum
      c) transmittance spectrum
      d) radiance spectrum
8) Low-resolution spectra simulation using a number of apparatus functions.
9) Possibility to extend with the user's functionality by adding custom lineshapes, 
   partitions sums and apparatus functions.

References:

[1] N.H. Ngo, D. Lisak, H. Tran, J.-M. Hartmann.
    An isolated line-shape model to go beyond the Voigt profile in 
    spectroscopic databases and radiative transfer codes.
    JQSRT, Volume 129, November 2013, Pages 89–100
    http://dx.doi.org/10.1016/j.jqsrt.2013.05.034

[2] A. L. Laraia, R. R. Gamache, J. Lamouroux, I. E. Gordon, L. S. Rothman.
    Total internal partition sums to support planetary remote sensing.
    Icarus, Volume 215, Issue 1, September 2011, Pages 391–400
    http://dx.doi.org/10.1016/j.icarus.2011.06.004

_______________________________________________________________________


This tutorial will give you an insight of how to use HAPI for Python.

First, let's choose a folder for our local database. Every time you start
your Python project, you have to specify explicitly the name of the 
database folder.

>>> db_begin('data')

So, let's download some data from the server and do some processing on it.
Suppose that we want to get line by line data on the main isotopologue of H2O.

For retrieving the data to the local database, user have to specify the following parameters:
1) Name of the local table which will store the downloaded data.
2) Either a pair of molecule and isotopologue HITRAN numbers (M and I), 
   or a "global" isotopologue ID (iso_id).
3) Wavenumber range (nu_min and nu_max)

N.B. If you specify the name which already exists in the database, 
the existing table with that name will be overrided. 

To get additional information on function fetch,
call getHelp:

>>> getHelp(fetch)
...

To download the data, simply call the function "fetch".
This will establish a connection with the main server and get the data using
the parameters listed above:

>>> fetch('H2O',1,1,3400,4100)
BEGIN DOWNLOAD: H2O
  65536 bytes written to data/H2O.data
  65536 bytes written to data/H2O.data
  65536 bytes written to data/H2O.data
...
  65536 bytes written to data/H2O.data
  65536 bytes written to data/H2O.data
  65536 bytes written to data/H2O.data
Header written to data/H2O.header
END DOWNLOAD
                     Lines parsed: 7524
PROCESSED

The output is shown right after the console line ">>>".
To check the file that you've just downloaded you can open the database
folder. The new plain text file should have a name "H2O.data" and
it should contain line-by-line data in HITRAN format.

N.B. If we want several isotopologues in one table, we should
use fetch_by_ids instead of just fetch. Fetch_by_ids takes a "global" 
isotopologue ID numbers as an input instead of HITRAN's "local" identification.
See getHelp(fetch_by_ids) to get more information on this.

To get a list of tables which are already in the database,
use tableList() function (it takes no arguments):
>>> tableList()

To learn about the table we just downloaded, let's use a function "describeTable".

>>> describeTable('H2O')
-----------------------------------------
H2O summary:
-----------------------------------------
Comment: 
Contains lines for H2(16O)
 in 3400.000-4100.000 wavenumber range
Number of rows: 7524
Table type: column-fixed
-----------------------------------------
            PAR_NAME           PAR_FORMAT

            molec_id                  %2d
        local_iso_id                  %1d
                  nu               %12.6f
                  sw               %10.3E
                   a               %10.3E
           gamma_air                %5.4f
          gamma_self                %5.3f
              elower               %10.4f
               n_air                %4.2f
           delta_air                %8.6f
 global_upper_quanta                 %15s
 global_lower_quanta                 %15s
  local_upper_quanta                 %15s
  local_lower_quanta                 %15s
                ierr                  %6s
                iref                 %12s
    line_mixing_flag                  %1s
                  gp                %7.1f
                 gpp                %7.1f
-----------------------------------------

This output tells how many rows are currenty in the table H2O, which 
wavenumber range was used by fetch(). Also this gives a basic information 
about parameters stored in the table.

So, having the table downloaded, one can perform different operations on it
using API.

Here is a list of operations currently available with API:
1) FILTERING 
2) OUTPUTTING
3) SORTING
4) GROUPING


  ////////////////////////////////
 /// FILTERING AND OUTPUTTING ///
////////////////////////////////

The table data can be filtered with the help of select() function.

Use simple select() call to output the table content:

>>> select('H2O')
MI          nu         S         A gair gsel        E_nair    dair  ...
11 1000.288940 1.957E-24 2.335E-02.07100.350 1813.22270.680.008260  ...
11 1000.532321 2.190E-28 1.305E-05.04630.281 2144.04590.39-.011030  ...
...

This will display the list of line parameters containing in the table "H2O".

That's the simplest way of using the function select(). Full information
on control parameters can be obtained via getHelp(select) statement.

Suppose that we need a lines from a table within some wavenumber range. 
That's what filtering is for. Let's apply a simple range filter on a table.

>>> select('H2O',Conditions=('between','nu',4000,4100))
MI          nu         S         A gair gsel        E_nair    dair     
 11 4000.188800 1.513E-25 1.105E-02.03340.298 1581.33570.51-.013910 ...
 11 4000.204070 3.482E-24 8.479E-03.08600.454  586.47920.61-.007000 ...
 11 4000.469910 3.268E-23 1.627E+00.05410.375 1255.91150.56-.013050 ...
......

As a result of this operation, we see a list of lines of H2O table,
whose wavenumbers lie between 4000 cm-1 and 4100 cm-1.
The condition is taken as an input parameter to API function "select".

To specify a subset of columns to display, use another control parameter - 
ParameterNames:

>>> select('H2O',ParameterNames=('nu','sw'),Conditions=('between','nu',4000,4100))

The usage of ParameterNames is outlined below in the section "Specifying a list 
of parameters". So far it worth mentioning that this parameter is a part 
of a powerful tool for displaying and processing tables from database.

In the next section we will show how to create quieries 
with more complex conditions.


  ////////////////////////////
 /// FILTERING CONDITIONS ///
////////////////////////////

Let's analyze the last example of filtering. Condition input variable is
as follows:

                    ('between','nu',4000,4100)

Thus, this is a python list (or tuple), containing logical expressions
defined under column names of the table. For example, 'nu' is a name of 
the column in 'H2O' table, and this column contains a transition wavenumber.
The structure of a simple condition is as follows:

                    (OPERATION,ARG1,ARG2,...)
                    
Where OPERATION must be in a set of predefined operations (see below),
and ARG1,ARG2 etc. are the arguments for this operation.
Conditions can be nested, i.e. ARG can itself be a condition (see examples).
The following operations are available in select (case insensitive):


DESCRIPTION                   LITERAL                     EXAMPLE
---------------------------------------------------------------------------------
Range:               'RANGE','BETWEEN':         ('BETWEEN','nu',0,1000)
Subset:              'IN','SUBSET':             ('IN','local_iso_id',[1,2,3,4])
And:                 '&','&&','AND':            ('AND',('<','nu',1000),('>','nu',10))
Or:                  '|','||','OR':             ('OR',('>','nu',1000),('<','nu',10))
Not:                 '!','NOT':                 ('NOT',('IN','local_iso_id',[1,2,3]))
Less than:           '<','LESS','LT':                 ('<','nu',1000)
More than:           '>','MORE','MT':                 ('>','sw',1.0e-20)
Less or equal than:  '<=','LESSOREQUAL','LTE':        ('<=','local_iso_id',10)
More or equal than   '>=','MOREOREQUAL','MTE':        ('>=','sw',1e-20)
Equal:               '=','==','EQ','EQUAL','EQUALS':  ('<=','local_iso_id',10)
Not equal:           '!=','<>','~=','NE','NOTEQUAL':  ('!=','local_iso_id',1)
Summation:           '+','SUM':                 ('+','v1','v2','v3')
Difference:          '-','DIFF':                ('-','nu','elow')
Multiplication:      '*','MUL':                 ('*','sw',0.98)
Division:            '/','DIV':                 ('/','A',2)
Cast to string:      'STR','STRING':            ('STR','some_string')
Cast to Python list  'LIST':                    ('LIST',[1,2,3,4,5])
Match regexp         'MATCH','LIKE':            ('MATCH','\w+','some string')
Search single match: 'SEARCH':                  ('SEARCH','\d \d \d','1 2 3 4')
Search all matches:  'FINDALL':                 ('FINDALL','\d','1 2 3 4 5')
Count within group:  'COUNT' :                  ('COUNT','local_iso_id')
---------------------------------------------------------------------------------
   
Let's create a query with more complex condition. Suppese that we are 
interested in all lines between 3500 and 4000 with 1e-19 intensity cutoff.
The query will look like this:

>>> Cond = ('AND',('BETWEEN','nu',3500,4000),('>=','Sw',1e-19))
>>> select('H2O',Conditions=Cond,DestinationTableName='tmp')

Here, apart from other parameters, we have used a new parameter 
DestinationTableName. This parameter contains a name of the table
where we want to put a result of the query. Thus we have chosen 
a name 'tmp' for a new table.


  ////////////////////////////////////
 /// ACCESSING COLUMNS IN A TABLE ///
////////////////////////////////////

To get an access to particular table column (or columns) all we need
is to get a column from a table and put it to Python variable.

For this purpose, there exist two functions:

  getColumn(...)
  getColumns(...)

The first one returns just one column at a time. The second one returns
a list of solumns.

So, here are some examples of how to use both:

>>> nu1 = getColumn('H2O','nu')
>>> nu2,sw2 = getColumns('H2O',['nu','sw'])

N.B. If you don't remember exact names of columns in a particular table,
use describeTable to get an info on it's structure!


  ///////////////////////////////////////
 /// SPECIFYING A LIST OF PARAMETERS ///
///////////////////////////////////////

Suppose that we want not only select a set of parameters/columns
from a table, but do a certain transformations with them (for example,
multiply column on a coefficient, or add one column to another etc...).
We can make it in two ways. First, we can extract a column from table
using one of the functions (getColumn or getColumns) and do the rest 
in Python. The second way is to do it on the level of select.
The select function has a control parameter "ParameterNames", which 
makes it possible to specify parameters we want to be selected, 
and evaluate some simple arithmetic expressions with them.

Assume that we need only wavenumber and intensity from H2O table.
Also we need to scale an intensity to the unitary abundance. To do so,
we must divide an 'sw' parameter by it's natural abundance (0.99731) for 
principal isotopologue of water).

Thus, we have to select two columns:  
wavenumber (nu) and scaled intensity (sw/0.99731)
>>> select('H2O',)


  ////////////////////////////
 /// SAVING QUERY TO DISK ///
////////////////////////////

To quickly save a result of a query to disk, the user can take an 
advantage of an additional parameter "File".
If this parameter is presented in function call, then the query is 
saved to file with the name which was specified in "File".

For example, select all lines from H2O and save the result in file 'H2O.txt':
>>> select('H2O',File='H2O.txt')


  ////////////////////////////////////////////
 /// GETTING INFORMATION ON ISOTOPOLOGUES ///
////////////////////////////////////////////

API provides the following auxillary information about isotopologues
present in HITRAN. Corresponding functions use the standard HITRAN
molecule-isotopologue notation:

1) Natural abundances
>>> abundance(mol_id,iso_id)

2) Molecular masses
>>> molecularMass(mol_id,iso_id)

3) Molecule names
>>> moleculeName(mol_id,iso_id)

4) Isotopologue names
>>> isotopologueName(mol_id,iso_id)

5) ISO_ID
>>> getHelp(ISO_ID)

The latter is a dictionary, which contain all information about 
isotopologues concentrated in one place.

"""


def print_data_tutorial():
    pydoc.pager(data_tutorial_text)


spectra_tutorial_text = """

CALCULATE YOUR SPECTRA!

Welcome to tutorial on calculating a spectra from line-by-line data.


  ///////////////
 /// PREFACE ///
///////////////

This tutorial will demonstrate how to use different lineshapes and partition
functions, and how to calculate synthetic spectra with respect to different 
instruments. It will be shown how to combine different parameters of spectral 
calculation to achieve better precision and performance for cross sections.

API provides a powerful tool to calculate cross-sections based on line-by-line
data containing in HITRAN. This features:

*) Python implementation of an HT (Hartmann-Tran [1]) lineshape which is used in 
   spectra simulations. This lineshape can also be reduced to a number of 
   conventional    line profiles such as Gaussian (Doppler), Lorentzian, Voigt, 
   Rautian, Speed-dependent Voigt and Rautian.
*) Python implementation of total internal partition sums (TIPS-2011 [2]) 
   which is used in spectra simulations.
*) High-resolution spectra simulation accounting pressure, 
   temperature and optical path length. The following spectral functions 
   can be calculated:
      a) absorption coefficient
      b) absorption spectrum
      c) transmittance spectrum
      d) radiance spectrum
*) Low-resolution spectra simulation using a number of apparatus functions.
*) Possibility to extend with the user's functionality by adding custom lineshapes, 
   partitions sums and apparatus functions.
*) An approach to function code is aimed to be flexible enough yet hopefully 
   intuitive.

References:

[1] N.H. Ngo, D. Lisak, H. Tran, J.-M. Hartmann.
    An isolated line-shape model to go beyond the Voigt profile in 
    spectroscopic databases and radiative transfer codes.
    JQSRT, Volume 129, November 2013, Pages 89–100
    http://dx.doi.org/10.1016/j.jqsrt.2013.05.034

[2] A. L. Laraia, R. R. Gamache, J. Lamouroux, I. E. Gordon, L. S. Rothman.
    Total internal partition sums to support planetary remote sensing.
    Icarus, Volume 215, Issue 1, September 2011, Pages 391–400
    http://dx.doi.org/10.1016/j.icarus.2011.06.004

            
  ///////////////////////////
 /// USING LINE PROFILES ///
///////////////////////////

Several lineshape (line profile) families are currently available:
1) Gaussian (Doppler) profile
2) Lorentzian profile
3) Voigt profile
4) Speed-dependent Voigt profile
5) Rautian profile
6) Speed-dependent Rautian profile
7) HT profile (Hartmann-Tran)

Each profile has it's own uniwue set of parameters. Normally one should
use profile parameters only in conjunction with their "native" profiles.

So, let's start exploring the available profiles using getHelp:
>>> getHelp(profiles)
Profiles available:
  HTP       : PROFILE_HT
  SDRautian : PROFILE_SDRAUTIAN
  Rautian   : PROFILE_RAUTIAN
  SDVoigt   : PROFILE_SDVOIGT
  Voigt     : PROFILE_VOIGT
  Lorentz   : PROFILE_LORENTZ
  Doppler   : PROFILE_DOPPLER

Output gives all available profiles. We can get additional info on each
of them just by calling getHelp(ProfileName):
>>> getHelp(PROFILE_HT)

Line profiles, adapted for using with HAPI, are written in Python and
heavily using the numerical library "Numpy". This means that the user
can calculate multiple values of particular profile at once having just
pasted a numpy array as a wavenumber grid (array). Let's give a short 
example of how to calculate HT profile on a numpy array.

>>> from numpy import arange
    w0 = 1000.
    GammaD = 0.005
    Gamma0 = 0.2
    Gamma2 = 0.01 * Gamma0
    Delta0 = 0.002
    Delta2 = 0.001 * Delta0
    nuVC = 0.2
    eta = 0.5
    Dw = 1.
    ww = arange(w0-Dw, w0+Dw, 0.01)  # GRID WITH THE STEP 0.01 
    l1 = PROFILE_HT(w0,GammaD,Gamma0,Gamma2,Delta0,Delta2,nuVC,eta,ww)[0]
    # now l1 contains values of HT profile calculates on the grid ww
    
On additional information about parameters see getHelp(PROFILE_HT).

It worth noting that PROFILE_HT returns 2 entities: real and imaginary part
of lineshape (as it described in the article given in preface). Apart from
HT, all other profiles return just one entity (the real part).


  ////////////////////////////
 /// USING PARTITION SUMS ///
////////////////////////////

As it was mentioned in the preface to this tutorial, the partition sums
are taken from the TIPS-2011 (the link is given above). Partition sums 
are taken for those isotopologues, which are present in HITRAN and in
TIPS-2011 simultaneousely.

N.B. Partition sums are omitted for the following isotopologues which
are in HITRAN at the moment:

ID       M     I         ISO                MOL
--------------------------------------------------
117      12    2     H(15N)(16O)3           HNO3
110      14    2     D(19F)                 HF
107      15    3     D(35Cl)                HCl
108      15    4     D(37Cl)                HCl
111      16    3     D(79Br)                HBr
112      16    4     D(81Br)                HBr
113      17    2     D(127I)                HI
118      22    2     (14N)(15N)             N2
119      29    2     (13C)(16O)(19F)2       COF2
 86      34    1     (16O)                  O
 92      39    1     (12C)H3(16O)H          CH3OH
114      47    1     (32S)(16O)3            SO3
--------------------------------------------------

The data on these isotopologues is not present in TIPS-2011 but is 
present in HITRAN. We're planning to add these molecules after TIPS-2013
is released.

To calculate a partition sum for most of the isotopologues in HITRAN,
we will use a function partitionSum (use getHelp for detailed info).
Let's just mention that 
The syntax is as follows: partitionSum(M,I,T), where M,I - standard 
HITRAN molecule-isotopologue notation, T - definition of temperature
range.

Usecase 1: temperatuer is defined by a list:
>>> Q = partitionSum(1,1,[70,80,90])

Usecase 2: temperature is defined by bounds and the step:
>>> T,Q = partiionSum(1,1,[70,3000],step=1.0)

In the latter example we calculate a partition sum on a range of
temperatures from 70K to 3000K using a step 1.0 K, and having arrays 
of temperature (T) and partition sum (Q) at the output.


  ///////////////////////////////////////////
 /// CALCULATING ABSORPTION COEFFICIENTS ///
///////////////////////////////////////////

Currently API can calculate the following spectral function at arbitrary
thermodynamic parameters:

1) Absorption coefficient
2) Absorption spectrum
3) Transmittance spectrum
4) Radiance spectrum

All these functions can be calculated with or without accounting of 
an instrument properties (apparatus function, resolution, path length etc...)

As it well known, the spectral functions such as absorption,
transmittance, and radiance spectra, are calculated on the basis
of the absorption coefficient. By that resaon, absorption coefficient
is the most important part of simulating a cross section. This part of
tutorial is devoted to demonstration how to calculate absorption 
coefficient from the HITRAN line-by-line data. Here we give a brief 
insight on basic parameters of calculation procedure, talk about some 
useful practices and precautions.

To calculate an absorption coefficient, we can use one of the following
functions:

-> absorptionCoefficient_HT
-> absorptionCoefficient_Voigt
-> absorptionCoefficient_Lorentz
-> absorptionCoefficient_Doppler

Each of these function calculates cross sections using different
lineshapes (the names a quite self-explanatory).
You can get detailed information on using each of these functions
by calling getHelp(function_name).

Let's look more closely to the cross sections based on the Lorentz profile.
For doing that, let's have a table downloaded from HITRANonline.

# get data on CO2 main isotopologue in the range 2000-2100 cm-1
>>> fetch('CO2',2,1,2000,2100)

OK, now we're ready to run a fast example of how to calculate an
absorption coefficient cross section:

>>> nu,coef = absorptionCoefficient_Lorentz(SourceTables='CO2')

This example calculates a Lorentz cross section using the whole set of 
lines in the "co2" table. This is the simplest possible way to use these
functions, because major part of parameters bound to their default values.

If we have matplotlib installed, then we can visualize it using a plotter:
>>> from pylab import plot
>>> plot(nu,coef) 

API provides a flexible control over a calculation procedure. This control
can be achieved by using a number of input parameters. So, let's dig 
into the depth of the settings.

The input parameters of absorptionCoefficient_Lorentz are as follows:

Name                          Default value
-------------------------------------------------------------------
SourceTables                  '__BUFFER__'
Components                    All isotopologues in SourceTables 
partitionFunction             PYTIPS
Environment                   {'T':296.,'p':1.}
WavenumberRange               depends on Components
WavenumberStep                0.01 cm-1
WavenumberWing                10 cm-1
WavenumberWingHW              50 HWHMs
IntensityThreshold            0 cm/molec
GammaL                        'gamma_air'
HITRAN_units                  True 
File                          None
Format                        '%e %e'
-------------------------------------------------------------------

Newt we'll give a brief explanation for each parameter. After each description
we'll make some notes about the usage of the correspondent parameter.


SourceTables:     (required parameter)
   
  List of source tables to take line-by-line data from.
  NOTE: User must provide at least one table in the list.

Components:    (optional parameter)

  List of tuples (M,I,D) to consider in cross section calculation.
  M here is a molecule number, I is an isotopologue number, 
  D is an abundance of the component.
  NOTE: If this input contains more than one tuple, then the output 
        is an absorption coefficient for mixture of corresponding gases.
  NOTE2: If omitted, then all data from the source tables is involved.

partitionFunction:    (optional parameter)

  Instance of partition function of the following format:
  Func(M,I,T), where Func - numae of function, (M,I) - HITRAN numbers
  for molecule and isotopologue, T - temperature.
  Function must return only one output - value of partition sum.
  NOTE: Deafult value is PYTIPS - python version of TIPS-2011

Environment:    (optional parameter)

  Python dictionary containing value of pressure and temperature.
  The format is as follows: Environment = {'p':pval,'T':tval}, 
  where "pval" and "tval" are corresponding values in atm and K 
  respectively.
  NOTE: Default value is {'p':1.0,'T':296.0}

WavenumberRange:    (optional parameter)

  List containing minimum and maximum value of wavenumber to consider
  in cross-section calculation. All lines that are out of htese bounds
  will be skipped. The firmat is as follows: WavenumberRange=[wn_low,wn_high]
  NOTE: If this parameter os skipped, then min and max are taken 
  from the data from SourceTables. Deprecated name is OmegaRange.

WavenumberStep:    (optional parameter)

  Value for the wavenumber step. 
  NOTE: Default value is 0.01 cm-1.
  NOTE2: Normally user would want to take the step under 0.001 when
         calculating absorption coefficient with Doppler profile 
         because of very narrow spectral lines. Deprecated name is OmegaStep.

WavenumberWing:    (optional parameter)

  Absolute value of the line wing in cm-1, i.e. distance from the center 
  of each line to the most far point where the profile is considered 
  to be non zero. Deprecated name is OmegaStep.
  NOTE: if omitted, then only OmegaWingHW is taken into account.

WavenumberWingHW:    (optional parameter)

  Relative value of the line wing in halfwidths. Deprecated name is OmegaWingHW.
  NOTE: The resulting wing is a maximum value from both OmegaWing and
  OmegaWingHW.

IntensityThreshold:    (optional parameter)

  Absolute value of minimum intensity in cm/molec to consider.
  NOTE: default value is 0.

GammaL:    (optional parameter)

  This is the name of broadening parameter to consider a "Lorentzian"
  part in the Voigt profile. In the current 160-char format there is 
  a choise between "gamma_air" and "gamma_self".
  NOTE: If the table has custom columns with a broadening coefficients,
        the user can specify the name of this column in GammaL. This
        would let the function calculate an absorption with custom
        broadening parameter.

HITRAN_units:    (optional parameter)

  Logical flag for units, in which the absorption coefficient shoould be 
  calculated. Currently, the choises are: cm^2/molec (if True) and
  cm-1 (if False).
  NOTE: to calculate other spectral functions like transmitance,
  radiance and absorption spectra, user should set HITRAN_units to False.

File:    (optional parameter)

  The name of the file to save the calculated absorption coefficient.
  The file is saved only if this parameter is specified.

Format:    (optional parameter)

  C-style format for the text data to be saved. Default value is "%e %e".
  NOTE: C-style output format specification (which are mostly valid for Python) 
        can be found, for instance, by the link: 
  http://www.gnu.org/software/libc/manual/html_node/Formatted-Output.html


N.B. Other functions such as absorptionCoefficient_Voigt(_HT,_Doppler) have
identical parameter sets so the description is the same for each function.


  ///////////////////////////////////////////////////////////////////
 /// CALCULATING ABSORPTION, TRANSMITTANCE, AND RADIANCE SPECTRA ///
///////////////////////////////////////////////////////////////////

Let's calculate an absorption, transmittance, and radiance
spectra on the basis of apsorption coefficient. In order to be consistent
with internal API's units, we need to have an absorption coefficient cm-1:

>>> nu,coef = absorptionCoefficient_Lorentz(SourceTables='CO2',HITRAN_units=False)

To calculate absorption spectrum, use the function absorptionSpectrum():
>>> nu,absorp = absorptionSpectrum(nu,coef) 

To calculate transmittance spectrum, use function transmittanceSpectrum():
>>> nu,trans = transmittanceSpectrum(nu,coef) 

To calculate radiance spectrum, use function radianceSpectrum():
>>> nu,radi = radianceSpectrum(nu,coef) 


The last three commands used a default path length (1 m).
To see complete info on all three functions, look for section
"calculating spectra" in getHelp()

Generally, all these three functions use similar set of parameters:

Wavenumber:       (required parameter) 

  Wavenumber grid to for spectrum. Deprecated name is Omegas.

AbsorptionCoefficient        (optional parameter)

  Absorption coefficient as input.

Environment={'T': 296.0, 'l': 100.0}       (optional parameter) 

  Environmental parameters for calculating  spectrum.
  This parameter is a bit specific for each of functions:
  For absorptionSpectrum() and transmittanceSpectrum() the default
  value is as follows: Environment={'l': 100.0}
  For transmittanceSpectrum() the default value, besides path length,
  contains a temperature: Environment={'T': 296.0, 'l': 100.0}
  NOTE: temperature must be equal to that which was used in 
  absorptionCoefficient_ routine!

File         (optional parameter)

  Filename of output file for calculated spectrum.
  If omitted, then the file is not created.

Format        (optional parameter)

  C-style format for spectra output file.
  NOTE: Default value is as follows: Format='%e %e'


  ///////////////////////////////////////
 /// APPLYING INSTRUMENTAL FUNCTIONS ///
///////////////////////////////////////

For comparison of the theoretical spectra with the real-world 
instruments output it's necessary to take into account instrumental resolution.
For this purpose HAPI has a function convolveSpectrum() which can emulate
spectra with lower resolution using custom instrumental functions.

The following instrumental functions are available:
1) Rectangular
2) Triangular
3) Gaussian
4) Diffraction
5) Michelson
6) Dispersion
7) Lorentz

To get a description of each instrumental function we can use getHelp():
>>> getHelp(slit_functions)
  RECTANGULAR : SLIT_RECTANGULAR
  TRIANGULAR  : SLIT_TRIANGULAR
  GAUSSIAN    : SLIT_GAUSSIAN
  DIFFRACTION : SLIT_DIFFRACTION
  MICHELSON   : SLIT_MICHELSON
  DISPERSION/LORENTZ : SLIT_DISPERSION
  
For instance,
>>> getHelp(SLIT_MICHELSON)
... will give a datailed info about Michelson's instrumental function.


The function convolveSpectrum() convolutes a high-resulution spectrum
with one of supplied instrumental (slit) functions. The folowing 
parameters of this function are provided:

Wavenumber     (required parameter)
  
  Array of wavenumbers in high-resolution input spectrum.
  Deprecated name is Omega.

CrossSection     (required parameter)

  Values of high-resolution input spectrum.

Resolution     (optional parameter)

  This parameter is passed to the slit function. It represents
  the resolution of corresponding instrument.
  NOTE: default value is 0.1 cm-1

AF_wing     (optional parameter)

  Width of an instrument function where it is considered non-zero.
  NOTE: default value is 10.0 cm-1

SlitFunction     (optional parameter)

  Custom instrumental function to convolve with spectrum.
  Format of the instrumental function must be as follows:
  Func(x,g), where Func - function name, x - wavenumber,
  g - resolution.
  NOTE: if omitted, then the default value is SLIT_RECTANGULAR


Before using the convolution procedure it worth giving some practical 
advices and remarks: 
1) Quality of a convolution depends on many things: quality of calculated 
spectra, width of AF_wing and WavenumberRange, Resolution, WavenumberStep etc ...
Most of these factors are taken from previus stages of spectral calculation.
Right choise of all these factors is crucial for the correct computation.
2) Dispersion, Diffraction and Michelson AF's don't work well in narrow 
wavenumber range because of their broad wings.
3) Generally one must consider WavenumberRange and AF_wing as wide as possible.
4) After applying a convolution, the resulting spectral range for 
the lower-resolution spectra is reduced by the doubled value of AF_wing.
For this reason, try to make an initial spectral range for high-resolution
spectrum (absorption, transmittance, radiance) sufficiently broad.

The following command will calculate a lower-resolution spectra from 
the CO2 transmittance, which was calculated in a previous section. 
The Spectral resolution is 1 cm-1, 

>>> nu_,trans_,i1,i2,slit = convolveSpectrum(nu,trans)

The outputs are: 

nu_, trans_ - wavenumbers and transmittance for the resulting 
              low-resolution spectrum.

i1,i2 - indexes for initial nu,trans spectrum denoting the part of 
        wavenumber range which was taken for lower resolution spectrum.
        => Low-res spectrum is calculated on nu[i1:i2]

Note, than to achieve more flexibility, one have to specify most of 
the optional parameters. For instance, more complete call is as follows:
>>> nu_,trans_,i1,i2,slit = convolveSpectrum(nu,trans,SlitFunction=SLIT_MICHELSON,Resolution=1.0,AF_wing=20.0)

"""


def print_spectra_tutorial():
    pydoc.pager(spectra_tutorial_text)


plotting_tutorial_text = """

PLOTTING THE SPECTRA WITH MATPLOTLIB

This tutorial briefly explains how to make plots using
the Matplotlib - Python library for plotting.

Prerequisites:
   To tun through this tutorial, user must have the following
   Python libraries installed:
   1) Matplotlib
       Matplotlib can be obtained by the link http://matplotlib.org/ 
   2) Numpy  (required by HAPI itself)
       Numpy can be obtained via pip:  
          sudo pip install numpy (under Linux and Mac)
          pip install numpy (under Windows)
       Or by the link http://www.numpy.org/
       
As an option, user can download one of the many scientific Python
distributions, such as Anaconda, Canopy etc...

So, let's calculate plot the basic entities which ar provided by HAPI.
To do so, we will do all necessary steps to download, filter and 
calculate cross sections "from scratch". To demonstrate the different
possibilities of matplotlib, we will mostly use Pylab - a part of 
Matplotlib with the interface similar to Matlab. Please note, that it's 
not the only way to use Matplotlib. More information can be found on it's site.

The next part is a step-by-step guide, demonstrating basic possilities
of HITRANonline API in conjunction with Matplotlib.

First, do some preliminary imports:
>>> from hapi import *
>>> from pylab import show,plot,subplot,xlim,ylim,title,legend,xlabel,ylabel,hold

Start the database 'data':
>>> db_begin('data') 

Download lines for main isotopologue of ozone in [3900,4050] range:
>>> fetch('O3',3,1,3900,4050)

PLot a sick spectrum using the function getStickXY()
>>> x,y = getStickXY('O3')
>>> plot(x,y); show()

Zoom in spectral region [4020,4035] cm-1:
>>> plot(x,y); xlim([4020,4035]); show()

Calculate and plot difference between Voigt and Lorentzian lineshape:
>>> wn = arange(3002,3008,0.01) # get wavenumber range of interest
>>> voi = PROFILE_VOIGT(3005,0.1,0.3,wn)[0]   # calc Voigt
>>> lor = PROFILE_LORENTZ(3005,0.3,wn)   # calc Lorentz
>>> diff = voi-lor    # calc difference
>>> subplot(2,1,1)   # upper panel
>>> plot(wn,voi,'red',wn,lor,'blue')  # plot both profiles
>>> legend(['Voigt','Lorentz'])   # show legend
>>> title('Voigt and Lorentz profiles')   # show title
>>> subplot(2,1,2)   # lower panel
>>> plot(wn,diff)   # plot diffenence
>>> title('Voigt-Lorentz residual')   # show title
>>> show()   # show all figures

Calculate and plot absorption coefficients for ozone using Voigt 
profile. Spectra are calculated for 4 cases of thermodynamic parameters: 
(1 atm, 296 K), (5 atm, 296 K), (1 atm, 500 K), and (5 atm, 500 K)
>>> nu1,coef1 = absorptionCoefficient_Voigt(((3,1),),'O3',
        WavenumberStep=0.01,HITRAN_units=False,GammaL='gamma_self',
        Environment={'p':1,'T':296.})
>>> nu2,coef2 = absorptionCoefficient_Voigt(((3,1),),'O3',
        WavenumberStep=0.01,HITRAN_units=False,GammaL='gamma_self',
        Environment={'p':5,'T':296.})
>>> nu3,coef3 = absorptionCoefficient_Voigt(((3,1),),'O3',
        WavenumberStep=0.01,HITRAN_units=False,GammaL='gamma_self',
        Environment={'p':1,'T':500.})
>>> nu4,coef4 = absorptionCoefficient_Voigt(((3,1),),'O3',
        WavenumberStep=0.01,HITRAN_units=False,GammaL='gamma_self',
        Environment={'p':5,'T':500.})
>>> subplot(2,2,1); plot(nu1,coef1); title('O3 k(w): p=1 atm, T=296K')
>>> subplot(2,2,2); plot(nu2,coef2); title('O3 k(w): p=5 atm, T=296K')
>>> subplot(2,2,3); plot(nu3,coef3); title('O3 k(w): p=1 atm, T=500K')
>>> subplot(2,2,4); plot(nu4,coef4); title('O3 k(w): p=5 atm, T=500K')
>>> show()

Calculate and plot absorption, transmittance and radiance spectra for 1 atm 
and 296K. Path length is set to 10 m.
>>> nu,absorp = absorptionSpectrum(nu1,coef1,Environment={'l':1000.})
>>> nu,transm = transmittanceSpectrum(nu1,coef1,Environment={'l':1000.})
>>> nu,radian = radianceSpectrum(nu1,coef1,Environment={'l':1000.,'T':296.})
>>> subplot(2,2,1); plot(nu1,coef1,'r'); title('O3 k(w): p=1 atm, T=296K')
>>> subplot(2,2,2); plot(nu,absorp,'g'); title('O3 absorption: p=1 atm, T=296K')
>>> subplot(2,2,3); plot(nu,transm,'b'); title('O3 transmittance: p=1 atm, T=296K')
>>> subplot(2,2,4); plot(nu,radian,'y'); title('O3 radiance: p=1 atm, T=296K')
>>> show()

Calculate and compare high resolution spectrum for O3 with lower resolution
spectrum convoluted with an instrumental function of ideal Michelson interferometer.
>>> nu_,trans_,i1,i2,slit = convolveSpectrum(nu,transm,SlitFunction=SLIT_MICHELSON,Resolution=1.0,AF_wing=20.0)
>>> plot(nu,transm,'red',nu_,trans_,'blue'); legend(['HI-RES','Michelson']); show()

"""


def print_plotting_tutorial():
    pydoc.pager(plotting_tutorial_text)


def getHelp(arg=None):
    """
    This function provides interactive manuals and tutorials.
    """
    if arg == None:
        print("--------------------------------------------------------------")
        print("Hello, this is an interactive help system of HITRANonline API.")
        print("--------------------------------------------------------------")
        print("Run getHelp(.) with one of the following arguments:")
        print("    tutorial  -  interactive tutorials on HAPI")
        print("    units     -  units used in calculations")
        print("    index     -  index of available HAPI functions")
    elif arg == "tutorial":
        print("-----------------------------------")
        print("This is a tutorial section of help.")
        print("-----------------------------------")
        print("Please choose the subject of tutorial:")
        print("    data      -  downloading the data and working with it")
        print("    spectra   -  calculating spectral functions")
        print("    plotting  -  visualizing data with matplotlib")
        print("    python    -  Python quick start guide")
    elif arg == "python":
        print_python_tutorial()
    elif arg == "data":
        print_data_tutorial()
    elif arg == "spectra":
        print_spectra_tutorial()
    elif arg == "plotting":
        print_plotting_tutorial()
    elif arg == "index":
        print("------------------------------")
        print("FETCHING DATA:")
        print("------------------------------")
        print("  fetch")
        print("  fetch_by_ids")
        print("")
        print("------------------------------")
        print("WORKING WITH DATA:")
        print("------------------------------")
        print("  db_begin")
        print("  db_commit")
        print("  tableList")
        print("  describe")
        print("  select")
        print("  sort")
        print("  extractColumns")
        print("  getColumn")
        print("  getColumns")
        print("  dropTable")
        print("")
        print("------------------------------")
        print("CALCULATING SPECTRA:")
        print("------------------------------")
        print("  profiles")
        print("  partitionSum")
        print("  absorptionCoefficient_HT")
        print("  absorptionCoefficient_Voigt")
        print("  absorptionCoefficient_SDVoigt")
        print("  absorptionCoefficient_Lorentz")
        print("  absorptionCoefficient_Doppler")
        print("  transmittanceSpectrum")
        print("  absorptionSpectrum")
        print("  radianceSpectrum")
        print("")
        print("------------------------------")
        print("CONVOLVING SPECTRA:")
        print("------------------------------")
        print("  convolveSpectrum")
        print("  slit_functions")
        print("")
        print("------------------------------")
        print("INFO ON ISOTOPOLOGUES:")
        print("------------------------------")
        print("  ISO_ID")
        print("  abundance")
        print("  molecularMass")
        print("  moleculeName")
        print("  isotopologueName")
        print("")
        print("------------------------------")
        print("MISCELLANEOUS:")
        print("------------------------------")
        print("  getStickXY")
        print("  read_hotw")
    elif arg == ISO:
        print_iso()
    elif arg == ISO_ID:
        print_iso_id()
    elif arg == profiles:
        print_profiles()
    elif arg == slit_functions:
        print_slit_functions()
    else:
        help(arg)


# Get atmospheric (natural) abundance
# for a specified isotopologue
# M - molecule number
# I - isotopologue number
def abundance(M, I):
    """
    INPUT PARAMETERS: 
        M: HITRAN molecule number
        I: HITRAN isotopologue number
    OUTPUT PARAMETERS: 
        Abbundance: natural abundance
    ---
    DESCRIPTION:
        Return natural (Earth) abundance of HITRAN isotolopogue.
    ---
    EXAMPLE OF USAGE:
        ab = abundance(1,1) # H2O
    ---
    """
    return ISO[(M, I)][ISO_INDEX["abundance"]]


# Get molecular mass
# for a specified isotopologue
# M - molecule number
# I - isotopologue number


def molecularMass(M, I):
    """
    INPUT PARAMETERS: 
        M: HITRAN molecule number
        I: HITRAN isotopologue number
    OUTPUT PARAMETERS: 
        MolMass: molecular mass
    ---
    DESCRIPTION:
        Return molecular mass of HITRAN isotolopogue.
    ---
    EXAMPLE OF USAGE:
        mass = molecularMass(1,1) # H2O
    ---
    """
    return ISO[(M, I)][ISO_INDEX["mass"]]


# Get molecule name
# for a specified isotopologue
# M - molecule number
# I - isotopologue number


def moleculeName(M):
    """
    INPUT PARAMETERS: 
        M: HITRAN molecule number
    OUTPUT PARAMETERS: 
        MolName: molecular name
    ---
    DESCRIPTION:
        Return name of HITRAN molecule.
    ---
    EXAMPLE OF USAGE:
        molname = moleculeName(1) # H2O
    ---
    """
    return ISO[(M, 1)][ISO_INDEX["mol_name"]]


# Get isotopologue name
# for a specified isotopologue
# M - molecule number
# I - isotopologue number


def isotopologueName(M, I):
    """
    INPUT PARAMETERS: 
        M: HITRAN molecule number
        I: HITRAN isotopologue number
    OUTPUT PARAMETERS: 
        IsoMass: isotopologue mass
    ---
    DESCRIPTION:
        Return name of HITRAN isotolopogue.
    ---
    EXAMPLE OF USAGE:
        isoname = isotopologueName(1,1) # H2O
    ---
    """
    return ISO[(M, I)][ISO_INDEX["iso_name"]]


# ----------------------- table list ----------------------------------


def tableList():
    """
    INPUT PARAMETERS: 
        none
    OUTPUT PARAMETERS: 
        TableList: a list of available tables
    ---
    DESCRIPTION:
        Return a list of tables present in database.
    ---
    EXAMPLE OF USAGE:
        lst = tableList()
    ---
    """

    return getTableList()


# ----------------------- describe ----------------------------------


def describe(TableName):
    """
    INPUT PARAMETERS: 
        TableName: name of the table to describe
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Print information about table, including 
        parameter names, formats and wavenumber range.
    ---
    EXAMPLE OF USAGE:
        describe('sampletab')
    ---
    """
    describeTable(TableName)


# ---------------------- /ISO.PY ---------------------------------------


def db_begin(db=None):
    """
    INPUT PARAMETERS: 
        db: database name (optional)
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Open a database connection. A database is stored 
        in a folder given in db input parameter.
        Default=data
    ---
    EXAMPLE OF USAGE:
        db_begin('bar')
    ---
    """
    databaseBegin(db)


def db_commit():
    """
    INPUT PARAMETERS: 
        none
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Commit all changes made to opened database.
        All tables will be saved in corresponding files.
    ---
    EXAMPLE OF USAGE:
        db_commit()
    ---
    """
    databaseCommit()


# ------------------ QUERY HITRAN ---------------------------------------


def comment(TableName, Comment):
    LOCAL_TABLE_CACHE[TableName]["header"]["comment"] = Comment


def fetch_by_ids(
    TableName, iso_id_list, numin, numax, ParameterGroups=[], Parameters=[]
):
    """
    INPUT PARAMETERS: 
        TableName:   local table name to fetch in (required)
        iso_id_list: list of isotopologue id's    (required)
        numin:       lower wavenumber bound       (required)
        numax:       upper wavenumber bound       (required)
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Download line-by-line data from HITRANonline server
        and save it to local table. The input parameter iso_id_list
        contains list of "global" isotopologue Ids (see help on ISO_ID).
        Note: this function is required if user wants to download
        multiple species into single table.
    ---
    EXAMPLE OF USAGE:
        fetch_by_ids('water',[1,2,3,4],4000,4100)
    ---
    """
    if type(iso_id_list) not in set([list, tuple]):
        iso_id_list = [iso_id_list]
    queryHITRAN(
        TableName,
        iso_id_list,
        numin,
        numax,
        pargroups=ParameterGroups,
        params=Parameters,
    )
    iso_names = [ISO_ID[i][ISO_ID_INDEX["iso_name"]] for i in iso_id_list]
    Comment = "Contains lines for " + ",".join(iso_names)
    Comment += "\n in %.3f-%.3f wavenumber range" % (numin, numax)
    comment(TableName, Comment)


# def queryHITRAN(TableName,iso_id_list,numin,numax):


def fetch(TableName, M, I, numin, numax, ParameterGroups=[], Parameters=[]):
    """
    INPUT PARAMETERS: 
        TableName:   local table name to fetch in (required)
        M:           HITRAN molecule number       (required)
        I:           HITRAN isotopologue number   (required)
        numin:       lower wavenumber bound       (required)
        numax:       upper wavenumber bound       (required)
    OUTPUT PARAMETERS: 
        none
    ---
    DESCRIPTION:
        Download line-by-line data from HITRANonline server
        and save it to local table. The input parameters M and I
        are the HITRAN molecule and isotopologue numbers.
        This function results in a table containing single 
        isotopologue specie. To have multiple species in a 
        single table use fetch_by_ids instead.
    ---
    EXAMPLE OF USAGE:
        fetch('HOH',1,1,4000,4100)
    ---
    """
    queryHITRAN(
        TableName,
        [ISO[(M, I)][ISO_INDEX["id"]]],
        numin,
        numax,
        pargroups=ParameterGroups,
        params=Parameters,
    )
    iso_name = ISO[(M, I)][ISO_INDEX["iso_name"]]
    Comment = "Contains lines for " + iso_name
    Comment += "\n in %.3f-%.3f wavenumber range" % (numin, numax)
    comment(TableName, Comment)


# ------------------ partition sum --------------------------------------

# ------------------- LAGRANGE INTERPOLATION ----------------------

# def AtoB(aa,bb,A,B,npt)


def AtoB(aa, A, B, npt):
    # ***************************
    # ...LaGrange 3- and 4-point interpolation
    # ...arrays A and B are the npt data points,  given aa, a value of the
    # ...A variable, the routine will find the corresponding bb value
    #
    # ...input:  aa
    # ...output: bb
    for I in range(2, npt + 1):
        if A[I - 1] >= aa:
            if I < 3 or I == npt:
                J = I
                if I < 3:
                    J = 3
                if I == npt:
                    J = npt
                J = J - 1  # zero index correction
                A0D1 = A[J - 2] - A[J - 1]
                if A0D1 == 0.0:
                    A0D1 = 0.0001
                A0D2 = A[J - 2] - A[J]
                if A0D2 == 0.0:
                    A0D2 = 0.0000
                A1D1 = A[J - 1] - A[J - 2]
                if A1D1 == 0.0:
                    A1D1 = 0.0001
                A1D2 = A[J - 1] - A[J]
                if A1D2 == 0.0:
                    A1D2 = 0.0001
                A2D1 = A[J] - A[J - 2]
                if A2D1 == 0.0:
                    A2D1 = 0.0001
                A2D2 = A[J] - A[J - 1]
                if A2D2 == 0.0:
                    A2D2 = 0.0001

                A0 = (aa - A[J - 1]) * (aa - A[J]) / (A0D1 * A0D2)
                A1 = (aa - A[J - 2]) * (aa - A[J]) / (A1D1 * A1D2)
                A2 = (aa - A[J - 2]) * (aa - A[J - 1]) / (A2D1 * A2D2)

                bb = A0 * B[J - 2] + A1 * B[J - 1] + A2 * B[J]

            else:
                J = I
                J = J - 1  # zero index correction
                A0D1 = A[J - 2] - A[J - 1]
                if A0D1 == 0.0:
                    A0D1 = 0.0001
                A0D2 = A[J - 2] - A[J]
                if A0D2 == 0.0:
                    A0D2 = 0.0001
                A0D3 = A[J - 2] - A[J + 1]
                if A0D3 == 0.0:
                    A0D3 = 0.0001
                A1D1 = A[J - 1] - A[J - 2]
                if A1D1 == 0.0:
                    A1D1 = 0.0001
                A1D2 = A[J - 1] - A[J]
                if A1D2 == 0.0:
                    A1D2 = 0.0001
                A1D3 = A[J - 1] - A[J + 1]
                if A1D3 == 0.0:
                    A1D3 = 0.0001

                A2D1 = A[J] - A[J - 2]
                if A2D1 == 0.0:
                    A2D1 = 0.0001
                A2D2 = A[J] - A[J - 1]
                if A2D2 == 0.0:
                    A2D2 = 0.0001
                A2D3 = A[J] - A[J + 1]
                if A2D3 == 0.0:
                    A2D3 = 0.0001

                A3D1 = A[J + 1] - A[J - 2]
                if A3D1 == 0.0:
                    A3D1 = 0.0001
                A3D2 = A[J + 1] - A[J - 1]
                if A3D2 == 0.0:
                    A3D2 = 0.0001
                A3D3 = A[J + 1] - A[J]
                if A3D3 == 0.0:
                    A3D3 = 0.0001

                A0 = (aa - A[J - 1]) * (aa - A[J]) * (aa - A[J + 1])
                A0 = A0 / (A0D1 * A0D2 * A0D3)
                A1 = (aa - A[J - 2]) * (aa - A[J]) * (aa - A[J + 1])
                A1 = A1 / (A1D1 * A1D2 * A1D3)
                A2 = (aa - A[J - 2]) * (aa - A[J - 1]) * (aa - A[J + 1])
                A2 = A2 / (A2D1 * A2D2 * A2D3)
                A3 = (aa - A[J - 2]) * (aa - A[J - 1]) * (aa - A[J])
                A3 = A3 / (A3D1 * A3D2 * A3D3)

                bb = A0 * B[J - 2] + A1 * B[J - 1] + A2 * B[J] + A3 * B[J + 1]

            break

    return bb


#  --------------- ISOTOPOLOGUE HASH ----------------------

TIPS_ISO_HASH = {}

#  --------------- STATISTICAL WEIGHT HASH ----------------------

TIPS_GSI_HASH = {}

#  --------------- INTERPOLATION NODES ----------------------

Tdat = __FloatType__(
    [
        60.0,
        85.0,
        110.0,
        135.0,
        160.0,
        185.0,
        210.0,
        235.0,
        260.0,
        285.0,
        310.0,
        335.0,
        360.0,
        385.0,
        410.0,
        435.0,
        460.0,
        485.0,
        510.0,
        535.0,
        560.0,
        585.0,
        610.0,
        635.0,
        660.0,
        685.0,
        710.0,
        735.0,
        760.0,
        785.0,
        810.0,
        835.0,
        860.0,
        885.0,
        910.0,
        935.0,
        960.0,
        985.0,
        1010.0,
        1035.0,
        1060.0,
        1085.0,
        1110.0,
        1135.0,
        1160.0,
        1185.0,
        1210.0,
        1235.0,
        1260.0,
        1285.0,
        1310.0,
        1335.0,
        1360.0,
        1385.0,
        1410.0,
        1435.0,
        1460.0,
        1485.0,
        1510.0,
        1535.0,
        1560.0,
        1585.0,
        1610.0,
        1635.0,
        1660.0,
        1685.0,
        1710.0,
        1735.0,
        1760.0,
        1785.0,
        1810.0,
        1835.0,
        1860.0,
        1885.0,
        1910.0,
        1935.0,
        1960.0,
        1985.0,
        2010.0,
        2035.0,
        2060.0,
        2085.0,
        2110.0,
        2135.0,
        2160.0,
        2185.0,
        2210.0,
        2235.0,
        2260.0,
        2285.0,
        2310.0,
        2335.0,
        2360.0,
        2385.0,
        2410.0,
        2435.0,
        2460.0,
        2485.0,
        2510.0,
        2535.0,
        2560.0,
        2585.0,
        2610.0,
        2635.0,
        2660.0,
        2685.0,
        2710.0,
        2735.0,
        2760.0,
        2785.0,
        2810.0,
        2835.0,
        2860.0,
        2885.0,
        2910.0,
        2935.0,
        2960.0,
        2985.0,
        3010.0,
    ]
)

TIPS_NPT = len(Tdat)


# REMARK
# float32 gives exactly the same results as fortran TIPS, because
# all constants in the fortran code given as xx.xxE+-XX, i.e.
# in single precision. By this fact all unsignificant figures
# over single precision are filled with digital garbage


#  --------------- H2O 161: M = 1, I = 1 ---------------------
M = 1
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.16824e02,
        0.27771e02,
        0.40408e02,
        0.54549e02,
        0.70054e02,
        0.86817e02,
        0.10475e03,
        0.12380e03,
        0.14391e03,
        0.16503e03,
        0.18714e03,
        0.21021e03,
        0.23425e03,
        0.25924e03,
        0.28518e03,
        0.31209e03,
        0.33997e03,
        0.36883e03,
        0.39870e03,
        0.42959e03,
        0.46152e03,
        0.49452e03,
        0.52860e03,
        0.56380e03,
        0.60015e03,
        0.63766e03,
        0.67637e03,
        0.71631e03,
        0.75750e03,
        0.79999e03,
        0.84380e03,
        0.88897e03,
        0.93553e03,
        0.98353e03,
        0.10330e04,
        0.10840e04,
        0.11365e04,
        0.11906e04,
        0.12463e04,
        0.13037e04,
        0.13628e04,
        0.14237e04,
        0.14863e04,
        0.15509e04,
        0.16173e04,
        0.16856e04,
        0.17559e04,
        0.18283e04,
        0.19028e04,
        0.19793e04,
        0.20581e04,
        0.21391e04,
        0.22224e04,
        0.23080e04,
        0.24067e04,
        0.24975e04,
        0.25908e04,
        0.26867e04,
        0.27853e04,
        0.28865e04,
        0.29904e04,
        0.30972e04,
        0.32068e04,
        0.33194e04,
        0.34349e04,
        0.35535e04,
        0.36752e04,
        0.38001e04,
        0.39282e04,
        0.40597e04,
        0.41945e04,
        0.43327e04,
        0.44745e04,
        0.46199e04,
        0.47688e04,
        0.49215e04,
        0.50780e04,
        0.52384e04,
        0.54027e04,
        0.55710e04,
        0.57434e04,
        0.59200e04,
        0.61008e04,
        0.62859e04,
        0.64754e04,
        0.66693e04,
        0.68679e04,
        0.70710e04,
        0.72788e04,
        0.74915e04,
        0.77090e04,
        0.79315e04,
        0.81590e04,
        0.83917e04,
        0.86296e04,
        0.88728e04,
        0.91214e04,
        0.93755e04,
        0.96351e04,
        0.99005e04,
        0.10171e05,
        0.10448e05,
        0.10731e05,
        0.11020e05,
        0.11315e05,
        0.11617e05,
        0.11924e05,
        0.12238e05,
        0.12559e05,
        0.12886e05,
        0.13220e05,
        0.13561e05,
        0.13909e05,
        0.14263e05,
        0.14625e05,
        0.14995e05,
        0.15371e05,
        0.15755e05,
        0.16147e05,
    ]
)


#  --------------- H2O 181: M = 1, I = 2 ---------------------
M = 1
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.15960e02,
        0.26999e02,
        0.39743e02,
        0.54003e02,
        0.69639e02,
        0.86543e02,
        0.10463e03,
        0.12384e03,
        0.14412e03,
        0.16542e03,
        0.18773e03,
        0.21103e03,
        0.23531e03,
        0.26057e03,
        0.28681e03,
        0.31406e03,
        0.34226e03,
        0.37130e03,
        0.40135e03,
        0.43243e03,
        0.46456e03,
        0.49777e03,
        0.53206e03,
        0.56748e03,
        0.60405e03,
        0.64179e03,
        0.68074e03,
        0.72093e03,
        0.76238e03,
        0.80513e03,
        0.84922e03,
        0.89467e03,
        0.94152e03,
        0.98982e03,
        0.10396e04,
        0.10909e04,
        0.11437e04,
        0.11982e04,
        0.12543e04,
        0.13120e04,
        0.13715e04,
        0.14328e04,
        0.14959e04,
        0.15608e04,
        0.16276e04,
        0.16964e04,
        0.17672e04,
        0.18401e04,
        0.19151e04,
        0.19922e04,
        0.20715e04,
        0.21531e04,
        0.22370e04,
        0.23232e04,
        0.24118e04,
        0.25030e04,
        0.25967e04,
        0.26929e04,
        0.27918e04,
        0.28934e04,
        0.29978e04,
        0.31050e04,
        0.32151e04,
        0.33281e04,
        0.34441e04,
        0.35632e04,
        0.36854e04,
        0.38108e04,
        0.39395e04,
        0.40715e04,
        0.42070e04,
        0.43459e04,
        0.44883e04,
        0.46343e04,
        0.47840e04,
        0.49374e04,
        0.50946e04,
        0.52558e04,
        0.54209e04,
        0.55900e04,
        0.57632e04,
        0.59407e04,
        0.61224e04,
        0.63084e04,
        0.64988e04,
        0.66938e04,
        0.68933e04,
        0.70975e04,
        0.73064e04,
        0.75202e04,
        0.77389e04,
        0.79625e04,
        0.81913e04,
        0.84252e04,
        0.86644e04,
        0.89089e04,
        0.91588e04,
        0.94143e04,
        0.96754e04,
        0.99422e04,
        0.10215e05,
        0.10493e05,
        0.10778e05,
        0.11068e05,
        0.11365e05,
        0.11668e05,
        0.11977e05,
        0.12293e05,
        0.12616e05,
        0.12945e05,
        0.13281e05,
        0.13624e05,
        0.13973e05,
        0.14330e05,
        0.14694e05,
        0.15066e05,
        0.15445e05,
        0.15831e05,
        0.16225e05,
    ]
)


#  --------------- H2O 171: M = 1, I = 3 ---------------------
M = 1
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.95371e02,
        0.16134e03,
        0.23750e03,
        0.32273e03,
        0.41617e03,
        0.51722e03,
        0.62540e03,
        0.74036e03,
        0.86185e03,
        0.98970e03,
        0.11238e04,
        0.12642e04,
        0.14097e04,
        0.15599e04,
        0.17159e04,
        0.18777e04,
        0.20453e04,
        0.22188e04,
        0.23983e04,
        0.25840e04,
        0.27760e04,
        0.29743e04,
        0.31792e04,
        0.33907e04,
        0.36091e04,
        0.38346e04,
        0.40672e04,
        0.43072e04,
        0.45547e04,
        0.48100e04,
        0.50732e04,
        0.53446e04,
        0.56244e04,
        0.59128e04,
        0.62100e04,
        0.65162e04,
        0.68317e04,
        0.71567e04,
        0.74915e04,
        0.78363e04,
        0.81914e04,
        0.85571e04,
        0.89335e04,
        0.93211e04,
        0.97200e04,
        0.10131e05,
        0.10553e05,
        0.10988e05,
        0.11435e05,
        0.11895e05,
        0.12368e05,
        0.12855e05,
        0.13356e05,
        0.13870e05,
        0.14399e05,
        0.14943e05,
        0.15502e05,
        0.16076e05,
        0.16666e05,
        0.17272e05,
        0.17895e05,
        0.18534e05,
        0.19191e05,
        0.19865e05,
        0.20557e05,
        0.21267e05,
        0.21996e05,
        0.22744e05,
        0.23512e05,
        0.24299e05,
        0.25106e05,
        0.25935e05,
        0.26784e05,
        0.27655e05,
        0.28547e05,
        0.29462e05,
        0.30400e05,
        0.31361e05,
        0.32345e05,
        0.33353e05,
        0.34386e05,
        0.35444e05,
        0.36527e05,
        0.37637e05,
        0.38772e05,
        0.39934e05,
        0.41124e05,
        0.42341e05,
        0.43587e05,
        0.44861e05,
        0.46165e05,
        0.47498e05,
        0.48862e05,
        0.50256e05,
        0.51682e05,
        0.53139e05,
        0.54629e05,
        0.56152e05,
        0.57708e05,
        0.59299e05,
        0.60923e05,
        0.62583e05,
        0.64279e05,
        0.66011e05,
        0.67779e05,
        0.69585e05,
        0.71429e05,
        0.73312e05,
        0.75234e05,
        0.77195e05,
        0.79197e05,
        0.81240e05,
        0.83325e05,
        0.85452e05,
        0.87622e05,
        0.89835e05,
        0.92093e05,
        0.94395e05,
        0.96743e05,
    ]
)


#  --------------- H2O 162: M = 1, I = 4 ---------------------
M = 1
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.75792e02,
        0.12986e03,
        0.19244e03,
        0.26253e03,
        0.33942e03,
        0.42259e03,
        0.51161e03,
        0.60619e03,
        0.70609e03,
        0.81117e03,
        0.92132e03,
        0.10365e04,
        0.11567e04,
        0.12820e04,
        0.14124e04,
        0.15481e04,
        0.16891e04,
        0.18355e04,
        0.19876e04,
        0.21455e04,
        0.23092e04,
        0.24791e04,
        0.26551e04,
        0.28376e04,
        0.30268e04,
        0.32258e04,
        0.34288e04,
        0.36392e04,
        0.38571e04,
        0.40828e04,
        0.43165e04,
        0.45584e04,
        0.48089e04,
        0.50681e04,
        0.53363e04,
        0.56139e04,
        0.59009e04,
        0.61979e04,
        0.65049e04,
        0.68224e04,
        0.71506e04,
        0.74898e04,
        0.78403e04,
        0.82024e04,
        0.85765e04,
        0.89628e04,
        0.93618e04,
        0.97736e04,
        0.10199e05,
        0.10637e05,
        0.11090e05,
        0.11557e05,
        0.12039e05,
        0.12535e05,
        0.13047e05,
        0.13575e05,
        0.14119e05,
        0.14679e05,
        0.15257e05,
        0.15851e05,
        0.16464e05,
        0.17094e05,
        0.17743e05,
        0.18411e05,
        0.19098e05,
        0.19805e05,
        0.20532e05,
        0.21280e05,
        0.22049e05,
        0.22840e05,
        0.23652e05,
        0.24487e05,
        0.25345e05,
        0.26227e05,
        0.27132e05,
        0.28062e05,
        0.29016e05,
        0.29997e05,
        0.31002e05,
        0.32035e05,
        0.33094e05,
        0.34180e05,
        0.35295e05,
        0.36438e05,
        0.37610e05,
        0.38812e05,
        0.40044e05,
        0.41306e05,
        0.42600e05,
        0.43926e05,
        0.45284e05,
        0.46675e05,
        0.48100e05,
        0.49559e05,
        0.51053e05,
        0.52583e05,
        0.54148e05,
        0.55750e05,
        0.57390e05,
        0.59067e05,
        0.60783e05,
        0.62539e05,
        0.64334e05,
        0.66170e05,
        0.68047e05,
        0.69967e05,
        0.71929e05,
        0.73934e05,
        0.75983e05,
        0.78078e05,
        0.80217e05,
        0.82403e05,
        0.84636e05,
        0.86917e05,
        0.89246e05,
        0.91625e05,
        0.94053e05,
        0.96533e05,
        0.99064e05,
    ]
)


#  --------------- H2O 182: M = 1, I = 5 ---------------------
M = 1
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.82770e02,
        0.13749e03,
        0.20083e03,
        0.27176e03,
        0.34955e03,
        0.43370e03,
        0.52376e03,
        0.61944e03,
        0.72050e03,
        0.82679e03,
        0.93821e03,
        0.10547e04,
        0.11763e04,
        0.13031e04,
        0.14350e04,
        0.15723e04,
        0.17150e04,
        0.18633e04,
        0.20172e04,
        0.21770e04,
        0.23429e04,
        0.25149e04,
        0.26934e04,
        0.28784e04,
        0.30702e04,
        0.32690e04,
        0.34750e04,
        0.36885e04,
        0.39096e04,
        0.41386e04,
        0.43758e04,
        0.46213e04,
        0.48755e04,
        0.51386e04,
        0.54109e04,
        0.56927e04,
        0.59841e04,
        0.62856e04,
        0.65973e04,
        0.69197e04,
        0.72529e04,
        0.75973e04,
        0.79533e04,
        0.83210e04,
        0.87009e04,
        0.90933e04,
        0.94985e04,
        0.99168e04,
        0.10348e05,
        0.10794e05,
        0.11254e05,
        0.11728e05,
        0.12217e05,
        0.12722e05,
        0.13242e05,
        0.13778e05,
        0.14331e05,
        0.14900e05,
        0.15486e05,
        0.16091e05,
        0.16713e05,
        0.17353e05,
        0.18012e05,
        0.18691e05,
        0.19389e05,
        0.20108e05,
        0.20847e05,
        0.21607e05,
        0.22388e05,
        0.23191e05,
        0.24017e05,
        0.24866e05,
        0.25738e05,
        0.26633e05,
        0.27553e05,
        0.28498e05,
        0.29468e05,
        0.30464e05,
        0.31486e05,
        0.32536e05,
        0.33612e05,
        0.34716e05,
        0.35849e05,
        0.37011e05,
        0.38202e05,
        0.39424e05,
        0.40676e05,
        0.41959e05,
        0.43274e05,
        0.44622e05,
        0.46002e05,
        0.47416e05,
        0.48864e05,
        0.50348e05,
        0.51866e05,
        0.53421e05,
        0.55012e05,
        0.56640e05,
        0.58307e05,
        0.60012e05,
        0.61757e05,
        0.63541e05,
        0.65366e05,
        0.67233e05,
        0.69141e05,
        0.71092e05,
        0.73087e05,
        0.75125e05,
        0.77209e05,
        0.79338e05,
        0.81513e05,
        0.83736e05,
        0.86006e05,
        0.88324e05,
        0.90693e05,
        0.93111e05,
        0.95580e05,
        0.98100e05,
        0.10067e06,
    ]
)


#  --------------- H2O 172: M = 1, I = 6 ---------------------
M = 1
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(36.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.49379e03,
        0.82021e03,
        0.11980e04,
        0.16211e04,
        0.20851e04,
        0.25870e04,
        0.31242e04,
        0.36949e04,
        0.42977e04,
        0.49317e04,
        0.55963e04,
        0.62911e04,
        0.70164e04,
        0.77722e04,
        0.85591e04,
        0.93777e04,
        0.10228e05,
        0.11112e05,
        0.12030e05,
        0.12983e05,
        0.13971e05,
        0.14997e05,
        0.16061e05,
        0.17163e05,
        0.18306e05,
        0.19491e05,
        0.20719e05,
        0.21991e05,
        0.23309e05,
        0.24673e05,
        0.26086e05,
        0.27549e05,
        0.29064e05,
        0.30631e05,
        0.32254e05,
        0.33932e05,
        0.35669e05,
        0.37464e05,
        0.39321e05,
        0.41242e05,
        0.43227e05,
        0.45279e05,
        0.47399e05,
        0.49589e05,
        0.51852e05,
        0.54189e05,
        0.56602e05,
        0.59094e05,
        0.61666e05,
        0.64320e05,
        0.67058e05,
        0.69883e05,
        0.72796e05,
        0.75801e05,
        0.78899e05,
        0.82092e05,
        0.85382e05,
        0.88773e05,
        0.92266e05,
        0.95863e05,
        0.99568e05,
        0.10338e06,
        0.10731e06,
        0.11135e06,
        0.11551e06,
        0.11979e06,
        0.12419e06,
        0.12871e06,
        0.13337e06,
        0.13815e06,
        0.14307e06,
        0.14812e06,
        0.15331e06,
        0.15865e06,
        0.16412e06,
        0.16975e06,
        0.17553e06,
        0.18146e06,
        0.18754e06,
        0.19379e06,
        0.20020e06,
        0.20678e06,
        0.21352e06,
        0.22044e06,
        0.22753e06,
        0.23480e06,
        0.24226e06,
        0.24990e06,
        0.25773e06,
        0.26575e06,
        0.27397e06,
        0.28239e06,
        0.29102e06,
        0.29985e06,
        0.30889e06,
        0.31814e06,
        0.32762e06,
        0.33731e06,
        0.34724e06,
        0.35739e06,
        0.36777e06,
        0.37840e06,
        0.38926e06,
        0.40038e06,
        0.41174e06,
        0.42335e06,
        0.43523e06,
        0.44737e06,
        0.45977e06,
        0.47245e06,
        0.48540e06,
        0.49863e06,
        0.51214e06,
        0.52595e06,
        0.54005e06,
        0.55444e06,
        0.56914e06,
        0.58415e06,
        0.59947e06,
    ]
)


#  --------------- CO2 626: M = 2, I = 1 ---------------------
M = 2
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.53642e02,
        0.75947e02,
        0.98292e02,
        0.12078e03,
        0.14364e03,
        0.16714e03,
        0.19160e03,
        0.21731e03,
        0.24454e03,
        0.27355e03,
        0.30456e03,
        0.33778e03,
        0.37343e03,
        0.41170e03,
        0.45280e03,
        0.49692e03,
        0.54427e03,
        0.59505e03,
        0.64948e03,
        0.70779e03,
        0.77019e03,
        0.83693e03,
        0.90825e03,
        0.98440e03,
        0.10656e04,
        0.11522e04,
        0.12445e04,
        0.13427e04,
        0.14471e04,
        0.15580e04,
        0.16759e04,
        0.18009e04,
        0.19334e04,
        0.20739e04,
        0.22225e04,
        0.23798e04,
        0.25462e04,
        0.27219e04,
        0.29074e04,
        0.31032e04,
        0.33097e04,
        0.35272e04,
        0.37564e04,
        0.39976e04,
        0.42514e04,
        0.45181e04,
        0.47985e04,
        0.50929e04,
        0.54019e04,
        0.57260e04,
        0.60659e04,
        0.64221e04,
        0.67952e04,
        0.71859e04,
        0.75946e04,
        0.80222e04,
        0.84691e04,
        0.89362e04,
        0.94241e04,
        0.99335e04,
        0.10465e05,
        0.11020e05,
        0.11598e05,
        0.12201e05,
        0.12828e05,
        0.13482e05,
        0.14163e05,
        0.14872e05,
        0.15609e05,
        0.16376e05,
        0.17173e05,
        0.18001e05,
        0.18861e05,
        0.19754e05,
        0.20682e05,
        0.21644e05,
        0.22643e05,
        0.23678e05,
        0.24752e05,
        0.25865e05,
        0.27018e05,
        0.28212e05,
        0.29449e05,
        0.30730e05,
        0.32055e05,
        0.33426e05,
        0.34845e05,
        0.36312e05,
        0.37828e05,
        0.39395e05,
        0.41015e05,
        0.42688e05,
        0.44416e05,
        0.46199e05,
        0.48041e05,
        0.49942e05,
        0.51902e05,
        0.53925e05,
        0.56011e05,
        0.58162e05,
        0.60379e05,
        0.62664e05,
        0.65019e05,
        0.67444e05,
        0.69942e05,
        0.72515e05,
        0.75163e05,
        0.77890e05,
        0.80695e05,
        0.83582e05,
        0.86551e05,
        0.89605e05,
        0.92746e05,
        0.95975e05,
        0.99294e05,
        0.10271e06,
        0.10621e06,
        0.10981e06,
        0.11351e06,
    ]
)


#  --------------- CO2 636: M = 2, I = 2 ---------------------
M = 2
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.10728e03,
        0.15189e03,
        0.19659e03,
        0.24164e03,
        0.28753e03,
        0.33486e03,
        0.38429e03,
        0.43643e03,
        0.49184e03,
        0.55104e03,
        0.61449e03,
        0.68263e03,
        0.75589e03,
        0.83468e03,
        0.91943e03,
        0.10106e04,
        0.11085e04,
        0.12137e04,
        0.13266e04,
        0.14477e04,
        0.15774e04,
        0.17163e04,
        0.18649e04,
        0.20237e04,
        0.21933e04,
        0.23743e04,
        0.25673e04,
        0.27729e04,
        0.29917e04,
        0.32245e04,
        0.34718e04,
        0.37345e04,
        0.40132e04,
        0.43087e04,
        0.46218e04,
        0.49533e04,
        0.53041e04,
        0.56749e04,
        0.60668e04,
        0.64805e04,
        0.69171e04,
        0.73774e04,
        0.78626e04,
        0.83736e04,
        0.89114e04,
        0.94772e04,
        0.10072e05,
        0.10697e05,
        0.11353e05,
        0.12042e05,
        0.12765e05,
        0.13523e05,
        0.14317e05,
        0.15148e05,
        0.16019e05,
        0.16930e05,
        0.17883e05,
        0.18879e05,
        0.19920e05,
        0.21008e05,
        0.22143e05,
        0.23328e05,
        0.24563e05,
        0.25852e05,
        0.27195e05,
        0.28594e05,
        0.30051e05,
        0.31568e05,
        0.33146e05,
        0.34788e05,
        0.36496e05,
        0.38271e05,
        0.40115e05,
        0.42031e05,
        0.44021e05,
        0.46086e05,
        0.48230e05,
        0.50453e05,
        0.52759e05,
        0.55150e05,
        0.57628e05,
        0.60195e05,
        0.62854e05,
        0.65608e05,
        0.68459e05,
        0.71409e05,
        0.74461e05,
        0.77618e05,
        0.80883e05,
        0.84258e05,
        0.87746e05,
        0.91350e05,
        0.95073e05,
        0.98918e05,
        0.10289e06,
        0.10698e06,
        0.11121e06,
        0.11558e06,
        0.12008e06,
        0.12472e06,
        0.12950e06,
        0.13443e06,
        0.13952e06,
        0.14475e06,
        0.15015e06,
        0.15571e06,
        0.16143e06,
        0.16732e06,
        0.17338e06,
        0.17962e06,
        0.18604e06,
        0.19264e06,
        0.19943e06,
        0.20642e06,
        0.21360e06,
        0.22098e06,
        0.22856e06,
        0.23636e06,
        0.24436e06,
    ]
)


#  --------------- CO2 628: M = 2, I = 3 ---------------------
M = 2
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.11368e03,
        0.16096e03,
        0.20833e03,
        0.25603e03,
        0.30452e03,
        0.35442e03,
        0.40640e03,
        0.46110e03,
        0.51910e03,
        0.58093e03,
        0.64709e03,
        0.71804e03,
        0.79422e03,
        0.87607e03,
        0.96402e03,
        0.10585e04,
        0.11600e04,
        0.12689e04,
        0.13857e04,
        0.15108e04,
        0.16449e04,
        0.17883e04,
        0.19416e04,
        0.21054e04,
        0.22803e04,
        0.24668e04,
        0.26655e04,
        0.28770e04,
        0.31021e04,
        0.33414e04,
        0.35956e04,
        0.38654e04,
        0.41516e04,
        0.44549e04,
        0.47761e04,
        0.51160e04,
        0.54755e04,
        0.58555e04,
        0.62568e04,
        0.66804e04,
        0.71273e04,
        0.75982e04,
        0.80944e04,
        0.86169e04,
        0.91666e04,
        0.97446e04,
        0.10352e05,
        0.10990e05,
        0.11660e05,
        0.12363e05,
        0.13101e05,
        0.13874e05,
        0.14683e05,
        0.15531e05,
        0.16418e05,
        0.17347e05,
        0.18317e05,
        0.19332e05,
        0.20392e05,
        0.21499e05,
        0.22654e05,
        0.23859e05,
        0.25116e05,
        0.26426e05,
        0.27792e05,
        0.29214e05,
        0.30695e05,
        0.32236e05,
        0.33840e05,
        0.35508e05,
        0.37242e05,
        0.39045e05,
        0.40917e05,
        0.42862e05,
        0.44881e05,
        0.46977e05,
        0.49152e05,
        0.51407e05,
        0.53746e05,
        0.56171e05,
        0.58683e05,
        0.61286e05,
        0.63981e05,
        0.66772e05,
        0.69661e05,
        0.72650e05,
        0.75742e05,
        0.78940e05,
        0.82246e05,
        0.85664e05,
        0.89196e05,
        0.92845e05,
        0.96613e05,
        0.10050e06,
        0.10452e06,
        0.10867e06,
        0.11295e06,
        0.11736e06,
        0.12191e06,
        0.12661e06,
        0.13145e06,
        0.13643e06,
        0.14157e06,
        0.14687e06,
        0.15232e06,
        0.15794e06,
        0.16372e06,
        0.16968e06,
        0.17580e06,
        0.18211e06,
        0.18859e06,
        0.19526e06,
        0.20213e06,
        0.20918e06,
        0.21643e06,
        0.22388e06,
        0.23154e06,
        0.23941e06,
        0.24750e06,
    ]
)


#  --------------- CO2 627: M = 2, I = 4 ---------------------
M = 2
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.66338e03,
        0.93923e03,
        0.12156e04,
        0.14938e04,
        0.17766e04,
        0.20676e04,
        0.23705e04,
        0.26891e04,
        0.30267e04,
        0.33866e04,
        0.37714e04,
        0.41839e04,
        0.46267e04,
        0.51023e04,
        0.56132e04,
        0.61618e04,
        0.67508e04,
        0.73827e04,
        0.80603e04,
        0.87863e04,
        0.95636e04,
        0.10395e05,
        0.11284e05,
        0.12233e05,
        0.13246e05,
        0.14326e05,
        0.15477e05,
        0.16702e05,
        0.18005e05,
        0.19390e05,
        0.20861e05,
        0.22422e05,
        0.24077e05,
        0.25832e05,
        0.27689e05,
        0.29655e05,
        0.31734e05,
        0.33931e05,
        0.36250e05,
        0.38698e05,
        0.41280e05,
        0.44002e05,
        0.46869e05,
        0.49886e05,
        0.53062e05,
        0.56400e05,
        0.59909e05,
        0.63594e05,
        0.67462e05,
        0.71521e05,
        0.75777e05,
        0.80238e05,
        0.84911e05,
        0.89804e05,
        0.94925e05,
        0.10028e06,
        0.10588e06,
        0.11173e06,
        0.11785e06,
        0.12423e06,
        0.13090e06,
        0.13785e06,
        0.14510e06,
        0.15265e06,
        0.16053e06,
        0.16873e06,
        0.17727e06,
        0.18615e06,
        0.19540e06,
        0.20501e06,
        0.21501e06,
        0.22540e06,
        0.23619e06,
        0.24740e06,
        0.25904e06,
        0.27112e06,
        0.28365e06,
        0.29664e06,
        0.31012e06,
        0.32409e06,
        0.33856e06,
        0.35356e06,
        0.36908e06,
        0.38516e06,
        0.40180e06,
        0.41902e06,
        0.43683e06,
        0.45525e06,
        0.47429e06,
        0.49397e06,
        0.51431e06,
        0.53532e06,
        0.55702e06,
        0.57943e06,
        0.60256e06,
        0.62644e06,
        0.65107e06,
        0.67648e06,
        0.70269e06,
        0.72972e06,
        0.75758e06,
        0.78629e06,
        0.81588e06,
        0.84636e06,
        0.87775e06,
        0.91008e06,
        0.94337e06,
        0.97763e06,
        0.10129e07,
        0.10492e07,
        0.10865e07,
        0.11249e07,
        0.11644e07,
        0.12050e07,
        0.12467e07,
        0.12896e07,
        0.13337e07,
        0.13789e07,
        0.14255e07,
    ]
)


#  --------------- CO2 638: M = 2, I = 5 ---------------------
M = 2
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.22737e03,
        0.32194e03,
        0.41671e03,
        0.51226e03,
        0.60963e03,
        0.71017e03,
        0.81528e03,
        0.92628e03,
        0.10444e04,
        0.11707e04,
        0.13061e04,
        0.14518e04,
        0.16085e04,
        0.17772e04,
        0.19588e04,
        0.21542e04,
        0.23644e04,
        0.25903e04,
        0.28330e04,
        0.30934e04,
        0.33726e04,
        0.36717e04,
        0.39918e04,
        0.43342e04,
        0.47001e04,
        0.50907e04,
        0.55074e04,
        0.59515e04,
        0.64244e04,
        0.69276e04,
        0.74626e04,
        0.80310e04,
        0.86344e04,
        0.92744e04,
        0.99528e04,
        0.10671e05,
        0.11432e05,
        0.12236e05,
        0.13086e05,
        0.13984e05,
        0.14932e05,
        0.15932e05,
        0.16985e05,
        0.18096e05,
        0.19265e05,
        0.20495e05,
        0.21788e05,
        0.23148e05,
        0.24576e05,
        0.26075e05,
        0.27648e05,
        0.29298e05,
        0.31027e05,
        0.32839e05,
        0.34736e05,
        0.36721e05,
        0.38798e05,
        0.40970e05,
        0.43240e05,
        0.45611e05,
        0.48087e05,
        0.50671e05,
        0.53368e05,
        0.56180e05,
        0.59111e05,
        0.62165e05,
        0.65347e05,
        0.68659e05,
        0.72107e05,
        0.75694e05,
        0.79425e05,
        0.83303e05,
        0.87334e05,
        0.91522e05,
        0.95872e05,
        0.10039e06,
        0.10507e06,
        0.10994e06,
        0.11498e06,
        0.12021e06,
        0.12563e06,
        0.13125e06,
        0.13707e06,
        0.14309e06,
        0.14933e06,
        0.15579e06,
        0.16247e06,
        0.16938e06,
        0.17653e06,
        0.18392e06,
        0.19156e06,
        0.19946e06,
        0.20761e06,
        0.21604e06,
        0.22473e06,
        0.23371e06,
        0.24298e06,
        0.25254e06,
        0.26240e06,
        0.27258e06,
        0.28307e06,
        0.29388e06,
        0.30502e06,
        0.31651e06,
        0.32834e06,
        0.34052e06,
        0.35307e06,
        0.36599e06,
        0.37929e06,
        0.39298e06,
        0.40706e06,
        0.42155e06,
        0.43645e06,
        0.45178e06,
        0.46753e06,
        0.48373e06,
        0.50038e06,
        0.51748e06,
        0.53506e06,
    ]
)


#  --------------- CO2 637: M = 2, I = 6 ---------------------
M = 2
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.13267e04,
        0.18785e04,
        0.24314e04,
        0.29888e04,
        0.35566e04,
        0.41426e04,
        0.47550e04,
        0.54013e04,
        0.60886e04,
        0.68232e04,
        0.76109e04,
        0.84574e04,
        0.93678e04,
        0.10348e05,
        0.11402e05,
        0.12536e05,
        0.13755e05,
        0.15065e05,
        0.16471e05,
        0.17980e05,
        0.19598e05,
        0.21330e05,
        0.23184e05,
        0.25166e05,
        0.27283e05,
        0.29543e05,
        0.31953e05,
        0.34521e05,
        0.37256e05,
        0.40164e05,
        0.43256e05,
        0.46541e05,
        0.50026e05,
        0.53723e05,
        0.57641e05,
        0.61790e05,
        0.66180e05,
        0.70823e05,
        0.75729e05,
        0.80910e05,
        0.86378e05,
        0.92145e05,
        0.98224e05,
        0.10463e06,
        0.11137e06,
        0.11846e06,
        0.12592e06,
        0.13375e06,
        0.14198e06,
        0.15062e06,
        0.15969e06,
        0.16920e06,
        0.17916e06,
        0.18959e06,
        0.20052e06,
        0.21196e06,
        0.22392e06,
        0.23642e06,
        0.24949e06,
        0.26314e06,
        0.27740e06,
        0.29227e06,
        0.30779e06,
        0.32398e06,
        0.34085e06,
        0.35842e06,
        0.37673e06,
        0.39579e06,
        0.41563e06,
        0.43626e06,
        0.45772e06,
        0.48003e06,
        0.50322e06,
        0.52730e06,
        0.55232e06,
        0.57829e06,
        0.60524e06,
        0.63320e06,
        0.66219e06,
        0.69226e06,
        0.72342e06,
        0.75571e06,
        0.78916e06,
        0.82380e06,
        0.85966e06,
        0.89678e06,
        0.93518e06,
        0.97490e06,
        0.10160e07,
        0.10585e07,
        0.11023e07,
        0.11477e07,
        0.11946e07,
        0.12430e07,
        0.12929e07,
        0.13445e07,
        0.13977e07,
        0.14526e07,
        0.15093e07,
        0.15677e07,
        0.16280e07,
        0.16901e07,
        0.17541e07,
        0.18200e07,
        0.18880e07,
        0.19579e07,
        0.20300e07,
        0.21042e07,
        0.21805e07,
        0.22591e07,
        0.23400e07,
        0.24232e07,
        0.25087e07,
        0.25967e07,
        0.26871e07,
        0.27801e07,
        0.28757e07,
        0.29739e07,
        0.30747e07,
    ]
)


#  --------------- CO2 828: M = 2, I = 7 ---------------------
M = 2
I = 7
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.60334e02,
        0.85430e02,
        0.11058e03,
        0.13590e03,
        0.16167e03,
        0.18821e03,
        0.21588e03,
        0.24502e03,
        0.27595e03,
        0.30896e03,
        0.34431e03,
        0.38225e03,
        0.42301e03,
        0.46684e03,
        0.51397e03,
        0.56464e03,
        0.61907e03,
        0.67753e03,
        0.74027e03,
        0.80753e03,
        0.87961e03,
        0.95676e03,
        0.10393e04,
        0.11275e04,
        0.12217e04,
        0.13222e04,
        0.14293e04,
        0.15434e04,
        0.16648e04,
        0.17940e04,
        0.19312e04,
        0.20769e04,
        0.22315e04,
        0.23954e04,
        0.25691e04,
        0.27529e04,
        0.29474e04,
        0.31530e04,
        0.33702e04,
        0.35995e04,
        0.38414e04,
        0.40965e04,
        0.43654e04,
        0.46484e04,
        0.49464e04,
        0.52598e04,
        0.55892e04,
        0.59353e04,
        0.62988e04,
        0.66803e04,
        0.70804e04,
        0.74998e04,
        0.79394e04,
        0.83998e04,
        0.88817e04,
        0.93859e04,
        0.99132e04,
        0.10464e05,
        0.11040e05,
        0.11642e05,
        0.12270e05,
        0.12925e05,
        0.13609e05,
        0.14321e05,
        0.15064e05,
        0.15838e05,
        0.16643e05,
        0.17482e05,
        0.18355e05,
        0.19263e05,
        0.20207e05,
        0.21188e05,
        0.22208e05,
        0.23267e05,
        0.24366e05,
        0.25508e05,
        0.26692e05,
        0.27921e05,
        0.29195e05,
        0.30516e05,
        0.31886e05,
        0.33304e05,
        0.34773e05,
        0.36294e05,
        0.37869e05,
        0.39499e05,
        0.41185e05,
        0.42929e05,
        0.44732e05,
        0.46596e05,
        0.48522e05,
        0.50513e05,
        0.52569e05,
        0.54692e05,
        0.56884e05,
        0.59146e05,
        0.61481e05,
        0.63890e05,
        0.66375e05,
        0.68937e05,
        0.71578e05,
        0.74301e05,
        0.77107e05,
        0.79998e05,
        0.82976e05,
        0.86043e05,
        0.89201e05,
        0.92452e05,
        0.95799e05,
        0.99242e05,
        0.10278e06,
        0.10643e06,
        0.11018e06,
        0.11403e06,
        0.11799e06,
        0.12206e06,
        0.12625e06,
        0.13055e06,
        0.13497e06,
    ]
)


#  --------------- CO2 728: M = 2, I = 8 ---------------------
M = 2
I = 8
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.70354e03,
        0.99615e03,
        0.12893e04,
        0.15846e04,
        0.18848e04,
        0.21940e04,
        0.25162e04,
        0.28554e04,
        0.32152e04,
        0.35991e04,
        0.40099e04,
        0.44507e04,
        0.49242e04,
        0.54332e04,
        0.59802e04,
        0.65681e04,
        0.71996e04,
        0.78776e04,
        0.86050e04,
        0.93847e04,
        0.10220e05,
        0.11114e05,
        0.12070e05,
        0.13091e05,
        0.14182e05,
        0.15345e05,
        0.16585e05,
        0.17906e05,
        0.19311e05,
        0.20805e05,
        0.22393e05,
        0.24078e05,
        0.25865e05,
        0.27760e05,
        0.29768e05,
        0.31893e05,
        0.34140e05,
        0.36516e05,
        0.39025e05,
        0.41674e05,
        0.44469e05,
        0.47416e05,
        0.50520e05,
        0.53789e05,
        0.57229e05,
        0.60847e05,
        0.64650e05,
        0.68645e05,
        0.72840e05,
        0.77242e05,
        0.81859e05,
        0.86699e05,
        0.91770e05,
        0.97081e05,
        0.10264e06,
        0.10846e06,
        0.11454e06,
        0.12090e06,
        0.12754e06,
        0.13447e06,
        0.14171e06,
        0.14927e06,
        0.15715e06,
        0.16536e06,
        0.17392e06,
        0.18284e06,
        0.19213e06,
        0.20179e06,
        0.21185e06,
        0.22231e06,
        0.23319e06,
        0.24450e06,
        0.25625e06,
        0.26845e06,
        0.28112e06,
        0.29427e06,
        0.30791e06,
        0.32206e06,
        0.33674e06,
        0.35196e06,
        0.36772e06,
        0.38406e06,
        0.40098e06,
        0.41850e06,
        0.43663e06,
        0.45539e06,
        0.47480e06,
        0.49488e06,
        0.51564e06,
        0.53710e06,
        0.55928e06,
        0.58219e06,
        0.60586e06,
        0.63029e06,
        0.65553e06,
        0.68157e06,
        0.70844e06,
        0.73616e06,
        0.76476e06,
        0.79424e06,
        0.82464e06,
        0.85597e06,
        0.88826e06,
        0.92153e06,
        0.95580e06,
        0.99108e06,
        0.10274e07,
        0.10648e07,
        0.11033e07,
        0.11429e07,
        0.11837e07,
        0.12256e07,
        0.12687e07,
        0.13131e07,
        0.13586e07,
        0.14055e07,
        0.14536e07,
        0.15031e07,
        0.15539e07,
    ]
)


#  --------------- CO2 727: M = 2, I = 9 ---------------------
M = 2
I = 9
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.20518e04,
        0.29051e04,
        0.37601e04,
        0.46209e04,
        0.54961e04,
        0.63969e04,
        0.73353e04,
        0.83227e04,
        0.93698e04,
        0.10486e05,
        0.11681e05,
        0.12962e05,
        0.14337e05,
        0.15815e05,
        0.17403e05,
        0.19110e05,
        0.20942e05,
        0.22909e05,
        0.25018e05,
        0.27278e05,
        0.29699e05,
        0.32290e05,
        0.35060e05,
        0.38019e05,
        0.41177e05,
        0.44545e05,
        0.48135e05,
        0.51957e05,
        0.56023e05,
        0.60346e05,
        0.64938e05,
        0.69812e05,
        0.74981e05,
        0.80461e05,
        0.86264e05,
        0.92406e05,
        0.98902e05,
        0.10577e06,
        0.11302e06,
        0.12067e06,
        0.12875e06,
        0.13726e06,
        0.14622e06,
        0.15566e06,
        0.16559e06,
        0.17604e06,
        0.18702e06,
        0.19855e06,
        0.21066e06,
        0.22336e06,
        0.23669e06,
        0.25065e06,
        0.26528e06,
        0.28061e06,
        0.29664e06,
        0.31342e06,
        0.33096e06,
        0.34930e06,
        0.36845e06,
        0.38845e06,
        0.40933e06,
        0.43111e06,
        0.45383e06,
        0.47751e06,
        0.50219e06,
        0.52790e06,
        0.55466e06,
        0.58252e06,
        0.61151e06,
        0.64166e06,
        0.67300e06,
        0.70558e06,
        0.73943e06,
        0.77458e06,
        0.81108e06,
        0.84896e06,
        0.88827e06,
        0.92904e06,
        0.97131e06,
        0.10151e07,
        0.10605e07,
        0.11076e07,
        0.11563e07,
        0.12068e07,
        0.12590e07,
        0.13130e07,
        0.13689e07,
        0.14267e07,
        0.14865e07,
        0.15483e07,
        0.16121e07,
        0.16781e07,
        0.17462e07,
        0.18165e07,
        0.18892e07,
        0.19641e07,
        0.20415e07,
        0.21213e07,
        0.22036e07,
        0.22884e07,
        0.23759e07,
        0.24661e07,
        0.25590e07,
        0.26547e07,
        0.27533e07,
        0.28549e07,
        0.29594e07,
        0.30670e07,
        0.31778e07,
        0.32918e07,
        0.34090e07,
        0.35296e07,
        0.36536e07,
        0.37812e07,
        0.39123e07,
        0.40470e07,
        0.41855e07,
        0.43278e07,
        0.44739e07,
    ]
)


#  --------------- CO2 838: M = 2, I = 10 ---------------------
M = 2
I = 10
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.12066e03,
        0.17085e03,
        0.22116e03,
        0.27190e03,
        0.32364e03,
        0.37711e03,
        0.43305e03,
        0.49219e03,
        0.55516e03,
        0.62256e03,
        0.69492e03,
        0.77276e03,
        0.85657e03,
        0.94685e03,
        0.10441e04,
        0.11488e04,
        0.12614e04,
        0.13826e04,
        0.15127e04,
        0.16525e04,
        0.18024e04,
        0.19630e04,
        0.21351e04,
        0.23191e04,
        0.25158e04,
        0.27260e04,
        0.29502e04,
        0.31892e04,
        0.34438e04,
        0.37148e04,
        0.40031e04,
        0.43094e04,
        0.46346e04,
        0.49797e04,
        0.53455e04,
        0.57331e04,
        0.61434e04,
        0.65775e04,
        0.70364e04,
        0.75212e04,
        0.80330e04,
        0.85730e04,
        0.91424e04,
        0.97423e04,
        0.10374e05,
        0.11039e05,
        0.11738e05,
        0.12474e05,
        0.13246e05,
        0.14057e05,
        0.14908e05,
        0.15801e05,
        0.16737e05,
        0.17717e05,
        0.18744e05,
        0.19819e05,
        0.20944e05,
        0.22120e05,
        0.23349e05,
        0.24634e05,
        0.25975e05,
        0.27376e05,
        0.28837e05,
        0.30361e05,
        0.31950e05,
        0.33605e05,
        0.35330e05,
        0.37126e05,
        0.38996e05,
        0.40942e05,
        0.42965e05,
        0.45069e05,
        0.47256e05,
        0.49528e05,
        0.51888e05,
        0.54338e05,
        0.56882e05,
        0.59521e05,
        0.62259e05,
        0.65097e05,
        0.68040e05,
        0.71090e05,
        0.74249e05,
        0.77522e05,
        0.80910e05,
        0.84417e05,
        0.88046e05,
        0.91801e05,
        0.95684e05,
        0.99699e05,
        0.10385e06,
        0.10814e06,
        0.11257e06,
        0.11715e06,
        0.12187e06,
        0.12675e06,
        0.13179e06,
        0.13699e06,
        0.14235e06,
        0.14788e06,
        0.15358e06,
        0.15946e06,
        0.16552e06,
        0.17176e06,
        0.17819e06,
        0.18482e06,
        0.19164e06,
        0.19867e06,
        0.20590e06,
        0.21335e06,
        0.22101e06,
        0.22889e06,
        0.23699e06,
        0.24533e06,
        0.25390e06,
        0.26271e06,
        0.27177e06,
        0.28108e06,
        0.29064e06,
    ]
)

#  --------------- CO2 838: M = 2, I = 0 ALIAS-----------------
TIPS_GSI_HASH[(M, 0)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, 0)] = TIPS_ISO_HASH[(M, I)]

#  --------------- CO2 837: M = 2, I = 11 ---------------------
M = 2
I = 11
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.14071e04,
        0.19923e04,
        0.25789e04,
        0.31704e04,
        0.37733e04,
        0.43962e04,
        0.50477e04,
        0.57360e04,
        0.64687e04,
        0.72525e04,
        0.80938e04,
        0.89984e04,
        0.99723e04,
        0.11021e05,
        0.12150e05,
        0.13366e05,
        0.14673e05,
        0.16079e05,
        0.17589e05,
        0.19211e05,
        0.20949e05,
        0.22812e05,
        0.24807e05,
        0.26940e05,
        0.29221e05,
        0.31656e05,
        0.34254e05,
        0.37023e05,
        0.39972e05,
        0.43111e05,
        0.46449e05,
        0.49996e05,
        0.53762e05,
        0.57756e05,
        0.61991e05,
        0.66477e05,
        0.71226e05,
        0.76249e05,
        0.81558e05,
        0.87167e05,
        0.93088e05,
        0.99334e05,
        0.10592e06,
        0.11286e06,
        0.12016e06,
        0.12785e06,
        0.13594e06,
        0.14444e06,
        0.15337e06,
        0.16274e06,
        0.17258e06,
        0.18290e06,
        0.19371e06,
        0.20504e06,
        0.21691e06,
        0.22933e06,
        0.24233e06,
        0.25592e06,
        0.27012e06,
        0.28496e06,
        0.30046e06,
        0.31663e06,
        0.33351e06,
        0.35111e06,
        0.36946e06,
        0.38858e06,
        0.40850e06,
        0.42924e06,
        0.45083e06,
        0.47329e06,
        0.49666e06,
        0.52095e06,
        0.54620e06,
        0.57243e06,
        0.59967e06,
        0.62796e06,
        0.65732e06,
        0.68778e06,
        0.71938e06,
        0.75214e06,
        0.78611e06,
        0.82131e06,
        0.85777e06,
        0.89553e06,
        0.93463e06,
        0.97511e06,
        0.10170e07,
        0.10603e07,
        0.11051e07,
        0.11514e07,
        0.11993e07,
        0.12488e07,
        0.12999e07,
        0.13527e07,
        0.14073e07,
        0.14636e07,
        0.15217e07,
        0.15816e07,
        0.16435e07,
        0.17072e07,
        0.17730e07,
        0.18408e07,
        0.19107e07,
        0.19827e07,
        0.20569e07,
        0.21334e07,
        0.22121e07,
        0.22931e07,
        0.23765e07,
        0.24624e07,
        0.25507e07,
        0.26416e07,
        0.27351e07,
        0.28312e07,
        0.29301e07,
        0.30317e07,
        0.31361e07,
        0.32434e07,
        0.33537e07,
    ]
)


#  --------------- O3 666: M = 3, I = 1 ---------------------
M = 3
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.30333e03,
        0.51126e03,
        0.75274e03,
        0.10241e04,
        0.13236e04,
        0.16508e04,
        0.20068e04,
        0.23935e04,
        0.28136e04,
        0.32703e04,
        0.37672e04,
        0.43082e04,
        0.48975e04,
        0.55395e04,
        0.62386e04,
        0.69996e04,
        0.78272e04,
        0.87264e04,
        0.97026e04,
        0.10761e05,
        0.11907e05,
        0.13146e05,
        0.14485e05,
        0.15929e05,
        0.17484e05,
        0.19158e05,
        0.20957e05,
        0.22887e05,
        0.24956e05,
        0.27172e05,
        0.29541e05,
        0.32072e05,
        0.34773e05,
        0.37652e05,
        0.40718e05,
        0.43979e05,
        0.47444e05,
        0.51123e05,
        0.55026e05,
        0.59161e05,
        0.63540e05,
        0.68172e05,
        0.73069e05,
        0.78240e05,
        0.83698e05,
        0.89453e05,
        0.95517e05,
        0.10190e06,
        0.10862e06,
        0.11569e06,
        0.12311e06,
        0.13091e06,
        0.13909e06,
        0.14767e06,
        0.15666e06,
        0.16608e06,
        0.17594e06,
        0.18626e06,
        0.19706e06,
        0.20834e06,
        0.22012e06,
        0.23242e06,
        0.24526e06,
        0.25866e06,
        0.27262e06,
        0.28717e06,
        0.30233e06,
        0.31811e06,
        0.33453e06,
        0.35161e06,
        0.36937e06,
        0.38784e06,
        0.40702e06,
        0.42694e06,
        0.44762e06,
        0.46909e06,
        0.49135e06,
        0.51444e06,
        0.53838e06,
        0.56318e06,
        0.58887e06,
        0.61548e06,
        0.64303e06,
        0.67153e06,
        0.70102e06,
        0.73153e06,
        0.76306e06,
        0.79566e06,
        0.82934e06,
        0.86413e06,
        0.90006e06,
        0.93716e06,
        0.97545e06,
        0.10150e07,
        0.10557e07,
        0.10977e07,
        0.11411e07,
        0.11858e07,
        0.12318e07,
        0.12792e07,
        0.13281e07,
        0.13784e07,
        0.14302e07,
        0.14835e07,
        0.15384e07,
        0.15948e07,
        0.16529e07,
        0.17126e07,
        0.17740e07,
        0.18371e07,
        0.19020e07,
        0.19686e07,
        0.20371e07,
        0.21074e07,
        0.21797e07,
        0.22538e07,
        0.23300e07,
        0.24081e07,
        0.24883e07,
    ]
)


#  --------------- O3 668: M = 3, I = 2 ---------------------
M = 3
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.64763e03,
        0.10916e04,
        0.16073e04,
        0.21870e04,
        0.28271e04,
        0.35272e04,
        0.42900e04,
        0.51197e04,
        0.60225e04,
        0.70057e04,
        0.80771e04,
        0.92455e04,
        0.10520e05,
        0.11911e05,
        0.13427e05,
        0.15079e05,
        0.16878e05,
        0.18834e05,
        0.20960e05,
        0.23267e05,
        0.25767e05,
        0.28472e05,
        0.31397e05,
        0.34553e05,
        0.37957e05,
        0.41620e05,
        0.45559e05,
        0.49790e05,
        0.54327e05,
        0.59187e05,
        0.64387e05,
        0.69944e05,
        0.75877e05,
        0.82203e05,
        0.88943e05,
        0.96114e05,
        0.10374e06,
        0.11184e06,
        0.12043e06,
        0.12954e06,
        0.13918e06,
        0.14939e06,
        0.16018e06,
        0.17159e06,
        0.18362e06,
        0.19632e06,
        0.20970e06,
        0.22380e06,
        0.23863e06,
        0.25423e06,
        0.27063e06,
        0.28786e06,
        0.30594e06,
        0.32490e06,
        0.34478e06,
        0.36561e06,
        0.38743e06,
        0.41026e06,
        0.43413e06,
        0.45909e06,
        0.48517e06,
        0.51241e06,
        0.54084e06,
        0.57049e06,
        0.60141e06,
        0.63365e06,
        0.66722e06,
        0.70219e06,
        0.73858e06,
        0.77644e06,
        0.81581e06,
        0.85674e06,
        0.89927e06,
        0.94345e06,
        0.98932e06,
        0.10369e07,
        0.10863e07,
        0.11375e07,
        0.11906e07,
        0.12457e07,
        0.13027e07,
        0.13618e07,
        0.14229e07,
        0.14862e07,
        0.15517e07,
        0.16194e07,
        0.16894e07,
        0.17618e07,
        0.18366e07,
        0.19139e07,
        0.19937e07,
        0.20761e07,
        0.21612e07,
        0.22490e07,
        0.23395e07,
        0.24330e07,
        0.25293e07,
        0.26286e07,
        0.27309e07,
        0.28363e07,
        0.29449e07,
        0.30568e07,
        0.31720e07,
        0.32905e07,
        0.34125e07,
        0.35381e07,
        0.36672e07,
        0.38000e07,
        0.39366e07,
        0.40770e07,
        0.42213e07,
        0.43696e07,
        0.45220e07,
        0.46785e07,
        0.48392e07,
        0.50043e07,
        0.51737e07,
        0.53476e07,
        0.55261e07,
    ]
)


#  --------------- O3 686: M = 3, I = 3 ---------------------
M = 3
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.31656e03,
        0.53355e03,
        0.78557e03,
        0.10688e04,
        0.13815e04,
        0.17235e04,
        0.20960e04,
        0.25011e04,
        0.29420e04,
        0.34223e04,
        0.39459e04,
        0.45172e04,
        0.51408e04,
        0.58213e04,
        0.65639e04,
        0.73735e04,
        0.82555e04,
        0.92152e04,
        0.10259e05,
        0.11391e05,
        0.12619e05,
        0.13949e05,
        0.15387e05,
        0.16940e05,
        0.18614e05,
        0.20417e05,
        0.22357e05,
        0.24440e05,
        0.26675e05,
        0.29070e05,
        0.31633e05,
        0.34374e05,
        0.37299e05,
        0.40420e05,
        0.43746e05,
        0.47285e05,
        0.51049e05,
        0.55047e05,
        0.59289e05,
        0.63788e05,
        0.68554e05,
        0.73598e05,
        0.78932e05,
        0.84568e05,
        0.90519e05,
        0.96796e05,
        0.10341e06,
        0.11039e06,
        0.11772e06,
        0.12544e06,
        0.13356e06,
        0.14208e06,
        0.15103e06,
        0.16041e06,
        0.17026e06,
        0.18057e06,
        0.19137e06,
        0.20268e06,
        0.21450e06,
        0.22687e06,
        0.23979e06,
        0.25328e06,
        0.26736e06,
        0.28206e06,
        0.29738e06,
        0.31336e06,
        0.33000e06,
        0.34733e06,
        0.36537e06,
        0.38414e06,
        0.40366e06,
        0.42396e06,
        0.44505e06,
        0.46696e06,
        0.48971e06,
        0.51332e06,
        0.53782e06,
        0.56323e06,
        0.58958e06,
        0.61689e06,
        0.64518e06,
        0.67448e06,
        0.70482e06,
        0.73623e06,
        0.76872e06,
        0.80234e06,
        0.83710e06,
        0.87303e06,
        0.91017e06,
        0.94853e06,
        0.98816e06,
        0.10291e07,
        0.10713e07,
        0.11149e07,
        0.11599e07,
        0.12063e07,
        0.12541e07,
        0.13034e07,
        0.13542e07,
        0.14066e07,
        0.14606e07,
        0.15161e07,
        0.15733e07,
        0.16322e07,
        0.16928e07,
        0.17552e07,
        0.18194e07,
        0.18854e07,
        0.19532e07,
        0.20230e07,
        0.20947e07,
        0.21684e07,
        0.22441e07,
        0.23219e07,
        0.24018e07,
        0.24838e07,
        0.25680e07,
        0.26545e07,
        0.27432e07,
    ]
)


#  --------------- O3 667: M = 3, I = 4 ---------------------
M = 3
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.37657e04,
        0.63472e04,
        0.93454e04,
        0.12715e05,
        0.16435e05,
        0.20502e05,
        0.24929e05,
        0.29742e05,
        0.34975e05,
        0.40668e05,
        0.46868e05,
        0.53624e05,
        0.60990e05,
        0.69018e05,
        0.77768e05,
        0.87296e05,
        0.97666e05,
        0.10894e06,
        0.12118e06,
        0.13446e06,
        0.14885e06,
        0.16441e06,
        0.18123e06,
        0.19938e06,
        0.21894e06,
        0.23998e06,
        0.26261e06,
        0.28690e06,
        0.31295e06,
        0.34084e06,
        0.37068e06,
        0.40256e06,
        0.43659e06,
        0.47287e06,
        0.51151e06,
        0.55262e06,
        0.59632e06,
        0.64272e06,
        0.69194e06,
        0.74412e06,
        0.79937e06,
        0.85783e06,
        0.91963e06,
        0.98492e06,
        0.10538e07,
        0.11265e07,
        0.12031e07,
        0.12837e07,
        0.13686e07,
        0.14579e07,
        0.15517e07,
        0.16502e07,
        0.17536e07,
        0.18621e07,
        0.19758e07,
        0.20949e07,
        0.22196e07,
        0.23501e07,
        0.24866e07,
        0.26292e07,
        0.27783e07,
        0.29339e07,
        0.30963e07,
        0.32658e07,
        0.34425e07,
        0.36266e07,
        0.38184e07,
        0.40181e07,
        0.42260e07,
        0.44422e07,
        0.46671e07,
        0.49008e07,
        0.51437e07,
        0.53959e07,
        0.56578e07,
        0.59296e07,
        0.62116e07,
        0.65040e07,
        0.68071e07,
        0.71213e07,
        0.74468e07,
        0.77838e07,
        0.81328e07,
        0.84939e07,
        0.88676e07,
        0.92541e07,
        0.96536e07,
        0.10067e08,
        0.10493e08,
        0.10934e08,
        0.11390e08,
        0.11860e08,
        0.12345e08,
        0.12846e08,
        0.13363e08,
        0.13895e08,
        0.14445e08,
        0.15011e08,
        0.15595e08,
        0.16196e08,
        0.16815e08,
        0.17453e08,
        0.18110e08,
        0.18786e08,
        0.19482e08,
        0.20198e08,
        0.20934e08,
        0.21691e08,
        0.22470e08,
        0.23270e08,
        0.24093e08,
        0.24939e08,
        0.25807e08,
        0.26699e08,
        0.27616e08,
        0.28556e08,
        0.29522e08,
        0.30514e08,
        0.31531e08,
    ]
)


#  --------------- O3 676: M = 3, I = 5 ---------------------
M = 3
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.18608e04,
        0.31363e04,
        0.46177e04,
        0.62826e04,
        0.81202e04,
        0.10129e05,
        0.12316e05,
        0.14693e05,
        0.17277e05,
        0.20089e05,
        0.23153e05,
        0.26492e05,
        0.30133e05,
        0.34103e05,
        0.38430e05,
        0.43145e05,
        0.48277e05,
        0.53858e05,
        0.59920e05,
        0.66497e05,
        0.73624e05,
        0.81336e05,
        0.89671e05,
        0.98668e05,
        0.10836e06,
        0.11880e06,
        0.13002e06,
        0.14207e06,
        0.15500e06,
        0.16884e06,
        0.18365e06,
        0.19947e06,
        0.21636e06,
        0.23438e06,
        0.25356e06,
        0.27398e06,
        0.29568e06,
        0.31873e06,
        0.34318e06,
        0.36911e06,
        0.39656e06,
        0.42561e06,
        0.45632e06,
        0.48877e06,
        0.52302e06,
        0.55914e06,
        0.59722e06,
        0.63732e06,
        0.67952e06,
        0.72390e06,
        0.77055e06,
        0.81954e06,
        0.87097e06,
        0.92491e06,
        0.98146e06,
        0.10407e07,
        0.11027e07,
        0.11677e07,
        0.12356e07,
        0.13066e07,
        0.13807e07,
        0.14582e07,
        0.15390e07,
        0.16233e07,
        0.17113e07,
        0.18029e07,
        0.18984e07,
        0.19978e07,
        0.21012e07,
        0.22089e07,
        0.23208e07,
        0.24372e07,
        0.25581e07,
        0.26837e07,
        0.28141e07,
        0.29494e07,
        0.30898e07,
        0.32354e07,
        0.33864e07,
        0.35428e07,
        0.37049e07,
        0.38728e07,
        0.40466e07,
        0.42264e07,
        0.44125e07,
        0.46050e07,
        0.48040e07,
        0.50098e07,
        0.52224e07,
        0.54420e07,
        0.56689e07,
        0.59031e07,
        0.61449e07,
        0.63943e07,
        0.66517e07,
        0.69172e07,
        0.71909e07,
        0.74731e07,
        0.77639e07,
        0.80635e07,
        0.83721e07,
        0.86900e07,
        0.90172e07,
        0.93541e07,
        0.97008e07,
        0.10058e08,
        0.10424e08,
        0.10802e08,
        0.11190e08,
        0.11589e08,
        0.11999e08,
        0.12420e08,
        0.12853e08,
        0.13298e08,
        0.13755e08,
        0.14223e08,
        0.14705e08,
        0.15199e08,
        0.15706e08,
    ]
)


#  --------------- O3 886: M = 3, I = 6 ---------------------
M = 3
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.67639e03,
        0.11401e04,
        0.16787e04,
        0.22843e04,
        0.29532e04,
        0.36856e04,
        0.44842e04,
        0.53545e04,
        0.63030e04,
        0.73381e04,
        0.84686e04,
        0.97040e04,
        0.11054e05,
        0.12530e05,
        0.14143e05,
        0.15903e05,
        0.17823e05,
        0.19915e05,
        0.22190e05,
        0.24663e05,
        0.27346e05,
        0.30254e05,
        0.33400e05,
        0.36800e05,
        0.40469e05,
        0.44423e05,
        0.48678e05,
        0.53251e05,
        0.58160e05,
        0.63423e05,
        0.69058e05,
        0.75085e05,
        0.81524e05,
        0.88395e05,
        0.95719e05,
        0.10352e06,
        0.11181e06,
        0.12063e06,
        0.12999e06,
        0.13991e06,
        0.15043e06,
        0.16157e06,
        0.17335e06,
        0.18580e06,
        0.19895e06,
        0.21283e06,
        0.22746e06,
        0.24288e06,
        0.25911e06,
        0.27619e06,
        0.29415e06,
        0.31301e06,
        0.33283e06,
        0.35362e06,
        0.37542e06,
        0.39827e06,
        0.42221e06,
        0.44726e06,
        0.47348e06,
        0.50089e06,
        0.52954e06,
        0.55947e06,
        0.59072e06,
        0.62332e06,
        0.65733e06,
        0.69279e06,
        0.72973e06,
        0.76821e06,
        0.80827e06,
        0.84996e06,
        0.89332e06,
        0.93840e06,
        0.98526e06,
        0.10339e07,
        0.10845e07,
        0.11370e07,
        0.11914e07,
        0.12479e07,
        0.13065e07,
        0.13672e07,
        0.14302e07,
        0.14953e07,
        0.15628e07,
        0.16327e07,
        0.17050e07,
        0.17798e07,
        0.18571e07,
        0.19371e07,
        0.20197e07,
        0.21051e07,
        0.21933e07,
        0.22844e07,
        0.23785e07,
        0.24755e07,
        0.25757e07,
        0.26790e07,
        0.27855e07,
        0.28954e07,
        0.30086e07,
        0.31253e07,
        0.32455e07,
        0.33693e07,
        0.34967e07,
        0.36280e07,
        0.37631e07,
        0.39021e07,
        0.40451e07,
        0.41922e07,
        0.43435e07,
        0.44990e07,
        0.46589e07,
        0.48232e07,
        0.49920e07,
        0.51654e07,
        0.53436e07,
        0.55265e07,
        0.57143e07,
        0.59071e07,
        0.61050e07,
    ]
)


#  --------------- O3 868: M = 3, I = 7 ---------------------
M = 3
I = 7
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.34615e03,
        0.58348e03,
        0.85915e03,
        0.11692e04,
        0.15117e04,
        0.18868e04,
        0.22960e04,
        0.27419e04,
        0.32278e04,
        0.37579e04,
        0.43366e04,
        0.49686e04,
        0.56591e04,
        0.64134e04,
        0.72369e04,
        0.81354e04,
        0.91148e04,
        0.10181e05,
        0.11341e05,
        0.12600e05,
        0.13966e05,
        0.15446e05,
        0.17046e05,
        0.18775e05,
        0.20640e05,
        0.22649e05,
        0.24810e05,
        0.27132e05,
        0.29624e05,
        0.32295e05,
        0.35154e05,
        0.38211e05,
        0.41475e05,
        0.44958e05,
        0.48670e05,
        0.52621e05,
        0.56823e05,
        0.61288e05,
        0.66026e05,
        0.71052e05,
        0.76376e05,
        0.82011e05,
        0.87972e05,
        0.94271e05,
        0.10092e06,
        0.10794e06,
        0.11534e06,
        0.12313e06,
        0.13134e06,
        0.13997e06,
        0.14905e06,
        0.15858e06,
        0.16859e06,
        0.17909e06,
        0.19010e06,
        0.20164e06,
        0.21373e06,
        0.22638e06,
        0.23962e06,
        0.25346e06,
        0.26792e06,
        0.28302e06,
        0.29879e06,
        0.31524e06,
        0.33240e06,
        0.35029e06,
        0.36892e06,
        0.38833e06,
        0.40853e06,
        0.42956e06,
        0.45142e06,
        0.47416e06,
        0.49778e06,
        0.52233e06,
        0.54781e06,
        0.57427e06,
        0.60172e06,
        0.63019e06,
        0.65971e06,
        0.69031e06,
        0.72201e06,
        0.75485e06,
        0.78886e06,
        0.82405e06,
        0.86048e06,
        0.89815e06,
        0.93711e06,
        0.97739e06,
        0.10190e07,
        0.10620e07,
        0.11065e07,
        0.11523e07,
        0.11997e07,
        0.12485e07,
        0.12990e07,
        0.13510e07,
        0.14046e07,
        0.14599e07,
        0.15169e07,
        0.15756e07,
        0.16361e07,
        0.16984e07,
        0.17626e07,
        0.18287e07,
        0.18966e07,
        0.19666e07,
        0.20386e07,
        0.21126e07,
        0.21887e07,
        0.22669e07,
        0.23474e07,
        0.24300e07,
        0.25150e07,
        0.26022e07,
        0.26919e07,
        0.27839e07,
        0.28784e07,
        0.29753e07,
        0.30749e07,
    ]
)


#  --------------- O3 678: M = 3, I = 8 ---------------------
M = 3
I = 8
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.39745e04,
        0.66993e04,
        0.98642e04,
        0.13422e05,
        0.17352e05,
        0.21652e05,
        0.26339e05,
        0.31442e05,
        0.37000e05,
        0.43058e05,
        0.49669e05,
        0.56885e05,
        0.64766e05,
        0.73372e05,
        0.82765e05,
        0.93011e05,
        0.10418e06,
        0.11633e06,
        0.12955e06,
        0.14390e06,
        0.15946e06,
        0.17632e06,
        0.19455e06,
        0.21424e06,
        0.23547e06,
        0.25835e06,
        0.28296e06,
        0.30939e06,
        0.33776e06,
        0.36816e06,
        0.40070e06,
        0.43549e06,
        0.47264e06,
        0.51228e06,
        0.55451e06,
        0.59947e06,
        0.64728e06,
        0.69807e06,
        0.75198e06,
        0.80915e06,
        0.86971e06,
        0.93381e06,
        0.10016e07,
        0.10733e07,
        0.11489e07,
        0.12287e07,
        0.13128e07,
        0.14015e07,
        0.14948e07,
        0.15930e07,
        0.16961e07,
        0.18045e07,
        0.19183e07,
        0.20378e07,
        0.21629e07,
        0.22942e07,
        0.24316e07,
        0.25754e07,
        0.27258e07,
        0.28831e07,
        0.30475e07,
        0.32192e07,
        0.33984e07,
        0.35855e07,
        0.37805e07,
        0.39838e07,
        0.41956e07,
        0.44162e07,
        0.46458e07,
        0.48847e07,
        0.51332e07,
        0.53916e07,
        0.56601e07,
        0.59390e07,
        0.62286e07,
        0.65292e07,
        0.68412e07,
        0.71647e07,
        0.75002e07,
        0.78479e07,
        0.82081e07,
        0.85813e07,
        0.89676e07,
        0.93676e07,
        0.97814e07,
        0.10209e08,
        0.10652e08,
        0.11110e08,
        0.11583e08,
        0.12071e08,
        0.12576e08,
        0.13097e08,
        0.13635e08,
        0.14190e08,
        0.14763e08,
        0.15354e08,
        0.15963e08,
        0.16592e08,
        0.17239e08,
        0.17906e08,
        0.18593e08,
        0.19301e08,
        0.20030e08,
        0.20780e08,
        0.21553e08,
        0.22347e08,
        0.23165e08,
        0.24006e08,
        0.24870e08,
        0.25759e08,
        0.26673e08,
        0.27612e08,
        0.28577e08,
        0.29568e08,
        0.30585e08,
        0.31631e08,
        0.32704e08,
        0.33805e08,
        0.34936e08,
    ]
)


#  --------------- O3 768: M = 3, I = 9 ---------------------
M = 3
I = 9
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.40228e04,
        0.67808e04,
        0.99842e04,
        0.13586e05,
        0.17564e05,
        0.21919e05,
        0.26665e05,
        0.31833e05,
        0.37461e05,
        0.43596e05,
        0.50286e05,
        0.57589e05,
        0.65562e05,
        0.74264e05,
        0.83761e05,
        0.94115e05,
        0.10540e06,
        0.11767e06,
        0.13102e06,
        0.14550e06,
        0.16121e06,
        0.17822e06,
        0.19661e06,
        0.21646e06,
        0.23788e06,
        0.26094e06,
        0.28574e06,
        0.31239e06,
        0.34097e06,
        0.37160e06,
        0.40437e06,
        0.43941e06,
        0.47683e06,
        0.51673e06,
        0.55925e06,
        0.60451e06,
        0.65262e06,
        0.70374e06,
        0.75799e06,
        0.81550e06,
        0.87643e06,
        0.94092e06,
        0.10091e07,
        0.10812e07,
        0.11572e07,
        0.12375e07,
        0.13221e07,
        0.14112e07,
        0.15050e07,
        0.16037e07,
        0.17074e07,
        0.18164e07,
        0.19307e07,
        0.20507e07,
        0.21765e07,
        0.23084e07,
        0.24464e07,
        0.25909e07,
        0.27421e07,
        0.29001e07,
        0.30652e07,
        0.32377e07,
        0.34177e07,
        0.36055e07,
        0.38014e07,
        0.40055e07,
        0.42182e07,
        0.44397e07,
        0.46703e07,
        0.49102e07,
        0.51597e07,
        0.54191e07,
        0.56886e07,
        0.59686e07,
        0.62593e07,
        0.65611e07,
        0.68742e07,
        0.71989e07,
        0.75356e07,
        0.78846e07,
        0.82461e07,
        0.86206e07,
        0.90083e07,
        0.94097e07,
        0.98249e07,
        0.10254e08,
        0.10699e08,
        0.11158e08,
        0.11632e08,
        0.12123e08,
        0.12629e08,
        0.13152e08,
        0.13691e08,
        0.14248e08,
        0.14823e08,
        0.15416e08,
        0.16027e08,
        0.16657e08,
        0.17307e08,
        0.17976e08,
        0.18665e08,
        0.19375e08,
        0.20106e08,
        0.20858e08,
        0.21633e08,
        0.22430e08,
        0.23250e08,
        0.24093e08,
        0.24960e08,
        0.25851e08,
        0.26767e08,
        0.27709e08,
        0.28676e08,
        0.29670e08,
        0.30691e08,
        0.31739e08,
        0.32815e08,
        0.33919e08,
        0.35053e08,
    ]
)


#  --------------- O3 786: M = 3, I = 10 ---------------------
M = 3
I = 10
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.39315e04,
        0.66267e04,
        0.97569e04,
        0.13276e05,
        0.17162e05,
        0.21414e05,
        0.26048e05,
        0.31094e05,
        0.36590e05,
        0.42581e05,
        0.49120e05,
        0.56260e05,
        0.64061e05,
        0.72580e05,
        0.81882e05,
        0.92031e05,
        0.10309e06,
        0.11514e06,
        0.12824e06,
        0.14247e06,
        0.15791e06,
        0.17463e06,
        0.19272e06,
        0.21226e06,
        0.23333e06,
        0.25604e06,
        0.28047e06,
        0.30673e06,
        0.33490e06,
        0.36510e06,
        0.39743e06,
        0.43200e06,
        0.46892e06,
        0.50831e06,
        0.55029e06,
        0.59498e06,
        0.64251e06,
        0.69301e06,
        0.74662e06,
        0.80347e06,
        0.86370e06,
        0.92747e06,
        0.99491e06,
        0.10662e07,
        0.11414e07,
        0.12208e07,
        0.13046e07,
        0.13928e07,
        0.14856e07,
        0.15833e07,
        0.16860e07,
        0.17939e07,
        0.19072e07,
        0.20261e07,
        0.21508e07,
        0.22814e07,
        0.24182e07,
        0.25614e07,
        0.27112e07,
        0.28679e07,
        0.30316e07,
        0.32026e07,
        0.33811e07,
        0.35674e07,
        0.37617e07,
        0.39642e07,
        0.41752e07,
        0.43950e07,
        0.46237e07,
        0.48618e07,
        0.51094e07,
        0.53668e07,
        0.56343e07,
        0.59123e07,
        0.62009e07,
        0.65005e07,
        0.68113e07,
        0.71338e07,
        0.74681e07,
        0.78147e07,
        0.81737e07,
        0.85457e07,
        0.89308e07,
        0.93295e07,
        0.97420e07,
        0.10169e08,
        0.10610e08,
        0.11066e08,
        0.11538e08,
        0.12025e08,
        0.12528e08,
        0.13048e08,
        0.13584e08,
        0.14138e08,
        0.14709e08,
        0.15298e08,
        0.15906e08,
        0.16532e08,
        0.17178e08,
        0.17843e08,
        0.18528e08,
        0.19234e08,
        0.19961e08,
        0.20710e08,
        0.21480e08,
        0.22272e08,
        0.23088e08,
        0.23926e08,
        0.24789e08,
        0.25675e08,
        0.26587e08,
        0.27523e08,
        0.28485e08,
        0.29474e08,
        0.30489e08,
        0.31532e08,
        0.32603e08,
        0.33701e08,
        0.34829e08,
    ]
)


#  --------------- O3 776: M = 3, I = 11 ---------------------
M = 3
I = 11
TIPS_GSI_HASH[(M, I)] = __FloatType__(36.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.23106e05,
        0.38945e05,
        0.57342e05,
        0.78021e05,
        0.10085e06,
        0.12582e06,
        0.15302e06,
        0.18262e06,
        0.21482e06,
        0.24989e06,
        0.28812e06,
        0.32983e06,
        0.37535e06,
        0.42501e06,
        0.47919e06,
        0.53825e06,
        0.60258e06,
        0.67256e06,
        0.74862e06,
        0.83118e06,
        0.92069e06,
        0.10176e07,
        0.11223e07,
        0.12354e07,
        0.13574e07,
        0.14887e07,
        0.16299e07,
        0.17816e07,
        0.19443e07,
        0.21187e07,
        0.23052e07,
        0.25047e07,
        0.27176e07,
        0.29447e07,
        0.31866e07,
        0.34441e07,
        0.37179e07,
        0.40087e07,
        0.43173e07,
        0.46444e07,
        0.49910e07,
        0.53578e07,
        0.57456e07,
        0.61554e07,
        0.65880e07,
        0.70444e07,
        0.75255e07,
        0.80322e07,
        0.85656e07,
        0.91266e07,
        0.97163e07,
        0.10336e08,
        0.10986e08,
        0.11668e08,
        0.12383e08,
        0.13133e08,
        0.13918e08,
        0.14739e08,
        0.15598e08,
        0.16496e08,
        0.17435e08,
        0.18415e08,
        0.19438e08,
        0.20505e08,
        0.21619e08,
        0.22779e08,
        0.23987e08,
        0.25246e08,
        0.26556e08,
        0.27920e08,
        0.29337e08,
        0.30811e08,
        0.32343e08,
        0.33934e08,
        0.35585e08,
        0.37300e08,
        0.39079e08,
        0.40924e08,
        0.42837e08,
        0.44819e08,
        0.46873e08,
        0.49001e08,
        0.51203e08,
        0.53483e08,
        0.55842e08,
        0.58282e08,
        0.60805e08,
        0.63414e08,
        0.66109e08,
        0.68894e08,
        0.71770e08,
        0.74740e08,
        0.77806e08,
        0.80970e08,
        0.84234e08,
        0.87600e08,
        0.91072e08,
        0.94651e08,
        0.98339e08,
        0.10214e09,
        0.10605e09,
        0.11009e09,
        0.11424e09,
        0.11851e09,
        0.12291e09,
        0.12744e09,
        0.13209e09,
        0.13688e09,
        0.14180e09,
        0.14687e09,
        0.15207e09,
        0.15742e09,
        0.16291e09,
        0.16855e09,
        0.17435e09,
        0.18030e09,
        0.18641e09,
        0.19268e09,
        0.19912e09,
    ]
)


#  --------------- O3 767: M = 3, I = 12 ---------------------
M = 3
I = 12
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.11692e05,
        0.19707e05,
        0.29017e05,
        0.39482e05,
        0.51038e05,
        0.63680e05,
        0.77450e05,
        0.92432e05,
        0.10873e06,
        0.12649e06,
        0.14584e06,
        0.16694e06,
        0.18996e06,
        0.21507e06,
        0.24245e06,
        0.27229e06,
        0.30478e06,
        0.34013e06,
        0.37853e06,
        0.42020e06,
        0.46536e06,
        0.51424e06,
        0.56708e06,
        0.62411e06,
        0.68559e06,
        0.75178e06,
        0.82296e06,
        0.89939e06,
        0.98137e06,
        0.10692e07,
        0.11631e07,
        0.12636e07,
        0.13708e07,
        0.14851e07,
        0.16069e07,
        0.17365e07,
        0.18742e07,
        0.20206e07,
        0.21758e07,
        0.23404e07,
        0.25148e07,
        0.26992e07,
        0.28943e07,
        0.31004e07,
        0.33179e07,
        0.35474e07,
        0.37892e07,
        0.40440e07,
        0.43121e07,
        0.45940e07,
        0.48904e07,
        0.52017e07,
        0.55285e07,
        0.58713e07,
        0.62306e07,
        0.66071e07,
        0.70014e07,
        0.74140e07,
        0.78456e07,
        0.82967e07,
        0.87681e07,
        0.92604e07,
        0.97742e07,
        0.10310e08,
        0.10869e08,
        0.11452e08,
        0.12059e08,
        0.12691e08,
        0.13348e08,
        0.14033e08,
        0.14745e08,
        0.15484e08,
        0.16253e08,
        0.17052e08,
        0.17881e08,
        0.18741e08,
        0.19634e08,
        0.20560e08,
        0.21520e08,
        0.22515e08,
        0.23546e08,
        0.24613e08,
        0.25718e08,
        0.26862e08,
        0.28046e08,
        0.29270e08,
        0.30536e08,
        0.31845e08,
        0.33197e08,
        0.34594e08,
        0.36037e08,
        0.37527e08,
        0.39065e08,
        0.40652e08,
        0.42289e08,
        0.43977e08,
        0.45719e08,
        0.47514e08,
        0.49363e08,
        0.51270e08,
        0.53233e08,
        0.55255e08,
        0.57337e08,
        0.59480e08,
        0.61686e08,
        0.63956e08,
        0.66290e08,
        0.68691e08,
        0.71160e08,
        0.73699e08,
        0.76307e08,
        0.78988e08,
        0.81743e08,
        0.84572e08,
        0.87478e08,
        0.90462e08,
        0.93525e08,
        0.96669e08,
        0.99896e08,
    ]
)


#  --------------- O3 888: M = 3, I = 13 ---------------------
M = 3
I = 13
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.36175e03,
        0.60978e03,
        0.89790e03,
        0.12219e04,
        0.15802e04,
        0.19728e04,
        0.24016e04,
        0.28696e04,
        0.33807e04,
        0.39394e04,
        0.45506e04,
        0.52196e04,
        0.59521e04,
        0.67538e04,
        0.76308e04,
        0.85894e04,
        0.96361e04,
        0.10777e05,
        0.12021e05,
        0.13373e05,
        0.14841e05,
        0.16434e05,
        0.18158e05,
        0.20023e05,
        0.22037e05,
        0.24208e05,
        0.26547e05,
        0.29061e05,
        0.31762e05,
        0.34659e05,
        0.37762e05,
        0.41083e05,
        0.44632e05,
        0.48421e05,
        0.52462e05,
        0.56766e05,
        0.61346e05,
        0.66215e05,
        0.71386e05,
        0.76873e05,
        0.82688e05,
        0.88848e05,
        0.95365e05,
        0.10226e06,
        0.10954e06,
        0.11722e06,
        0.12532e06,
        0.13387e06,
        0.14286e06,
        0.15233e06,
        0.16229e06,
        0.17275e06,
        0.18374e06,
        0.19528e06,
        0.20737e06,
        0.22006e06,
        0.23335e06,
        0.24726e06,
        0.26182e06,
        0.27705e06,
        0.29297e06,
        0.30960e06,
        0.32696e06,
        0.34509e06,
        0.36399e06,
        0.38371e06,
        0.40425e06,
        0.42566e06,
        0.44794e06,
        0.47114e06,
        0.49527e06,
        0.52036e06,
        0.54644e06,
        0.57354e06,
        0.60169e06,
        0.63091e06,
        0.66124e06,
        0.69270e06,
        0.72533e06,
        0.75916e06,
        0.79421e06,
        0.83053e06,
        0.86814e06,
        0.90708e06,
        0.94737e06,
        0.98907e06,
        0.10322e07,
        0.10768e07,
        0.11229e07,
        0.11705e07,
        0.12197e07,
        0.12705e07,
        0.13230e07,
        0.13771e07,
        0.14330e07,
        0.14906e07,
        0.15501e07,
        0.16114e07,
        0.16745e07,
        0.17397e07,
        0.18067e07,
        0.18759e07,
        0.19470e07,
        0.20203e07,
        0.20957e07,
        0.21733e07,
        0.22532e07,
        0.23353e07,
        0.24198e07,
        0.25067e07,
        0.25960e07,
        0.26878e07,
        0.27821e07,
        0.28790e07,
        0.29785e07,
        0.30807e07,
        0.31857e07,
        0.32934e07,
        0.34040e07,
    ]
)


#  --------------- O3 887: M = 3, I = 14 ---------------------
M = 3
I = 14
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.42000e04,
        0.70796e04,
        0.10424e05,
        0.14186e05,
        0.18342e05,
        0.22896e05,
        0.27866e05,
        0.33285e05,
        0.39199e05,
        0.45659e05,
        0.52720e05,
        0.60444e05,
        0.68895e05,
        0.78139e05,
        0.88246e05,
        0.99288e05,
        0.11134e06,
        0.12447e06,
        0.13877e06,
        0.15431e06,
        0.17119e06,
        0.18949e06,
        0.20930e06,
        0.23071e06,
        0.25383e06,
        0.27875e06,
        0.30558e06,
        0.33442e06,
        0.36539e06,
        0.39861e06,
        0.43418e06,
        0.47224e06,
        0.51291e06,
        0.55632e06,
        0.60260e06,
        0.65189e06,
        0.70434e06,
        0.76008e06,
        0.81927e06,
        0.88206e06,
        0.94862e06,
        0.10191e07,
        0.10937e07,
        0.11725e07,
        0.12558e07,
        0.13436e07,
        0.14363e07,
        0.15340e07,
        0.16368e07,
        0.17450e07,
        0.18588e07,
        0.19784e07,
        0.21040e07,
        0.22358e07,
        0.23741e07,
        0.25190e07,
        0.26708e07,
        0.28297e07,
        0.29961e07,
        0.31700e07,
        0.33518e07,
        0.35417e07,
        0.37400e07,
        0.39469e07,
        0.41628e07,
        0.43878e07,
        0.46224e07,
        0.48667e07,
        0.51210e07,
        0.53858e07,
        0.56611e07,
        0.59475e07,
        0.62451e07,
        0.65544e07,
        0.68755e07,
        0.72089e07,
        0.75550e07,
        0.79139e07,
        0.82861e07,
        0.86720e07,
        0.90719e07,
        0.94861e07,
        0.99151e07,
        0.10359e08,
        0.10819e08,
        0.11294e08,
        0.11786e08,
        0.12294e08,
        0.12820e08,
        0.13363e08,
        0.13924e08,
        0.14503e08,
        0.15101e08,
        0.15719e08,
        0.16356e08,
        0.17013e08,
        0.17690e08,
        0.18389e08,
        0.19109e08,
        0.19851e08,
        0.20616e08,
        0.21404e08,
        0.22215e08,
        0.23050e08,
        0.23910e08,
        0.24794e08,
        0.25704e08,
        0.26640e08,
        0.27603e08,
        0.28593e08,
        0.29610e08,
        0.30656e08,
        0.31731e08,
        0.32835e08,
        0.33969e08,
        0.35133e08,
        0.36329e08,
        0.37556e08,
        0.38816e08,
    ]
)


#  --------------- O3 878: M = 3, I = 15 ---------------------
M = 3
I = 15
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.21250e04,
        0.35820e04,
        0.52744e04,
        0.71778e04,
        0.92814e04,
        0.11586e05,
        0.14102e05,
        0.16845e05,
        0.19839e05,
        0.23108e05,
        0.26680e05,
        0.30588e05,
        0.34861e05,
        0.39534e05,
        0.44642e05,
        0.50219e05,
        0.56305e05,
        0.62937e05,
        0.70155e05,
        0.78001e05,
        0.86516e05,
        0.95747e05,
        0.10574e06,
        0.11653e06,
        0.12819e06,
        0.14075e06,
        0.15427e06,
        0.16881e06,
        0.18441e06,
        0.20114e06,
        0.21906e06,
        0.23823e06,
        0.25871e06,
        0.28056e06,
        0.30386e06,
        0.32867e06,
        0.35507e06,
        0.38312e06,
        0.41291e06,
        0.44450e06,
        0.47799e06,
        0.51344e06,
        0.55095e06,
        0.59060e06,
        0.63248e06,
        0.67667e06,
        0.72327e06,
        0.77238e06,
        0.82409e06,
        0.87850e06,
        0.93571e06,
        0.99583e06,
        0.10590e07,
        0.11252e07,
        0.11947e07,
        0.12675e07,
        0.13438e07,
        0.14237e07,
        0.15072e07,
        0.15946e07,
        0.16859e07,
        0.17814e07,
        0.18810e07,
        0.19849e07,
        0.20934e07,
        0.22064e07,
        0.23242e07,
        0.24469e07,
        0.25747e07,
        0.27076e07,
        0.28459e07,
        0.29897e07,
        0.31391e07,
        0.32944e07,
        0.34557e07,
        0.36231e07,
        0.37968e07,
        0.39770e07,
        0.41639e07,
        0.43576e07,
        0.45583e07,
        0.47663e07,
        0.49816e07,
        0.52045e07,
        0.54352e07,
        0.56739e07,
        0.59207e07,
        0.61759e07,
        0.64396e07,
        0.67121e07,
        0.69936e07,
        0.72844e07,
        0.75845e07,
        0.78943e07,
        0.82139e07,
        0.85436e07,
        0.88837e07,
        0.92342e07,
        0.95956e07,
        0.99680e07,
        0.10352e08,
        0.10747e08,
        0.11154e08,
        0.11573e08,
        0.12004e08,
        0.12448e08,
        0.12904e08,
        0.13374e08,
        0.13857e08,
        0.14353e08,
        0.14864e08,
        0.15388e08,
        0.15927e08,
        0.16481e08,
        0.17050e08,
        0.17634e08,
        0.18234e08,
        0.18849e08,
        0.19481e08,
    ]
)


#  --------------- O3 778: M = 3, I = 16 ---------------------
M = 3
I = 16
TIPS_GSI_HASH[(M, I)] = __FloatType__(36.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.24692e05,
        0.41621e05,
        0.61284e05,
        0.83394e05,
        0.10782e06,
        0.13457e06,
        0.16375e06,
        0.19554e06,
        0.23020e06,
        0.26801e06,
        0.30930e06,
        0.35443e06,
        0.40375e06,
        0.45763e06,
        0.51650e06,
        0.58075e06,
        0.65080e06,
        0.72711e06,
        0.81012e06,
        0.90030e06,
        0.99815e06,
        0.11042e07,
        0.12189e07,
        0.13428e07,
        0.14765e07,
        0.16206e07,
        0.17757e07,
        0.19423e07,
        0.21212e07,
        0.23129e07,
        0.25181e07,
        0.27377e07,
        0.29721e07,
        0.32223e07,
        0.34890e07,
        0.37729e07,
        0.40750e07,
        0.43959e07,
        0.47365e07,
        0.50978e07,
        0.54807e07,
        0.58860e07,
        0.63147e07,
        0.67678e07,
        0.72463e07,
        0.77512e07,
        0.82836e07,
        0.88445e07,
        0.94351e07,
        0.10056e08,
        0.10710e08,
        0.11396e08,
        0.12117e08,
        0.12873e08,
        0.13666e08,
        0.14497e08,
        0.15367e08,
        0.16279e08,
        0.17232e08,
        0.18229e08,
        0.19271e08,
        0.20359e08,
        0.21495e08,
        0.22681e08,
        0.23917e08,
        0.25206e08,
        0.26549e08,
        0.27948e08,
        0.29404e08,
        0.30920e08,
        0.32496e08,
        0.34135e08,
        0.35838e08,
        0.37608e08,
        0.39445e08,
        0.41353e08,
        0.43332e08,
        0.45385e08,
        0.47514e08,
        0.49721e08,
        0.52007e08,
        0.54376e08,
        0.56829e08,
        0.59367e08,
        0.61995e08,
        0.64712e08,
        0.67523e08,
        0.70429e08,
        0.73432e08,
        0.76535e08,
        0.79740e08,
        0.83050e08,
        0.86467e08,
        0.89993e08,
        0.93632e08,
        0.97385e08,
        0.10126e09,
        0.10525e09,
        0.10936e09,
        0.11360e09,
        0.11796e09,
        0.12246e09,
        0.12709e09,
        0.13186e09,
        0.13677e09,
        0.14182e09,
        0.14701e09,
        0.15236e09,
        0.15785e09,
        0.16350e09,
        0.16931e09,
        0.17528e09,
        0.18141e09,
        0.18771e09,
        0.19418e09,
        0.20082e09,
        0.20764e09,
        0.21465e09,
        0.22183e09,
    ]
)


#  --------------- O3 787: M = 3, I = 17 ---------------------
M = 3
I = 17
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.12211e05,
        0.20582e05,
        0.30305e05,
        0.41237e05,
        0.53314e05,
        0.66536e05,
        0.80957e05,
        0.96672e05,
        0.11380e06,
        0.13250e06,
        0.15292e06,
        0.17524e06,
        0.19965e06,
        0.22632e06,
        0.25546e06,
        0.28728e06,
        0.32199e06,
        0.35980e06,
        0.40094e06,
        0.44565e06,
        0.49417e06,
        0.54676e06,
        0.60366e06,
        0.66516e06,
        0.73152e06,
        0.80305e06,
        0.88002e06,
        0.96276e06,
        0.10516e07,
        0.11468e07,
        0.12488e07,
        0.13578e07,
        0.14743e07,
        0.15987e07,
        0.17312e07,
        0.18723e07,
        0.20225e07,
        0.21820e07,
        0.23514e07,
        0.25310e07,
        0.27214e07,
        0.29230e07,
        0.31362e07,
        0.33616e07,
        0.35997e07,
        0.38509e07,
        0.41158e07,
        0.43949e07,
        0.46887e07,
        0.49980e07,
        0.53231e07,
        0.56647e07,
        0.60234e07,
        0.63998e07,
        0.67946e07,
        0.72084e07,
        0.76418e07,
        0.80955e07,
        0.85702e07,
        0.90666e07,
        0.95854e07,
        0.10127e08,
        0.10693e08,
        0.11284e08,
        0.11900e08,
        0.12542e08,
        0.13211e08,
        0.13907e08,
        0.14633e08,
        0.15388e08,
        0.16173e08,
        0.16990e08,
        0.17838e08,
        0.18720e08,
        0.19636e08,
        0.20586e08,
        0.21573e08,
        0.22596e08,
        0.23657e08,
        0.24757e08,
        0.25896e08,
        0.27077e08,
        0.28299e08,
        0.29565e08,
        0.30874e08,
        0.32229e08,
        0.33630e08,
        0.35079e08,
        0.36576e08,
        0.38123e08,
        0.39721e08,
        0.41371e08,
        0.43075e08,
        0.44833e08,
        0.46647e08,
        0.48518e08,
        0.50448e08,
        0.52438e08,
        0.54489e08,
        0.56603e08,
        0.58780e08,
        0.61023e08,
        0.63332e08,
        0.65710e08,
        0.68157e08,
        0.70676e08,
        0.73266e08,
        0.75931e08,
        0.78672e08,
        0.81490e08,
        0.84386e08,
        0.87363e08,
        0.90422e08,
        0.93564e08,
        0.96791e08,
        0.10011e09,
        0.10351e09,
        0.10700e09,
        0.11059e09,
    ]
)


#  --------------- O3 777: M = 3, I = 18 ---------------------
M = 3
I = 18
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.71750e05,
        0.12094e06,
        0.17807e06,
        0.24230e06,
        0.31324e06,
        0.39088e06,
        0.47550e06,
        0.56764e06,
        0.66800e06,
        0.77740e06,
        0.89677e06,
        0.10271e07,
        0.11694e07,
        0.13249e07,
        0.14945e07,
        0.16796e07,
        0.18813e07,
        0.21009e07,
        0.23396e07,
        0.25989e07,
        0.28801e07,
        0.31847e07,
        0.35140e07,
        0.38698e07,
        0.42535e07,
        0.46669e07,
        0.51115e07,
        0.55893e07,
        0.61019e07,
        0.66513e07,
        0.72393e07,
        0.78680e07,
        0.85395e07,
        0.92558e07,
        0.10019e08,
        0.10832e08,
        0.11696e08,
        0.12614e08,
        0.13588e08,
        0.14621e08,
        0.15716e08,
        0.16875e08,
        0.18100e08,
        0.19395e08,
        0.20762e08,
        0.22205e08,
        0.23726e08,
        0.25328e08,
        0.27015e08,
        0.28789e08,
        0.30654e08,
        0.32614e08,
        0.34671e08,
        0.36830e08,
        0.39093e08,
        0.41465e08,
        0.43949e08,
        0.46549e08,
        0.49269e08,
        0.52112e08,
        0.55084e08,
        0.58188e08,
        0.61428e08,
        0.64809e08,
        0.68335e08,
        0.72010e08,
        0.75840e08,
        0.79828e08,
        0.83979e08,
        0.88299e08,
        0.92792e08,
        0.97463e08,
        0.10232e09,
        0.10736e09,
        0.11260e09,
        0.11803e09,
        0.12367e09,
        0.12952e09,
        0.13559e09,
        0.14187e09,
        0.14839e09,
        0.15513e09,
        0.16212e09,
        0.16935e09,
        0.17683e09,
        0.18457e09,
        0.19257e09,
        0.20085e09,
        0.20940e09,
        0.21824e09,
        0.22736e09,
        0.23678e09,
        0.24651e09,
        0.25655e09,
        0.26691e09,
        0.27759e09,
        0.28861e09,
        0.29997e09,
        0.31167e09,
        0.32374e09,
        0.33616e09,
        0.34896e09,
        0.36214e09,
        0.37571e09,
        0.38967e09,
        0.40404e09,
        0.41882e09,
        0.43403e09,
        0.44966e09,
        0.46573e09,
        0.48226e09,
        0.49923e09,
        0.51668e09,
        0.53460e09,
        0.55301e09,
        0.57191e09,
        0.59131e09,
        0.61123e09,
        0.63167e09,
    ]
)


#  --------------- N2O 446: M = 4, I = 1 ---------------------
M = 4
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(9.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.89943e03,
        0.12734e04,
        0.16489e04,
        0.20293e04,
        0.24205e04,
        0.28289e04,
        0.32609e04,
        0.37222e04,
        0.42180e04,
        0.47529e04,
        0.53312e04,
        0.59572e04,
        0.66348e04,
        0.73683e04,
        0.81616e04,
        0.90190e04,
        0.99450e04,
        0.10944e05,
        0.12021e05,
        0.13180e05,
        0.14426e05,
        0.15766e05,
        0.17203e05,
        0.18745e05,
        0.20396e05,
        0.22162e05,
        0.24051e05,
        0.26069e05,
        0.28222e05,
        0.30517e05,
        0.32962e05,
        0.35564e05,
        0.38331e05,
        0.41271e05,
        0.44393e05,
        0.47704e05,
        0.51214e05,
        0.54932e05,
        0.58868e05,
        0.63030e05,
        0.67429e05,
        0.72075e05,
        0.76979e05,
        0.82151e05,
        0.87604e05,
        0.93348e05,
        0.99395e05,
        0.10576e06,
        0.11245e06,
        0.11948e06,
        0.12686e06,
        0.13461e06,
        0.14275e06,
        0.15128e06,
        0.16021e06,
        0.16958e06,
        0.17938e06,
        0.18964e06,
        0.20037e06,
        0.21159e06,
        0.22331e06,
        0.23556e06,
        0.24834e06,
        0.26169e06,
        0.27561e06,
        0.29012e06,
        0.30525e06,
        0.32101e06,
        0.33743e06,
        0.35452e06,
        0.37230e06,
        0.39080e06,
        0.41004e06,
        0.43004e06,
        0.45082e06,
        0.47241e06,
        0.49483e06,
        0.51810e06,
        0.54225e06,
        0.56730e06,
        0.59329e06,
        0.62022e06,
        0.64814e06,
        0.67707e06,
        0.70703e06,
        0.73806e06,
        0.77018e06,
        0.80342e06,
        0.83781e06,
        0.87338e06,
        0.91016e06,
        0.94818e06,
        0.98748e06,
        0.10281e07,
        0.10700e07,
        0.11133e07,
        0.11581e07,
        0.12042e07,
        0.12519e07,
        0.13010e07,
        0.13517e07,
        0.14040e07,
        0.14579e07,
        0.15134e07,
        0.15707e07,
        0.16297e07,
        0.16905e07,
        0.17530e07,
        0.18175e07,
        0.18838e07,
        0.19521e07,
        0.20224e07,
        0.20947e07,
        0.21690e07,
        0.22455e07,
        0.23242e07,
        0.24050e07,
        0.24881e07,
        0.25735e07,
    ]
)


#  --------------- N2O 456: M = 4, I = 2 ---------------------
M = 4
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.59966e03,
        0.84903e03,
        0.10995e04,
        0.13538e04,
        0.16158e04,
        0.18903e04,
        0.21815e04,
        0.24934e04,
        0.28295e04,
        0.31927e04,
        0.35862e04,
        0.40128e04,
        0.44752e04,
        0.49763e04,
        0.55189e04,
        0.61059e04,
        0.67404e04,
        0.74256e04,
        0.81646e04,
        0.89609e04,
        0.98180e04,
        0.10740e05,
        0.11729e05,
        0.12791e05,
        0.13930e05,
        0.15149e05,
        0.16453e05,
        0.17847e05,
        0.19335e05,
        0.20922e05,
        0.22614e05,
        0.24416e05,
        0.26333e05,
        0.28371e05,
        0.30535e05,
        0.32833e05,
        0.35269e05,
        0.37851e05,
        0.40585e05,
        0.43478e05,
        0.46537e05,
        0.49769e05,
        0.53182e05,
        0.56783e05,
        0.60580e05,
        0.64582e05,
        0.68796e05,
        0.73232e05,
        0.77898e05,
        0.82803e05,
        0.87957e05,
        0.93369e05,
        0.99048e05,
        0.10501e06,
        0.11125e06,
        0.11780e06,
        0.12465e06,
        0.13182e06,
        0.13933e06,
        0.14718e06,
        0.15539e06,
        0.16396e06,
        0.17291e06,
        0.18226e06,
        0.19201e06,
        0.20218e06,
        0.21278e06,
        0.22383e06,
        0.23534e06,
        0.24733e06,
        0.25980e06,
        0.27278e06,
        0.28628e06,
        0.30032e06,
        0.31491e06,
        0.33007e06,
        0.34581e06,
        0.36216e06,
        0.37912e06,
        0.39673e06,
        0.41499e06,
        0.43392e06,
        0.45355e06,
        0.47389e06,
        0.49496e06,
        0.51678e06,
        0.53937e06,
        0.56276e06,
        0.58695e06,
        0.61199e06,
        0.63788e06,
        0.66464e06,
        0.69231e06,
        0.72090e06,
        0.75044e06,
        0.78094e06,
        0.81244e06,
        0.84496e06,
        0.87853e06,
        0.91316e06,
        0.94889e06,
        0.98573e06,
        0.10237e07,
        0.10629e07,
        0.11033e07,
        0.11449e07,
        0.11877e07,
        0.12319e07,
        0.12773e07,
        0.13241e07,
        0.13723e07,
        0.14219e07,
        0.14729e07,
        0.15254e07,
        0.15793e07,
        0.16349e07,
        0.16919e07,
        0.17506e07,
        0.18109e07,
    ]
)


#  --------------- N2O 546: M = 4, I = 3 ---------------------
M = 4
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.62051e03,
        0.87856e03,
        0.11377e04,
        0.14003e04,
        0.16705e04,
        0.19529e04,
        0.22518e04,
        0.25713e04,
        0.29149e04,
        0.32859e04,
        0.36873e04,
        0.41220e04,
        0.45929e04,
        0.51028e04,
        0.56547e04,
        0.62515e04,
        0.68963e04,
        0.75923e04,
        0.83428e04,
        0.91511e04,
        0.10021e05,
        0.10956e05,
        0.11960e05,
        0.13036e05,
        0.14190e05,
        0.15425e05,
        0.16746e05,
        0.18158e05,
        0.19664e05,
        0.21271e05,
        0.22984e05,
        0.24806e05,
        0.26745e05,
        0.28806e05,
        0.30995e05,
        0.33317e05,
        0.35780e05,
        0.38389e05,
        0.41151e05,
        0.44073e05,
        0.47162e05,
        0.50425e05,
        0.53871e05,
        0.57505e05,
        0.61338e05,
        0.65375e05,
        0.69628e05,
        0.74102e05,
        0.78808e05,
        0.83755e05,
        0.88951e05,
        0.94407e05,
        0.10013e06,
        0.10614e06,
        0.11243e06,
        0.11902e06,
        0.12593e06,
        0.13316e06,
        0.14072e06,
        0.14862e06,
        0.15689e06,
        0.16552e06,
        0.17453e06,
        0.18394e06,
        0.19376e06,
        0.20399e06,
        0.21466e06,
        0.22578e06,
        0.23737e06,
        0.24942e06,
        0.26198e06,
        0.27503e06,
        0.28861e06,
        0.30273e06,
        0.31741e06,
        0.33265e06,
        0.34848e06,
        0.36492e06,
        0.38197e06,
        0.39967e06,
        0.41803e06,
        0.43706e06,
        0.45679e06,
        0.47723e06,
        0.49840e06,
        0.52033e06,
        0.54303e06,
        0.56653e06,
        0.59084e06,
        0.61599e06,
        0.64200e06,
        0.66888e06,
        0.69667e06,
        0.72539e06,
        0.75506e06,
        0.78569e06,
        0.81733e06,
        0.84998e06,
        0.88369e06,
        0.91846e06,
        0.95433e06,
        0.99132e06,
        0.10295e07,
        0.10688e07,
        0.11093e07,
        0.11511e07,
        0.11941e07,
        0.12384e07,
        0.12840e07,
        0.13310e07,
        0.13793e07,
        0.14291e07,
        0.14803e07,
        0.15329e07,
        0.15871e07,
        0.16428e07,
        0.17000e07,
        0.17589e07,
        0.18194e07,
    ]
)


#  --------------- N2O 448: M = 4, I = 4 ---------------------
M = 4
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(9.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.95253e03,
        0.13487e04,
        0.17465e04,
        0.21498e04,
        0.25648e04,
        0.29986e04,
        0.34580e04,
        0.39493e04,
        0.44779e04,
        0.50488e04,
        0.56669e04,
        0.63366e04,
        0.70625e04,
        0.78488e04,
        0.87003e04,
        0.96216e04,
        0.10617e05,
        0.11692e05,
        0.12852e05,
        0.14102e05,
        0.15447e05,
        0.16893e05,
        0.18446e05,
        0.20112e05,
        0.21898e05,
        0.23811e05,
        0.25856e05,
        0.28042e05,
        0.30377e05,
        0.32866e05,
        0.35520e05,
        0.38345e05,
        0.41351e05,
        0.44545e05,
        0.47939e05,
        0.51540e05,
        0.55359e05,
        0.59405e05,
        0.63689e05,
        0.68222e05,
        0.73015e05,
        0.78078e05,
        0.83424e05,
        0.89064e05,
        0.95012e05,
        0.10128e06,
        0.10788e06,
        0.11482e06,
        0.12213e06,
        0.12981e06,
        0.13788e06,
        0.14635e06,
        0.15524e06,
        0.16456e06,
        0.17433e06,
        0.18457e06,
        0.19530e06,
        0.20652e06,
        0.21827e06,
        0.23055e06,
        0.24338e06,
        0.25679e06,
        0.27079e06,
        0.28541e06,
        0.30066e06,
        0.31656e06,
        0.33314e06,
        0.35042e06,
        0.36841e06,
        0.38715e06,
        0.40666e06,
        0.42695e06,
        0.44805e06,
        0.46999e06,
        0.49279e06,
        0.51649e06,
        0.54109e06,
        0.56664e06,
        0.59315e06,
        0.62066e06,
        0.64919e06,
        0.67877e06,
        0.70943e06,
        0.74121e06,
        0.77413e06,
        0.80822e06,
        0.84351e06,
        0.88004e06,
        0.91783e06,
        0.95693e06,
        0.99737e06,
        0.10392e07,
        0.10824e07,
        0.11270e07,
        0.11732e07,
        0.12208e07,
        0.12700e07,
        0.13208e07,
        0.13732e07,
        0.14272e07,
        0.14830e07,
        0.15405e07,
        0.15999e07,
        0.16610e07,
        0.17240e07,
        0.17890e07,
        0.18559e07,
        0.19248e07,
        0.19957e07,
        0.20687e07,
        0.21439e07,
        0.22213e07,
        0.23009e07,
        0.23828e07,
        0.24671e07,
        0.25537e07,
        0.26428e07,
        0.27343e07,
        0.28284e07,
    ]
)


#  --------------- N2O 447: M = 4, I = 5 ---------------------
M = 4
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(54.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.55598e04,
        0.78718e04,
        0.10193e05,
        0.12546e05,
        0.14966e05,
        0.17495e05,
        0.20171e05,
        0.23031e05,
        0.26106e05,
        0.29426e05,
        0.33018e05,
        0.36908e05,
        0.41121e05,
        0.45684e05,
        0.50622e05,
        0.55962e05,
        0.61731e05,
        0.67958e05,
        0.74671e05,
        0.81902e05,
        0.89681e05,
        0.98043e05,
        0.10702e06,
        0.11665e06,
        0.12697e06,
        0.13801e06,
        0.14983e06,
        0.16244e06,
        0.17591e06,
        0.19028e06,
        0.20558e06,
        0.22188e06,
        0.23920e06,
        0.25762e06,
        0.27718e06,
        0.29793e06,
        0.31993e06,
        0.34323e06,
        0.36791e06,
        0.39401e06,
        0.42160e06,
        0.45074e06,
        0.48151e06,
        0.51397e06,
        0.54819e06,
        0.58424e06,
        0.62221e06,
        0.66215e06,
        0.70416e06,
        0.74832e06,
        0.79470e06,
        0.84340e06,
        0.89450e06,
        0.94808e06,
        0.10042e07,
        0.10631e07,
        0.11247e07,
        0.11892e07,
        0.12567e07,
        0.13272e07,
        0.14009e07,
        0.14779e07,
        0.15583e07,
        0.16422e07,
        0.17298e07,
        0.18211e07,
        0.19163e07,
        0.20154e07,
        0.21187e07,
        0.22263e07,
        0.23382e07,
        0.24546e07,
        0.25757e07,
        0.27016e07,
        0.28324e07,
        0.29683e07,
        0.31095e07,
        0.32560e07,
        0.34081e07,
        0.35659e07,
        0.37295e07,
        0.38991e07,
        0.40750e07,
        0.42572e07,
        0.44459e07,
        0.46414e07,
        0.48437e07,
        0.50531e07,
        0.52698e07,
        0.54939e07,
        0.57257e07,
        0.59653e07,
        0.62129e07,
        0.64688e07,
        0.67331e07,
        0.70061e07,
        0.72880e07,
        0.75790e07,
        0.78792e07,
        0.81891e07,
        0.85086e07,
        0.88382e07,
        0.91780e07,
        0.95283e07,
        0.98893e07,
        0.10261e08,
        0.10644e08,
        0.11039e08,
        0.11445e08,
        0.11864e08,
        0.12294e08,
        0.12738e08,
        0.13194e08,
        0.13663e08,
        0.14145e08,
        0.14641e08,
        0.15151e08,
        0.15675e08,
        0.16214e08,
    ]
)


#  --------------- CO 26: M = 5, I = 1 ---------------------
M = 5
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.21948e02,
        0.30961e02,
        0.39980e02,
        0.49004e02,
        0.58035e02,
        0.67071e02,
        0.76112e02,
        0.85160e02,
        0.94213e02,
        0.10327e03,
        0.11234e03,
        0.12142e03,
        0.13050e03,
        0.13960e03,
        0.14872e03,
        0.15787e03,
        0.16704e03,
        0.17624e03,
        0.18548e03,
        0.19477e03,
        0.20411e03,
        0.21350e03,
        0.22295e03,
        0.23248e03,
        0.24207e03,
        0.25175e03,
        0.26151e03,
        0.27136e03,
        0.28130e03,
        0.29134e03,
        0.30148e03,
        0.31172e03,
        0.32207e03,
        0.33253e03,
        0.34312e03,
        0.35381e03,
        0.36463e03,
        0.37557e03,
        0.38663e03,
        0.39782e03,
        0.40914e03,
        0.42060e03,
        0.43218e03,
        0.44389e03,
        0.45575e03,
        0.46774e03,
        0.47987e03,
        0.49213e03,
        0.50454e03,
        0.51708e03,
        0.52978e03,
        0.54261e03,
        0.55559e03,
        0.56871e03,
        0.58198e03,
        0.59540e03,
        0.60896e03,
        0.62267e03,
        0.63653e03,
        0.65055e03,
        0.66470e03,
        0.67901e03,
        0.69347e03,
        0.70808e03,
        0.72284e03,
        0.73776e03,
        0.75283e03,
        0.76805e03,
        0.78342e03,
        0.79895e03,
        0.81463e03,
        0.83047e03,
        0.84646e03,
        0.86260e03,
        0.87891e03,
        0.89536e03,
        0.91197e03,
        0.92874e03,
        0.94566e03,
        0.96275e03,
        0.97998e03,
        0.99738e03,
        0.10149e04,
        0.10326e04,
        0.10505e04,
        0.10685e04,
        0.10867e04,
        0.11051e04,
        0.11236e04,
        0.11422e04,
        0.11611e04,
        0.11800e04,
        0.11992e04,
        0.12185e04,
        0.12380e04,
        0.12576e04,
        0.12774e04,
        0.12973e04,
        0.13174e04,
        0.13377e04,
        0.13581e04,
        0.13787e04,
        0.13994e04,
        0.14203e04,
        0.14414e04,
        0.14627e04,
        0.14841e04,
        0.15056e04,
        0.15273e04,
        0.15492e04,
        0.15713e04,
        0.15935e04,
        0.16159e04,
        0.16384e04,
        0.16611e04,
        0.16840e04,
        0.17070e04,
        0.17302e04,
        0.17536e04,
    ]
)


#  --------------- CO 36: M = 5, I = 2 ---------------------
M = 5
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.45888e02,
        0.64745e02,
        0.83615e02,
        0.10250e03,
        0.12139e03,
        0.14030e03,
        0.15921e03,
        0.17814e03,
        0.19708e03,
        0.21604e03,
        0.23501e03,
        0.25400e03,
        0.27302e03,
        0.29207e03,
        0.31117e03,
        0.33031e03,
        0.34952e03,
        0.36880e03,
        0.38817e03,
        0.40764e03,
        0.42723e03,
        0.44694e03,
        0.46679e03,
        0.48679e03,
        0.50696e03,
        0.52730e03,
        0.54783e03,
        0.56855e03,
        0.58948e03,
        0.61061e03,
        0.63198e03,
        0.65357e03,
        0.67539e03,
        0.69747e03,
        0.71979e03,
        0.74237e03,
        0.76521e03,
        0.78832e03,
        0.81169e03,
        0.83534e03,
        0.85927e03,
        0.88348e03,
        0.90798e03,
        0.93277e03,
        0.95784e03,
        0.98322e03,
        0.10089e04,
        0.10349e04,
        0.10611e04,
        0.10877e04,
        0.11146e04,
        0.11418e04,
        0.11693e04,
        0.11971e04,
        0.12253e04,
        0.12537e04,
        0.12825e04,
        0.13115e04,
        0.13409e04,
        0.13707e04,
        0.14007e04,
        0.14311e04,
        0.14617e04,
        0.14928e04,
        0.15241e04,
        0.15558e04,
        0.15877e04,
        0.16200e04,
        0.16527e04,
        0.16857e04,
        0.17190e04,
        0.17526e04,
        0.17866e04,
        0.18209e04,
        0.18555e04,
        0.18905e04,
        0.19258e04,
        0.19614e04,
        0.19974e04,
        0.20337e04,
        0.20703e04,
        0.21073e04,
        0.21446e04,
        0.21823e04,
        0.22203e04,
        0.22586e04,
        0.22973e04,
        0.23363e04,
        0.23756e04,
        0.24153e04,
        0.24553e04,
        0.24957e04,
        0.25364e04,
        0.25775e04,
        0.26189e04,
        0.26606e04,
        0.27027e04,
        0.27451e04,
        0.27879e04,
        0.28310e04,
        0.28745e04,
        0.29183e04,
        0.29625e04,
        0.30070e04,
        0.30518e04,
        0.30970e04,
        0.31425e04,
        0.31885e04,
        0.32347e04,
        0.32813e04,
        0.33282e04,
        0.33755e04,
        0.34231e04,
        0.34711e04,
        0.35194e04,
        0.35681e04,
        0.36172e04,
        0.36666e04,
        0.37163e04,
    ]
)


#  --------------- CO 28: M = 5, I = 3 ---------------------
M = 5
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.23030e02,
        0.32495e02,
        0.41966e02,
        0.51443e02,
        0.60926e02,
        0.70415e02,
        0.79910e02,
        0.89410e02,
        0.98918e02,
        0.10843e03,
        0.11795e03,
        0.12749e03,
        0.13703e03,
        0.14659e03,
        0.15618e03,
        0.16579e03,
        0.17543e03,
        0.18511e03,
        0.19483e03,
        0.20461e03,
        0.21444e03,
        0.22434e03,
        0.23430e03,
        0.24435e03,
        0.25447e03,
        0.26468e03,
        0.27499e03,
        0.28540e03,
        0.29591e03,
        0.30652e03,
        0.31725e03,
        0.32810e03,
        0.33906e03,
        0.35014e03,
        0.36136e03,
        0.37270e03,
        0.38417e03,
        0.39577e03,
        0.40752e03,
        0.41940e03,
        0.43142e03,
        0.44358e03,
        0.45589e03,
        0.46834e03,
        0.48094e03,
        0.49369e03,
        0.50659e03,
        0.51964e03,
        0.53284e03,
        0.54619e03,
        0.55971e03,
        0.57337e03,
        0.58719e03,
        0.60117e03,
        0.61530e03,
        0.62959e03,
        0.64405e03,
        0.65866e03,
        0.67343e03,
        0.68837e03,
        0.70346e03,
        0.71872e03,
        0.73414e03,
        0.74972e03,
        0.76547e03,
        0.78138e03,
        0.79745e03,
        0.81369e03,
        0.83010e03,
        0.84667e03,
        0.86341e03,
        0.88031e03,
        0.89738e03,
        0.91462e03,
        0.93202e03,
        0.94960e03,
        0.96734e03,
        0.98524e03,
        0.10033e04,
        0.10216e04,
        0.10400e04,
        0.10586e04,
        0.10773e04,
        0.10962e04,
        0.11153e04,
        0.11346e04,
        0.11540e04,
        0.11737e04,
        0.11934e04,
        0.12134e04,
        0.12335e04,
        0.12538e04,
        0.12743e04,
        0.12949e04,
        0.13157e04,
        0.13367e04,
        0.13578e04,
        0.13792e04,
        0.14007e04,
        0.14223e04,
        0.14442e04,
        0.14662e04,
        0.14884e04,
        0.15108e04,
        0.15333e04,
        0.15560e04,
        0.15789e04,
        0.16020e04,
        0.16252e04,
        0.16486e04,
        0.16722e04,
        0.16960e04,
        0.17199e04,
        0.17441e04,
        0.17684e04,
        0.17928e04,
        0.18175e04,
        0.18423e04,
        0.18673e04,
    ]
)


#  --------------- CO 27: M = 5, I = 4 ---------------------
M = 5
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.13505e03,
        0.19054e03,
        0.24606e03,
        0.30161e03,
        0.35720e03,
        0.41283e03,
        0.46848e03,
        0.52418e03,
        0.57991e03,
        0.63568e03,
        0.69149e03,
        0.74737e03,
        0.80332e03,
        0.85937e03,
        0.91553e03,
        0.97183e03,
        0.10283e04,
        0.10850e04,
        0.11420e04,
        0.11992e04,
        0.12568e04,
        0.13147e04,
        0.13730e04,
        0.14318e04,
        0.14910e04,
        0.15507e04,
        0.16110e04,
        0.16718e04,
        0.17332e04,
        0.17952e04,
        0.18579e04,
        0.19212e04,
        0.19852e04,
        0.20499e04,
        0.21153e04,
        0.21815e04,
        0.22484e04,
        0.23161e04,
        0.23846e04,
        0.24539e04,
        0.25240e04,
        0.25949e04,
        0.26666e04,
        0.27392e04,
        0.28127e04,
        0.28869e04,
        0.29621e04,
        0.30381e04,
        0.31150e04,
        0.31928e04,
        0.32715e04,
        0.33511e04,
        0.34316e04,
        0.35129e04,
        0.35952e04,
        0.36785e04,
        0.37626e04,
        0.38477e04,
        0.39336e04,
        0.40206e04,
        0.41084e04,
        0.41972e04,
        0.42869e04,
        0.43776e04,
        0.44692e04,
        0.45618e04,
        0.46553e04,
        0.47498e04,
        0.48452e04,
        0.49416e04,
        0.50390e04,
        0.51373e04,
        0.52366e04,
        0.53368e04,
        0.54381e04,
        0.55403e04,
        0.56435e04,
        0.57476e04,
        0.58527e04,
        0.59588e04,
        0.60659e04,
        0.61739e04,
        0.62829e04,
        0.63930e04,
        0.65040e04,
        0.66160e04,
        0.67290e04,
        0.68429e04,
        0.69579e04,
        0.70739e04,
        0.71908e04,
        0.73088e04,
        0.74277e04,
        0.75477e04,
        0.76686e04,
        0.77905e04,
        0.79135e04,
        0.80374e04,
        0.81624e04,
        0.82883e04,
        0.84153e04,
        0.85432e04,
        0.86722e04,
        0.88022e04,
        0.89331e04,
        0.90651e04,
        0.91982e04,
        0.93322e04,
        0.94672e04,
        0.96033e04,
        0.97404e04,
        0.98785e04,
        0.10018e05,
        0.10158e05,
        0.10299e05,
        0.10441e05,
        0.10584e05,
        0.10728e05,
        0.10874e05,
    ]
)


#  --------------- CO 38: M = 5, I = 5 ---------------------
M = 5
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.48264e02,
        0.68112e02,
        0.87974e02,
        0.10785e03,
        0.12773e03,
        0.14763e03,
        0.16754e03,
        0.18747e03,
        0.20741e03,
        0.22736e03,
        0.24733e03,
        0.26732e03,
        0.28735e03,
        0.30741e03,
        0.32752e03,
        0.34770e03,
        0.36794e03,
        0.38828e03,
        0.40871e03,
        0.42926e03,
        0.44994e03,
        0.47077e03,
        0.49175e03,
        0.51290e03,
        0.53424e03,
        0.55578e03,
        0.57752e03,
        0.59948e03,
        0.62166e03,
        0.64409e03,
        0.66676e03,
        0.68969e03,
        0.71287e03,
        0.73633e03,
        0.76006e03,
        0.78407e03,
        0.80836e03,
        0.83295e03,
        0.85784e03,
        0.88302e03,
        0.90851e03,
        0.93431e03,
        0.96042e03,
        0.98686e03,
        0.10136e04,
        0.10407e04,
        0.10681e04,
        0.10958e04,
        0.11238e04,
        0.11522e04,
        0.11809e04,
        0.12100e04,
        0.12393e04,
        0.12691e04,
        0.12991e04,
        0.13295e04,
        0.13603e04,
        0.13914e04,
        0.14228e04,
        0.14546e04,
        0.14867e04,
        0.15192e04,
        0.15520e04,
        0.15852e04,
        0.16187e04,
        0.16526e04,
        0.16869e04,
        0.17215e04,
        0.17564e04,
        0.17917e04,
        0.18274e04,
        0.18634e04,
        0.18998e04,
        0.19365e04,
        0.19736e04,
        0.20111e04,
        0.20489e04,
        0.20871e04,
        0.21256e04,
        0.21645e04,
        0.22038e04,
        0.22434e04,
        0.22834e04,
        0.23238e04,
        0.23645e04,
        0.24056e04,
        0.24471e04,
        0.24889e04,
        0.25311e04,
        0.25736e04,
        0.26166e04,
        0.26599e04,
        0.27035e04,
        0.27476e04,
        0.27920e04,
        0.28368e04,
        0.28819e04,
        0.29275e04,
        0.29733e04,
        0.30196e04,
        0.30662e04,
        0.31133e04,
        0.31606e04,
        0.32084e04,
        0.32565e04,
        0.33050e04,
        0.33539e04,
        0.34032e04,
        0.34528e04,
        0.35028e04,
        0.35532e04,
        0.36040e04,
        0.36551e04,
        0.37067e04,
        0.37586e04,
        0.38108e04,
        0.38635e04,
        0.39165e04,
        0.39699e04,
    ]
)


#  --------------- CO 37: M = 5, I = 6 ---------------------
M = 5
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.28271e03,
        0.39894e03,
        0.51524e03,
        0.63162e03,
        0.74807e03,
        0.86459e03,
        0.98119e03,
        0.10979e04,
        0.12146e04,
        0.13314e04,
        0.14484e04,
        0.15654e04,
        0.16826e04,
        0.18000e04,
        0.19176e04,
        0.20355e04,
        0.21538e04,
        0.22725e04,
        0.23916e04,
        0.25114e04,
        0.26318e04,
        0.27529e04,
        0.28749e04,
        0.29977e04,
        0.31215e04,
        0.32463e04,
        0.33721e04,
        0.34991e04,
        0.36274e04,
        0.37568e04,
        0.38876e04,
        0.40197e04,
        0.41533e04,
        0.42882e04,
        0.44247e04,
        0.45626e04,
        0.47022e04,
        0.48433e04,
        0.49860e04,
        0.51304e04,
        0.52763e04,
        0.54240e04,
        0.55735e04,
        0.57246e04,
        0.58775e04,
        0.60321e04,
        0.61886e04,
        0.63468e04,
        0.65068e04,
        0.66687e04,
        0.68324e04,
        0.69980e04,
        0.71654e04,
        0.73347e04,
        0.75058e04,
        0.76789e04,
        0.78539e04,
        0.80307e04,
        0.82096e04,
        0.83903e04,
        0.85729e04,
        0.87576e04,
        0.89441e04,
        0.91326e04,
        0.93230e04,
        0.95154e04,
        0.97098e04,
        0.99061e04,
        0.10104e05,
        0.10305e05,
        0.10507e05,
        0.10711e05,
        0.10918e05,
        0.11126e05,
        0.11336e05,
        0.11549e05,
        0.11763e05,
        0.11979e05,
        0.12198e05,
        0.12418e05,
        0.12640e05,
        0.12865e05,
        0.13091e05,
        0.13320e05,
        0.13550e05,
        0.13783e05,
        0.14018e05,
        0.14254e05,
        0.14493e05,
        0.14734e05,
        0.14977e05,
        0.15221e05,
        0.15468e05,
        0.15718e05,
        0.15969e05,
        0.16222e05,
        0.16477e05,
        0.16734e05,
        0.16994e05,
        0.17255e05,
        0.17519e05,
        0.17784e05,
        0.18052e05,
        0.18322e05,
        0.18594e05,
        0.18868e05,
        0.19144e05,
        0.19422e05,
        0.19703e05,
        0.19985e05,
        0.20270e05,
        0.20556e05,
        0.20845e05,
        0.21136e05,
        0.21429e05,
        0.21724e05,
        0.22021e05,
        0.22320e05,
        0.22622e05,
    ]
)


#  --------------- CH4 211: M = 6, I = 1 ---------------------
M = 6
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.54800e02,
        0.91500e02,
        0.13410e03,
        0.18180e03,
        0.23410e03,
        0.29070e03,
        0.35140e03,
        0.41600e03,
        0.48450e03,
        0.55720e03,
        0.63420e03,
        0.71600e03,
        0.80310e03,
        0.89590e03,
        0.99520e03,
        0.11017e04,
        0.12161e04,
        0.13393e04,
        0.14721e04,
        0.16155e04,
        0.17706e04,
        0.19384e04,
        0.21202e04,
        0.23172e04,
        0.25307e04,
        0.27624e04,
        0.30137e04,
        0.32864e04,
        0.35823e04,
        0.39034e04,
        0.42519e04,
        0.46300e04,
        0.50402e04,
        0.54853e04,
        0.59679e04,
        0.64913e04,
        0.70588e04,
        0.76739e04,
        0.83404e04,
        0.90625e04,
        0.98446e04,
        0.10691e05,
        0.11608e05,
        0.12600e05,
        0.13674e05,
        0.14835e05,
        0.16090e05,
        0.17447e05,
        0.18914e05,
        0.20500e05,
        0.22212e05,
        0.24063e05,
        0.26061e05,
        0.28218e05,
        0.30548e05,
        0.33063e05,
        0.35778e05,
        0.38708e05,
        0.41871e05,
        0.45284e05,
        0.48970e05,
        0.52940e05,
        0.57230e05,
        0.61860e05,
        0.66860e05,
        0.72250e05,
        0.78070e05,
        0.84350e05,
        0.91130e05,
        0.98450e05,
        0.10635e06,
        0.11488e06,
        0.12408e06,
        0.13403e06,
        0.14480e06,
        0.15640e06,
        0.16890e06,
        0.18240e06,
        0.19700e06,
        0.21280e06,
        0.22980e06,
        0.24830e06,
        0.26820e06,
        0.28970e06,
        0.31290e06,
        0.33800e06,
        0.36520e06,
        0.39450e06,
        0.42600e06,
        0.46000e06,
        0.49700e06,
        0.53700e06,
        0.58100e06,
        0.62700e06,
        0.67800e06,
        0.73300e06,
        0.79200e06,
        0.85600e06,
        0.92500e06,
        0.10000e07,
        0.10800e07,
        0.11670e07,
        0.12610e07,
        0.13620e07,
        0.14720e07,
        0.15910e07,
        0.17190e07,
        0.18600e07,
        0.20100e07,
        0.21700e07,
        0.23400e07,
        0.25300e07,
        0.27300e07,
        0.29500e07,
        0.31800e07,
        0.34300e07,
        0.37000e07,
        0.39900e07,
        0.42856e07,
    ]
)


#  --------------- CH4 311: M = 6, I = 2 ---------------------
M = 6
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.10958e03,
        0.18304e03,
        0.26818e03,
        0.36356e03,
        0.46820e03,
        0.58141e03,
        0.70270e03,
        0.83186e03,
        0.96893e03,
        0.11142e04,
        0.12682e04,
        0.14316e04,
        0.16055e04,
        0.17909e04,
        0.19891e04,
        0.22016e04,
        0.24297e04,
        0.26752e04,
        0.29399e04,
        0.32255e04,
        0.35342e04,
        0.38680e04,
        0.42294e04,
        0.46208e04,
        0.50449e04,
        0.55046e04,
        0.60030e04,
        0.65434e04,
        0.71293e04,
        0.77646e04,
        0.84535e04,
        0.92004e04,
        0.10010e05,
        0.10888e05,
        0.11838e05,
        0.12869e05,
        0.13984e05,
        0.15193e05,
        0.16501e05,
        0.17916e05,
        0.19448e05,
        0.21104e05,
        0.22895e05,
        0.24830e05,
        0.26921e05,
        0.29180e05,
        0.31618e05,
        0.34250e05,
        0.37090e05,
        0.40152e05,
        0.43454e05,
        0.47012e05,
        0.50845e05,
        0.54973e05,
        0.59416e05,
        0.64197e05,
        0.69340e05,
        0.74870e05,
        0.80813e05,
        0.87198e05,
        0.94055e05,
        0.10142e06,
        0.10932e06,
        0.11779e06,
        0.12688e06,
        0.13662e06,
        0.14706e06,
        0.15824e06,
        0.17021e06,
        0.18302e06,
        0.19673e06,
        0.21139e06,
        0.22706e06,
        0.24381e06,
        0.26171e06,
        0.28082e06,
        0.30122e06,
        0.32299e06,
        0.34621e06,
        0.37097e06,
        0.39737e06,
        0.42551e06,
        0.45548e06,
        0.48739e06,
        0.52136e06,
        0.55752e06,
        0.59598e06,
        0.63688e06,
        0.68036e06,
        0.72657e06,
        0.77566e06,
        0.82780e06,
        0.88316e06,
        0.94191e06,
        0.10043e07,
        0.10704e07,
        0.11405e07,
        0.12148e07,
        0.12936e07,
        0.13770e07,
        0.14654e07,
        0.15589e07,
        0.16579e07,
        0.17627e07,
        0.18736e07,
        0.19908e07,
        0.21147e07,
        0.22456e07,
        0.23840e07,
        0.25301e07,
        0.26844e07,
        0.28474e07,
        0.30193e07,
        0.32007e07,
        0.33921e07,
        0.35939e07,
        0.38067e07,
        0.40310e07,
        0.42673e07,
    ]
)


#  --------------- CH4 212: M = 6, I = 3 ---------------------
M = 6
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.44079e03,
        0.73786e03,
        0.10822e04,
        0.14679e04,
        0.18913e04,
        0.23497e04,
        0.28415e04,
        0.33665e04,
        0.39257e04,
        0.45211e04,
        0.51562e04,
        0.58349e04,
        0.65624e04,
        0.73445e04,
        0.81872e04,
        0.90978e04,
        0.10084e05,
        0.11153e05,
        0.12315e05,
        0.13579e05,
        0.14955e05,
        0.16455e05,
        0.18089e05,
        0.19871e05,
        0.21816e05,
        0.23937e05,
        0.26251e05,
        0.28776e05,
        0.31531e05,
        0.34535e05,
        0.37811e05,
        0.41384e05,
        0.45278e05,
        0.49521e05,
        0.54144e05,
        0.59178e05,
        0.64657e05,
        0.70621e05,
        0.77108e05,
        0.84161e05,
        0.91828e05,
        0.10016e06,
        0.10921e06,
        0.11903e06,
        0.12968e06,
        0.14124e06,
        0.15378e06,
        0.16736e06,
        0.18207e06,
        0.19800e06,
        0.21524e06,
        0.23389e06,
        0.25405e06,
        0.27585e06,
        0.29939e06,
        0.32482e06,
        0.35226e06,
        0.38186e06,
        0.41379e06,
        0.44821e06,
        0.48529e06,
        0.52522e06,
        0.56821e06,
        0.61447e06,
        0.66422e06,
        0.71771e06,
        0.77519e06,
        0.83693e06,
        0.90323e06,
        0.97438e06,
        0.10507e07,
        0.11326e07,
        0.12203e07,
        0.13143e07,
        0.14150e07,
        0.15228e07,
        0.16382e07,
        0.17616e07,
        0.18935e07,
        0.20346e07,
        0.21853e07,
        0.23463e07,
        0.25181e07,
        0.27016e07,
        0.28973e07,
        0.31060e07,
        0.33284e07,
        0.35655e07,
        0.38181e07,
        0.40870e07,
        0.43733e07,
        0.46780e07,
        0.50020e07,
        0.53467e07,
        0.57130e07,
        0.61023e07,
        0.65158e07,
        0.69549e07,
        0.74211e07,
        0.79158e07,
        0.84407e07,
        0.89973e07,
        0.95874e07,
        0.10213e08,
        0.10875e08,
        0.11577e08,
        0.12320e08,
        0.13107e08,
        0.13940e08,
        0.14820e08,
        0.15752e08,
        0.16736e08,
        0.17777e08,
        0.18877e08,
        0.20038e08,
        0.21265e08,
        0.22560e08,
        0.23927e08,
        0.25369e08,
    ]
)


#  --------------- CH4 312: M = 6, I = 4 ---------------------
M = 6
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.88231e03,
        0.14770e04,
        0.21661e04,
        0.29384e04,
        0.37859e04,
        0.47034e04,
        0.56879e04,
        0.67388e04,
        0.78581e04,
        0.90501e04,
        0.10321e05,
        0.11680e05,
        0.13136e05,
        0.14702e05,
        0.16389e05,
        0.18212e05,
        0.20186e05,
        0.22328e05,
        0.24654e05,
        0.27185e05,
        0.29941e05,
        0.32943e05,
        0.36216e05,
        0.39786e05,
        0.43681e05,
        0.47930e05,
        0.52567e05,
        0.57625e05,
        0.63144e05,
        0.69164e05,
        0.75730e05,
        0.82890e05,
        0.90693e05,
        0.99198e05,
        0.10846e06,
        0.11855e06,
        0.12954e06,
        0.14149e06,
        0.15450e06,
        0.16864e06,
        0.18402e06,
        0.20072e06,
        0.21886e06,
        0.23856e06,
        0.25993e06,
        0.28312e06,
        0.30825e06,
        0.33550e06,
        0.36501e06,
        0.39696e06,
        0.43155e06,
        0.46896e06,
        0.50942e06,
        0.55315e06,
        0.60039e06,
        0.65141e06,
        0.70648e06,
        0.76589e06,
        0.82997e06,
        0.89904e06,
        0.97346e06,
        0.10536e07,
        0.11399e07,
        0.12327e07,
        0.13326e07,
        0.14400e07,
        0.15554e07,
        0.16793e07,
        0.18124e07,
        0.19553e07,
        0.21085e07,
        0.22729e07,
        0.24490e07,
        0.26378e07,
        0.28400e07,
        0.30565e07,
        0.32881e07,
        0.35360e07,
        0.38010e07,
        0.40843e07,
        0.43870e07,
        0.47103e07,
        0.50555e07,
        0.54239e07,
        0.58169e07,
        0.62361e07,
        0.66830e07,
        0.71592e07,
        0.76666e07,
        0.82069e07,
        0.87820e07,
        0.93940e07,
        0.10045e08,
        0.10737e08,
        0.11473e08,
        0.12256e08,
        0.13086e08,
        0.13969e08,
        0.14905e08,
        0.15899e08,
        0.16954e08,
        0.18072e08,
        0.19258e08,
        0.20515e08,
        0.21847e08,
        0.23257e08,
        0.24750e08,
        0.26331e08,
        0.28004e08,
        0.29774e08,
        0.31646e08,
        0.33625e08,
        0.35716e08,
        0.37926e08,
        0.40261e08,
        0.42726e08,
        0.45329e08,
        0.48077e08,
        0.50975e08,
    ]
)


#  --------------- O2 66: M = 7, I = 1 ---------------------
M = 7
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.44334e02,
        0.62460e02,
        0.80596e02,
        0.98738e02,
        0.11688e03,
        0.13503e03,
        0.15319e03,
        0.17136e03,
        0.18954e03,
        0.20775e03,
        0.22600e03,
        0.24431e03,
        0.26270e03,
        0.28119e03,
        0.29981e03,
        0.31857e03,
        0.33750e03,
        0.35662e03,
        0.37594e03,
        0.39550e03,
        0.41529e03,
        0.43535e03,
        0.45568e03,
        0.47630e03,
        0.49722e03,
        0.51844e03,
        0.53998e03,
        0.56185e03,
        0.58406e03,
        0.60660e03,
        0.62949e03,
        0.65274e03,
        0.67635e03,
        0.70031e03,
        0.72465e03,
        0.74936e03,
        0.77444e03,
        0.79990e03,
        0.82574e03,
        0.85197e03,
        0.87858e03,
        0.90558e03,
        0.93297e03,
        0.96076e03,
        0.98895e03,
        0.10175e04,
        0.10465e04,
        0.10759e04,
        0.11057e04,
        0.11359e04,
        0.11665e04,
        0.11976e04,
        0.12290e04,
        0.12609e04,
        0.12931e04,
        0.13258e04,
        0.13590e04,
        0.13925e04,
        0.14265e04,
        0.14609e04,
        0.14958e04,
        0.15311e04,
        0.15669e04,
        0.16031e04,
        0.16397e04,
        0.16768e04,
        0.17144e04,
        0.17524e04,
        0.17909e04,
        0.18298e04,
        0.18692e04,
        0.19091e04,
        0.19495e04,
        0.19904e04,
        0.20318e04,
        0.20736e04,
        0.21160e04,
        0.21588e04,
        0.22022e04,
        0.22461e04,
        0.22905e04,
        0.23354e04,
        0.23809e04,
        0.24268e04,
        0.24734e04,
        0.25204e04,
        0.25680e04,
        0.26162e04,
        0.26649e04,
        0.27142e04,
        0.27641e04,
        0.28145e04,
        0.28655e04,
        0.29171e04,
        0.29693e04,
        0.30221e04,
        0.30755e04,
        0.31295e04,
        0.31841e04,
        0.32393e04,
        0.32951e04,
        0.33516e04,
        0.34087e04,
        0.34665e04,
        0.35249e04,
        0.35839e04,
        0.36436e04,
        0.37040e04,
        0.37650e04,
        0.38267e04,
        0.38891e04,
        0.39522e04,
        0.40159e04,
        0.40804e04,
        0.41455e04,
        0.42114e04,
        0.42780e04,
        0.43452e04,
        0.44132e04,
    ]
)


#  --------------- O2 68: M = 7, I = 2 ---------------------
M = 7
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.89206e02,
        0.12759e03,
        0.16600e03,
        0.20442e03,
        0.24285e03,
        0.28128e03,
        0.31973e03,
        0.35821e03,
        0.39672e03,
        0.43530e03,
        0.47398e03,
        0.51281e03,
        0.55183e03,
        0.59108e03,
        0.63062e03,
        0.67051e03,
        0.71078e03,
        0.75148e03,
        0.79265e03,
        0.83435e03,
        0.87659e03,
        0.91941e03,
        0.96285e03,
        0.10069e04,
        0.10517e04,
        0.10971e04,
        0.11432e04,
        0.11901e04,
        0.12377e04,
        0.12861e04,
        0.13352e04,
        0.13851e04,
        0.14358e04,
        0.14872e04,
        0.15395e04,
        0.15926e04,
        0.16466e04,
        0.17013e04,
        0.17569e04,
        0.18134e04,
        0.18706e04,
        0.19288e04,
        0.19877e04,
        0.20476e04,
        0.21083e04,
        0.21698e04,
        0.22323e04,
        0.22956e04,
        0.23598e04,
        0.24248e04,
        0.24908e04,
        0.25576e04,
        0.26253e04,
        0.26940e04,
        0.27635e04,
        0.28339e04,
        0.29052e04,
        0.29775e04,
        0.30506e04,
        0.31247e04,
        0.31997e04,
        0.32756e04,
        0.33524e04,
        0.34302e04,
        0.35089e04,
        0.35885e04,
        0.36691e04,
        0.37506e04,
        0.38331e04,
        0.39166e04,
        0.40010e04,
        0.40864e04,
        0.41727e04,
        0.42601e04,
        0.43484e04,
        0.44377e04,
        0.45280e04,
        0.46193e04,
        0.47116e04,
        0.48049e04,
        0.48992e04,
        0.49946e04,
        0.50909e04,
        0.51883e04,
        0.52868e04,
        0.53863e04,
        0.54868e04,
        0.55884e04,
        0.56911e04,
        0.57949e04,
        0.58997e04,
        0.60056e04,
        0.61126e04,
        0.62207e04,
        0.63298e04,
        0.64401e04,
        0.65516e04,
        0.66641e04,
        0.67778e04,
        0.68926e04,
        0.70085e04,
        0.71256e04,
        0.72439e04,
        0.73633e04,
        0.74839e04,
        0.76056e04,
        0.77286e04,
        0.78527e04,
        0.79781e04,
        0.81046e04,
        0.82324e04,
        0.83613e04,
        0.84915e04,
        0.86229e04,
        0.87556e04,
        0.88895e04,
        0.90247e04,
        0.91611e04,
        0.92988e04,
    ]
)


#  --------------- O2 67: M = 7, I = 3 ---------------------
M = 7
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.52071e03,
        0.74484e03,
        0.96908e03,
        0.11934e04,
        0.14177e04,
        0.16422e04,
        0.18667e04,
        0.20913e04,
        0.23161e04,
        0.25413e04,
        0.27671e04,
        0.29936e04,
        0.32212e04,
        0.34501e04,
        0.36806e04,
        0.39130e04,
        0.41476e04,
        0.43846e04,
        0.46242e04,
        0.48668e04,
        0.51125e04,
        0.53615e04,
        0.56140e04,
        0.58701e04,
        0.61300e04,
        0.63938e04,
        0.66617e04,
        0.69337e04,
        0.72099e04,
        0.74904e04,
        0.77754e04,
        0.80647e04,
        0.83586e04,
        0.86571e04,
        0.89602e04,
        0.92680e04,
        0.95805e04,
        0.98977e04,
        0.10220e05,
        0.10547e05,
        0.10878e05,
        0.11215e05,
        0.11556e05,
        0.11903e05,
        0.12254e05,
        0.12611e05,
        0.12972e05,
        0.13338e05,
        0.13710e05,
        0.14086e05,
        0.14468e05,
        0.14855e05,
        0.15247e05,
        0.15644e05,
        0.16046e05,
        0.16453e05,
        0.16866e05,
        0.17283e05,
        0.17706e05,
        0.18135e05,
        0.18568e05,
        0.19007e05,
        0.19452e05,
        0.19901e05,
        0.20356e05,
        0.20817e05,
        0.21283e05,
        0.21754e05,
        0.22231e05,
        0.22713e05,
        0.23201e05,
        0.23695e05,
        0.24194e05,
        0.24699e05,
        0.25209e05,
        0.25725e05,
        0.26247e05,
        0.26775e05,
        0.27308e05,
        0.27847e05,
        0.28393e05,
        0.28944e05,
        0.29500e05,
        0.30063e05,
        0.30632e05,
        0.31207e05,
        0.31788e05,
        0.32375e05,
        0.32968e05,
        0.33568e05,
        0.34173e05,
        0.34785e05,
        0.35403e05,
        0.36028e05,
        0.36659e05,
        0.37296e05,
        0.37939e05,
        0.38590e05,
        0.39246e05,
        0.39909e05,
        0.40579e05,
        0.41256e05,
        0.41939e05,
        0.42629e05,
        0.43325e05,
        0.44029e05,
        0.44739e05,
        0.45456e05,
        0.46180e05,
        0.46911e05,
        0.47649e05,
        0.48394e05,
        0.49146e05,
        0.49905e05,
        0.50671e05,
        0.51445e05,
        0.52226e05,
        0.53014e05,
        0.53809e05,
    ]
)


#  --------------- NO 46: M = 8, I = 1 ---------------------
M = 8
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.15840e03,
        0.23971e03,
        0.33080e03,
        0.42907e03,
        0.53251e03,
        0.63972e03,
        0.74975e03,
        0.86195e03,
        0.97582e03,
        0.10911e04,
        0.12074e04,
        0.13248e04,
        0.14430e04,
        0.15621e04,
        0.16820e04,
        0.18027e04,
        0.19243e04,
        0.20468e04,
        0.21703e04,
        0.22948e04,
        0.24204e04,
        0.25472e04,
        0.26753e04,
        0.28046e04,
        0.29354e04,
        0.30676e04,
        0.32013e04,
        0.33365e04,
        0.34734e04,
        0.36120e04,
        0.37522e04,
        0.38942e04,
        0.40379e04,
        0.41835e04,
        0.43310e04,
        0.44803e04,
        0.46316e04,
        0.47849e04,
        0.49400e04,
        0.50972e04,
        0.52564e04,
        0.54176e04,
        0.55809e04,
        0.57462e04,
        0.59137e04,
        0.60832e04,
        0.62548e04,
        0.64286e04,
        0.66045e04,
        0.67825e04,
        0.69628e04,
        0.71451e04,
        0.73297e04,
        0.75164e04,
        0.77053e04,
        0.78964e04,
        0.80897e04,
        0.82853e04,
        0.84830e04,
        0.86830e04,
        0.88852e04,
        0.90896e04,
        0.92963e04,
        0.95052e04,
        0.97164e04,
        0.99297e04,
        0.10145e05,
        0.10363e05,
        0.10583e05,
        0.10806e05,
        0.11031e05,
        0.11258e05,
        0.11487e05,
        0.11718e05,
        0.11952e05,
        0.12188e05,
        0.12426e05,
        0.12667e05,
        0.12910e05,
        0.13155e05,
        0.13403e05,
        0.13652e05,
        0.13905e05,
        0.14159e05,
        0.14416e05,
        0.14675e05,
        0.14936e05,
        0.15199e05,
        0.15465e05,
        0.15733e05,
        0.16004e05,
        0.16277e05,
        0.16552e05,
        0.16829e05,
        0.17109e05,
        0.17391e05,
        0.17675e05,
        0.17962e05,
        0.18251e05,
        0.18542e05,
        0.18836e05,
        0.19131e05,
        0.19430e05,
        0.19730e05,
        0.20033e05,
        0.20338e05,
        0.20646e05,
        0.20955e05,
        0.21268e05,
        0.21582e05,
        0.21899e05,
        0.22218e05,
        0.22539e05,
        0.22863e05,
        0.23189e05,
        0.23518e05,
        0.23848e05,
        0.24181e05,
        0.24517e05,
    ]
)


#  --------------- NO 56: M = 8, I = 2 ---------------------
M = 8
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.10942e03,
        0.16560e03,
        0.22856e03,
        0.29647e03,
        0.36795e03,
        0.44204e03,
        0.51808e03,
        0.59561e03,
        0.67432e03,
        0.75396e03,
        0.83439e03,
        0.91551e03,
        0.99725e03,
        0.10796e04,
        0.11625e04,
        0.12460e04,
        0.13302e04,
        0.14150e04,
        0.15005e04,
        0.15868e04,
        0.16739e04,
        0.17618e04,
        0.18506e04,
        0.19404e04,
        0.20311e04,
        0.21229e04,
        0.22158e04,
        0.23098e04,
        0.24050e04,
        0.25013e04,
        0.25989e04,
        0.26976e04,
        0.27977e04,
        0.28991e04,
        0.30018e04,
        0.31058e04,
        0.32112e04,
        0.33180e04,
        0.34262e04,
        0.35358e04,
        0.36468e04,
        0.37593e04,
        0.38732e04,
        0.39885e04,
        0.41054e04,
        0.42237e04,
        0.43436e04,
        0.44649e04,
        0.45877e04,
        0.47121e04,
        0.48379e04,
        0.49654e04,
        0.50943e04,
        0.52248e04,
        0.53568e04,
        0.54904e04,
        0.56255e04,
        0.57622e04,
        0.59004e04,
        0.60403e04,
        0.61816e04,
        0.63246e04,
        0.64692e04,
        0.66152e04,
        0.67630e04,
        0.69123e04,
        0.70631e04,
        0.72156e04,
        0.73696e04,
        0.75253e04,
        0.76825e04,
        0.78414e04,
        0.80018e04,
        0.81638e04,
        0.83275e04,
        0.84927e04,
        0.86596e04,
        0.88280e04,
        0.89981e04,
        0.91698e04,
        0.93430e04,
        0.95180e04,
        0.96945e04,
        0.98726e04,
        0.10052e05,
        0.10234e05,
        0.10417e05,
        0.10601e05,
        0.10788e05,
        0.10975e05,
        0.11165e05,
        0.11356e05,
        0.11549e05,
        0.11743e05,
        0.11939e05,
        0.12137e05,
        0.12336e05,
        0.12537e05,
        0.12739e05,
        0.12943e05,
        0.13149e05,
        0.13356e05,
        0.13565e05,
        0.13776e05,
        0.13988e05,
        0.14202e05,
        0.14418e05,
        0.14635e05,
        0.14853e05,
        0.15074e05,
        0.15296e05,
        0.15520e05,
        0.15745e05,
        0.15972e05,
        0.16200e05,
        0.16431e05,
        0.16663e05,
        0.16896e05,
        0.17131e05,
    ]
)


#  --------------- NO 48: M = 8, I = 3 ---------------------
M = 8
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.16695e03,
        0.25269e03,
        0.34876e03,
        0.45239e03,
        0.56148e03,
        0.67455e03,
        0.79059e03,
        0.90891e03,
        0.10290e04,
        0.11506e04,
        0.12733e04,
        0.13971e04,
        0.15219e04,
        0.16476e04,
        0.17742e04,
        0.19017e04,
        0.20302e04,
        0.21598e04,
        0.22904e04,
        0.24223e04,
        0.25553e04,
        0.26897e04,
        0.28255e04,
        0.29628e04,
        0.31016e04,
        0.32420e04,
        0.33842e04,
        0.35280e04,
        0.36736e04,
        0.38211e04,
        0.39704e04,
        0.41217e04,
        0.42750e04,
        0.44302e04,
        0.45876e04,
        0.47469e04,
        0.49084e04,
        0.50720e04,
        0.52378e04,
        0.54058e04,
        0.55759e04,
        0.57483e04,
        0.59230e04,
        0.60999e04,
        0.62791e04,
        0.64605e04,
        0.66443e04,
        0.68304e04,
        0.70187e04,
        0.72095e04,
        0.74026e04,
        0.75980e04,
        0.77958e04,
        0.79960e04,
        0.81986e04,
        0.84036e04,
        0.86109e04,
        0.88207e04,
        0.90328e04,
        0.92474e04,
        0.94644e04,
        0.96839e04,
        0.99057e04,
        0.10130e05,
        0.10357e05,
        0.10586e05,
        0.10817e05,
        0.11052e05,
        0.11288e05,
        0.11527e05,
        0.11768e05,
        0.12012e05,
        0.12259e05,
        0.12507e05,
        0.12759e05,
        0.13012e05,
        0.13269e05,
        0.13527e05,
        0.13788e05,
        0.14052e05,
        0.14318e05,
        0.14587e05,
        0.14858e05,
        0.15131e05,
        0.15408e05,
        0.15686e05,
        0.15967e05,
        0.16251e05,
        0.16537e05,
        0.16825e05,
        0.17116e05,
        0.17410e05,
        0.17706e05,
        0.18004e05,
        0.18305e05,
        0.18609e05,
        0.18915e05,
        0.19224e05,
        0.19535e05,
        0.19848e05,
        0.20164e05,
        0.20483e05,
        0.20804e05,
        0.21127e05,
        0.21453e05,
        0.21782e05,
        0.22113e05,
        0.22447e05,
        0.22783e05,
        0.23122e05,
        0.23463e05,
        0.23807e05,
        0.24153e05,
        0.24502e05,
        0.24853e05,
        0.25207e05,
        0.25563e05,
        0.25922e05,
        0.26283e05,
    ]
)


#  --------------- SO2 626: M = 9, I = 1 ---------------------
M = 9
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.52899e03,
        0.89171e03,
        0.13139e04,
        0.17915e04,
        0.23246e04,
        0.29155e04,
        0.35675e04,
        0.42848e04,
        0.50723e04,
        0.59352e04,
        0.68794e04,
        0.79109e04,
        0.90366e04,
        0.10264e05,
        0.11599e05,
        0.13052e05,
        0.14629e05,
        0.16340e05,
        0.18193e05,
        0.20199e05,
        0.22366e05,
        0.24704e05,
        0.27225e05,
        0.29938e05,
        0.32855e05,
        0.35987e05,
        0.39346e05,
        0.42944e05,
        0.46794e05,
        0.50909e05,
        0.55302e05,
        0.59986e05,
        0.64977e05,
        0.70288e05,
        0.75934e05,
        0.81931e05,
        0.88294e05,
        0.95040e05,
        0.10219e06,
        0.10975e06,
        0.11774e06,
        0.12619e06,
        0.13511e06,
        0.14452e06,
        0.15443e06,
        0.16487e06,
        0.17586e06,
        0.18742e06,
        0.19957e06,
        0.21234e06,
        0.22573e06,
        0.23978e06,
        0.25451e06,
        0.26995e06,
        0.28611e06,
        0.30302e06,
        0.32071e06,
        0.33920e06,
        0.35852e06,
        0.37869e06,
        0.39974e06,
        0.42171e06,
        0.44461e06,
        0.46848e06,
        0.49334e06,
        0.51922e06,
        0.54617e06,
        0.57419e06,
        0.60334e06,
        0.63363e06,
        0.66511e06,
        0.69780e06,
        0.73174e06,
        0.76696e06,
        0.80349e06,
        0.84138e06,
        0.88066e06,
        0.92136e06,
        0.96352e06,
        0.10072e07,
        0.10524e07,
        0.10992e07,
        0.11475e07,
        0.11976e07,
        0.12493e07,
        0.13028e07,
        0.13580e07,
        0.14151e07,
        0.14741e07,
        0.15349e07,
        0.15977e07,
        0.16625e07,
        0.17293e07,
        0.17982e07,
        0.18693e07,
        0.19425e07,
        0.20180e07,
        0.20958e07,
        0.21758e07,
        0.22583e07,
        0.23432e07,
        0.24305e07,
        0.25204e07,
        0.26129e07,
        0.27080e07,
        0.28058e07,
        0.29064e07,
        0.30097e07,
        0.31159e07,
        0.32250e07,
        0.33371e07,
        0.34522e07,
        0.35705e07,
        0.36918e07,
        0.38164e07,
        0.39442e07,
        0.40754e07,
        0.42099e07,
        0.43479e07,
    ]
)


#  --------------- SO2 646: M = 9, I = 2 ---------------------
M = 9
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.53140e03,
        0.89578e03,
        0.13199e04,
        0.17997e04,
        0.23353e04,
        0.29288e04,
        0.35837e04,
        0.43043e04,
        0.50953e04,
        0.59621e04,
        0.69104e04,
        0.79465e04,
        0.90772e04,
        0.10310e05,
        0.11651e05,
        0.13110e05,
        0.14694e05,
        0.16413e05,
        0.18274e05,
        0.20289e05,
        0.22465e05,
        0.24814e05,
        0.27345e05,
        0.30070e05,
        0.33000e05,
        0.36145e05,
        0.39519e05,
        0.43133e05,
        0.46999e05,
        0.51132e05,
        0.55544e05,
        0.60248e05,
        0.65260e05,
        0.70594e05,
        0.76264e05,
        0.82287e05,
        0.88678e05,
        0.95453e05,
        0.10263e06,
        0.11022e06,
        0.11825e06,
        0.12674e06,
        0.13569e06,
        0.14514e06,
        0.15510e06,
        0.16558e06,
        0.17662e06,
        0.18823e06,
        0.20043e06,
        0.21325e06,
        0.22670e06,
        0.24081e06,
        0.25561e06,
        0.27111e06,
        0.28733e06,
        0.30432e06,
        0.32208e06,
        0.34065e06,
        0.36005e06,
        0.38031e06,
        0.40145e06,
        0.42351e06,
        0.44651e06,
        0.47047e06,
        0.49544e06,
        0.52144e06,
        0.54849e06,
        0.57664e06,
        0.60591e06,
        0.63633e06,
        0.66794e06,
        0.70077e06,
        0.73485e06,
        0.77022e06,
        0.80691e06,
        0.84496e06,
        0.88440e06,
        0.92527e06,
        0.96761e06,
        0.10115e07,
        0.10568e07,
        0.11038e07,
        0.11524e07,
        0.12027e07,
        0.12546e07,
        0.13083e07,
        0.13638e07,
        0.14211e07,
        0.14803e07,
        0.15414e07,
        0.16045e07,
        0.16695e07,
        0.17366e07,
        0.18059e07,
        0.18772e07,
        0.19507e07,
        0.20265e07,
        0.21046e07,
        0.21850e07,
        0.22678e07,
        0.23531e07,
        0.24408e07,
        0.25310e07,
        0.26239e07,
        0.27194e07,
        0.28176e07,
        0.29186e07,
        0.30224e07,
        0.31290e07,
        0.32386e07,
        0.33512e07,
        0.34668e07,
        0.35855e07,
        0.37074e07,
        0.38324e07,
        0.39608e07,
        0.40925e07,
        0.42276e07,
        0.43662e07,
    ]
)

#  --------------- NO2 646: M = 10, I = 1 ---------------------
M = 10
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.12046e04,
        0.20297e04,
        0.29875e04,
        0.40626e04,
        0.52463e04,
        0.65350e04,
        0.79286e04,
        0.94298e04,
        0.11043e05,
        0.12776e05,
        0.14634e05,
        0.16627e05,
        0.18765e05,
        0.21056e05,
        0.23511e05,
        0.26143e05,
        0.28961e05,
        0.31979e05,
        0.35209e05,
        0.38663e05,
        0.42355e05,
        0.46300e05,
        0.50510e05,
        0.55001e05,
        0.59787e05,
        0.64884e05,
        0.70308e05,
        0.76075e05,
        0.82201e05,
        0.88704e05,
        0.95602e05,
        0.10291e06,
        0.11065e06,
        0.11884e06,
        0.12750e06,
        0.13665e06,
        0.14631e06,
        0.15650e06,
        0.16724e06,
        0.17856e06,
        0.19047e06,
        0.20301e06,
        0.21618e06,
        0.23002e06,
        0.24456e06,
        0.25981e06,
        0.27580e06,
        0.29256e06,
        0.31012e06,
        0.32850e06,
        0.34773e06,
        0.36784e06,
        0.38886e06,
        0.41082e06,
        0.43374e06,
        0.45766e06,
        0.48262e06,
        0.50863e06,
        0.53574e06,
        0.56398e06,
        0.59339e06,
        0.62398e06,
        0.65581e06,
        0.68891e06,
        0.72331e06,
        0.75905e06,
        0.79617e06,
        0.83470e06,
        0.87469e06,
        0.91617e06,
        0.95919e06,
        0.10038e07,
        0.10500e07,
        0.10979e07,
        0.11474e07,
        0.11988e07,
        0.12519e07,
        0.13068e07,
        0.13636e07,
        0.14224e07,
        0.14831e07,
        0.15459e07,
        0.16107e07,
        0.16776e07,
        0.17467e07,
        0.18180e07,
        0.18916e07,
        0.19675e07,
        0.20458e07,
        0.21265e07,
        0.22097e07,
        0.22954e07,
        0.23837e07,
        0.24747e07,
        0.25684e07,
        0.26648e07,
        0.27641e07,
        0.28662e07,
        0.29713e07,
        0.30794e07,
        0.31905e07,
        0.33048e07,
        0.34223e07,
        0.35430e07,
        0.36670e07,
        0.37944e07,
        0.39253e07,
        0.40597e07,
        0.41976e07,
        0.43393e07,
        0.44846e07,
        0.46337e07,
        0.47867e07,
        0.49437e07,
        0.51046e07,
        0.52696e07,
        0.54388e07,
        0.56122e07,
        0.57900e07,
    ]
)


#  --------------- NH3 4111: M = 11, I = 1 ---------------------
M = 11
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.16013e03,
        0.26692e03,
        0.39067e03,
        0.52933e03,
        0.68153e03,
        0.84641e03,
        0.10234e04,
        0.12125e04,
        0.14136e04,
        0.16272e04,
        0.18537e04,
        0.20937e04,
        0.23481e04,
        0.26177e04,
        0.29035e04,
        0.32065e04,
        0.35279e04,
        0.38688e04,
        0.42304e04,
        0.46141e04,
        0.50212e04,
        0.54531e04,
        0.59114e04,
        0.63976e04,
        0.69133e04,
        0.74602e04,
        0.80401e04,
        0.86549e04,
        0.93066e04,
        0.99971e04,
        0.10729e05,
        0.11504e05,
        0.12324e05,
        0.13193e05,
        0.14112e05,
        0.15085e05,
        0.16114e05,
        0.17201e05,
        0.18352e05,
        0.19567e05,
        0.20851e05,
        0.22208e05,
        0.23640e05,
        0.25152e05,
        0.26747e05,
        0.28430e05,
        0.30205e05,
        0.32077e05,
        0.34050e05,
        0.36128e05,
        0.38317e05,
        0.40623e05,
        0.43050e05,
        0.45605e05,
        0.48292e05,
        0.51119e05,
        0.54091e05,
        0.57215e05,
        0.60498e05,
        0.63947e05,
        0.67569e05,
        0.71372e05,
        0.75364e05,
        0.79552e05,
        0.83946e05,
        0.88553e05,
        0.93384e05,
        0.98447e05,
        0.10375e06,
        0.10931e06,
        0.11513e06,
        0.12122e06,
        0.12760e06,
        0.13427e06,
        0.14125e06,
        0.14855e06,
        0.15619e06,
        0.16417e06,
        0.17250e06,
        0.18121e06,
        0.19031e06,
        0.19981e06,
        0.20973e06,
        0.22008e06,
        0.23088e06,
        0.24215e06,
        0.25390e06,
        0.26615e06,
        0.27892e06,
        0.29223e06,
        0.30610e06,
        0.32055e06,
        0.33559e06,
        0.35125e06,
        0.36756e06,
        0.38453e06,
        0.40219e06,
        0.42056e06,
        0.43967e06,
        0.45953e06,
        0.48019e06,
        0.50165e06,
        0.52396e06,
        0.54714e06,
        0.57122e06,
        0.59622e06,
        0.62218e06,
        0.64913e06,
        0.67710e06,
        0.70613e06,
        0.73624e06,
        0.76748e06,
        0.79988e06,
        0.83347e06,
        0.86829e06,
        0.90439e06,
        0.94180e06,
        0.98056e06,
        0.10207e07,
    ]
)


#  --------------- NH3 5111: M = 11, I = 2 ---------------------
M = 11
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.10697e03,
        0.17832e03,
        0.26100e03,
        0.35364e03,
        0.45533e03,
        0.56549e03,
        0.68377e03,
        0.81007e03,
        0.94447e03,
        0.10872e04,
        0.12385e04,
        0.13988e04,
        0.15688e04,
        0.17490e04,
        0.19399e04,
        0.21424e04,
        0.23571e04,
        0.25848e04,
        0.28264e04,
        0.30828e04,
        0.33548e04,
        0.36434e04,
        0.39496e04,
        0.42745e04,
        0.46190e04,
        0.49845e04,
        0.53720e04,
        0.57828e04,
        0.62182e04,
        0.66796e04,
        0.71684e04,
        0.76862e04,
        0.82344e04,
        0.88149e04,
        0.94292e04,
        0.10079e05,
        0.10767e05,
        0.11494e05,
        0.12262e05,
        0.13074e05,
        0.13932e05,
        0.14839e05,
        0.15796e05,
        0.16806e05,
        0.17872e05,
        0.18997e05,
        0.20183e05,
        0.21434e05,
        0.22752e05,
        0.24141e05,
        0.25604e05,
        0.27145e05,
        0.28767e05,
        0.30475e05,
        0.32271e05,
        0.34160e05,
        0.36146e05,
        0.38234e05,
        0.40428e05,
        0.42733e05,
        0.45154e05,
        0.47696e05,
        0.50364e05,
        0.53163e05,
        0.56100e05,
        0.59180e05,
        0.62408e05,
        0.65792e05,
        0.69339e05,
        0.73053e05,
        0.76943e05,
        0.81016e05,
        0.85279e05,
        0.89740e05,
        0.94406e05,
        0.99287e05,
        0.10439e06,
        0.10972e06,
        0.11530e06,
        0.12112e06,
        0.12720e06,
        0.13355e06,
        0.14018e06,
        0.14711e06,
        0.15433e06,
        0.16186e06,
        0.16971e06,
        0.17791e06,
        0.18645e06,
        0.19534e06,
        0.20462e06,
        0.21428e06,
        0.22434e06,
        0.23481e06,
        0.24572e06,
        0.25706e06,
        0.26887e06,
        0.28116e06,
        0.29393e06,
        0.30722e06,
        0.32103e06,
        0.33539e06,
        0.35031e06,
        0.36581e06,
        0.38191e06,
        0.39864e06,
        0.41600e06,
        0.43403e06,
        0.45274e06,
        0.47215e06,
        0.49230e06,
        0.51319e06,
        0.53487e06,
        0.55734e06,
        0.58064e06,
        0.60478e06,
        0.62981e06,
        0.65574e06,
        0.68260e06,
    ]
)


#  --------------- HNO3 146: M = 12, I = 1 ---------------------
M = 12
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.15010e05,
        0.25316e05,
        0.37374e05,
        0.51216e05,
        0.67105e05,
        0.85473e05,
        0.10688e06,
        0.13201e06,
        0.16165e06,
        0.19671e06,
        0.23825e06,
        0.28749e06,
        0.34583e06,
        0.41490e06,
        0.49657e06,
        0.59302e06,
        0.70673e06,
        0.84054e06,
        0.99775e06,
        0.11821e07,
        0.13978e07,
        0.16498e07,
        0.19436e07,
        0.22855e07,
        0.26825e07,
        0.31428e07,
        0.36753e07,
        0.42903e07,
        0.49993e07,
        0.58151e07,
        0.67523e07,
        0.78269e07,
        0.90572e07,
        0.10463e08,
        0.12067e08,
        0.13895e08,
        0.15973e08,
        0.18333e08,
        0.21009e08,
        0.24039e08,
        0.27464e08,
        0.31331e08,
        0.35690e08,
        0.40597e08,
        0.46115e08,
        0.52310e08,
        0.59257e08,
        0.67037e08,
        0.75739e08,
        0.85461e08,
        0.96310e08,
        0.10840e09,
        0.12186e09,
        0.13683e09,
        0.15346e09,
        0.17191e09,
        0.19236e09,
        0.21501e09,
        0.24006e09,
        0.26774e09,
        0.29830e09,
        0.33200e09,
        0.36914e09,
        0.41002e09,
        0.45498e09,
        0.50438e09,
        0.55862e09,
        0.61812e09,
        0.68332e09,
        0.75473e09,
        0.83286e09,
        0.91828e09,
        0.10116e10,
        0.11134e10,
        0.12245e10,
        0.13456e10,
        0.14775e10,
        0.16210e10,
        0.17771e10,
        0.19467e10,
        0.21309e10,
        0.23309e10,
        0.25477e10,
        0.27827e10,
        0.30372e10,
        0.33127e10,
        0.36107e10,
        0.39329e10,
        0.42809e10,
        0.46567e10,
        0.50623e10,
        0.54997e10,
        0.59711e10,
        0.64789e10,
        0.70257e10,
        0.76140e10,
        0.82468e10,
        0.89269e10,
        0.96575e10,
        0.10442e11,
        0.11284e11,
        0.12187e11,
        0.13155e11,
        0.14193e11,
        0.15304e11,
        0.16494e11,
        0.17767e11,
        0.19129e11,
        0.20585e11,
        0.22140e11,
        0.23802e11,
        0.25576e11,
        0.27469e11,
        0.29489e11,
        0.31642e11,
        0.33937e11,
        0.36382e11,
        0.38985e11,
        0.41757e11,
    ]
)


#  --------------- HNO3 156: M = 12, I = 2 --------------------- NOT IN TIPS-2011
M = 12
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- OH 61: M = 13, I = 1 ---------------------
M = 13
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.20066e02,
        0.24774e02,
        0.30309e02,
        0.36357e02,
        0.42745e02,
        0.49371e02,
        0.56168e02,
        0.63093e02,
        0.70116e02,
        0.77217e02,
        0.84380e02,
        0.91594e02,
        0.98850e02,
        0.10614e03,
        0.11346e03,
        0.12081e03,
        0.12818e03,
        0.13557e03,
        0.14298e03,
        0.15041e03,
        0.15785e03,
        0.16531e03,
        0.17278e03,
        0.18027e03,
        0.18778e03,
        0.19530e03,
        0.20284e03,
        0.21040e03,
        0.21797e03,
        0.22556e03,
        0.23318e03,
        0.24082e03,
        0.24848e03,
        0.25617e03,
        0.26389e03,
        0.27163e03,
        0.27941e03,
        0.28721e03,
        0.29505e03,
        0.30292e03,
        0.31084e03,
        0.31878e03,
        0.32677e03,
        0.33480e03,
        0.34287e03,
        0.35099e03,
        0.35915e03,
        0.36736e03,
        0.37561e03,
        0.38391e03,
        0.39227e03,
        0.40067e03,
        0.40913e03,
        0.41764e03,
        0.42620e03,
        0.43482e03,
        0.44350e03,
        0.45223e03,
        0.46102e03,
        0.46987e03,
        0.47878e03,
        0.48775e03,
        0.49679e03,
        0.50588e03,
        0.51503e03,
        0.52425e03,
        0.53354e03,
        0.54288e03,
        0.55229e03,
        0.56177e03,
        0.57132e03,
        0.58092e03,
        0.59060e03,
        0.60035e03,
        0.61016e03,
        0.62004e03,
        0.62999e03,
        0.64001e03,
        0.65010e03,
        0.66025e03,
        0.67049e03,
        0.68078e03,
        0.69115e03,
        0.70160e03,
        0.71211e03,
        0.72269e03,
        0.73335e03,
        0.74408e03,
        0.75488e03,
        0.76576e03,
        0.77671e03,
        0.78773e03,
        0.79883e03,
        0.81000e03,
        0.82124e03,
        0.83256e03,
        0.84396e03,
        0.85542e03,
        0.86696e03,
        0.87858e03,
        0.89027e03,
        0.90204e03,
        0.91389e03,
        0.92580e03,
        0.93781e03,
        0.94988e03,
        0.96203e03,
        0.97425e03,
        0.98656e03,
        0.99893e03,
        0.10114e04,
        0.10239e04,
        0.10365e04,
        0.10492e04,
        0.10620e04,
        0.10748e04,
        0.10878e04,
        0.11007e04,
        0.11138e04,
    ]
)

#  --------------- OH 81: M = 13, I = 2 ---------------------
M = 13
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.20124e02,
        0.24876e02,
        0.30457e02,
        0.36553e02,
        0.42991e02,
        0.49666e02,
        0.56513e02,
        0.63489e02,
        0.70563e02,
        0.77715e02,
        0.84929e02,
        0.92195e02,
        0.99504e02,
        0.10685e03,
        0.11423e03,
        0.12164e03,
        0.12907e03,
        0.13654e03,
        0.14403e03,
        0.15154e03,
        0.15909e03,
        0.16666e03,
        0.17427e03,
        0.18191e03,
        0.18959e03,
        0.19731e03,
        0.20507e03,
        0.21287e03,
        0.22073e03,
        0.22863e03,
        0.23658e03,
        0.24459e03,
        0.25266e03,
        0.26078e03,
        0.26897e03,
        0.27722e03,
        0.28554e03,
        0.29393e03,
        0.30238e03,
        0.31091e03,
        0.31952e03,
        0.32820e03,
        0.33696e03,
        0.34579e03,
        0.35471e03,
        0.36371e03,
        0.37279e03,
        0.38196e03,
        0.39121e03,
        0.40055e03,
        0.40998e03,
        0.41949e03,
        0.42910e03,
        0.43879e03,
        0.44858e03,
        0.45845e03,
        0.46843e03,
        0.47849e03,
        0.48865e03,
        0.49890e03,
        0.50924e03,
        0.51969e03,
        0.53022e03,
        0.54086e03,
        0.55159e03,
        0.56242e03,
        0.57335e03,
        0.58437e03,
        0.59550e03,
        0.60673e03,
        0.61805e03,
        0.62947e03,
        0.64100e03,
        0.65263e03,
        0.66435e03,
        0.67618e03,
        0.68811e03,
        0.70014e03,
        0.71228e03,
        0.72451e03,
        0.73685e03,
        0.74929e03,
        0.76184e03,
        0.77449e03,
        0.78724e03,
        0.80009e03,
        0.81306e03,
        0.82612e03,
        0.83929e03,
        0.85256e03,
        0.86594e03,
        0.87942e03,
        0.89301e03,
        0.90670e03,
        0.92050e03,
        0.93440e03,
        0.94841e03,
        0.96253e03,
        0.97675e03,
        0.99108e03,
        0.10055e04,
        0.10201e04,
        0.10347e04,
        0.10495e04,
        0.10643e04,
        0.10793e04,
        0.10944e04,
        0.11096e04,
        0.11248e04,
        0.11402e04,
        0.11558e04,
        0.11714e04,
        0.11871e04,
        0.12029e04,
        0.12189e04,
        0.12349e04,
        0.12511e04,
        0.12673e04,
        0.12837e04,
    ]
)

#  --------------- OH 62: M = 13, I = 3 ---------------------
M = 13
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.41032e02,
        0.54704e02,
        0.70201e02,
        0.86985e02,
        0.10469e03,
        0.12306e03,
        0.14194e03,
        0.16119e03,
        0.18075e03,
        0.20054e03,
        0.22053e03,
        0.24068e03,
        0.26096e03,
        0.28135e03,
        0.30183e03,
        0.32241e03,
        0.34305e03,
        0.36376e03,
        0.38453e03,
        0.40535e03,
        0.42622e03,
        0.44714e03,
        0.46811e03,
        0.48913e03,
        0.51019e03,
        0.53131e03,
        0.55246e03,
        0.57368e03,
        0.59495e03,
        0.61627e03,
        0.63766e03,
        0.65912e03,
        0.68064e03,
        0.70223e03,
        0.72390e03,
        0.74565e03,
        0.76749e03,
        0.78941e03,
        0.81143e03,
        0.83355e03,
        0.85578e03,
        0.87810e03,
        0.90054e03,
        0.92310e03,
        0.94577e03,
        0.96857e03,
        0.99149e03,
        0.10145e04,
        0.10377e04,
        0.10611e04,
        0.10845e04,
        0.11081e04,
        0.11319e04,
        0.11558e04,
        0.11798e04,
        0.12040e04,
        0.12284e04,
        0.12529e04,
        0.12776e04,
        0.13025e04,
        0.13275e04,
        0.13527e04,
        0.13781e04,
        0.14036e04,
        0.14293e04,
        0.14552e04,
        0.14813e04,
        0.15076e04,
        0.15340e04,
        0.15606e04,
        0.15874e04,
        0.16144e04,
        0.16416e04,
        0.16690e04,
        0.16965e04,
        0.17243e04,
        0.17522e04,
        0.17804e04,
        0.18087e04,
        0.18373e04,
        0.18660e04,
        0.18949e04,
        0.19241e04,
        0.19534e04,
        0.19829e04,
        0.20127e04,
        0.20426e04,
        0.20727e04,
        0.21031e04,
        0.21336e04,
        0.21644e04,
        0.21954e04,
        0.22266e04,
        0.22579e04,
        0.22895e04,
        0.23213e04,
        0.23534e04,
        0.23856e04,
        0.24180e04,
        0.24506e04,
        0.24835e04,
        0.25166e04,
        0.25499e04,
        0.25834e04,
        0.26171e04,
        0.26510e04,
        0.26852e04,
        0.27195e04,
        0.27541e04,
        0.27889e04,
        0.28239e04,
        0.28592e04,
        0.28946e04,
        0.29303e04,
        0.29661e04,
        0.30023e04,
        0.30386e04,
        0.30751e04,
        0.31119e04,
    ]
)


#  --------------- HF 19: M = 14, I = 1 ---------------------
M = 14
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.95958e01,
        0.12933e02,
        0.16295e02,
        0.19666e02,
        0.23043e02,
        0.26425e02,
        0.29809e02,
        0.33195e02,
        0.36584e02,
        0.39974e02,
        0.43366e02,
        0.46759e02,
        0.50154e02,
        0.53550e02,
        0.56947e02,
        0.60346e02,
        0.63746e02,
        0.67148e02,
        0.70550e02,
        0.73955e02,
        0.77361e02,
        0.80769e02,
        0.84179e02,
        0.87591e02,
        0.91006e02,
        0.94424e02,
        0.97846e02,
        0.10127e03,
        0.10470e03,
        0.10813e03,
        0.11157e03,
        0.11502e03,
        0.11847e03,
        0.12193e03,
        0.12540e03,
        0.12888e03,
        0.13236e03,
        0.13586e03,
        0.13936e03,
        0.14288e03,
        0.14641e03,
        0.14995e03,
        0.15351e03,
        0.15708e03,
        0.16066e03,
        0.16426e03,
        0.16788e03,
        0.17151e03,
        0.17516e03,
        0.17882e03,
        0.18251e03,
        0.18621e03,
        0.18994e03,
        0.19368e03,
        0.19745e03,
        0.20123e03,
        0.20504e03,
        0.20887e03,
        0.21272e03,
        0.21659e03,
        0.22049e03,
        0.22441e03,
        0.22836e03,
        0.23233e03,
        0.23632e03,
        0.24034e03,
        0.24439e03,
        0.24846e03,
        0.25255e03,
        0.25668e03,
        0.26083e03,
        0.26501e03,
        0.26921e03,
        0.27344e03,
        0.27770e03,
        0.28199e03,
        0.28631e03,
        0.29066e03,
        0.29503e03,
        0.29944e03,
        0.30387e03,
        0.30833e03,
        0.31282e03,
        0.31735e03,
        0.32190e03,
        0.32648e03,
        0.33110e03,
        0.33574e03,
        0.34042e03,
        0.34512e03,
        0.34986e03,
        0.35463e03,
        0.35943e03,
        0.36426e03,
        0.36913e03,
        0.37402e03,
        0.37895e03,
        0.38391e03,
        0.38891e03,
        0.39393e03,
        0.39899e03,
        0.40408e03,
        0.40921e03,
        0.41436e03,
        0.41955e03,
        0.42478e03,
        0.43004e03,
        0.43533e03,
        0.44065e03,
        0.44601e03,
        0.45140e03,
        0.45683e03,
        0.46229e03,
        0.46779e03,
        0.47332e03,
        0.47888e03,
        0.48448e03,
        0.49011e03,
        0.49578e03,
    ]
)

#  --------------- HF 29: M = 14, I = 2 --------------------- not in TIPS-2011
M = 14
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- HСl 15: M = 15, I = 1 --------------------
M = 15
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.34775e02,
        0.48060e02,
        0.61370e02,
        0.74692e02,
        0.88024e02,
        0.10136e03,
        0.11471e03,
        0.12806e03,
        0.14141e03,
        0.15478e03,
        0.16814e03,
        0.18151e03,
        0.19489e03,
        0.20827e03,
        0.22166e03,
        0.23506e03,
        0.24847e03,
        0.26189e03,
        0.27533e03,
        0.28878e03,
        0.30225e03,
        0.31575e03,
        0.32928e03,
        0.34284e03,
        0.35645e03,
        0.37009e03,
        0.38378e03,
        0.39753e03,
        0.41134e03,
        0.42521e03,
        0.43914e03,
        0.45316e03,
        0.46725e03,
        0.48142e03,
        0.49568e03,
        0.51003e03,
        0.52448e03,
        0.53902e03,
        0.55368e03,
        0.56843e03,
        0.58330e03,
        0.59829e03,
        0.61339e03,
        0.62862e03,
        0.64396e03,
        0.65944e03,
        0.67504e03,
        0.69078e03,
        0.70665e03,
        0.72265e03,
        0.73880e03,
        0.75508e03,
        0.77151e03,
        0.78809e03,
        0.80481e03,
        0.82168e03,
        0.83870e03,
        0.85587e03,
        0.87320e03,
        0.89068e03,
        0.90832e03,
        0.92611e03,
        0.94407e03,
        0.96218e03,
        0.98046e03,
        0.99889e03,
        0.10175e04,
        0.10363e04,
        0.10552e04,
        0.10743e04,
        0.10936e04,
        0.11130e04,
        0.11326e04,
        0.11524e04,
        0.11723e04,
        0.11924e04,
        0.12127e04,
        0.12332e04,
        0.12538e04,
        0.12746e04,
        0.12956e04,
        0.13168e04,
        0.13381e04,
        0.13597e04,
        0.13814e04,
        0.14032e04,
        0.14253e04,
        0.14475e04,
        0.14700e04,
        0.14926e04,
        0.15153e04,
        0.15383e04,
        0.15615e04,
        0.15848e04,
        0.16083e04,
        0.16320e04,
        0.16559e04,
        0.16800e04,
        0.17043e04,
        0.17287e04,
        0.17533e04,
        0.17782e04,
        0.18032e04,
        0.18284e04,
        0.18538e04,
        0.18794e04,
        0.19051e04,
        0.19311e04,
        0.19573e04,
        0.19836e04,
        0.20102e04,
        0.20369e04,
        0.20638e04,
        0.20910e04,
        0.21183e04,
        0.21458e04,
        0.21735e04,
        0.22014e04,
        0.22295e04,
    ]
)


#  --------------- HСl 17: M = 15, I = 2 ---------------------
M = 15
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.34823e02,
        0.48128e02,
        0.61458e02,
        0.74801e02,
        0.88152e02,
        0.10151e03,
        0.11488e03,
        0.12825e03,
        0.14162e03,
        0.15500e03,
        0.16839e03,
        0.18178e03,
        0.19518e03,
        0.20858e03,
        0.22199e03,
        0.23541e03,
        0.24884e03,
        0.26228e03,
        0.27574e03,
        0.28921e03,
        0.30270e03,
        0.31622e03,
        0.32977e03,
        0.34336e03,
        0.35698e03,
        0.37065e03,
        0.38436e03,
        0.39813e03,
        0.41196e03,
        0.42585e03,
        0.43981e03,
        0.45384e03,
        0.46796e03,
        0.48215e03,
        0.49644e03,
        0.51081e03,
        0.52528e03,
        0.53986e03,
        0.55453e03,
        0.56932e03,
        0.58421e03,
        0.59922e03,
        0.61435e03,
        0.62960e03,
        0.64498e03,
        0.66048e03,
        0.67611e03,
        0.69187e03,
        0.70777e03,
        0.72381e03,
        0.73998e03,
        0.75630e03,
        0.77276e03,
        0.78936e03,
        0.80612e03,
        0.82302e03,
        0.84007e03,
        0.85727e03,
        0.87463e03,
        0.89215e03,
        0.90982e03,
        0.92765e03,
        0.94563e03,
        0.96378e03,
        0.98209e03,
        0.10006e04,
        0.10192e04,
        0.10380e04,
        0.10570e04,
        0.10761e04,
        0.10954e04,
        0.11149e04,
        0.11345e04,
        0.11543e04,
        0.11743e04,
        0.11945e04,
        0.12148e04,
        0.12353e04,
        0.12560e04,
        0.12768e04,
        0.12979e04,
        0.13191e04,
        0.13405e04,
        0.13620e04,
        0.13838e04,
        0.14057e04,
        0.14278e04,
        0.14501e04,
        0.14726e04,
        0.14952e04,
        0.15180e04,
        0.15410e04,
        0.15642e04,
        0.15876e04,
        0.16112e04,
        0.16349e04,
        0.16589e04,
        0.16830e04,
        0.17073e04,
        0.17318e04,
        0.17565e04,
        0.17814e04,
        0.18064e04,
        0.18317e04,
        0.18572e04,
        0.18828e04,
        0.19086e04,
        0.19346e04,
        0.19609e04,
        0.19873e04,
        0.20139e04,
        0.20406e04,
        0.20676e04,
        0.20948e04,
        0.21222e04,
        0.21498e04,
        0.21775e04,
        0.22055e04,
        0.22337e04,
    ]
)


#  --------------- HСl 25: M = 15, I = 3 --------------------- not in TIPS-2011
M = 15
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- HСl 27: M = 15, I = 4 --------------------- not in TIPS-2011
M = 15
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- HBr 19: M = 16, I = 1 ---------------------
M = 16
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.42744e02,
        0.59373e02,
        0.76023e02,
        0.92685e02,
        0.10936e03,
        0.12604e03,
        0.14272e03,
        0.15942e03,
        0.17612e03,
        0.19282e03,
        0.20954e03,
        0.22626e03,
        0.24299e03,
        0.25973e03,
        0.27648e03,
        0.29325e03,
        0.31004e03,
        0.32686e03,
        0.34371e03,
        0.36060e03,
        0.37753e03,
        0.39451e03,
        0.41156e03,
        0.42868e03,
        0.44587e03,
        0.46314e03,
        0.48051e03,
        0.49798e03,
        0.51556e03,
        0.53325e03,
        0.55106e03,
        0.56900e03,
        0.58708e03,
        0.60530e03,
        0.62367e03,
        0.64219e03,
        0.66088e03,
        0.67972e03,
        0.69874e03,
        0.71793e03,
        0.73730e03,
        0.75685e03,
        0.77659e03,
        0.79652e03,
        0.81664e03,
        0.83696e03,
        0.85748e03,
        0.87820e03,
        0.89914e03,
        0.92028e03,
        0.94163e03,
        0.96319e03,
        0.98498e03,
        0.10070e04,
        0.10292e04,
        0.10516e04,
        0.10743e04,
        0.10972e04,
        0.11203e04,
        0.11437e04,
        0.11673e04,
        0.11911e04,
        0.12151e04,
        0.12394e04,
        0.12640e04,
        0.12887e04,
        0.13137e04,
        0.13390e04,
        0.13645e04,
        0.13902e04,
        0.14162e04,
        0.14424e04,
        0.14689e04,
        0.14956e04,
        0.15226e04,
        0.15498e04,
        0.15773e04,
        0.16050e04,
        0.16330e04,
        0.16612e04,
        0.16897e04,
        0.17185e04,
        0.17475e04,
        0.17767e04,
        0.18062e04,
        0.18360e04,
        0.18660e04,
        0.18963e04,
        0.19269e04,
        0.19577e04,
        0.19888e04,
        0.20202e04,
        0.20518e04,
        0.20837e04,
        0.21158e04,
        0.21482e04,
        0.21809e04,
        0.22139e04,
        0.22471e04,
        0.22806e04,
        0.23143e04,
        0.23484e04,
        0.23827e04,
        0.24173e04,
        0.24521e04,
        0.24873e04,
        0.25227e04,
        0.25584e04,
        0.25943e04,
        0.26306e04,
        0.26671e04,
        0.27039e04,
        0.27409e04,
        0.27783e04,
        0.28159e04,
        0.28538e04,
        0.28920e04,
        0.29305e04,
        0.29693e04,
    ]
)


#  --------------- HBr 11: M = 16, I = 2 ---------------------
M = 16
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.42756e02,
        0.59390e02,
        0.76045e02,
        0.92713e02,
        0.10939e03,
        0.12607e03,
        0.14277e03,
        0.15947e03,
        0.17617e03,
        0.19288e03,
        0.20960e03,
        0.22633e03,
        0.24306e03,
        0.25981e03,
        0.27656e03,
        0.29334e03,
        0.31014e03,
        0.32696e03,
        0.34381e03,
        0.36071e03,
        0.37764e03,
        0.39464e03,
        0.41169e03,
        0.42881e03,
        0.44601e03,
        0.46329e03,
        0.48066e03,
        0.49813e03,
        0.51572e03,
        0.53341e03,
        0.55123e03,
        0.56918e03,
        0.58727e03,
        0.60549e03,
        0.62387e03,
        0.64240e03,
        0.66109e03,
        0.67994e03,
        0.69896e03,
        0.71816e03,
        0.73754e03,
        0.75710e03,
        0.77684e03,
        0.79678e03,
        0.81691e03,
        0.83724e03,
        0.85776e03,
        0.87850e03,
        0.89943e03,
        0.92058e03,
        0.94194e03,
        0.96352e03,
        0.98531e03,
        0.10073e04,
        0.10295e04,
        0.10520e04,
        0.10747e04,
        0.10976e04,
        0.11207e04,
        0.11441e04,
        0.11677e04,
        0.11915e04,
        0.12156e04,
        0.12399e04,
        0.12644e04,
        0.12892e04,
        0.13142e04,
        0.13395e04,
        0.13650e04,
        0.13907e04,
        0.14167e04,
        0.14429e04,
        0.14694e04,
        0.14961e04,
        0.15231e04,
        0.15504e04,
        0.15778e04,
        0.16056e04,
        0.16336e04,
        0.16618e04,
        0.16903e04,
        0.17191e04,
        0.17481e04,
        0.17773e04,
        0.18069e04,
        0.18367e04,
        0.18667e04,
        0.18970e04,
        0.19276e04,
        0.19584e04,
        0.19895e04,
        0.20209e04,
        0.20525e04,
        0.20844e04,
        0.21166e04,
        0.21490e04,
        0.21817e04,
        0.22147e04,
        0.22479e04,
        0.22814e04,
        0.23152e04,
        0.23492e04,
        0.23835e04,
        0.24181e04,
        0.24530e04,
        0.24882e04,
        0.25236e04,
        0.25593e04,
        0.25952e04,
        0.26315e04,
        0.26680e04,
        0.27048e04,
        0.27419e04,
        0.27793e04,
        0.28169e04,
        0.28549e04,
        0.28931e04,
        0.29316e04,
        0.29703e04,
    ]
)


#  --------------- HBr 29: M = 16, I = 3 --------------------- not in TIPS-2011
M = 16
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- HBr 21: M = 16, I = 4 --------------------- not in TIPS-2011
M = 16
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- HI 17: M = 17, I = 1 ---------------------
M = 17
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.82031e02,
        0.11447e03,
        0.14694e03,
        0.17943e03,
        0.21194e03,
        0.24445e03,
        0.27699e03,
        0.30953e03,
        0.34209e03,
        0.37466e03,
        0.40725e03,
        0.43986e03,
        0.47249e03,
        0.50517e03,
        0.53789e03,
        0.57068e03,
        0.60354e03,
        0.63650e03,
        0.66957e03,
        0.70278e03,
        0.73614e03,
        0.76967e03,
        0.80340e03,
        0.83735e03,
        0.87153e03,
        0.90596e03,
        0.94067e03,
        0.97566e03,
        0.10110e04,
        0.10466e04,
        0.10826e04,
        0.11189e04,
        0.11555e04,
        0.11926e04,
        0.12300e04,
        0.12679e04,
        0.13061e04,
        0.13448e04,
        0.13839e04,
        0.14235e04,
        0.14635e04,
        0.15039e04,
        0.15448e04,
        0.15862e04,
        0.16280e04,
        0.16704e04,
        0.17132e04,
        0.17565e04,
        0.18003e04,
        0.18446e04,
        0.18894e04,
        0.19347e04,
        0.19806e04,
        0.20269e04,
        0.20738e04,
        0.21212e04,
        0.21691e04,
        0.22176e04,
        0.22666e04,
        0.23162e04,
        0.23662e04,
        0.24169e04,
        0.24680e04,
        0.25198e04,
        0.25720e04,
        0.26249e04,
        0.26783e04,
        0.27322e04,
        0.27867e04,
        0.28418e04,
        0.28975e04,
        0.29537e04,
        0.30105e04,
        0.30678e04,
        0.31258e04,
        0.31843e04,
        0.32434e04,
        0.33031e04,
        0.33633e04,
        0.34242e04,
        0.34856e04,
        0.35477e04,
        0.36103e04,
        0.36735e04,
        0.37373e04,
        0.38018e04,
        0.38668e04,
        0.39324e04,
        0.39986e04,
        0.40654e04,
        0.41329e04,
        0.42009e04,
        0.42696e04,
        0.43388e04,
        0.44087e04,
        0.44792e04,
        0.45503e04,
        0.46221e04,
        0.46944e04,
        0.47674e04,
        0.48410e04,
        0.49152e04,
        0.49901e04,
        0.50656e04,
        0.51417e04,
        0.52185e04,
        0.52959e04,
        0.53739e04,
        0.54526e04,
        0.55319e04,
        0.56118e04,
        0.56924e04,
        0.57736e04,
        0.58555e04,
        0.59380e04,
        0.60212e04,
        0.61050e04,
        0.61895e04,
        0.62746e04,
    ]
)


#  --------------- HI 27: M = 17, I = 2 --------------------- not in TIPS-2011
M = 17
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- ClO 56: M = 18, I = 1 ---------------------
M = 18
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.53847e03,
        0.76580e03,
        0.10017e04,
        0.12511e04,
        0.15168e04,
        0.18001e04,
        0.21014e04,
        0.24206e04,
        0.27577e04,
        0.31127e04,
        0.34857e04,
        0.38765e04,
        0.42854e04,
        0.47124e04,
        0.51575e04,
        0.56208e04,
        0.61025e04,
        0.66026e04,
        0.71211e04,
        0.76582e04,
        0.82138e04,
        0.87882e04,
        0.93813e04,
        0.99932e04,
        0.10624e05,
        0.11273e05,
        0.11942e05,
        0.12629e05,
        0.13336e05,
        0.14061e05,
        0.14806e05,
        0.15570e05,
        0.16353e05,
        0.17155e05,
        0.17976e05,
        0.18816e05,
        0.19676e05,
        0.20555e05,
        0.21453e05,
        0.22371e05,
        0.23308e05,
        0.24264e05,
        0.25240e05,
        0.26236e05,
        0.27250e05,
        0.28284e05,
        0.29338e05,
        0.30412e05,
        0.31505e05,
        0.32617e05,
        0.33749e05,
        0.34901e05,
        0.36072e05,
        0.37263e05,
        0.38474e05,
        0.39705e05,
        0.40955e05,
        0.42225e05,
        0.43515e05,
        0.44825e05,
        0.46154e05,
        0.47504e05,
        0.48873e05,
        0.50262e05,
        0.51672e05,
        0.53101e05,
        0.54549e05,
        0.56019e05,
        0.57508e05,
        0.59017e05,
        0.60546e05,
        0.62095e05,
        0.63665e05,
        0.65254e05,
        0.66864e05,
        0.68494e05,
        0.70144e05,
        0.71814e05,
        0.73504e05,
        0.75215e05,
        0.76946e05,
        0.78698e05,
        0.80470e05,
        0.82261e05,
        0.84074e05,
        0.85907e05,
        0.87760e05,
        0.89633e05,
        0.91527e05,
        0.93442e05,
        0.95377e05,
        0.97333e05,
        0.99309e05,
        0.10131e06,
        0.10332e06,
        0.10536e06,
        0.10742e06,
        0.10950e06,
        0.11160e06,
        0.11372e06,
        0.11586e06,
        0.11802e06,
        0.12020e06,
        0.12241e06,
        0.12463e06,
        0.12688e06,
        0.12914e06,
        0.13143e06,
        0.13374e06,
        0.13607e06,
        0.13842e06,
        0.14079e06,
        0.14318e06,
        0.14559e06,
        0.14802e06,
        0.15048e06,
        0.15295e06,
        0.15545e06,
        0.15797e06,
    ]
)


#  --------------- ClO 76: M = 18, I = 2 ---------------------
M = 18
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.54775e03,
        0.77899e03,
        0.10189e04,
        0.12726e04,
        0.15430e04,
        0.18313e04,
        0.21378e04,
        0.24627e04,
        0.28059e04,
        0.31674e04,
        0.35472e04,
        0.39454e04,
        0.43621e04,
        0.47972e04,
        0.52508e04,
        0.57232e04,
        0.62143e04,
        0.67242e04,
        0.72531e04,
        0.78010e04,
        0.83678e04,
        0.89537e04,
        0.95589e04,
        0.10183e05,
        0.10827e05,
        0.11490e05,
        0.12172e05,
        0.12874e05,
        0.13595e05,
        0.14335e05,
        0.15095e05,
        0.15875e05,
        0.16674e05,
        0.17493e05,
        0.18332e05,
        0.19190e05,
        0.20068e05,
        0.20965e05,
        0.21882e05,
        0.22820e05,
        0.23776e05,
        0.24753e05,
        0.25750e05,
        0.26766e05,
        0.27803e05,
        0.28859e05,
        0.29935e05,
        0.31032e05,
        0.32148e05,
        0.33284e05,
        0.34441e05,
        0.35617e05,
        0.36814e05,
        0.38031e05,
        0.39267e05,
        0.40524e05,
        0.41802e05,
        0.43099e05,
        0.44417e05,
        0.45755e05,
        0.47113e05,
        0.48492e05,
        0.49891e05,
        0.51310e05,
        0.52750e05,
        0.54210e05,
        0.55690e05,
        0.57191e05,
        0.58713e05,
        0.60255e05,
        0.61817e05,
        0.63400e05,
        0.65004e05,
        0.66628e05,
        0.68272e05,
        0.69938e05,
        0.71624e05,
        0.73331e05,
        0.75058e05,
        0.76806e05,
        0.78575e05,
        0.80364e05,
        0.82175e05,
        0.84006e05,
        0.85858e05,
        0.87731e05,
        0.89625e05,
        0.91539e05,
        0.93475e05,
        0.95431e05,
        0.97409e05,
        0.99407e05,
        0.10143e06,
        0.10347e06,
        0.10553e06,
        0.10761e06,
        0.10972e06,
        0.11184e06,
        0.11399e06,
        0.11615e06,
        0.11834e06,
        0.12055e06,
        0.12278e06,
        0.12503e06,
        0.12731e06,
        0.12960e06,
        0.13192e06,
        0.13425e06,
        0.13661e06,
        0.13899e06,
        0.14139e06,
        0.14382e06,
        0.14626e06,
        0.14873e06,
        0.15121e06,
        0.15372e06,
        0.15625e06,
        0.15880e06,
        0.16138e06,
    ]
)


#  --------------- OCS 622: M = 19, I = 1 ---------------------
M = 19
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.20609e03,
        0.29199e03,
        0.37861e03,
        0.46737e03,
        0.56024e03,
        0.65929e03,
        0.76649e03,
        0.88361e03,
        0.10123e04,
        0.11541e04,
        0.13105e04,
        0.14829e04,
        0.16728e04,
        0.18818e04,
        0.21113e04,
        0.23629e04,
        0.26383e04,
        0.29391e04,
        0.32672e04,
        0.36245e04,
        0.40128e04,
        0.44343e04,
        0.48911e04,
        0.53853e04,
        0.59193e04,
        0.64956e04,
        0.71166e04,
        0.77849e04,
        0.85033e04,
        0.92746e04,
        0.10102e05,
        0.10988e05,
        0.11936e05,
        0.12949e05,
        0.14032e05,
        0.15186e05,
        0.16416e05,
        0.17726e05,
        0.19120e05,
        0.20601e05,
        0.22173e05,
        0.23842e05,
        0.25611e05,
        0.27484e05,
        0.29468e05,
        0.31566e05,
        0.33783e05,
        0.36124e05,
        0.38595e05,
        0.41202e05,
        0.43949e05,
        0.46842e05,
        0.49888e05,
        0.53092e05,
        0.56460e05,
        0.59999e05,
        0.63716e05,
        0.67616e05,
        0.71708e05,
        0.75997e05,
        0.80491e05,
        0.85197e05,
        0.90124e05,
        0.95278e05,
        0.10067e06,
        0.10630e06,
        0.11219e06,
        0.11833e06,
        0.12475e06,
        0.13144e06,
        0.13842e06,
        0.14570e06,
        0.15328e06,
        0.16117e06,
        0.16940e06,
        0.17795e06,
        0.18686e06,
        0.19611e06,
        0.20574e06,
        0.21574e06,
        0.22613e06,
        0.23692e06,
        0.24813e06,
        0.25975e06,
        0.27182e06,
        0.28433e06,
        0.29730e06,
        0.31074e06,
        0.32467e06,
        0.33909e06,
        0.35403e06,
        0.36950e06,
        0.38551e06,
        0.40207e06,
        0.41920e06,
        0.43691e06,
        0.45522e06,
        0.47415e06,
        0.49370e06,
        0.51390e06,
        0.53476e06,
        0.55629e06,
        0.57852e06,
        0.60146e06,
        0.62513e06,
        0.64954e06,
        0.67471e06,
        0.70067e06,
        0.72742e06,
        0.75499e06,
        0.78339e06,
        0.81265e06,
        0.84279e06,
        0.87381e06,
        0.90576e06,
        0.93863e06,
        0.97246e06,
        0.10073e07,
        0.10431e07,
    ]
)


#  --------------- OCS 624: M = 19, I = 2 ---------------------
M = 19
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.21125e03,
        0.29930e03,
        0.38809e03,
        0.47911e03,
        0.57437e03,
        0.67603e03,
        0.78610e03,
        0.90643e03,
        0.10387e04,
        0.11846e04,
        0.13456e04,
        0.15231e04,
        0.17188e04,
        0.19342e04,
        0.21709e04,
        0.24304e04,
        0.27145e04,
        0.30250e04,
        0.33638e04,
        0.37328e04,
        0.41339e04,
        0.45694e04,
        0.50415e04,
        0.55524e04,
        0.61045e04,
        0.67004e04,
        0.73427e04,
        0.80340e04,
        0.87773e04,
        0.95755e04,
        0.10432e05,
        0.11349e05,
        0.12330e05,
        0.13380e05,
        0.14500e05,
        0.15696e05,
        0.16970e05,
        0.18327e05,
        0.19770e05,
        0.21305e05,
        0.22934e05,
        0.24663e05,
        0.26497e05,
        0.28439e05,
        0.30495e05,
        0.32669e05,
        0.34968e05,
        0.37396e05,
        0.39958e05,
        0.42661e05,
        0.45510e05,
        0.48511e05,
        0.51669e05,
        0.54993e05,
        0.58487e05,
        0.62159e05,
        0.66014e05,
        0.70061e05,
        0.74306e05,
        0.78757e05,
        0.83421e05,
        0.88305e05,
        0.93418e05,
        0.98767e05,
        0.10436e06,
        0.11021e06,
        0.11632e06,
        0.12270e06,
        0.12936e06,
        0.13631e06,
        0.14355e06,
        0.15111e06,
        0.15898e06,
        0.16718e06,
        0.17572e06,
        0.18460e06,
        0.19385e06,
        0.20346e06,
        0.21346e06,
        0.22385e06,
        0.23464e06,
        0.24585e06,
        0.25748e06,
        0.26956e06,
        0.28209e06,
        0.29509e06,
        0.30856e06,
        0.32252e06,
        0.33699e06,
        0.35198e06,
        0.36750e06,
        0.38357e06,
        0.40020e06,
        0.41741e06,
        0.43521e06,
        0.45362e06,
        0.47264e06,
        0.49231e06,
        0.51263e06,
        0.53362e06,
        0.55529e06,
        0.57768e06,
        0.60078e06,
        0.62462e06,
        0.64922e06,
        0.67459e06,
        0.70075e06,
        0.72773e06,
        0.75554e06,
        0.78419e06,
        0.81372e06,
        0.84413e06,
        0.87546e06,
        0.90771e06,
        0.94092e06,
        0.97509e06,
        0.10103e07,
        0.10464e07,
        0.10837e07,
    ]
)


#  --------------- OCS 632: M = 19, I = 3 ---------------------
M = 19
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.41351e03,
        0.58591e03,
        0.76004e03,
        0.93907e03,
        0.11273e04,
        0.13289e04,
        0.15481e04,
        0.17884e04,
        0.20533e04,
        0.23459e04,
        0.26692e04,
        0.30264e04,
        0.34205e04,
        0.38547e04,
        0.43323e04,
        0.48565e04,
        0.54309e04,
        0.60592e04,
        0.67451e04,
        0.74928e04,
        0.83064e04,
        0.91903e04,
        0.10149e05,
        0.11187e05,
        0.12310e05,
        0.13523e05,
        0.14831e05,
        0.16240e05,
        0.17756e05,
        0.19384e05,
        0.21132e05,
        0.23005e05,
        0.25011e05,
        0.27157e05,
        0.29449e05,
        0.31896e05,
        0.34506e05,
        0.37286e05,
        0.40245e05,
        0.43392e05,
        0.46735e05,
        0.50284e05,
        0.54048e05,
        0.58038e05,
        0.62263e05,
        0.66733e05,
        0.71460e05,
        0.76455e05,
        0.81728e05,
        0.87292e05,
        0.93159e05,
        0.99341e05,
        0.10585e06,
        0.11270e06,
        0.11991e06,
        0.12748e06,
        0.13543e06,
        0.14378e06,
        0.15255e06,
        0.16174e06,
        0.17137e06,
        0.18146e06,
        0.19202e06,
        0.20308e06,
        0.21465e06,
        0.22674e06,
        0.23937e06,
        0.25257e06,
        0.26635e06,
        0.28073e06,
        0.29573e06,
        0.31137e06,
        0.32767e06,
        0.34466e06,
        0.36235e06,
        0.38076e06,
        0.39992e06,
        0.41985e06,
        0.44057e06,
        0.46211e06,
        0.48450e06,
        0.50775e06,
        0.53189e06,
        0.55695e06,
        0.58295e06,
        0.60992e06,
        0.63789e06,
        0.66688e06,
        0.69693e06,
        0.72806e06,
        0.76030e06,
        0.79368e06,
        0.82823e06,
        0.86399e06,
        0.90097e06,
        0.93923e06,
        0.97878e06,
        0.10197e07,
        0.10619e07,
        0.11056e07,
        0.11506e07,
        0.11972e07,
        0.12453e07,
        0.12949e07,
        0.13460e07,
        0.13988e07,
        0.14533e07,
        0.15094e07,
        0.15673e07,
        0.16270e07,
        0.16884e07,
        0.17518e07,
        0.18170e07,
        0.18842e07,
        0.19533e07,
        0.20245e07,
        0.20978e07,
        0.21732e07,
        0.22507e07,
    ]
)


#  --------------- OCS 623: M = 19, I = 4 ---------------------
M = 19
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.83485e03,
        0.11828e04,
        0.15337e04,
        0.18934e04,
        0.22697e04,
        0.26712e04,
        0.31059e04,
        0.35809e04,
        0.41030e04,
        0.46785e04,
        0.53133e04,
        0.60135e04,
        0.67850e04,
        0.76338e04,
        0.85663e04,
        0.95888e04,
        0.10708e05,
        0.11931e05,
        0.13265e05,
        0.14718e05,
        0.16298e05,
        0.18012e05,
        0.19870e05,
        0.21881e05,
        0.24054e05,
        0.26399e05,
        0.28926e05,
        0.31646e05,
        0.34570e05,
        0.37710e05,
        0.41077e05,
        0.44685e05,
        0.48545e05,
        0.52672e05,
        0.57078e05,
        0.61780e05,
        0.66790e05,
        0.72125e05,
        0.77801e05,
        0.83833e05,
        0.90239e05,
        0.97036e05,
        0.10424e06,
        0.11188e06,
        0.11996e06,
        0.12850e06,
        0.13754e06,
        0.14708e06,
        0.15715e06,
        0.16777e06,
        0.17896e06,
        0.19076e06,
        0.20317e06,
        0.21623e06,
        0.22996e06,
        0.24438e06,
        0.25953e06,
        0.27543e06,
        0.29211e06,
        0.30959e06,
        0.32791e06,
        0.34710e06,
        0.36718e06,
        0.38820e06,
        0.41017e06,
        0.43314e06,
        0.45713e06,
        0.48219e06,
        0.50835e06,
        0.53564e06,
        0.56409e06,
        0.59376e06,
        0.62468e06,
        0.65688e06,
        0.69041e06,
        0.72530e06,
        0.76161e06,
        0.79937e06,
        0.83862e06,
        0.87941e06,
        0.92179e06,
        0.96581e06,
        0.10115e07,
        0.10589e07,
        0.11081e07,
        0.11591e07,
        0.12120e07,
        0.12669e07,
        0.13237e07,
        0.13825e07,
        0.14435e07,
        0.15066e07,
        0.15718e07,
        0.16394e07,
        0.17093e07,
        0.17815e07,
        0.18562e07,
        0.19334e07,
        0.20132e07,
        0.20956e07,
        0.21807e07,
        0.22685e07,
        0.23592e07,
        0.24528e07,
        0.25494e07,
        0.26490e07,
        0.27517e07,
        0.28576e07,
        0.29667e07,
        0.30792e07,
        0.31951e07,
        0.33145e07,
        0.34374e07,
        0.35640e07,
        0.36943e07,
        0.38285e07,
        0.39665e07,
        0.41085e07,
        0.42546e07,
    ]
)


#  --------------- OCS 822: M = 19, I = 5 ---------------------
M = 19
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.21967e03,
        0.31126e03,
        0.40370e03,
        0.49862e03,
        0.59823e03,
        0.70481e03,
        0.82050e03,
        0.94724e03,
        0.10868e04,
        0.12409e04,
        0.14112e04,
        0.15993e04,
        0.18067e04,
        0.20353e04,
        0.22866e04,
        0.25624e04,
        0.28645e04,
        0.31950e04,
        0.35558e04,
        0.39490e04,
        0.43767e04,
        0.48413e04,
        0.53452e04,
        0.58909e04,
        0.64810e04,
        0.71182e04,
        0.78053e04,
        0.85454e04,
        0.93413e04,
        0.10196e05,
        0.11114e05,
        0.12098e05,
        0.13151e05,
        0.14277e05,
        0.15480e05,
        0.16764e05,
        0.18133e05,
        0.19592e05,
        0.21144e05,
        0.22794e05,
        0.24548e05,
        0.26409e05,
        0.28383e05,
        0.30475e05,
        0.32689e05,
        0.35033e05,
        0.37511e05,
        0.40128e05,
        0.42892e05,
        0.45808e05,
        0.48882e05,
        0.52121e05,
        0.55532e05,
        0.59121e05,
        0.62895e05,
        0.66861e05,
        0.71028e05,
        0.75402e05,
        0.79991e05,
        0.84803e05,
        0.89847e05,
        0.95130e05,
        0.10066e06,
        0.10645e06,
        0.11251e06,
        0.11883e06,
        0.12545e06,
        0.13236e06,
        0.13957e06,
        0.14710e06,
        0.15495e06,
        0.16313e06,
        0.17166e06,
        0.18055e06,
        0.18980e06,
        0.19944e06,
        0.20946e06,
        0.21989e06,
        0.23073e06,
        0.24200e06,
        0.25371e06,
        0.26587e06,
        0.27850e06,
        0.29161e06,
        0.30521e06,
        0.31931e06,
        0.33394e06,
        0.34910e06,
        0.36482e06,
        0.38109e06,
        0.39795e06,
        0.41541e06,
        0.43348e06,
        0.45217e06,
        0.47151e06,
        0.49151e06,
        0.51219e06,
        0.53356e06,
        0.55565e06,
        0.57847e06,
        0.60204e06,
        0.62637e06,
        0.65149e06,
        0.67742e06,
        0.70417e06,
        0.73176e06,
        0.76023e06,
        0.78957e06,
        0.81982e06,
        0.85100e06,
        0.88313e06,
        0.91622e06,
        0.95031e06,
        0.98541e06,
        0.10216e07,
        0.10587e07,
        0.10970e07,
        0.11364e07,
        0.11769e07,
    ]
)


#  --------------- H2CO 126: M = 20, I = 2 ---------------------
M = 20
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.25934e03,
        0.43623e03,
        0.64143e03,
        0.87152e03,
        0.11241e04,
        0.13975e04,
        0.16906e04,
        0.20029e04,
        0.23344e04,
        0.26857e04,
        0.30577e04,
        0.34518e04,
        0.38698e04,
        0.43138e04,
        0.47860e04,
        0.52890e04,
        0.58256e04,
        0.63985e04,
        0.70109e04,
        0.76660e04,
        0.83673e04,
        0.91184e04,
        0.99230e04,
        0.10785e05,
        0.11710e05,
        0.12700e05,
        0.13762e05,
        0.14900e05,
        0.16119e05,
        0.17425e05,
        0.18823e05,
        0.20320e05,
        0.21923e05,
        0.23637e05,
        0.25471e05,
        0.27432e05,
        0.29527e05,
        0.31765e05,
        0.34155e05,
        0.36706e05,
        0.39428e05,
        0.42330e05,
        0.45424e05,
        0.48720e05,
        0.52231e05,
        0.55968e05,
        0.59945e05,
        0.64175e05,
        0.68672e05,
        0.73450e05,
        0.78526e05,
        0.83915e05,
        0.89634e05,
        0.95701e05,
        0.10213e06,
        0.10895e06,
        0.11618e06,
        0.12383e06,
        0.13193e06,
        0.14049e06,
        0.14956e06,
        0.15914e06,
        0.16927e06,
        0.17997e06,
        0.19127e06,
        0.20320e06,
        0.21578e06,
        0.22906e06,
        0.24306e06,
        0.25782e06,
        0.27336e06,
        0.28974e06,
        0.30698e06,
        0.32513e06,
        0.34422e06,
        0.36430e06,
        0.38542e06,
        0.40761e06,
        0.43093e06,
        0.45542e06,
        0.48114e06,
        0.50813e06,
        0.53646e06,
        0.56617e06,
        0.59733e06,
        0.63000e06,
        0.66423e06,
        0.70010e06,
        0.73767e06,
        0.77701e06,
        0.81818e06,
        0.86127e06,
        0.90635e06,
        0.95349e06,
        0.10028e07,
        0.10543e07,
        0.11082e07,
        0.11644e07,
        0.12232e07,
        0.12845e07,
        0.13485e07,
        0.14154e07,
        0.14851e07,
        0.15578e07,
        0.16337e07,
        0.17127e07,
        0.17952e07,
        0.18810e07,
        0.19705e07,
        0.20637e07,
        0.21607e07,
        0.22617e07,
        0.23669e07,
        0.24763e07,
        0.25901e07,
        0.27085e07,
        0.28316e07,
        0.29596e07,
        0.30926e07,
    ]
)


#  --------------- H2CO 136: M = 20, I = 2 ---------------------
M = 20
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.53173e03,
        0.89447e03,
        0.13153e04,
        0.17871e04,
        0.23051e04,
        0.28658e04,
        0.34669e04,
        0.41073e04,
        0.47872e04,
        0.55074e04,
        0.62702e04,
        0.70785e04,
        0.79357e04,
        0.88462e04,
        0.98147e04,
        0.10846e05,
        0.11946e05,
        0.13121e05,
        0.14377e05,
        0.15721e05,
        0.17159e05,
        0.18699e05,
        0.20349e05,
        0.22118e05,
        0.24013e05,
        0.26045e05,
        0.28222e05,
        0.30555e05,
        0.33055e05,
        0.35733e05,
        0.38601e05,
        0.41671e05,
        0.44958e05,
        0.48474e05,
        0.52235e05,
        0.56255e05,
        0.60552e05,
        0.65142e05,
        0.70043e05,
        0.75275e05,
        0.80856e05,
        0.86808e05,
        0.93152e05,
        0.99913e05,
        0.10711e06,
        0.11478e06,
        0.12293e06,
        0.13161e06,
        0.14083e06,
        0.15063e06,
        0.16104e06,
        0.17209e06,
        0.18382e06,
        0.19626e06,
        0.20945e06,
        0.22343e06,
        0.23825e06,
        0.25394e06,
        0.27054e06,
        0.28812e06,
        0.30671e06,
        0.32636e06,
        0.34713e06,
        0.36907e06,
        0.39224e06,
        0.41671e06,
        0.44252e06,
        0.46975e06,
        0.49845e06,
        0.52872e06,
        0.56060e06,
        0.59418e06,
        0.62954e06,
        0.66676e06,
        0.70591e06,
        0.74710e06,
        0.79040e06,
        0.83591e06,
        0.88373e06,
        0.93395e06,
        0.98669e06,
        0.10421e07,
        0.11001e07,
        0.11611e07,
        0.12250e07,
        0.12920e07,
        0.13622e07,
        0.14357e07,
        0.15128e07,
        0.15934e07,
        0.16779e07,
        0.17662e07,
        0.18587e07,
        0.19554e07,
        0.20565e07,
        0.21621e07,
        0.22725e07,
        0.23879e07,
        0.25084e07,
        0.26342e07,
        0.27655e07,
        0.29026e07,
        0.30456e07,
        0.31947e07,
        0.33502e07,
        0.35124e07,
        0.36814e07,
        0.38575e07,
        0.40410e07,
        0.42321e07,
        0.44311e07,
        0.46382e07,
        0.48538e07,
        0.50782e07,
        0.53116e07,
        0.55544e07,
        0.58068e07,
        0.60693e07,
        0.63421e07,
    ]
)


#  --------------- H2CO 128: M = 20, I = 3 ---------------------
M = 20
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.27198e03,
        0.45755e03,
        0.67282e03,
        0.91421e03,
        0.11792e04,
        0.14660e04,
        0.17735e04,
        0.21012e04,
        0.24490e04,
        0.28175e04,
        0.32077e04,
        0.36212e04,
        0.40598e04,
        0.45256e04,
        0.50211e04,
        0.55488e04,
        0.61116e04,
        0.67127e04,
        0.73552e04,
        0.80426e04,
        0.87783e04,
        0.95663e04,
        0.10410e05,
        0.11315e05,
        0.12285e05,
        0.13324e05,
        0.14438e05,
        0.15632e05,
        0.16911e05,
        0.18281e05,
        0.19748e05,
        0.21319e05,
        0.23000e05,
        0.24799e05,
        0.26723e05,
        0.28780e05,
        0.30978e05,
        0.33326e05,
        0.35834e05,
        0.38510e05,
        0.41365e05,
        0.44410e05,
        0.47656e05,
        0.51115e05,
        0.54798e05,
        0.58719e05,
        0.62891e05,
        0.67329e05,
        0.72047e05,
        0.77060e05,
        0.82385e05,
        0.88039e05,
        0.94039e05,
        0.10040e06,
        0.10715e06,
        0.11431e06,
        0.12189e06,
        0.12991e06,
        0.13841e06,
        0.14740e06,
        0.15691e06,
        0.16696e06,
        0.17759e06,
        0.18882e06,
        0.20067e06,
        0.21318e06,
        0.22639e06,
        0.24032e06,
        0.25501e06,
        0.27049e06,
        0.28680e06,
        0.30398e06,
        0.32207e06,
        0.34111e06,
        0.36114e06,
        0.38221e06,
        0.40436e06,
        0.42765e06,
        0.45211e06,
        0.47781e06,
        0.50479e06,
        0.53311e06,
        0.56283e06,
        0.59400e06,
        0.62669e06,
        0.66097e06,
        0.69688e06,
        0.73451e06,
        0.77393e06,
        0.81520e06,
        0.85840e06,
        0.90360e06,
        0.95090e06,
        0.10004e07,
        0.10521e07,
        0.11061e07,
        0.11626e07,
        0.12216e07,
        0.12833e07,
        0.13476e07,
        0.14148e07,
        0.14849e07,
        0.15581e07,
        0.16344e07,
        0.17140e07,
        0.17969e07,
        0.18834e07,
        0.19735e07,
        0.20674e07,
        0.21651e07,
        0.22669e07,
        0.23729e07,
        0.24832e07,
        0.25980e07,
        0.27174e07,
        0.28416e07,
        0.29708e07,
        0.31050e07,
        0.32446e07,
    ]
)


#  --------------- HOCl 165: M = 21, I = 1 ---------------------
M = 21
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.17041e04,
        0.28708e04,
        0.42250e04,
        0.57456e04,
        0.74211e04,
        0.92470e04,
        0.11225e05,
        0.13359e05,
        0.15657e05,
        0.18129e05,
        0.20785e05,
        0.23637e05,
        0.26696e05,
        0.29974e05,
        0.33484e05,
        0.37239e05,
        0.41252e05,
        0.45536e05,
        0.50105e05,
        0.54973e05,
        0.60152e05,
        0.65659e05,
        0.71507e05,
        0.77711e05,
        0.84286e05,
        0.91249e05,
        0.98614e05,
        0.10640e06,
        0.11462e06,
        0.12330e06,
        0.13244e06,
        0.14208e06,
        0.15222e06,
        0.16289e06,
        0.17411e06,
        0.18589e06,
        0.19825e06,
        0.21123e06,
        0.22483e06,
        0.23908e06,
        0.25400e06,
        0.26962e06,
        0.28596e06,
        0.30303e06,
        0.32087e06,
        0.33950e06,
        0.35895e06,
        0.37923e06,
        0.40038e06,
        0.42243e06,
        0.44539e06,
        0.46930e06,
        0.49419e06,
        0.52008e06,
        0.54700e06,
        0.57498e06,
        0.60406e06,
        0.63426e06,
        0.66562e06,
        0.69816e06,
        0.73192e06,
        0.76692e06,
        0.80322e06,
        0.84083e06,
        0.87979e06,
        0.92014e06,
        0.96192e06,
        0.10052e07,
        0.10499e07,
        0.10961e07,
        0.11440e07,
        0.11934e07,
        0.12445e07,
        0.12973e07,
        0.13518e07,
        0.14081e07,
        0.14661e07,
        0.15261e07,
        0.15879e07,
        0.16516e07,
        0.17174e07,
        0.17851e07,
        0.18550e07,
        0.19269e07,
        0.20010e07,
        0.20773e07,
        0.21559e07,
        0.22367e07,
        0.23200e07,
        0.24056e07,
        0.24936e07,
        0.25842e07,
        0.26773e07,
        0.27730e07,
        0.28714e07,
        0.29724e07,
        0.30763e07,
        0.31829e07,
        0.32924e07,
        0.34049e07,
        0.35203e07,
        0.36387e07,
        0.37603e07,
        0.38850e07,
        0.40129e07,
        0.41441e07,
        0.42786e07,
        0.44165e07,
        0.45579e07,
        0.47028e07,
        0.48512e07,
        0.50033e07,
        0.51592e07,
        0.53187e07,
        0.54822e07,
        0.56495e07,
        0.58208e07,
        0.59961e07,
        0.61755e07,
    ]
)


#  --------------- HOCl 167: M = 21, I = 2 ---------------------
M = 21
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.17342e04,
        0.29215e04,
        0.42998e04,
        0.58473e04,
        0.75524e04,
        0.94107e04,
        0.11423e05,
        0.13595e05,
        0.15935e05,
        0.18450e05,
        0.21154e05,
        0.24056e05,
        0.27168e05,
        0.30505e05,
        0.34077e05,
        0.37899e05,
        0.41983e05,
        0.46343e05,
        0.50993e05,
        0.55947e05,
        0.61218e05,
        0.66822e05,
        0.72774e05,
        0.79088e05,
        0.85780e05,
        0.92866e05,
        0.10036e06,
        0.10829e06,
        0.11665e06,
        0.12548e06,
        0.13479e06,
        0.14460e06,
        0.15492e06,
        0.16578e06,
        0.17719e06,
        0.18918e06,
        0.20177e06,
        0.21497e06,
        0.22881e06,
        0.24332e06,
        0.25851e06,
        0.27440e06,
        0.29102e06,
        0.30840e06,
        0.32656e06,
        0.34552e06,
        0.36531e06,
        0.38595e06,
        0.40748e06,
        0.42991e06,
        0.45328e06,
        0.47762e06,
        0.50295e06,
        0.52929e06,
        0.55669e06,
        0.58517e06,
        0.61477e06,
        0.64550e06,
        0.67741e06,
        0.71053e06,
        0.74489e06,
        0.78052e06,
        0.81745e06,
        0.85573e06,
        0.89539e06,
        0.93645e06,
        0.97897e06,
        0.10230e07,
        0.10685e07,
        0.11156e07,
        0.11643e07,
        0.12146e07,
        0.12666e07,
        0.13203e07,
        0.13757e07,
        0.14330e07,
        0.14921e07,
        0.15531e07,
        0.16160e07,
        0.16809e07,
        0.17478e07,
        0.18168e07,
        0.18878e07,
        0.19611e07,
        0.20365e07,
        0.21141e07,
        0.21941e07,
        0.22764e07,
        0.23611e07,
        0.24482e07,
        0.25378e07,
        0.26300e07,
        0.27248e07,
        0.28222e07,
        0.29223e07,
        0.30251e07,
        0.31308e07,
        0.32393e07,
        0.33508e07,
        0.34652e07,
        0.35827e07,
        0.37032e07,
        0.38269e07,
        0.39539e07,
        0.40840e07,
        0.42176e07,
        0.43545e07,
        0.44948e07,
        0.46387e07,
        0.47861e07,
        0.49372e07,
        0.50920e07,
        0.52506e07,
        0.54130e07,
        0.55793e07,
        0.57496e07,
        0.59239e07,
        0.61024e07,
        0.62850e07,
    ]
)


#  --------------- N2 44: M = 22, I = 1 ---------------------
M = 22
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.95487e02,
        0.13466e03,
        0.17386e03,
        0.21307e03,
        0.25230e03,
        0.29154e03,
        0.33080e03,
        0.37008e03,
        0.40937e03,
        0.44868e03,
        0.48800e03,
        0.52736e03,
        0.56674e03,
        0.60616e03,
        0.64562e03,
        0.68515e03,
        0.72475e03,
        0.76445e03,
        0.80426e03,
        0.84420e03,
        0.88430e03,
        0.92457e03,
        0.96505e03,
        0.10057e04,
        0.10467e04,
        0.10879e04,
        0.11293e04,
        0.11711e04,
        0.12132e04,
        0.12556e04,
        0.12984e04,
        0.13416e04,
        0.13851e04,
        0.14291e04,
        0.14734e04,
        0.15182e04,
        0.15635e04,
        0.16091e04,
        0.16553e04,
        0.17019e04,
        0.17490e04,
        0.17965e04,
        0.18446e04,
        0.18932e04,
        0.19422e04,
        0.19918e04,
        0.20419e04,
        0.20926e04,
        0.21437e04,
        0.21954e04,
        0.22477e04,
        0.23004e04,
        0.23538e04,
        0.24077e04,
        0.24621e04,
        0.25171e04,
        0.25727e04,
        0.26288e04,
        0.26856e04,
        0.27428e04,
        0.28007e04,
        0.28591e04,
        0.29181e04,
        0.29777e04,
        0.30379e04,
        0.30986e04,
        0.31600e04,
        0.32219e04,
        0.32844e04,
        0.33475e04,
        0.34112e04,
        0.34755e04,
        0.35404e04,
        0.36059e04,
        0.36720e04,
        0.37387e04,
        0.38060e04,
        0.38739e04,
        0.39424e04,
        0.40115e04,
        0.40812e04,
        0.41515e04,
        0.42224e04,
        0.42939e04,
        0.43661e04,
        0.44388e04,
        0.45122e04,
        0.45861e04,
        0.46607e04,
        0.47359e04,
        0.48117e04,
        0.48882e04,
        0.49652e04,
        0.50428e04,
        0.51211e04,
        0.52000e04,
        0.52795e04,
        0.53596e04,
        0.54404e04,
        0.55217e04,
        0.56037e04,
        0.56863e04,
        0.57695e04,
        0.58533e04,
        0.59378e04,
        0.60229e04,
        0.61086e04,
        0.61950e04,
        0.62819e04,
        0.63695e04,
        0.64577e04,
        0.65465e04,
        0.66360e04,
        0.67261e04,
        0.68168e04,
        0.69081e04,
        0.70001e04,
        0.70927e04,
        0.71859e04,
    ]
)


#  --------------- N2 45: M = 22, I = 2 --------------------- not in TIPS-2011
M = 22
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- HCN 124: M = 23, I = 1 ---------------------
M = 23
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.17143e03,
        0.24209e03,
        0.31285e03,
        0.38392e03,
        0.45582e03,
        0.52929e03,
        0.60515e03,
        0.68424e03,
        0.76731e03,
        0.85505e03,
        0.94805e03,
        0.10468e04,
        0.11519e04,
        0.12637e04,
        0.13826e04,
        0.15090e04,
        0.16435e04,
        0.17863e04,
        0.19378e04,
        0.20985e04,
        0.22689e04,
        0.24492e04,
        0.26401e04,
        0.28418e04,
        0.30550e04,
        0.32801e04,
        0.35176e04,
        0.37680e04,
        0.40318e04,
        0.43097e04,
        0.46021e04,
        0.49097e04,
        0.52330e04,
        0.55727e04,
        0.59294e04,
        0.63038e04,
        0.66964e04,
        0.71081e04,
        0.75396e04,
        0.79915e04,
        0.84646e04,
        0.89596e04,
        0.94774e04,
        0.10019e05,
        0.10585e05,
        0.11176e05,
        0.11793e05,
        0.12437e05,
        0.13108e05,
        0.13809e05,
        0.14540e05,
        0.15301e05,
        0.16094e05,
        0.16919e05,
        0.17779e05,
        0.18673e05,
        0.19603e05,
        0.20570e05,
        0.21575e05,
        0.22619e05,
        0.23704e05,
        0.24831e05,
        0.26000e05,
        0.27213e05,
        0.28472e05,
        0.29778e05,
        0.31131e05,
        0.32534e05,
        0.33987e05,
        0.35493e05,
        0.37052e05,
        0.38666e05,
        0.40336e05,
        0.42064e05,
        0.43852e05,
        0.45701e05,
        0.47612e05,
        0.49587e05,
        0.51629e05,
        0.53738e05,
        0.55916e05,
        0.58165e05,
        0.60486e05,
        0.62883e05,
        0.65355e05,
        0.67905e05,
        0.70536e05,
        0.73249e05,
        0.76045e05,
        0.78927e05,
        0.81897e05,
        0.84957e05,
        0.88108e05,
        0.91354e05,
        0.94696e05,
        0.98136e05,
        0.10168e06,
        0.10532e06,
        0.10907e06,
        0.11292e06,
        0.11689e06,
        0.12096e06,
        0.12516e06,
        0.12946e06,
        0.13389e06,
        0.13844e06,
        0.14311e06,
        0.14791e06,
        0.15284e06,
        0.15790e06,
        0.16310e06,
        0.16843e06,
        0.17391e06,
        0.17953e06,
        0.18529e06,
        0.19120e06,
        0.19726e06,
        0.20348e06,
        0.20986e06,
    ]
)


#  --------------- HCN 134: M = 23, I = 2 ---------------------
M = 23
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.35186e03,
        0.49693e03,
        0.64221e03,
        0.78815e03,
        0.93585e03,
        0.10868e04,
        0.12428e04,
        0.14056e04,
        0.15766e04,
        0.17574e04,
        0.19491e04,
        0.21528e04,
        0.23695e04,
        0.26002e04,
        0.28457e04,
        0.31068e04,
        0.33845e04,
        0.36795e04,
        0.39926e04,
        0.43249e04,
        0.46770e04,
        0.50500e04,
        0.54447e04,
        0.58621e04,
        0.63032e04,
        0.67690e04,
        0.72606e04,
        0.77789e04,
        0.83252e04,
        0.89005e04,
        0.95062e04,
        0.10143e05,
        0.10813e05,
        0.11517e05,
        0.12256e05,
        0.13032e05,
        0.13846e05,
        0.14699e05,
        0.15593e05,
        0.16530e05,
        0.17511e05,
        0.18538e05,
        0.19612e05,
        0.20734e05,
        0.21908e05,
        0.23134e05,
        0.24414e05,
        0.25750e05,
        0.27145e05,
        0.28599e05,
        0.30115e05,
        0.31694e05,
        0.33340e05,
        0.35054e05,
        0.36838e05,
        0.38694e05,
        0.40625e05,
        0.42633e05,
        0.44720e05,
        0.46889e05,
        0.49142e05,
        0.51481e05,
        0.53910e05,
        0.56430e05,
        0.59045e05,
        0.61757e05,
        0.64568e05,
        0.67482e05,
        0.70502e05,
        0.73630e05,
        0.76869e05,
        0.80223e05,
        0.83694e05,
        0.87285e05,
        0.91000e05,
        0.94843e05,
        0.98815e05,
        0.10292e06,
        0.10716e06,
        0.11155e06,
        0.11608e06,
        0.12075e06,
        0.12558e06,
        0.13056e06,
        0.13570e06,
        0.14100e06,
        0.14647e06,
        0.15211e06,
        0.15793e06,
        0.16392e06,
        0.17009e06,
        0.17646e06,
        0.18301e06,
        0.18976e06,
        0.19671e06,
        0.20387e06,
        0.21123e06,
        0.21881e06,
        0.22660e06,
        0.23462e06,
        0.24287e06,
        0.25135e06,
        0.26007e06,
        0.26903e06,
        0.27824e06,
        0.28771e06,
        0.29743e06,
        0.30742e06,
        0.31767e06,
        0.32820e06,
        0.33901e06,
        0.35011e06,
        0.36150e06,
        0.37319e06,
        0.38518e06,
        0.39749e06,
        0.41010e06,
        0.42304e06,
        0.43631e06,
    ]
)


#  --------------- HCN 135: M = 23, I = 3 ---------------------
M = 23
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.11863e03,
        0.16755e03,
        0.21653e03,
        0.26576e03,
        0.31559e03,
        0.36656e03,
        0.41926e03,
        0.47428e03,
        0.53214e03,
        0.59333e03,
        0.65824e03,
        0.72727e03,
        0.80074e03,
        0.87898e03,
        0.96227e03,
        0.10509e04,
        0.11452e04,
        0.12454e04,
        0.13518e04,
        0.14647e04,
        0.15844e04,
        0.17112e04,
        0.18455e04,
        0.19875e04,
        0.21377e04,
        0.22962e04,
        0.24636e04,
        0.26402e04,
        0.28263e04,
        0.30224e04,
        0.32289e04,
        0.34461e04,
        0.36745e04,
        0.39145e04,
        0.41667e04,
        0.44314e04,
        0.47092e04,
        0.50005e04,
        0.53059e04,
        0.56259e04,
        0.59609e04,
        0.63116e04,
        0.66785e04,
        0.70622e04,
        0.74633e04,
        0.78823e04,
        0.83200e04,
        0.87769e04,
        0.92536e04,
        0.97509e04,
        0.10269e05,
        0.10810e05,
        0.11373e05,
        0.11959e05,
        0.12570e05,
        0.13205e05,
        0.13866e05,
        0.14554e05,
        0.15268e05,
        0.16011e05,
        0.16782e05,
        0.17583e05,
        0.18415e05,
        0.19279e05,
        0.20174e05,
        0.21103e05,
        0.22067e05,
        0.23065e05,
        0.24100e05,
        0.25172e05,
        0.26282e05,
        0.27432e05,
        0.28622e05,
        0.29853e05,
        0.31127e05,
        0.32445e05,
        0.33807e05,
        0.35215e05,
        0.36670e05,
        0.38174e05,
        0.39727e05,
        0.41330e05,
        0.42986e05,
        0.44695e05,
        0.46459e05,
        0.48278e05,
        0.50155e05,
        0.52091e05,
        0.54086e05,
        0.56143e05,
        0.58263e05,
        0.60447e05,
        0.62696e05,
        0.65013e05,
        0.67399e05,
        0.69856e05,
        0.72384e05,
        0.74986e05,
        0.77663e05,
        0.80416e05,
        0.83249e05,
        0.86161e05,
        0.89156e05,
        0.92233e05,
        0.95397e05,
        0.98648e05,
        0.10199e06,
        0.10542e06,
        0.10894e06,
        0.11256e06,
        0.11627e06,
        0.12009e06,
        0.12400e06,
        0.12802e06,
        0.13214e06,
        0.13636e06,
        0.14070e06,
        0.14515e06,
        0.14971e06,
    ]
)


#  --------------- CH3Cl 215: M = 24, I = 1 ---------------------
M = 24
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.50529e04,
        0.85123e04,
        0.12528e05,
        0.17036e05,
        0.22005e05,
        0.27429e05,
        0.33325e05,
        0.39734e05,
        0.46713e05,
        0.54336e05,
        0.62690e05,
        0.71876e05,
        0.82006e05,
        0.93204e05,
        0.10560e06,
        0.11936e06,
        0.13463e06,
        0.15158e06,
        0.17043e06,
        0.19137e06,
        0.21464e06,
        0.24049e06,
        0.26920e06,
        0.30107e06,
        0.33642e06,
        0.37563e06,
        0.41907e06,
        0.46719e06,
        0.52045e06,
        0.57936e06,
        0.64448e06,
        0.71641e06,
        0.79582e06,
        0.88341e06,
        0.97997e06,
        0.10863e07,
        0.12034e07,
        0.13323e07,
        0.14739e07,
        0.16295e07,
        0.18003e07,
        0.19877e07,
        0.21932e07,
        0.24183e07,
        0.26649e07,
        0.29346e07,
        0.32296e07,
        0.35519e07,
        0.39039e07,
        0.42881e07,
        0.47072e07,
        0.51639e07,
        0.56615e07,
        0.62032e07,
        0.67926e07,
        0.74335e07,
        0.81299e07,
        0.88862e07,
        0.97071e07,
        0.10598e08,
        0.11563e08,
        0.12609e08,
        0.13742e08,
        0.14968e08,
        0.16294e08,
        0.17728e08,
        0.19277e08,
        0.20950e08,
        0.22756e08,
        0.24704e08,
        0.26805e08,
        0.29069e08,
        0.31507e08,
        0.34132e08,
        0.36957e08,
        0.39995e08,
        0.43260e08,
        0.46769e08,
        0.50538e08,
        0.54583e08,
        0.58923e08,
        0.63578e08,
        0.68568e08,
        0.73914e08,
        0.79640e08,
        0.85770e08,
        0.92329e08,
        0.99345e08,
        0.10685e09,
        0.11486e09,
        0.12342e09,
        0.13257e09,
        0.14233e09,
        0.15274e09,
        0.16384e09,
        0.17568e09,
        0.18829e09,
        0.20173e09,
        0.21604e09,
        0.23127e09,
        0.24748e09,
        0.26471e09,
        0.28304e09,
        0.30252e09,
        0.32322e09,
        0.34520e09,
        0.36853e09,
        0.39330e09,
        0.41958e09,
        0.44745e09,
        0.47701e09,
        0.50833e09,
        0.54151e09,
        0.57667e09,
        0.61389e09,
        0.65329e09,
        0.69498e09,
        0.73909e09,
        0.78573e09,
    ]
)


#  --------------- CH3Cl 217: M = 24, I = 2 ---------------------
M = 24
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.51327e04,
        0.86469e04,
        0.12726e05,
        0.17306e05,
        0.22354e05,
        0.27863e05,
        0.33853e05,
        0.40364e05,
        0.47453e05,
        0.55197e05,
        0.63684e05,
        0.73016e05,
        0.83306e05,
        0.94681e05,
        0.10728e06,
        0.12125e06,
        0.13676e06,
        0.15399e06,
        0.17313e06,
        0.19441e06,
        0.21804e06,
        0.24430e06,
        0.27347e06,
        0.30584e06,
        0.34176e06,
        0.38158e06,
        0.42572e06,
        0.47460e06,
        0.52871e06,
        0.58855e06,
        0.65471e06,
        0.72778e06,
        0.80844e06,
        0.89743e06,
        0.99552e06,
        0.11036e07,
        0.12225e07,
        0.13534e07,
        0.14973e07,
        0.16553e07,
        0.18289e07,
        0.20193e07,
        0.22280e07,
        0.24567e07,
        0.27072e07,
        0.29812e07,
        0.32808e07,
        0.36083e07,
        0.39659e07,
        0.43562e07,
        0.47819e07,
        0.52459e07,
        0.57514e07,
        0.63017e07,
        0.69005e07,
        0.75515e07,
        0.82590e07,
        0.90273e07,
        0.98613e07,
        0.10766e08,
        0.11747e08,
        0.12809e08,
        0.13960e08,
        0.15206e08,
        0.16553e08,
        0.18010e08,
        0.19584e08,
        0.21283e08,
        0.23118e08,
        0.25097e08,
        0.27231e08,
        0.29531e08,
        0.32008e08,
        0.34674e08,
        0.37544e08,
        0.40630e08,
        0.43948e08,
        0.47513e08,
        0.51341e08,
        0.55451e08,
        0.59860e08,
        0.64589e08,
        0.69658e08,
        0.75089e08,
        0.80906e08,
        0.87134e08,
        0.93797e08,
        0.10092e09,
        0.10854e09,
        0.11669e09,
        0.12539e09,
        0.13467e09,
        0.14459e09,
        0.15517e09,
        0.16645e09,
        0.17847e09,
        0.19129e09,
        0.20494e09,
        0.21948e09,
        0.23495e09,
        0.25141e09,
        0.26893e09,
        0.28754e09,
        0.30733e09,
        0.32836e09,
        0.35069e09,
        0.37440e09,
        0.39956e09,
        0.42626e09,
        0.45457e09,
        0.48460e09,
        0.51642e09,
        0.55013e09,
        0.58585e09,
        0.62366e09,
        0.66369e09,
        0.70605e09,
        0.75085e09,
        0.79824e09,
    ]
)


#  --------------- H2O2 1661: M = 25, I = 1 ---------------------
M = 25
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.62392e03,
        0.10958e04,
        0.16692e04,
        0.23492e04,
        0.31427e04,
        0.40574e04,
        0.51014e04,
        0.62840e04,
        0.76157e04,
        0.91085e04,
        0.10776e05,
        0.12633e05,
        0.14696e05,
        0.16983e05,
        0.19515e05,
        0.22312e05,
        0.25396e05,
        0.28792e05,
        0.32526e05,
        0.36625e05,
        0.41118e05,
        0.46036e05,
        0.51410e05,
        0.57275e05,
        0.63667e05,
        0.70623e05,
        0.78185e05,
        0.86394e05,
        0.95295e05,
        0.10493e06,
        0.11536e06,
        0.12662e06,
        0.13878e06,
        0.15188e06,
        0.16600e06,
        0.18118e06,
        0.19750e06,
        0.21503e06,
        0.23383e06,
        0.25398e06,
        0.27556e06,
        0.29864e06,
        0.32333e06,
        0.34970e06,
        0.37784e06,
        0.40786e06,
        0.43985e06,
        0.47392e06,
        0.51018e06,
        0.54874e06,
        0.58972e06,
        0.63324e06,
        0.67943e06,
        0.72843e06,
        0.78037e06,
        0.83540e06,
        0.89366e06,
        0.95530e06,
        0.10205e07,
        0.10894e07,
        0.11622e07,
        0.12391e07,
        0.13202e07,
        0.14057e07,
        0.14959e07,
        0.15909e07,
        0.16910e07,
        0.17963e07,
        0.19072e07,
        0.20237e07,
        0.21463e07,
        0.22750e07,
        0.24102e07,
        0.25522e07,
        0.27012e07,
        0.28575e07,
        0.30213e07,
        0.31931e07,
        0.33730e07,
        0.35615e07,
        0.37588e07,
        0.39653e07,
        0.41813e07,
        0.44072e07,
        0.46433e07,
        0.48901e07,
        0.51479e07,
        0.54171e07,
        0.56982e07,
        0.59915e07,
        0.62976e07,
        0.66167e07,
        0.69495e07,
        0.72963e07,
        0.76577e07,
        0.80342e07,
        0.84262e07,
        0.88343e07,
        0.92591e07,
        0.97011e07,
        0.10161e08,
        0.10639e08,
        0.11136e08,
        0.11652e08,
        0.12189e08,
        0.12746e08,
        0.13325e08,
        0.13926e08,
        0.14550e08,
        0.15198e08,
        0.15870e08,
        0.16566e08,
        0.17289e08,
        0.18038e08,
        0.18814e08,
        0.19619e08,
        0.20452e08,
        0.21315e08,
        0.22209e08,
    ]
)


#  --------------- C2H2 1221: M = 26, I = 1 ---------------------
M = 26
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.71617e02,
        0.10121e03,
        0.13092e03,
        0.16104e03,
        0.19218e03,
        0.22509e03,
        0.26062e03,
        0.29959e03,
        0.34281e03,
        0.39103e03,
        0.44503e03,
        0.50558e03,
        0.57346e03,
        0.64950e03,
        0.73457e03,
        0.82960e03,
        0.93557e03,
        0.10535e04,
        0.11846e04,
        0.13301e04,
        0.14911e04,
        0.16692e04,
        0.18658e04,
        0.20825e04,
        0.23211e04,
        0.25833e04,
        0.28711e04,
        0.31867e04,
        0.35323e04,
        0.39102e04,
        0.43230e04,
        0.47735e04,
        0.52645e04,
        0.57991e04,
        0.63807e04,
        0.70127e04,
        0.76988e04,
        0.84430e04,
        0.92495e04,
        0.10123e05,
        0.11067e05,
        0.12088e05,
        0.13191e05,
        0.14381e05,
        0.15664e05,
        0.17047e05,
        0.18536e05,
        0.20137e05,
        0.21859e05,
        0.23710e05,
        0.25696e05,
        0.27827e05,
        0.30112e05,
        0.32561e05,
        0.35183e05,
        0.37990e05,
        0.40991e05,
        0.44199e05,
        0.47626e05,
        0.51285e05,
        0.55189e05,
        0.59353e05,
        0.63791e05,
        0.68518e05,
        0.73551e05,
        0.78908e05,
        0.84604e05,
        0.90661e05,
        0.97095e05,
        0.10393e06,
        0.11118e06,
        0.11888e06,
        0.12704e06,
        0.13569e06,
        0.14486e06,
        0.15457e06,
        0.16485e06,
        0.17572e06,
        0.18722e06,
        0.19938e06,
        0.21223e06,
        0.22581e06,
        0.24014e06,
        0.25527e06,
        0.27123e06,
        0.28807e06,
        0.30582e06,
        0.32452e06,
        0.34423e06,
        0.36498e06,
        0.38683e06,
        0.40982e06,
        0.43401e06,
        0.45944e06,
        0.48618e06,
        0.51428e06,
        0.54380e06,
        0.57480e06,
        0.60735e06,
        0.64151e06,
        0.67735e06,
        0.71494e06,
        0.75436e06,
        0.79568e06,
        0.83898e06,
        0.88434e06,
        0.93184e06,
        0.98158e06,
        0.10336e07,
        0.10881e07,
        0.11451e07,
        0.12047e07,
        0.12670e07,
        0.13321e07,
        0.14002e07,
        0.14713e07,
        0.15455e07,
        0.16231e07,
        0.17040e07,
    ]
)


#  --------------- C2H2 1231: M = 26, I = 2 ---------------------
M = 26
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.28647e03,
        0.40486e03,
        0.52369e03,
        0.64419e03,
        0.76874e03,
        0.90040e03,
        0.10425e04,
        0.11984e04,
        0.13713e04,
        0.15642e04,
        0.17802e04,
        0.20223e04,
        0.22939e04,
        0.25981e04,
        0.29384e04,
        0.33185e04,
        0.37424e04,
        0.42142e04,
        0.47386e04,
        0.53203e04,
        0.59646e04,
        0.66769e04,
        0.74633e04,
        0.83302e04,
        0.92845e04,
        0.10333e05,
        0.11485e05,
        0.12747e05,
        0.14129e05,
        0.15641e05,
        0.17292e05,
        0.19094e05,
        0.21058e05,
        0.23197e05,
        0.25523e05,
        0.28051e05,
        0.30796e05,
        0.33773e05,
        0.36999e05,
        0.40492e05,
        0.44270e05,
        0.48354e05,
        0.52765e05,
        0.57525e05,
        0.62658e05,
        0.68189e05,
        0.74144e05,
        0.80551e05,
        0.87439e05,
        0.94840e05,
        0.10279e06,
        0.11131e06,
        0.12045e06,
        0.13025e06,
        0.14074e06,
        0.15196e06,
        0.16397e06,
        0.17680e06,
        0.19051e06,
        0.20514e06,
        0.22076e06,
        0.23742e06,
        0.25517e06,
        0.27408e06,
        0.29421e06,
        0.31564e06,
        0.33842e06,
        0.36265e06,
        0.38839e06,
        0.41572e06,
        0.44474e06,
        0.47553e06,
        0.50818e06,
        0.54278e06,
        0.57945e06,
        0.61829e06,
        0.65940e06,
        0.70289e06,
        0.74890e06,
        0.79754e06,
        0.84894e06,
        0.90324e06,
        0.96057e06,
        0.10211e07,
        0.10849e07,
        0.11523e07,
        0.12233e07,
        0.12981e07,
        0.13769e07,
        0.14599e07,
        0.15473e07,
        0.16393e07,
        0.17361e07,
        0.18378e07,
        0.19447e07,
        0.20571e07,
        0.21752e07,
        0.22992e07,
        0.24294e07,
        0.25661e07,
        0.27094e07,
        0.28598e07,
        0.30175e07,
        0.31828e07,
        0.33560e07,
        0.35374e07,
        0.37274e07,
        0.39264e07,
        0.41346e07,
        0.43525e07,
        0.45805e07,
        0.48188e07,
        0.50681e07,
        0.53286e07,
        0.56008e07,
        0.58852e07,
        0.61823e07,
        0.64924e07,
        0.68162e07,
    ]
)


#  --------------- C2H2 1222: M = 26, I = 3 ---------------------
M = 26
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.24843e03,
        0.35373e03,
        0.45997e03,
        0.56930e03,
        0.68497e03,
        0.81065e03,
        0.94999e03,
        0.11065e04,
        0.12837e04,
        0.14848e04,
        0.17135e04,
        0.19731e04,
        0.22675e04,
        0.26205e04,
        0.29999e04,
        0.34276e04,
        0.39086e04,
        0.44486e04,
        0.50533e04,
        0.57294e04,
        0.64837e04,
        0.73237e04,
        0.82576e04,
        0.92941e04,
        0.10443e05,
        0.11714e05,
        0.13117e05,
        0.14666e05,
        0.16373e05,
        0.18250e05,
        0.20313e05,
        0.22578e05,
        0.25060e05,
        0.27777e05,
        0.30750e05,
        0.33997e05,
        0.37541e05,
        0.41405e05,
        0.45614e05,
        0.50192e05,
        0.55170e05,
        0.60576e05,
        0.66441e05,
        0.72799e05,
        0.79686e05,
        0.87140e05,
        0.95199e05,
        0.10391e06,
        0.11331e06,
        0.12345e06,
        0.13438e06,
        0.14615e06,
        0.15882e06,
        0.17245e06,
        0.18710e06,
        0.20283e06,
        0.21972e06,
        0.23783e06,
        0.25724e06,
        0.27804e06,
        0.30030e06,
        0.32411e06,
        0.34958e06,
        0.37679e06,
        0.40585e06,
        0.43686e06,
        0.46994e06,
        0.50521e06,
        0.54280e06,
        0.58282e06,
        0.62542e06,
        0.67074e06,
        0.71892e06,
        0.77013e06,
        0.82453e06,
        0.88228e06,
        0.94356e06,
        0.10086e07,
        0.10775e07,
        0.11505e07,
        0.12279e07,
        0.13098e07,
        0.13964e07,
        0.14881e07,
        0.15850e07,
        0.16875e07,
        0.17957e07,
        0.19100e07,
        0.20307e07,
        0.21580e07,
        0.22923e07,
        0.24339e07,
        0.25831e07,
        0.27404e07,
        0.29060e07,
        0.30803e07,
        0.32638e07,
        0.34568e07,
        0.36598e07,
        0.38733e07,
        0.40976e07,
        0.43332e07,
        0.45807e07,
        0.48406e07,
        0.51133e07,
        0.53995e07,
        0.56997e07,
        0.60144e07,
        0.63444e07,
        0.66901e07,
        0.70524e07,
        0.74317e07,
        0.78289e07,
        0.82447e07,
        0.86797e07,
        0.91348e07,
        0.96108e07,
        0.10108e08,
        0.10629e08,
    ]
)


#  --------------- C2H6 1221: M = 27, I = 1 ---------------------
M = 27
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.47267e04,
        0.80011e04,
        0.11928e05,
        0.16564e05,
        0.21985e05,
        0.28287e05,
        0.35590e05,
        0.44049e05,
        0.53862e05,
        0.65277e05,
        0.78597e05,
        0.94191e05,
        0.11250e06,
        0.13407e06,
        0.15952e06,
        0.18962e06,
        0.22526e06,
        0.26751e06,
        0.31763e06,
        0.37714e06,
        0.44780e06,
        0.53174e06,
        0.63145e06,
        0.74989e06,
        0.89056e06,
        0.10576e07,
        0.12559e07,
        0.14912e07,
        0.17704e07,
        0.21013e07,
        0.24936e07,
        0.29582e07,
        0.35083e07,
        0.41591e07,
        0.49286e07,
        0.58379e07,
        0.69116e07,
        0.81787e07,
        0.96728e07,
        0.11433e08,
        0.13506e08,
        0.15945e08,
        0.18812e08,
        0.22180e08,
        0.26134e08,
        0.30770e08,
        0.36204e08,
        0.42565e08,
        0.50008e08,
        0.58708e08,
        0.68868e08,
        0.80725e08,
        0.94548e08,
        0.11065e09,
        0.12940e09,
        0.15119e09,
        0.17652e09,
        0.20593e09,
        0.24003e09,
        0.27956e09,
        0.32533e09,
        0.37829e09,
        0.43951e09,
        0.51021e09,
        0.59180e09,
        0.68588e09,
        0.79427e09,
        0.91904e09,
        0.10625e10,
        0.12275e10,
        0.14168e10,
        0.16341e10,
        0.18831e10,
        0.21684e10,
        0.24949e10,
        0.28684e10,
        0.32951e10,
        0.37823e10,
        0.43382e10,
        0.49719e10,
        0.56938e10,
        0.65156e10,
        0.74502e10,
        0.85125e10,
        0.97190e10,
        0.11088e11,
        0.12641e11,
        0.14401e11,
        0.16393e11,
        0.18648e11,
        0.21198e11,
        0.24079e11,
        0.27332e11,
        0.31003e11,
        0.35142e11,
        0.39807e11,
        0.45060e11,
        0.50972e11,
        0.57620e11,
        0.65091e11,
        0.73483e11,
        0.82902e11,
        0.93467e11,
        0.10531e12,
        0.11858e12,
        0.13343e12,
        0.15005e12,
        0.16864e12,
        0.18941e12,
        0.21260e12,
        0.23849e12,
        0.26737e12,
        0.29957e12,
        0.33545e12,
        0.37541e12,
        0.41987e12,
        0.46934e12,
        0.52432e12,
        0.58542e12,
    ]
)


#  --------------- C2H6 1231: M = 27, I = 2 ---------------------
M = 27
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.24128e04,
        0.40845e04,
        0.60896e04,
        0.84564e04,
        0.11224e05,
        0.14442e05,
        0.18170e05,
        0.22490e05,
        0.27501e05,
        0.33329e05,
        0.40131e05,
        0.48094e05,
        0.57446e05,
        0.68459e05,
        0.81458e05,
        0.96828e05,
        0.11503e06,
        0.13661e06,
        0.16221e06,
        0.19260e06,
        0.22869e06,
        0.27156e06,
        0.32249e06,
        0.38298e06,
        0.45483e06,
        0.54015e06,
        0.64144e06,
        0.76164e06,
        0.90423e06,
        0.10733e07,
        0.12737e07,
        0.15110e07,
        0.17920e07,
        0.21245e07,
        0.25176e07,
        0.29821e07,
        0.35307e07,
        0.41780e07,
        0.49414e07,
        0.58408e07,
        0.68999e07,
        0.81461e07,
        0.96110e07,
        0.11332e08,
        0.13352e08,
        0.15721e08,
        0.18497e08,
        0.21748e08,
        0.25551e08,
        0.29997e08,
        0.35189e08,
        0.41248e08,
        0.48313e08,
        0.56542e08,
        0.66122e08,
        0.77262e08,
        0.90206e08,
        0.10523e09,
        0.12267e09,
        0.14287e09,
        0.16626e09,
        0.19333e09,
        0.22462e09,
        0.26076e09,
        0.30247e09,
        0.35056e09,
        0.40596e09,
        0.46974e09,
        0.54310e09,
        0.62740e09,
        0.72420e09,
        0.83527e09,
        0.96260e09,
        0.11084e10,
        0.12754e10,
        0.14663e10,
        0.16845e10,
        0.19336e10,
        0.22178e10,
        0.25418e10,
        0.29109e10,
        0.33311e10,
        0.38090e10,
        0.43522e10,
        0.49691e10,
        0.56693e10,
        0.64633e10,
        0.73631e10,
        0.83821e10,
        0.95352e10,
        0.10839e11,
        0.12312e11,
        0.13976e11,
        0.15854e11,
        0.17971e11,
        0.20357e11,
        0.23043e11,
        0.26067e11,
        0.29467e11,
        0.33289e11,
        0.37581e11,
        0.42399e11,
        0.47804e11,
        0.53862e11,
        0.60649e11,
        0.68247e11,
        0.76750e11,
        0.86257e11,
        0.96882e11,
        0.10875e12,
        0.12199e12,
        0.13677e12,
        0.15325e12,
        0.17160e12,
        0.19204e12,
        0.21480e12,
        0.24010e12,
        0.26824e12,
        0.29950e12,
    ]
)


#  --------------- PH3 1111: M = 28, I = 1 ---------------------
M = 28
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.29652e03,
        0.49643e03,
        0.72810e03,
        0.98777e03,
        0.12729e04,
        0.15820e04,
        0.19145e04,
        0.22708e04,
        0.26520e04,
        0.30600e04,
        0.34971e04,
        0.39662e04,
        0.44702e04,
        0.50126e04,
        0.55970e04,
        0.62273e04,
        0.69075e04,
        0.76421e04,
        0.84357e04,
        0.92933e04,
        0.10220e05,
        0.11222e05,
        0.12304e05,
        0.13473e05,
        0.14736e05,
        0.16099e05,
        0.17571e05,
        0.19160e05,
        0.20873e05,
        0.22720e05,
        0.24710e05,
        0.26854e05,
        0.29162e05,
        0.31646e05,
        0.34317e05,
        0.37188e05,
        0.40273e05,
        0.43585e05,
        0.47140e05,
        0.50953e05,
        0.55040e05,
        0.59419e05,
        0.64108e05,
        0.69127e05,
        0.74496e05,
        0.80236e05,
        0.86369e05,
        0.92918e05,
        0.99909e05,
        0.10737e06,
        0.11532e06,
        0.12380e06,
        0.13282e06,
        0.14244e06,
        0.15266e06,
        0.16354e06,
        0.17511e06,
        0.18739e06,
        0.20044e06,
        0.21430e06,
        0.22900e06,
        0.24459e06,
        0.26111e06,
        0.27862e06,
        0.29716e06,
        0.31680e06,
        0.33757e06,
        0.35954e06,
        0.38277e06,
        0.40733e06,
        0.43326e06,
        0.46065e06,
        0.48955e06,
        0.52005e06,
        0.55222e06,
        0.58614e06,
        0.62188e06,
        0.65953e06,
        0.69917e06,
        0.74091e06,
        0.78483e06,
        0.83103e06,
        0.87960e06,
        0.93067e06,
        0.98432e06,
        0.10407e07,
        0.10999e07,
        0.11620e07,
        0.12272e07,
        0.12956e07,
        0.13673e07,
        0.14425e07,
        0.15212e07,
        0.16038e07,
        0.16902e07,
        0.17808e07,
        0.18755e07,
        0.19746e07,
        0.20784e07,
        0.21868e07,
        0.23002e07,
        0.24187e07,
        0.25425e07,
        0.26719e07,
        0.28070e07,
        0.29480e07,
        0.30952e07,
        0.32488e07,
        0.34091e07,
        0.35762e07,
        0.37504e07,
        0.39320e07,
        0.41213e07,
        0.43185e07,
        0.45239e07,
        0.47378e07,
        0.49605e07,
        0.51923e07,
        0.54335e07,
    ]
)


#  --------------- COF2 269: M = 29, I = 1 ---------------------
M = 29
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.54999e04,
        0.92749e04,
        0.13668e05,
        0.18643e05,
        0.24224e05,
        0.30487e05,
        0.37547e05,
        0.45543e05,
        0.54639e05,
        0.65019e05,
        0.76886e05,
        0.90462e05,
        0.10600e06,
        0.12377e06,
        0.14407e06,
        0.16723e06,
        0.19363e06,
        0.22367e06,
        0.25780e06,
        0.29650e06,
        0.34031e06,
        0.38982e06,
        0.44568e06,
        0.50859e06,
        0.57932e06,
        0.65872e06,
        0.74770e06,
        0.84724e06,
        0.95844e06,
        0.10825e07,
        0.12205e07,
        0.13741e07,
        0.15446e07,
        0.17336e07,
        0.19428e07,
        0.21742e07,
        0.24296e07,
        0.27113e07,
        0.30214e07,
        0.33626e07,
        0.37373e07,
        0.41484e07,
        0.45989e07,
        0.50921e07,
        0.56313e07,
        0.62202e07,
        0.68626e07,
        0.75628e07,
        0.83251e07,
        0.91542e07,
        0.10055e08,
        0.11033e08,
        0.12093e08,
        0.13242e08,
        0.14486e08,
        0.15831e08,
        0.17284e08,
        0.18853e08,
        0.20546e08,
        0.22371e08,
        0.24335e08,
        0.26450e08,
        0.28724e08,
        0.31167e08,
        0.33790e08,
        0.36605e08,
        0.39623e08,
        0.42856e08,
        0.46318e08,
        0.50022e08,
        0.53983e08,
        0.58215e08,
        0.62735e08,
        0.67558e08,
        0.72702e08,
        0.78186e08,
        0.84028e08,
        0.90247e08,
        0.96865e08,
        0.10390e09,
        0.11138e09,
        0.11933e09,
        0.12777e09,
        0.13672e09,
        0.14622e09,
        0.15629e09,
        0.16695e09,
        0.17825e09,
        0.19021e09,
        0.20287e09,
        0.21625e09,
        0.23039e09,
        0.24534e09,
        0.26113e09,
        0.27779e09,
        0.29538e09,
        0.31392e09,
        0.33348e09,
        0.35409e09,
        0.37580e09,
        0.39867e09,
        0.42274e09,
        0.44806e09,
        0.47470e09,
        0.50271e09,
        0.53215e09,
        0.56308e09,
        0.59557e09,
        0.62968e09,
        0.66548e09,
        0.70304e09,
        0.74243e09,
        0.78374e09,
        0.82703e09,
        0.87240e09,
        0.91992e09,
        0.96967e09,
        0.10218e10,
        0.10763e10,
    ]
)


#  --------------- COF2 369: M = 29, I = 2 --------------------- not in TIPS-2011
M = 29
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- SF6 29: M = 30, I = 1 ---------------------
M = 30
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.46373e05,
        0.78844e05,
        0.11939e06,
        0.17183e06,
        0.24247e06,
        0.34059e06,
        0.47963e06,
        0.67906e06,
        0.96713e06,
        0.13848e07,
        0.19911e07,
        0.28714e07,
        0.41481e07,
        0.59956e07,
        0.86617e07,
        0.12496e08,
        0.17991e08,
        0.25832e08,
        0.36971e08,
        0.52724e08,
        0.74895e08,
        0.10595e09,
        0.14923e09,
        0.20925e09,
        0.29208e09,
        0.40582e09,
        0.56124e09,
        0.77259e09,
        0.10586e10,
        0.14439e10,
        0.19605e10,
        0.26500e10,
        0.35662e10,
        0.47781e10,
        0.63747e10,
        0.84689e10,
        0.11205e11,
        0.14765e11,
        0.19378e11,
        0.25336e11,
        0.32998e11,
        0.42819e11,
        0.55361e11,
        0.71323e11,
        0.91569e11,
        0.11716e12,
        0.14941e12,
        0.18992e12,
        0.24065e12,
        0.30398e12,
        0.38283e12,
        0.48069e12,
        0.60182e12,
        0.75136e12,
        0.93546e12,
        0.11615e13,
        0.14384e13,
        0.17767e13,
        0.21890e13,
        0.26903e13,
        0.32984e13,
        0.40344e13,
        0.49232e13,
        0.59942e13,
        0.72819e13,
        0.88272e13,
        0.10678e14,
        0.12889e14,
        0.15527e14,
        0.18666e14,
        0.22397e14,
        0.26823e14,
        0.32062e14,
        0.38253e14,
        0.45558e14,
        0.54161e14,
        0.64277e14,
        0.76153e14,
        0.90072e14,
        0.10636e15,
        0.12539e15,
        0.14759e15,
        0.17345e15,
        0.20354e15,
        0.23848e15,
        0.27902e15,
        0.32597e15,
        0.38028e15,
        0.44303e15,
        0.51542e15,
        0.59883e15,
        0.69482e15,
        0.80516e15,
        0.93182e15,
        0.10770e16,
        0.12434e16,
        0.14336e16,
        0.16511e16,
        0.18992e16,
        0.21821e16,
        0.25043e16,
        0.28709e16,
        0.32875e16,
        0.37604e16,
        0.42968e16,
        0.49046e16,
        0.55925e16,
        0.63704e16,
        0.72492e16,
        0.82411e16,
        0.93596e16,
        0.10620e17,
        0.12038e17,
        0.13633e17,
        0.15425e17,
        0.17438e17,
        0.19694e17,
        0.22224e17,
        0.25057e17,
    ]
)


#  --------------- H2S 121: M = 31, I = 1 ---------------------
M = 31
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.47192e02,
        0.78671e02,
        0.11510e03,
        0.15589e03,
        0.20061e03,
        0.24896e03,
        0.30070e03,
        0.35571e03,
        0.41386e03,
        0.47513e03,
        0.53951e03,
        0.60703e03,
        0.67772e03,
        0.75167e03,
        0.82896e03,
        0.90969e03,
        0.99396e03,
        0.10819e04,
        0.11736e04,
        0.12692e04,
        0.13689e04,
        0.14727e04,
        0.15809e04,
        0.16937e04,
        0.18111e04,
        0.19333e04,
        0.20606e04,
        0.21931e04,
        0.23309e04,
        0.24744e04,
        0.26236e04,
        0.27788e04,
        0.29403e04,
        0.31081e04,
        0.32825e04,
        0.34638e04,
        0.36522e04,
        0.38478e04,
        0.40510e04,
        0.42619e04,
        0.44808e04,
        0.47080e04,
        0.49437e04,
        0.51881e04,
        0.54415e04,
        0.57042e04,
        0.59764e04,
        0.62584e04,
        0.65505e04,
        0.68529e04,
        0.71660e04,
        0.74899e04,
        0.78251e04,
        0.81718e04,
        0.85303e04,
        0.89008e04,
        0.92838e04,
        0.96795e04,
        0.10088e05,
        0.10510e05,
        0.10946e05,
        0.11396e05,
        0.11860e05,
        0.12339e05,
        0.12833e05,
        0.13342e05,
        0.13867e05,
        0.14408e05,
        0.14966e05,
        0.15540e05,
        0.16132e05,
        0.16741e05,
        0.17368e05,
        0.18013e05,
        0.18677e05,
        0.19361e05,
        0.20064e05,
        0.20786e05,
        0.21529e05,
        0.22293e05,
        0.23078e05,
        0.23885e05,
        0.24714e05,
        0.25565e05,
        0.26439e05,
        0.27337e05,
        0.28258e05,
        0.29204e05,
        0.30174e05,
        0.31170e05,
        0.32191e05,
        0.33239e05,
        0.34313e05,
        0.35414e05,
        0.36543e05,
        0.37700e05,
        0.38886e05,
        0.40101e05,
        0.41346e05,
        0.42621e05,
        0.43926e05,
        0.45263e05,
        0.46631e05,
        0.48033e05,
        0.49466e05,
        0.50934e05,
        0.52435e05,
        0.53971e05,
        0.55542e05,
        0.57149e05,
        0.58792e05,
        0.60472e05,
        0.62190e05,
        0.63946e05,
        0.65740e05,
        0.67574e05,
        0.69448e05,
        0.71362e05,
        0.73318e05,
    ]
)


#  --------------- H2S 141: M = 31, I = 2 ---------------------
M = 31
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.47310e02,
        0.78869e02,
        0.11539e03,
        0.15628e03,
        0.20112e03,
        0.24959e03,
        0.30147e03,
        0.35661e03,
        0.41491e03,
        0.47634e03,
        0.54088e03,
        0.60857e03,
        0.67945e03,
        0.75359e03,
        0.83107e03,
        0.91201e03,
        0.99649e03,
        0.10846e04,
        0.11766e04,
        0.12724e04,
        0.13724e04,
        0.14765e04,
        0.15850e04,
        0.16980e04,
        0.18157e04,
        0.19382e04,
        0.20658e04,
        0.21987e04,
        0.23369e04,
        0.24807e04,
        0.26303e04,
        0.27859e04,
        0.29478e04,
        0.31160e04,
        0.32909e04,
        0.34727e04,
        0.36615e04,
        0.38576e04,
        0.40613e04,
        0.42728e04,
        0.44923e04,
        0.47200e04,
        0.49563e04,
        0.52013e04,
        0.54554e04,
        0.57188e04,
        0.59917e04,
        0.62744e04,
        0.65672e04,
        0.68704e04,
        0.71843e04,
        0.75090e04,
        0.78451e04,
        0.81926e04,
        0.85520e04,
        0.89236e04,
        0.93075e04,
        0.97042e04,
        0.10114e05,
        0.10537e05,
        0.10974e05,
        0.11425e05,
        0.11890e05,
        0.12370e05,
        0.12866e05,
        0.13376e05,
        0.13903e05,
        0.14445e05,
        0.15004e05,
        0.15580e05,
        0.16173e05,
        0.16784e05,
        0.17412e05,
        0.18059e05,
        0.18725e05,
        0.19410e05,
        0.20115e05,
        0.20839e05,
        0.21584e05,
        0.22350e05,
        0.23137e05,
        0.23946e05,
        0.24777e05,
        0.25630e05,
        0.26507e05,
        0.27407e05,
        0.28330e05,
        0.29278e05,
        0.30251e05,
        0.31249e05,
        0.32273e05,
        0.33324e05,
        0.34401e05,
        0.35505e05,
        0.36637e05,
        0.37797e05,
        0.38985e05,
        0.40204e05,
        0.41451e05,
        0.42729e05,
        0.44038e05,
        0.45379e05,
        0.46751e05,
        0.48155e05,
        0.49593e05,
        0.51064e05,
        0.52569e05,
        0.54109e05,
        0.55684e05,
        0.57295e05,
        0.58943e05,
        0.60627e05,
        0.62349e05,
        0.64109e05,
        0.65908e05,
        0.67747e05,
        0.69625e05,
        0.71544e05,
        0.73505e05,
    ]
)


#  --------------- H2S 131: M = 30, I = 3 ---------------------
M = 31
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.18901e03,
        0.31509e03,
        0.46102e03,
        0.62437e03,
        0.80349e03,
        0.99713e03,
        0.12044e04,
        0.14247e04,
        0.16576e04,
        0.19030e04,
        0.21609e04,
        0.24313e04,
        0.27145e04,
        0.30106e04,
        0.33202e04,
        0.36436e04,
        0.39811e04,
        0.43332e04,
        0.47005e04,
        0.50835e04,
        0.54827e04,
        0.58987e04,
        0.63321e04,
        0.67836e04,
        0.72538e04,
        0.77434e04,
        0.82532e04,
        0.87838e04,
        0.93360e04,
        0.99106e04,
        0.10508e05,
        0.11130e05,
        0.11777e05,
        0.12449e05,
        0.13147e05,
        0.13874e05,
        0.14628e05,
        0.15412e05,
        0.16225e05,
        0.17070e05,
        0.17947e05,
        0.18857e05,
        0.19801e05,
        0.20780e05,
        0.21795e05,
        0.22847e05,
        0.23937e05,
        0.25067e05,
        0.26236e05,
        0.27448e05,
        0.28702e05,
        0.29999e05,
        0.31342e05,
        0.32730e05,
        0.34166e05,
        0.35650e05,
        0.37184e05,
        0.38769e05,
        0.40406e05,
        0.42097e05,
        0.43842e05,
        0.45644e05,
        0.47503e05,
        0.49421e05,
        0.51399e05,
        0.53439e05,
        0.55542e05,
        0.57709e05,
        0.59942e05,
        0.62242e05,
        0.64611e05,
        0.67051e05,
        0.69563e05,
        0.72148e05,
        0.74808e05,
        0.77545e05,
        0.80360e05,
        0.83255e05,
        0.86232e05,
        0.89291e05,
        0.92435e05,
        0.95667e05,
        0.98986e05,
        0.10240e06,
        0.10590e06,
        0.10949e06,
        0.11318e06,
        0.11697e06,
        0.12086e06,
        0.12484e06,
        0.12893e06,
        0.13313e06,
        0.13743e06,
        0.14184e06,
        0.14637e06,
        0.15100e06,
        0.15575e06,
        0.16062e06,
        0.16560e06,
        0.17071e06,
        0.17594e06,
        0.18129e06,
        0.18677e06,
        0.19238e06,
        0.19813e06,
        0.20400e06,
        0.21002e06,
        0.21617e06,
        0.22246e06,
        0.22890e06,
        0.23548e06,
        0.24221e06,
        0.24909e06,
        0.25612e06,
        0.26331e06,
        0.27065e06,
        0.27816e06,
        0.28583e06,
        0.29366e06,
    ]
)


#  --------------- HCOOH 126: M = 32, I = 1 ---------------------
M = 32
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.31899e04,
        0.53773e04,
        0.79205e04,
        0.10792e05,
        0.13993e05,
        0.17550e05,
        0.21509e05,
        0.25930e05,
        0.30885e05,
        0.36460e05,
        0.42750e05,
        0.49864e05,
        0.57926e05,
        0.67071e05,
        0.77453e05,
        0.89243e05,
        0.10263e06,
        0.11783e06,
        0.13507e06,
        0.15462e06,
        0.17676e06,
        0.20183e06,
        0.23018e06,
        0.26221e06,
        0.29836e06,
        0.33911e06,
        0.38501e06,
        0.43664e06,
        0.49467e06,
        0.55981e06,
        0.63286e06,
        0.71470e06,
        0.80628e06,
        0.90865e06,
        0.10230e07,
        0.11505e07,
        0.12927e07,
        0.14509e07,
        0.16269e07,
        0.18225e07,
        0.20396e07,
        0.22804e07,
        0.25472e07,
        0.28425e07,
        0.31692e07,
        0.35301e07,
        0.39285e07,
        0.43681e07,
        0.48525e07,
        0.53858e07,
        0.59727e07,
        0.66178e07,
        0.73265e07,
        0.81042e07,
        0.89571e07,
        0.98918e07,
        0.10915e08,
        0.12035e08,
        0.13259e08,
        0.14597e08,
        0.16057e08,
        0.17650e08,
        0.19387e08,
        0.21279e08,
        0.23339e08,
        0.25579e08,
        0.28016e08,
        0.30663e08,
        0.33536e08,
        0.36655e08,
        0.40037e08,
        0.43701e08,
        0.47671e08,
        0.51967e08,
        0.56614e08,
        0.61639e08,
        0.67068e08,
        0.72930e08,
        0.79257e08,
        0.86082e08,
        0.93439e08,
        0.10137e09,
        0.10990e09,
        0.11909e09,
        0.12898e09,
        0.13960e09,
        0.15102e09,
        0.16329e09,
        0.17646e09,
        0.19059e09,
        0.20575e09,
        0.22200e09,
        0.23941e09,
        0.25806e09,
        0.27802e09,
        0.29938e09,
        0.32223e09,
        0.34666e09,
        0.37276e09,
        0.40064e09,
        0.43041e09,
        0.46218e09,
        0.49607e09,
        0.53221e09,
        0.57074e09,
        0.61179e09,
        0.65551e09,
        0.70206e09,
        0.75159e09,
        0.80430e09,
        0.86034e09,
        0.91992e09,
        0.98324e09,
        0.10505e10,
        0.11219e10,
        0.11977e10,
        0.12782e10,
        0.13635e10,
        0.14540e10,
    ]
)


#  --------------- HO2 166: M = 33, I = 1 ---------------------
M = 33
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.39277e03,
        0.66062e03,
        0.97123e03,
        0.13194e04,
        0.17014e04,
        0.21148e04,
        0.25578e04,
        0.30296e04,
        0.35297e04,
        0.40585e04,
        0.46167e04,
        0.52055e04,
        0.58264e04,
        0.64809e04,
        0.71707e04,
        0.78978e04,
        0.86641e04,
        0.94715e04,
        0.10322e05,
        0.11218e05,
        0.12161e05,
        0.13154e05,
        0.14198e05,
        0.15296e05,
        0.16449e05,
        0.17661e05,
        0.18933e05,
        0.20267e05,
        0.21666e05,
        0.23133e05,
        0.24669e05,
        0.26277e05,
        0.27960e05,
        0.29720e05,
        0.31560e05,
        0.33482e05,
        0.35489e05,
        0.37584e05,
        0.39769e05,
        0.42048e05,
        0.44423e05,
        0.46898e05,
        0.49475e05,
        0.52157e05,
        0.54948e05,
        0.57850e05,
        0.60868e05,
        0.64003e05,
        0.67261e05,
        0.70643e05,
        0.74154e05,
        0.77797e05,
        0.81575e05,
        0.85492e05,
        0.89553e05,
        0.93760e05,
        0.98118e05,
        0.10263e06,
        0.10730e06,
        0.11213e06,
        0.11713e06,
        0.12230e06,
        0.12765e06,
        0.13317e06,
        0.13888e06,
        0.14478e06,
        0.15086e06,
        0.15715e06,
        0.16363e06,
        0.17032e06,
        0.17723e06,
        0.18434e06,
        0.19168e06,
        0.19924e06,
        0.20704e06,
        0.21506e06,
        0.22333e06,
        0.23185e06,
        0.24061e06,
        0.24963e06,
        0.25891e06,
        0.26846e06,
        0.27828e06,
        0.28838e06,
        0.29876e06,
        0.30943e06,
        0.32039e06,
        0.33166e06,
        0.34323e06,
        0.35512e06,
        0.36732e06,
        0.37985e06,
        0.39271e06,
        0.40590e06,
        0.41944e06,
        0.43333e06,
        0.44758e06,
        0.46219e06,
        0.47717e06,
        0.49252e06,
        0.50826e06,
        0.52439e06,
        0.54091e06,
        0.55784e06,
        0.57518e06,
        0.59293e06,
        0.61112e06,
        0.62973e06,
        0.64878e06,
        0.66828e06,
        0.68824e06,
        0.70866e06,
        0.72955e06,
        0.75091e06,
        0.77276e06,
        0.79511e06,
        0.81795e06,
        0.84131e06,
        0.86518e06,
    ]
)


#  --------------- O 6: M = 34, I = 1 --------------------- not in TIPS-2011
M = 34
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- ClONO2 5646: M = 35, I = 1 ---------------------
M = 35
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.11444e06,
        0.21121e06,
        0.34858e06,
        0.53934e06,
        0.80041e06,
        0.11539e07,
        0.16286e07,
        0.22614e07,
        0.30992e07,
        0.42015e07,
        0.56426e07,
        0.75152e07,
        0.99344e07,
        0.13042e08,
        0.17012e08,
        0.22058e08,
        0.28437e08,
        0.36463e08,
        0.46514e08,
        0.59042e08,
        0.74589e08,
        0.93801e08,
        0.11744e09,
        0.14643e09,
        0.18181e09,
        0.22486e09,
        0.27705e09,
        0.34009e09,
        0.41598e09,
        0.50705e09,
        0.61599e09,
        0.74590e09,
        0.90037e09,
        0.10835e10,
        0.13001e10,
        0.15554e10,
        0.18556e10,
        0.22079e10,
        0.26200e10,
        0.31012e10,
        0.36615e10,
        0.43126e10,
        0.50675e10,
        0.59409e10,
        0.69492e10,
        0.81110e10,
        0.94469e10,
        0.10980e11,
        0.12736e11,
        0.14745e11,
        0.17037e11,
        0.19649e11,
        0.22620e11,
        0.25994e11,
        0.29819e11,
        0.34150e11,
        0.39044e11,
        0.44568e11,
        0.50794e11,
        0.57799e11,
        0.65672e11,
        0.74506e11,
        0.84408e11,
        0.95490e11,
        0.10788e12,
        0.12171e12,
        0.13713e12,
        0.15431e12,
        0.17342e12,
        0.19465e12,
        0.21822e12,
        0.24435e12,
        0.27329e12,
        0.30530e12,
        0.34069e12,
        0.37976e12,
        0.42286e12,
        0.47034e12,
        0.52262e12,
        0.58012e12,
        0.64330e12,
        0.71267e12,
        0.78875e12,
        0.87214e12,
        0.96344e12,
        0.10633e13,
        0.11725e13,
        0.12918e13,
        0.14220e13,
        0.15640e13,
        0.17188e13,
        0.18873e13,
        0.20706e13,
        0.22700e13,
        0.24866e13,
        0.27218e13,
        0.29771e13,
        0.32538e13,
        0.35537e13,
        0.38784e13,
        0.42299e13,
        0.46100e13,
        0.50208e13,
        0.54645e13,
        0.59435e13,
        0.64603e13,
        0.70175e13,
        0.76180e13,
        0.82647e13,
        0.89608e13,
        0.97097e13,
        0.10515e14,
        0.11380e14,
        0.12310e14,
        0.13307e14,
        0.14378e14,
        0.15526e14,
        0.16756e14,
        0.18075e14,
    ]
)


#  --------------- ClONO2 7646: M = 35, I = 2 ---------------------
M = 35
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.11735e06,
        0.21659e06,
        0.35745e06,
        0.55307e06,
        0.82078e06,
        0.11833e07,
        0.16700e07,
        0.23189e07,
        0.31781e07,
        0.43084e07,
        0.57862e07,
        0.77065e07,
        0.10187e08,
        0.13374e08,
        0.17445e08,
        0.22619e08,
        0.29161e08,
        0.37391e08,
        0.47698e08,
        0.60545e08,
        0.76487e08,
        0.96188e08,
        0.12043e09,
        0.15015e09,
        0.18644e09,
        0.23059e09,
        0.28410e09,
        0.34874e09,
        0.42657e09,
        0.51995e09,
        0.63167e09,
        0.76489e09,
        0.92329e09,
        0.11111e10,
        0.13331e10,
        0.15950e10,
        0.19029e10,
        0.22641e10,
        0.26867e10,
        0.31801e10,
        0.37547e10,
        0.44224e10,
        0.51965e10,
        0.60921e10,
        0.71261e10,
        0.83174e10,
        0.96873e10,
        0.11260e11,
        0.13061e11,
        0.15120e11,
        0.17471e11,
        0.20149e11,
        0.23196e11,
        0.26656e11,
        0.30578e11,
        0.35019e11,
        0.40038e11,
        0.45703e11,
        0.52087e11,
        0.59270e11,
        0.67343e11,
        0.76403e11,
        0.86556e11,
        0.97921e11,
        0.11062e12,
        0.12481e12,
        0.14062e12,
        0.15824e12,
        0.17783e12,
        0.19961e12,
        0.22377e12,
        0.25057e12,
        0.28024e12,
        0.31308e12,
        0.34936e12,
        0.38943e12,
        0.43362e12,
        0.48232e12,
        0.53593e12,
        0.59489e12,
        0.65968e12,
        0.73081e12,
        0.80883e12,
        0.89434e12,
        0.98797e12,
        0.10904e13,
        0.12024e13,
        0.13247e13,
        0.14582e13,
        0.16038e13,
        0.17625e13,
        0.19353e13,
        0.21233e13,
        0.23278e13,
        0.25499e13,
        0.27911e13,
        0.30528e13,
        0.33366e13,
        0.36442e13,
        0.39772e13,
        0.43376e13,
        0.47273e13,
        0.51486e13,
        0.56036e13,
        0.60948e13,
        0.66248e13,
        0.71962e13,
        0.78119e13,
        0.84751e13,
        0.91889e13,
        0.99569e13,
        0.10783e14,
        0.11670e14,
        0.12623e14,
        0.13646e14,
        0.14744e14,
        0.15921e14,
        0.17183e14,
        0.18535e14,
    ]
)


#  --------------- NOp 46: M = 36, I = 1 ---------------------
M = 36
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.63956e02,
        0.90185e02,
        0.11642e03,
        0.14265e03,
        0.16889e03,
        0.19513e03,
        0.22138e03,
        0.24763e03,
        0.27388e03,
        0.30013e03,
        0.32639e03,
        0.35266e03,
        0.37894e03,
        0.40523e03,
        0.43155e03,
        0.45790e03,
        0.48429e03,
        0.51074e03,
        0.53725e03,
        0.56383e03,
        0.59052e03,
        0.61731e03,
        0.64422e03,
        0.67127e03,
        0.69846e03,
        0.72582e03,
        0.75335e03,
        0.78108e03,
        0.80901e03,
        0.83715e03,
        0.86552e03,
        0.89413e03,
        0.92298e03,
        0.95208e03,
        0.98144e03,
        0.10111e04,
        0.10410e04,
        0.10712e04,
        0.11017e04,
        0.11325e04,
        0.11636e04,
        0.11950e04,
        0.12268e04,
        0.12588e04,
        0.12912e04,
        0.13239e04,
        0.13570e04,
        0.13903e04,
        0.14241e04,
        0.14581e04,
        0.14926e04,
        0.15273e04,
        0.15624e04,
        0.15979e04,
        0.16337e04,
        0.16699e04,
        0.17065e04,
        0.17434e04,
        0.17806e04,
        0.18183e04,
        0.18563e04,
        0.18947e04,
        0.19334e04,
        0.19725e04,
        0.20120e04,
        0.20519e04,
        0.20921e04,
        0.21327e04,
        0.21737e04,
        0.22151e04,
        0.22568e04,
        0.22990e04,
        0.23415e04,
        0.23844e04,
        0.24276e04,
        0.24713e04,
        0.25153e04,
        0.25598e04,
        0.26046e04,
        0.26497e04,
        0.26953e04,
        0.27413e04,
        0.27876e04,
        0.28343e04,
        0.28815e04,
        0.29290e04,
        0.29769e04,
        0.30251e04,
        0.30738e04,
        0.31229e04,
        0.31723e04,
        0.32222e04,
        0.32724e04,
        0.33230e04,
        0.33740e04,
        0.34254e04,
        0.34772e04,
        0.35294e04,
        0.35819e04,
        0.36349e04,
        0.36883e04,
        0.37420e04,
        0.37961e04,
        0.38507e04,
        0.39056e04,
        0.39609e04,
        0.40166e04,
        0.40727e04,
        0.41292e04,
        0.41861e04,
        0.42434e04,
        0.43010e04,
        0.43591e04,
        0.44176e04,
        0.44764e04,
        0.45357e04,
        0.45953e04,
        0.46554e04,
        0.47158e04,
    ]
)


#  --------------- HOBr 169: M = 37, I = 1 ---------------------
M = 37
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.24445e04,
        0.41206e04,
        0.60683e04,
        0.82610e04,
        0.10689e05,
        0.13352e05,
        0.16261e05,
        0.19427e05,
        0.22867e05,
        0.26600e05,
        0.30643e05,
        0.35018e05,
        0.39745e05,
        0.44844e05,
        0.50338e05,
        0.56249e05,
        0.62599e05,
        0.69410e05,
        0.76706e05,
        0.84509e05,
        0.92845e05,
        0.10174e06,
        0.11121e06,
        0.12128e06,
        0.13199e06,
        0.14335e06,
        0.15540e06,
        0.16815e06,
        0.18165e06,
        0.19591e06,
        0.21096e06,
        0.22684e06,
        0.24358e06,
        0.26120e06,
        0.27974e06,
        0.29922e06,
        0.31969e06,
        0.34118e06,
        0.36372e06,
        0.38735e06,
        0.41210e06,
        0.43800e06,
        0.46511e06,
        0.49345e06,
        0.52307e06,
        0.55400e06,
        0.58628e06,
        0.61997e06,
        0.65509e06,
        0.69170e06,
        0.72984e06,
        0.76954e06,
        0.81087e06,
        0.85386e06,
        0.89856e06,
        0.94502e06,
        0.99329e06,
        0.10434e07,
        0.10955e07,
        0.11495e07,
        0.12055e07,
        0.12636e07,
        0.13238e07,
        0.13862e07,
        0.14508e07,
        0.15177e07,
        0.15870e07,
        0.16587e07,
        0.17328e07,
        0.18095e07,
        0.18888e07,
        0.19707e07,
        0.20554e07,
        0.21428e07,
        0.22331e07,
        0.23263e07,
        0.24225e07,
        0.25217e07,
        0.26241e07,
        0.27296e07,
        0.28385e07,
        0.29506e07,
        0.30662e07,
        0.31853e07,
        0.33079e07,
        0.34341e07,
        0.35641e07,
        0.36979e07,
        0.38355e07,
        0.39771e07,
        0.41228e07,
        0.42725e07,
        0.44265e07,
        0.45848e07,
        0.47474e07,
        0.49145e07,
        0.50862e07,
        0.52624e07,
        0.54435e07,
        0.56293e07,
        0.58201e07,
        0.60159e07,
        0.62168e07,
        0.64229e07,
        0.66343e07,
        0.68511e07,
        0.70734e07,
        0.73013e07,
        0.75349e07,
        0.77742e07,
        0.80196e07,
        0.82709e07,
        0.85283e07,
        0.87920e07,
        0.90620e07,
        0.93385e07,
        0.96215e07,
        0.99112e07,
        0.10208e08,
    ]
)


#  --------------- HOBr 161: M = 37, I = 2 ---------------------
M = 37
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.24350e04,
        0.41047e04,
        0.60448e04,
        0.82291e04,
        0.10648e05,
        0.13301e05,
        0.16200e05,
        0.19355e05,
        0.22784e05,
        0.26504e05,
        0.30534e05,
        0.34895e05,
        0.39607e05,
        0.44691e05,
        0.50169e05,
        0.56063e05,
        0.62394e05,
        0.69186e05,
        0.76461e05,
        0.84243e05,
        0.92555e05,
        0.10142e06,
        0.11087e06,
        0.12091e06,
        0.13159e06,
        0.14292e06,
        0.15494e06,
        0.16766e06,
        0.18112e06,
        0.19534e06,
        0.21036e06,
        0.22620e06,
        0.24289e06,
        0.26047e06,
        0.27896e06,
        0.29840e06,
        0.31882e06,
        0.34025e06,
        0.36274e06,
        0.38630e06,
        0.41099e06,
        0.43683e06,
        0.46387e06,
        0.49215e06,
        0.52169e06,
        0.55255e06,
        0.58475e06,
        0.61836e06,
        0.65340e06,
        0.68992e06,
        0.72796e06,
        0.76757e06,
        0.80880e06,
        0.85169e06,
        0.89628e06,
        0.94263e06,
        0.99079e06,
        0.10408e07,
        0.10927e07,
        0.11466e07,
        0.12025e07,
        0.12605e07,
        0.13205e07,
        0.13828e07,
        0.14472e07,
        0.15140e07,
        0.15831e07,
        0.16546e07,
        0.17286e07,
        0.18051e07,
        0.18842e07,
        0.19660e07,
        0.20504e07,
        0.21377e07,
        0.22277e07,
        0.23207e07,
        0.24167e07,
        0.25157e07,
        0.26178e07,
        0.27231e07,
        0.28317e07,
        0.29436e07,
        0.30589e07,
        0.31777e07,
        0.33001e07,
        0.34260e07,
        0.35557e07,
        0.36892e07,
        0.38265e07,
        0.39678e07,
        0.41131e07,
        0.42626e07,
        0.44162e07,
        0.45741e07,
        0.47364e07,
        0.49031e07,
        0.50744e07,
        0.52503e07,
        0.54309e07,
        0.56164e07,
        0.58067e07,
        0.60021e07,
        0.62025e07,
        0.64081e07,
        0.66191e07,
        0.68354e07,
        0.70572e07,
        0.72846e07,
        0.75177e07,
        0.77565e07,
        0.80013e07,
        0.82521e07,
        0.85090e07,
        0.87721e07,
        0.90415e07,
        0.93173e07,
        0.95997e07,
        0.98888e07,
        0.10185e08,
    ]
)


#  --------------- C2H4 221: M = 38, I = 1 ---------------------
M = 38
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.95843e03,
        0.16137e04,
        0.23744e04,
        0.32285e04,
        0.41694e04,
        0.51963e04,
        0.63143e04,
        0.75337e04,
        0.88702e04,
        0.10344e05,
        0.11978e05,
        0.13802e05,
        0.15846e05,
        0.18145e05,
        0.20740e05,
        0.23675e05,
        0.27000e05,
        0.30770e05,
        0.35048e05,
        0.39905e05,
        0.45420e05,
        0.51680e05,
        0.58786e05,
        0.66850e05,
        0.75997e05,
        0.86369e05,
        0.98123e05,
        0.11144e06,
        0.12651e06,
        0.14356e06,
        0.16284e06,
        0.18463e06,
        0.20923e06,
        0.23699e06,
        0.26831e06,
        0.30360e06,
        0.34334e06,
        0.38808e06,
        0.43840e06,
        0.49495e06,
        0.55847e06,
        0.62976e06,
        0.70973e06,
        0.79935e06,
        0.89973e06,
        0.10121e07,
        0.11378e07,
        0.12782e07,
        0.14351e07,
        0.16102e07,
        0.18055e07,
        0.20231e07,
        0.22656e07,
        0.25354e07,
        0.28356e07,
        0.31692e07,
        0.35398e07,
        0.39511e07,
        0.44074e07,
        0.49132e07,
        0.54736e07,
        0.60940e07,
        0.67803e07,
        0.75392e07,
        0.83776e07,
        0.93035e07,
        0.10325e08,
        0.11452e08,
        0.12694e08,
        0.14062e08,
        0.15567e08,
        0.17224e08,
        0.19045e08,
        0.21046e08,
        0.23243e08,
        0.25655e08,
        0.28300e08,
        0.31200e08,
        0.34377e08,
        0.37856e08,
        0.41662e08,
        0.45826e08,
        0.50378e08,
        0.55351e08,
        0.60781e08,
        0.66707e08,
        0.73172e08,
        0.80219e08,
        0.87899e08,
        0.96262e08,
        0.10537e09,
        0.11527e09,
        0.12604e09,
        0.13775e09,
        0.15047e09,
        0.16428e09,
        0.17927e09,
        0.19553e09,
        0.21316e09,
        0.23226e09,
        0.25296e09,
        0.27537e09,
        0.29963e09,
        0.32587e09,
        0.35425e09,
        0.38492e09,
        0.41805e09,
        0.45383e09,
        0.49246e09,
        0.53413e09,
        0.57908e09,
        0.62754e09,
        0.67977e09,
        0.73602e09,
        0.79660e09,
        0.86179e09,
        0.93194e09,
        0.10074e10,
        0.10885e10,
    ]
)


#  --------------- C2H4 231: M = 38, I = 2 ---------------------
M = 38
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.39228e04,
        0.66051e04,
        0.97190e04,
        0.13215e05,
        0.17066e05,
        0.21270e05,
        0.25846e05,
        0.30838e05,
        0.36309e05,
        0.42341e05,
        0.49032e05,
        0.56496e05,
        0.64862e05,
        0.74275e05,
        0.84897e05,
        0.96912e05,
        0.11052e06,
        0.12595e06,
        0.14347e06,
        0.16335e06,
        0.18592e06,
        0.21155e06,
        0.24064e06,
        0.27365e06,
        0.31109e06,
        0.35354e06,
        0.40166e06,
        0.45615e06,
        0.51785e06,
        0.58765e06,
        0.66657e06,
        0.75575e06,
        0.85646e06,
        0.97011e06,
        0.10983e07,
        0.12428e07,
        0.14055e07,
        0.15886e07,
        0.17945e07,
        0.20260e07,
        0.22861e07,
        0.25779e07,
        0.29052e07,
        0.32721e07,
        0.36830e07,
        0.41429e07,
        0.46573e07,
        0.52323e07,
        0.58744e07,
        0.65912e07,
        0.73906e07,
        0.82816e07,
        0.92740e07,
        0.10379e08,
        0.11607e08,
        0.12973e08,
        0.14490e08,
        0.16174e08,
        0.18042e08,
        0.20112e08,
        0.22406e08,
        0.24945e08,
        0.27755e08,
        0.30861e08,
        0.34293e08,
        0.38083e08,
        0.42266e08,
        0.46878e08,
        0.51961e08,
        0.57560e08,
        0.63724e08,
        0.70504e08,
        0.77959e08,
        0.86150e08,
        0.95145e08,
        0.10502e09,
        0.11585e09,
        0.12772e09,
        0.14072e09,
        0.15496e09,
        0.17054e09,
        0.18759e09,
        0.20622e09,
        0.22658e09,
        0.24880e09,
        0.27306e09,
        0.29952e09,
        0.32837e09,
        0.35981e09,
        0.39404e09,
        0.43131e09,
        0.47186e09,
        0.51595e09,
        0.56387e09,
        0.61594e09,
        0.67247e09,
        0.73382e09,
        0.80038e09,
        0.87255e09,
        0.95076e09,
        0.10355e10,
        0.11272e10,
        0.12265e10,
        0.13339e10,
        0.14501e10,
        0.15756e10,
        0.17113e10,
        0.18577e10,
        0.20159e10,
        0.21865e10,
        0.23705e10,
        0.25688e10,
        0.27826e10,
        0.30129e10,
        0.32608e10,
        0.35277e10,
        0.38149e10,
        0.41237e10,
        0.44557e10,
    ]
)


#  --------------- CH3OH 2161: M = 39, I = 1 --------------------- not in TIPS-2011
M = 39
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


#  --------------- CH3Br 219: M = 40, I = 1 ---------------------
M = 40
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.70299e04,
        0.11847e05,
        0.17442e05,
        0.23741e05,
        0.30723e05,
        0.38408e05,
        0.46851e05,
        0.56138e05,
        0.66375e05,
        0.77692e05,
        0.90239e05,
        0.10418e06,
        0.11972e06,
        0.13704e06,
        0.15639e06,
        0.17801e06,
        0.20218e06,
        0.22920e06,
        0.25940e06,
        0.29316e06,
        0.33087e06,
        0.37296e06,
        0.41992e06,
        0.47229e06,
        0.53062e06,
        0.59557e06,
        0.66781e06,
        0.74812e06,
        0.83731e06,
        0.93629e06,
        0.10461e07,
        0.11677e07,
        0.13023e07,
        0.14513e07,
        0.16159e07,
        0.17978e07,
        0.19985e07,
        0.22199e07,
        0.24638e07,
        0.27324e07,
        0.30280e07,
        0.33529e07,
        0.37099e07,
        0.41019e07,
        0.45319e07,
        0.50034e07,
        0.55199e07,
        0.60853e07,
        0.67039e07,
        0.73801e07,
        0.81189e07,
        0.89255e07,
        0.98056e07,
        0.10765e08,
        0.11811e08,
        0.12949e08,
        0.14188e08,
        0.15535e08,
        0.17000e08,
        0.18590e08,
        0.20317e08,
        0.22190e08,
        0.24220e08,
        0.26421e08,
        0.28804e08,
        0.31383e08,
        0.34173e08,
        0.37189e08,
        0.40448e08,
        0.43967e08,
        0.47765e08,
        0.51862e08,
        0.56280e08,
        0.61040e08,
        0.66167e08,
        0.71686e08,
        0.77624e08,
        0.84009e08,
        0.90873e08,
        0.98247e08,
        0.10616e09,
        0.11466e09,
        0.12378e09,
        0.13356e09,
        0.14403e09,
        0.15526e09,
        0.16728e09,
        0.18014e09,
        0.19391e09,
        0.20863e09,
        0.22436e09,
        0.24117e09,
        0.25913e09,
        0.27830e09,
        0.29875e09,
        0.32057e09,
        0.34384e09,
        0.36864e09,
        0.39506e09,
        0.42320e09,
        0.45316e09,
        0.48504e09,
        0.51896e09,
        0.55502e09,
        0.59336e09,
        0.63410e09,
        0.67738e09,
        0.72334e09,
        0.77212e09,
        0.82388e09,
        0.87879e09,
        0.93701e09,
        0.99873e09,
        0.10641e10,
        0.11334e10,
        0.12068e10,
        0.12845e10,
        0.13667e10,
        0.14536e10,
    ]
)


#  --------------- CH3Br 211: M = 40, I = 2 ---------------------
M = 40
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.70566e04,
        0.11892e05,
        0.17508e05,
        0.23832e05,
        0.30841e05,
        0.38557e05,
        0.47036e05,
        0.56362e05,
        0.66644e05,
        0.78011e05,
        0.90615e05,
        0.10462e06,
        0.12023e06,
        0.13763e06,
        0.15707e06,
        0.17880e06,
        0.20308e06,
        0.23023e06,
        0.26059e06,
        0.29451e06,
        0.33240e06,
        0.37471e06,
        0.42191e06,
        0.47453e06,
        0.53316e06,
        0.59843e06,
        0.67104e06,
        0.75176e06,
        0.84141e06,
        0.94090e06,
        0.10512e07,
        0.11735e07,
        0.13088e07,
        0.14585e07,
        0.16241e07,
        0.18069e07,
        0.20086e07,
        0.22312e07,
        0.24764e07,
        0.27464e07,
        0.30435e07,
        0.33702e07,
        0.37291e07,
        0.41231e07,
        0.45554e07,
        0.50294e07,
        0.55486e07,
        0.61171e07,
        0.67389e07,
        0.74188e07,
        0.81616e07,
        0.89725e07,
        0.98573e07,
        0.10822e08,
        0.11873e08,
        0.13018e08,
        0.14263e08,
        0.15618e08,
        0.17090e08,
        0.18689e08,
        0.20425e08,
        0.22308e08,
        0.24350e08,
        0.26563e08,
        0.28959e08,
        0.31552e08,
        0.34357e08,
        0.37389e08,
        0.40666e08,
        0.44204e08,
        0.48023e08,
        0.52143e08,
        0.56585e08,
        0.61371e08,
        0.66526e08,
        0.72076e08,
        0.78046e08,
        0.84467e08,
        0.91369e08,
        0.98783e08,
        0.10674e09,
        0.11529e09,
        0.12446e09,
        0.13429e09,
        0.14482e09,
        0.15611e09,
        0.16820e09,
        0.18113e09,
        0.19497e09,
        0.20978e09,
        0.22560e09,
        0.24250e09,
        0.26056e09,
        0.27983e09,
        0.30040e09,
        0.32234e09,
        0.34574e09,
        0.37068e09,
        0.39725e09,
        0.42555e09,
        0.45567e09,
        0.48773e09,
        0.52184e09,
        0.55811e09,
        0.59666e09,
        0.63763e09,
        0.68115e09,
        0.72736e09,
        0.77642e09,
        0.82847e09,
        0.88368e09,
        0.94223e09,
        0.10043e10,
        0.10701e10,
        0.11397e10,
        0.12135e10,
        0.12916e10,
        0.13743e10,
        0.14618e10,
    ]
)


#  --------------- CH3CN 2124: M = 41, I = 1 ---------------------
M = 41
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.54361e04,
        0.91953e04,
        0.13708e05,
        0.19097e05,
        0.25531e05,
        0.33206e05,
        0.42337e05,
        0.53173e05,
        0.66002e05,
        0.81163e05,
        0.99053e05,
        0.12014e06,
        0.14496e06,
        0.17414e06,
        0.20843e06,
        0.24866e06,
        0.29580e06,
        0.35099e06,
        0.41551e06,
        0.49085e06,
        0.57871e06,
        0.68104e06,
        0.80008e06,
        0.93836e06,
        0.10988e07,
        0.12848e07,
        0.14999e07,
        0.17487e07,
        0.20359e07,
        0.23670e07,
        0.27484e07,
        0.31871e07,
        0.36912e07,
        0.42697e07,
        0.49328e07,
        0.56921e07,
        0.65605e07,
        0.75526e07,
        0.86847e07,
        0.99753e07,
        0.11445e08,
        0.13116e08,
        0.15016e08,
        0.17172e08,
        0.19617e08,
        0.22386e08,
        0.25520e08,
        0.29063e08,
        0.33064e08,
        0.37578e08,
        0.42667e08,
        0.48397e08,
        0.54844e08,
        0.62090e08,
        0.70228e08,
        0.79358e08,
        0.89592e08,
        0.10105e09,
        0.11388e09,
        0.12822e09,
        0.14424e09,
        0.16212e09,
        0.18205e09,
        0.20427e09,
        0.22900e09,
        0.25652e09,
        0.28710e09,
        0.32107e09,
        0.35877e09,
        0.40059e09,
        0.44692e09,
        0.49822e09,
        0.55500e09,
        0.61777e09,
        0.68712e09,
        0.76370e09,
        0.84819e09,
        0.94135e09,
        0.10440e10,
        0.11570e10,
        0.12814e10,
        0.14181e10,
        0.15684e10,
        0.17334e10,
        0.19145e10,
        0.21131e10,
        0.23308e10,
        0.25693e10,
        0.28304e10,
        0.31161e10,
        0.34285e10,
        0.37698e10,
        0.41426e10,
        0.45496e10,
        0.49935e10,
        0.54776e10,
        0.60051e10,
        0.65796e10,
        0.72049e10,
        0.78853e10,
        0.86251e10,
        0.94291e10,
        0.10303e11,
        0.11251e11,
        0.12280e11,
        0.13396e11,
        0.14606e11,
        0.15916e11,
        0.17336e11,
        0.18873e11,
        0.20536e11,
        0.22334e11,
        0.24278e11,
        0.26379e11,
        0.28647e11,
        0.31096e11,
        0.33739e11,
        0.36589e11,
        0.39661e11,
    ]
)


#  --------------- CH3CN 2134: M = 41, I = 2 --------------------- not in HITRAN-2012
M = 41
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.10906e05,
        0.18458e05,
        0.27552e05,
        0.38455e05,
        0.51523e05,
        0.67161e05,
        0.85818e05,
        0.10801e06,
        0.13434e06,
        0.16550e06,
        0.20234e06,
        0.24581e06,
        0.29705e06,
        0.35737e06,
        0.42831e06,
        0.51162e06,
        0.60936e06,
        0.72387e06,
        0.85786e06,
        0.10145e07,
        0.11972e07,
        0.14102e07,
        0.16582e07,
        0.19465e07,
        0.22813e07,
        0.26695e07,
        0.31190e07,
        0.36390e07,
        0.42397e07,
        0.49328e07,
        0.57314e07,
        0.66507e07,
        0.77076e07,
        0.89211e07,
        0.10313e08,
        0.11907e08,
        0.13732e08,
        0.15817e08,
        0.18198e08,
        0.20914e08,
        0.24007e08,
        0.27527e08,
        0.31529e08,
        0.36073e08,
        0.41228e08,
        0.47070e08,
        0.53683e08,
        0.61162e08,
        0.69612e08,
        0.79149e08,
        0.89903e08,
        0.10202e09,
        0.11565e09,
        0.13098e09,
        0.14820e09,
        0.16753e09,
        0.18921e09,
        0.21349e09,
        0.24066e09,
        0.27106e09,
        0.30502e09,
        0.34293e09,
        0.38523e09,
        0.43237e09,
        0.48486e09,
        0.54328e09,
        0.60823e09,
        0.68039e09,
        0.76049e09,
        0.84935e09,
        0.94784e09,
        0.10569e10,
        0.11777e10,
        0.13112e10,
        0.14588e10,
        0.16217e10,
        0.18016e10,
        0.19999e10,
        0.22185e10,
        0.24592e10,
        0.27241e10,
        0.30155e10,
        0.33357e10,
        0.36875e10,
        0.40736e10,
        0.44971e10,
        0.49615e10,
        0.54702e10,
        0.60273e10,
        0.66369e10,
        0.73035e10,
        0.80322e10,
        0.88282e10,
        0.96972e10,
        0.10645e11,
        0.11679e11,
        0.12806e11,
        0.14034e11,
        0.15370e11,
        0.16824e11,
        0.18406e11,
        0.20125e11,
        0.21992e11,
        0.24020e11,
        0.26221e11,
        0.28608e11,
        0.31197e11,
        0.34002e11,
        0.37040e11,
        0.40330e11,
        0.43889e11,
        0.47739e11,
        0.51902e11,
        0.56400e11,
        0.61259e11,
        0.66504e11,
        0.72165e11,
        0.78272e11,
        0.84856e11,
    ]
)


#  --------------- CH3CN 3124: M = 41, I = 3 --------------------- not in HITRAN-2012
M = 41
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.11223e05,
        0.18985e05,
        0.28307e05,
        0.39441e05,
        0.52744e05,
        0.68620e05,
        0.87523e05,
        0.10997e06,
        0.13658e06,
        0.16806e06,
        0.20524e06,
        0.24910e06,
        0.30080e06,
        0.36165e06,
        0.43319e06,
        0.51722e06,
        0.61579e06,
        0.73127e06,
        0.86640e06,
        0.10243e07,
        0.12086e07,
        0.14234e07,
        0.16735e07,
        0.19642e07,
        0.23017e07,
        0.26931e07,
        0.31464e07,
        0.36706e07,
        0.42762e07,
        0.49749e07,
        0.57801e07,
        0.67069e07,
        0.77722e07,
        0.89955e07,
        0.10398e08,
        0.12006e08,
        0.13845e08,
        0.15947e08,
        0.18346e08,
        0.21083e08,
        0.24201e08,
        0.27748e08,
        0.31781e08,
        0.36361e08,
        0.41556e08,
        0.47442e08,
        0.54106e08,
        0.61643e08,
        0.70157e08,
        0.79767e08,
        0.90604e08,
        0.10281e09,
        0.11655e09,
        0.13199e09,
        0.14935e09,
        0.16882e09,
        0.19065e09,
        0.21512e09,
        0.24250e09,
        0.27312e09,
        0.30733e09,
        0.34553e09,
        0.38814e09,
        0.43562e09,
        0.48851e09,
        0.54736e09,
        0.61279e09,
        0.68548e09,
        0.76617e09,
        0.85568e09,
        0.95489e09,
        0.10648e10,
        0.11864e10,
        0.13209e10,
        0.14695e10,
        0.16337e10,
        0.18148e10,
        0.20146e10,
        0.22348e10,
        0.24772e10,
        0.27441e10,
        0.30375e10,
        0.33601e10,
        0.37143e10,
        0.41032e10,
        0.45298e10,
        0.49975e10,
        0.55099e10,
        0.60709e10,
        0.66849e10,
        0.73563e10,
        0.80902e10,
        0.88918e10,
        0.97670e10,
        0.10722e11,
        0.11763e11,
        0.12898e11,
        0.14134e11,
        0.15480e11,
        0.16945e11,
        0.18537e11,
        0.20269e11,
        0.22149e11,
        0.24191e11,
        0.26408e11,
        0.28812e11,
        0.31419e11,
        0.34244e11,
        0.37303e11,
        0.40616e11,
        0.44201e11,
        0.48078e11,
        0.52269e11,
        0.56799e11,
        0.61692e11,
        0.66974e11,
        0.72675e11,
        0.78824e11,
        0.85454e11,
    ]
)


#  --------------- CH3CN 3134: M = 41, I = 4 --------------------- not in HITRAN-2012
M = 41
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.22522e05,
        0.38117e05,
        0.56899e05,
        0.79412e05,
        0.10640e06,
        0.13870e06,
        0.17726e06,
        0.22314e06,
        0.27761e06,
        0.34214e06,
        0.41847e06,
        0.50862e06,
        0.61497e06,
        0.74028e06,
        0.88774e06,
        0.10611e07,
        0.12646e07,
        0.15031e07,
        0.17825e07,
        0.21092e07,
        0.24908e07,
        0.29358e07,
        0.34541e07,
        0.40571e07,
        0.47576e07,
        0.55703e07,
        0.65120e07,
        0.76018e07,
        0.88614e07,
        0.10315e08,
        0.11992e08,
        0.13922e08,
        0.16142e08,
        0.18693e08,
        0.21619e08,
        0.24973e08,
        0.28812e08,
        0.33202e08,
        0.38216e08,
        0.43936e08,
        0.50455e08,
        0.57876e08,
        0.66315e08,
        0.75901e08,
        0.86779e08,
        0.99110e08,
        0.11307e09,
        0.12887e09,
        0.14672e09,
        0.16688e09,
        0.18961e09,
        0.21523e09,
        0.24407e09,
        0.27651e09,
        0.31295e09,
        0.35387e09,
        0.39975e09,
        0.45118e09,
        0.50875e09,
        0.57315e09,
        0.64512e09,
        0.72549e09,
        0.81517e09,
        0.91514e09,
        0.10265e10,
        0.11504e10,
        0.12883e10,
        0.14414e10,
        0.16115e10,
        0.18001e10,
        0.20093e10,
        0.22410e10,
        0.24975e10,
        0.27812e10,
        0.30948e10,
        0.34412e10,
        0.38235e10,
        0.42452e10,
        0.47101e10,
        0.52220e10,
        0.57856e10,
        0.64055e10,
        0.70869e10,
        0.78355e10,
        0.86574e10,
        0.95591e10,
        0.10548e11,
        0.11631e11,
        0.12817e11,
        0.14116e11,
        0.15536e11,
        0.17088e11,
        0.18785e11,
        0.20636e11,
        0.22657e11,
        0.24861e11,
        0.27264e11,
        0.29881e11,
        0.32730e11,
        0.35832e11,
        0.39205e11,
        0.42871e11,
        0.46855e11,
        0.51182e11,
        0.55878e11,
        0.60973e11,
        0.66497e11,
        0.72484e11,
        0.78970e11,
        0.85992e11,
        0.93592e11,
        0.10181e12,
        0.11070e12,
        0.12031e12,
        0.13069e12,
        0.14189e12,
        0.15398e12,
        0.16703e12,
        0.18110e12,
    ]
)


#  --------------- CF4 29: M = 42, I = 1 ---------------------
M = 42
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.76233e04,
        0.12867e05,
        0.19059e05,
        0.26316e05,
        0.34895e05,
        0.45145e05,
        0.57461e05,
        0.72259e05,
        0.89950e05,
        0.11092e06,
        0.13550e06,
        0.16399e06,
        0.19658e06,
        0.23341e06,
        0.27457e06,
        0.32004e06,
        0.36978e06,
        0.42369e06,
        0.48161e06,
        0.54338e06,
        0.60880e06,
        0.67764e06,
        0.55684e07,
        0.71250e07,
        0.90615e07,
        0.11458e08,
        0.14407e08,
        0.18021e08,
        0.22428e08,
        0.27778e08,
        0.34247e08,
        0.42038e08,
        0.51386e08,
        0.62559e08,
        0.75869e08,
        0.91670e08,
        0.11037e09,
        0.13242e09,
        0.15836e09,
        0.18878e09,
        0.22436e09,
        0.26584e09,
        0.31410e09,
        0.37008e09,
        0.43488e09,
        0.50970e09,
        0.59589e09,
        0.69496e09,
        0.80858e09,
        0.93863e09,
        0.10872e10,
        0.12565e10,
        0.14491e10,
        0.16679e10,
        0.19159e10,
        0.21966e10,
        0.25136e10,
        0.28711e10,
        0.32740e10,
        0.37260e10,
        0.42340e10,
        0.48030e10,
        0.54400e10,
        0.61520e10,
        0.69470e10,
        0.78320e10,
        0.88170e10,
        0.99120e10,
        0.11130e11,
        0.12470e11,
        0.13970e11,
        0.15620e11,
        0.17440e11,
        0.19450e11,
        0.21670e11,
        0.24100e11,
        0.26790e11,
        0.29730e11,
        0.33000e11,
        0.36500e11,
        0.40400e11,
        0.44600e11,
        0.49300e11,
        0.54300e11,
        0.59800e11,
        0.65800e11,
        0.72400e11,
        0.79500e11,
        0.87200e11,
        0.95500e11,
        0.10500e12,
        0.11400e12,
        0.12500e12,
        0.13600e12,
        0.14900e12,
        0.16200e12,
        0.17700e12,
        0.19200e12,
        0.21000e12,
        0.23000e12,
        0.25000e12,
        0.27000e12,
        0.29000e12,
        0.31000e12,
        0.34000e12,
        0.36000e12,
        0.39000e12,
        0.42000e12,
        0.46000e12,
        0.49000e12,
        0.53000e12,
        0.57000e12,
        0.61000e12,
        0.66000e12,
        0.70000e12,
        0.75000e12,
        0.81000e12,
        0.86000e12,
        0.93000e12,
    ]
)


#  --------------- C4H2 1221: M = 43, I = 1 ---------------------
M = 43
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.57628e03,
        0.84874e03,
        0.11789e04,
        0.15952e04,
        0.21317e04,
        0.28324e04,
        0.37543e04,
        0.49705e04,
        0.65754e04,
        0.86894e04,
        0.11466e05,
        0.15099e05,
        0.19834e05,
        0.25980e05,
        0.33920e05,
        0.44132e05,
        0.57210e05,
        0.73884e05,
        0.95049e05,
        0.12180e06,
        0.15548e06,
        0.19771e06,
        0.25045e06,
        0.31606e06,
        0.39739e06,
        0.49786e06,
        0.62152e06,
        0.77324e06,
        0.95878e06,
        0.11850e07,
        0.14599e07,
        0.17930e07,
        0.21956e07,
        0.26807e07,
        0.32637e07,
        0.39626e07,
        0.47983e07,
        0.57951e07,
        0.69813e07,
        0.83896e07,
        0.10058e08,
        0.12030e08,
        0.14356e08,
        0.17093e08,
        0.20309e08,
        0.24079e08,
        0.28491e08,
        0.33644e08,
        0.39651e08,
        0.46642e08,
        0.54764e08,
        0.64184e08,
        0.75091e08,
        0.87699e08,
        0.10225e09,
        0.11902e09,
        0.13832e09,
        0.16049e09,
        0.18593e09,
        0.21507e09,
        0.24841e09,
        0.28650e09,
        0.32996e09,
        0.37949e09,
        0.43586e09,
        0.49993e09,
        0.57266e09,
        0.65513e09,
        0.74852e09,
        0.85418e09,
        0.97356e09,
        0.11083e10,
        0.12602e10,
        0.14313e10,
        0.16238e10,
        0.18401e10,
        0.20829e10,
        0.23553e10,
        0.26605e10,
        0.30021e10,
        0.33841e10,
        0.38109e10,
        0.42874e10,
        0.48187e10,
        0.54107e10,
        0.60698e10,
        0.68029e10,
        0.76176e10,
        0.85223e10,
        0.95260e10,
        0.10639e11,
        0.11871e11,
        0.13236e11,
        0.14744e11,
        0.16412e11,
        0.18253e11,
        0.20285e11,
        0.22526e11,
        0.24995e11,
        0.27714e11,
        0.30705e11,
        0.33995e11,
        0.37609e11,
        0.41579e11,
        0.45934e11,
        0.50711e11,
        0.55947e11,
        0.61681e11,
        0.67957e11,
        0.74824e11,
        0.82330e11,
        0.90532e11,
        0.99487e11,
        0.10926e12,
        0.11992e12,
        0.13154e12,
        0.14420e12,
        0.15799e12,
        0.17299e12,
    ]
)


#  --------------- HC3N 12224: M = 44, I = 1 --------------------- 1224 in HITRAN, 12224 in TIPS
M = 44
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.16683e04,
        0.24538e04,
        0.33995e04,
        0.45769e04,
        0.60637e04,
        0.79533e04,
        0.10360e05,
        0.13422e05,
        0.17311e05,
        0.22232e05,
        0.28434e05,
        0.36215e05,
        0.45932e05,
        0.58011e05,
        0.72958e05,
        0.91370e05,
        0.11395e06,
        0.14153e06,
        0.17507e06,
        0.21570e06,
        0.26475e06,
        0.32372e06,
        0.39440e06,
        0.47881e06,
        0.57930e06,
        0.69856e06,
        0.83968e06,
        0.10062e07,
        0.12021e07,
        0.14320e07,
        0.17011e07,
        0.20153e07,
        0.23812e07,
        0.28065e07,
        0.32996e07,
        0.38701e07,
        0.45287e07,
        0.52876e07,
        0.61602e07,
        0.71616e07,
        0.83088e07,
        0.96206e07,
        0.11118e08,
        0.12824e08,
        0.14765e08,
        0.16969e08,
        0.19469e08,
        0.22299e08,
        0.25498e08,
        0.29110e08,
        0.33181e08,
        0.37763e08,
        0.42914e08,
        0.48697e08,
        0.55180e08,
        0.62440e08,
        0.70558e08,
        0.79627e08,
        0.89743e08,
        0.10102e09,
        0.11356e09,
        0.12752e09,
        0.14301e09,
        0.16020e09,
        0.17925e09,
        0.20035e09,
        0.22367e09,
        0.24945e09,
        0.27790e09,
        0.30928e09,
        0.34385e09,
        0.38191e09,
        0.42376e09,
        0.46975e09,
        0.52023e09,
        0.57562e09,
        0.63632e09,
        0.70279e09,
        0.77553e09,
        0.85506e09,
        0.94195e09,
        0.10368e10,
        0.11403e10,
        0.12531e10,
        0.13759e10,
        0.15097e10,
        0.16552e10,
        0.18133e10,
        0.19851e10,
        0.21715e10,
        0.23738e10,
        0.25931e10,
        0.28307e10,
        0.30879e10,
        0.33662e10,
        0.36672e10,
        0.39926e10,
        0.43439e10,
        0.47233e10,
        0.51325e10,
        0.55738e10,
        0.60493e10,
        0.65615e10,
        0.71129e10,
        0.77061e10,
        0.83441e10,
        0.90298e10,
        0.97664e10,
        0.10557e11,
        0.11406e11,
        0.12317e11,
        0.13293e11,
        0.14339e11,
        0.15459e11,
        0.16659e11,
        0.17942e11,
        0.19316e11,
        0.20784e11,
        0.22353e11,
    ]
)


#  --------------- HC3N 12234: M = 44, I = 2 --------------------- see above
M = 44
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.33507e04,
        0.49290e04,
        0.68293e04,
        0.91959e04,
        0.12185e05,
        0.15986e05,
        0.20828e05,
        0.26993e05,
        0.34824e05,
        0.44739e05,
        0.57239e05,
        0.72931e05,
        0.92539e05,
        0.11693e06,
        0.14713e06,
        0.18435e06,
        0.23004e06,
        0.28588e06,
        0.35384e06,
        0.43625e06,
        0.53580e06,
        0.65562e06,
        0.79933e06,
        0.97115e06,
        0.11759e07,
        0.14191e07,
        0.17073e07,
        0.20476e07,
        0.24486e07,
        0.29196e07,
        0.34716e07,
        0.41169e07,
        0.48696e07,
        0.57453e07,
        0.67621e07,
        0.79402e07,
        0.93022e07,
        0.10874e08,
        0.12684e08,
        0.14764e08,
        0.17150e08,
        0.19884e08,
        0.23009e08,
        0.26576e08,
        0.30641e08,
        0.35265e08,
        0.40518e08,
        0.46477e08,
        0.53225e08,
        0.60856e08,
        0.69475e08,
        0.79195e08,
        0.90143e08,
        0.10246e09,
        0.11629e09,
        0.13182e09,
        0.14921e09,
        0.16868e09,
        0.19045e09,
        0.21477e09,
        0.24189e09,
        0.27211e09,
        0.30575e09,
        0.34316e09,
        0.38471e09,
        0.43083e09,
        0.48196e09,
        0.53858e09,
        0.60125e09,
        0.67052e09,
        0.74704e09,
        0.83148e09,
        0.92459e09,
        0.10272e10,
        0.11401e10,
        0.12643e10,
        0.14007e10,
        0.15506e10,
        0.17150e10,
        0.18953e10,
        0.20928e10,
        0.23090e10,
        0.25456e10,
        0.28042e10,
        0.30867e10,
        0.33951e10,
        0.37316e10,
        0.40984e10,
        0.44981e10,
        0.49332e10,
        0.54067e10,
        0.59216e10,
        0.64812e10,
        0.70890e10,
        0.77488e10,
        0.84645e10,
        0.92405e10,
        0.10081e11,
        0.10992e11,
        0.11978e11,
        0.13044e11,
        0.14197e11,
        0.15443e11,
        0.16789e11,
        0.18243e11,
        0.19810e11,
        0.21501e11,
        0.23324e11,
        0.25288e11,
        0.27403e11,
        0.29680e11,
        0.32130e11,
        0.34764e11,
        0.37596e11,
        0.40639e11,
        0.43907e11,
        0.47416e11,
        0.51181e11,
        0.55220e11,
    ]
)


#  --------------- HC3N 12324: M = 44, I = 3 --------------------- see above
M = 44
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.33506e04,
        0.49280e04,
        0.68267e04,
        0.91901e04,
        0.12174e05,
        0.15966e05,
        0.20793e05,
        0.26936e05,
        0.34734e05,
        0.44598e05,
        0.57026e05,
        0.72612e05,
        0.92071e05,
        0.11625e06,
        0.14616e06,
        0.18298e06,
        0.22813e06,
        0.28323e06,
        0.35022e06,
        0.43133e06,
        0.52918e06,
        0.64677e06,
        0.78761e06,
        0.95571e06,
        0.11557e07,
        0.13929e07,
        0.16734e07,
        0.20041e07,
        0.23929e07,
        0.28488e07,
        0.33820e07,
        0.40040e07,
        0.47280e07,
        0.55686e07,
        0.65423e07,
        0.76678e07,
        0.89661e07,
        0.10460e08,
        0.12177e08,
        0.14145e08,
        0.16397e08,
        0.18970e08,
        0.21903e08,
        0.25242e08,
        0.29036e08,
        0.33339e08,
        0.38214e08,
        0.43726e08,
        0.49949e08,
        0.56965e08,
        0.64864e08,
        0.73743e08,
        0.83711e08,
        0.94886e08,
        0.10740e09,
        0.12139e09,
        0.13701e09,
        0.15443e09,
        0.17384e09,
        0.19543e09,
        0.21943e09,
        0.24607e09,
        0.27561e09,
        0.30832e09,
        0.34452e09,
        0.38453e09,
        0.42870e09,
        0.47742e09,
        0.53110e09,
        0.59020e09,
        0.65518e09,
        0.72659e09,
        0.80496e09,
        0.89092e09,
        0.98510e09,
        0.10882e10,
        0.12010e10,
        0.13242e10,
        0.14588e10,
        0.16056e10,
        0.17657e10,
        0.19401e10,
        0.21299e10,
        0.23363e10,
        0.25606e10,
        0.28043e10,
        0.30687e10,
        0.33553e10,
        0.36660e10,
        0.40024e10,
        0.43665e10,
        0.47601e10,
        0.51856e10,
        0.56450e10,
        0.61408e10,
        0.66756e10,
        0.72520e10,
        0.78729e10,
        0.85413e10,
        0.92604e10,
        0.10034e11,
        0.10864e11,
        0.11757e11,
        0.12714e11,
        0.13742e11,
        0.14843e11,
        0.16023e11,
        0.17287e11,
        0.18640e11,
        0.20087e11,
        0.21634e11,
        0.23288e11,
        0.25054e11,
        0.26939e11,
        0.28950e11,
        0.31096e11,
        0.33382e11,
        0.35819e11,
        0.38413e11,
    ]
)


#  --------------- HC3N 13224: M = 44, I = 4 --------------------- see above
M = 44
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.34439e04,
        0.50672e04,
        0.70230e04,
        0.94603e04,
        0.12542e05,
        0.16462e05,
        0.21461e05,
        0.27833e05,
        0.35935e05,
        0.46204e05,
        0.59168e05,
        0.75463e05,
        0.95854e05,
        0.12126e06,
        0.15276e06,
        0.19165e06,
        0.23947e06,
        0.29802e06,
        0.36943e06,
        0.45619e06,
        0.56121e06,
        0.68789e06,
        0.84018e06,
        0.10227e07,
        0.12407e07,
        0.15003e07,
        0.18086e07,
        0.21738e07,
        0.26052e07,
        0.31134e07,
        0.37106e07,
        0.44109e07,
        0.52300e07,
        0.61861e07,
        0.72996e07,
        0.85939e07,
        0.10095e08,
        0.11833e08,
        0.13841e08,
        0.16158e08,
        0.18825e08,
        0.21890e08,
        0.25407e08,
        0.29436e08,
        0.34045e08,
        0.39308e08,
        0.45309e08,
        0.52143e08,
        0.59912e08,
        0.68734e08,
        0.78737e08,
        0.90065e08,
        0.10288e09,
        0.11735e09,
        0.13367e09,
        0.15206e09,
        0.17277e09,
        0.19604e09,
        0.22217e09,
        0.25148e09,
        0.28432e09,
        0.32108e09,
        0.36218e09,
        0.40809e09,
        0.45932e09,
        0.51644e09,
        0.58004e09,
        0.65082e09,
        0.72950e09,
        0.81690e09,
        0.91388e09,
        0.10214e10,
        0.11405e10,
        0.12724e10,
        0.14182e10,
        0.15794e10,
        0.17573e10,
        0.19536e10,
        0.21701e10,
        0.24086e10,
        0.26711e10,
        0.29599e10,
        0.32774e10,
        0.36262e10,
        0.40090e10,
        0.44290e10,
        0.48895e10,
        0.53939e10,
        0.59462e10,
        0.65504e10,
        0.72111e10,
        0.79332e10,
        0.87217e10,
        0.95823e10,
        0.10521e11,
        0.11544e11,
        0.12659e11,
        0.13874e11,
        0.15195e11,
        0.16632e11,
        0.18194e11,
        0.19892e11,
        0.21735e11,
        0.23736e11,
        0.25907e11,
        0.28260e11,
        0.30810e11,
        0.33572e11,
        0.36563e11,
        0.39799e11,
        0.43299e11,
        0.47083e11,
        0.51172e11,
        0.55588e11,
        0.60355e11,
        0.65500e11,
        0.71049e11,
        0.77031e11,
        0.83478e11,
    ]
)


#  --------------- HC3N 12225: M = 44, I = 5 --------------------- see above
M = 44
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.11455e04,
        0.16850e04,
        0.23345e04,
        0.31432e04,
        0.41647e04,
        0.54630e04,
        0.71168e04,
        0.92219e04,
        0.11895e05,
        0.15279e05,
        0.19545e05,
        0.24897e05,
        0.31584e05,
        0.39899e05,
        0.50190e05,
        0.62871e05,
        0.78428e05,
        0.97434e05,
        0.12056e06,
        0.14859e06,
        0.18243e06,
        0.22314e06,
        0.27194e06,
        0.33026e06,
        0.39972e06,
        0.48219e06,
        0.57983e06,
        0.69509e06,
        0.83077e06,
        0.99009e06,
        0.11767e07,
        0.13946e07,
        0.16487e07,
        0.19441e07,
        0.22868e07,
        0.26836e07,
        0.31420e07,
        0.36704e07,
        0.42786e07,
        0.49770e07,
        0.57776e07,
        0.66938e07,
        0.77404e07,
        0.89339e07,
        0.10293e08,
        0.11837e08,
        0.13590e08,
        0.15576e08,
        0.17823e08,
        0.20362e08,
        0.23227e08,
        0.26454e08,
        0.30085e08,
        0.34166e08,
        0.38745e08,
        0.43877e08,
        0.49622e08,
        0.56046e08,
        0.63219e08,
        0.71222e08,
        0.80138e08,
        0.90062e08,
        0.10110e09,
        0.11335e09,
        0.12695e09,
        0.14202e09,
        0.15870e09,
        0.17716e09,
        0.19756e09,
        0.22009e09,
        0.24493e09,
        0.27232e09,
        0.30247e09,
        0.33565e09,
        0.37211e09,
        0.41217e09,
        0.45613e09,
        0.50433e09,
        0.55714e09,
        0.61497e09,
        0.67823e09,
        0.74739e09,
        0.82293e09,
        0.90540e09,
        0.99536e09,
        0.10934e10,
        0.12002e10,
        0.13165e10,
        0.14430e10,
        0.15805e10,
        0.17299e10,
        0.18922e10,
        0.20682e10,
        0.22591e10,
        0.24660e10,
        0.26901e10,
        0.29326e10,
        0.31951e10,
        0.34788e10,
        0.37854e10,
        0.41166e10,
        0.44741e10,
        0.48598e10,
        0.52758e10,
        0.57240e10,
        0.62069e10,
        0.67269e10,
        0.72864e10,
        0.78882e10,
        0.85352e10,
        0.92305e10,
        0.99773e10,
        0.10779e11,
        0.11639e11,
        0.12562e11,
        0.13552e11,
        0.14612e11,
        0.15748e11,
        0.16964e11,
    ]
)


#  --------------- HC3N 22224: M = 44, I = 6 --------------------- see above
M = 44
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(9.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.27029e04,
        0.39999e04,
        0.55894e04,
        0.76092e04,
        0.10219e05,
        0.13616e05,
        0.18042e05,
        0.23798e05,
        0.31255e05,
        0.40867e05,
        0.53189e05,
        0.68897e05,
        0.88807e05,
        0.11390e06,
        0.14537e06,
        0.18461e06,
        0.23330e06,
        0.29342e06,
        0.36733e06,
        0.45779e06,
        0.56802e06,
        0.70182e06,
        0.86361e06,
        0.10585e07,
        0.12925e07,
        0.15725e07,
        0.19064e07,
        0.23034e07,
        0.27739e07,
        0.33302e07,
        0.39858e07,
        0.47566e07,
        0.56604e07,
        0.67176e07,
        0.79511e07,
        0.93872e07,
        0.11055e08,
        0.12989e08,
        0.15225e08,
        0.17806e08,
        0.20779e08,
        0.24197e08,
        0.28119e08,
        0.32612e08,
        0.37749e08,
        0.43612e08,
        0.50294e08,
        0.57895e08,
        0.66528e08,
        0.76318e08,
        0.87403e08,
        0.99937e08,
        0.11409e09,
        0.13004e09,
        0.14800e09,
        0.16819e09,
        0.19086e09,
        0.21629e09,
        0.24476e09,
        0.27661e09,
        0.31219e09,
        0.35189e09,
        0.39615e09,
        0.44542e09,
        0.50021e09,
        0.56108e09,
        0.62862e09,
        0.70350e09,
        0.78641e09,
        0.87814e09,
        0.97952e09,
        0.10915e10,
        0.12149e10,
        0.13510e10,
        0.15008e10,
        0.16656e10,
        0.18468e10,
        0.20457e10,
        0.22640e10,
        0.25032e10,
        0.27653e10,
        0.30522e10,
        0.33659e10,
        0.37088e10,
        0.40832e10,
        0.44917e10,
        0.49371e10,
        0.54224e10,
        0.59508e10,
        0.65256e10,
        0.71507e10,
        0.78298e10,
        0.85671e10,
        0.93672e10,
        0.10235e11,
        0.11175e11,
        0.12193e11,
        0.13295e11,
        0.14487e11,
        0.15776e11,
        0.17168e11,
        0.18671e11,
        0.20293e11,
        0.22043e11,
        0.23929e11,
        0.25960e11,
        0.28148e11,
        0.30502e11,
        0.33034e11,
        0.35756e11,
        0.38681e11,
        0.41823e11,
        0.45195e11,
        0.48812e11,
        0.52692e11,
        0.56850e11,
        0.61306e11,
        0.66076e11,
        0.71183e11,
    ]
)


#  --------------- H2 11: M = 45, I = 1 ---------------------
M = 45
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.15265e01,
        0.22243e01,
        0.29619e01,
        0.36724e01,
        0.43456e01,
        0.49880e01,
        0.56090e01,
        0.62165e01,
        0.68161e01,
        0.74113e01,
        0.80044e01,
        0.85966e01,
        0.91887e01,
        0.97810e01,
        0.10374e02,
        0.10967e02,
        0.11561e02,
        0.12156e02,
        0.12751e02,
        0.13347e02,
        0.13944e02,
        0.14541e02,
        0.15139e02,
        0.15738e02,
        0.16337e02,
        0.16937e02,
        0.17538e02,
        0.18140e02,
        0.18743e02,
        0.19346e02,
        0.19951e02,
        0.20556e02,
        0.21163e02,
        0.21771e02,
        0.22379e02,
        0.22990e02,
        0.23601e02,
        0.24214e02,
        0.24829e02,
        0.25445e02,
        0.26063e02,
        0.26683e02,
        0.27304e02,
        0.27928e02,
        0.28553e02,
        0.29181e02,
        0.29811e02,
        0.30443e02,
        0.31078e02,
        0.31715e02,
        0.32355e02,
        0.32997e02,
        0.33643e02,
        0.34291e02,
        0.34942e02,
        0.35596e02,
        0.36253e02,
        0.36914e02,
        0.37578e02,
        0.38245e02,
        0.38916e02,
        0.39590e02,
        0.40268e02,
        0.40949e02,
        0.41635e02,
        0.42324e02,
        0.43017e02,
        0.43715e02,
        0.44416e02,
        0.45122e02,
        0.45831e02,
        0.46546e02,
        0.47264e02,
        0.47987e02,
        0.48714e02,
        0.49446e02,
        0.50183e02,
        0.50925e02,
        0.51671e02,
        0.52422e02,
        0.53178e02,
        0.53939e02,
        0.54705e02,
        0.55476e02,
        0.56252e02,
        0.57033e02,
        0.57820e02,
        0.58612e02,
        0.59409e02,
        0.60212e02,
        0.61020e02,
        0.61833e02,
        0.62652e02,
        0.63477e02,
        0.64308e02,
        0.65144e02,
        0.65986e02,
        0.66833e02,
        0.67687e02,
        0.68546e02,
        0.69411e02,
        0.70283e02,
        0.71160e02,
        0.72043e02,
        0.72933e02,
        0.73829e02,
        0.74730e02,
        0.75638e02,
        0.76553e02,
        0.77473e02,
        0.78400e02,
        0.79333e02,
        0.80273e02,
        0.81219e02,
        0.82172e02,
        0.83131e02,
        0.84097e02,
        0.85069e02,
        0.86048e02,
    ]
)


#  --------------- H2 12: M = 45, I = 2 ---------------------
M = 45
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.81692e01,
        0.10308e02,
        0.12557e02,
        0.14848e02,
        0.17159e02,
        0.19482e02,
        0.21815e02,
        0.24153e02,
        0.26497e02,
        0.28845e02,
        0.31197e02,
        0.33552e02,
        0.35910e02,
        0.38272e02,
        0.40636e02,
        0.43002e02,
        0.45372e02,
        0.47744e02,
        0.50119e02,
        0.52496e02,
        0.54877e02,
        0.57261e02,
        0.59649e02,
        0.62040e02,
        0.64435e02,
        0.66835e02,
        0.69240e02,
        0.71650e02,
        0.74066e02,
        0.76489e02,
        0.78918e02,
        0.81354e02,
        0.83799e02,
        0.86252e02,
        0.88715e02,
        0.91187e02,
        0.93669e02,
        0.96163e02,
        0.98668e02,
        0.10118e03,
        0.10371e03,
        0.10626e03,
        0.10881e03,
        0.11138e03,
        0.11397e03,
        0.11657e03,
        0.11919e03,
        0.12182e03,
        0.12447e03,
        0.12714e03,
        0.12982e03,
        0.13252e03,
        0.13524e03,
        0.13798e03,
        0.14074e03,
        0.14352e03,
        0.14632e03,
        0.14914e03,
        0.15198e03,
        0.15484e03,
        0.15772e03,
        0.16062e03,
        0.16355e03,
        0.16649e03,
        0.16946e03,
        0.17246e03,
        0.17547e03,
        0.17851e03,
        0.18157e03,
        0.18466e03,
        0.18777e03,
        0.19090e03,
        0.19406e03,
        0.19725e03,
        0.20045e03,
        0.20369e03,
        0.20695e03,
        0.21023e03,
        0.21354e03,
        0.21687e03,
        0.22024e03,
        0.22362e03,
        0.22704e03,
        0.23048e03,
        0.23394e03,
        0.23744e03,
        0.24096e03,
        0.24451e03,
        0.24808e03,
        0.25169e03,
        0.25532e03,
        0.25897e03,
        0.26266e03,
        0.26638e03,
        0.27012e03,
        0.27389e03,
        0.27769e03,
        0.28152e03,
        0.28537e03,
        0.28926e03,
        0.29317e03,
        0.29712e03,
        0.30109e03,
        0.30509e03,
        0.30913e03,
        0.31319e03,
        0.31728e03,
        0.32140e03,
        0.32555e03,
        0.32974e03,
        0.33395e03,
        0.33819e03,
        0.34246e03,
        0.34677e03,
        0.35110e03,
        0.35547e03,
        0.35987e03,
        0.36429e03,
        0.36875e03,
    ]
)


#  --------------- CS 22: M = 46, I = 1 ---------------------
M = 46
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.51416e02,
        0.72723e02,
        0.94044e02,
        0.11538e03,
        0.13673e03,
        0.15810e03,
        0.17949e03,
        0.20093e03,
        0.22245e03,
        0.24407e03,
        0.26582e03,
        0.28776e03,
        0.30992e03,
        0.33233e03,
        0.35504e03,
        0.37807e03,
        0.40147e03,
        0.42525e03,
        0.44944e03,
        0.47406e03,
        0.49914e03,
        0.52468e03,
        0.55071e03,
        0.57723e03,
        0.60427e03,
        0.63183e03,
        0.65991e03,
        0.68854e03,
        0.71771e03,
        0.74743e03,
        0.77771e03,
        0.80855e03,
        0.83996e03,
        0.87193e03,
        0.90449e03,
        0.93762e03,
        0.97134e03,
        0.10056e04,
        0.10405e04,
        0.10760e04,
        0.11121e04,
        0.11487e04,
        0.11860e04,
        0.12239e04,
        0.12623e04,
        0.13014e04,
        0.13410e04,
        0.13813e04,
        0.14222e04,
        0.14637e04,
        0.15057e04,
        0.15484e04,
        0.15917e04,
        0.16357e04,
        0.16802e04,
        0.17253e04,
        0.17711e04,
        0.18175e04,
        0.18645e04,
        0.19121e04,
        0.19603e04,
        0.20091e04,
        0.20586e04,
        0.21087e04,
        0.21594e04,
        0.22107e04,
        0.22626e04,
        0.23152e04,
        0.23684e04,
        0.24222e04,
        0.24767e04,
        0.25317e04,
        0.25874e04,
        0.26438e04,
        0.27007e04,
        0.27583e04,
        0.28165e04,
        0.28754e04,
        0.29348e04,
        0.29949e04,
        0.30557e04,
        0.31170e04,
        0.31790e04,
        0.32417e04,
        0.33049e04,
        0.33688e04,
        0.34334e04,
        0.34986e04,
        0.35644e04,
        0.36308e04,
        0.36979e04,
        0.37656e04,
        0.38340e04,
        0.39030e04,
        0.39727e04,
        0.40430e04,
        0.41139e04,
        0.41855e04,
        0.42577e04,
        0.43306e04,
        0.44041e04,
        0.44782e04,
        0.45530e04,
        0.46284e04,
        0.47045e04,
        0.47813e04,
        0.48587e04,
        0.49367e04,
        0.50154e04,
        0.50947e04,
        0.51747e04,
        0.52553e04,
        0.53366e04,
        0.54185e04,
        0.55011e04,
        0.55844e04,
        0.56683e04,
        0.57528e04,
        0.58380e04,
    ]
)


#  --------------- CS 24: M = 46, I = 2 ---------------------
M = 46
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.52247e02,
        0.73900e02,
        0.95568e02,
        0.11725e03,
        0.13895e03,
        0.16066e03,
        0.18241e03,
        0.20420e03,
        0.22607e03,
        0.24805e03,
        0.27018e03,
        0.29249e03,
        0.31503e03,
        0.33784e03,
        0.36096e03,
        0.38442e03,
        0.40824e03,
        0.43247e03,
        0.45712e03,
        0.48221e03,
        0.50778e03,
        0.53382e03,
        0.56037e03,
        0.58743e03,
        0.61501e03,
        0.64312e03,
        0.67179e03,
        0.70100e03,
        0.73077e03,
        0.76111e03,
        0.79202e03,
        0.82351e03,
        0.85559e03,
        0.88824e03,
        0.92149e03,
        0.95533e03,
        0.98977e03,
        0.10248e04,
        0.10605e04,
        0.10967e04,
        0.11336e04,
        0.11710e04,
        0.12091e04,
        0.12478e04,
        0.12871e04,
        0.13270e04,
        0.13675e04,
        0.14087e04,
        0.14505e04,
        0.14929e04,
        0.15359e04,
        0.15795e04,
        0.16238e04,
        0.16687e04,
        0.17142e04,
        0.17604e04,
        0.18071e04,
        0.18546e04,
        0.19026e04,
        0.19513e04,
        0.20006e04,
        0.20505e04,
        0.21011e04,
        0.21523e04,
        0.22042e04,
        0.22566e04,
        0.23098e04,
        0.23635e04,
        0.24179e04,
        0.24730e04,
        0.25286e04,
        0.25850e04,
        0.26419e04,
        0.26995e04,
        0.27578e04,
        0.28167e04,
        0.28762e04,
        0.29364e04,
        0.29972e04,
        0.30587e04,
        0.31208e04,
        0.31836e04,
        0.32470e04,
        0.33111e04,
        0.33758e04,
        0.34412e04,
        0.35072e04,
        0.35739e04,
        0.36412e04,
        0.37092e04,
        0.37778e04,
        0.38471e04,
        0.39171e04,
        0.39877e04,
        0.40589e04,
        0.41309e04,
        0.42034e04,
        0.42767e04,
        0.43505e04,
        0.44251e04,
        0.45003e04,
        0.45762e04,
        0.46527e04,
        0.47299e04,
        0.48077e04,
        0.48863e04,
        0.49654e04,
        0.50453e04,
        0.51258e04,
        0.52070e04,
        0.52888e04,
        0.53713e04,
        0.54545e04,
        0.55383e04,
        0.56229e04,
        0.57080e04,
        0.57939e04,
        0.58804e04,
        0.59676e04,
    ]
)


#  --------------- CS 32: M = 46, I = 3 ---------------------
M = 46
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.10889e03,
        0.15403e03,
        0.19920e03,
        0.24440e03,
        0.28964e03,
        0.33491e03,
        0.38026e03,
        0.42571e03,
        0.47134e03,
        0.51722e03,
        0.56342e03,
        0.61005e03,
        0.65719e03,
        0.70493e03,
        0.75334e03,
        0.80249e03,
        0.85245e03,
        0.90329e03,
        0.95504e03,
        0.10078e04,
        0.10615e04,
        0.11163e04,
        0.11721e04,
        0.12291e04,
        0.12872e04,
        0.13464e04,
        0.14068e04,
        0.14684e04,
        0.15311e04,
        0.15951e04,
        0.16604e04,
        0.17268e04,
        0.17945e04,
        0.18635e04,
        0.19337e04,
        0.20051e04,
        0.20779e04,
        0.21519e04,
        0.22272e04,
        0.23038e04,
        0.23817e04,
        0.24609e04,
        0.25414e04,
        0.26232e04,
        0.27064e04,
        0.27908e04,
        0.28765e04,
        0.29636e04,
        0.30520e04,
        0.31417e04,
        0.32327e04,
        0.33251e04,
        0.34188e04,
        0.35138e04,
        0.36102e04,
        0.37079e04,
        0.38070e04,
        0.39074e04,
        0.40091e04,
        0.41122e04,
        0.42166e04,
        0.43224e04,
        0.44295e04,
        0.45380e04,
        0.46478e04,
        0.47590e04,
        0.48715e04,
        0.49854e04,
        0.51007e04,
        0.52173e04,
        0.53353e04,
        0.54547e04,
        0.55754e04,
        0.56975e04,
        0.58210e04,
        0.59458e04,
        0.60720e04,
        0.61996e04,
        0.63285e04,
        0.64589e04,
        0.65906e04,
        0.67236e04,
        0.68581e04,
        0.69940e04,
        0.71312e04,
        0.72698e04,
        0.74098e04,
        0.75512e04,
        0.76940e04,
        0.78381e04,
        0.79837e04,
        0.81307e04,
        0.82790e04,
        0.84287e04,
        0.85799e04,
        0.87324e04,
        0.88864e04,
        0.90417e04,
        0.91984e04,
        0.93566e04,
        0.95161e04,
        0.96771e04,
        0.98394e04,
        0.10003e05,
        0.10168e05,
        0.10335e05,
        0.10503e05,
        0.10672e05,
        0.10843e05,
        0.11015e05,
        0.11189e05,
        0.11364e05,
        0.11541e05,
        0.11719e05,
        0.11898e05,
        0.12079e05,
        0.12261e05,
        0.12444e05,
        0.12630e05,
    ]
)


#  --------------- CS 23: M = 46, I = 4 ---------------------
M = 46
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.0)
TIPS_ISO_HASH[(M, I)] = float32(
    [
        0.20737e03,
        0.29330e03,
        0.37930e03,
        0.46535e03,
        0.55145e03,
        0.63764e03,
        0.72394e03,
        0.81043e03,
        0.89722e03,
        0.98443e03,
        0.10722e04,
        0.11607e04,
        0.12501e04,
        0.13406e04,
        0.14323e04,
        0.15253e04,
        0.16197e04,
        0.17158e04,
        0.18135e04,
        0.19129e04,
        0.20142e04,
        0.21174e04,
        0.22226e04,
        0.23298e04,
        0.24391e04,
        0.25504e04,
        0.26639e04,
        0.27796e04,
        0.28976e04,
        0.30177e04,
        0.31401e04,
        0.32648e04,
        0.33918e04,
        0.35211e04,
        0.36527e04,
        0.37867e04,
        0.39231e04,
        0.40618e04,
        0.42029e04,
        0.43463e04,
        0.44922e04,
        0.46405e04,
        0.47912e04,
        0.49443e04,
        0.50999e04,
        0.52579e04,
        0.54183e04,
        0.55812e04,
        0.57465e04,
        0.59143e04,
        0.60846e04,
        0.62573e04,
        0.64325e04,
        0.66102e04,
        0.67903e04,
        0.69729e04,
        0.71581e04,
        0.73457e04,
        0.75358e04,
        0.77284e04,
        0.79235e04,
        0.81211e04,
        0.83212e04,
        0.85239e04,
        0.87290e04,
        0.89367e04,
        0.91469e04,
        0.93596e04,
        0.95748e04,
        0.97926e04,
        0.10013e05,
        0.10236e05,
        0.10461e05,
        0.10689e05,
        0.10920e05,
        0.11153e05,
        0.11388e05,
        0.11626e05,
        0.11867e05,
        0.12110e05,
        0.12356e05,
        0.12604e05,
        0.12855e05,
        0.13109e05,
        0.13365e05,
        0.13623e05,
        0.13884e05,
        0.14148e05,
        0.14415e05,
        0.14683e05,
        0.14955e05,
        0.15229e05,
        0.15506e05,
        0.15785e05,
        0.16067e05,
        0.16351e05,
        0.16638e05,
        0.16928e05,
        0.17220e05,
        0.17515e05,
        0.17813e05,
        0.18113e05,
        0.18416e05,
        0.18721e05,
        0.19029e05,
        0.19340e05,
        0.19653e05,
        0.19969e05,
        0.20287e05,
        0.20608e05,
        0.20932e05,
        0.21258e05,
        0.21587e05,
        0.21919e05,
        0.22253e05,
        0.22590e05,
        0.22930e05,
        0.23272e05,
        0.23617e05,
    ]
)


#  --------------- SO3 26: M = 46, I = 1 --------------------- not in TIPS-2011
M = 47
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.0)
TIPS_ISO_HASH[(M, I)] = float32([0.0])


# NOT IN HITRAN, BUT PRESENT IN TIPS-2011
#                                                        ... extracted from iso_comparison
#
# id       M    I     COMMENT          TIPS_M    TIPS_I        iso_name                 abundance          mass       mol_name
# 101    1001    1    not in HITRAN      45                     H                        \N                 \N              H
#
# 102    1002    1    not in HITRAN      45                     He                       \N                 \N              He
#
# 104    1018    1    not in HITRAN      45                     Ar                       \N                 \N              Ar
#
#                    not in HITRAN      45         4224                                                                    C2N2
#                    not in HITRAN      45         5225                                                                    C2N2
#
#                    not in HITRAN      48         26                                                                      SO
#                    not in HITRAN      48         46                                                                      SO
#                    not in HITRAN      48         28                                                                      SO
#
#                    not in HITRAN      49         1221                                                                    C3H4
#
#                    not in HITRAN      50         2111                                                                    CH3
#
#                    not in HITRAN      51         222                                                                     CS2
#                    not in HITRAN      51         224                                                                     CS2
#                    not in HITRAN      51         223                                                                     CS2
#                    not in HITRAN      51         232                                                                     CS2


#  --------------- TIPS IMPLEMENTATION ----------------------


def BD_TIPS_2011_PYTHON(M, I, T):

    # out of temperature range
    if T < 70.0 or T > 3000.0:
        # Qt = -1.
        # gi = 0.
        # return gi,Qt
        raise Exception("TIPS: T must be between 70K and 3000K.")

    try:
        # get statistical weight for specified isotopologue
        gi = TIPS_GSI_HASH[(M, I)]
        # interpolate partition sum for specified isotopologue
        Qt = AtoB(T, Tdat, TIPS_ISO_HASH[(M, I)], TIPS_NPT)
    except KeyError:
        raise Exception("TIPS: no data for M,I = %d,%d." % (M, I))

    return gi, Qt


# Total internal partition sum
# M - molecule number
# I - isotopologue number
# T - temperature (K)
# returns (StatWeight,PartitionSum)


def partitionSum(M, I, T, step=None):
    """
    INPUT PARAMETERS: 
        M: HITRAN molecule number              (required)
        I: HITRAN isotopologue number          (required)
        T: temperature conditions              (required)
        step:  step to calculate temperatures  (optional)
    OUTPUT PARAMETERS:
        TT: list of temperatures (present only if T is a list)
        PartSum: partition sums calculated on a list of temperatures
    ---
    DESCRIPTION:
        Calculate range of partition sums at different temperatures.
        This function uses a python implementation of TIPS-2011 code:

        Reference:
            A. L. Laraia, R. R. Gamache, J. Lamouroux, I. E. Gordon, L. S. Rothman.
            Total internal partition sums to support planetary remote sensing.
            Icarus, Volume 215, Issue 1, September 2011, Pages 391–400
            http://dx.doi.org/10.1016/j.icarus.2011.06.004

        Output depends on a structure of input parameter T so that:
            1) If T is a scalar/list and step IS NOT provided,
                then calculate partition sums over each value of T.
            2) If T is a list and step parameter IS provided,
                then calculate partition sums between T[0] and T[1]
                with a given step.
    ---
    EXAMPLE OF USAGE:
        PartSum = partitionSum(1,1,[296,1000])
        TT,PartSum = partitionSum(1,1,[296,1000],step=0.1)
    ---
    """
    # partitionSum
    if not step:
        if type(T) not in set([list, tuple]):
            return BD_TIPS_2011_PYTHON(M, I, T)[1]
        else:
            return [BD_TIPS_2011_PYTHON(M, I, temp)[1] for temp in T]
    else:
        # n = (T[1]-T[0])/step
        # TT = linspace(T[0],T[1],n)
        TT = arange(T[0], T[1], step)
        return TT, array([BD_TIPS_2011_PYTHON(M, I, temp)[1] for temp in TT])


# ------------------ partition sum --------------------------------------


# ------------------ LINESHAPES -----------------------------------------

# ------------------ complex probability function -----------------------
# define static data
zone = __ComplexType__(1.0e0 + 0.0e0j)
zi = __ComplexType__(0.0e0 + 1.0e0j)
tt = __FloatType__(
    [
        0.5e0,
        1.5e0,
        2.5e0,
        3.5e0,
        4.5e0,
        5.5e0,
        6.5e0,
        7.5e0,
        8.5e0,
        9.5e0,
        10.5e0,
        11.5e0,
        12.5e0,
        13.5e0,
        14.5e0,
    ]
)
pipwoeronehalf = __FloatType__(0.564189583547756e0)

# "naive" implementation for benchmarks


def cpf3(X, Y):

    # X,Y,WR,WI - numpy arrays
    if type(X) != ndarray:
        if type(X) not in set([list, tuple]):
            X = array([X])
        else:
            X = array(X)
    if type(Y) != ndarray:
        if type(Y) not in set([list, tuple]):
            Y = array([Y])
        else:
            Y = array(Y)

    zm1 = zone / __ComplexType__(X + zi * Y)  # maybe redundant
    zm2 = zm1 ** 2
    zsum = zone
    zterm = zone

    for tt_i in tt:
        zterm *= zm2 * tt_i
        zsum += zterm

    zsum *= zi * zm1 * pipwoeronehalf

    return zsum.real, zsum.imag


T = __FloatType__(
    [
        0.314240376e0,
        0.947788391e0,
        1.59768264e0,
        2.27950708e0,
        3.02063703e0,
        3.8897249e0,
    ]
)
U = __FloatType__(
    [
        1.01172805e0,
        -0.75197147e0,
        1.2557727e-2,
        1.00220082e-2,
        -2.42068135e-4,
        5.00848061e-7,
    ]
)
S = __FloatType__(
    [
        1.393237e0,
        0.231152406e0,
        -0.155351466e0,
        6.21836624e-3,
        9.19082986e-5,
        -6.27525958e-7,
    ]
)

# Complex probability function implementation (Humlicek)


def cpf(X, Y):

    # X,Y,WR,WI - numpy arrays
    if type(X) != ndarray:
        if type(X) not in set([list, tuple]):
            X = array([X])
        else:
            X = array(X)
    if type(Y) != ndarray:
        if type(Y) not in set([list, tuple]):
            Y = array([Y])
        else:
            Y = array(Y)

    # REGION3
    index_REGION3 = where(sqrt(X ** 2 + Y ** 2) > __FloatType__(8.0e0))
    X_REGION3 = X[index_REGION3]
    Y_REGION3 = Y[index_REGION3]
    zm1 = zone / __ComplexType__(X_REGION3 + zi * Y_REGION3)
    zm2 = zm1 ** 2
    zsum_REGION3 = zone
    zterm = zone
    for tt_i in tt:
        zterm *= zm2 * tt_i
        zsum_REGION3 += zterm
    zsum_REGION3 *= zi * zm1 * pipwoeronehalf

    index_REGION12 = setdiff1d(array(arange(len(X))), array(index_REGION3))
    X_REGION12 = X[index_REGION12]
    Y_REGION12 = Y[index_REGION12]

    WR = __FloatType__(0.0e0)
    WI = __FloatType__(0.0e0)

    # REGION12
    Y1_REGION12 = Y_REGION12 + __FloatType__(1.5e0)
    Y2_REGION12 = Y1_REGION12 ** 2

    # REGION2
    subindex_REGION2 = where(
        (Y_REGION12 <= 0.85e0) & (abs(X_REGION12) >= (18.1e0 * Y_REGION12 + 1.65e0))
    )

    index_REGION2 = index_REGION12[subindex_REGION2]

    X_REGION2 = X[index_REGION2]
    Y_REGION2 = Y[index_REGION2]
    Y1_REGION2 = Y1_REGION12[subindex_REGION2]
    Y2_REGION2 = Y2_REGION12[subindex_REGION2]
    Y3_REGION2 = Y_REGION2 + __FloatType__(3.0e0)

    WR_REGION2 = WR
    WI_REGION2 = WI

    WR_REGION2 = zeros(len(X_REGION2))
    ii = abs(X_REGION2) < __FloatType__(12.0e0)
    WR_REGION2[ii] = exp(-X_REGION2[ii] ** 2)
    WR_REGION2[~ii] = WR

    for I in range(6):
        R_REGION2 = X_REGION2 - T[I]
        R2_REGION2 = R_REGION2 ** 2
        D_REGION2 = __FloatType__(1.0e0) / (R2_REGION2 + Y2_REGION2)
        D1_REGION2 = Y1_REGION2 * D_REGION2
        D2_REGION2 = R_REGION2 * D_REGION2
        WR_REGION2 = WR_REGION2 + Y_REGION2 * (
            U[I] * (R_REGION2 * D2_REGION2 - 1.5e0 * D1_REGION2)
            + S[I] * Y3_REGION2 * D2_REGION2
        ) / (R2_REGION2 + 2.25e0)
        R_REGION2 = X_REGION2 + T[I]
        R2_REGION2 = R_REGION2 ** 2
        D_REGION2 = __FloatType__(1.0e0) / (R2_REGION2 + Y2_REGION2)
        D3_REGION2 = Y1_REGION2 * D_REGION2
        D4_REGION2 = R_REGION2 * D_REGION2
        WR_REGION2 = WR_REGION2 + Y_REGION2 * (
            U[I] * (R_REGION2 * D4_REGION2 - 1.5e0 * D3_REGION2)
            - S[I] * Y3_REGION2 * D4_REGION2
        ) / (R2_REGION2 + 2.25e0)
        WI_REGION2 = (
            WI_REGION2
            + U[I] * (D2_REGION2 + D4_REGION2)
            + S[I] * (D1_REGION2 - D3_REGION2)
        )

    # REGION3
    index_REGION1 = setdiff1d(array(index_REGION12), array(index_REGION2))
    X_REGION1 = X[index_REGION1]
    Y_REGION1 = X[index_REGION1]

    subindex_REGION1 = setdiff1d(
        array(arange(len(index_REGION12))), array(subindex_REGION2)
    )
    Y1_REGION1 = Y1_REGION12[subindex_REGION1]
    Y2_REGION1 = Y2_REGION12[subindex_REGION1]

    WR_REGION1 = WR
    WI_REGION1 = WI

    for I in range(6):
        R_REGION1 = X_REGION1 - T[I]
        D_REGION1 = __FloatType__(1.0e0) / (R_REGION1 ** 2 + Y2_REGION1)
        D1_REGION1 = Y1_REGION1 * D_REGION1
        D2_REGION1 = R_REGION1 * D_REGION1
        R_REGION1 = X_REGION1 + T[I]
        D_REGION1 = __FloatType__(1.0e0) / (R_REGION1 ** 2 + Y2_REGION1)
        D3_REGION1 = Y1_REGION1 * D_REGION1
        D4_REGION1 = R_REGION1 * D_REGION1

        WR_REGION1 = (
            WR_REGION1
            + U[I] * (D1_REGION1 + D3_REGION1)
            - S[I] * (D2_REGION1 - D4_REGION1)
        )
        WI_REGION1 = (
            WI_REGION1
            + U[I] * (D2_REGION1 + D4_REGION1)
            + S[I] * (D1_REGION1 - D3_REGION1)
        )

    # total result
    WR_TOTAL = zeros(len(X))
    WI_TOTAL = zeros(len(X))
    # REGION3
    WR_TOTAL[index_REGION3] = zsum_REGION3.real
    WI_TOTAL[index_REGION3] = zsum_REGION3.imag
    # REGION2
    WR_TOTAL[index_REGION2] = WR_REGION2
    WI_TOTAL[index_REGION2] = WI_REGION2
    # REGION1
    WR_TOTAL[index_REGION1] = WR_REGION1
    WI_TOTAL[index_REGION1] = WI_REGION1

    return WR_TOTAL, WI_TOTAL


hcpf = cpf  # stub for initial cpf

# ------------------ Schreier CPF ------------------------

# "Optimized implementations of rational approximations
#  for the Voigt and complex error function".
# Franz Schreier. JQSRT 112 (2011) 1010-10250
# doi:10.1016/j.jqsrt.2010.12.010

# Enable this if numpy.polyval doesn't perform well.
"""    
def polyval(p, x):
    y = zeros(x.shape, dtype=float)
    for i, v in enumerate(p):
        y *= x
        y += v
    return y
"""


def cef(x, y, N):
    # Computes the function w(z) = exp(-zA2) erfc(-iz) using a rational
    # series with N terms. It is assumed that Im(z) > 0 or Im(z) = 0.
    z = x + 1.0j * y
    M = 2 * N
    M2 = 2 * M
    k = arange(-M + 1, M)  # '; # M2 = no. of sampling points.
    L = sqrt(N / sqrt(2))  # Optimal choice of L.
    theta = k * pi / M
    t = L * tan(theta / 2)  # Variables theta and t.
    # f = exp(-t.A2)*(LA2+t.A2); f = [0; f]; # Function to be transformed.
    f = zeros(len(t) + 1)
    f[0] = 0
    f[1:] = exp(-(t ** 2)) * (L ** 2 + t ** 2)
    # f = insert(exp(-t**2)*(L**2+t**2),0,0)
    a = real(fft(fftshift(f))) / M2  # Coefficients of transform.
    a = flipud(a[1 : N + 1])  # Reorder coefficients.
    Z = (L + 1.0j * z) / (L - 1.0j * z)
    p = polyval(a, Z)  # Polynomial evaluation.
    w = 2 * p / (L - 1.0j * z) ** 2 + (1 / sqrt(pi)) / (L - 1.0j * z)  # Evaluate w(z).
    return w


# weideman24 by default
# weideman24 = lambda x,y: cef(x,y,24)
def weideman(x, y, n):
    return cef(x, y, n)


def hum1_wei(x, y, n=24):
    t = y - 1.0j * x
    cerf = 1 / sqrt(pi) * t / (0.5 + t ** 2)
    """
    z = x+1j*y
    cerf = 1j*z/sqrt(pi)/(z**2-0.5)
    """
    mask = abs(x) + y < 15.0
    if any(mask):
        w24 = weideman(x[mask], y[mask], n)
        place(cerf, mask, w24)
    return cerf.real, cerf.imag


VARIABLES["CPF"] = hum1_wei
# VARIABLES['CPF'] = cpf

# ------------------ Hartmann-Tran Profile (HTP) ------------------------


def pcqsdhc(sg0, GamD, Gam0, Gam2, Shift0, Shift2, anuVC, eta, sg):
    # -------------------------------------------------
    #      "pCqSDHC": partially-Correlated quadratic-Speed-Dependent Hard-Collision
    #      Subroutine to Compute the complex normalized spectral shape of an
    #      isolated line by the pCqSDHC model
    #
    #      Reference:
    #      H. Tran, N.H. Ngo, J.-M. Hartmann.
    #      Efficient computation of some speed-dependent isolated line profiles.
    #      JQSRT, Volume 129, November 2013, Pages 199–203
    #      http://dx.doi.org/10.1016/j.jqsrt.2013.06.015
    #
    #      Input/Output Parameters of Routine (Arguments or Common)
    #      ---------------------------------
    #      T          : Temperature in Kelvin (Input).
    #      amM1       : Molar mass of the absorber in g/mol(Input).
    #      sg0        : Unperturbed line position in cm-1 (Input).
    #      GamD       : Doppler HWHM in cm-1 (Input)
    #      Gam0       : Speed-averaged line-width in cm-1 (Input).
    #      Gam2       : Speed dependence of the line-width in cm-1 (Input).
    #      anuVC      : Velocity-changing frequency in cm-1 (Input).
    #      eta        : Correlation parameter, No unit (Input).
    #      Shift0     : Speed-averaged line-shift in cm-1 (Input).
    #      Shift2     : Speed dependence of the line-shift in cm-1 (Input)
    #      sg         : Current WaveNumber of the Computation in cm-1 (Input).
    #
    #      Output Quantities (through Common Statements)
    #      -----------------
    #      LS_pCqSDHC_R: Real part of the normalized spectral shape (cm)
    #      LS_pCqSDHC_I: Imaginary part of the normalized spectral shape (cm)
    #
    #      Called Routines: 'CPF'      (Complex Probability Function)
    #      ---------------  'CPF3'      (Complex Probability Function for the region 3)
    #
    #      Called By: Main Program
    #      ---------
    #
    #     Double Precision Version
    #
    # -------------------------------------------------

    # sg is the only vector argument which is passed to fusnction

    if type(sg) not in set([array, ndarray, list, tuple]):
        sg = array([sg])

    number_of_points = len(sg)
    Aterm_GLOBAL = zeros(number_of_points, dtype=__ComplexType__)
    Bterm_GLOBAL = zeros(number_of_points, dtype=__ComplexType__)

    cte = sqrt(log(2.0e0)) / GamD
    rpi = sqrt(pi)
    iz = __ComplexType__(0.0e0 + 1.0e0j)

    c0 = __ComplexType__(Gam0 + 1.0e0j * Shift0)
    c2 = __ComplexType__(Gam2 + 1.0e0j * Shift2)
    c0t = __ComplexType__((1.0e0 - eta) * (c0 - 1.5e0 * c2) + anuVC)
    c2t = __ComplexType__((1.0e0 - eta) * c2)

    # PART1
    if abs(c2t) == 0.0e0:
        # print('PART1') # DEBUG
        Z1 = (iz * (sg0 - sg) + c0t) * cte
        xZ1 = -Z1.imag
        yZ1 = Z1.real
        # WR1,WI1 = cpf(xZ1,yZ1)
        WR1, WI1 = VARIABLES["CPF"](xZ1, yZ1)
        Aterm_GLOBAL = rpi * cte * __ComplexType__(WR1 + 1.0e0j * WI1)
        index_Z1 = abs(Z1) <= 4.0e3
        index_NOT_Z1 = ~index_Z1
        if any(index_Z1):
            # print('PART1/Z1') # DEBUG
            Bterm_GLOBAL = (
                rpi
                * cte
                * ((1.0e0 - Z1 ** 2) * __ComplexType__(WR1 + 1.0e0j * WI1) + Z1 / rpi)
            )
        if any(index_NOT_Z1):
            # print('PART1/~Z1') # DEBUG
            Bterm_GLOBAL = cte * (
                rpi * __ComplexType__(WR1 + 1.0e0j * WI1)
                + 0.5e0 / Z1
                - 0.75e0 / (Z1 ** 3)
            )
    else:
        # PART2, PART3 AND PART4   (PART4 IS A MAIN PART)

        # X - vector, Y - scalar
        X = (iz * (sg0 - sg) + c0t) / c2t
        Y = __ComplexType__(1.0e0 / ((2.0e0 * cte * c2t)) ** 2)
        csqrtY = (Gam2 - iz * Shift2) / (
            2.0e0 * cte * (1.0e0 - eta) * (Gam2 ** 2 + Shift2 ** 2)
        )

        index_PART2 = abs(X) <= 3.0e-8 * abs(Y)
        index_PART3 = (abs(Y) <= 1.0e-15 * abs(X)) & ~index_PART2
        index_PART4 = ~(index_PART2 | index_PART3)

        # PART4
        if any(index_PART4):
            # print('PART4') # DEBUG
            X_TMP = X[index_PART4]
            Z1 = sqrt(X_TMP + Y) - csqrtY
            Z2 = Z1 + __FloatType__(2.0e0) * csqrtY
            xZ1 = -Z1.imag
            yZ1 = Z1.real
            xZ2 = -Z2.imag
            yZ2 = Z2.real
            SZ1 = sqrt(xZ1 ** 2 + yZ1 ** 2)
            SZ2 = sqrt(xZ2 ** 2 + yZ2 ** 2)
            DSZ = abs(SZ1 - SZ2)
            SZmx = maximum(SZ1, SZ2)
            SZmn = minimum(SZ1, SZ2)
            length_PART4 = len(index_PART4)
            WR1_PART4 = zeros(length_PART4)
            WI1_PART4 = zeros(length_PART4)
            WR2_PART4 = zeros(length_PART4)
            WI2_PART4 = zeros(length_PART4)
            index_CPF3 = (DSZ <= 1.0e0) & (SZmx > 8.0e0) & (SZmn <= 8.0e0)
            index_CPF = ~index_CPF3  # can be removed
            if any(index_CPF3):
                # print('PART4/CPF3') # DEBUG
                WR1, WI1 = cpf3(xZ1[index_CPF3], yZ1[index_CPF3])
                WR2, WI2 = cpf3(xZ2[index_CPF3], yZ2[index_CPF3])
                WR1_PART4[index_CPF3] = WR1
                WI1_PART4[index_CPF3] = WI1
                WR2_PART4[index_CPF3] = WR2
                WI2_PART4[index_CPF3] = WI2
            if any(index_CPF):
                # print('PART4/CPF') # DEBUG
                # print(VARIABLES['CPF'])
                # WR1,WI1 = cpf(xZ1[index_CPF],yZ1[index_CPF])
                # WR2,WI2 = cpf(xZ2[index_CPF],yZ2[index_CPF])
                WR1, WI1 = VARIABLES["CPF"](xZ1[index_CPF], yZ1[index_CPF])
                WR2, WI2 = VARIABLES["CPF"](xZ2[index_CPF], yZ2[index_CPF])
                WR1_PART4[index_CPF] = WR1
                WI1_PART4[index_CPF] = WI1
                WR2_PART4[index_CPF] = WR2
                WI2_PART4[index_CPF] = WI2

            Aterm = (
                rpi
                * cte
                * (
                    __ComplexType__(WR1_PART4 + 1.0e0j * WI1_PART4)
                    - __ComplexType__(WR2_PART4 + 1.0e0j * WI2_PART4)
                )
            )
            Bterm = (
                -1.0e0
                + rpi
                / (2.0e0 * csqrtY)
                * (1.0e0 - Z1 ** 2)
                * __ComplexType__(WR1_PART4 + 1.0e0j * WI1_PART4)
                - rpi
                / (2.0e0 * csqrtY)
                * (1.0e0 - Z2 ** 2)
                * __ComplexType__(WR2_PART4 + 1.0e0j * WI2_PART4)
            ) / c2t
            Aterm_GLOBAL[index_PART4] = Aterm
            Bterm_GLOBAL[index_PART4] = Bterm

        # PART2
        if any(index_PART2):
            # print('PART2') # DEBUG
            X_TMP = X[index_PART2]
            Z1 = (iz * (sg0 - sg[index_PART2]) + c0t) * cte
            Z2 = sqrt(X_TMP + Y) + csqrtY
            xZ1 = -Z1.imag
            yZ1 = Z1.real
            xZ2 = -Z2.imag
            yZ2 = Z2.real
            # WR1_PART2,WI1_PART2 = cpf(xZ1,yZ1)
            # WR2_PART2,WI2_PART2 = cpf(xZ2,yZ2)
            WR1_PART2, WI1_PART2 = VARIABLES["CPF"](xZ1, yZ1)
            WR2_PART2, WI2_PART2 = VARIABLES["CPF"](xZ2, yZ2)
            Aterm = (
                rpi
                * cte
                * (
                    __ComplexType__(WR1_PART2 + 1.0e0j * WI1_PART2)
                    - __ComplexType__(WR2_PART2 + 1.0e0j * WI2_PART2)
                )
            )
            Bterm = (
                -1.0e0
                + rpi
                / (2.0e0 * csqrtY)
                * (1.0e0 - Z1 ** 2)
                * __ComplexType__(WR1_PART2 + 1.0e0j * WI1_PART2)
                - rpi
                / (2.0e0 * csqrtY)
                * (1.0e0 - Z2 ** 2)
                * __ComplexType__(WR2_PART2 + 1.0e0j * WI2_PART2)
            ) / c2t
            Aterm_GLOBAL[index_PART2] = Aterm
            Bterm_GLOBAL[index_PART2] = Bterm

        # PART3
        if any(index_PART3):
            # print('PART3') # DEBUG
            X_TMP = X[index_PART3]
            xZ1 = -sqrt(X_TMP + Y).imag
            yZ1 = sqrt(X_TMP + Y).real
            # WR1_PART3,WI1_PART3 =  cpf(xZ1,yZ1)
            WR1_PART3, WI1_PART3 = VARIABLES["CPF"](xZ1, yZ1)
            index_ABS = abs(sqrt(X_TMP)) <= 4.0e3
            index_NOT_ABS = ~index_ABS
            Aterm = zeros(len(index_PART3), dtype=__ComplexType__)
            Bterm = zeros(len(index_PART3), dtype=__ComplexType__)
            if any(index_ABS):
                xXb = -sqrt(X).imag
                yXb = sqrt(X).real
                # WRb,WIb = cpf(xXb,yXb)
                WRb, WIb = VARIABLES["CPF"](xXb, yXb)
                Aterm[index_ABS] = (2.0e0 * rpi / c2t) * (
                    1.0e0 / rpi
                    - sqrt(X_TMP[index_ABS]) * __ComplexType__(WRb + 1.0e0j * WIb)
                )
                Bterm[index_ABS] = (1.0e0 / c2t) * (
                    -1.0e0
                    + 2.0e0
                    * rpi
                    * (1.0e0 - X_TMP[index_ABS] - 2.0e0 * Y)
                    * (
                        1.0e0 / rpi
                        - sqrt(X_TMP[index_ABS]) * __ComplexType__(WRb + 1.0e0j * WIb)
                    )
                    + 2.0e0
                    * rpi
                    * sqrt(X_TMP[index_ABS] + Y)
                    * __ComplexType__(WR1_PART3 + 1.0e0j * WI1_PART3)
                )
            if any(index_NOT_ABS):
                Aterm[index_NOT_ABS] = (1.0e0 / c2t) * (
                    1.0e0 / X_TMP[index_NOT_ABS] - 1.5e0 / (X_TMP[index_NOT_ABS] ** 2)
                )
                Bterm[index_NOT_ABS] = (1.0e0 / c2t) * (
                    -1.0e0
                    + (1.0e0 - X_TMP[index_NOT_ABS] - 2.0e0 * Y)
                    * (
                        1.0e0 / X_TMP[index_NOT_ABS]
                        - 1.5e0 / (X_TMP[index_NOT_ABS] ** 2)
                    )
                    + 2.0e0
                    * rpi
                    * sqrt(X_TMP[index_NOT_ABS] + Y)
                    * __ComplexType__(WR1 + 1.0e0j * WI1)
                )
            Aterm_GLOBAL[index_PART3] = Aterm
            Bterm_GLOBAL[index_PART3] = Bterm

    # common part
    LS_pCqSDHC = (1.0e0 / pi) * (
        Aterm_GLOBAL
        / (
            1.0e0
            - (anuVC - eta * (c0 - 1.5e0 * c2)) * Aterm_GLOBAL
            + eta * c2 * Bterm_GLOBAL
        )
    )
    # print('pcqsdc_end',sum(LS_pCqSDHC.real),sum(LS_pCqSDHC.imag))
    return LS_pCqSDHC.real, LS_pCqSDHC.imag


# ------------------  CROSS-SECTIONS, XSECT.PY --------------------------------

# set interfaces for TIPS(M,I,T)
def PYTIPS(M, I, T):
    return BD_TIPS_2011_PYTHON(M, I, T)[1]


# set interfaces for profiles
# PYHTP = pcqsdhc
# PROFILE_HTP = PYHTP
# PROFILE_VOIGT = lambda sg0,GamD,Gam0,sg: PROFILE_HTP(sg0,GamD,Gam0,cZero,cZero,cZero,cZero,cZero,sg)
# PROFILE_LORENTZ = lambda sg0,Gam0,sg: Gam0/(pi*(Gam0**2+(sg-sg0)**2))
# PROFILE_DOPPLER = lambda sg0,GamD,sg: cSqrtLn2divSqrtPi*exp(-cLn2*((sg-sg0)/GamD)**2)/GamD


def PROFILE_HT(sg0, GamD, Gam0, Gam2, Shift0, Shift2, anuVC, eta, sg):
    """
    #-------------------------------------------------
    #      "pCqSDHC": partially-Correlated quadratic-Speed-Dependent Hard-Collision
    #      Subroutine to Compute the complex normalized spectral shape of an 
    #      isolated line by the pCqSDHC model
    #
    #      References:
    #
    #      1) N.H. Ngo, D. Lisak, H. Tran, J.-M. Hartmann.
    #         An isolated line-shape model to go beyond the Voigt profile in 
    #         spectroscopic databases and radiative transfer codes.
    #         JQSRT, Volume 129, November 2013, Pages 89–100
    #         http://dx.doi.org/10.1016/j.jqsrt.2013.05.034
    #
    #      2) H. Tran, N.H. Ngo, J.-M. Hartmann.
    #         Efficient computation of some speed-dependent isolated line profiles.
    #         JQSRT, Volume 129, November 2013, Pages 199–203
    #         http://dx.doi.org/10.1016/j.jqsrt.2013.06.015
    #
    #      3) H. Tran, N.H. Ngo, J.-M. Hartmann.
    #         Erratum to “Efficient computation of some speed-dependent isolated line profiles”.
    #         JQSRT, Volume 134, February 2014, Pages 104
    #         http://dx.doi.org/10.1016/j.jqsrt.2013.10.015
    #
    #      Input/Output Parameters of Routine (Arguments or Common)
    #      ---------------------------------
    #      T       : Temperature in Kelvin (Input).
    #      amM1    : Molar mass of the absorber in g/mol(Input).
    #      sg0     : Unperturbed line position in cm-1 (Input).
    #      GamD    : Doppler HWHM in cm-1 (Input)
    #      Gam0    : Speed-averaged line-width in cm-1 (Input).       
    #      Gam2    : Speed dependence of the line-width in cm-1 (Input).
    #      anuVC   : Velocity-changing frequency in cm-1 (Input).
    #      eta     : Correlation parameter, No unit (Input).
    #      Shift0  : Speed-averaged line-shift in cm-1 (Input).
    #      Shift2  : Speed dependence of the line-shift in cm-1 (Input)       
    #      sg      : Current WaveNumber of the Computation in cm-1 (Input).
    #
    #      The function has two outputs:
    #      -----------------
    #      (1): Real part of the normalized spectral shape (cm)
    #      (2): Imaginary part of the normalized spectral shape (cm)
    #
    #      Called Routines: 'CPF'       (Complex Probability Function)
    #      ---------------  'CPF3'      (Complex Probability Function for the region 3)
    #
    #      Based on a double precision Fortran version
    #
    #-------------------------------------------------
    """
    return pcqsdhc(sg0, GamD, Gam0, Gam2, Shift0, Shift2, anuVC, eta, sg)


PROFILE_HTP = PROFILE_HT  # stub for backwards compatibility


def PROFILE_SDRAUTIAN(sg0, GamD, Gam0, Gam2, Shift0, Shift2, anuVC, sg):
    """
    # Speed dependent Rautian profile based on HTP.
    # Input parameters:
    #      sg0     : Unperturbed line position in cm-1 (Input).
    #      GamD    : Doppler HWHM in cm-1 (Input)
    #      Gam0    : Speed-averaged line-width in cm-1 (Input).       
    #      Gam2    : Speed dependence of the line-width in cm-1 (Input).
    #      anuVC   : Velocity-changing frequency in cm-1 (Input).
    #      Shift0  : Speed-averaged line-shift in cm-1 (Input).
    #      Shift2  : Speed dependence of the line-shift in cm-1 (Input)       
    #      sg      : Current WaveNumber of the Computation in cm-1 (Input).
    """
    return pcqsdhc(sg0, GamD, Gam0, Gam2, Shift0, Shift2, anuVC, cZero, sg)


def PROFILE_RAUTIAN(sg0, GamD, Gam0, Shift0, anuVC, eta, sg):
    """
    # Rautian profile based on HTP.
    # Input parameters:
    #      sg0     : Unperturbed line position in cm-1 (Input).
    #      GamD    : Doppler HWHM in cm-1 (Input)
    #      Gam0    : Speed-averaged line-width in cm-1 (Input).       
    #      anuVC   : Velocity-changing frequency in cm-1 (Input).
    #      Shift0  : Speed-averaged line-shift in cm-1 (Input).
    #      sg      : Current WaveNumber of the Computation in cm-1 (Input).
    """
    return pcqsdhc(sg0, GamD, Gam0, cZero, Shift0, cZero, anuVC, cZero, sg)


def PROFILE_SDVOIGT(sg0, GamD, Gam0, Gam2, Shift0, Shift2, sg):
    """
    # Speed dependent Voigt profile based on HTP.
    # Input parameters:
    #      sg0     : Unperturbed line position in cm-1 (Input).
    #      GamD    : Doppler HWHM in cm-1 (Input)
    #      Gam0    : Speed-averaged line-width in cm-1 (Input).       
    #      Gam2    : Speed dependence of the line-width in cm-1 (Input).
    #      Shift0  : Speed-averaged line-shift in cm-1 (Input).
    #      Shift2  : Speed dependence of the line-shift in cm-1 (Input)       
    #      sg      : Current WaveNumber of the Computation in cm-1 (Input).
    """
    return pcqsdhc(sg0, GamD, Gam0, Gam2, Shift0, Shift2, cZero, cZero, sg)


def PROFILE_VOIGT(sg0, GamD, Gam0, sg):
    """
    # Voigt profile based on HTP.
    # Input parameters:
    #   sg0: Unperturbed line position in cm-1 (Input).
    #   GamD: Doppler HWHM in cm-1 (Input)
    #   Gam0: Speed-averaged line-width in cm-1 (Input).       
    #   sg: Current WaveNumber of the Computation in cm-1 (Input).
    """
    return PROFILE_HTP(sg0, GamD, Gam0, cZero, cZero, cZero, cZero, cZero, sg)


def PROFILE_LORENTZ(sg0, Gam0, sg):
    """
    # Lorentz profile.
    # Input parameters:
    #   sg0: Unperturbed line position in cm-1 (Input).
    #   Gam0: Speed-averaged line-width in cm-1 (Input).       
    #   sg: Current WaveNumber of the Computation in cm-1 (Input).
    """
    return Gam0 / (pi * (Gam0 ** 2 + (sg - sg0) ** 2))


def PROFILE_DOPPLER(sg0, GamD, sg):
    """
    # Doppler profile.
    # Input parameters:
    #   sg0: Unperturbed line position in cm-1 (Input).
    #   GamD: Doppler HWHM in cm-1 (Input)
    #   sg: Current WaveNumber of the Computation in cm-1 (Input).
    """
    return cSqrtLn2divSqrtPi * exp(-cLn2 * ((sg - sg0) / GamD) ** 2) / GamD


# Volume concentration of all gas molecules at the pressure p and temperature T


def volumeConcentration(p, T):
    return (p / 9.869233e-7) / (cBolts * T)  # CGS


# ------------------------------- PARAMETER DEPENDENCIES --------------------------------

# temperature dependence for intencities (HITRAN)


def EnvironmentDependency_Intensity(
    LineIntensityRef, T, Tref, SigmaT, SigmaTref, LowerStateEnergy, LineCenter
):
    const = __FloatType__(1.4388028496642257)
    ch = exp(-const * LowerStateEnergy / T) * (1 - exp(-const * LineCenter / T))
    zn = exp(-const * LowerStateEnergy / Tref) * (1 - exp(-const * LineCenter / Tref))
    LineIntensity = LineIntensityRef * SigmaTref / SigmaT * ch / zn
    return LineIntensity


# environmental dependence for GammaD (HTP, Voigt)    # Tref/T ????


def EnvironmentDependency_GammaD(GammaD_ref, T, Tref):
    # Doppler parameters do not depend on pressure!
    return GammaD_ref * sqrt(T / Tref)


# environmental dependence for Gamma0 (HTP, Voigt)


def EnvironmentDependency_Gamma0(Gamma0_ref, T, Tref, p, pref, TempRatioPower):
    return Gamma0_ref * p / pref * (Tref / T) ** TempRatioPower


# environmental dependence for Gamma2 (HTP)


def EnvironmentDependency_Gamma2(Gamma2_ref, T, Tref, p, pref, TempRatioPower):
    return Gamma2_ref * p / pref * (Tref / T) ** TempRatioPower


# environmental dependence for Delta0 (HTP)


def EnvironmentDependency_Delta0(Delta0_ref, p, pref):
    return Delta0_ref * p / pref


# environmental dependence for Delta2 (HTP)


def EnvironmentDependency_Delta2(Delta2_ref, p, pref):
    return Delta2_ref * p / pref


# environmental dependence for anuVC (HTP)


def EnvironmentDependency_anuVC(anuVC_ref, T, Tref, p, pref):
    return anuVC_ref * Tref / T * p / pref


# ------------------------------- /PARAMETER DEPENDENCIES --------------------------------

# ------------------------------- BINGINGS --------------------------------


# default parameter bindings
DefaultParameterBindings = {}

# default temperature dependencies
DefaultEnvironmentDependencyBindings = {}

# ------------------------------- /BINGINGS --------------------------------

# default values for intensity threshold
DefaultIntensityThreshold = 0.0  # cm*molec

# default value for omega wing in halfwidths (from center)
DefaultOmegaWingHW = 50.0  # cm-1    HOTW default


# check and argument for being a tuple or list
# this is connected with a "bug" that in Python
# (val) is not a tuple, but (val,) is a tuple
def listOfTuples(a):
    if type(a) not in set([list, tuple]):
        a = [a]
    return a


# determine default parameters from those which are passed to absorptionCoefficient_...
def getDefaultValuesForXsect(
    Components,
    SourceTables,
    Environment,
    OmegaRange,
    OmegaStep,
    OmegaWing,
    IntensityThreshold,
    Format,
):
    if SourceTables[0] == None:
        SourceTables = [
            "__BUFFER__",
        ]
    if Environment == None:
        Environment = {"T": 296.0, "p": 1.0}
    if Components == [None]:
        CompDict = {}
        for TableName in SourceTables:
            # check table existance
            if TableName not in list(LOCAL_TABLE_CACHE.keys()):
                raise Exception(
                    "%s: no such table. Check tableList() for more info." % TableName
                )
            mol_ids = LOCAL_TABLE_CACHE[TableName]["data"]["molec_id"]
            iso_ids = LOCAL_TABLE_CACHE[TableName]["data"]["local_iso_id"]
            if len(mol_ids) != len(iso_ids):
                raise Exception("Lengths if mol_ids and iso_ids differ!")
            MI_zip = list(zip(mol_ids, iso_ids))
            MI_zip = set(MI_zip)
            for mol_id, iso_id in MI_zip:
                CompDict[(mol_id, iso_id)] = None
        Components = list(CompDict.keys())
    if OmegaRange == None:
        omega_min = float("inf")
        omega_max = float("-inf")
        for TableName in SourceTables:
            nu = LOCAL_TABLE_CACHE[TableName]["data"]["nu"]
            numin = min(nu)
            numax = max(nu)
            if omega_min > numin:
                omega_min = numin
            if omega_max < numax:
                omega_max = numax
        OmegaRange = (omega_min, omega_max)
    OmegaDelta = OmegaRange[1] - OmegaRange[0]
    if OmegaStep == None:
        # OmegaStep = OmegaDelta/100.
        OmegaStep = 0.01  # cm-1
    if OmegaWing == None:
        # OmegaWing = OmegaDelta/10.
        OmegaWing = 0.0  # cm-1
    if not Format:
        """
        Infinitesimal = 1e-14 # put this to header in next version!
        min_number_of_digits = 4 # minimal number of digits after dec. pnt.
        last_digit_pos = 0
        while modf(OmegaStep * 10**last_digit_pos)[0] > Infinitesimal:
            last_digit_pos += 1
        actual_number_of_digits = max(min_number_of_digits,last_digit_pos)
        Format = '%%.%df %%e' % actual_number_of_digits
        """
        Format = "%.12f %e"
    return (
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    )


# save numpy arrays to file
# arrays must have same dimensions
def save_to_file(fname, fformat, *arg):
    f = open(fname, "w")
    for i in range(len(arg[0])):
        argline = []
        for j in range(len(arg)):
            argline.append(arg[j][i])
        f.write((fformat + "\n") % tuple(argline))
    f.close()


# ==========================================================================================
# =========================== NEW ABSORPTION COEFFICIENT ===================================
# ==========================================================================================


def absorptionCoefficient_HT(
    Components=None,
    SourceTables=None,
    partitionFunction=PYTIPS,
    Environment=None,
    OmegaRange=None,
    OmegaStep=None,
    OmegaWing=None,
    IntensityThreshold=DefaultIntensityThreshold,
    OmegaWingHW=DefaultOmegaWingHW,
    GammaL="gamma_air",
    HITRAN_units=True,
    LineShift=True,
    File=None,
    Format=None,
    OmegaGrid=None,
    WavenumberRange=None,
    WavenumberStep=None,
    WavenumberWing=None,
    WavenumberWingHW=None,
    WavenumberGrid=None,
    Diluent={},
    EnvDependences=None,
):
    """
    INPUT PARAMETERS: 
        Components:  list of tuples [(M,I,D)], where
                        M - HITRAN molecule number,
                        I - HITRAN isotopologue number,
                        D - relative abundance (optional)
        SourceTables:  list of tables from which to calculate cross-section   (optional)
        partitionFunction:  pointer to partition function (default is PYTIPS) (optional)
        Environment:  dictionary containing thermodynamic parameters.
                        'p' - pressure in atmospheres,
                        'T' - temperature in Kelvin
                        Default={'p':1.,'T':296.}
        WavenumberRange:  wavenumber range to consider.
        WavenumberStep:   wavenumber step to consider. 
        WavenumberWing:   absolute wing for calculating a lineshape (in cm-1) 
        WavenumberWingHW:  relative wing for calculating a lineshape (in halfwidths)
        IntensityThreshold:  threshold for intensities
        GammaL:  specifies broadening parameter ('gamma_air' or 'gamma_self')
        HITRAN_units:  use cm2/molecule (True) or cm-1 (False) for absorption coefficient
        File:   write output to file (if specified)
        Format:  c-format of file output (accounts for significant digits in WavenumberStep)
    OUTPUT PARAMETERS: 
        Wavenum: wavenumber grid with respect to parameters WavenumberRange and WavenumberStep
        Xsect: absorption coefficient calculated on the grid
    ---
    DESCRIPTION:
        Calculate absorption coefficient using HT profile.
        Absorption coefficient is calculated at arbitrary temperature and pressure.
        User can vary a wide range of parameters to control a process of calculation.
        The choise of these parameters depends on properties of a particular linelist.
        Default values are a sort of guess which gives a decent precision (on average) 
        for a reasonable amount of cpu time. To increase calculation accuracy,
        user should use a trial and error method.
    ---
    EXAMPLE OF USAGE:
        nu,coef = absorptionCoefficient_HT(((2,1),),'co2',WavenumberStep=0.01,
                                              HITRAN_units=False,GammaL='gamma_self')
    ---
    """

    warn("To get the most up-to-date version please check http://hitran.org/hapi")

    # Paremeters OmegaRange,OmegaStep,OmegaWing,OmegaWingHW, and OmegaGrid
    # are deprecated and given for backward compatibility with the older versions.
    if WavenumberRange:
        OmegaRange = WavenumberRange
    if WavenumberStep:
        OmegaStep = WavenumberStep
    if WavenumberWing:
        OmegaWing = WavenumberWing
    if WavenumberWingHW:
        OmegaWingHW = WavenumberWingHW
    if WavenumberGrid:
        OmegaGrid = WavenumberGrid

    # "bug" with 1-element list
    Components = listOfTuples(Components)
    SourceTables = listOfTuples(SourceTables)

    # determine final input values
    (
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    ) = getDefaultValuesForXsect(
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    )

    # warn user about too large omega step
    if OmegaStep > 0.1:
        warn("Big wavenumber step: possible accuracy decline")

    # get uniform linespace for cross-section
    # number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    # Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        # Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.0)  # K
    pref = __FloatType__(1.0)  # atm

    # actual temperature and pressure
    T = Environment["T"]  # K
    p = Environment["p"]  # atm

    # Find reference temperature
    TRanges = [(0, 100), (100, 200), (200, 400), (400, float("inf"))]
    Trefs = [50.0, 150.0, 296.0, 700.0]
    for TRange, TrefHT in zip(TRanges, Trefs):
        if T >= TRange[0] and T < TRange[1]:
            break

    # create dictionary from Components
    ABUNDANCES = {}
    NATURAL_ABUNDANCES = {}
    for Component in Components:
        M = Component[0]
        I = Component[1]
        if len(Component) >= 3:
            ni = Component[2]
        else:
            try:
                ni = ISO[(M, I)][ISO_INDEX["abundance"]]
            except KeyError:
                raise Exception("cannot find component M,I = %d,%d." % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX["abundance"]]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # setup the default empty environment dependence function
    if not EnvDependences:

        def EnvDependences(ENV, LINE):
            return {}

    Env = Environment.copy()
    Env["Tref"] = Tref
    Env["pref"] = pref

    # setup the Diluent variable
    GammaL = GammaL.lower()
    if not Diluent:
        if GammaL == "gamma_air":
            Diluent = {"air": 1.0}
        elif GammaL == "gamma_self":
            Diluent = {"self": 1.0}
        else:
            raise Exception("Unknown GammaL value: %s" % GammaL)

    # Simple check
    for key in Diluent:
        val = Diluent[key]
        if val < 0 and val > 1:
            raise Exception("Diluent fraction must be in [0,1]")

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get the number of rows
        nline = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]

        # get parameter names for each table
        parnames = list(LOCAL_TABLE_CACHE[TableName]["data"].keys())

        # loop through line centers (single stream)
        for RowID in range(nline):

            # Get the custom environment dependences
            Line = {
                parname: LOCAL_TABLE_CACHE[TableName]["data"][parname][RowID]
                for parname in parnames
            }
            CustomEnvDependences = EnvDependences(Env, Line)

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]["data"]["nu"][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]["data"]["sw"][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]["data"]["elower"][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["molec_id"][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["local_iso_id"][RowID]

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES:
                continue

            # partition functions for T and Tref
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            if "sw" in CustomEnvDependences:
                LineIntensity = CustomEnvDependences["sw"]
            else:
                LineIntensity = EnvironmentDependency_Intensity(
                    LineIntensityDB,
                    T,
                    Tref,
                    SigmaT,
                    SigmaTref,
                    LowerStateEnergyDB,
                    LineCenterDB,
                )

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            if LineIntensity < IntensityThreshold:
                continue

            #   doppler broadening coefficient (GammaD)
            cMassMol = 1.66053873e-27  # hapi
            m = molecularMass(MoleculeNumberDB, IsoNumberDB) * cMassMol * 1000
            GammaD = sqrt(2 * cBolts * T * log(2) / m / cc ** 2) * LineCenterDB

            #   pressure broadening coefficient
            Gamma0 = 0.0
            Shift0 = 0.0
            Gamma2 = 0.0
            Shift2 = 0.0
            Delta2 = 0.0
            NuVC = 0.0
            EtaNumer = 0.0
            for species in Diluent:
                species_lower = species.lower()

                abun = Diluent[species]

                # Search for broadening HWHM.
                try:
                    # search for HT-style name
                    Gamma0DB = LOCAL_TABLE_CACHE[TableName]["data"][
                        "gamma_HT_0_%s_%d" % (species_lower, TrefHT)
                    ][RowID]
                    if Gamma0DB == 0.0:
                        raise Exception
                except:
                    try:
                        # search for Voigt-style name
                        Gamma0DB = LOCAL_TABLE_CACHE[TableName]["data"][
                            "gamma_%s" % species_lower
                        ][RowID]
                    except:
                        Gamma0DB = 0.0

                # Search for temperature exponent for broadening HWHM.
                try:
                    # search for HT-style name
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][
                        "n_HT_%s_%d" % (species_lower, TrefHT)
                    ][RowID]
                    if TempRatioPowerDB == 0.0:
                        raise Exception
                    Tref = TrefHT
                except:
                    Tref = 296.0
                    try:
                        # search for Voigt-style name
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][
                            "n_%s" % species_lower
                        ][RowID]
                        if species_lower == "self" and TempRatioPowerDB == 0.0:
                            # same for self as for air
                            TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][
                                "n_air"
                            ][RowID]
                    except:
                        # print('TempRatioPowerDB is set to zero')
                        # TempRatioPowerDB = 0
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][
                            "n_air"
                        ][RowID]

                # Add to the final Gamma0
                Gamma0T = CustomEnvDependences.get(
                    "gamma_HT_0_%s_%d" % (species_lower, TrefHT),
                    CustomEnvDependences.get(
                        "gamma_%s" % species_lower,
                        EnvironmentDependency_Gamma0(
                            Gamma0DB, T, Tref, p, pref, TempRatioPowerDB
                        ),
                    ),
                )
                Gamma0 += abun * Gamma0T

                # Search for shift.
                try:
                    # search for HT-style name
                    Shift0DB = LOCAL_TABLE_CACHE[TableName]["data"][
                        "delta_HT_0_%s_%d" % (species_lower, TrefHT)
                    ][RowID]
                    if Shift0DB == 0.0:
                        raise Exception
                except:
                    try:
                        # search for Voigt-style name
                        Shift0DB = LOCAL_TABLE_CACHE[TableName]["data"][
                            "delta_%s" % species_lower
                        ][RowID]
                    except:
                        Shift0DB = 0.0

                # Search for temperature dependence for shift.
                try:
                    # search for HT-style name
                    deltap = LOCAL_TABLE_CACHE[TableName]["data"][
                        "deltap_HT_%s_%d" % (species_lower, TrefHT)
                    ][RowID]
                    if deltap == 0.0:
                        raise Exception
                    Tref = TrefHT
                except:
                    Tref = 296.0
                    try:
                        # search for Voigt-style name
                        deltap = LOCAL_TABLE_CACHE[TableName]["data"][
                            "deltap_%s" % species_lower
                        ][RowID]
                    except:
                        deltap = 0.0

                Shift0T = CustomEnvDependences.get(
                    "deltap_HT_%s_%d" % (species_lower, TrefHT),
                    CustomEnvDependences.get(
                        "deltap_%s" % species_lower,
                        ((Shift0DB + deltap * (T - Tref)) * p / pref),
                    ),
                )
                Shift0 += abun * Shift0T

                # Search for speed dependence for HWHM.
                try:
                    Gamma2DB = LOCAL_TABLE_CACHE[TableName]["data"][
                        "gamma_HT_2_%s_%d" % (species_lower, TrefHT)
                    ][RowID]
                    if Gamma2DB == 0.0:
                        raise Exception
                except:
                    try:
                        SDDB = LOCAL_TABLE_CACHE[TableName]["data"][
                            "SD_%s" % species_lower
                        ][RowID]
                        Gamma2DB = SDDB * Gamma0DB
                    except:
                        Gamma2DB = 0.0

                Gamma2 += abun * CustomEnvDependences.get(
                    "gamma_HT_2_%s_%d" % (species_lower, TrefHT), Gamma2DB * (p / pref)
                )

                # Search for speed dependence for shift.
                try:
                    Delta2DB = LOCAL_TABLE_CACHE[TableName]["data"][
                        "delta_HT_2_%s_%d" % (species_lower, TrefHT)
                    ][RowID]
                except:
                    Delta2DB = 0.0

                Delta2 += abun * CustomEnvDependences.get(
                    "delta_HT_2_%s_%d" % (species_lower, TrefHT), Delta2DB * p / pref
                )

                # Search for frequency of VC
                try:
                    NuVCDB = LOCAL_TABLE_CACHE[TableName]["data"][
                        "nu_HT_%s" % species_lower
                    ][RowID]
                except:
                    NuVCDB = 0.0

                # Search for temperature exponent for frequency of VC
                try:
                    KappaDB = LOCAL_TABLE_CACHE[TableName]["data"][
                        "kappa_HT_%s" % species_lower
                    ][RowID]
                except:
                    KappaDB = 0.0

                NuVC += abun * CustomEnvDependences.get(
                    "nu_HT_%s" % species_lower, NuVCDB * (Tref / T) ** KappaDB * p
                )

                # Setup correlation parameter
                try:
                    EtaDB = LOCAL_TABLE_CACHE[TableName]["data"][
                        "eta_HT_%s" % species_lower
                    ][RowID]
                except:
                    EtaDB = 0.0

                EtaNumer += EtaDB * abun * (Gamma0T + 1j * Shift0T)

            Eta = EtaNumer / (Gamma0 + 1j * Shift0)

            #   get final wing of the line according to Gamma0, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW * Gamma0, OmegaWingHW * GammaD)

            #   shift coefficient
            # Shift0 = Shift0DB*p/pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB - OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB + OmegaWingF)
            lineshape_vals = PROFILE_HT(
                LineCenterDB,
                GammaD,
                Gamma0,
                Gamma2,
                Shift0,
                Shift2,
                NuVC,
                Eta,
                Omegas[BoundIndexLower:BoundIndexUpper],
            )[0]
            Xsect[BoundIndexLower:BoundIndexUpper] += (
                factor
                / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * LineIntensity
                * lineshape_vals
            )
            # print(LineCenterDB,GammaD,Gamma0,Gamma2,Shift0,sum(lineshape_vals),sum(Xsect))
            # raise Exception

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def absorptionCoefficient_SDVoigt(
    Components=None,
    SourceTables=None,
    partitionFunction=PYTIPS,
    Environment=None,
    OmegaRange=None,
    OmegaStep=None,
    OmegaWing=None,
    IntensityThreshold=DefaultIntensityThreshold,
    OmegaWingHW=DefaultOmegaWingHW,
    GammaL="gamma_air",
    HITRAN_units=True,
    LineShift=True,
    File=None,
    Format=None,
    OmegaGrid=None,
    WavenumberRange=None,
    WavenumberStep=None,
    WavenumberWing=None,
    WavenumberWingHW=None,
    WavenumberGrid=None,
    Diluent={},
    EnvDependences=None,
):
    """
    INPUT PARAMETERS: 
        Components:  list of tuples [(M,I,D)], where
                        M - HITRAN molecule number,
                        I - HITRAN isotopologue number,
                        D - relative abundance (optional)
        SourceTables:  list of tables from which to calculate cross-section   (optional)
        partitionFunction:  pointer to partition function (default is PYTIPS) (optional)
        Environment:  dictionary containing thermodynamic parameters.
                        'p' - pressure in atmospheres,
                        'T' - temperature in Kelvin
                        Default={'p':1.,'T':296.}
        WavenumberRange:  wavenumber range to consider.
        WavenumberStep:   wavenumber step to consider. 
        WavenumberWing:   absolute wing for calculating a lineshape (in cm-1) 
        WavenumberWingHW:  relative wing for calculating a lineshape (in halfwidths)
        IntensityThreshold:  threshold for intensities
        GammaL:  specifies broadening parameter ('gamma_air' or 'gamma_self')
        HITRAN_units:  use cm2/molecule (True) or cm-1 (False) for absorption coefficient
        File:   write output to file (if specified)
        Format:  c-format of file output (accounts for significant digits in WavenumberStep)
    OUTPUT PARAMETERS: 
        Wavenum: wavenumber grid with respect to parameters WavenumberRange and WavenumberStep
        Xsect: absorption coefficient calculated on the grid
    ---
    DESCRIPTION:
        Calculate absorption coefficient using SDVoigt profile.
        Absorption coefficient is calculated at arbitrary temperature and pressure.
        User can vary a wide range of parameters to control a process of calculation.
        The choise of these parameters depends on properties of a particular linelist.
        Default values are a sort of guess which gives a decent precision (on average) 
        for a reasonable amount of cpu time. To increase calculation accuracy,
        user should use a trial and error method.
    ---
    EXAMPLE OF USAGE:
        nu,coef = absorptionCoefficient_SDVoigt(((2,1),),'co2',WavenumberStep=0.01,
                                              HITRAN_units=False,GammaL='gamma_self')
    ---
    """

    warn("To get the most up-to-date version please check http://hitran.org/hapi")

    # Paremeters OmegaRange,OmegaStep,OmegaWing,OmegaWingHW, and OmegaGrid
    # are deprecated and given for backward compatibility with the older versions.
    if WavenumberRange:
        OmegaRange = WavenumberRange
    if WavenumberStep:
        OmegaStep = WavenumberStep
    if WavenumberWing:
        OmegaWing = WavenumberWing
    if WavenumberWingHW:
        OmegaWingHW = WavenumberWingHW
    if WavenumberGrid:
        OmegaGrid = WavenumberGrid

    # "bug" with 1-element list
    Components = listOfTuples(Components)
    SourceTables = listOfTuples(SourceTables)

    # determine final input values
    (
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    ) = getDefaultValuesForXsect(
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    )

    # warn user about too large omega step
    if OmegaStep > 0.1:
        warn("Big wavenumber step: possible accuracy decline")

    # get uniform linespace for cross-section
    # number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    # Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        # Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.0)  # K
    pref = __FloatType__(1.0)  # atm

    # actual temperature and pressure
    T = Environment["T"]  # K
    p = Environment["p"]  # atm

    # create dictionary from Components
    ABUNDANCES = {}
    NATURAL_ABUNDANCES = {}
    for Component in Components:
        M = Component[0]
        I = Component[1]
        if len(Component) >= 3:
            ni = Component[2]
        else:
            try:
                ni = ISO[(M, I)][ISO_INDEX["abundance"]]
            except KeyError:
                raise Exception("cannot find component M,I = %d,%d." % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX["abundance"]]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # setup the default empty environment dependence function
    if not EnvDependences:

        def EnvDependences(ENV, LINE):
            return {}

    Env = Environment.copy()
    Env["Tref"] = Tref
    Env["pref"] = pref

    # setup the Diluent variable
    GammaL = GammaL.lower()
    if not Diluent:
        if GammaL == "gamma_air":
            Diluent = {"air": 1.0}
        elif GammaL == "gamma_self":
            Diluent = {"self": 1.0}
        else:
            raise Exception("Unknown GammaL value: %s" % GammaL)

    # Simple check
    for key in Diluent:
        val = Diluent[key]
        if val < 0 and val > 1:
            raise Exception("Diluent fraction must be in [0,1]")

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get the number of rows
        nline = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]

        # get parameter names for each table
        parnames = list(LOCAL_TABLE_CACHE[TableName]["data"].keys())

        # loop through line centers (single stream)
        for RowID in range(nline):

            # Get the custom environment dependences
            Line = {
                parname: LOCAL_TABLE_CACHE[TableName]["data"][parname][RowID]
                for parname in parnames
            }
            CustomEnvDependences = EnvDependences(Env, Line)

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]["data"]["nu"][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]["data"]["sw"][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]["data"]["elower"][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["molec_id"][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["local_iso_id"][RowID]

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES:
                continue

            # partition functions for T and Tref
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            if "sw" in CustomEnvDependences:
                LineIntensity = CustomEnvDependences["sw"]
            else:
                LineIntensity = EnvironmentDependency_Intensity(
                    LineIntensityDB,
                    T,
                    Tref,
                    SigmaT,
                    SigmaTref,
                    LowerStateEnergyDB,
                    LineCenterDB,
                )

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            if LineIntensity < IntensityThreshold:
                continue

            #   doppler broadening coefficient (GammaD)
            cMassMol = 1.66053873e-27  # hapi
            m = molecularMass(MoleculeNumberDB, IsoNumberDB) * cMassMol * 1000
            GammaD = sqrt(2 * cBolts * T * log(2) / m / cc ** 2) * LineCenterDB

            #   pressure broadening coefficient
            Gamma0 = 0.0
            Shift0 = 0.0
            Gamma2 = 0.0
            Shift2 = 0.0
            for species in Diluent:
                species_lower = species.lower()

                abun = Diluent[species]

                gamma_name = "gamma_" + species_lower
                try:
                    Gamma0DB = LOCAL_TABLE_CACHE[TableName]["data"][gamma_name][RowID]
                except:
                    Gamma0DB = 0.0

                n_name = "n_" + species_lower
                try:
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][n_name][
                        RowID
                    ]
                    if species_lower == "self" and TempRatioPowerDB == 0.0:
                        # same for self as for air
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][
                            "n_air"
                        ][RowID]
                except:
                    # TempRatioPowerDB = 0
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"]["n_air"][
                        RowID
                    ]

                # Add to the final Gamma0
                Gamma0 += abun * CustomEnvDependences.get(
                    gamma_name,  # default ->
                    EnvironmentDependency_Gamma0(
                        Gamma0DB, T, Tref, p, pref, TempRatioPowerDB
                    ),
                )

                delta_name = "delta_" + species_lower
                try:
                    Shift0DB = LOCAL_TABLE_CACHE[TableName]["data"][delta_name][RowID]
                except:
                    Shift0DB = 0.0

                deltap_name = "deltap_" + species_lower
                try:
                    deltap = LOCAL_TABLE_CACHE[TableName]["data"][deltap_name][RowID]
                except:
                    deltap = 0.0

                Shift0 += abun * CustomEnvDependences.get(
                    delta_name,  # default ->
                    ((Shift0DB + deltap * (T - Tref)) * p / pref),
                )

                SD_name = "SD_" + species_lower
                try:
                    SDDB = LOCAL_TABLE_CACHE[TableName]["data"][SD_name][RowID]
                except:
                    SDDB = 0.0

                Gamma2 += (
                    abun
                    * CustomEnvDependences.get(SD_name, SDDB * p / pref)  # default ->
                    * Gamma0DB
                )

                #   get final wing of the line according to Gamma0, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW * Gamma0, OmegaWingHW * GammaD)

            #   shift coefficient
            # Shift0 = Shift0DB*p/pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB - OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB + OmegaWingF)
            lineshape_vals = PROFILE_SDVOIGT(
                LineCenterDB,
                GammaD,
                Gamma0,
                Gamma2,
                Shift0,
                Shift2,
                Omegas[BoundIndexLower:BoundIndexUpper],
            )[0]
            Xsect[BoundIndexLower:BoundIndexUpper] += (
                factor
                / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * LineIntensity
                * lineshape_vals
            )
            # print(LineCenterDB,GammaD,Gamma0,Gamma2,Shift0,sum(lineshape_vals),sum(Xsect))
            # raise Exception

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def absorptionCoefficient_Voigt(
    Components=None,
    SourceTables=None,
    partitionFunction=PYTIPS,
    Environment=None,
    OmegaRange=None,
    OmegaStep=None,
    OmegaWing=None,
    IntensityThreshold=DefaultIntensityThreshold,
    OmegaWingHW=DefaultOmegaWingHW,
    GammaL="gamma_air",
    HITRAN_units=True,
    LineShift=True,
    File=None,
    Format=None,
    OmegaGrid=None,
    WavenumberRange=None,
    WavenumberStep=None,
    WavenumberWing=None,
    WavenumberWingHW=None,
    WavenumberGrid=None,
    Diluent={},
    EnvDependences=None,
):
    """
    INPUT PARAMETERS: 
        Components:  list of tuples [(M,I,D)], where
                        M - HITRAN molecule number,
                        I - HITRAN isotopologue number,
                        D - relative abundance (optional)
        SourceTables:  list of tables from which to calculate cross-section   (optional)
        partitionFunction:  pointer to partition function (default is PYTIPS) (optional)
        Environment:  dictionary containing thermodynamic parameters.
                        'p' - pressure in atmospheres,
                        'T' - temperature in Kelvin
                        Default={'p':1.,'T':296.}
        WavenumberRange:  wavenumber range to consider.
        WavenumberStep:   wavenumber step to consider. 
        WavenumberWing:   absolute wing for calculating a lineshape (in cm-1) 
        WavenumberWingHW:  relative wing for calculating a lineshape (in halfwidths)
        IntensityThreshold:  threshold for intensities
        GammaL:  specifies broadening parameter ('gamma_air' or 'gamma_self')
        HITRAN_units:  use cm2/molecule (True) or cm-1 (False) for absorption coefficient
        File:   write output to file (if specified)
        Format:  c-format of file output (accounts for significant digits in WavenumberStep)
    OUTPUT PARAMETERS: 
        Wavenum: wavenumber grid with respect to parameters WavenumberRange and WavenumberStep
        Xsect: absorption coefficient calculated on the grid
    ---
    DESCRIPTION:
        Calculate absorption coefficient using Voigt profile.
        Absorption coefficient is calculated at arbitrary temperature and pressure.
        User can vary a wide range of parameters to control a process of calculation.
        The choise of these parameters depends on properties of a particular linelist.
        Default values are a sort of guess which gives a decent precision (on average) 
        for a reasonable amount of cpu time. To increase calculation accuracy,
        user should use a trial and error method.
    ---
    EXAMPLE OF USAGE:
        nu,coef = absorptionCoefficient_Voigt(((2,1),),'co2',WavenumberStep=0.01,
                                              HITRAN_units=False,GammaL='gamma_self')
    ---
    """

    # Paremeters OmegaRange,OmegaStep,OmegaWing,OmegaWingHW, and OmegaGrid
    # are deprecated and given for backward compatibility with the older versions.
    if WavenumberRange:
        OmegaRange = WavenumberRange
    if WavenumberStep:
        OmegaStep = WavenumberStep
    if WavenumberWing:
        OmegaWing = WavenumberWing
    if WavenumberWingHW:
        OmegaWingHW = WavenumberWingHW
    if WavenumberGrid:
        OmegaGrid = WavenumberGrid

    # "bug" with 1-element list
    Components = listOfTuples(Components)
    SourceTables = listOfTuples(SourceTables)

    # determine final input values
    (
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    ) = getDefaultValuesForXsect(
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    )

    # warn user about too large omega step
    if OmegaStep > 0.1:
        warn("Big wavenumber step: possible accuracy decline")

    # get uniform linespace for cross-section
    # number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    # Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        # Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.0)  # K
    pref = __FloatType__(1.0)  # atm

    # actual temperature and pressure
    T = Environment["T"]  # K
    p = Environment["p"]  # atm

    # create dictionary from Components
    ABUNDANCES = {}
    NATURAL_ABUNDANCES = {}
    for Component in Components:
        M = Component[0]
        I = Component[1]
        if len(Component) >= 3:
            ni = Component[2]
        else:
            try:
                ni = ISO[(M, I)][ISO_INDEX["abundance"]]
            except KeyError:
                raise Exception("cannot find component M,I = %d,%d." % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX["abundance"]]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # setup the default empty environment dependence function
    if not EnvDependences:

        def EnvDependences(ENV, LINE):
            return {}

    Env = Environment.copy()
    Env["Tref"] = Tref
    Env["pref"] = pref

    # setup the Diluent variable
    GammaL = GammaL.lower()
    if not Diluent:
        if GammaL == "gamma_air":
            Diluent = {"air": 1.0}
        elif GammaL == "gamma_self":
            Diluent = {"self": 1.0}
        else:
            raise Exception("Unknown GammaL value: %s" % GammaL)

    # Simple check
    for key in Diluent:
        val = Diluent[key]
        if val < 0 and val > 1:
            raise Exception("Diluent fraction must be in [0,1]")

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get the number of rows
        nline = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]

        # get parameter names for each table
        parnames = list(LOCAL_TABLE_CACHE[TableName]["data"].keys())

        # loop through line centers (single stream)
        for RowID in range(nline):

            # Get the custom environment dependences
            Line = {
                parname: LOCAL_TABLE_CACHE[TableName]["data"][parname][RowID]
                for parname in parnames
            }
            CustomEnvDependences = EnvDependences(Env, Line)

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]["data"]["nu"][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]["data"]["sw"][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]["data"]["elower"][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["molec_id"][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["local_iso_id"][RowID]

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES:
                continue

            # partition functions for T and Tref
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            if "sw" in CustomEnvDependences:
                LineIntensity = CustomEnvDependences["sw"]
            else:
                LineIntensity = EnvironmentDependency_Intensity(
                    LineIntensityDB,
                    T,
                    Tref,
                    SigmaT,
                    SigmaTref,
                    LowerStateEnergyDB,
                    LineCenterDB,
                )

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            if LineIntensity < IntensityThreshold:
                continue

            #   doppler broadening coefficient (GammaD)
            cMassMol = 1.66053873e-27  # hapi
            m = molecularMass(MoleculeNumberDB, IsoNumberDB) * cMassMol * 1000
            GammaD = sqrt(2 * cBolts * T * log(2) / m / cc ** 2) * LineCenterDB

            #   pressure broadening coefficient
            Gamma0 = 0.0
            Shift0 = 0.0
            for species in Diluent:
                species_lower = species.lower()

                abun = Diluent[species]

                gamma_name = "gamma_" + species_lower
                try:
                    Gamma0DB = LOCAL_TABLE_CACHE[TableName]["data"][gamma_name][RowID]
                except:
                    Gamma0DB = 0.0

                n_name = "n_" + species_lower
                try:
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][n_name][
                        RowID
                    ]
                    if species_lower == "self" and TempRatioPowerDB == 0.0:
                        # same for self as for air
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][
                            "n_air"
                        ][RowID]
                except:
                    # TempRatioPowerDB = 0
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"]["n_air"][
                        RowID
                    ]

                # Add to the final Gamma0
                Gamma0 += abun * CustomEnvDependences.get(
                    gamma_name,  # default ->
                    EnvironmentDependency_Gamma0(
                        Gamma0DB, T, Tref, p, pref, TempRatioPowerDB
                    ),
                )

                delta_name = "delta_" + species_lower
                try:
                    Shift0DB = LOCAL_TABLE_CACHE[TableName]["data"][delta_name][RowID]
                except:
                    Shift0DB = 0.0

                deltap_name = "deltap_" + species_lower
                try:
                    deltap = LOCAL_TABLE_CACHE[TableName]["data"][deltap_name][RowID]
                except:
                    deltap = 0.0

                Shift0 += abun * CustomEnvDependences.get(
                    delta_name,  # default ->
                    ((Shift0DB + deltap * (T - Tref)) * p / pref),
                )

            #   get final wing of the line according to Gamma0, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW * Gamma0, OmegaWingHW * GammaD)

            #   shift coefficient
            # Shift0 = Shift0DB*p/pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB - OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB + OmegaWingF)
            lineshape_vals = PROFILE_VOIGT(
                LineCenterDB + Shift0,
                GammaD,
                Gamma0,
                Omegas[BoundIndexLower:BoundIndexUpper],
            )[0]
            Xsect[BoundIndexLower:BoundIndexUpper] += (
                factor
                / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * LineIntensity
                * lineshape_vals
            )

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def absorptionCoefficient_Lorentz(
    Components=None,
    SourceTables=None,
    partitionFunction=PYTIPS,
    Environment=None,
    OmegaRange=None,
    OmegaStep=None,
    OmegaWing=None,
    IntensityThreshold=DefaultIntensityThreshold,
    OmegaWingHW=DefaultOmegaWingHW,
    GammaL="gamma_air",
    HITRAN_units=True,
    LineShift=True,
    File=None,
    Format=None,
    OmegaGrid=None,
    WavenumberRange=None,
    WavenumberStep=None,
    WavenumberWing=None,
    WavenumberWingHW=None,
    WavenumberGrid=None,
    Diluent={},
    EnvDependences=None,
):
    """
    INPUT PARAMETERS: 
        Components:  list of tuples [(M,I,D)], where
                        M - HITRAN molecule number,
                        I - HITRAN isotopologue number,
                        D - relative abundance (optional)
        SourceTables:  list of tables from which to calculate cross-section   (optional)
        partitionFunction:  pointer to partition function (default is PYTIPS) (optional)
        Environment:  dictionary containing thermodynamic parameters.
                        'p' - pressure in atmospheres,
                        'T' - temperature in Kelvin
                        Default={'p':1.,'T':296.}
        WavenumberRange:  wavenumber range to consider.
        WavenumberStep:   wavenumber step to consider. 
        WavenumberWing:   absolute wing for calculating a lineshape (in cm-1) 
        WavenumberWingHW:  relative wing for calculating a lineshape (in halfwidths)
        IntensityThreshold:  threshold for intensities
        GammaL:  specifies broadening parameter ('gamma_air' or 'gamma_self')
        HITRAN_units:  use cm2/molecule (True) or cm-1 (False) for absorption coefficient
        File:   write output to file (if specified)
        Format:  c-format of file output (accounts for significant digits in WavenumberStep)
    OUTPUT PARAMETERS: 
        Wavenum: wavenumber grid with respect to parameters WavenumberRange and WavenumberStep
        Xsect: absorption coefficient calculated on the grid
    ---
    DESCRIPTION:
        Calculate absorption coefficient using Lorentz profile.
        Absorption coefficient is calculated at arbitrary temperature and pressure.
        User can vary a wide range of parameters to control a process of calculation.
        The choise of these parameters depends on properties of a particular linelist.
        Default values are a sort of guess which gives a decent precision (on average) 
        for a reasonable amount of cpu time. To increase calculation accuracy,
        user should use a trial and error method.
    ---
    EXAMPLE OF USAGE:
        nu,coef = absorptionCoefficient_Lorentz(((2,1),),'co2',WavenumberStep=0.01,
                                              HITRAN_units=False,GammaL='gamma_self')
    ---
    """

    # Paremeters OmegaRange,OmegaStep,OmegaWing,OmegaWingHW, and OmegaGrid
    # are deprecated and given for backward compatibility with the older versions.
    if WavenumberRange:
        OmegaRange = WavenumberRange
    if WavenumberStep:
        OmegaStep = WavenumberStep
    if WavenumberWing:
        OmegaWing = WavenumberWing
    if WavenumberWingHW:
        OmegaWingHW = WavenumberWingHW
    if WavenumberGrid:
        OmegaGrid = WavenumberGrid

    # "bug" with 1-element list
    Components = listOfTuples(Components)
    SourceTables = listOfTuples(SourceTables)

    # determine final input values
    (
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    ) = getDefaultValuesForXsect(
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    )

    # warn user about too large omega step
    if OmegaStep > 0.1:
        warn("Big wavenumber step: possible accuracy decline")

    # get uniform linespace for cross-section
    # number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    # Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        # Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.0)  # K
    pref = __FloatType__(1.0)  # atm

    # actual temperature and pressure
    T = Environment["T"]  # K
    p = Environment["p"]  # atm

    # create dictionary from Components
    ABUNDANCES = {}
    NATURAL_ABUNDANCES = {}
    for Component in Components:
        M = Component[0]
        I = Component[1]
        if len(Component) >= 3:
            ni = Component[2]
        else:
            try:
                ni = ISO[(M, I)][ISO_INDEX["abundance"]]
            except KeyError:
                raise Exception("cannot find component M,I = %d,%d." % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX["abundance"]]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # setup the default empty environment dependence function
    if not EnvDependences:

        def EnvDependences(ENV, LINE):
            return {}

    Env = Environment.copy()
    Env["Tref"] = Tref
    Env["pref"] = pref

    # setup the Diluent variable
    GammaL = GammaL.lower()
    if not Diluent:
        if GammaL == "gamma_air":
            Diluent = {"air": 1.0}
        elif GammaL == "gamma_self":
            Diluent = {"self": 1.0}
        else:
            raise Exception("Unknown GammaL value: %s" % GammaL)

    # Simple check
    for key in Diluent:
        val = Diluent[key]
        if val < 0 and val > 1:
            raise Exception("Diluent fraction must be in [0,1]")

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get the number of rows
        nline = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]

        # get parameter names for each table
        parnames = list(LOCAL_TABLE_CACHE[TableName]["data"].keys())

        # loop through line centers (single stream)
        for RowID in range(nline):

            # Get the custom environment dependences
            Line = {
                parname: LOCAL_TABLE_CACHE[TableName]["data"][parname][RowID]
                for parname in parnames
            }
            CustomEnvDependences = EnvDependences(Env, Line)

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]["data"]["nu"][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]["data"]["sw"][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]["data"]["elower"][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["molec_id"][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["local_iso_id"][RowID]

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES:
                continue

            # partition functions for T and Tref
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            if "sw" in CustomEnvDependences:
                LineIntensity = CustomEnvDependences["sw"]
            else:
                LineIntensity = EnvironmentDependency_Intensity(
                    LineIntensityDB,
                    T,
                    Tref,
                    SigmaT,
                    SigmaTref,
                    LowerStateEnergyDB,
                    LineCenterDB,
                )

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            if LineIntensity < IntensityThreshold:
                continue

            #   pressure broadening coefficient
            Gamma0 = 0.0
            Shift0 = 0.0
            for species in Diluent:
                species_lower = species.lower()

                abun = Diluent[species]

                gamma_name = "gamma_" + species_lower
                try:
                    Gamma0DB = LOCAL_TABLE_CACHE[TableName]["data"][gamma_name][RowID]
                except:
                    Gamma0DB = 0.0

                n_name = "n_" + species_lower
                try:
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][n_name][
                        RowID
                    ]
                    if species_lower == "self" and TempRatioPowerDB == 0.0:
                        # same for self as for air
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"][
                            "n_air"
                        ][RowID]
                except:
                    # TempRatioPowerDB = 0
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]["data"]["n_air"][
                        RowID
                    ]

                # Add to the final Gamma0
                Gamma0 += abun * CustomEnvDependences.get(
                    gamma_name,  # default ->
                    EnvironmentDependency_Gamma0(
                        Gamma0DB, T, Tref, p, pref, TempRatioPowerDB
                    ),
                )

                delta_name = "delta_" + species_lower
                try:
                    Shift0DB = LOCAL_TABLE_CACHE[TableName]["data"][delta_name][RowID]
                except:
                    Shift0DB = 0.0

                deltap_name = "deltap_" + species_lower
                try:
                    deltap = LOCAL_TABLE_CACHE[TableName]["data"][deltap_name][RowID]
                except:
                    deltap = 0.0

                Shift0 += abun * CustomEnvDependences.get(
                    delta_name,  # default ->
                    ((Shift0DB + deltap * (T - Tref)) * p / pref),
                )

            #   get final wing of the line according to Gamma0, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW * Gamma0)

            #   shift coefficient
            # Shift0 = Shift0DB*p/pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB - OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB + OmegaWingF)
            lineshape_vals = PROFILE_LORENTZ(
                LineCenterDB + Shift0, Gamma0, Omegas[BoundIndexLower:BoundIndexUpper]
            )
            Xsect[BoundIndexLower:BoundIndexUpper] += (
                factor
                / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * LineIntensity
                * lineshape_vals
            )

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


# Alias for a profile selector
absorptionCoefficient = absorptionCoefficient_HT

# ==========================================================================================
# =========================== /NEW ABSORPTION COEFFICIENT ===================================
# ==========================================================================================

# calculate apsorption for Doppler profile


def absorptionCoefficient_Doppler(
    Components=None,
    SourceTables=None,
    partitionFunction=PYTIPS,
    Environment=None,
    OmegaRange=None,
    OmegaStep=None,
    OmegaWing=None,
    IntensityThreshold=DefaultIntensityThreshold,
    OmegaWingHW=DefaultOmegaWingHW,
    ParameterBindings=DefaultParameterBindings,
    EnvironmentDependencyBindings=DefaultEnvironmentDependencyBindings,
    GammaL="dummy",
    HITRAN_units=True,
    LineShift=True,
    File=None,
    Format=None,
    OmegaGrid=None,
    WavenumberRange=None,
    WavenumberStep=None,
    WavenumberWing=None,
    WavenumberWingHW=None,
    WavenumberGrid=None,
):
    """
    INPUT PARAMETERS: 
        Components:  list of tuples [(M,I,D)], where
                        M - HITRAN molecule number,
                        I - HITRAN isotopologue number,
                        D - abundance (optional)
        SourceTables:  list of tables from which to calculate cross-section   (optional)
        partitionFunction:  pointer to partition function (default is PYTIPS) (optional)
        Environment:  dictionary containing thermodynamic parameters.
                        'p' - pressure in atmospheres,
                        'T' - temperature in Kelvin
                        Default={'p':1.,'T':296.}
        WavenumberRange:  wavenumber range to consider.
        WavenumberStep:   wavenumber step to consider. 
        WavenumberWing:   absolute wing for calculating a lineshape (in cm-1) 
        WavenumberWingHW:  relative wing for calculating a lineshape (in halfwidths)
        IntensityThreshold:  threshold for intensities
        GammaL:  specifies broadening parameter ('gamma_air' or 'gamma_self')
        HITRAN_units:  use cm2/molecule (True) or cm-1 (False) for absorption coefficient
        File:   write output to file (if specified)
        Format:  c-format of file output (accounts for significant digits in WavenumberStep)
    OUTPUT PARAMETERS: 
        Wavenum: wavenumber grid with respect to parameters OmegaRange and OmegaStep
        Xsect: absorption coefficient calculated on the grid
    ---
    DESCRIPTION:
        Calculate absorption coefficient using Doppler (Gauss) profile.
        Absorption coefficient is calculated at arbitrary temperature and pressure.
        User can vary a wide range of parameters to control a process of calculation.
        The choise of these parameters depends on properties of a particular linelist.
        Default values are a sort of guess which give a decent precision (on average) 
        for a reasonable amount of cpu time. To increase calculation accuracy,
        user should use a trial and error method.
    ---
    EXAMPLE OF USAGE:
        nu,coef = absorptionCoefficient_Doppler(((2,1),),'co2',WavenumberStep=0.01,
                                              HITRAN_units=False,GammaL='gamma_self')
    ---
    """

    if WavenumberRange:
        OmegaRange = WavenumberRange
    if WavenumberStep:
        OmegaStep = WavenumberStep
    if WavenumberWing:
        OmegaWing = WavenumberWing
    if WavenumberWingHW:
        OmegaWingHW = WavenumberWingHW
    if WavenumberGrid:
        OmegaGrid = WavenumberGrid

    # "bug" with 1-element list
    Components = listOfTuples(Components)
    SourceTables = listOfTuples(SourceTables)

    # determine final input values
    (
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    ) = getDefaultValuesForXsect(
        Components,
        SourceTables,
        Environment,
        OmegaRange,
        OmegaStep,
        OmegaWing,
        IntensityThreshold,
        Format,
    )
    # special for Doppler case: set OmegaStep to a smaller value
    if not OmegaStep:
        OmegaStep = 0.001

    # warn user about too large omega step
    if OmegaStep > 0.005:
        warn("Big wavenumber step: possible accuracy decline")

    # get uniform linespace for cross-section
    # number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    # Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        # Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.0)  # K
    pref = __FloatType__(1.0)  # atm

    # actual temperature and pressure
    T = Environment["T"]  # K
    p = Environment["p"]  # atm

    # create dictionary from Components
    ABUNDANCES = {}
    NATURAL_ABUNDANCES = {}
    for Component in Components:
        M = Component[0]
        I = Component[1]
        if len(Component) >= 3:
            ni = Component[2]
        else:
            try:
                ni = ISO[(M, I)][ISO_INDEX["abundance"]]
            except KeyError:
                raise Exception("cannot find component M,I = %d,%d." % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX["abundance"]]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get line centers
        nline = LOCAL_TABLE_CACHE[TableName]["header"]["number_of_rows"]

        # loop through line centers (single stream)
        for RowID in range(nline):

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]["data"]["nu"][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]["data"]["sw"][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]["data"]["elower"][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["molec_id"][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]["data"]["local_iso_id"][RowID]
            if LineShift:
                Shift0DB = LOCAL_TABLE_CACHE[TableName]["data"]["delta_air"][RowID]
            else:
                Shift0DB = 0

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES:
                continue

            # partition functions for T and Tref
            # TODO: optimize
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            LineIntensity = EnvironmentDependency_Intensity(
                LineIntensityDB,
                T,
                Tref,
                SigmaT,
                SigmaTref,
                LowerStateEnergyDB,
                LineCenterDB,
            )

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            # TODO: apply wing narrowing instead of filtering, this would be more appropriate
            if LineIntensity < IntensityThreshold:
                continue

            #   doppler broadening coefficient (GammaD)
            # GammaDDB = cSqrtLn2*LineCenterDB/cc*sqrt(2*cBolts*T/molecularMass(MoleculeNumberDB,IsoNumberDB))
            # GammaD = EnvironmentDependency_GammaD(GammaDDB,T,Tref)
            # print(GammaD)

            cMassMol = 1.66053873e-27
            # cSqrt2Ln2 = 1.1774100225
            fSqrtMass = sqrt(molecularMass(MoleculeNumberDB, IsoNumberDB))
            # fSqrtMass = sqrt(32831.2508809)
            cc_ = 2.99792458e8
            cBolts_ = 1.3806503e-23
            # cBolts_ = 1.3806488E-23
            GammaD = (
                (cSqrt2Ln2 / cc_)
                * sqrt(cBolts_ / cMassMol)
                * sqrt(T)
                * LineCenterDB
                / fSqrtMass
            )

            # GammaD = 4.30140e-7*LineCenterDB*sqrt(T/molecularMass(MoleculeNumberDB,IsoNumberDB))

            # cc_ = 2.99792458e8 # 2.99792458e10 # 2.99792458e8
            # cBolts_ = 1.3806503e-23 #1.3806488E-16 # 1.380648813E-16 # 1.3806503e-23 # 1.3806488E-23
            # GammaD = sqrt(log(2))*LineCenterDB*sqrt(2*cBolts_*T/(cMassMol*molecularMass(MoleculeNumberDB,IsoNumberDB)*cc_**2))
            # print(GammaD)

            #   get final wing of the line according to GammaD, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW * GammaD)

            #   shift coefficient
            Shift0 = Shift0DB * p / pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB - OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB + OmegaWingF)
            lineshape_vals = PROFILE_DOPPLER(
                LineCenterDB + Shift0, GammaD, Omegas[BoundIndexLower:BoundIndexUpper]
            )
            # lineshape_vals = PROFILE_VOIGT(LineCenterDB,GammaD,cZero,Omegas[BoundIndexLower:BoundIndexUpper])[0]
            # Xsect[BoundIndexLower:BoundIndexUpper] += lineshape_vals # DEBUG
            Xsect[BoundIndexLower:BoundIndexUpper] += (
                factor
                / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]
                * LineIntensity
                * lineshape_vals
            )

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


# ---------------------------------------------------------------------------
# SHORTCUTS AND ALIASES FOR ABSORPTION COEFFICIENTS
# ---------------------------------------------------------------------------

absorptionCoefficient_Gauss = absorptionCoefficient_Doppler


def abscoef_HT(table=None, step=None, grid=None, env={"T": 296.0, "p": 1.0}, file=None):
    return absorptionCoefficient_HT(
        SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file
    )


def abscoef_Voigt(
    table=None, step=None, grid=None, env={"T": 296.0, "p": 1.0}, file=None
):
    return absorptionCoefficient_Voigt(
        SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file
    )


def abscoef_Lorentz(
    table=None, step=None, grid=None, env={"T": 296.0, "p": 1.0}, file=None
):
    return absorptionCoefficient_Lorentz(
        SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file
    )


def abscoef_Doppler(
    table=None, step=None, grid=None, env={"T": 296.0, "p": 1.0}, file=None
):
    return absorptionCoefficient_Doppler(
        SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file
    )


abscoef_Gauss = abscoef_Doppler


# default
def abscoef(table=None, step=None, grid=None, env={"T": 296.0, "p": 1.0}, file=None):
    return absorptionCoefficient_Lorentz(
        SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file
    )


# ---------------------------------------------------------------------------


def transmittanceSpectrum(
    Omegas,
    AbsorptionCoefficient,
    Environment={"l": 100.0},
    File=None,
    Format="%e %e",
    Wavenumber=None,
):
    """
    INPUT PARAMETERS: 
        Wavenumber/Omegas:   wavenumber grid                    (required)
        AbsorptionCoefficient:  absorption coefficient on grid  (required)
        Environment:  dictionary containing path length in cm.
                      Default={'l':100.}
        File:         name of the output file                 (optional) 
        Format: c format used in file output, default '%e %e' (optional)
    OUTPUT PARAMETERS: 
        Wavenum: wavenumber grid
        Xsect:  transmittance spectrum calculated on the grid
    ---
    DESCRIPTION:
        Calculate a transmittance spectrum (dimensionless) based
        on previously calculated absorption coefficient.
        Transmittance spectrum is calculated at an arbitrary
        optical path length 'l' (1 m by default)
    ---
    EXAMPLE OF USAGE:
        nu,trans = transmittanceSpectrum(nu,coef)
    ---
    """
    # compatibility with older versions
    if Wavenumber:
        Omegas = Wavenumber
    l = Environment["l"]
    Xsect = exp(-AbsorptionCoefficient * l)
    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def absorptionSpectrum(
    Omegas,
    AbsorptionCoefficient,
    Environment={"l": 100.0},
    File=None,
    Format="%e %e",
    Wavenumber=None,
):
    """
    INPUT PARAMETERS: 
        Wavenumber/Omegas:   wavenumber grid                    (required)
        AbsorptionCoefficient:  absorption coefficient on grid  (required)
        Environment:  dictionary containing path length in cm.
                      Default={'l':100.}
        File:         name of the output file                 (optional) 
        Format: c format used in file output, default '%e %e' (optional)
    OUTPUT PARAMETERS: 
        Wavenum: wavenumber grid
        Xsect:  transmittance spectrum calculated on the grid
    ---
    DESCRIPTION:
        Calculate an absorption spectrum (dimensionless) based
        on previously calculated absorption coefficient.
        Absorption spectrum is calculated at an arbitrary
        optical path length 'l' (1 m by default)
    ---
    EXAMPLE OF USAGE:
        nu,absorp = absorptionSpectrum(nu,coef)
    ---
    """
    # compatibility with older versions
    if Wavenumber:
        Omegas = Wavenumber
    l = Environment["l"]
    Xsect = 1 - exp(-AbsorptionCoefficient * l)
    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def radianceSpectrum(
    Omegas,
    AbsorptionCoefficient,
    Environment={"l": 100.0, "T": 296.0},
    File=None,
    Format="%e %e",
    Wavenumber=None,
):
    """
    INPUT PARAMETERS: 
        Wavenumber/Omegas:   wavenumber grid                   (required)
        AbsorptionCoefficient:  absorption coefficient on grid (required)
        Environment:  dictionary containing path length in cm.
                      and temperature in Kelvin.
                      Default={'l':100.,'T':296.}
        File:         name of the output file                 (optional) 
        Format: c format used in file output, default '%e %e' (optional)
    OUTPUT PARAMETERS: 
        Wavenum: wavenumber grid
        Xsect:  radiance spectrum calculated on the grid
    ---
    DESCRIPTION:
        Calculate a radiance spectrum (in W/sr/cm^2/cm-1) based
        on previously calculated absorption coefficient.
        Radiance spectrum is calculated at an arbitrary
        optical path length 'l' (1 m by default) and 
        temperature 'T' (296 K by default). For obtaining a
        physically meaningful result 'T' must be the same 
        as a temperature which was used in absorption coefficient.
    ---
    EXAMPLE OF USAGE:
        nu,radi = radianceSpectrum(nu,coef)
    ---
    """
    # compatibility with older versions
    if Wavenumber:
        Omegas = Wavenumber
    l = Environment["l"]
    T = Environment["T"]
    Alw = 1 - exp(-AbsorptionCoefficient * l)
    LBBTw = (
        2
        * hh
        * cc ** 2
        * Omegas ** 3
        / (exp(hh * cc * Omegas / (cBolts * T)) - 1)
        * 1.0e-7
    )
    Xsect = Alw * LBBTw  # W/sr/cm**2/cm**-1
    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


# GET X,Y FOR FINE PLOTTING OF A STICK SPECTRUM
def getStickXY(TableName):
    """
    Get X and Y for fine plotting of a stick spectrum.
    Usage: X,Y = getStickXY(TableName).
    """
    cent, intens = getColumns(TableName, ("nu", "sw"))
    n = len(cent)
    cent_ = zeros(n * 3)
    intens_ = zeros(n * 3)
    for i in range(n):
        intens_[3 * i] = 0
        intens_[3 * i + 1] = intens[i]
        intens_[3 * i + 2] = 0
        cent_[(3 * i) : (3 * i + 3)] = cent[i]
    return cent_, intens_


# /GET X,Y FOR FINE PLOTTING OF A STICK SPECTRUM


# LOW-RES SPECTRA (CONVOLUTION WITH APPARATUS FUNCTION)

# /LOW-RES SPECTRA (CONVOLUTION WITH APPARATUS FUNCTION)

# /----------------------------------------------------------------------------


# ------------------  HITRAN-ON-THE-WEB COMPATIBILITY -------------------------


def read_hotw(filename):
    """
    Read cross-section file fetched from HITRAN-on-the-Web.
    The format of the file line must be as follows: 
      nu, coef
    Other lines are omitted.
    """
    import sys

    f = open(filename, "r")
    nu = []
    coef = []
    for line in f:
        pars = line.split()
        try:
            nu.append(float(pars[0]))
            coef.append(float(pars[1]))
        except:
            if False:
                print(sys.exc_info())
            else:
                pass
    return array(nu), array(coef)


# alias for read_hotw for backwards compatibility
read_xsect = read_hotw

# /----------------------------------------------------------------------------

# ------------------  SPECTRAL CONVOLUTION -------------------------

# rectangular slit function


def SLIT_RECTANGULAR(x, g):
    """
    Instrumental (slit) function.
    B(x) = 1/γ , if |x| ≤ γ/2 & B(x) = 0, if |x| > γ/2,
    where γ is a slit width or the instrumental resolution.
    """
    index_inner = abs(x) <= g / 2
    index_outer = ~index_inner
    y = zeros(len(x))
    y[index_inner] = 1 / g
    y[index_outer] = 0
    return y


# triangular slit function


def SLIT_TRIANGULAR(x, g):
    """
    Instrumental (slit) function.
    B(x) = 1/γ*(1-|x|/γ), if |x| ≤ γ & B(x) = 0, if |x| > γ,
    where γ is the line width equal to the half base of the triangle.
    """
    index_inner = abs(x) <= g
    index_outer = ~index_inner
    y = zeros(len(x))
    y[index_inner] = 1 / g * (1 - abs(x[index_inner]) / g)
    y[index_outer] = 0
    return y


# gaussian slit function


def SLIT_GAUSSIAN(x, g):
    """
    Instrumental (slit) function.
    B(x) = sqrt(ln(2)/pi)/γ*exp(-ln(2)*(x/γ)**2),
    where γ/2 is a gaussian half-width at half-maximum.
    """
    g /= 2
    return sqrt(log(2)) / (sqrt(pi) * g) * exp(-log(2) * (x / g) ** 2)


# dispersion slit function


def SLIT_DISPERSION(x, g):
    """
    Instrumental (slit) function.
    B(x) = γ/pi/(x**2+γ**2),
    where γ/2 is a lorentzian half-width at half-maximum.
    """
    g /= 2
    return g / pi / (x ** 2 + g ** 2)


# cosinus slit function


def SLIT_COSINUS(x, g):
    return (cos(pi / g * x) + 1) / (2 * g)


# diffraction slit function


def SLIT_DIFFRACTION(x, g):
    """
    Instrumental (slit) function.
    """
    y = zeros(len(x))
    index_zero = x == 0
    index_nonzero = ~index_zero
    dk_ = pi / g
    x_ = dk_ * x[index_nonzero]
    w_ = sin(x_)
    r_ = w_ ** 2 / x_ ** 2
    y[index_zero] = 1
    y[index_nonzero] = r_ / g
    return y


# apparatus function of the ideal Michelson interferometer


def SLIT_MICHELSON(x, g):
    """
    Instrumental (slit) function.
    B(x) = 2/γ*sin(2pi*x/γ)/(2pi*x/γ) if x!=0 else 1,
    where 1/γ is the maximum optical path difference.
    """
    y = zeros(len(x))
    index_zero = x == 0
    index_nonzero = ~index_zero
    dk_ = 2 * pi / g
    x_ = dk_ * x[index_nonzero]
    y[index_zero] = 1
    y[index_nonzero] = 2 / g * sin(x_) / x_
    return y


# spectral convolution with an apparatus (slit) function


def convolveSpectrum(
    Omega,
    CrossSection,
    Resolution=0.1,
    AF_wing=10.0,
    SlitFunction=SLIT_RECTANGULAR,
    Wavenumber=None,
):
    """
    INPUT PARAMETERS: 
        Wavenumber/Omega:    wavenumber grid                     (required)
        CrossSection:  high-res cross section calculated on grid (required)
        Resolution:    instrumental resolution γ                 (optional)
        AF_wing:       instrumental function wing                (optional)
        SlitFunction:  instrumental function for low-res spectra calculation (optional)
    OUTPUT PARAMETERS: 
        Wavenum: wavenumber grid
        CrossSection: low-res cross section calculated on grid
        i1: lower index in Omega input
        i2: higher index in Omega input
        slit: slit function calculated over grid [-AF_wing; AF_wing]
                with the step equal to instrumental resolution. 
    ---
    DESCRIPTION:
        Produce a simulation of experimental spectrum via the convolution 
        of a “dry” spectrum with an instrumental function.
        Instrumental function is provided as a parameter and
        is calculated in a grid with the width=AF_wing and step=Resolution.
    ---
    EXAMPLE OF USAGE:
        nu_,radi_,i,j,slit = convolveSpectrum(nu,radi,Resolution=2.0,AF_wing=10.0,
                                                SlitFunction=SLIT_MICHELSON)
    ---
    """
    # compatibility with older versions
    if Wavenumber:
        Omega = Wavenumber
    step = Omega[1] - Omega[0]
    if step >= Resolution:
        raise Exception("step must be less than resolution")
    # x = arange(-AF_wing,AF_wing+step,step)
    x = arange_(-AF_wing, AF_wing + step, step)  # fix
    slit = SlitFunction(x, Resolution)
    # FIXING THE BUG: normalize slit function
    slit /= sum(slit) * step  # simple normalization
    left_bnd = len(slit) / 2
    right_bnd = len(Omega) - len(slit) / 2
    # CrossSectionLowRes = convolve(CrossSection,slit,mode='valid')*step
    CrossSectionLowRes = convolve(CrossSection, slit, mode="same") * step
    # return Omega[left_bnd:right_bnd],CrossSectionLowRes,left_bnd,right_bnd,slit
    return (
        Omega[left_bnd:right_bnd],
        CrossSectionLowRes[left_bnd:right_bnd],
        left_bnd,
        right_bnd,
        slit,
    )


# DEBUG
# spectral convolution with an apparatus (slit) function


def convolveSpectrumSame(
    Omega, CrossSection, Resolution=0.1, AF_wing=10.0, SlitFunction=SLIT_RECTANGULAR
):
    """
    Convolves cross section with a slit function with given parameters.
    """
    step = Omega[1] - Omega[0]
    x = arange(-AF_wing, AF_wing + step, step)
    slit = SlitFunction(x, Resolution)
    print("step=")
    print(step)
    print("x=")
    print(x)
    print("slitfunc=")
    print(SlitFunction)
    CrossSectionLowRes = convolve(CrossSection, slit, mode="same") * step
    return Omega, CrossSectionLowRes, None, None, slit


# DEBUG


def convolveSpectrumFull(
    Omega, CrossSection, Resolution=0.1, AF_wing=10.0, SlitFunction=SLIT_RECTANGULAR
):
    """
    Convolves cross section with a slit function with given parameters.
    """
    step = Omega[1] - Omega[0]
    x = arange(-AF_wing, AF_wing + step, step)
    slit = SlitFunction(x, Resolution)
    print("step=")
    print(step)
    print("x=")
    print(x)
    print("slitfunc=")
    print(SlitFunction)
    CrossSectionLowRes = convolve(CrossSection, slit, mode="full") * step
    return Omega, CrossSectionLowRes, None, None


# ------------------------------------------------------------------
