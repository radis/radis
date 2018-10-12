# -*- coding: utf-8 -*-

'''
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
'''

from __future__ import print_function, absolute_import, division, unicode_literals

#import httplib
#import urllib2
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
#from collections import OrderedDict
from warnings import warn, simplefilter
import pydoc
import six
from six.moves import range
from six.moves import zip

# Enable warning repetitions
simplefilter('always', UserWarning)

# Python 3 compatibility
try:
    import urllib.request as urllib2
except ImportError:
    import six.moves.urllib.request, six.moves.urllib.error, six.moves.urllib.parse

HAPI_VERSION = '1.1.0.6'
# CHANGES:
# FIXED GRID BUG (ver. 1.1.0.1)
# FIXED OUTPUT FORMAT FOR CROSS-SECTIONS (ver. 1.1.0.1)
# ADDED CPF BY SCHREIER (JQSRT_112_2011) (ver. 1.1.0.2)
# OPTIMIZED EXPRESSION EVALUATIONS FOR SELECT (ver. 1.1.0.3)
# ADDED SUPPORT FOR MIXTURES (ver. 1.1.0.4)
# ADDED SUPPORT FOR USER-DEFINED ENV DEPENDENCES (ver. 1.1.0.5)
# ADDED PROFILE SELECTION (ALPHA) (ver. 1.1.0.6)

# version header
print('HAPI version: %s' % HAPI_VERSION)
print('To get the most up-to-date version please check http://hitran.org/hapi')

# define precision
__ComplexType__ = complex128
__IntegerType__ = int64
__FloatType__ = float64

# define zero
cZero = __FloatType__(0.)

# physical constants
cBolts = 1.380648813E-16  # erg/K, CGS
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
    warn('GLOBAL_DEBUG is set to True!')

GLOBAL_CURRENT_DIR = '.'

GLOBAL_HITRAN_APIKEY = 'e20e4bd3-e12c-4931-99e0-4c06e88536bd'

GLOBAL_USER = 'user'
GLOBAL_REQUISITES = []

GLOBAL_CONNECTION = []
GLOBAL_DATABASE = 'hitran'

LOCAL_HOST = 'http://localhost'

# DEBUG switch
if GLOBAL_DEBUG:
    GLOBAL_HOST = LOCAL_HOST+':8000'  # localhost
else:
    GLOBAL_HOST = 'http://hitran.org'

# this is a backup url in the case GLOBAL_HOST does not work
GLOBAL_HOST_BACKUP = 'http://hitranazure.cloudapp.net/'

# In this "robust" version of arange the grid doesn't suffer
# from the shift of the nodes due to error accumulation.
# This effect is pronounced only if the step is sufficiently small.


def arange_(lower, upper, step):
    npnt = floor((upper-lower)/step)+1
    upper_new = lower + step*(npnt-1)
    if abs((upper-upper_new)-step) < 1e-10:
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
GlobalQueryString = ''

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
        raise Exception('can''t setup connection')

# interface for HTTP-get method
# Connection must be established before use


def httpGet(URL, Connection=GLOBAL_CONNECTION):
    Method = 'get'
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
    Host = getAttribute(Connection, 'host')
    HostQuery = parseToFrontend(Query)
    URL = Host+HostQuery
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


NODE_NAME = 'local'

GLOBAL_NODENAMES = {
    0: 'hitran-main',
    1: 'local'
}

GLOBAL_NODELIST = {
    0: {  # main HITRAN node
        'host': GLOBAL_HOST,
        'ACCESS_KEY': '9b6a7975-2a84-43d8-920e-f4dea9db6805'  # guest
    },
    1: {  # local node prototype
        'host': LOCAL_HOST,
        'ACCESS_KEY': '6cfd7040-24a6-4197-81f9-6e25e50005b2',  # admin
    }
}


def createNode(NodeID, NodeList=GLOBAL_NODELIST):
    # create a node, throw if exists
    node = NodeList.get(NodeID)
    if node:
        raise Exception('node %s already exists' % NodeName)
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
            raise Exception('node %s doesn''t have property %s' %
                            (ModeName, Propname))
    else:
        raise Exception('no such node %s' % Nodename)


def setNodeProperty(NodeID, PropName, PropValue, NodeList=GLOBAL_NODELIST):
    # set a property for certain node
    # throw exception if node not found
    # if the property doesn't exist it will appear
    node = NodeList.get(NodeID)
    if not node:
        raise Exception('no such node %s ' % NodeName)
    NodeList[PropName] = PropValue
    return


def resolveNodeID(NodeName, NodeNames=GLOBAL_NODENAMES):
    for NodeID in NodeNames.keys():
        if NodeNames[NodeID] == NodeName:
            return NodeID


def checkAccess(DBName, TableName, NodeName, UserName, Requisites, NodeList=GLOBAL_NODELIST, NodeNames=GLOBAL_NODENAMES):
    # simple node-level authentication (bridge to AUTH system)
    NodeID = resolveNodeID(NodeName, NodeNames)
    Node = NodeList[NodeID]
    if Requisites.key in Node['keys_allowed']:
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
    'admin': {
        'ACCESS_KEY': '6cfd7040-24a6-4197-81f9-6e25e50005b2',
        'LEVEL': 'ADMIN'
    },
    'guest': {
        'ACCESS_KEY': '9b6a7975-2a84-43d8-920e-f4dea9db6805',
        'LEVEL': 'USER'
    }
}


def addUser():
    pass


def deleteUser():
    pass


def authenticate(UserName, Requisites, Privileges=GLOBAL_PRIVILEGES):
    # Authentication
    key_list = [Privileges[User]['ACCESS_KEY'] for User in Privileges.keys]
    return True if Requisites.AccessKey in key_list else False


def checkPrivileges(Path, UserName=GLOBAL_USER, Requisites=GLOBAL_REQUISITES,
                    Privileges=GLOBAL_PRIVILEGES, NodeList=GLOBAL_NODELIST, Nodenames=GLOBAL_NODENAMES):
    # Privileges are checked before executing every query (needs optimization)
    # Path example: SOME_DB::SOME_TABLE::SOME_NODE
    if not authenticate(UserName, Requisites, Privileges):
        return False
    (DBName, TableName, NodeName) = Path.split('::')
    # loop on all nodes , use NODE_MANAGER's functions instead of
    #   working with GLOBAL_NODELIST directly
    if not checkAccess(DBName, TableName, NodeName, UserName, Requisites, NodeList, NodeNames):
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

BACKEND_DATABASE_NAME_DEFAULT = '.'

VARIABLES = {}
VARIABLES['BACKEND_DATABASE_NAME'] = BACKEND_DATABASE_NAME_DEFAULT

# For this node local DB is schema-dependent!
LOCAL_TABLE_CACHE = {
    'sampletab': {  # table
        'header': {  # header
            'order': ('column1', 'column2', 'column3'),
            'format': {
                'column1': '%10d',
                'column2': '%20f',
                'column3': '%30s'
            },
            'default': {
                'column1': 0,
                'column2': 0.0,
                'column3': ''
            },
            'number_of_rows': 3,
            'size_in_bytes': None,
            'table_name': 'sampletab',
            'table_type': 'strict'
        },  # /header
        'data': {
            'column1': [1, 2, 3],
            'column2': [10.5, 11.5, 12.5],
            'column3': ['one', 'two', 'three']
        },  # /data
    }  # /table
}  # hash-map of tables

# FORMAT CONVERSION LAYER

# converts between TRANSPORT_FORMAT and OBJECT_FORMAT
HITRAN_FORMAT_160 = {
    'M': {'pos':   1,   'len':  2,   'format': '%2d'},
    'I': {'pos':   3,   'len':  1,   'format': '%1d'},
    'nu': {'pos':   4,   'len': 12,   'format': '%12f'},
    'S': {'pos':  16,   'len': 10,   'format': '%10f'},
    'R': {'pos':  26,   'len':  0,   'format': '%0f'},
    'A': {'pos':  26,   'len': 10,   'format': '%10f'},
    'gamma_air': {'pos':  36,   'len':  5,   'format': '%5f'},
    'gamma_self': {'pos':  41,   'len':  5,   'format': '%5f'},
    'E_': {'pos':  46,   'len': 10,   'format': '%10f'},
    'n_air': {'pos':  56,   'len':  4,   'format': '%4f'},
    'delta_air': {'pos':  60,   'len':  8,   'format': '%8f'},
    'V': {'pos':  68,   'len': 15,   'format': '%15s'},
    'V_': {'pos':  83,   'len': 15,   'format': '%15s'},
    'Q': {'pos':  98,   'len': 15,   'format': '%15s'},
    'Q_': {'pos': 113,   'len': 15,   'format': '%15s'},
    'Ierr': {'pos': 128,   'len':  6,   'format': '%6s'},
    'Iref': {'pos': 134,   'len': 12,   'format': '%12s'},
    'flag': {'pos': 146,   'len':  1,   'format': '%1s'},
    'g': {'pos': 147,   'len':  7,   'format': '%7f'},
    'g_': {'pos': 154,   'len':  7,   'format': '%7f'}
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
        "gpp"
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
        "global_lower_quanta": "%15s"
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
        "global_lower_quanta": "000"
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
        "global_lower_quanta": "Electronic and vibrational quantum numbers and labels for the lower state of a transition"
    },
}

PARAMETER_META = \
    {
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
            "has_error": 0
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
            "has_error": 0
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
            "has_error": 0
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 0
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
        },
        "elower": {
            "id": 11,
            "name": "elower",
            "name_html": "<em>E\"</em>",
            "table_name": "",
            "description": "Lower-state energy",
            "description_html": "Lower-state energy",
            "default_fmt": "%10.4f",
            "default_units": "cm-1",
            "data_type": "float",
            "selectable": 1,
            "has_reference": 0,
            "has_error": 0
        },
        "gp": {
            "id": 12,
            "name": "gp",
            "name_html": "<em>g</em>\'",
            "table_name": "",
            "description": "Upper state degeneracy",
            "description_html": "Upper state degeneracy",
            "default_fmt": "%5d",
            "default_units": "",
            "data_type": "int",
            "selectable": 1,
            "has_reference": 0,
            "has_error": 0
        },
        "gpp": {
            "id": 13,
            "name": "gpp",
            "name_html": "<em>g</em>\"",
            "table_name": "",
            "description": "Lower state degeneracy",
            "description_html": "Lower state degeneracy",
            "default_fmt": "%5d",
            "default_units": "",
            "data_type": "int",
            "selectable": 1,
            "has_reference": 0,
            "has_error": 0
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
            "has_error": 0
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
            "has_error": 0
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
            "has_error": 0
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
            "has_error": 0
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
            "has_error": 0
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
            "has_error": 0
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
            "has_error": 0
        },
        "deltap_air": {
            "id": 21,
            "name": "deltap_air",
            "name_html": "<em>&delta;\'</em><sub>air</sub>",
            "table_name": "prm_deltap_air",
            "description": "Linear temperature dependence coefficient for air-induced pressure shift",
            "description_html": "Linear temperature dependence coefficient for air-induced pressure shift",
            "default_fmt": "%10.3e",
            "default_units": "",
            "data_type": "float",
            "selectable": 1,
            "has_reference": 1,
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
        },
        "deltap_self": {
            "id": 24,
            "name": "deltap_self",
            "name_html": "<em>&delta;\'</em><sub>self</sub>",
            "table_name": "prm_deltap_self",
            "description": "Linear temperature dependence coefficient for self-induced pressure shift",
            "description_html": "Linear temperature dependence coefficient for self-induced pressure shift",
            "default_fmt": "%10.3e",
            "default_units": "",
            "data_type": "float",
            "selectable": 1,
            "has_reference": 1,
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
        },
        "statep": {
            "id": 33,
            "name": "statep",
            "name_html": "qns\'",
            "table_name": "",
            "description": "Upper state quantum numbers",
            "description_html": "Upper state quantum numbers",
            "default_fmt": "%256s",
            "default_units": "",
            "data_type": "str",
            "selectable": 1,
            "has_reference": 0,
            "has_error": 0
        },
        "statepp": {
            "id": 34,
            "name": "statepp",
            "name_html": "qns\"",
            "table_name": "",
            "description": "Lower state quantum numbers",
            "description_html": "Lower state quantum numbers",
            "default_fmt": "%256s",
            "default_units": "",
            "data_type": "str",
            "selectable": 1,
            "has_reference": 0,
            "has_error": 0
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
            "has_error": 1
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
            "has_error": 0
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
            "has_error": 0
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
        },
        "deltap_h2": {
            "id": 41,
            "name": "deltap_H2",
            "name_html": "<em>&delta;\'</em><sub>H2</sub>",
            "table_name": "prm_deltap_H2",
            "description": "Linear temperature dependence coefficient for H2-induced pressure shift",
            "description_html": "Linear temperature dependence coefficient for H<sub>2</sub>-induced pressure shift",
            "default_fmt": "%10.3e",
            "default_units": "",
            "data_type": "float",
            "selectable": 1,
            "has_reference": 1,
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
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
            "has_error": 1
        },
        "gamma_HT_0_self_50": {
            "default_fmt": "%6.4f",
        },
        "n_HT_self_50": {
            "default_fmt": "%9.6f",
        },
        "gamma_HT_2_self_50": {
            "default_fmt": "%6.4f",
        },
        "delta_HT_0_self_50": {
            "default_fmt": "%9.6f",
        },
        "deltap_HT_self_50": {
            "default_fmt": "%9.6f",
        },
        "delta_HT_2_self_50": {
            "default_fmt": "%9.6f",
        },
        "gamma_HT_0_self_150": {
            "default_fmt": "%6.4f",
        },
        "n_HT_self_150": {
            "default_fmt": "%9.6f",
        },
        "gamma_HT_2_self_150": {
            "default_fmt": "%6.4f",
        },
        "delta_HT_0_self_150": {
            "default_fmt": "%9.6f",
        },
        "deltap_HT_self_150": {
            "default_fmt": "%9.6f",
        },
        "delta_HT_2_self_150": {
            "default_fmt": "%9.6f",
        },
        "gamma_HT_0_self_296": {
            "default_fmt": "%6.4f",
        },
        "n_HT_self_296": {
            "default_fmt": "%9.6f",
        },
        "gamma_HT_2_self_296": {
            "default_fmt": "%6.4f",
        },
        "delta_HT_0_self_296": {
            "default_fmt": "%9.6f",
        },
        "deltap_HT_self_296": {
            "default_fmt": "%9.6f",
        },
        "delta_HT_2_self_296": {
            "default_fmt": "%9.6f",
        },
        "gamma_HT_0_self_700": {
            "default_fmt": "%6.4f",
        },
        "n_HT_self_700": {
            "default_fmt": "%9.6f",
        },
        "gamma_HT_2_self_700": {
            "default_fmt": "%6.4f",
        },
        "delta_HT_0_self_700": {
            "default_fmt": "%9.6f",
        },
        "deltap_HT_self_700": {
            "default_fmt": "%9.6f",
        },
        "delta_HT_2_self_700": {
            "default_fmt": "%9.6f",
        },
        "nu_HT_self": {
            "default_fmt": "%6.4f",
        },
        "kappa_HT_self": {
            "default_fmt": "%9.6f",
        },
        "eta_HT_self": {
            "default_fmt": "%9.6f",
        },
    }


def transport2object(TransportData):
    pass


def object2transport(ObjectData):
    pass


def getFullTableAndHeaderName(TableName):
    # print('TableName=',TableName)
    fullpath_data = VARIABLES['BACKEND_DATABASE_NAME'] + \
        '/' + TableName + '.data'
    if not os.path.isfile(fullpath_data):
        fullpath_data = VARIABLES['BACKEND_DATABASE_NAME'] + \
            '/' + TableName + '.par'
        if not os.path.isfile(fullpath_data) and TableName != 'sampletab':
            raise Exception('Lonely header \"%s\"' % fullpath_data)
    fullpath_header = VARIABLES['BACKEND_DATABASE_NAME'] + \
        '/' + TableName + '.header'
    return fullpath_data, fullpath_header


def getParameterFormat(ParameterName, TableName):
    return LOCAL_TABLE_CACHE[TableName]['header']['format']


def getTableHeader(TableName):
    return LOCAL_TABLE_CACHE[TableName]['header']

# RowObject = list of tuples like (name,value,format)


def addRowObject(RowObject, TableName):
    # add RowObject to TableObject in CACHE
    # check consistency first
    if [p[0] for p in RowObject] != LOCAL_TABLE_CACHE[TableName]['header']['order']:
        raise Exception('The row is not consistent with the table')
    for par_name, par_value, par_format in RowObject:
        LOCAL_TABLE_CACHE[TableName]['data'][par_name] += par_value
    pass


def getRowObject(RowID, TableName):
    # return RowObject from TableObject in CACHE
    RowObject = []
    for par_name in LOCAL_TABLE_CACHE[TableName]['header']['order']:
        par_value = LOCAL_TABLE_CACHE[TableName]['data'][par_name][RowID]
        par_format = LOCAL_TABLE_CACHE[TableName]['header']['format'][par_name]
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
        #LOCAL_TABLE_CACHE[TableName]['data'][par_name] += [par_value]
        LOCAL_TABLE_CACHE[TableName]['data'][par_name].append(par_value)


def setRowObject(RowID, RowObject, TableName):
    number_of_rows = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    if RowID >= 0 and RowID < number_of_rows:
        for par_name, par_value, par_format in RowObject:
            LOCAL_TABLE_CACHE[TableName]['data'][par_name][RowID] = par_value
    else:
        # !!! XXX ATTENTION: THIS IS A TEMPORARY INSERTION XXX !!!
        LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows'] += 1
        addRowObject(RowObject, TableName)


def getDefaultRowObject(TableName):
    # get a default RowObject from a table
    RowObject = []
    for par_name in LOCAL_TABLE_CACHE[TableName]['header']['order']:
        par_value = LOCAL_TABLE_CACHE[TableName]['header']['default'][par_name]
        par_format = LOCAL_TABLE_CACHE[TableName]['header']['format'][par_name]
        RowObject.append((par_name, par_value, par_format))
    return RowObject


def subsetOfRowObject(ParameterNames, RowObject):
    # return a subset of RowObject according to
    #RowObjectNew = []
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


#FORMAT_PYTHON_REGEX = '^\%([0-9]*)\.?([0-9]*)([dfs])$'
FORMAT_PYTHON_REGEX = '^\%(\d*)(\.(\d*))?([edfsEDFS])$'

# Fortran string formatting
#  based on a pythonic format string


def formatString(par_format, par_value, lang='FORTRAN'):
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
    if ty.lower() in set(['f', 'e']):
        lng = int(lng) if lng else 0
        lngpnt = int(lngpnt) if lngpnt else 0
        result = par_format % par_value
        res = result.strip()
        if lng == lngpnt+1:
            if res[0:1] == '0':
                result = '%%%ds' % lng % res[1:]
        if par_value < 0:
            if res[1:2] == '0':
                result = '%%%ds' % lng % (res[0:1]+res[2:])
    return result


def formatGetLength(fmt, lang='FORTRAN'):
    regex = FORMAT_PYTHON_REGEX


def putRowObjectToString(RowObject):
    # serialize RowObject to string
    # TODO: support different languages (C,Fortran)
    output_string = ''
    for par_name, par_value, par_format in RowObject:
        # Python formatting
        #output_string += par_format % par_value
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
    "global_lower_quanta": "V_"
}


def putTableHeaderToString(TableName):
    output_string = ''
    regex = FORMAT_PYTHON_REGEX
    for par_name in LOCAL_TABLE_CACHE[TableName]['header']['order']:
        par_format = LOCAL_TABLE_CACHE[TableName]['header']['format'][par_name]
        (lng, trail, lngpnt, ty) = re.search(regex, par_format).groups()
        fmt = '%%%ss' % lng
        try:
            par_name_short = PARAMETER_NICKNAMES[par_name]
        except:
            par_name_short = par_name
        #output_string += fmt % par_name
        output_string += (fmt % par_name_short)[:int(lng)]
    return output_string


def getRowObjectFromString(input_string, TableName):
    # restore RowObject from string, get formats and names in TableName
    # print 'getRowObjectFromString:'
    pos = 0
    RowObject = []
    # print 'Header: '+str(LOCAL_TABLE_CACHE[TableName]['header'])
    for par_name in LOCAL_TABLE_CACHE[TableName]['header']['order']:
        #print 'ITERATION\npos: '+str(pos) #
        #print 'par_name: '+par_name #
        par_format = LOCAL_TABLE_CACHE[TableName]['header']['format'][par_name]
        #print 'par_format: '+par_format #
        regex = '^\%([0-9]+)\.?[0-9]*([dfs])$'
        regex = FORMAT_PYTHON_REGEX
        #print 'par_name: '+par_name #
        (lng, trail, lngpnt, ty) = re.search(regex, par_format).groups()
        lng = int(lng)
        #print 'lng,ty:'+str((lng,ty)) #
        par_value = input_string[pos:(pos+lng)]
        #print 'par_value: '+par_value #
        if ty == 'd':  # integer value
            par_value = int(par_value)
        elif ty.lower() in set(['e', 'f']):  # float value
            par_value = float(par_value)
        elif ty == 's':  # string value
            # par_value = par_value.strip() # strip spaces and tabs
            pass  # don't strip string value
        else:
            print('err1')
            raise Exception('Format \"%s\" is unknown' % par_format)
        RowObject.append((par_name, par_value, par_format))
        pos += lng
    # Do the same but now for extra (comma-separated) parameters
    if 'extra' in set(LOCAL_TABLE_CACHE[TableName]['header']):
        csv_chunks = input_string.split(LOCAL_TABLE_CACHE[TableName]['header'].
                                        get('extra_separator', ','))
        # Disregard the first "column-fixed" container if it presents:
        if LOCAL_TABLE_CACHE[TableName]['header'].get('order', []):
            pos = 1
        else:
            pos = 0
        for par_name in LOCAL_TABLE_CACHE[TableName]['header']['extra']:
            par_format = LOCAL_TABLE_CACHE[TableName]['header']['extra_format'][par_name]
            regex = '^\%([0-9]+)\.?[0-9]*([dfs])$'
            regex = FORMAT_PYTHON_REGEX
            (lng, trail, lngpnt, ty) = re.search(regex, par_format).groups()
            lng = int(lng)
            par_value = csv_chunks[pos]
            if ty == 'd':  # integer value
                try:
                    par_value = int(par_value)
                except:
                    par_value = 0
            elif ty.lower() in set(['e', 'f']):  # float value
                try:
                    par_value = float(par_value)
                except:
                    par_value = 0.0
            elif ty == 's':  # string value
                # par_value = par_value.strip() # strip spaces and tabs
                pass  # don't strip string value
            else:
                print('err')
                raise Exception('Format \"%s\" is unknown' % par_format)
            RowObject.append((par_name, par_value, par_format))
            pos += 1
    return RowObject
    # LOCAL_TABLE_CACHE[TableName]['data'][par_name] += par_value # or append()?

# Conversion between OBJECT_FORMAT and STORAGE_FORMAT
# This will substitute putTableToStorage and getTableFromStorage


def cache2storage(TableName):
    # print 'cache2storage:'
    try:
        os.mkdir(VARIABLES['BACKEND_DATABASE_NAME'])
    except:
        pass
    fullpath_data, fullpath_header = getFullTableAndHeaderName(TableName)
    # print 'fullpath_data:'+fullpath_data
    # print 'fullpath_header'+fullpath_header
    # check if file exists and throw an exception
    #if isfile(fullpath_data): raise Exception('Table \"%s\" already exists',NewTableName)
    #if isfile(fullpath_header): raise Exception('SCHEMA IS BROKEN')
    OutfileData = open(fullpath_data, 'w')
    OutfileHeader = open(fullpath_header, 'w')
    # write table data
    line_count = 1
    line_number = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    for RowID in range(0, LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']):
        # print '%d line from %d' % (line_count,line_number)
        line_count += 1
        RowObject = getRowObject(RowID, TableName)
        # print 'RowObject:'+str(RowObject)
        raw_string = putRowObjectToString(RowObject)
        # print 'RowObject_string:'+raw_string
        OutfileData.write(raw_string+'\n')
    # write table header
    TableHeader = getTableHeader(TableName)
    OutfileHeader.write(json.dumps(TableHeader, indent=2))


def storage2cache(TableName):
    # print 'storage2cache:'
    # print('TableName',TableName)
    fullpath_data, fullpath_header = getFullTableAndHeaderName(TableName)
    InfileData = open(fullpath_data, 'r')
    InfileHeader = open(fullpath_header, 'r')
    # try:
    header_text = InfileHeader.read()
    try:
        Header = json.loads(header_text)
    except:
        print('HEADER:')
        print(header_text)
        raise Exception('Invalid header')
    # print 'Header:'+str(Header)
    LOCAL_TABLE_CACHE[TableName] = {}
    LOCAL_TABLE_CACHE[TableName]['header'] = Header
    LOCAL_TABLE_CACHE[TableName]['data'] = {}
    # Check if Header['order'] and Header['extra'] contain
    #  parameters with same names, raise exception if true.
    #intersct = set(Header['order']).intersection(set(Header.get('extra',[])))
    intersct = set(Header.get('order', [])).intersection(
        set(Header.get('extra', [])))
    if intersct:
        raise Exception('Parameters with the same names: {}'.format(intersct))
    # initialize empty data to avoid problems
    glob_order = []
    glob_format = {}
    glob_default = {}
    if "order" in list(LOCAL_TABLE_CACHE[TableName]['header'].keys()):
        glob_order += LOCAL_TABLE_CACHE[TableName]['header']['order']
        glob_format.update(LOCAL_TABLE_CACHE[TableName]['header']['format'])
        glob_default.update(LOCAL_TABLE_CACHE[TableName]['header']['default'])
        for par_name in LOCAL_TABLE_CACHE[TableName]['header']['order']:
            LOCAL_TABLE_CACHE[TableName]['data'][par_name] = []
    if "extra" in list(LOCAL_TABLE_CACHE[TableName]['header'].keys()):
        glob_order += LOCAL_TABLE_CACHE[TableName]['header']['extra']
        glob_format.update(
            LOCAL_TABLE_CACHE[TableName]['header']['extra_format'])
        for par_name in LOCAL_TABLE_CACHE[TableName]['header']['extra']:
            glob_default[par_name] = PARAMETER_META[par_name]['default_fmt']
            LOCAL_TABLE_CACHE[TableName]['data'][par_name] = []
    line_count = 0
    #line_number = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    for line in InfileData:
        # print '%d line from %d' % (line_count,line_number)
        #print 'line: '+line #
        try:
            RowObject = getRowObjectFromString(line, TableName)
            line_count += 1
        except:
            continue
        # print 'RowObject: '+str(RowObject)
        addRowObject(RowObject, TableName)
    # except:
    #    raise Exception('TABLE FETCHING ERROR')
    LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows'] = line_count
    # Delete all character-separated values, treat them as column-fixed.
    try:
        del LOCAL_TABLE_CACHE[TableName]['header']['extra']
        del LOCAL_TABLE_CACHE[TableName]['header']['extra_format']
        del LOCAL_TABLE_CACHE[TableName]['header']['extra_separator']
    except:
        pass
    # Update header.order/format with header.extra/format if exist.
    LOCAL_TABLE_CACHE[TableName]['header']['order'] = glob_order
    LOCAL_TABLE_CACHE[TableName]['header']['format'] = glob_format
    LOCAL_TABLE_CACHE[TableName]['header']['default'] = glob_default
    InfileData.close()
    InfileHeader.close()
    print('                     Lines parsed: %d' % line_count)
    pass

# / FORMAT CONVERSION LAYER


def getTableNamesFromStorage(StorageName):
    file_names = listdir(StorageName)
    table_names = []
    for file_name in file_names:
        # search all files with "header" extensions
        #matchObject = re.search('(\w+)\.header$',file_name)
        matchObject = re.search('(.+)\.header$', file_name)
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
            #fname,fext = re.search('(\w+)\.(\w+)',file_name).groups()
            fname, fext = re.search('(.+)\.(\w+)', file_name).groups()
        except:
            continue
        if fext == 'header':
            headers[fname] = True
    for file_name in file_names:
        # check if extension is 'par' and the header is absent
        try:
            #fname,fext = re.search('(\w+)\.(\w+)',file_name).groups()
            fname, fext = re.search('(.+)\.(\w+)', file_name).groups()
        except:
            continue
        if fext == 'par' and fname not in headers:
            parfiles_without_header.append(fname)
    return parfiles_without_header


def createHeader(TableName):
    fname = TableName+'.header'
    fp = open(VARIABLES['BACKEND_DATABASE_NAME']+'/'+fname, 'w')
    if os.path.isfile(TableName):
        raise Exception('File \"%s\" already exists!' % fname)
    fp.write(json.dumps(HITRAN_DEFAULT_HEADER, indent=2))
    fp.close()


def loadCache():
    # print 'loadCache:'
    print('Using '+VARIABLES['BACKEND_DATABASE_NAME']+'\n')
    LOCAL_TABLE_CACHE = {}  # ?????
    table_names = getTableNamesFromStorage(VARIABLES['BACKEND_DATABASE_NAME'])
    # print('table_names=',table_names)
    parfiles_without_header = scanForNewParfiles(
        VARIABLES['BACKEND_DATABASE_NAME'])
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
        VARIABLES['BACKEND_DATABASE_NAME'] = db
    else:
        VARIABLES['BACKEND_DATABASE_NAME'] = BACKEND_DATABASE_NAME_DEFAULT
    # print 'databaseBegin:'
    # print(os.path.isdir("/home/el"))
    # print(os.path.exists("/home/el/myfile.txt"))
    if not os.path.exists(VARIABLES['BACKEND_DATABASE_NAME']):
        os.mkdir(VARIABLES['BACKEND_DATABASE_NAME'])
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
SAMPLE_CONDITIONS = ('AND', ('SET', 'internal_iso_id', [
                     1, 2, 3, 4, 5, 6]), ('>=', 'nu', 0), ('<=', 'nu', 100))

# sample hitranonline protocol
# http://hitran.cloudapp.net/lbl/5?output_format_id=1&iso_ids_list=5&numin=0&numax=100&access=api&key=e20e4bd3-e12c-4931-99e0-4c06e88536bd

CONDITION_OPERATIONS = set(['AND', 'OR', 'NOT', 'RANGE', 'IN', '<', '>', '<=', '>=',
                            '==', '!=', 'LIKE', 'STR', '+', '-', '*', '/', 'MATCH', 'SEARCH', 'FINDALL'])

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
        if args[i-1] >= args[i]:
            return False
    return True


def operationMORE(args):
    # any number of args
    for i in range(1, len(args)):
        if args[i-1] <= args[i]:
            return False
    return True


def operationLESSOREQUAL(args):
    # any number of args
    for i in range(1, len(args)):
        if args[i-1] > args[i]:
            return False
    return True


def operationMOREOREQUAL(args):
    # any number of args
    for i in range(1, len(args)):
        if args[i-1] < args[i]:
            return False
    return True


def operationEQUAL(args):
    # any number of args
    for i in range(1, len(args)):
        if args[i] != args[i-1]:
            return False
    return True


def operationNOTEQUAL(arg1, arg2):
    return arg1 != arg2


def operationSUM(args):
    # any numbers of arguments
    if type(args[0]) in set([int, float]):
        result = 0
    elif type(args[0]) in set([str, six.text_type]):
        result = ''
    else:
        raise Exception('SUM error: unknown arg type')
    for arg in args:
        result += arg
    return result


def operationDIFF(arg1, arg2):
    return arg1-arg2


def operationMUL(args):
    # any numbers of arguments
    if type(args[0]) in set([int, float]):
        result = 1
    else:
        raise Exception('MUL error: unknown arg type')
    for arg in args:
        result *= arg
    return result


def operationDIV(arg1, arg2):
    return arg1/arg2


def operationSTR(arg):
    # transform arg to str
    if type(arg) != str:
        raise Exception('Type mismatch: STR')
    return arg


def operationSET(arg):
    # transform arg to list
    if type(arg) not in set([list, tuple, set]):
        raise Exception('Type mismatch: SET')
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
        result.append(('STR', item))
    return result


def operationFINDALL(arg1, arg2):
    # Search all groups of a regex
    # Output a list of groups of entries
    # XXX: If a group has more than 1 entry,
    #    there could be potential problems
    list_of_groups = re.findall(arg1, arg2)
    result = []
    for item in list_of_groups:
        result.append(('STR', item))
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
GROUP_FUNCTION_NAMES = {'COUNT':  0,
                        'SUM':  0,
                        'MUL':  1,
                        'AVG':  0,
                        'MIN': +1e100,
                        'MAX': -1e100,
                        'SSQ': 0,
                        }


def clearGroupIndex():
    #GROUP_INDEX = {}
    for key in GROUP_INDEX.keys():
        del GROUP_INDEX[key]


def getValueFromGroupIndex(GroupIndexKey, FunctionName):
    # If no such index_key, create it and return a value
    if FunctionName not in GROUP_FUNCTION_NAMES:
        raise Exception('No such function \"%s\"' % FunctionName)
    # In the case if NewRowObjectDefault is requested
    if not GroupIndexKey:
        return GROUP_FUNCTION_NAMES[FunctionName]
    if FunctionName not in GROUP_INDEX[GroupIndexKey]['FUNCTIONS']:
        GROUP_INDEX[GroupIndexKey]['FUNCTIONS'][FunctionName] = {}
        GROUP_INDEX[GroupIndexKey]['FUNCTIONS'][FunctionName]['FLAG'] = True
        GROUP_INDEX[GroupIndexKey]['FUNCTIONS'][FunctionName]['VALUE'] = \
            GROUP_FUNCTION_NAMES[FunctionName]
    return GROUP_INDEX[GroupIndexKey]['FUNCTIONS'][FunctionName]['VALUE']


def setValueToGroupIndex(GroupIndexKey, FunctionName, Value):
    GROUP_INDEX[GroupIndexKey]['FUNCTIONS'][FunctionName]['VALUE'] = Value


def initializeGroup(GroupIndexKey):
    if GroupIndexKey not in GROUP_INDEX:
        print('GROUP_DESC[COUNT]='+str(GROUP_DESC['COUNT']))
        GROUP_INDEX[GroupIndexKey] = {}
        GROUP_INDEX[GroupIndexKey]['FUNCTIONS'] = {}
        GROUP_INDEX[GroupIndexKey]['ROWID'] = len(GROUP_INDEX) - 1
    for FunctionName in GROUP_FUNCTION_NAMES:
        # initialize function flags (UpdateFlag)
        if FunctionName in GROUP_INDEX[GroupIndexKey]['FUNCTIONS']:
            GROUP_INDEX[GroupIndexKey]['FUNCTIONS'][FunctionName]['FLAG'] = True
    print('initializeGroup: GROUP_INDEX='+str(GROUP_INDEX))


def groupCOUNT(GroupIndexKey):
    FunctionName = 'COUNT'
    Value = getValueFromGroupIndex(GroupIndexKey, FunctionName)
    if GroupIndexKey:
        if GROUP_INDEX[GroupIndexKey]['FUNCTIONS'][FunctionName]['FLAG']:
            GROUP_INDEX[GroupIndexKey]['FUNCTIONS'][FunctionName]['FLAG'] = False
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


OPERATORS = {\
    # List
    'LIST': lambda args: operationLIST(args),
    # And
    '&': lambda args: operationAND(args),
    '&&': lambda args: operationAND(args),
    'AND': lambda args: operationAND(args),
    # Or
    '|': lambda args: operationOR(args),
    '||': lambda args: operationOR(args),
    'OR': lambda args: operationOR(args),
    # Not
    '!': lambda args: operationNOT(args[0]),
    'NOT': lambda args: operationNOT(args[0]),
    # Between
    'RANGE': lambda args: operationRANGE(args[0], args[1], args[2]),
    'BETWEEN': lambda args: operationRANGE(args[0], args[1], args[2]),
    # Subset
    'IN': lambda args: operationSUBSET(args[0], args[1]),
    'SUBSET': lambda args: operationSUBSET(args[0], args[1]),
    # Less
    '<': lambda args: operationLESS(args),
    'LESS': lambda args: operationLESS(args),
    'LT': lambda args: operationLESS(args),
    # More
    '>': lambda args: operationMORE(args),
    'MORE': lambda args: operationMORE(args),
    'MT': lambda args: operationMORE(args),
    # Less or equal
    '<=': lambda args: operationLESSOREQUAL(args),
    'LESSOREQUAL': lambda args: operationLESSOREQUAL(args),
    'LTE': lambda args: operationLESSOREQUAL(args),
    # More or equal
    '>=': lambda args: operationMOREOREQUAL(args),
    'MOREOREQUAL': lambda args: operationMOREOREQUAL(args),
    'MTE': lambda args: operationMOREOREQUAL(args),
    # Equal
    '=': lambda args: operationEQUAL(args),
    '==': lambda args: operationEQUAL(args),
    'EQ': lambda args: operationEQUAL(args),
    'EQUAL': lambda args: operationEQUAL(args),
    'EQUALS': lambda args: operationEQUAL(args),
    # Not equal
    '!=': lambda args: operationNOTEQUAL(args[0], args[1]),
    '<>': lambda args: operationNOTEQUAL(args[0], args[1]),
    '~=': lambda args: operationNOTEQUAL(args[0], args[1]),
    'NE': lambda args: operationNOTEQUAL(args[0], args[1]),
    'NOTEQUAL': lambda args: operationNOTEQUAL(args[0], args[1]),
    # Plus
    '+': lambda args: operationSUM(args),
    'SUM': lambda args: operationSUM(args),
    # Minus
    '-': lambda args: operationDIFF(args[0], args[1]),
    'DIFF': lambda args: operationDIFF(args[0], args[1]),
    # Mul
    '*': lambda args: operationMUL(args),
    'MUL': lambda args: operationMUL(args),
    # Div
    '/': lambda args: operationDIV(args[0], args[1]),
    'DIV': lambda args: operationDIV(args[0], args[1]),
    # Regexp match
    'MATCH': lambda args: operationMATCH(args[0], args[1]),
    'LIKE': lambda args: operationMATCH(args[0], args[1]),
    # Regexp search
    'SEARCH': lambda args: operationSEARCH(args[0], args[1]),
    # Regexp findal
    'FINDALL': lambda args: operationFINDALL(args[0], args[1]),
    # Group count
    'COUNT': lambda args: groupCOUNT(GroupIndexKey),
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
        if head in set(['STR', 'STRING']):  # one arg
            return operationSTR(root[1])
        elif head in set(['SET']):
            return operationSET(root[1])
        tail = root[1:]
        args = []
        # evaluate arguments recursively
        for element in tail:  # resolve tree by recursion
            args.append(evaluateExpression(
                element, VarDictionary, GroupIndexKey))
        # call functions with evaluated arguments
        try:
            return OPERATORS[head](args)
        except KeyError:
            raise Exception('Unknown operator: %s' % head)
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
    #VarDictionary = getVarDictionary(RowObject)
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
        return '%10d'
    elif Type is float:
        return '%25.15E'
    elif Type is str:
        return '%20s'
    elif Type is bool:
        return '%2d'
    else:
        raise Exception('Unknown type')


def getDefaultValue(Type):
    if Type is int:
        return 0
    elif Type is float:
        return 0.0
    elif Type is str:
        return ''
    elif Type is bool:
        return False
    else:
        raise Exception('Unknown type')

# VarDictionary = Context (this name is more suitable)

# GroupIndexKey is a key to special structure/dictionary GROUP_INDEX.
# GROUP_INDEX contains information needed to calculate streamed group functions
#  such as COUNT, AVG, MIN, MAX etc...


def newRowObject(ParameterNames, RowObject, VarDictionary, ContextFormat, GroupIndexKey=None):
    # Return a subset of RowObject according to
    # ParameterNames include either parnames
    #  or expressions containing parnames literals
    # ContextFormat contains format for ParNames
    anoncount = 0
    RowObjectNew = []
    for expr in ParameterNames:
        if type(expr) in set([list, tuple]):  # bind
            head = expr[0]
            if head in set(['let', 'bind', 'LET', 'BIND']):
                par_name = expr[1]
                par_expr = expr[2]
            else:
                par_name = "#%d" % anoncount
                anoncount += 1
                par_expr = expr
            par_value = evaluateExpression(
                par_expr, VarDictionary, GroupIndexKey)
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

QUERY_BUFFER = '__BUFFER__'


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
    print('-----------------------------------------')
    print(TableName+' summary:')
    try:
        print('-----------------------------------------')
        print('Comment: \n'+LOCAL_TABLE_CACHE[TableName]['header']['comment'])
    except:
        pass
    print('Number of rows: ' +
          str(LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']))
    print('Table type: ' +
          str(LOCAL_TABLE_CACHE[TableName]['header']['table_type']))
    print('-----------------------------------------')
    print('            PAR_NAME           PAR_FORMAT')
    print('')
    for par_name in LOCAL_TABLE_CACHE[TableName]['header']['order']:
        par_format = LOCAL_TABLE_CACHE[TableName]['header']['format'][par_name]
        print('%20s %20s' % (par_name, par_format))
    print('-----------------------------------------')

# Write a table to File or STDOUT


def outputTable(TableName, Conditions=None, File=None, Header=True):
    # Display or record table with condition checking
    if File:
        Header = False
        OutputFile = open(File, 'w')
    if Header:
        headstr = putTableHeaderToString(TableName)
        if File:
            OutputFile.write(headstr)
        else:
            print(headstr)
    for RowID in range(0, LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']):
        RowObject = getRowObject(RowID, TableName)
        VarDictionary = getVarDictionary(RowObject)
        VarDictionary['LineNumber'] = RowID
        if not checkRowObject(RowObject, Conditions, VarDictionary):
            continue
        raw_string = putRowObjectToString(RowObject)
        if File:
            OutputFile.write(raw_string+'\n')
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
    LOCAL_TABLE_CACHE[TableName]['header'] = {}
    LOCAL_TABLE_CACHE[TableName]['header']['order'] = header_order
    LOCAL_TABLE_CACHE[TableName]['header']['format'] = header_format
    LOCAL_TABLE_CACHE[TableName]['header']['default'] = header_default
    LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows'] = 0
    LOCAL_TABLE_CACHE[TableName]['header']['size_in_bytes'] = 0
    LOCAL_TABLE_CACHE[TableName]['header']['table_name'] = TableName
    LOCAL_TABLE_CACHE[TableName]['header']['table_type'] = 'column-fixed'
    LOCAL_TABLE_CACHE[TableName]['data'] = data


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
        #LOCAL_TABLE_CACHE[TableName] = {}
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
    return LOCAL_TABLE_CACHE[TableName]['data'][ParameterName]

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
        Columns.append(LOCAL_TABLE_CACHE[TableName]['data'][par_name])
    return Columns


def addColumn(TableName, ParameterName, Before=None, Expression=None, Type=None, Default=None, Format=None):
    if ParameterName in LOCAL_TABLE_CACHE[TableName]['header']['format']:
        raise Exception('Column \"%s\" already exists' % ParameterName)
    if not Type:
        Type = float
    if not Default:
        Default = getDefaultValue(Type)
    if not Format:
        Format = getDefaultFormat(Type)
    number_of_rows = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    # Mess with data
    if not Expression:
        LOCAL_TABLE_CACHE[TableName]['data'][ParameterName] = [
            Default for i in range(0, number_of_rows)]
    else:
        data = []
        for RowID in range(0, number_of_rows):
            RowObject = getRowObject(RowID, TableName)
            VarDictionary = getVarDictionary(RowObject)
            VarDictionary['LineNumber'] = RowID
            par_value = evaluateExpression(Expression, VarDictionary)
            data.append(par_value)
            LOCAL_TABLE_CACHE[TableName]['data'][ParameterName] = data
    # Mess with header
    header_order = LOCAL_TABLE_CACHE[TableName]['header']['order']
    if not Before:
        header_order.append(ParameterName)
    else:
        #i = 0
        # for par_name in header_order:
        #    if par_name == Before: break
        #    i += 1
        i = header_order.index(Before)
        header_order = header_order[:i] + [ParameterName, ] + header_order[i:]
    LOCAL_TABLE_CACHE[TableName]['header']['order'] = header_order
    LOCAL_TABLE_CACHE[TableName]['header']['format'][ParameterName] = Format
    LOCAL_TABLE_CACHE[TableName]['header']['default'][ParameterName] = Default


def deleteColumn(TableName, ParameterName):
    if ParameterName not in LOCAL_TABLE_CACHE[TableName]['header']['format']:
        raise Exception('No such column \"%s\"' % ParameterName)
    # Mess with data
    i = LOCAL_TABLE_CACHE[TableName]['header']['order'].index(ParameterName)
    del LOCAL_TABLE_CACHE[TableName]['header']['order'][i]
    del LOCAL_TABLE_CACHE[TableName]['header']['format'][ParameterName]
    del LOCAL_TABLE_CACHE[TableName]['header']['default'][ParameterName]
    if not LOCAL_TABLE_CACHE[TableName]['header']['order']:
        LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows'] = 0
    # Mess with header
    del LOCAL_TABLE_CACHE[TableName]['data'][ParameterName]


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
    #table_columns = LOCAL_TABLE_CACHE[TableName]['data'].keys()
    #table_length = len(TableObject['header']['number_of_rows'])
    # if ParameterNames=='*':
    #   ParameterNames = table_columns
    # check if Conditions contain elements which are not in the TableObject
    #condition_variables = getConditionVariables(Conditions)
    #strange_pars = set(condition_variables)-set(table_variables)
    # if strange_pars:
    #   raise Exception('The following parameters are not in the table \"%s\"' % (TableName,list(strange_pars)))
    # do full scan each time
    if DestinationTableName == TableName:
        raise Exception('Selecting into source table is forbidden')
    table_length = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    row_count = 0
    for RowID in range(0, table_length):
        RowObject = getRowObject(RowID, TableName)
        VarDictionary = getVarDictionary(RowObject)
        VarDictionary['LineNumber'] = RowID
        ContextFormat = getContextFormat(RowObject)
        RowObjectNew = newRowObject(
            ParameterNames, RowObject, VarDictionary, ContextFormat)
        if checkRowObject(RowObject, Conditions, VarDictionary):
            addRowObject(RowObjectNew, DestinationTableName)
            row_count += 1
    LOCAL_TABLE_CACHE[DestinationTableName]['header']['number_of_rows'] += row_count


def length(TableName):
    tab_len = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    #print(str(tab_len)+' rows in '+TableName)
    return tab_len

# Select parameters from a table with certain conditions.
# Parameters can be the names or expressions.
# Conditions contain a list of expressions in a special language.
# Set Output to False to suppress output
# Set File=FileName to redirect output to a file.


def select(TableName, DestinationTableName=QUERY_BUFFER, ParameterNames=None, Conditions=None, Output=True, File=None):
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
            '%s: no such table. Check tableList() for more info.' % TableName)
    if not ParameterNames:
        ParameterNames = LOCAL_TABLE_CACHE[TableName]['header']['order']
    # clear QUERY_BUFFER for the new result
    LOCAL_TABLE_CACHE[DestinationTableName] = {}
    RowObjectDefault = getDefaultRowObject(TableName)
    VarDictionary = getVarDictionary(RowObjectDefault)
    ContextFormat = getContextFormat(RowObjectDefault)
    RowObjectDefaultNew = newRowObject(
        ParameterNames, RowObjectDefault, VarDictionary, ContextFormat)
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
        LOCAL_TABLE_CACHE[DestinationTableName]['header'] = LOCAL_TABLE_CACHE[TableName]['header']
        LOCAL_TABLE_CACHE[DestinationTableName]['data'] = {}
    LOCAL_TABLE_CACHE[DestinationTableName]['header']['number_of_rows'] = len(
        RowIDList)
    # print 'AT: RowIDList = '+str(RowIDList)
    for par_name in LOCAL_TABLE_CACHE[DestinationTableName]['header']['order']:
        par_data = LOCAL_TABLE_CACHE[TableName]['data'][par_name]
        LOCAL_TABLE_CACHE[DestinationTableName]['data'][par_name] = [
            par_data[i] for i in RowIDList]


def compareLESS(RowObject1, RowObject2, ParameterNames):
    # print 'CL/'
    # arg1 and arg2 are RowObjects
    # Compare them according to ParameterNames
    # Simple validity check:
    # if len(arg1) != len(arg2):
    #   raise Exception('Arguments have different lengths')
    #RowObject1Subset = subsetOfRowObject(ParameterNames,RowObject1)
    #RowObject2Subset = subsetOfRowObject(ParameterNames,RowObject2)
    # return RowObject1Subset < RowObject2Subset
    row1 = []
    row2 = []
    #n = len(RowObject1)
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
        #pivot = lst[0]
        #lesser = quickSort([x for x in lst[1:] if x < pivot])
        #greater = quickSort([x for x in lst[1:] if x >= pivot])
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
        greater = quickSort(greater_index, TableName,
                            ParameterNames, Accending)
        # return lesser + [pivot_index] + greater
        if Accending:
            return lesser + [PivotID] + greater
        else:
            return greater + [PivotID] + lesser

# Sorting must work well on the table itself!


def sort(TableName, DestinationTableName=None, ParameterNames=None, Accending=True, Output=False, File=None):
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
    number_of_rows = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    index = list(range(0, number_of_rows))
    # print 'num = '+str(number_of_rows)
    if not DestinationTableName:
        DestinationTableName = TableName
    # if names are not provided use all parameters in sorting
    if not ParameterNames:
        ParameterNames = LOCAL_TABLE_CACHE[TableName]['header']['order']
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


def group(TableName, DestinationTableName=QUERY_BUFFER, ParameterNames=None, GroupParameterNames=None, Output=True):
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
        raise Exception('TableName and DestinationTableName must be different')
    #if not ParameterNames: ParameterNames=LOCAL_TABLE_CACHE[TableName]['header']['order']
    # Prepare the new DestinationTable
    RowObjectDefault = getDefaultRowObject(TableName)
    VarDictionary = getVarDictionary(RowObjectDefault)
    ContextFormat = getContextFormat(RowObjectDefault)
    RowObjectDefaultNew = newRowObject(
        ParameterNames, RowObjectDefault, VarDictionary, ContextFormat)
    dropTable(DestinationTableName)  # redundant
    createTable(DestinationTableName, RowObjectDefaultNew)
    # Loop through rows of source Table
    # On each iteration group functions update GROUP_INDEX (see description above)
    number_of_rows = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    # STAGE 1: CREATE GROUPS
    print('LOOP:')
    for RowID in range(0, number_of_rows):
        print('--------------------------------')
        print('RowID='+str(RowID))
        # RowObject from source table
        RowObject = getRowObject(RowID, TableName)
        VarDictionary = getVarDictionary(RowObject)
        print('VarDictionary='+str(VarDictionary))
        # This is a trick which makes evaluateExpression function
        #   not consider first expression as an operation
        GroupParameterNames_ = ['LIST'] + list(GroupParameterNames)
        GroupIndexKey = evaluateExpression(GroupParameterNames_, VarDictionary)
        # List is an unhashable type in Python!
        GroupIndexKey = tuple(GroupIndexKey)
        initializeGroup(GroupIndexKey)
        print('GROUP_INDEX='+str(GROUP_INDEX))
        ContextFormat = getContextFormat(RowObject)
        RowObjectNew = newRowObject(
            ParameterNames, RowObject, VarDictionary, ContextFormat, GroupIndexKey)
        RowIDGroup = GROUP_INDEX[GroupIndexKey]['ROWID']
        setRowObject(RowIDGroup, RowObjectNew, DestinationTableName)
    # Output result if required
    if Output and DestinationTableName == QUERY_BUFFER:
        outputTable(DestinationTableName, File=File)

# /GROUPING =========================================================

# EXTRACTING ========================================================


REGEX_INTEGER = '[+-]?\d+'
REGEX_STRING = '[^\s]+'
REGEX_FLOAT_F = '[+-]?\d*\.?\d+'
REGEX_FLOAT_E = '[+-]?\d*\.?\d+[eEfF]?[+-]?\d+'


def REGEX_INTEGER_FIXCOL(n): return '\d{%d}' % n


def REGEX_STRING_FIXCOL(n): return '[^\s]{%d}' % n


def REGEX_FLOAT_F_FIXCOL(n): return '[\+\-\.\d]{%d}' % n


def REGEX_FLOAT_E_FIXCOL(n): return '[\+\-\.\deEfF]{%d}' % n

# Extract sub-columns from string column


def extractColumns(TableName, SourceParameterName, ParameterFormats, ParameterNames=None, FixCol=False):
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
    if type(LOCAL_TABLE_CACHE[TableName]['header']['default'][SourceParameterName]) not in set([str, six.text_type]):
        raise Exception('Source parameter must be a string')
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
                par_name = '#%d' % i
                fmt = LOCAL_TABLE_CACHE[TableName]['header']['format'].get(
                    par_name, None)
                if not fmt:
                    break
            ParameterNames.append(par_name)
    # check if ParameterNames are valid
    Intersection = set(ParameterNames).intersection(
        LOCAL_TABLE_CACHE[TableName]['header']['order'])
    if Intersection:
        raise Exception('Parameters %s already exist' %
                        str(list(Intersection)))
    # loop over ParameterNames to prepare LOCAL_TABLE_CACHE
    i = 0
    for par_name in ParameterNames:
        par_format = ParameterFormats[i]
        LOCAL_TABLE_CACHE[TableName]['header']['format'][par_name] = par_format
        LOCAL_TABLE_CACHE[TableName]['data'][par_name] = []
        i += 1
    # append new parameters in order list
    LOCAL_TABLE_CACHE[TableName]['header']['order'] += ParameterNames
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
        if ty == 'd':
            par_type = int
            if FixCol:
                format_regex_part = REGEX_INTEGER_FIXCOL(lng)
            else:
                format_regex_part = REGEX_INTEGER
        elif ty == 's':
            par_type = str
            if FixCol:
                format_regex_part = REGEX_STRING_FIXCOL(lng)
            else:
                format_regex_part = REGEX_STRING
        elif ty == 'f':
            par_type = float
            if FixCol:
                format_regex_part = REGEX_FLOAT_F_FIXCOL(lng)
            else:
                format_regex_part = REGEX_FLOAT_F
        elif ty == 'e':
            par_type = float
            if FixCol:
                format_regex_part = REGEX_FLOAT_E_FIXCOL(lng)
            else:
                format_regex_part = REGEX_FLOAT_E
        else:
            raise Exception('Unknown data type')
        format_regex.append('('+format_regex_part+')')
        format_types.append(par_type)
        def_val = getDefaultValue(par_type)
        LOCAL_TABLE_CACHE[TableName]['header']['default'][par_name] = def_val
        i += 1
    format_regex = '\s*'.join(format_regex)
    # print 'format_regex='+str(format_regex)
    # return format_regex
    # loop through values of SourceParameter
    for SourceParameterString in LOCAL_TABLE_CACHE[TableName]['data'][SourceParameterName]:
        try:
            ExtractedValues = list(
                re.search(format_regex, SourceParameterString).groups())
        except:
            raise Exception('Error with line \"%s\"' % SourceParameterString)
        i = 0
        # loop through all parameters which are supposed to be extracted
        for par_name in ParameterNames:
            # print 'ExtractedValues[i]='+ExtractedValues[i]
            # print 'par_name='+par_name
            par_value = format_types[i](ExtractedValues[i])
            LOCAL_TABLE_CACHE[TableName]['data'][par_name].append(par_value)
            i += 1
    # explicitly check that number of rows are equal
    number_of_rows = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
    number_of_rows2 = len(
        LOCAL_TABLE_CACHE[TableName]['data'][SourceParameterName])
    number_of_rows3 = len(
        LOCAL_TABLE_CACHE[TableName]['data'][ParameterNames[0]])
    if not (number_of_rows == number_of_rows2 == number_of_rows3):
        raise Exception('Error while extracting parameters: check your regexp')

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
PARLIST_DOTPAR = ['par_line', ]
PARLIST_ID = ['trans_id', ]
PARLIST_STANDARD = ['molec_id', 'local_iso_id', 'nu', 'sw', 'a', 'elower', 'gamma_air',
                    'delta_air', 'gamma_self', 'n_air', 'n_self', 'gp', 'gpp']
PARLIST_LABELS = ['statep', 'statepp']
PARLIST_LINEMIXING = ['y_air', 'y_self']

PARLIST_VOIGT_AIR = ['gamma_air', 'delta_air', 'deltap_air', 'n_air']
PARLIST_VOIGT_SELF = ['gamma_self', 'delta_self', 'deltap_self', 'n_self']
PARLIST_VOIGT_H2 = ['gamma_H2', 'delta_H2', 'deltap_H2', 'n_H2']
PARLIST_VOIGT_CO2 = ['gamma_CO2', 'delta_CO2', 'n_CO2']
PARLIST_VOIGT_HE = ['gamma_He', 'delta_He', 'n_He']
PARLIST_VOIGT_ALL = mergeParlist(PARLIST_VOIGT_AIR, PARLIST_VOIGT_SELF,
                                 PARLIST_VOIGT_H2, PARLIST_VOIGT_CO2,
                                 PARLIST_VOIGT_HE)

PARLIST_SDVOIGT_AIR = ['gamma_air', 'delta_air',
                       'deltap_air', 'n_air', 'SD_air']
PARLIST_SDVOIGT_SELF = ['gamma_self', 'delta_self',
                        'deltap_self', 'n_self', 'SD_self']
PARLIST_SDVOIGT_H2 = []
PARLIST_SDVOIGT_CO2 = []
PARLIST_SDVOIGT_HE = []
PARLIST_SDVOIGT_ALL = mergeParlist(PARLIST_SDVOIGT_AIR, PARLIST_SDVOIGT_SELF,
                                   PARLIST_SDVOIGT_H2, PARLIST_SDVOIGT_CO2,
                                   PARLIST_SDVOIGT_HE)

PARLIST_GALATRY_AIR = ['gamma_air', 'delta_air',
                       'deltap_air', 'n_air', 'beta_g_air']
PARLIST_GALATRY_SELF = ['gamma_self', 'delta_self',
                        'deltap_self', 'n_self', 'beta_g_self']
PARLIST_GALATRY_H2 = []
PARLIST_GALATRY_CO2 = []
PARLIST_GALATRY_HE = []
PARLIST_GALATRY_ALL = mergeParlist(PARLIST_GALATRY_AIR, PARLIST_GALATRY_SELF,
                                   PARLIST_GALATRY_H2, PARLIST_GALATRY_CO2,
                                   PARLIST_GALATRY_HE)

PARLIST_ALL = mergeParlist(PARLIST_ID, PARLIST_DOTPAR, PARLIST_STANDARD,
                           PARLIST_LABELS, PARLIST_LINEMIXING, PARLIST_VOIGT_ALL,
                           PARLIST_SDVOIGT_ALL, PARLIST_GALATRY_ALL)

PARAMETER_GROUPS = {
    'par_line': PARLIST_DOTPAR,
    '160-char': PARLIST_DOTPAR,
    '.par': PARLIST_DOTPAR,
    'id': PARLIST_ID,
    'standard': PARLIST_STANDARD,
    'labels': PARLIST_LABELS,
    'linemixing': PARLIST_LINEMIXING,
    'voigt_air': PARLIST_VOIGT_AIR,
    'voigt_self': PARLIST_VOIGT_SELF,
    'voigt_h2': PARLIST_VOIGT_H2,
    'voigt_co2': PARLIST_VOIGT_CO2,
    'voigt_he': PARLIST_VOIGT_HE,
    'voigt': PARLIST_VOIGT_ALL,
    'sdvoigt_air': PARLIST_SDVOIGT_AIR,
    'sdvoigt_self': PARLIST_SDVOIGT_SELF,
    'sdvoigt_h2': PARLIST_SDVOIGT_H2,
    'sdvoigt_co2': PARLIST_SDVOIGT_CO2,
    'sdvoigt_he': PARLIST_SDVOIGT_HE,
    'sdvoigt': PARLIST_SDVOIGT_ALL,
    'galatry_air': PARLIST_GALATRY_AIR,
    'galatry_self': PARLIST_GALATRY_SELF,
    'galatry_h2': PARLIST_GALATRY_H2,
    'galatry_co2': PARLIST_GALATRY_CO2,
    'galatry_he': PARLIST_GALATRY_HE,
    'galatry': PARLIST_GALATRY_ALL,
    'all': PARLIST_ALL
}


def prepareParlist(pargroups=[], params=[], dotpar=True):
    # Apply defaults
    parlist_default = []
    if dotpar:
        parlist_default += ['par_line']
    #parlist_default += PARAMETER_GROUPS['id']

    # Make a dictionary of "assumed" parameters.
    ASSUMED_PARAMS = {}
    if 'par_line' in set(parlist_default):
        ASSUMED_PARAMS = HITRAN_DEFAULT_HEADER['format']

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
    HEADER = {'table_name': '', 'number_of_rows': -1, 'format': {},
              'default': {}, 'table_type': 'column-fixed',
              'size_in_bytes': -1, 'order': [], 'description': {}}

    # Add column-fixed 160-character part, if specified in parlist.
    if 'par_line' in set(parlist):
        HEADER['order'] = HITRAN_DEFAULT_HEADER['order']
        HEADER['format'] = HITRAN_DEFAULT_HEADER['format']
        HEADER['default'] = HITRAN_DEFAULT_HEADER['default']
        HEADER['description'] = HITRAN_DEFAULT_HEADER['description']

    # Insert all other parameters in the "extra" section of the header.
    #while 'par_line' in parlist: parlist.remove('par_line')
    plist = [v for v in parlist if v != 'par_line']
    HEADER['extra'] = []
    HEADER['extra_format'] = {}
    HEADER['extra_separator'] = ','
    for param in plist:
        param = param.lower()
        HEADER['extra'].append(param)
        HEADER['extra_format'][param] = PARAMETER_META[param]['default_fmt']

    return HEADER


def queryHITRAN(TableName, iso_id_list, numin, numax, pargroups=[], params=[], dotpar=True, head=False):
    #import httplib
    #conn = httplib.HTTPConnection('hitranazure.cloudapp.com')
    # conn.Request('')
    #r = conn.getresponse()
    # print r.status, r.reason
    #data1 = data1.read
    # TableHeader = HITRAN_DEFAULT_HEADER # deprecated
    ParameterList = prepareParlist(
        pargroups=pargroups, params=params, dotpar=dotpar)
    TableHeader = prepareHeader(ParameterList)
    TableHeader['table_name'] = TableName
    DataFileName = VARIABLES['BACKEND_DATABASE_NAME'] + \
        '/' + TableName + '.data'
    HeaderFileName = VARIABLES['BACKEND_DATABASE_NAME'] + \
        '/' + TableName + '.header'
    # if TableName in LOCAL_TABLE_CACHE.keys():
    #   raise Exception('Table \"%s\" exists' % TableName)
    # if os.path.isfile(DataFileName):
    #   raise Exception('File \"%s\" exists' % DataFileName)
    # if os.path.isfile(HeaderFileName):
    #   raise Exception('!!File \"%s\" exists' % HeaderFileName)
    # create URL
    iso_id_list_str = [str(iso_id) for iso_id in iso_id_list]
    iso_id_list_str = ','.join(iso_id_list_str)
    # url = 'http://hitran.cloudapp.net' + '/lbl/5?' + \
    # url = 'http://hitranazure.cloudapp.net' + '/lbl/5?' + \
    # 'iso_ids_list=' + iso_id_list_str + '&' + \
    # 'numin=' + str(numin) + '&' + \
    # 'numax=' + str(numax) + '&' + \
    # 'access=api' + '&' + \
    #'key=' + GLOBAL_HITRAN_APIKEY
    if pargroups or params:  # custom par search
        url = GLOBAL_HOST + '/lbl/api?' + \
            'iso_ids_list=' + iso_id_list_str + '&' + \
            'numin=' + str(numin) + '&' + \
            'numax=' + str(numax) + '&' + \
            'head=' + str(head) + '&' + \
            'fixwidth=0&sep=[comma]&' +\
            'request_params=' + ','.join(ParameterList)
    else:  # old-fashioned .par search
        url = GLOBAL_HOST + '/lbl/api?' + \
            'iso_ids_list=' + iso_id_list_str + '&' + \
            'numin=' + str(numin) + '&' + \
            'numax=' + str(numax)
    #raise Exception(url)
    # Download data by chunks.
    try:
        req = six.moves.urllib.request.urlopen(url)
    except six.moves.urllib.error.HTTPError:
        raise Exception('Failed to retrieve data for given parameters.')
    except six.moves.urllib.error.URLError:
        raise Exception(
            'Cannot connect to %s. Try again or edit GLOBAL_HOST variable.' % GLOBAL_HOST)
    # CHUNK = 16 * 1024 # default value
    CHUNK = 64 * 1024
    print('BEGIN DOWNLOAD: '+TableName)
    with open(DataFileName, 'w') as fp:
        while True:
            chunk = req.read(CHUNK)
            if not chunk:
                break
            fp.write(chunk.decode('utf-8'))
            print('  %d bytes written to %s' % (CHUNK, DataFileName))
    with open(HeaderFileName, 'w') as fp:
        fp.write(json.dumps(TableHeader, indent=2))
        print('Header written to %s' % HeaderFileName)
    print('END DOWNLOAD')
    # Set comment
    # Get this table to LOCAL_TABLE_CACHE
    storage2cache(TableName)
    print('PROCESSED')


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

    dbname, tablename, nodename = NewTablePath.split('::')
    dbname1, tablename1, nodename1 = SourceTablePath.split('::')

    if not NODE_READY:
        raise Exception(
            'Node \"%s\" is not ready. Call nodeInit()' % NODE_NAME)

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


ISO_ID_INDEX = {

    'M': 0,
    'I': 1,
    'iso_name': 2,
    'abundance': 3,
    'mass': 4,
    'mol_name': 5

}

#    id           M    I    iso_name                    abundance           mass        mol_name

ISO_ID = {

    1: [1,   1,  'H2(16O)',                   0.997317,           18.010565,      'H2O'],
    2: [1,   2,  'H2(18O)',                   0.00199983,         20.014811,      'H2O'],
    3: [1,   3,  'H2(17O)',                   0.000372,           19.01478,       'H2O'],
    4: [1,   4,  'HD(16O)',                   0.00031069,         19.01674,       'H2O'],
    5: [1,   5,  'HD(18O)',                   0.000000623,        21.020985,      'H2O'],
    6: [1,   6,  'HD(17O)',                   0.000000116,        20.020956,      'H2O'],
    7: [2,   1,  '(12C)(16O)2',               0.9842,             43.98983,       'CO2'],
    8: [2,   2,  '(13C)(16O)2',               0.01106,            44.993185,      'CO2'],
    9: [2,   3,  '(16O)(12C)(18O)',           0.0039471,          45.994076,      'CO2'],
    10: [2,   4,  '(16O)(12C)(17O)',           0.000734,           44.994045,      'CO2'],
    11: [2,   5,  '(16O)(13C)(18O)',           0.00004434,         46.997431,      'CO2'],
    12: [2,   6,  '(16O)(13C)(17O)',           0.00000825,         45.9974,        'CO2'],
    13: [2,   7,  '(12C)(18O)2',               0.0000039573,       47.998322,      'CO2'],
    14: [2,   8,  '(17O)(12C)(18O)',           0.00000147,         46.998291,      'CO2'],
    15: [2,   0,  '(13C)(18O)2',               0.000000044967,     49.001675,      'CO2'],
    120: [2,  11,  '(18O)(13C)(17O)',           0.00000001654,      48.00165,       'CO2'],
    121: [2,   9,  '(12C)(17O)2',               0.0000001368,       45.998262,      'CO2'],
    16: [3,   1,  '(16O)3',                    0.992901,           47.984745,      'O3'],
    17: [3,   2,  '(16O)(16O)(18O)',           0.00398194,         49.988991,      'O3'],
    18: [3,   3,  '(16O)(18O)(16O)',           0.00199097,         49.988991,      'O3'],
    19: [3,   4,  '(16O)(16O)(17O)',           0.00074,            48.98896,       'O3'],
    20: [3,   5,  '(16O)(17O)(16O)',           0.00037,            48.98896,       'O3'],
    21: [4,   1,  '(14N)2(16O)',               0.990333,           44.001062,      'N2O'],
    22: [4,   2,  '(14N)(15N)(16O)',           0.0036409,          44.998096,      'N2O'],
    23: [4,   3,  '(15N)(14N)(16O)',           0.0036409,          44.998096,      'N2O'],
    24: [4,   4,  '(14N)2(18O)',               0.00198582,         46.005308,      'N2O'],
    25: [4,   5,  '(14N)2(17O)',               0.000369,           45.005278,      'N2O'],
    26: [5,   1,  '(12C)(16O)',                0.98654,            27.994915,      'CO'],
    27: [5,   2,  '(13C)(16O)',                0.01108,            28.99827,       'CO'],
    28: [5,   3,  '(12C)(18O)',                0.0019782,          29.999161,      'CO'],
    29: [5,   4,  '(12C)(17O)',                0.000368,           28.99913,       'CO'],
    30: [5,   5,  '(13C)(18O)',                0.00002222,         31.002516,      'CO'],
    31: [5,   6,  '(13C)(17O)',                0.00000413,         30.002485,      'CO'],
    32: [6,   1,  '(12C)H4',                   0.98827,            16.0313,        'CH4'],
    33: [6,   2,  '(13C)H4',                   0.0111,             17.034655,      'CH4'],
    34: [6,   3,  '(12C)H3D',                  0.00061575,         17.037475,      'CH4'],
    35: [6,   4,  '(13C)H3D',                  0.0000049203,       18.04083,       'CH4'],
    36: [7,   1,  '(16O)2',                    0.995262,           31.98983,       'O2'],
    37: [7,   2,  '(16O)(18O)',                0.00399141,         33.994076,      'O2'],
    38: [7,   3,  '(16O)(17O)',                0.000742,           32.994045,      'O2'],
    39: [8,   1,  '(14N)(16O)',                0.993974,           29.997989,      'NO'],
    40: [8,   2,  '(15N)(16O)',                0.0036543,          30.995023,      'NO'],
    41: [8,   3,  '(14N)(18O)',                0.00199312,         32.002234,      'NO'],
    42: [9,   1,  '(32S)(16O)2',               0.94568,            63.961901,      'SO2'],
    43: [9,   2,  '(34S)(16O)2',               0.04195,            65.957695,      'SO2'],
    44: [10,   1,  '(14N)(16O)2',               0.991616,           45.992904,      'NO2'],
    45: [11,   1,  '(14N)H3',                   0.9958715,          17.026549,      'NH3'],
    46: [11,   2,  '(15N)H3',                   0.0036613,          18.023583,      'NH3'],
    47: [12,   1,  'H(14N)(16O)3',              0.98911,            62.995644,      'HNO3'],
    117: [12,   2,  'H(15N)(16O)3',              0.003636,           63.99268,      'HNO3'],
    48: [13,   1,  '(16O)H',                    0.997473,           17.00274,       'OH'],
    49: [13,   2,  '(18O)H',                    0.00200014,         19.006986,      'OH'],
    50: [13,   3,  '(16O)D',                    0.00015537,         18.008915,      'OH'],
    51: [14,   1,  'H(19F)',                    0.99984425,         20.006229,      'HF'],
    110: [14,   2,  'D(19F)',                    0.000115,           21.0125049978,  'HF'],
    52: [15,   1,  'H(35Cl)',                   0.757587,           35.976678,      'HCl'],
    53: [15,   2,  'H(37Cl)',                   0.242257,           37.973729,      'HCl'],
    107: [15,   3,  'D(35Cl)',                   0.000118005,        36.9829544578,  'HCl'],
    108: [15,   4,  'D(37Cl)',                   0.000037735,        38.9800043678,  'HCl'],
    54: [16,   1,  'H(79Br)',                   0.50678,            79.92616,       'HBr'],
    55: [16,   2,  'H(81Br)',                   0.49306,            81.924115,      'HBr'],
    111: [16,   3,  'D(79Br)',                   0.0000582935,       80.9324388778,  'HBr'],
    112: [16,   4,  'D(81Br)',                   0.0000567065,       82.9303923778,  'HBr'],
    56: [17,   1,  'H(127I)',                   0.99984425,         127.912297,     'HI'],
    113: [17,   2,  'D(127I)',                   0.000115,           128.918574778,  'HI'],
    57: [18,   1,  '(35Cl)(16O)',               0.75591,            50.963768,      'ClO'],
    58: [18,   2,  '(37Cl)(16O)',               0.24172,            52.960819,      'ClO'],
    59: [19,   1,  '(16O)(12C)(32S)',           0.93739,            59.966986,      'OCS'],
    60: [19,   2,  '(16O)(12C)(34S)',           0.04158,            61.96278,       'OCS'],
    61: [19,   3,  '(16O)(13C)(32S)',           0.01053,            60.970341,      'OCS'],
    62: [19,   4,  '(16O)(12C)(33S)',           0.01053,            60.966371,      'OCS'],
    63: [19,   5,  '(18O)(12C)(32S)',           0.00188,            61.971231,      'OCS'],
    64: [20,   1,  'H2(12C)(16O)',              0.98624,            30.010565,      'H2CO'],
    65: [20,   2,  'H2(13C)(16O)',              0.01108,            31.01392,       'H2CO'],
    66: [20,   3,  'H2(12C)(18O)',              0.0019776,          32.014811,      'H2CO'],
    67: [21,   1,  'H(16O)(35Cl)',              0.75579,            51.971593,      'HOCl'],
    68: [21,   2,  'H(16O)(37Cl)',              0.24168,            53.968644,      'HOCl'],
    69: [22,   1,  '(14N)2',                    0.9926874,          28.006147,      'N2'],
    118: [22,   2,  '(14N)(15N)',                0.0072535,          29.997989,      'N2'],
    70: [23,   1,  'H(12C)(14N)',               0.98511,            27.010899,      'HCN'],
    71: [23,   2,  'H(13C)(14N)',               0.01107,            28.014254,      'HCN'],
    72: [23,   3,  'H(12C)(15N)',               0.0036217,          28.007933,      'HCN'],
    73: [24,   1,  '(12C)H3(35Cl)',             0.74894,            49.992328,      'CH3Cl'],
    74: [24,   2,  '(12C)H3(37Cl)',             0.23949,            51.989379,      'CH3Cl'],
    75: [25,   1,  'H2(16O)2',                  0.994952,           34.00548,       'H2O2'],
    76: [26,   1,  '(12C)2H2',                  0.9776,             26.01565,       'C2H2'],
    77: [26,   2,  '(12C)(13C)H2',              0.02197,            27.019005,      'C2H2'],
    105: [26,   3,  '(12C)2HD',                  0.00030455,         27.021825,      'C2H2'],
    78: [27,   1,  '(12C)2H6',                  0.97699,            30.04695,       'C2H6'],
    106: [27,   2,  '(12C)H3(13C)H3',            0.021952611,        31.050305,      'C2H6'],
    79: [28,   1,  '(31P)H3',                   0.99953283,         33.997238,      'PH3'],
    80: [29,   1,  '(12C)(16O)(19F)2',          0.98654,            65.991722,      'COF2'],
    119: [29,   2,  '(13C)(16O)(19F)2',          0.0110834,          66.995083,      'COF2'],
    81: [31,   1,  'H2(32S)',                   0.94988,            33.987721,      'H2S'],
    82: [31,   2,  'H2(34S)',                   0.04214,            35.983515,      'H2S'],
    83: [31,   3,  'H2(33S)',                   0.007498,           34.987105,      'H2S'],
    84: [32,   1,  'H(12C)(16O)(16O)H',         0.983898,           46.00548,       'HCOOH'],
    85: [33,   1,  'H(16O)2',                   0.995107,           32.997655,      'HO2'],
    86: [34,   1,  '(16O)',                     0.997628,           15.994915,      'O'],
    87: [36,   1,  '(14N)(16O)+',               0.993974,           29.997989,      'NOp'],
    88: [37,   1,  'H(16O)(79Br)',              0.5056,             95.921076,      'HOBr'],
    89: [37,   2,  'H(16O)(81Br)',              0.4919,             97.919027,      'HOBr'],
    90: [38,   1,  '(12C)2H4',                  0.9773,             28.0313,        'C2H4'],
    91: [38,   2,  '(12C)H2(13C)H2',            0.02196,            29.034655,      'C2H4'],
    92: [39,   1,  '(12C)H3(16O)H',             0.98593,            32.026215,      'CH3OH'],
    93: [40,   1,  '(12C)H3(79Br)',             0.5013,             93.941811,      'CH3Br'],
    94: [40,   2,  '(12C)H3(81Br)',             0.48766,            95.939764,      'CH3Br'],
    95: [41,   1,  '(12C)H3(12C)(14N)',         0.97482,            41.026549,      'CH3CN'],
    96: [42,   1,  '(12C)(19F)4',               0.9893,             87.993616,      'CF4'],
    116: [43,   1,  '(12C)4H2',                  0.955998,           50.01565,       'C4H2'],
    109: [44,   1,  'H(12C)3(14N)',              0.9646069,          51.01089903687, 'HC3N'],
    103: [45,   1,  'H2',                        0.999688,           2.01565,        'H2'],
    115: [45,   2,  'HD',                        0.00022997,         3.021825,       'H2'],
    97: [46,   1,  '(12C)(32S)',                0.939624,           43.971036,      'CS'],
    98: [46,   2,  '(12C)(34S)',                0.0416817,          45.966787,      'CS'],
    99: [46,   3,  '(13C)(32S)',                0.0105565,          44.974368,      'CS'],
    100: [46,   4,  '(12C)(33S)',                0.00741668,         44.970399,      'CS'],
    114: [47,   1,  '(32S)(16O)3',               0.9423964,          79.95682,       'SO3'],
    101: [1001,   1,  'H',                         None,               None,           'H'],
    102: [1002,   1,  'He',                        None,               None,           'He'],
    104: [1018,   1,  'Ar',                        None,               None,           'Ar'],

}

ISO_INDEX = {

    'id': 0,
    'iso_name': 1,
    'abundance': 2,
    'mass': 3,
    'mol_name': 4

}

#        M    I             id    iso_name                    abundance           mass        mol_name

ISO = {

    (1,   1): [1,  'H2(16O)',                   0.997317,           18.010565,      'H2O'],
    (1,   2): [2,  'H2(18O)',                   0.00199983,         20.014811,      'H2O'],
    (1,   3): [3,  'H2(17O)',                   0.000372,           19.01478,       'H2O'],
    (1,   4): [4,  'HD(16O)',                   0.00031069,         19.01674,       'H2O'],
    (1,   5): [5,  'HD(18O)',                   0.000000623,        21.020985,      'H2O'],
    (1,   6): [6,  'HD(17O)',                   0.000000116,        20.020956,      'H2O'],
    (2,   1): [7,  '(12C)(16O)2',               0.9842,             43.98983,       'CO2'],
    (2,   2): [8,  '(13C)(16O)2',               0.01106,            44.993185,      'CO2'],
    (2,   3): [9,  '(16O)(12C)(18O)',           0.0039471,          45.994076,      'CO2'],
    (2,   4): [10,  '(16O)(12C)(17O)',           0.000734,           44.994045,      'CO2'],
    (2,   5): [11,  '(16O)(13C)(18O)',           0.00004434,         46.997431,      'CO2'],
    (2,   6): [12,  '(16O)(13C)(17O)',           0.00000825,         45.9974,        'CO2'],
    (2,   7): [13,  '(12C)(18O)2',               0.0000039573,       47.998322,      'CO2'],
    (2,   8): [14,  '(17O)(12C)(18O)',           0.00000147,         46.998291,      'CO2'],
    (2,   0): [15,  '(13C)(18O)2',               0.000000044967,     49.001675,      'CO2'],
    (2,  11): [120,  '(18O)(13C)(17O)',           0.00000001654,      48.00165,       'CO2'],
    (2,   9): [121,  '(12C)(17O)2',               0.0000001368,       45.998262,      'CO2'],
    (3,   1): [16,  '(16O)3',                    0.992901,           47.984745,      'O3'],
    (3,   2): [17,  '(16O)(16O)(18O)',           0.00398194,         49.988991,      'O3'],
    (3,   3): [18,  '(16O)(18O)(16O)',           0.00199097,         49.988991,      'O3'],
    (3,   4): [19,  '(16O)(16O)(17O)',           0.00074,            48.98896,       'O3'],
    (3,   5): [20,  '(16O)(17O)(16O)',           0.00037,            48.98896,       'O3'],
    (4,   1): [21,  '(14N)2(16O)',               0.990333,           44.001062,      'N2O'],
    (4,   2): [22,  '(14N)(15N)(16O)',           0.0036409,          44.998096,      'N2O'],
    (4,   3): [23,  '(15N)(14N)(16O)',           0.0036409,          44.998096,      'N2O'],
    (4,   4): [24,  '(14N)2(18O)',               0.00198582,         46.005308,      'N2O'],
    (4,   5): [25,  '(14N)2(17O)',               0.000369,           45.005278,      'N2O'],
    (5,   1): [26,  '(12C)(16O)',                0.98654,            27.994915,      'CO'],
    (5,   2): [27,  '(13C)(16O)',                0.01108,            28.99827,       'CO'],
    (5,   3): [28,  '(12C)(18O)',                0.0019782,          29.999161,      'CO'],
    (5,   4): [29,  '(12C)(17O)',                0.000368,           28.99913,       'CO'],
    (5,   5): [30,  '(13C)(18O)',                0.00002222,         31.002516,      'CO'],
    (5,   6): [31,  '(13C)(17O)',                0.00000413,         30.002485,      'CO'],
    (6,   1): [32,  '(12C)H4',                   0.98827,            16.0313,        'CH4'],
    (6,   2): [33,  '(13C)H4',                   0.0111,             17.034655,      'CH4'],
    (6,   3): [34,  '(12C)H3D',                  0.00061575,         17.037475,      'CH4'],
    (6,   4): [35,  '(13C)H3D',                  0.0000049203,       18.04083,       'CH4'],
    (7,   1): [36,  '(16O)2',                    0.995262,           31.98983,       'O2'],
    (7,   2): [37,  '(16O)(18O)',                0.00399141,         33.994076,      'O2'],
    (7,   3): [38,  '(16O)(17O)',                0.000742,           32.994045,      'O2'],
    (8,   1): [39,  '(14N)(16O)',                0.993974,           29.997989,      'NO'],
    (8,   2): [40,  '(15N)(16O)',                0.0036543,          30.995023,      'NO'],
    (8,   3): [41,  '(14N)(18O)',                0.00199312,         32.002234,      'NO'],
    (9,   1): [42,  '(32S)(16O)2',               0.94568,            63.961901,      'SO2'],
    (9,   2): [43,  '(34S)(16O)2',               0.04195,            65.957695,      'SO2'],
    (10,   1): [44,  '(14N)(16O)2',               0.991616,           45.992904,      'NO2'],
    (11,   1): [45,  '(14N)H3',                   0.9958715,          17.026549,      'NH3'],
    (11,   2): [46,  '(15N)H3',                   0.0036613,          18.023583,      'NH3'],
    (12,   1): [47,  'H(14N)(16O)3',              0.98911,            62.995644,      'HNO3'],
    (12,   2): [117,  'H(15N)(16O)3',              0.003636,           63.99268,      'HNO3'],
    (13,   1): [48,  '(16O)H',                    0.997473,           17.00274,       'OH'],
    (13,   2): [49,  '(18O)H',                    0.00200014,         19.006986,      'OH'],
    (13,   3): [50,  '(16O)D',                    0.00015537,         18.008915,      'OH'],
    (14,   1): [51,  'H(19F)',                    0.99984425,         20.006229,      'HF'],
    (14,   2): [110,  'D(19F)',                    0.000115,           21.0125049978,  'HF'],
    (15,   1): [52,  'H(35Cl)',                   0.757587,           35.976678,      'HCl'],
    (15,   2): [53,  'H(37Cl)',                   0.242257,           37.973729,      'HCl'],
    (15,   3): [107,  'D(35Cl)',                   0.000118005,        36.9829544578,  'HCl'],
    (15,   4): [108,  'D(37Cl)',                   0.000037735,        38.9800043678,  'HCl'],
    (16,   1): [54,  'H(79Br)',                   0.50678,            79.92616,       'HBr'],
    (16,   2): [55,  'H(81Br)',                   0.49306,            81.924115,      'HBr'],
    (16,   3): [111,  'D(79Br)',                   0.0000582935,       80.9324388778,  'HBr'],
    (16,   4): [112,  'D(81Br)',                   0.0000567065,       82.9303923778,  'HBr'],
    (17,   1): [56,  'H(127I)',                   0.99984425,         127.912297,     'HI'],
    (17,   2): [113,  'D(127I)',                   0.000115,           128.918574778,  'HI'],
    (18,   1): [57,  '(35Cl)(16O)',               0.75591,            50.963768,      'ClO'],
    (18,   2): [58,  '(37Cl)(16O)',               0.24172,            52.960819,      'ClO'],
    (19,   1): [59,  '(16O)(12C)(32S)',           0.93739,            59.966986,      'OCS'],
    (19,   2): [60,  '(16O)(12C)(34S)',           0.04158,            61.96278,       'OCS'],
    (19,   3): [61,  '(16O)(13C)(32S)',           0.01053,            60.970341,      'OCS'],
    (19,   4): [62,  '(16O)(12C)(33S)',           0.01053,            60.966371,      'OCS'],
    (19,   5): [63,  '(18O)(12C)(32S)',           0.00188,            61.971231,      'OCS'],
    (20,   1): [64,  'H2(12C)(16O)',              0.98624,            30.010565,      'H2CO'],
    (20,   2): [65,  'H2(13C)(16O)',              0.01108,            31.01392,       'H2CO'],
    (20,   3): [66,  'H2(12C)(18O)',              0.0019776,          32.014811,      'H2CO'],
    (21,   1): [67,  'H(16O)(35Cl)',              0.75579,            51.971593,      'HOCl'],
    (21,   2): [68,  'H(16O)(37Cl)',              0.24168,            53.968644,      'HOCl'],
    (22,   1): [69,  '(14N)2',                    0.9926874,          28.006147,      'N2'],
    (22,   2): [118,  '(14N)(15N)',                0.0072535,          29.997989,      'N2'],
    (23,   1): [70,  'H(12C)(14N)',               0.98511,            27.010899,      'HCN'],
    (23,   2): [71,  'H(13C)(14N)',               0.01107,            28.014254,      'HCN'],
    (23,   3): [72,  'H(12C)(15N)',               0.0036217,          28.007933,      'HCN'],
    (24,   1): [73,  '(12C)H3(35Cl)',             0.74894,            49.992328,      'CH3Cl'],
    (24,   2): [74,  '(12C)H3(37Cl)',             0.23949,            51.989379,      'CH3Cl'],
    (25,   1): [75,  'H2(16O)2',                  0.994952,           34.00548,       'H2O2'],
    (26,   1): [76,  '(12C)2H2',                  0.9776,             26.01565,       'C2H2'],
    (26,   2): [77,  '(12C)(13C)H2',              0.02197,            27.019005,      'C2H2'],
    (26,   3): [105,  '(12C)2HD',                  0.00030455,         27.021825,      'C2H2'],
    (27,   1): [78,  '(12C)2H6',                  0.97699,            30.04695,       'C2H6'],
    (27,   2): [106,  '(12C)H3(13C)H3',            0.021952611,        31.050305,      'C2H6'],
    (28,   1): [79,  '(31P)H3',                   0.99953283,         33.997238,      'PH3'],
    (29,   1): [80,  '(12C)(16O)(19F)2',          0.98654,            65.991722,      'COF2'],
    (29,   2): [119,  '(13C)(16O)(19F)2',          0.0110834,          66.995083,      'COF2'],
    (31,   1): [81,  'H2(32S)',                   0.94988,            33.987721,      'H2S'],
    (31,   2): [82,  'H2(34S)',                   0.04214,            35.983515,      'H2S'],
    (31,   3): [83,  'H2(33S)',                   0.007498,           34.987105,      'H2S'],
    (32,   1): [84,  'H(12C)(16O)(16O)H',         0.983898,           46.00548,       'HCOOH'],
    (33,   1): [85,  'H(16O)2',                   0.995107,           32.997655,      'HO2'],
    (34,   1): [86,  '(16O)',                     0.997628,           15.994915,      'O'],
    (36,   1): [87,  '(14N)(16O)+',               0.993974,           29.997989,      'NOp'],
    (37,   1): [88,  'H(16O)(79Br)',              0.5056,             95.921076,      'HOBr'],
    (37,   2): [89,  'H(16O)(81Br)',              0.4919,             97.919027,      'HOBr'],
    (38,   1): [90,  '(12C)2H4',                  0.9773,             28.0313,        'C2H4'],
    (38,   2): [91,  '(12C)H2(13C)H2',            0.02196,            29.034655,      'C2H4'],
    (39,   1): [92,  '(12C)H3(16O)H',             0.98593,            32.026215,      'CH3OH'],
    (40,   1): [93,  '(12C)H3(79Br)',             0.5013,             93.941811,      'CH3Br'],
    (40,   2): [94,  '(12C)H3(81Br)',             0.48766,            95.939764,      'CH3Br'],
    (41,   1): [95,  '(12C)H3(12C)(14N)',         0.97482,            41.026549,      'CH3CN'],
    (42,   1): [96,  '(12C)(19F)4',               0.9893,             87.993616,      'CF4'],
    (43,   1): [116,  '(12C)4H2',                  0.955998,           50.01565,       'C4H2'],
    (44,   1): [109,  'H(12C)3(14N)',              0.9646069,          51.01089903687, 'HC3N'],
    (45,   1): [103,  'H2',                        0.999688,           2.01565,        'H2'],
    (45,   2): [115,  'HD',                        0.00022997,         3.021825,       'H2'],
    (46,   1): [97,  '(12C)(32S)',                0.939624,           43.971036,      'CS'],
    (46,   2): [98,  '(12C)(34S)',                0.0416817,          45.966787,      'CS'],
    (46,   3): [99,  '(13C)(32S)',                0.0105565,          44.974368,      'CS'],
    (46,   4): [100,  '(12C)(33S)',                0.00741668,         44.970399,      'CS'],
    (47,   1): [114,  '(32S)(16O)3',               0.9423964,          79.95682,       'SO3'],
    (1001,   1): [101,  'H',                         None,               None,           'H'],
    (1002,   1): [102,  'He',                        None,               None,           'He'],
    (1018,   1): [104,  'Ar',                        None,               None,           'Ar'],

}


def print_iso():
    print('The dictionary \"ISO\" contains information on isotopologues in HITRAN\n')
    print('   M    I          id                  iso_name   abundance      mass        mol_name')
    for i in ISO:
        ab = ISO[i][ISO_INDEX['abundance']]
        ma = ISO[i][ISO_INDEX['mass']]
        ab = ab if ab else -1
        ma = ma if ma else -1
        print('%4i %4i     : %5i %25s %10f %10f %15s' % (
            i[0], i[1], ISO[i][ISO_INDEX['id']], ISO[i][ISO_INDEX['iso_name']], ab, ma, ISO[i][ISO_INDEX['mol_name']]))


def print_iso_id():
    print('The dictionary \"ISO_ID\" contains information on \"global\" IDs of isotopologues in HITRAN\n')
    print('   id            M    I                    iso_name       abundance       mass        mol_name')
    for i in ISO_ID:
        ab = ISO_ID[i][ISO_ID_INDEX['abundance']]
        ma = ISO_ID[i][ISO_ID_INDEX['mass']]
        ab = ab if ab else -1
        ma = ma if ma else -1
        print('%5i     :   %4i %4i   %25s %15.10f %10f %15s' % (
            i, ISO_ID[i][ISO_ID_INDEX['M']], ISO_ID[i][ISO_ID_INDEX['I']], ISO_ID[i][ISO_ID_INDEX['iso_name']], ab, ma, ISO_ID[i][ISO_ID_INDEX['mol_name']]))


profiles = 'profiles'


def print_profiles():
    print('Profiles available:')
    print('  HT        : PROFILE_HT')
    print('  SDRautian : PROFILE_SDRAUTIAN')
    print('  Rautian   : PROFILE_RAUTIAN')
    print('  SDVoigt   : PROFILE_SDVOIGT')
    print('  Voigt     : PROFILE_VOIGT')
    print('  Lorentz   : PROFILE_LORENTZ')
    print('  Doppler   : PROFILE_DOPPLER')


slit_functions = 'slit_functions'


def print_slit_functions():
    print('  RECTANGULAR : SLIT_RECTANGULAR')
    print('  TRIANGULAR  : SLIT_TRIANGULAR')
    print('  GAUSSIAN    : SLIT_GAUSSIAN')
    print('  DIFFRACTION : SLIT_DIFFRACTION')
    print('  MICHELSON   : SLIT_MICHELSON')
    print('  DISPERSION/LORENTZ : SLIT_DISPERSION')


tutorial = 'tutorial'
units = 'units'
index = 'index'
data = 'data'
spectra = 'spectra'
plotting = 'plotting'
python = 'python'

python_tutorial_text = \
    """
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


data_tutorial_text = \
    """

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


spectra_tutorial_text = \
    """

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


plotting_tutorial_text = \
    """

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
        print('--------------------------------------------------------------')
        print('Hello, this is an interactive help system of HITRANonline API.')
        print('--------------------------------------------------------------')
        print('Run getHelp(.) with one of the following arguments:')
        print('    tutorial  -  interactive tutorials on HAPI')
        print('    units     -  units used in calculations')
        print('    index     -  index of available HAPI functions')
    elif arg == 'tutorial':
        print('-----------------------------------')
        print('This is a tutorial section of help.')
        print('-----------------------------------')
        print('Please choose the subject of tutorial:')
        print('    data      -  downloading the data and working with it')
        print('    spectra   -  calculating spectral functions')
        print('    plotting  -  visualizing data with matplotlib')
        print('    python    -  Python quick start guide')
    elif arg == 'python':
        print_python_tutorial()
    elif arg == 'data':
        print_data_tutorial()
    elif arg == 'spectra':
        print_spectra_tutorial()
    elif arg == 'plotting':
        print_plotting_tutorial()
    elif arg == 'index':
        print('------------------------------')
        print('FETCHING DATA:')
        print('------------------------------')
        print('  fetch')
        print('  fetch_by_ids')
        print('')
        print('------------------------------')
        print('WORKING WITH DATA:')
        print('------------------------------')
        print('  db_begin')
        print('  db_commit')
        print('  tableList')
        print('  describe')
        print('  select')
        print('  sort')
        print('  extractColumns')
        print('  getColumn')
        print('  getColumns')
        print('  dropTable')
        print('')
        print('------------------------------')
        print('CALCULATING SPECTRA:')
        print('------------------------------')
        print('  profiles')
        print('  partitionSum')
        print('  absorptionCoefficient_HT')
        print('  absorptionCoefficient_Voigt')
        print('  absorptionCoefficient_SDVoigt')
        print('  absorptionCoefficient_Lorentz')
        print('  absorptionCoefficient_Doppler')
        print('  transmittanceSpectrum')
        print('  absorptionSpectrum')
        print('  radianceSpectrum')
        print('')
        print('------------------------------')
        print('CONVOLVING SPECTRA:')
        print('------------------------------')
        print('  convolveSpectrum')
        print('  slit_functions')
        print('')
        print('------------------------------')
        print('INFO ON ISOTOPOLOGUES:')
        print('------------------------------')
        print('  ISO_ID')
        print('  abundance')
        print('  molecularMass')
        print('  moleculeName')
        print('  isotopologueName')
        print('')
        print('------------------------------')
        print('MISCELLANEOUS:')
        print('------------------------------')
        print('  getStickXY')
        print('  read_hotw')
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
    return ISO[(M, I)][ISO_INDEX['abundance']]

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
    return ISO[(M, I)][ISO_INDEX['mass']]

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
    return ISO[(M, 1)][ISO_INDEX['mol_name']]

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
    return ISO[(M, I)][ISO_INDEX['iso_name']]

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
    LOCAL_TABLE_CACHE[TableName]['header']['comment'] = Comment


def fetch_by_ids(TableName, iso_id_list, numin, numax, ParameterGroups=[], Parameters=[]):
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
    queryHITRAN(TableName, iso_id_list, numin, numax,
                pargroups=ParameterGroups, params=Parameters)
    iso_names = [ISO_ID[i][ISO_ID_INDEX['iso_name']] for i in iso_id_list]
    Comment = 'Contains lines for '+','.join(iso_names)
    Comment += ('\n in %.3f-%.3f wavenumber range' % (numin, numax))
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
    queryHITRAN(TableName, [ISO[(M, I)][ISO_INDEX['id']]], numin, numax,
                pargroups=ParameterGroups, params=Parameters)
    iso_name = ISO[(M, I)][ISO_INDEX['iso_name']]
    Comment = 'Contains lines for '+iso_name
    Comment += ('\n in %.3f-%.3f wavenumber range' % (numin, numax))
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
    for I in range(2, npt+1):
        if A[I-1] >= aa:
            if I < 3 or I == npt:
                J = I
                if I < 3:
                    J = 3
                if I == npt:
                    J = npt
                J = J-1   # zero index correction
                A0D1 = A[J-2]-A[J-1]
                if A0D1 == 0.0:
                    A0D1 = 0.0001
                A0D2 = A[J-2]-A[J]
                if A0D2 == 0.0:
                    A0D2 = 0.0000
                A1D1 = A[J-1]-A[J-2]
                if A1D1 == 0.0:
                    A1D1 = 0.0001
                A1D2 = A[J-1]-A[J]
                if A1D2 == 0.0:
                    A1D2 = 0.0001
                A2D1 = A[J]-A[J-2]
                if A2D1 == 0.0:
                    A2D1 = 0.0001
                A2D2 = A[J]-A[J-1]
                if A2D2 == 0.0:
                    A2D2 = 0.0001

                A0 = (aa-A[J-1])*(aa-A[J])/(A0D1*A0D2)
                A1 = (aa-A[J-2])*(aa-A[J])/(A1D1*A1D2)
                A2 = (aa-A[J-2])*(aa-A[J-1])/(A2D1*A2D2)

                bb = A0*B[J-2] + A1*B[J-1] + A2*B[J]

            else:
                J = I
                J = J-1   # zero index correction
                A0D1 = A[J-2]-A[J-1]
                if A0D1 == 0.0:
                    A0D1 = 0.0001
                A0D2 = A[J-2]-A[J]
                if A0D2 == 0.0:
                    A0D2 = 0.0001
                A0D3 = (A[J-2]-A[J+1])
                if A0D3 == 0.0:
                    A0D3 = 0.0001
                A1D1 = A[J-1]-A[J-2]
                if A1D1 == 0.0:
                    A1D1 = 0.0001
                A1D2 = A[J-1]-A[J]
                if A1D2 == 0.0:
                    A1D2 = 0.0001
                A1D3 = A[J-1]-A[J+1]
                if A1D3 == 0.0:
                    A1D3 = 0.0001

                A2D1 = A[J]-A[J-2]
                if A2D1 == 0.0:
                    A2D1 = 0.0001
                A2D2 = A[J]-A[J-1]
                if A2D2 == 0.0:
                    A2D2 = 0.0001
                A2D3 = A[J]-A[J+1]
                if A2D3 == 0.0:
                    A2D3 = 0.0001

                A3D1 = A[J+1]-A[J-2]
                if A3D1 == 0.0:
                    A3D1 = 0.0001
                A3D2 = A[J+1]-A[J-1]
                if A3D2 == 0.0:
                    A3D2 = 0.0001
                A3D3 = A[J+1]-A[J]
                if A3D3 == 0.0:
                    A3D3 = 0.0001

                A0 = (aa-A[J-1])*(aa-A[J])*(aa-A[J+1])
                A0 = A0/(A0D1*A0D2*A0D3)
                A1 = (aa-A[J-2])*(aa-A[J])*(aa-A[J+1])
                A1 = A1/(A1D1*A1D2*A1D3)
                A2 = (aa-A[J-2])*(aa-A[J-1])*(aa-A[J+1])
                A2 = A2/(A2D1*A2D2*A2D3)
                A3 = (aa-A[J-2])*(aa-A[J-1])*(aa-A[J])
                A3 = A3/(A3D1*A3D2*A3D3)

                bb = A0*B[J-2] + A1*B[J-1] + A2*B[J] + A3*B[J+1]

            break

    return bb


#  --------------- ISOTOPOLOGUE HASH ----------------------

TIPS_ISO_HASH = {}

#  --------------- STATISTICAL WEIGHT HASH ----------------------

TIPS_GSI_HASH = {}

#  --------------- INTERPOLATION NODES ----------------------

Tdat = __FloatType__([60.,  85., 110., 135., 160., 185., 210., 235.,
                      260., 285., 310., 335., 360., 385., 410., 435., 460., 485.,
                      510., 535., 560., 585., 610., 635., 660., 685., 710., 735.,
                      760., 785., 810., 835., 860., 885., 910., 935., 960., 985.,
                      1010., 1035., 1060., 1085., 1110., 1135., 1160., 1185., 1210., 1235.,
                      1260., 1285., 1310., 1335., 1360., 1385., 1410., 1435., 1460., 1485.,
                      1510., 1535., 1560., 1585., 1610., 1635., 1660., 1685., 1710., 1735.,
                      1760., 1785., 1810., 1835., 1860., 1885., 1910., 1935., 1960., 1985.,
                      2010., 2035., 2060., 2085., 2110., 2135., 2160., 2185., 2210., 2235.,
                      2260., 2285., 2310., 2335., 2360., 2385., 2410., 2435., 2460., 2485.,
                      2510., 2535., 2560., 2585., 2610., 2635., 2660., 2685., 2710., 2735.,
                      2760., 2785., 2810., 2835., 2860., 2885., 2910., 2935., 2960., 2985.,
                      3010.])

TIPS_NPT = len(Tdat)


# REMARK
# float32 gives exactly the same results as fortran TIPS, because
# all constants in the fortran code given as xx.xxE+-XX, i.e.
# in single precision. By this fact all unsignificant figures
# over single precision are filled with digital garbage


#  --------------- H2O 161: M = 1, I = 1 ---------------------
M = 1
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.16824E+02, 0.27771E+02, 0.40408E+02,
                                 0.54549E+02, 0.70054E+02, 0.86817E+02, 0.10475E+03, 0.12380E+03,
                                 0.14391E+03, 0.16503E+03, 0.18714E+03, 0.21021E+03, 0.23425E+03,
                                 0.25924E+03, 0.28518E+03, 0.31209E+03, 0.33997E+03, 0.36883E+03,
                                 0.39870E+03, 0.42959E+03, 0.46152E+03, 0.49452E+03, 0.52860E+03,
                                 0.56380E+03, 0.60015E+03, 0.63766E+03, 0.67637E+03, 0.71631E+03,
                                 0.75750E+03, 0.79999E+03, 0.84380E+03, 0.88897E+03, 0.93553E+03,
                                 0.98353E+03, 0.10330E+04, 0.10840E+04, 0.11365E+04, 0.11906E+04,
                                 0.12463E+04, 0.13037E+04, 0.13628E+04, 0.14237E+04, 0.14863E+04,
                                 0.15509E+04, 0.16173E+04, 0.16856E+04, 0.17559E+04, 0.18283E+04,
                                 0.19028E+04, 0.19793E+04, 0.20581E+04, 0.21391E+04, 0.22224E+04,
                                 0.23080E+04, 0.24067E+04, 0.24975E+04, 0.25908E+04, 0.26867E+04,
                                 0.27853E+04, 0.28865E+04, 0.29904E+04, 0.30972E+04, 0.32068E+04,
                                 0.33194E+04, 0.34349E+04, 0.35535E+04, 0.36752E+04, 0.38001E+04,
                                 0.39282E+04, 0.40597E+04, 0.41945E+04, 0.43327E+04, 0.44745E+04,
                                 0.46199E+04, 0.47688E+04, 0.49215E+04, 0.50780E+04, 0.52384E+04,
                                 0.54027E+04, 0.55710E+04, 0.57434E+04, 0.59200E+04, 0.61008E+04,
                                 0.62859E+04, 0.64754E+04, 0.66693E+04, 0.68679E+04, 0.70710E+04,
                                 0.72788E+04, 0.74915E+04, 0.77090E+04, 0.79315E+04, 0.81590E+04,
                                 0.83917E+04, 0.86296E+04, 0.88728E+04, 0.91214E+04, 0.93755E+04,
                                 0.96351E+04, 0.99005E+04, 0.10171E+05, 0.10448E+05, 0.10731E+05,
                                 0.11020E+05, 0.11315E+05, 0.11617E+05, 0.11924E+05, 0.12238E+05,
                                 0.12559E+05, 0.12886E+05, 0.13220E+05, 0.13561E+05, 0.13909E+05,
                                 0.14263E+05, 0.14625E+05, 0.14995E+05, 0.15371E+05, 0.15755E+05,
                                 0.16147E+05])


#  --------------- H2O 181: M = 1, I = 2 ---------------------
M = 1
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.15960E+02, 0.26999E+02, 0.39743E+02,
                                 0.54003E+02, 0.69639E+02, 0.86543E+02, 0.10463E+03, 0.12384E+03,
                                 0.14412E+03, 0.16542E+03, 0.18773E+03, 0.21103E+03, 0.23531E+03,
                                 0.26057E+03, 0.28681E+03, 0.31406E+03, 0.34226E+03, 0.37130E+03,
                                 0.40135E+03, 0.43243E+03, 0.46456E+03, 0.49777E+03, 0.53206E+03,
                                 0.56748E+03, 0.60405E+03, 0.64179E+03, 0.68074E+03, 0.72093E+03,
                                 0.76238E+03, 0.80513E+03, 0.84922E+03, 0.89467E+03, 0.94152E+03,
                                 0.98982E+03, 0.10396E+04, 0.10909E+04, 0.11437E+04, 0.11982E+04,
                                 0.12543E+04, 0.13120E+04, 0.13715E+04, 0.14328E+04, 0.14959E+04,
                                 0.15608E+04, 0.16276E+04, 0.16964E+04, 0.17672E+04, 0.18401E+04,
                                 0.19151E+04, 0.19922E+04, 0.20715E+04, 0.21531E+04, 0.22370E+04,
                                 0.23232E+04, 0.24118E+04, 0.25030E+04, 0.25967E+04, 0.26929E+04,
                                 0.27918E+04, 0.28934E+04, 0.29978E+04, 0.31050E+04, 0.32151E+04,
                                 0.33281E+04, 0.34441E+04, 0.35632E+04, 0.36854E+04, 0.38108E+04,
                                 0.39395E+04, 0.40715E+04, 0.42070E+04, 0.43459E+04, 0.44883E+04,
                                 0.46343E+04, 0.47840E+04, 0.49374E+04, 0.50946E+04, 0.52558E+04,
                                 0.54209E+04, 0.55900E+04, 0.57632E+04, 0.59407E+04, 0.61224E+04,
                                 0.63084E+04, 0.64988E+04, 0.66938E+04, 0.68933E+04, 0.70975E+04,
                                 0.73064E+04, 0.75202E+04, 0.77389E+04, 0.79625E+04, 0.81913E+04,
                                 0.84252E+04, 0.86644E+04, 0.89089E+04, 0.91588E+04, 0.94143E+04,
                                 0.96754E+04, 0.99422E+04, 0.10215E+05, 0.10493E+05, 0.10778E+05,
                                 0.11068E+05, 0.11365E+05, 0.11668E+05, 0.11977E+05, 0.12293E+05,
                                 0.12616E+05, 0.12945E+05, 0.13281E+05, 0.13624E+05, 0.13973E+05,
                                 0.14330E+05, 0.14694E+05, 0.15066E+05, 0.15445E+05, 0.15831E+05,
                                 0.16225E+05])


#  --------------- H2O 171: M = 1, I = 3 ---------------------
M = 1
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.95371E+02, 0.16134E+03, 0.23750E+03,
                                 0.32273E+03, 0.41617E+03, 0.51722E+03, 0.62540E+03, 0.74036E+03,
                                 0.86185E+03, 0.98970E+03, 0.11238E+04, 0.12642E+04, 0.14097E+04,
                                 0.15599E+04, 0.17159E+04, 0.18777E+04, 0.20453E+04, 0.22188E+04,
                                 0.23983E+04, 0.25840E+04, 0.27760E+04, 0.29743E+04, 0.31792E+04,
                                 0.33907E+04, 0.36091E+04, 0.38346E+04, 0.40672E+04, 0.43072E+04,
                                 0.45547E+04, 0.48100E+04, 0.50732E+04, 0.53446E+04, 0.56244E+04,
                                 0.59128E+04, 0.62100E+04, 0.65162E+04, 0.68317E+04, 0.71567E+04,
                                 0.74915E+04, 0.78363E+04, 0.81914E+04, 0.85571E+04, 0.89335E+04,
                                 0.93211E+04, 0.97200E+04, 0.10131E+05, 0.10553E+05, 0.10988E+05,
                                 0.11435E+05, 0.11895E+05, 0.12368E+05, 0.12855E+05, 0.13356E+05,
                                 0.13870E+05, 0.14399E+05, 0.14943E+05, 0.15502E+05, 0.16076E+05,
                                 0.16666E+05, 0.17272E+05, 0.17895E+05, 0.18534E+05, 0.19191E+05,
                                 0.19865E+05, 0.20557E+05, 0.21267E+05, 0.21996E+05, 0.22744E+05,
                                 0.23512E+05, 0.24299E+05, 0.25106E+05, 0.25935E+05, 0.26784E+05,
                                 0.27655E+05, 0.28547E+05, 0.29462E+05, 0.30400E+05, 0.31361E+05,
                                 0.32345E+05, 0.33353E+05, 0.34386E+05, 0.35444E+05, 0.36527E+05,
                                 0.37637E+05, 0.38772E+05, 0.39934E+05, 0.41124E+05, 0.42341E+05,
                                 0.43587E+05, 0.44861E+05, 0.46165E+05, 0.47498E+05, 0.48862E+05,
                                 0.50256E+05, 0.51682E+05, 0.53139E+05, 0.54629E+05, 0.56152E+05,
                                 0.57708E+05, 0.59299E+05, 0.60923E+05, 0.62583E+05, 0.64279E+05,
                                 0.66011E+05, 0.67779E+05, 0.69585E+05, 0.71429E+05, 0.73312E+05,
                                 0.75234E+05, 0.77195E+05, 0.79197E+05, 0.81240E+05, 0.83325E+05,
                                 0.85452E+05, 0.87622E+05, 0.89835E+05, 0.92093E+05, 0.94395E+05,
                                 0.96743E+05])


#  --------------- H2O 162: M = 1, I = 4 ---------------------
M = 1
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.75792E+02, 0.12986E+03, 0.19244E+03,
                                 0.26253E+03, 0.33942E+03, 0.42259E+03, 0.51161E+03, 0.60619E+03,
                                 0.70609E+03, 0.81117E+03, 0.92132E+03, 0.10365E+04, 0.11567E+04,
                                 0.12820E+04, 0.14124E+04, 0.15481E+04, 0.16891E+04, 0.18355E+04,
                                 0.19876E+04, 0.21455E+04, 0.23092E+04, 0.24791E+04, 0.26551E+04,
                                 0.28376E+04, 0.30268E+04, 0.32258E+04, 0.34288E+04, 0.36392E+04,
                                 0.38571E+04, 0.40828E+04, 0.43165E+04, 0.45584E+04, 0.48089E+04,
                                 0.50681E+04, 0.53363E+04, 0.56139E+04, 0.59009E+04, 0.61979E+04,
                                 0.65049E+04, 0.68224E+04, 0.71506E+04, 0.74898E+04, 0.78403E+04,
                                 0.82024E+04, 0.85765E+04, 0.89628E+04, 0.93618E+04, 0.97736E+04,
                                 0.10199E+05, 0.10637E+05, 0.11090E+05, 0.11557E+05, 0.12039E+05,
                                 0.12535E+05, 0.13047E+05, 0.13575E+05, 0.14119E+05, 0.14679E+05,
                                 0.15257E+05, 0.15851E+05, 0.16464E+05, 0.17094E+05, 0.17743E+05,
                                 0.18411E+05, 0.19098E+05, 0.19805E+05, 0.20532E+05, 0.21280E+05,
                                 0.22049E+05, 0.22840E+05, 0.23652E+05, 0.24487E+05, 0.25345E+05,
                                 0.26227E+05, 0.27132E+05, 0.28062E+05, 0.29016E+05, 0.29997E+05,
                                 0.31002E+05, 0.32035E+05, 0.33094E+05, 0.34180E+05, 0.35295E+05,
                                 0.36438E+05, 0.37610E+05, 0.38812E+05, 0.40044E+05, 0.41306E+05,
                                 0.42600E+05, 0.43926E+05, 0.45284E+05, 0.46675E+05, 0.48100E+05,
                                 0.49559E+05, 0.51053E+05, 0.52583E+05, 0.54148E+05, 0.55750E+05,
                                 0.57390E+05, 0.59067E+05, 0.60783E+05, 0.62539E+05, 0.64334E+05,
                                 0.66170E+05, 0.68047E+05, 0.69967E+05, 0.71929E+05, 0.73934E+05,
                                 0.75983E+05, 0.78078E+05, 0.80217E+05, 0.82403E+05, 0.84636E+05,
                                 0.86917E+05, 0.89246E+05, 0.91625E+05, 0.94053E+05, 0.96533E+05,
                                 0.99064E+05])


#  --------------- H2O 182: M = 1, I = 5 ---------------------
M = 1
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.82770E+02, 0.13749E+03, 0.20083E+03,
                                 0.27176E+03, 0.34955E+03, 0.43370E+03, 0.52376E+03, 0.61944E+03,
                                 0.72050E+03, 0.82679E+03, 0.93821E+03, 0.10547E+04, 0.11763E+04,
                                 0.13031E+04, 0.14350E+04, 0.15723E+04, 0.17150E+04, 0.18633E+04,
                                 0.20172E+04, 0.21770E+04, 0.23429E+04, 0.25149E+04, 0.26934E+04,
                                 0.28784E+04, 0.30702E+04, 0.32690E+04, 0.34750E+04, 0.36885E+04,
                                 0.39096E+04, 0.41386E+04, 0.43758E+04, 0.46213E+04, 0.48755E+04,
                                 0.51386E+04, 0.54109E+04, 0.56927E+04, 0.59841E+04, 0.62856E+04,
                                 0.65973E+04, 0.69197E+04, 0.72529E+04, 0.75973E+04, 0.79533E+04,
                                 0.83210E+04, 0.87009E+04, 0.90933E+04, 0.94985E+04, 0.99168E+04,
                                 0.10348E+05, 0.10794E+05, 0.11254E+05, 0.11728E+05, 0.12217E+05,
                                 0.12722E+05, 0.13242E+05, 0.13778E+05, 0.14331E+05, 0.14900E+05,
                                 0.15486E+05, 0.16091E+05, 0.16713E+05, 0.17353E+05, 0.18012E+05,
                                 0.18691E+05, 0.19389E+05, 0.20108E+05, 0.20847E+05, 0.21607E+05,
                                 0.22388E+05, 0.23191E+05, 0.24017E+05, 0.24866E+05, 0.25738E+05,
                                 0.26633E+05, 0.27553E+05, 0.28498E+05, 0.29468E+05, 0.30464E+05,
                                 0.31486E+05, 0.32536E+05, 0.33612E+05, 0.34716E+05, 0.35849E+05,
                                 0.37011E+05, 0.38202E+05, 0.39424E+05, 0.40676E+05, 0.41959E+05,
                                 0.43274E+05, 0.44622E+05, 0.46002E+05, 0.47416E+05, 0.48864E+05,
                                 0.50348E+05, 0.51866E+05, 0.53421E+05, 0.55012E+05, 0.56640E+05,
                                 0.58307E+05, 0.60012E+05, 0.61757E+05, 0.63541E+05, 0.65366E+05,
                                 0.67233E+05, 0.69141E+05, 0.71092E+05, 0.73087E+05, 0.75125E+05,
                                 0.77209E+05, 0.79338E+05, 0.81513E+05, 0.83736E+05, 0.86006E+05,
                                 0.88324E+05, 0.90693E+05, 0.93111E+05, 0.95580E+05, 0.98100E+05,
                                 0.10067E+06])


#  --------------- H2O 172: M = 1, I = 6 ---------------------
M = 1
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(36.)
TIPS_ISO_HASH[(M, I)] = float32([0.49379E+03, 0.82021E+03, 0.11980E+04,
                                 0.16211E+04, 0.20851E+04, 0.25870E+04, 0.31242E+04, 0.36949E+04,
                                 0.42977E+04, 0.49317E+04, 0.55963E+04, 0.62911E+04, 0.70164E+04,
                                 0.77722E+04, 0.85591E+04, 0.93777E+04, 0.10228E+05, 0.11112E+05,
                                 0.12030E+05, 0.12983E+05, 0.13971E+05, 0.14997E+05, 0.16061E+05,
                                 0.17163E+05, 0.18306E+05, 0.19491E+05, 0.20719E+05, 0.21991E+05,
                                 0.23309E+05, 0.24673E+05, 0.26086E+05, 0.27549E+05, 0.29064E+05,
                                 0.30631E+05, 0.32254E+05, 0.33932E+05, 0.35669E+05, 0.37464E+05,
                                 0.39321E+05, 0.41242E+05, 0.43227E+05, 0.45279E+05, 0.47399E+05,
                                 0.49589E+05, 0.51852E+05, 0.54189E+05, 0.56602E+05, 0.59094E+05,
                                 0.61666E+05, 0.64320E+05, 0.67058E+05, 0.69883E+05, 0.72796E+05,
                                 0.75801E+05, 0.78899E+05, 0.82092E+05, 0.85382E+05, 0.88773E+05,
                                 0.92266E+05, 0.95863E+05, 0.99568E+05, 0.10338E+06, 0.10731E+06,
                                 0.11135E+06, 0.11551E+06, 0.11979E+06, 0.12419E+06, 0.12871E+06,
                                 0.13337E+06, 0.13815E+06, 0.14307E+06, 0.14812E+06, 0.15331E+06,
                                 0.15865E+06, 0.16412E+06, 0.16975E+06, 0.17553E+06, 0.18146E+06,
                                 0.18754E+06, 0.19379E+06, 0.20020E+06, 0.20678E+06, 0.21352E+06,
                                 0.22044E+06, 0.22753E+06, 0.23480E+06, 0.24226E+06, 0.24990E+06,
                                 0.25773E+06, 0.26575E+06, 0.27397E+06, 0.28239E+06, 0.29102E+06,
                                 0.29985E+06, 0.30889E+06, 0.31814E+06, 0.32762E+06, 0.33731E+06,
                                 0.34724E+06, 0.35739E+06, 0.36777E+06, 0.37840E+06, 0.38926E+06,
                                 0.40038E+06, 0.41174E+06, 0.42335E+06, 0.43523E+06, 0.44737E+06,
                                 0.45977E+06, 0.47245E+06, 0.48540E+06, 0.49863E+06, 0.51214E+06,
                                 0.52595E+06, 0.54005E+06, 0.55444E+06, 0.56914E+06, 0.58415E+06,
                                 0.59947E+06])


#  --------------- CO2 626: M = 2, I = 1 ---------------------
M = 2
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.53642E+02, 0.75947E+02, 0.98292E+02,
                                 0.12078E+03, 0.14364E+03, 0.16714E+03, 0.19160E+03, 0.21731E+03,
                                 0.24454E+03, 0.27355E+03, 0.30456E+03, 0.33778E+03, 0.37343E+03,
                                 0.41170E+03, 0.45280E+03, 0.49692E+03, 0.54427E+03, 0.59505E+03,
                                 0.64948E+03, 0.70779E+03, 0.77019E+03, 0.83693E+03, 0.90825E+03,
                                 0.98440E+03, 0.10656E+04, 0.11522E+04, 0.12445E+04, 0.13427E+04,
                                 0.14471E+04, 0.15580E+04, 0.16759E+04, 0.18009E+04, 0.19334E+04,
                                 0.20739E+04, 0.22225E+04, 0.23798E+04, 0.25462E+04, 0.27219E+04,
                                 0.29074E+04, 0.31032E+04, 0.33097E+04, 0.35272E+04, 0.37564E+04,
                                 0.39976E+04, 0.42514E+04, 0.45181E+04, 0.47985E+04, 0.50929E+04,
                                 0.54019E+04, 0.57260E+04, 0.60659E+04, 0.64221E+04, 0.67952E+04,
                                 0.71859E+04, 0.75946E+04, 0.80222E+04, 0.84691E+04, 0.89362E+04,
                                 0.94241E+04, 0.99335E+04, 0.10465E+05, 0.11020E+05, 0.11598E+05,
                                 0.12201E+05, 0.12828E+05, 0.13482E+05, 0.14163E+05, 0.14872E+05,
                                 0.15609E+05, 0.16376E+05, 0.17173E+05, 0.18001E+05, 0.18861E+05,
                                 0.19754E+05, 0.20682E+05, 0.21644E+05, 0.22643E+05, 0.23678E+05,
                                 0.24752E+05, 0.25865E+05, 0.27018E+05, 0.28212E+05, 0.29449E+05,
                                 0.30730E+05, 0.32055E+05, 0.33426E+05, 0.34845E+05, 0.36312E+05,
                                 0.37828E+05, 0.39395E+05, 0.41015E+05, 0.42688E+05, 0.44416E+05,
                                 0.46199E+05, 0.48041E+05, 0.49942E+05, 0.51902E+05, 0.53925E+05,
                                 0.56011E+05, 0.58162E+05, 0.60379E+05, 0.62664E+05, 0.65019E+05,
                                 0.67444E+05, 0.69942E+05, 0.72515E+05, 0.75163E+05, 0.77890E+05,
                                 0.80695E+05, 0.83582E+05, 0.86551E+05, 0.89605E+05, 0.92746E+05,
                                 0.95975E+05, 0.99294E+05, 0.10271E+06, 0.10621E+06, 0.10981E+06,
                                 0.11351E+06])


#  --------------- CO2 636: M = 2, I = 2 ---------------------
M = 2
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.10728E+03, 0.15189E+03, 0.19659E+03,
                                 0.24164E+03, 0.28753E+03, 0.33486E+03, 0.38429E+03, 0.43643E+03,
                                 0.49184E+03, 0.55104E+03, 0.61449E+03, 0.68263E+03, 0.75589E+03,
                                 0.83468E+03, 0.91943E+03, 0.10106E+04, 0.11085E+04, 0.12137E+04,
                                 0.13266E+04, 0.14477E+04, 0.15774E+04, 0.17163E+04, 0.18649E+04,
                                 0.20237E+04, 0.21933E+04, 0.23743E+04, 0.25673E+04, 0.27729E+04,
                                 0.29917E+04, 0.32245E+04, 0.34718E+04, 0.37345E+04, 0.40132E+04,
                                 0.43087E+04, 0.46218E+04, 0.49533E+04, 0.53041E+04, 0.56749E+04,
                                 0.60668E+04, 0.64805E+04, 0.69171E+04, 0.73774E+04, 0.78626E+04,
                                 0.83736E+04, 0.89114E+04, 0.94772E+04, 0.10072E+05, 0.10697E+05,
                                 0.11353E+05, 0.12042E+05, 0.12765E+05, 0.13523E+05, 0.14317E+05,
                                 0.15148E+05, 0.16019E+05, 0.16930E+05, 0.17883E+05, 0.18879E+05,
                                 0.19920E+05, 0.21008E+05, 0.22143E+05, 0.23328E+05, 0.24563E+05,
                                 0.25852E+05, 0.27195E+05, 0.28594E+05, 0.30051E+05, 0.31568E+05,
                                 0.33146E+05, 0.34788E+05, 0.36496E+05, 0.38271E+05, 0.40115E+05,
                                 0.42031E+05, 0.44021E+05, 0.46086E+05, 0.48230E+05, 0.50453E+05,
                                 0.52759E+05, 0.55150E+05, 0.57628E+05, 0.60195E+05, 0.62854E+05,
                                 0.65608E+05, 0.68459E+05, 0.71409E+05, 0.74461E+05, 0.77618E+05,
                                 0.80883E+05, 0.84258E+05, 0.87746E+05, 0.91350E+05, 0.95073E+05,
                                 0.98918E+05, 0.10289E+06, 0.10698E+06, 0.11121E+06, 0.11558E+06,
                                 0.12008E+06, 0.12472E+06, 0.12950E+06, 0.13443E+06, 0.13952E+06,
                                 0.14475E+06, 0.15015E+06, 0.15571E+06, 0.16143E+06, 0.16732E+06,
                                 0.17338E+06, 0.17962E+06, 0.18604E+06, 0.19264E+06, 0.19943E+06,
                                 0.20642E+06, 0.21360E+06, 0.22098E+06, 0.22856E+06, 0.23636E+06,
                                 0.24436E+06])


#  --------------- CO2 628: M = 2, I = 3 ---------------------
M = 2
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.11368E+03, 0.16096E+03, 0.20833E+03,
                                 0.25603E+03, 0.30452E+03, 0.35442E+03, 0.40640E+03, 0.46110E+03,
                                 0.51910E+03, 0.58093E+03, 0.64709E+03, 0.71804E+03, 0.79422E+03,
                                 0.87607E+03, 0.96402E+03, 0.10585E+04, 0.11600E+04, 0.12689E+04,
                                 0.13857E+04, 0.15108E+04, 0.16449E+04, 0.17883E+04, 0.19416E+04,
                                 0.21054E+04, 0.22803E+04, 0.24668E+04, 0.26655E+04, 0.28770E+04,
                                 0.31021E+04, 0.33414E+04, 0.35956E+04, 0.38654E+04, 0.41516E+04,
                                 0.44549E+04, 0.47761E+04, 0.51160E+04, 0.54755E+04, 0.58555E+04,
                                 0.62568E+04, 0.66804E+04, 0.71273E+04, 0.75982E+04, 0.80944E+04,
                                 0.86169E+04, 0.91666E+04, 0.97446E+04, 0.10352E+05, 0.10990E+05,
                                 0.11660E+05, 0.12363E+05, 0.13101E+05, 0.13874E+05, 0.14683E+05,
                                 0.15531E+05, 0.16418E+05, 0.17347E+05, 0.18317E+05, 0.19332E+05,
                                 0.20392E+05, 0.21499E+05, 0.22654E+05, 0.23859E+05, 0.25116E+05,
                                 0.26426E+05, 0.27792E+05, 0.29214E+05, 0.30695E+05, 0.32236E+05,
                                 0.33840E+05, 0.35508E+05, 0.37242E+05, 0.39045E+05, 0.40917E+05,
                                 0.42862E+05, 0.44881E+05, 0.46977E+05, 0.49152E+05, 0.51407E+05,
                                 0.53746E+05, 0.56171E+05, 0.58683E+05, 0.61286E+05, 0.63981E+05,
                                 0.66772E+05, 0.69661E+05, 0.72650E+05, 0.75742E+05, 0.78940E+05,
                                 0.82246E+05, 0.85664E+05, 0.89196E+05, 0.92845E+05, 0.96613E+05,
                                 0.10050E+06, 0.10452E+06, 0.10867E+06, 0.11295E+06, 0.11736E+06,
                                 0.12191E+06, 0.12661E+06, 0.13145E+06, 0.13643E+06, 0.14157E+06,
                                 0.14687E+06, 0.15232E+06, 0.15794E+06, 0.16372E+06, 0.16968E+06,
                                 0.17580E+06, 0.18211E+06, 0.18859E+06, 0.19526E+06, 0.20213E+06,
                                 0.20918E+06, 0.21643E+06, 0.22388E+06, 0.23154E+06, 0.23941E+06,
                                 0.24750E+06])


#  --------------- CO2 627: M = 2, I = 4 ---------------------
M = 2
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.66338E+03, 0.93923E+03, 0.12156E+04,
                                 0.14938E+04, 0.17766E+04, 0.20676E+04, 0.23705E+04, 0.26891E+04,
                                 0.30267E+04, 0.33866E+04, 0.37714E+04, 0.41839E+04, 0.46267E+04,
                                 0.51023E+04, 0.56132E+04, 0.61618E+04, 0.67508E+04, 0.73827E+04,
                                 0.80603E+04, 0.87863E+04, 0.95636E+04, 0.10395E+05, 0.11284E+05,
                                 0.12233E+05, 0.13246E+05, 0.14326E+05, 0.15477E+05, 0.16702E+05,
                                 0.18005E+05, 0.19390E+05, 0.20861E+05, 0.22422E+05, 0.24077E+05,
                                 0.25832E+05, 0.27689E+05, 0.29655E+05, 0.31734E+05, 0.33931E+05,
                                 0.36250E+05, 0.38698E+05, 0.41280E+05, 0.44002E+05, 0.46869E+05,
                                 0.49886E+05, 0.53062E+05, 0.56400E+05, 0.59909E+05, 0.63594E+05,
                                 0.67462E+05, 0.71521E+05, 0.75777E+05, 0.80238E+05, 0.84911E+05,
                                 0.89804E+05, 0.94925E+05, 0.10028E+06, 0.10588E+06, 0.11173E+06,
                                 0.11785E+06, 0.12423E+06, 0.13090E+06, 0.13785E+06, 0.14510E+06,
                                 0.15265E+06, 0.16053E+06, 0.16873E+06, 0.17727E+06, 0.18615E+06,
                                 0.19540E+06, 0.20501E+06, 0.21501E+06, 0.22540E+06, 0.23619E+06,
                                 0.24740E+06, 0.25904E+06, 0.27112E+06, 0.28365E+06, 0.29664E+06,
                                 0.31012E+06, 0.32409E+06, 0.33856E+06, 0.35356E+06, 0.36908E+06,
                                 0.38516E+06, 0.40180E+06, 0.41902E+06, 0.43683E+06, 0.45525E+06,
                                 0.47429E+06, 0.49397E+06, 0.51431E+06, 0.53532E+06, 0.55702E+06,
                                 0.57943E+06, 0.60256E+06, 0.62644E+06, 0.65107E+06, 0.67648E+06,
                                 0.70269E+06, 0.72972E+06, 0.75758E+06, 0.78629E+06, 0.81588E+06,
                                 0.84636E+06, 0.87775E+06, 0.91008E+06, 0.94337E+06, 0.97763E+06,
                                 0.10129E+07, 0.10492E+07, 0.10865E+07, 0.11249E+07, 0.11644E+07,
                                 0.12050E+07, 0.12467E+07, 0.12896E+07, 0.13337E+07, 0.13789E+07,
                                 0.14255E+07])


#  --------------- CO2 638: M = 2, I = 5 ---------------------
M = 2
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.22737E+03, 0.32194E+03, 0.41671E+03,
                                 0.51226E+03, 0.60963E+03, 0.71017E+03, 0.81528E+03, 0.92628E+03,
                                 0.10444E+04, 0.11707E+04, 0.13061E+04, 0.14518E+04, 0.16085E+04,
                                 0.17772E+04, 0.19588E+04, 0.21542E+04, 0.23644E+04, 0.25903E+04,
                                 0.28330E+04, 0.30934E+04, 0.33726E+04, 0.36717E+04, 0.39918E+04,
                                 0.43342E+04, 0.47001E+04, 0.50907E+04, 0.55074E+04, 0.59515E+04,
                                 0.64244E+04, 0.69276E+04, 0.74626E+04, 0.80310E+04, 0.86344E+04,
                                 0.92744E+04, 0.99528E+04, 0.10671E+05, 0.11432E+05, 0.12236E+05,
                                 0.13086E+05, 0.13984E+05, 0.14932E+05, 0.15932E+05, 0.16985E+05,
                                 0.18096E+05, 0.19265E+05, 0.20495E+05, 0.21788E+05, 0.23148E+05,
                                 0.24576E+05, 0.26075E+05, 0.27648E+05, 0.29298E+05, 0.31027E+05,
                                 0.32839E+05, 0.34736E+05, 0.36721E+05, 0.38798E+05, 0.40970E+05,
                                 0.43240E+05, 0.45611E+05, 0.48087E+05, 0.50671E+05, 0.53368E+05,
                                 0.56180E+05, 0.59111E+05, 0.62165E+05, 0.65347E+05, 0.68659E+05,
                                 0.72107E+05, 0.75694E+05, 0.79425E+05, 0.83303E+05, 0.87334E+05,
                                 0.91522E+05, 0.95872E+05, 0.10039E+06, 0.10507E+06, 0.10994E+06,
                                 0.11498E+06, 0.12021E+06, 0.12563E+06, 0.13125E+06, 0.13707E+06,
                                 0.14309E+06, 0.14933E+06, 0.15579E+06, 0.16247E+06, 0.16938E+06,
                                 0.17653E+06, 0.18392E+06, 0.19156E+06, 0.19946E+06, 0.20761E+06,
                                 0.21604E+06, 0.22473E+06, 0.23371E+06, 0.24298E+06, 0.25254E+06,
                                 0.26240E+06, 0.27258E+06, 0.28307E+06, 0.29388E+06, 0.30502E+06,
                                 0.31651E+06, 0.32834E+06, 0.34052E+06, 0.35307E+06, 0.36599E+06,
                                 0.37929E+06, 0.39298E+06, 0.40706E+06, 0.42155E+06, 0.43645E+06,
                                 0.45178E+06, 0.46753E+06, 0.48373E+06, 0.50038E+06, 0.51748E+06,
                                 0.53506E+06])


#  --------------- CO2 637: M = 2, I = 6 ---------------------
M = 2
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.13267E+04, 0.18785E+04, 0.24314E+04,
                                 0.29888E+04, 0.35566E+04, 0.41426E+04, 0.47550E+04, 0.54013E+04,
                                 0.60886E+04, 0.68232E+04, 0.76109E+04, 0.84574E+04, 0.93678E+04,
                                 0.10348E+05, 0.11402E+05, 0.12536E+05, 0.13755E+05, 0.15065E+05,
                                 0.16471E+05, 0.17980E+05, 0.19598E+05, 0.21330E+05, 0.23184E+05,
                                 0.25166E+05, 0.27283E+05, 0.29543E+05, 0.31953E+05, 0.34521E+05,
                                 0.37256E+05, 0.40164E+05, 0.43256E+05, 0.46541E+05, 0.50026E+05,
                                 0.53723E+05, 0.57641E+05, 0.61790E+05, 0.66180E+05, 0.70823E+05,
                                 0.75729E+05, 0.80910E+05, 0.86378E+05, 0.92145E+05, 0.98224E+05,
                                 0.10463E+06, 0.11137E+06, 0.11846E+06, 0.12592E+06, 0.13375E+06,
                                 0.14198E+06, 0.15062E+06, 0.15969E+06, 0.16920E+06, 0.17916E+06,
                                 0.18959E+06, 0.20052E+06, 0.21196E+06, 0.22392E+06, 0.23642E+06,
                                 0.24949E+06, 0.26314E+06, 0.27740E+06, 0.29227E+06, 0.30779E+06,
                                 0.32398E+06, 0.34085E+06, 0.35842E+06, 0.37673E+06, 0.39579E+06,
                                 0.41563E+06, 0.43626E+06, 0.45772E+06, 0.48003E+06, 0.50322E+06,
                                 0.52730E+06, 0.55232E+06, 0.57829E+06, 0.60524E+06, 0.63320E+06,
                                 0.66219E+06, 0.69226E+06, 0.72342E+06, 0.75571E+06, 0.78916E+06,
                                 0.82380E+06, 0.85966E+06, 0.89678E+06, 0.93518E+06, 0.97490E+06,
                                 0.10160E+07, 0.10585E+07, 0.11023E+07, 0.11477E+07, 0.11946E+07,
                                 0.12430E+07, 0.12929E+07, 0.13445E+07, 0.13977E+07, 0.14526E+07,
                                 0.15093E+07, 0.15677E+07, 0.16280E+07, 0.16901E+07, 0.17541E+07,
                                 0.18200E+07, 0.18880E+07, 0.19579E+07, 0.20300E+07, 0.21042E+07,
                                 0.21805E+07, 0.22591E+07, 0.23400E+07, 0.24232E+07, 0.25087E+07,
                                 0.25967E+07, 0.26871E+07, 0.27801E+07, 0.28757E+07, 0.29739E+07,
                                 0.30747E+07])


#  --------------- CO2 828: M = 2, I = 7 ---------------------
M = 2
I = 7
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.60334E+02, 0.85430E+02, 0.11058E+03,
                                 0.13590E+03, 0.16167E+03, 0.18821E+03, 0.21588E+03, 0.24502E+03,
                                 0.27595E+03, 0.30896E+03, 0.34431E+03, 0.38225E+03, 0.42301E+03,
                                 0.46684E+03, 0.51397E+03, 0.56464E+03, 0.61907E+03, 0.67753E+03,
                                 0.74027E+03, 0.80753E+03, 0.87961E+03, 0.95676E+03, 0.10393E+04,
                                 0.11275E+04, 0.12217E+04, 0.13222E+04, 0.14293E+04, 0.15434E+04,
                                 0.16648E+04, 0.17940E+04, 0.19312E+04, 0.20769E+04, 0.22315E+04,
                                 0.23954E+04, 0.25691E+04, 0.27529E+04, 0.29474E+04, 0.31530E+04,
                                 0.33702E+04, 0.35995E+04, 0.38414E+04, 0.40965E+04, 0.43654E+04,
                                 0.46484E+04, 0.49464E+04, 0.52598E+04, 0.55892E+04, 0.59353E+04,
                                 0.62988E+04, 0.66803E+04, 0.70804E+04, 0.74998E+04, 0.79394E+04,
                                 0.83998E+04, 0.88817E+04, 0.93859E+04, 0.99132E+04, 0.10464E+05,
                                 0.11040E+05, 0.11642E+05, 0.12270E+05, 0.12925E+05, 0.13609E+05,
                                 0.14321E+05, 0.15064E+05, 0.15838E+05, 0.16643E+05, 0.17482E+05,
                                 0.18355E+05, 0.19263E+05, 0.20207E+05, 0.21188E+05, 0.22208E+05,
                                 0.23267E+05, 0.24366E+05, 0.25508E+05, 0.26692E+05, 0.27921E+05,
                                 0.29195E+05, 0.30516E+05, 0.31886E+05, 0.33304E+05, 0.34773E+05,
                                 0.36294E+05, 0.37869E+05, 0.39499E+05, 0.41185E+05, 0.42929E+05,
                                 0.44732E+05, 0.46596E+05, 0.48522E+05, 0.50513E+05, 0.52569E+05,
                                 0.54692E+05, 0.56884E+05, 0.59146E+05, 0.61481E+05, 0.63890E+05,
                                 0.66375E+05, 0.68937E+05, 0.71578E+05, 0.74301E+05, 0.77107E+05,
                                 0.79998E+05, 0.82976E+05, 0.86043E+05, 0.89201E+05, 0.92452E+05,
                                 0.95799E+05, 0.99242E+05, 0.10278E+06, 0.10643E+06, 0.11018E+06,
                                 0.11403E+06, 0.11799E+06, 0.12206E+06, 0.12625E+06, 0.13055E+06,
                                 0.13497E+06])


#  --------------- CO2 728: M = 2, I = 8 ---------------------
M = 2
I = 8
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.70354E+03, 0.99615E+03, 0.12893E+04,
                                 0.15846E+04, 0.18848E+04, 0.21940E+04, 0.25162E+04, 0.28554E+04,
                                 0.32152E+04, 0.35991E+04, 0.40099E+04, 0.44507E+04, 0.49242E+04,
                                 0.54332E+04, 0.59802E+04, 0.65681E+04, 0.71996E+04, 0.78776E+04,
                                 0.86050E+04, 0.93847E+04, 0.10220E+05, 0.11114E+05, 0.12070E+05,
                                 0.13091E+05, 0.14182E+05, 0.15345E+05, 0.16585E+05, 0.17906E+05,
                                 0.19311E+05, 0.20805E+05, 0.22393E+05, 0.24078E+05, 0.25865E+05,
                                 0.27760E+05, 0.29768E+05, 0.31893E+05, 0.34140E+05, 0.36516E+05,
                                 0.39025E+05, 0.41674E+05, 0.44469E+05, 0.47416E+05, 0.50520E+05,
                                 0.53789E+05, 0.57229E+05, 0.60847E+05, 0.64650E+05, 0.68645E+05,
                                 0.72840E+05, 0.77242E+05, 0.81859E+05, 0.86699E+05, 0.91770E+05,
                                 0.97081E+05, 0.10264E+06, 0.10846E+06, 0.11454E+06, 0.12090E+06,
                                 0.12754E+06, 0.13447E+06, 0.14171E+06, 0.14927E+06, 0.15715E+06,
                                 0.16536E+06, 0.17392E+06, 0.18284E+06, 0.19213E+06, 0.20179E+06,
                                 0.21185E+06, 0.22231E+06, 0.23319E+06, 0.24450E+06, 0.25625E+06,
                                 0.26845E+06, 0.28112E+06, 0.29427E+06, 0.30791E+06, 0.32206E+06,
                                 0.33674E+06, 0.35196E+06, 0.36772E+06, 0.38406E+06, 0.40098E+06,
                                 0.41850E+06, 0.43663E+06, 0.45539E+06, 0.47480E+06, 0.49488E+06,
                                 0.51564E+06, 0.53710E+06, 0.55928E+06, 0.58219E+06, 0.60586E+06,
                                 0.63029E+06, 0.65553E+06, 0.68157E+06, 0.70844E+06, 0.73616E+06,
                                 0.76476E+06, 0.79424E+06, 0.82464E+06, 0.85597E+06, 0.88826E+06,
                                 0.92153E+06, 0.95580E+06, 0.99108E+06, 0.10274E+07, 0.10648E+07,
                                 0.11033E+07, 0.11429E+07, 0.11837E+07, 0.12256E+07, 0.12687E+07,
                                 0.13131E+07, 0.13586E+07, 0.14055E+07, 0.14536E+07, 0.15031E+07,
                                 0.15539E+07])


#  --------------- CO2 727: M = 2, I = 9 ---------------------
M = 2
I = 9
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.20518E+04, 0.29051E+04, 0.37601E+04,
                                 0.46209E+04, 0.54961E+04, 0.63969E+04, 0.73353E+04, 0.83227E+04,
                                 0.93698E+04, 0.10486E+05, 0.11681E+05, 0.12962E+05, 0.14337E+05,
                                 0.15815E+05, 0.17403E+05, 0.19110E+05, 0.20942E+05, 0.22909E+05,
                                 0.25018E+05, 0.27278E+05, 0.29699E+05, 0.32290E+05, 0.35060E+05,
                                 0.38019E+05, 0.41177E+05, 0.44545E+05, 0.48135E+05, 0.51957E+05,
                                 0.56023E+05, 0.60346E+05, 0.64938E+05, 0.69812E+05, 0.74981E+05,
                                 0.80461E+05, 0.86264E+05, 0.92406E+05, 0.98902E+05, 0.10577E+06,
                                 0.11302E+06, 0.12067E+06, 0.12875E+06, 0.13726E+06, 0.14622E+06,
                                 0.15566E+06, 0.16559E+06, 0.17604E+06, 0.18702E+06, 0.19855E+06,
                                 0.21066E+06, 0.22336E+06, 0.23669E+06, 0.25065E+06, 0.26528E+06,
                                 0.28061E+06, 0.29664E+06, 0.31342E+06, 0.33096E+06, 0.34930E+06,
                                 0.36845E+06, 0.38845E+06, 0.40933E+06, 0.43111E+06, 0.45383E+06,
                                 0.47751E+06, 0.50219E+06, 0.52790E+06, 0.55466E+06, 0.58252E+06,
                                 0.61151E+06, 0.64166E+06, 0.67300E+06, 0.70558E+06, 0.73943E+06,
                                 0.77458E+06, 0.81108E+06, 0.84896E+06, 0.88827E+06, 0.92904E+06,
                                 0.97131E+06, 0.10151E+07, 0.10605E+07, 0.11076E+07, 0.11563E+07,
                                 0.12068E+07, 0.12590E+07, 0.13130E+07, 0.13689E+07, 0.14267E+07,
                                 0.14865E+07, 0.15483E+07, 0.16121E+07, 0.16781E+07, 0.17462E+07,
                                 0.18165E+07, 0.18892E+07, 0.19641E+07, 0.20415E+07, 0.21213E+07,
                                 0.22036E+07, 0.22884E+07, 0.23759E+07, 0.24661E+07, 0.25590E+07,
                                 0.26547E+07, 0.27533E+07, 0.28549E+07, 0.29594E+07, 0.30670E+07,
                                 0.31778E+07, 0.32918E+07, 0.34090E+07, 0.35296E+07, 0.36536E+07,
                                 0.37812E+07, 0.39123E+07, 0.40470E+07, 0.41855E+07, 0.43278E+07,
                                 0.44739E+07])


#  --------------- CO2 838: M = 2, I = 10 ---------------------
M = 2
I = 10
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.12066E+03, 0.17085E+03, 0.22116E+03,
                                 0.27190E+03, 0.32364E+03, 0.37711E+03, 0.43305E+03, 0.49219E+03,
                                 0.55516E+03, 0.62256E+03, 0.69492E+03, 0.77276E+03, 0.85657E+03,
                                 0.94685E+03, 0.10441E+04, 0.11488E+04, 0.12614E+04, 0.13826E+04,
                                 0.15127E+04, 0.16525E+04, 0.18024E+04, 0.19630E+04, 0.21351E+04,
                                 0.23191E+04, 0.25158E+04, 0.27260E+04, 0.29502E+04, 0.31892E+04,
                                 0.34438E+04, 0.37148E+04, 0.40031E+04, 0.43094E+04, 0.46346E+04,
                                 0.49797E+04, 0.53455E+04, 0.57331E+04, 0.61434E+04, 0.65775E+04,
                                 0.70364E+04, 0.75212E+04, 0.80330E+04, 0.85730E+04, 0.91424E+04,
                                 0.97423E+04, 0.10374E+05, 0.11039E+05, 0.11738E+05, 0.12474E+05,
                                 0.13246E+05, 0.14057E+05, 0.14908E+05, 0.15801E+05, 0.16737E+05,
                                 0.17717E+05, 0.18744E+05, 0.19819E+05, 0.20944E+05, 0.22120E+05,
                                 0.23349E+05, 0.24634E+05, 0.25975E+05, 0.27376E+05, 0.28837E+05,
                                 0.30361E+05, 0.31950E+05, 0.33605E+05, 0.35330E+05, 0.37126E+05,
                                 0.38996E+05, 0.40942E+05, 0.42965E+05, 0.45069E+05, 0.47256E+05,
                                 0.49528E+05, 0.51888E+05, 0.54338E+05, 0.56882E+05, 0.59521E+05,
                                 0.62259E+05, 0.65097E+05, 0.68040E+05, 0.71090E+05, 0.74249E+05,
                                 0.77522E+05, 0.80910E+05, 0.84417E+05, 0.88046E+05, 0.91801E+05,
                                 0.95684E+05, 0.99699E+05, 0.10385E+06, 0.10814E+06, 0.11257E+06,
                                 0.11715E+06, 0.12187E+06, 0.12675E+06, 0.13179E+06, 0.13699E+06,
                                 0.14235E+06, 0.14788E+06, 0.15358E+06, 0.15946E+06, 0.16552E+06,
                                 0.17176E+06, 0.17819E+06, 0.18482E+06, 0.19164E+06, 0.19867E+06,
                                 0.20590E+06, 0.21335E+06, 0.22101E+06, 0.22889E+06, 0.23699E+06,
                                 0.24533E+06, 0.25390E+06, 0.26271E+06, 0.27177E+06, 0.28108E+06,
                                 0.29064E+06])

#  --------------- CO2 838: M = 2, I = 0 ALIAS-----------------
TIPS_GSI_HASH[(M, 0)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, 0)] = TIPS_ISO_HASH[(M, I)]

#  --------------- CO2 837: M = 2, I = 11 ---------------------
M = 2
I = 11
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.14071E+04, 0.19923E+04, 0.25789E+04,
                                 0.31704E+04, 0.37733E+04, 0.43962E+04, 0.50477E+04, 0.57360E+04,
                                 0.64687E+04, 0.72525E+04, 0.80938E+04, 0.89984E+04, 0.99723E+04,
                                 0.11021E+05, 0.12150E+05, 0.13366E+05, 0.14673E+05, 0.16079E+05,
                                 0.17589E+05, 0.19211E+05, 0.20949E+05, 0.22812E+05, 0.24807E+05,
                                 0.26940E+05, 0.29221E+05, 0.31656E+05, 0.34254E+05, 0.37023E+05,
                                 0.39972E+05, 0.43111E+05, 0.46449E+05, 0.49996E+05, 0.53762E+05,
                                 0.57756E+05, 0.61991E+05, 0.66477E+05, 0.71226E+05, 0.76249E+05,
                                 0.81558E+05, 0.87167E+05, 0.93088E+05, 0.99334E+05, 0.10592E+06,
                                 0.11286E+06, 0.12016E+06, 0.12785E+06, 0.13594E+06, 0.14444E+06,
                                 0.15337E+06, 0.16274E+06, 0.17258E+06, 0.18290E+06, 0.19371E+06,
                                 0.20504E+06, 0.21691E+06, 0.22933E+06, 0.24233E+06, 0.25592E+06,
                                 0.27012E+06, 0.28496E+06, 0.30046E+06, 0.31663E+06, 0.33351E+06,
                                 0.35111E+06, 0.36946E+06, 0.38858E+06, 0.40850E+06, 0.42924E+06,
                                 0.45083E+06, 0.47329E+06, 0.49666E+06, 0.52095E+06, 0.54620E+06,
                                 0.57243E+06, 0.59967E+06, 0.62796E+06, 0.65732E+06, 0.68778E+06,
                                 0.71938E+06, 0.75214E+06, 0.78611E+06, 0.82131E+06, 0.85777E+06,
                                 0.89553E+06, 0.93463E+06, 0.97511E+06, 0.10170E+07, 0.10603E+07,
                                 0.11051E+07, 0.11514E+07, 0.11993E+07, 0.12488E+07, 0.12999E+07,
                                 0.13527E+07, 0.14073E+07, 0.14636E+07, 0.15217E+07, 0.15816E+07,
                                 0.16435E+07, 0.17072E+07, 0.17730E+07, 0.18408E+07, 0.19107E+07,
                                 0.19827E+07, 0.20569E+07, 0.21334E+07, 0.22121E+07, 0.22931E+07,
                                 0.23765E+07, 0.24624E+07, 0.25507E+07, 0.26416E+07, 0.27351E+07,
                                 0.28312E+07, 0.29301E+07, 0.30317E+07, 0.31361E+07, 0.32434E+07,
                                 0.33537E+07])


#  --------------- O3 666: M = 3, I = 1 ---------------------
M = 3
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.30333E+03, 0.51126E+03, 0.75274E+03,
                                 0.10241E+04, 0.13236E+04, 0.16508E+04, 0.20068E+04, 0.23935E+04,
                                 0.28136E+04, 0.32703E+04, 0.37672E+04, 0.43082E+04, 0.48975E+04,
                                 0.55395E+04, 0.62386E+04, 0.69996E+04, 0.78272E+04, 0.87264E+04,
                                 0.97026E+04, 0.10761E+05, 0.11907E+05, 0.13146E+05, 0.14485E+05,
                                 0.15929E+05, 0.17484E+05, 0.19158E+05, 0.20957E+05, 0.22887E+05,
                                 0.24956E+05, 0.27172E+05, 0.29541E+05, 0.32072E+05, 0.34773E+05,
                                 0.37652E+05, 0.40718E+05, 0.43979E+05, 0.47444E+05, 0.51123E+05,
                                 0.55026E+05, 0.59161E+05, 0.63540E+05, 0.68172E+05, 0.73069E+05,
                                 0.78240E+05, 0.83698E+05, 0.89453E+05, 0.95517E+05, 0.10190E+06,
                                 0.10862E+06, 0.11569E+06, 0.12311E+06, 0.13091E+06, 0.13909E+06,
                                 0.14767E+06, 0.15666E+06, 0.16608E+06, 0.17594E+06, 0.18626E+06,
                                 0.19706E+06, 0.20834E+06, 0.22012E+06, 0.23242E+06, 0.24526E+06,
                                 0.25866E+06, 0.27262E+06, 0.28717E+06, 0.30233E+06, 0.31811E+06,
                                 0.33453E+06, 0.35161E+06, 0.36937E+06, 0.38784E+06, 0.40702E+06,
                                 0.42694E+06, 0.44762E+06, 0.46909E+06, 0.49135E+06, 0.51444E+06,
                                 0.53838E+06, 0.56318E+06, 0.58887E+06, 0.61548E+06, 0.64303E+06,
                                 0.67153E+06, 0.70102E+06, 0.73153E+06, 0.76306E+06, 0.79566E+06,
                                 0.82934E+06, 0.86413E+06, 0.90006E+06, 0.93716E+06, 0.97545E+06,
                                 0.10150E+07, 0.10557E+07, 0.10977E+07, 0.11411E+07, 0.11858E+07,
                                 0.12318E+07, 0.12792E+07, 0.13281E+07, 0.13784E+07, 0.14302E+07,
                                 0.14835E+07, 0.15384E+07, 0.15948E+07, 0.16529E+07, 0.17126E+07,
                                 0.17740E+07, 0.18371E+07, 0.19020E+07, 0.19686E+07, 0.20371E+07,
                                 0.21074E+07, 0.21797E+07, 0.22538E+07, 0.23300E+07, 0.24081E+07,
                                 0.24883E+07])


#  --------------- O3 668: M = 3, I = 2 ---------------------
M = 3
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.64763E+03, 0.10916E+04, 0.16073E+04,
                                 0.21870E+04, 0.28271E+04, 0.35272E+04, 0.42900E+04, 0.51197E+04,
                                 0.60225E+04, 0.70057E+04, 0.80771E+04, 0.92455E+04, 0.10520E+05,
                                 0.11911E+05, 0.13427E+05, 0.15079E+05, 0.16878E+05, 0.18834E+05,
                                 0.20960E+05, 0.23267E+05, 0.25767E+05, 0.28472E+05, 0.31397E+05,
                                 0.34553E+05, 0.37957E+05, 0.41620E+05, 0.45559E+05, 0.49790E+05,
                                 0.54327E+05, 0.59187E+05, 0.64387E+05, 0.69944E+05, 0.75877E+05,
                                 0.82203E+05, 0.88943E+05, 0.96114E+05, 0.10374E+06, 0.11184E+06,
                                 0.12043E+06, 0.12954E+06, 0.13918E+06, 0.14939E+06, 0.16018E+06,
                                 0.17159E+06, 0.18362E+06, 0.19632E+06, 0.20970E+06, 0.22380E+06,
                                 0.23863E+06, 0.25423E+06, 0.27063E+06, 0.28786E+06, 0.30594E+06,
                                 0.32490E+06, 0.34478E+06, 0.36561E+06, 0.38743E+06, 0.41026E+06,
                                 0.43413E+06, 0.45909E+06, 0.48517E+06, 0.51241E+06, 0.54084E+06,
                                 0.57049E+06, 0.60141E+06, 0.63365E+06, 0.66722E+06, 0.70219E+06,
                                 0.73858E+06, 0.77644E+06, 0.81581E+06, 0.85674E+06, 0.89927E+06,
                                 0.94345E+06, 0.98932E+06, 0.10369E+07, 0.10863E+07, 0.11375E+07,
                                 0.11906E+07, 0.12457E+07, 0.13027E+07, 0.13618E+07, 0.14229E+07,
                                 0.14862E+07, 0.15517E+07, 0.16194E+07, 0.16894E+07, 0.17618E+07,
                                 0.18366E+07, 0.19139E+07, 0.19937E+07, 0.20761E+07, 0.21612E+07,
                                 0.22490E+07, 0.23395E+07, 0.24330E+07, 0.25293E+07, 0.26286E+07,
                                 0.27309E+07, 0.28363E+07, 0.29449E+07, 0.30568E+07, 0.31720E+07,
                                 0.32905E+07, 0.34125E+07, 0.35381E+07, 0.36672E+07, 0.38000E+07,
                                 0.39366E+07, 0.40770E+07, 0.42213E+07, 0.43696E+07, 0.45220E+07,
                                 0.46785E+07, 0.48392E+07, 0.50043E+07, 0.51737E+07, 0.53476E+07,
                                 0.55261E+07])


#  --------------- O3 686: M = 3, I = 3 ---------------------
M = 3
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.31656E+03, 0.53355E+03, 0.78557E+03,
                                 0.10688E+04, 0.13815E+04, 0.17235E+04, 0.20960E+04, 0.25011E+04,
                                 0.29420E+04, 0.34223E+04, 0.39459E+04, 0.45172E+04, 0.51408E+04,
                                 0.58213E+04, 0.65639E+04, 0.73735E+04, 0.82555E+04, 0.92152E+04,
                                 0.10259E+05, 0.11391E+05, 0.12619E+05, 0.13949E+05, 0.15387E+05,
                                 0.16940E+05, 0.18614E+05, 0.20417E+05, 0.22357E+05, 0.24440E+05,
                                 0.26675E+05, 0.29070E+05, 0.31633E+05, 0.34374E+05, 0.37299E+05,
                                 0.40420E+05, 0.43746E+05, 0.47285E+05, 0.51049E+05, 0.55047E+05,
                                 0.59289E+05, 0.63788E+05, 0.68554E+05, 0.73598E+05, 0.78932E+05,
                                 0.84568E+05, 0.90519E+05, 0.96796E+05, 0.10341E+06, 0.11039E+06,
                                 0.11772E+06, 0.12544E+06, 0.13356E+06, 0.14208E+06, 0.15103E+06,
                                 0.16041E+06, 0.17026E+06, 0.18057E+06, 0.19137E+06, 0.20268E+06,
                                 0.21450E+06, 0.22687E+06, 0.23979E+06, 0.25328E+06, 0.26736E+06,
                                 0.28206E+06, 0.29738E+06, 0.31336E+06, 0.33000E+06, 0.34733E+06,
                                 0.36537E+06, 0.38414E+06, 0.40366E+06, 0.42396E+06, 0.44505E+06,
                                 0.46696E+06, 0.48971E+06, 0.51332E+06, 0.53782E+06, 0.56323E+06,
                                 0.58958E+06, 0.61689E+06, 0.64518E+06, 0.67448E+06, 0.70482E+06,
                                 0.73623E+06, 0.76872E+06, 0.80234E+06, 0.83710E+06, 0.87303E+06,
                                 0.91017E+06, 0.94853E+06, 0.98816E+06, 0.10291E+07, 0.10713E+07,
                                 0.11149E+07, 0.11599E+07, 0.12063E+07, 0.12541E+07, 0.13034E+07,
                                 0.13542E+07, 0.14066E+07, 0.14606E+07, 0.15161E+07, 0.15733E+07,
                                 0.16322E+07, 0.16928E+07, 0.17552E+07, 0.18194E+07, 0.18854E+07,
                                 0.19532E+07, 0.20230E+07, 0.20947E+07, 0.21684E+07, 0.22441E+07,
                                 0.23219E+07, 0.24018E+07, 0.24838E+07, 0.25680E+07, 0.26545E+07,
                                 0.27432E+07])


#  --------------- O3 667: M = 3, I = 4 ---------------------
M = 3
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.37657E+04, 0.63472E+04, 0.93454E+04,
                                 0.12715E+05, 0.16435E+05, 0.20502E+05, 0.24929E+05, 0.29742E+05,
                                 0.34975E+05, 0.40668E+05, 0.46868E+05, 0.53624E+05, 0.60990E+05,
                                 0.69018E+05, 0.77768E+05, 0.87296E+05, 0.97666E+05, 0.10894E+06,
                                 0.12118E+06, 0.13446E+06, 0.14885E+06, 0.16441E+06, 0.18123E+06,
                                 0.19938E+06, 0.21894E+06, 0.23998E+06, 0.26261E+06, 0.28690E+06,
                                 0.31295E+06, 0.34084E+06, 0.37068E+06, 0.40256E+06, 0.43659E+06,
                                 0.47287E+06, 0.51151E+06, 0.55262E+06, 0.59632E+06, 0.64272E+06,
                                 0.69194E+06, 0.74412E+06, 0.79937E+06, 0.85783E+06, 0.91963E+06,
                                 0.98492E+06, 0.10538E+07, 0.11265E+07, 0.12031E+07, 0.12837E+07,
                                 0.13686E+07, 0.14579E+07, 0.15517E+07, 0.16502E+07, 0.17536E+07,
                                 0.18621E+07, 0.19758E+07, 0.20949E+07, 0.22196E+07, 0.23501E+07,
                                 0.24866E+07, 0.26292E+07, 0.27783E+07, 0.29339E+07, 0.30963E+07,
                                 0.32658E+07, 0.34425E+07, 0.36266E+07, 0.38184E+07, 0.40181E+07,
                                 0.42260E+07, 0.44422E+07, 0.46671E+07, 0.49008E+07, 0.51437E+07,
                                 0.53959E+07, 0.56578E+07, 0.59296E+07, 0.62116E+07, 0.65040E+07,
                                 0.68071E+07, 0.71213E+07, 0.74468E+07, 0.77838E+07, 0.81328E+07,
                                 0.84939E+07, 0.88676E+07, 0.92541E+07, 0.96536E+07, 0.10067E+08,
                                 0.10493E+08, 0.10934E+08, 0.11390E+08, 0.11860E+08, 0.12345E+08,
                                 0.12846E+08, 0.13363E+08, 0.13895E+08, 0.14445E+08, 0.15011E+08,
                                 0.15595E+08, 0.16196E+08, 0.16815E+08, 0.17453E+08, 0.18110E+08,
                                 0.18786E+08, 0.19482E+08, 0.20198E+08, 0.20934E+08, 0.21691E+08,
                                 0.22470E+08, 0.23270E+08, 0.24093E+08, 0.24939E+08, 0.25807E+08,
                                 0.26699E+08, 0.27616E+08, 0.28556E+08, 0.29522E+08, 0.30514E+08,
                                 0.31531E+08])


#  --------------- O3 676: M = 3, I = 5 ---------------------
M = 3
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.18608E+04, 0.31363E+04, 0.46177E+04,
                                 0.62826E+04, 0.81202E+04, 0.10129E+05, 0.12316E+05, 0.14693E+05,
                                 0.17277E+05, 0.20089E+05, 0.23153E+05, 0.26492E+05, 0.30133E+05,
                                 0.34103E+05, 0.38430E+05, 0.43145E+05, 0.48277E+05, 0.53858E+05,
                                 0.59920E+05, 0.66497E+05, 0.73624E+05, 0.81336E+05, 0.89671E+05,
                                 0.98668E+05, 0.10836E+06, 0.11880E+06, 0.13002E+06, 0.14207E+06,
                                 0.15500E+06, 0.16884E+06, 0.18365E+06, 0.19947E+06, 0.21636E+06,
                                 0.23438E+06, 0.25356E+06, 0.27398E+06, 0.29568E+06, 0.31873E+06,
                                 0.34318E+06, 0.36911E+06, 0.39656E+06, 0.42561E+06, 0.45632E+06,
                                 0.48877E+06, 0.52302E+06, 0.55914E+06, 0.59722E+06, 0.63732E+06,
                                 0.67952E+06, 0.72390E+06, 0.77055E+06, 0.81954E+06, 0.87097E+06,
                                 0.92491E+06, 0.98146E+06, 0.10407E+07, 0.11027E+07, 0.11677E+07,
                                 0.12356E+07, 0.13066E+07, 0.13807E+07, 0.14582E+07, 0.15390E+07,
                                 0.16233E+07, 0.17113E+07, 0.18029E+07, 0.18984E+07, 0.19978E+07,
                                 0.21012E+07, 0.22089E+07, 0.23208E+07, 0.24372E+07, 0.25581E+07,
                                 0.26837E+07, 0.28141E+07, 0.29494E+07, 0.30898E+07, 0.32354E+07,
                                 0.33864E+07, 0.35428E+07, 0.37049E+07, 0.38728E+07, 0.40466E+07,
                                 0.42264E+07, 0.44125E+07, 0.46050E+07, 0.48040E+07, 0.50098E+07,
                                 0.52224E+07, 0.54420E+07, 0.56689E+07, 0.59031E+07, 0.61449E+07,
                                 0.63943E+07, 0.66517E+07, 0.69172E+07, 0.71909E+07, 0.74731E+07,
                                 0.77639E+07, 0.80635E+07, 0.83721E+07, 0.86900E+07, 0.90172E+07,
                                 0.93541E+07, 0.97008E+07, 0.10058E+08, 0.10424E+08, 0.10802E+08,
                                 0.11190E+08, 0.11589E+08, 0.11999E+08, 0.12420E+08, 0.12853E+08,
                                 0.13298E+08, 0.13755E+08, 0.14223E+08, 0.14705E+08, 0.15199E+08,
                                 0.15706E+08])


#  --------------- O3 886: M = 3, I = 6 ---------------------
M = 3
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.67639E+03, 0.11401E+04, 0.16787E+04,
                                 0.22843E+04, 0.29532E+04, 0.36856E+04, 0.44842E+04, 0.53545E+04,
                                 0.63030E+04, 0.73381E+04, 0.84686E+04, 0.97040E+04, 0.11054E+05,
                                 0.12530E+05, 0.14143E+05, 0.15903E+05, 0.17823E+05, 0.19915E+05,
                                 0.22190E+05, 0.24663E+05, 0.27346E+05, 0.30254E+05, 0.33400E+05,
                                 0.36800E+05, 0.40469E+05, 0.44423E+05, 0.48678E+05, 0.53251E+05,
                                 0.58160E+05, 0.63423E+05, 0.69058E+05, 0.75085E+05, 0.81524E+05,
                                 0.88395E+05, 0.95719E+05, 0.10352E+06, 0.11181E+06, 0.12063E+06,
                                 0.12999E+06, 0.13991E+06, 0.15043E+06, 0.16157E+06, 0.17335E+06,
                                 0.18580E+06, 0.19895E+06, 0.21283E+06, 0.22746E+06, 0.24288E+06,
                                 0.25911E+06, 0.27619E+06, 0.29415E+06, 0.31301E+06, 0.33283E+06,
                                 0.35362E+06, 0.37542E+06, 0.39827E+06, 0.42221E+06, 0.44726E+06,
                                 0.47348E+06, 0.50089E+06, 0.52954E+06, 0.55947E+06, 0.59072E+06,
                                 0.62332E+06, 0.65733E+06, 0.69279E+06, 0.72973E+06, 0.76821E+06,
                                 0.80827E+06, 0.84996E+06, 0.89332E+06, 0.93840E+06, 0.98526E+06,
                                 0.10339E+07, 0.10845E+07, 0.11370E+07, 0.11914E+07, 0.12479E+07,
                                 0.13065E+07, 0.13672E+07, 0.14302E+07, 0.14953E+07, 0.15628E+07,
                                 0.16327E+07, 0.17050E+07, 0.17798E+07, 0.18571E+07, 0.19371E+07,
                                 0.20197E+07, 0.21051E+07, 0.21933E+07, 0.22844E+07, 0.23785E+07,
                                 0.24755E+07, 0.25757E+07, 0.26790E+07, 0.27855E+07, 0.28954E+07,
                                 0.30086E+07, 0.31253E+07, 0.32455E+07, 0.33693E+07, 0.34967E+07,
                                 0.36280E+07, 0.37631E+07, 0.39021E+07, 0.40451E+07, 0.41922E+07,
                                 0.43435E+07, 0.44990E+07, 0.46589E+07, 0.48232E+07, 0.49920E+07,
                                 0.51654E+07, 0.53436E+07, 0.55265E+07, 0.57143E+07, 0.59071E+07,
                                 0.61050E+07])


#  --------------- O3 868: M = 3, I = 7 ---------------------
M = 3
I = 7
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.34615E+03, 0.58348E+03, 0.85915E+03,
                                 0.11692E+04, 0.15117E+04, 0.18868E+04, 0.22960E+04, 0.27419E+04,
                                 0.32278E+04, 0.37579E+04, 0.43366E+04, 0.49686E+04, 0.56591E+04,
                                 0.64134E+04, 0.72369E+04, 0.81354E+04, 0.91148E+04, 0.10181E+05,
                                 0.11341E+05, 0.12600E+05, 0.13966E+05, 0.15446E+05, 0.17046E+05,
                                 0.18775E+05, 0.20640E+05, 0.22649E+05, 0.24810E+05, 0.27132E+05,
                                 0.29624E+05, 0.32295E+05, 0.35154E+05, 0.38211E+05, 0.41475E+05,
                                 0.44958E+05, 0.48670E+05, 0.52621E+05, 0.56823E+05, 0.61288E+05,
                                 0.66026E+05, 0.71052E+05, 0.76376E+05, 0.82011E+05, 0.87972E+05,
                                 0.94271E+05, 0.10092E+06, 0.10794E+06, 0.11534E+06, 0.12313E+06,
                                 0.13134E+06, 0.13997E+06, 0.14905E+06, 0.15858E+06, 0.16859E+06,
                                 0.17909E+06, 0.19010E+06, 0.20164E+06, 0.21373E+06, 0.22638E+06,
                                 0.23962E+06, 0.25346E+06, 0.26792E+06, 0.28302E+06, 0.29879E+06,
                                 0.31524E+06, 0.33240E+06, 0.35029E+06, 0.36892E+06, 0.38833E+06,
                                 0.40853E+06, 0.42956E+06, 0.45142E+06, 0.47416E+06, 0.49778E+06,
                                 0.52233E+06, 0.54781E+06, 0.57427E+06, 0.60172E+06, 0.63019E+06,
                                 0.65971E+06, 0.69031E+06, 0.72201E+06, 0.75485E+06, 0.78886E+06,
                                 0.82405E+06, 0.86048E+06, 0.89815E+06, 0.93711E+06, 0.97739E+06,
                                 0.10190E+07, 0.10620E+07, 0.11065E+07, 0.11523E+07, 0.11997E+07,
                                 0.12485E+07, 0.12990E+07, 0.13510E+07, 0.14046E+07, 0.14599E+07,
                                 0.15169E+07, 0.15756E+07, 0.16361E+07, 0.16984E+07, 0.17626E+07,
                                 0.18287E+07, 0.18966E+07, 0.19666E+07, 0.20386E+07, 0.21126E+07,
                                 0.21887E+07, 0.22669E+07, 0.23474E+07, 0.24300E+07, 0.25150E+07,
                                 0.26022E+07, 0.26919E+07, 0.27839E+07, 0.28784E+07, 0.29753E+07,
                                 0.30749E+07])


#  --------------- O3 678: M = 3, I = 8 ---------------------
M = 3
I = 8
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.39745E+04, 0.66993E+04, 0.98642E+04,
                                 0.13422E+05, 0.17352E+05, 0.21652E+05, 0.26339E+05, 0.31442E+05,
                                 0.37000E+05, 0.43058E+05, 0.49669E+05, 0.56885E+05, 0.64766E+05,
                                 0.73372E+05, 0.82765E+05, 0.93011E+05, 0.10418E+06, 0.11633E+06,
                                 0.12955E+06, 0.14390E+06, 0.15946E+06, 0.17632E+06, 0.19455E+06,
                                 0.21424E+06, 0.23547E+06, 0.25835E+06, 0.28296E+06, 0.30939E+06,
                                 0.33776E+06, 0.36816E+06, 0.40070E+06, 0.43549E+06, 0.47264E+06,
                                 0.51228E+06, 0.55451E+06, 0.59947E+06, 0.64728E+06, 0.69807E+06,
                                 0.75198E+06, 0.80915E+06, 0.86971E+06, 0.93381E+06, 0.10016E+07,
                                 0.10733E+07, 0.11489E+07, 0.12287E+07, 0.13128E+07, 0.14015E+07,
                                 0.14948E+07, 0.15930E+07, 0.16961E+07, 0.18045E+07, 0.19183E+07,
                                 0.20378E+07, 0.21629E+07, 0.22942E+07, 0.24316E+07, 0.25754E+07,
                                 0.27258E+07, 0.28831E+07, 0.30475E+07, 0.32192E+07, 0.33984E+07,
                                 0.35855E+07, 0.37805E+07, 0.39838E+07, 0.41956E+07, 0.44162E+07,
                                 0.46458E+07, 0.48847E+07, 0.51332E+07, 0.53916E+07, 0.56601E+07,
                                 0.59390E+07, 0.62286E+07, 0.65292E+07, 0.68412E+07, 0.71647E+07,
                                 0.75002E+07, 0.78479E+07, 0.82081E+07, 0.85813E+07, 0.89676E+07,
                                 0.93676E+07, 0.97814E+07, 0.10209E+08, 0.10652E+08, 0.11110E+08,
                                 0.11583E+08, 0.12071E+08, 0.12576E+08, 0.13097E+08, 0.13635E+08,
                                 0.14190E+08, 0.14763E+08, 0.15354E+08, 0.15963E+08, 0.16592E+08,
                                 0.17239E+08, 0.17906E+08, 0.18593E+08, 0.19301E+08, 0.20030E+08,
                                 0.20780E+08, 0.21553E+08, 0.22347E+08, 0.23165E+08, 0.24006E+08,
                                 0.24870E+08, 0.25759E+08, 0.26673E+08, 0.27612E+08, 0.28577E+08,
                                 0.29568E+08, 0.30585E+08, 0.31631E+08, 0.32704E+08, 0.33805E+08,
                                 0.34936E+08])


#  --------------- O3 768: M = 3, I = 9 ---------------------
M = 3
I = 9
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.40228E+04, 0.67808E+04, 0.99842E+04,
                                 0.13586E+05, 0.17564E+05, 0.21919E+05, 0.26665E+05, 0.31833E+05,
                                 0.37461E+05, 0.43596E+05, 0.50286E+05, 0.57589E+05, 0.65562E+05,
                                 0.74264E+05, 0.83761E+05, 0.94115E+05, 0.10540E+06, 0.11767E+06,
                                 0.13102E+06, 0.14550E+06, 0.16121E+06, 0.17822E+06, 0.19661E+06,
                                 0.21646E+06, 0.23788E+06, 0.26094E+06, 0.28574E+06, 0.31239E+06,
                                 0.34097E+06, 0.37160E+06, 0.40437E+06, 0.43941E+06, 0.47683E+06,
                                 0.51673E+06, 0.55925E+06, 0.60451E+06, 0.65262E+06, 0.70374E+06,
                                 0.75799E+06, 0.81550E+06, 0.87643E+06, 0.94092E+06, 0.10091E+07,
                                 0.10812E+07, 0.11572E+07, 0.12375E+07, 0.13221E+07, 0.14112E+07,
                                 0.15050E+07, 0.16037E+07, 0.17074E+07, 0.18164E+07, 0.19307E+07,
                                 0.20507E+07, 0.21765E+07, 0.23084E+07, 0.24464E+07, 0.25909E+07,
                                 0.27421E+07, 0.29001E+07, 0.30652E+07, 0.32377E+07, 0.34177E+07,
                                 0.36055E+07, 0.38014E+07, 0.40055E+07, 0.42182E+07, 0.44397E+07,
                                 0.46703E+07, 0.49102E+07, 0.51597E+07, 0.54191E+07, 0.56886E+07,
                                 0.59686E+07, 0.62593E+07, 0.65611E+07, 0.68742E+07, 0.71989E+07,
                                 0.75356E+07, 0.78846E+07, 0.82461E+07, 0.86206E+07, 0.90083E+07,
                                 0.94097E+07, 0.98249E+07, 0.10254E+08, 0.10699E+08, 0.11158E+08,
                                 0.11632E+08, 0.12123E+08, 0.12629E+08, 0.13152E+08, 0.13691E+08,
                                 0.14248E+08, 0.14823E+08, 0.15416E+08, 0.16027E+08, 0.16657E+08,
                                 0.17307E+08, 0.17976E+08, 0.18665E+08, 0.19375E+08, 0.20106E+08,
                                 0.20858E+08, 0.21633E+08, 0.22430E+08, 0.23250E+08, 0.24093E+08,
                                 0.24960E+08, 0.25851E+08, 0.26767E+08, 0.27709E+08, 0.28676E+08,
                                 0.29670E+08, 0.30691E+08, 0.31739E+08, 0.32815E+08, 0.33919E+08,
                                 0.35053E+08])


#  --------------- O3 786: M = 3, I = 10 ---------------------
M = 3
I = 10
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.39315E+04, 0.66267E+04, 0.97569E+04,
                                 0.13276E+05, 0.17162E+05, 0.21414E+05, 0.26048E+05, 0.31094E+05,
                                 0.36590E+05, 0.42581E+05, 0.49120E+05, 0.56260E+05, 0.64061E+05,
                                 0.72580E+05, 0.81882E+05, 0.92031E+05, 0.10309E+06, 0.11514E+06,
                                 0.12824E+06, 0.14247E+06, 0.15791E+06, 0.17463E+06, 0.19272E+06,
                                 0.21226E+06, 0.23333E+06, 0.25604E+06, 0.28047E+06, 0.30673E+06,
                                 0.33490E+06, 0.36510E+06, 0.39743E+06, 0.43200E+06, 0.46892E+06,
                                 0.50831E+06, 0.55029E+06, 0.59498E+06, 0.64251E+06, 0.69301E+06,
                                 0.74662E+06, 0.80347E+06, 0.86370E+06, 0.92747E+06, 0.99491E+06,
                                 0.10662E+07, 0.11414E+07, 0.12208E+07, 0.13046E+07, 0.13928E+07,
                                 0.14856E+07, 0.15833E+07, 0.16860E+07, 0.17939E+07, 0.19072E+07,
                                 0.20261E+07, 0.21508E+07, 0.22814E+07, 0.24182E+07, 0.25614E+07,
                                 0.27112E+07, 0.28679E+07, 0.30316E+07, 0.32026E+07, 0.33811E+07,
                                 0.35674E+07, 0.37617E+07, 0.39642E+07, 0.41752E+07, 0.43950E+07,
                                 0.46237E+07, 0.48618E+07, 0.51094E+07, 0.53668E+07, 0.56343E+07,
                                 0.59123E+07, 0.62009E+07, 0.65005E+07, 0.68113E+07, 0.71338E+07,
                                 0.74681E+07, 0.78147E+07, 0.81737E+07, 0.85457E+07, 0.89308E+07,
                                 0.93295E+07, 0.97420E+07, 0.10169E+08, 0.10610E+08, 0.11066E+08,
                                 0.11538E+08, 0.12025E+08, 0.12528E+08, 0.13048E+08, 0.13584E+08,
                                 0.14138E+08, 0.14709E+08, 0.15298E+08, 0.15906E+08, 0.16532E+08,
                                 0.17178E+08, 0.17843E+08, 0.18528E+08, 0.19234E+08, 0.19961E+08,
                                 0.20710E+08, 0.21480E+08, 0.22272E+08, 0.23088E+08, 0.23926E+08,
                                 0.24789E+08, 0.25675E+08, 0.26587E+08, 0.27523E+08, 0.28485E+08,
                                 0.29474E+08, 0.30489E+08, 0.31532E+08, 0.32603E+08, 0.33701E+08,
                                 0.34829E+08])


#  --------------- O3 776: M = 3, I = 11 ---------------------
M = 3
I = 11
TIPS_GSI_HASH[(M, I)] = __FloatType__(36.)
TIPS_ISO_HASH[(M, I)] = float32([0.23106E+05, 0.38945E+05, 0.57342E+05,
                                 0.78021E+05, 0.10085E+06, 0.12582E+06, 0.15302E+06, 0.18262E+06,
                                 0.21482E+06, 0.24989E+06, 0.28812E+06, 0.32983E+06, 0.37535E+06,
                                 0.42501E+06, 0.47919E+06, 0.53825E+06, 0.60258E+06, 0.67256E+06,
                                 0.74862E+06, 0.83118E+06, 0.92069E+06, 0.10176E+07, 0.11223E+07,
                                 0.12354E+07, 0.13574E+07, 0.14887E+07, 0.16299E+07, 0.17816E+07,
                                 0.19443E+07, 0.21187E+07, 0.23052E+07, 0.25047E+07, 0.27176E+07,
                                 0.29447E+07, 0.31866E+07, 0.34441E+07, 0.37179E+07, 0.40087E+07,
                                 0.43173E+07, 0.46444E+07, 0.49910E+07, 0.53578E+07, 0.57456E+07,
                                 0.61554E+07, 0.65880E+07, 0.70444E+07, 0.75255E+07, 0.80322E+07,
                                 0.85656E+07, 0.91266E+07, 0.97163E+07, 0.10336E+08, 0.10986E+08,
                                 0.11668E+08, 0.12383E+08, 0.13133E+08, 0.13918E+08, 0.14739E+08,
                                 0.15598E+08, 0.16496E+08, 0.17435E+08, 0.18415E+08, 0.19438E+08,
                                 0.20505E+08, 0.21619E+08, 0.22779E+08, 0.23987E+08, 0.25246E+08,
                                 0.26556E+08, 0.27920E+08, 0.29337E+08, 0.30811E+08, 0.32343E+08,
                                 0.33934E+08, 0.35585E+08, 0.37300E+08, 0.39079E+08, 0.40924E+08,
                                 0.42837E+08, 0.44819E+08, 0.46873E+08, 0.49001E+08, 0.51203E+08,
                                 0.53483E+08, 0.55842E+08, 0.58282E+08, 0.60805E+08, 0.63414E+08,
                                 0.66109E+08, 0.68894E+08, 0.71770E+08, 0.74740E+08, 0.77806E+08,
                                 0.80970E+08, 0.84234E+08, 0.87600E+08, 0.91072E+08, 0.94651E+08,
                                 0.98339E+08, 0.10214E+09, 0.10605E+09, 0.11009E+09, 0.11424E+09,
                                 0.11851E+09, 0.12291E+09, 0.12744E+09, 0.13209E+09, 0.13688E+09,
                                 0.14180E+09, 0.14687E+09, 0.15207E+09, 0.15742E+09, 0.16291E+09,
                                 0.16855E+09, 0.17435E+09, 0.18030E+09, 0.18641E+09, 0.19268E+09,
                                 0.19912E+09])


#  --------------- O3 767: M = 3, I = 12 ---------------------
M = 3
I = 12
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.11692E+05, 0.19707E+05, 0.29017E+05,
                                 0.39482E+05, 0.51038E+05, 0.63680E+05, 0.77450E+05, 0.92432E+05,
                                 0.10873E+06, 0.12649E+06, 0.14584E+06, 0.16694E+06, 0.18996E+06,
                                 0.21507E+06, 0.24245E+06, 0.27229E+06, 0.30478E+06, 0.34013E+06,
                                 0.37853E+06, 0.42020E+06, 0.46536E+06, 0.51424E+06, 0.56708E+06,
                                 0.62411E+06, 0.68559E+06, 0.75178E+06, 0.82296E+06, 0.89939E+06,
                                 0.98137E+06, 0.10692E+07, 0.11631E+07, 0.12636E+07, 0.13708E+07,
                                 0.14851E+07, 0.16069E+07, 0.17365E+07, 0.18742E+07, 0.20206E+07,
                                 0.21758E+07, 0.23404E+07, 0.25148E+07, 0.26992E+07, 0.28943E+07,
                                 0.31004E+07, 0.33179E+07, 0.35474E+07, 0.37892E+07, 0.40440E+07,
                                 0.43121E+07, 0.45940E+07, 0.48904E+07, 0.52017E+07, 0.55285E+07,
                                 0.58713E+07, 0.62306E+07, 0.66071E+07, 0.70014E+07, 0.74140E+07,
                                 0.78456E+07, 0.82967E+07, 0.87681E+07, 0.92604E+07, 0.97742E+07,
                                 0.10310E+08, 0.10869E+08, 0.11452E+08, 0.12059E+08, 0.12691E+08,
                                 0.13348E+08, 0.14033E+08, 0.14745E+08, 0.15484E+08, 0.16253E+08,
                                 0.17052E+08, 0.17881E+08, 0.18741E+08, 0.19634E+08, 0.20560E+08,
                                 0.21520E+08, 0.22515E+08, 0.23546E+08, 0.24613E+08, 0.25718E+08,
                                 0.26862E+08, 0.28046E+08, 0.29270E+08, 0.30536E+08, 0.31845E+08,
                                 0.33197E+08, 0.34594E+08, 0.36037E+08, 0.37527E+08, 0.39065E+08,
                                 0.40652E+08, 0.42289E+08, 0.43977E+08, 0.45719E+08, 0.47514E+08,
                                 0.49363E+08, 0.51270E+08, 0.53233E+08, 0.55255E+08, 0.57337E+08,
                                 0.59480E+08, 0.61686E+08, 0.63956E+08, 0.66290E+08, 0.68691E+08,
                                 0.71160E+08, 0.73699E+08, 0.76307E+08, 0.78988E+08, 0.81743E+08,
                                 0.84572E+08, 0.87478E+08, 0.90462E+08, 0.93525E+08, 0.96669E+08,
                                 0.99896E+08])


#  --------------- O3 888: M = 3, I = 13 ---------------------
M = 3
I = 13
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.36175E+03, 0.60978E+03, 0.89790E+03,
                                 0.12219E+04, 0.15802E+04, 0.19728E+04, 0.24016E+04, 0.28696E+04,
                                 0.33807E+04, 0.39394E+04, 0.45506E+04, 0.52196E+04, 0.59521E+04,
                                 0.67538E+04, 0.76308E+04, 0.85894E+04, 0.96361E+04, 0.10777E+05,
                                 0.12021E+05, 0.13373E+05, 0.14841E+05, 0.16434E+05, 0.18158E+05,
                                 0.20023E+05, 0.22037E+05, 0.24208E+05, 0.26547E+05, 0.29061E+05,
                                 0.31762E+05, 0.34659E+05, 0.37762E+05, 0.41083E+05, 0.44632E+05,
                                 0.48421E+05, 0.52462E+05, 0.56766E+05, 0.61346E+05, 0.66215E+05,
                                 0.71386E+05, 0.76873E+05, 0.82688E+05, 0.88848E+05, 0.95365E+05,
                                 0.10226E+06, 0.10954E+06, 0.11722E+06, 0.12532E+06, 0.13387E+06,
                                 0.14286E+06, 0.15233E+06, 0.16229E+06, 0.17275E+06, 0.18374E+06,
                                 0.19528E+06, 0.20737E+06, 0.22006E+06, 0.23335E+06, 0.24726E+06,
                                 0.26182E+06, 0.27705E+06, 0.29297E+06, 0.30960E+06, 0.32696E+06,
                                 0.34509E+06, 0.36399E+06, 0.38371E+06, 0.40425E+06, 0.42566E+06,
                                 0.44794E+06, 0.47114E+06, 0.49527E+06, 0.52036E+06, 0.54644E+06,
                                 0.57354E+06, 0.60169E+06, 0.63091E+06, 0.66124E+06, 0.69270E+06,
                                 0.72533E+06, 0.75916E+06, 0.79421E+06, 0.83053E+06, 0.86814E+06,
                                 0.90708E+06, 0.94737E+06, 0.98907E+06, 0.10322E+07, 0.10768E+07,
                                 0.11229E+07, 0.11705E+07, 0.12197E+07, 0.12705E+07, 0.13230E+07,
                                 0.13771E+07, 0.14330E+07, 0.14906E+07, 0.15501E+07, 0.16114E+07,
                                 0.16745E+07, 0.17397E+07, 0.18067E+07, 0.18759E+07, 0.19470E+07,
                                 0.20203E+07, 0.20957E+07, 0.21733E+07, 0.22532E+07, 0.23353E+07,
                                 0.24198E+07, 0.25067E+07, 0.25960E+07, 0.26878E+07, 0.27821E+07,
                                 0.28790E+07, 0.29785E+07, 0.30807E+07, 0.31857E+07, 0.32934E+07,
                                 0.34040E+07])


#  --------------- O3 887: M = 3, I = 14 ---------------------
M = 3
I = 14
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.42000E+04, 0.70796E+04, 0.10424E+05,
                                 0.14186E+05, 0.18342E+05, 0.22896E+05, 0.27866E+05, 0.33285E+05,
                                 0.39199E+05, 0.45659E+05, 0.52720E+05, 0.60444E+05, 0.68895E+05,
                                 0.78139E+05, 0.88246E+05, 0.99288E+05, 0.11134E+06, 0.12447E+06,
                                 0.13877E+06, 0.15431E+06, 0.17119E+06, 0.18949E+06, 0.20930E+06,
                                 0.23071E+06, 0.25383E+06, 0.27875E+06, 0.30558E+06, 0.33442E+06,
                                 0.36539E+06, 0.39861E+06, 0.43418E+06, 0.47224E+06, 0.51291E+06,
                                 0.55632E+06, 0.60260E+06, 0.65189E+06, 0.70434E+06, 0.76008E+06,
                                 0.81927E+06, 0.88206E+06, 0.94862E+06, 0.10191E+07, 0.10937E+07,
                                 0.11725E+07, 0.12558E+07, 0.13436E+07, 0.14363E+07, 0.15340E+07,
                                 0.16368E+07, 0.17450E+07, 0.18588E+07, 0.19784E+07, 0.21040E+07,
                                 0.22358E+07, 0.23741E+07, 0.25190E+07, 0.26708E+07, 0.28297E+07,
                                 0.29961E+07, 0.31700E+07, 0.33518E+07, 0.35417E+07, 0.37400E+07,
                                 0.39469E+07, 0.41628E+07, 0.43878E+07, 0.46224E+07, 0.48667E+07,
                                 0.51210E+07, 0.53858E+07, 0.56611E+07, 0.59475E+07, 0.62451E+07,
                                 0.65544E+07, 0.68755E+07, 0.72089E+07, 0.75550E+07, 0.79139E+07,
                                 0.82861E+07, 0.86720E+07, 0.90719E+07, 0.94861E+07, 0.99151E+07,
                                 0.10359E+08, 0.10819E+08, 0.11294E+08, 0.11786E+08, 0.12294E+08,
                                 0.12820E+08, 0.13363E+08, 0.13924E+08, 0.14503E+08, 0.15101E+08,
                                 0.15719E+08, 0.16356E+08, 0.17013E+08, 0.17690E+08, 0.18389E+08,
                                 0.19109E+08, 0.19851E+08, 0.20616E+08, 0.21404E+08, 0.22215E+08,
                                 0.23050E+08, 0.23910E+08, 0.24794E+08, 0.25704E+08, 0.26640E+08,
                                 0.27603E+08, 0.28593E+08, 0.29610E+08, 0.30656E+08, 0.31731E+08,
                                 0.32835E+08, 0.33969E+08, 0.35133E+08, 0.36329E+08, 0.37556E+08,
                                 0.38816E+08])


#  --------------- O3 878: M = 3, I = 15 ---------------------
M = 3
I = 15
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.21250E+04, 0.35820E+04, 0.52744E+04,
                                 0.71778E+04, 0.92814E+04, 0.11586E+05, 0.14102E+05, 0.16845E+05,
                                 0.19839E+05, 0.23108E+05, 0.26680E+05, 0.30588E+05, 0.34861E+05,
                                 0.39534E+05, 0.44642E+05, 0.50219E+05, 0.56305E+05, 0.62937E+05,
                                 0.70155E+05, 0.78001E+05, 0.86516E+05, 0.95747E+05, 0.10574E+06,
                                 0.11653E+06, 0.12819E+06, 0.14075E+06, 0.15427E+06, 0.16881E+06,
                                 0.18441E+06, 0.20114E+06, 0.21906E+06, 0.23823E+06, 0.25871E+06,
                                 0.28056E+06, 0.30386E+06, 0.32867E+06, 0.35507E+06, 0.38312E+06,
                                 0.41291E+06, 0.44450E+06, 0.47799E+06, 0.51344E+06, 0.55095E+06,
                                 0.59060E+06, 0.63248E+06, 0.67667E+06, 0.72327E+06, 0.77238E+06,
                                 0.82409E+06, 0.87850E+06, 0.93571E+06, 0.99583E+06, 0.10590E+07,
                                 0.11252E+07, 0.11947E+07, 0.12675E+07, 0.13438E+07, 0.14237E+07,
                                 0.15072E+07, 0.15946E+07, 0.16859E+07, 0.17814E+07, 0.18810E+07,
                                 0.19849E+07, 0.20934E+07, 0.22064E+07, 0.23242E+07, 0.24469E+07,
                                 0.25747E+07, 0.27076E+07, 0.28459E+07, 0.29897E+07, 0.31391E+07,
                                 0.32944E+07, 0.34557E+07, 0.36231E+07, 0.37968E+07, 0.39770E+07,
                                 0.41639E+07, 0.43576E+07, 0.45583E+07, 0.47663E+07, 0.49816E+07,
                                 0.52045E+07, 0.54352E+07, 0.56739E+07, 0.59207E+07, 0.61759E+07,
                                 0.64396E+07, 0.67121E+07, 0.69936E+07, 0.72844E+07, 0.75845E+07,
                                 0.78943E+07, 0.82139E+07, 0.85436E+07, 0.88837E+07, 0.92342E+07,
                                 0.95956E+07, 0.99680E+07, 0.10352E+08, 0.10747E+08, 0.11154E+08,
                                 0.11573E+08, 0.12004E+08, 0.12448E+08, 0.12904E+08, 0.13374E+08,
                                 0.13857E+08, 0.14353E+08, 0.14864E+08, 0.15388E+08, 0.15927E+08,
                                 0.16481E+08, 0.17050E+08, 0.17634E+08, 0.18234E+08, 0.18849E+08,
                                 0.19481E+08])


#  --------------- O3 778: M = 3, I = 16 ---------------------
M = 3
I = 16
TIPS_GSI_HASH[(M, I)] = __FloatType__(36.)
TIPS_ISO_HASH[(M, I)] = float32([0.24692E+05, 0.41621E+05, 0.61284E+05,
                                 0.83394E+05, 0.10782E+06, 0.13457E+06, 0.16375E+06, 0.19554E+06,
                                 0.23020E+06, 0.26801E+06, 0.30930E+06, 0.35443E+06, 0.40375E+06,
                                 0.45763E+06, 0.51650E+06, 0.58075E+06, 0.65080E+06, 0.72711E+06,
                                 0.81012E+06, 0.90030E+06, 0.99815E+06, 0.11042E+07, 0.12189E+07,
                                 0.13428E+07, 0.14765E+07, 0.16206E+07, 0.17757E+07, 0.19423E+07,
                                 0.21212E+07, 0.23129E+07, 0.25181E+07, 0.27377E+07, 0.29721E+07,
                                 0.32223E+07, 0.34890E+07, 0.37729E+07, 0.40750E+07, 0.43959E+07,
                                 0.47365E+07, 0.50978E+07, 0.54807E+07, 0.58860E+07, 0.63147E+07,
                                 0.67678E+07, 0.72463E+07, 0.77512E+07, 0.82836E+07, 0.88445E+07,
                                 0.94351E+07, 0.10056E+08, 0.10710E+08, 0.11396E+08, 0.12117E+08,
                                 0.12873E+08, 0.13666E+08, 0.14497E+08, 0.15367E+08, 0.16279E+08,
                                 0.17232E+08, 0.18229E+08, 0.19271E+08, 0.20359E+08, 0.21495E+08,
                                 0.22681E+08, 0.23917E+08, 0.25206E+08, 0.26549E+08, 0.27948E+08,
                                 0.29404E+08, 0.30920E+08, 0.32496E+08, 0.34135E+08, 0.35838E+08,
                                 0.37608E+08, 0.39445E+08, 0.41353E+08, 0.43332E+08, 0.45385E+08,
                                 0.47514E+08, 0.49721E+08, 0.52007E+08, 0.54376E+08, 0.56829E+08,
                                 0.59367E+08, 0.61995E+08, 0.64712E+08, 0.67523E+08, 0.70429E+08,
                                 0.73432E+08, 0.76535E+08, 0.79740E+08, 0.83050E+08, 0.86467E+08,
                                 0.89993E+08, 0.93632E+08, 0.97385E+08, 0.10126E+09, 0.10525E+09,
                                 0.10936E+09, 0.11360E+09, 0.11796E+09, 0.12246E+09, 0.12709E+09,
                                 0.13186E+09, 0.13677E+09, 0.14182E+09, 0.14701E+09, 0.15236E+09,
                                 0.15785E+09, 0.16350E+09, 0.16931E+09, 0.17528E+09, 0.18141E+09,
                                 0.18771E+09, 0.19418E+09, 0.20082E+09, 0.20764E+09, 0.21465E+09,
                                 0.22183E+09])


#  --------------- O3 787: M = 3, I = 17 ---------------------
M = 3
I = 17
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.12211E+05, 0.20582E+05, 0.30305E+05,
                                 0.41237E+05, 0.53314E+05, 0.66536E+05, 0.80957E+05, 0.96672E+05,
                                 0.11380E+06, 0.13250E+06, 0.15292E+06, 0.17524E+06, 0.19965E+06,
                                 0.22632E+06, 0.25546E+06, 0.28728E+06, 0.32199E+06, 0.35980E+06,
                                 0.40094E+06, 0.44565E+06, 0.49417E+06, 0.54676E+06, 0.60366E+06,
                                 0.66516E+06, 0.73152E+06, 0.80305E+06, 0.88002E+06, 0.96276E+06,
                                 0.10516E+07, 0.11468E+07, 0.12488E+07, 0.13578E+07, 0.14743E+07,
                                 0.15987E+07, 0.17312E+07, 0.18723E+07, 0.20225E+07, 0.21820E+07,
                                 0.23514E+07, 0.25310E+07, 0.27214E+07, 0.29230E+07, 0.31362E+07,
                                 0.33616E+07, 0.35997E+07, 0.38509E+07, 0.41158E+07, 0.43949E+07,
                                 0.46887E+07, 0.49980E+07, 0.53231E+07, 0.56647E+07, 0.60234E+07,
                                 0.63998E+07, 0.67946E+07, 0.72084E+07, 0.76418E+07, 0.80955E+07,
                                 0.85702E+07, 0.90666E+07, 0.95854E+07, 0.10127E+08, 0.10693E+08,
                                 0.11284E+08, 0.11900E+08, 0.12542E+08, 0.13211E+08, 0.13907E+08,
                                 0.14633E+08, 0.15388E+08, 0.16173E+08, 0.16990E+08, 0.17838E+08,
                                 0.18720E+08, 0.19636E+08, 0.20586E+08, 0.21573E+08, 0.22596E+08,
                                 0.23657E+08, 0.24757E+08, 0.25896E+08, 0.27077E+08, 0.28299E+08,
                                 0.29565E+08, 0.30874E+08, 0.32229E+08, 0.33630E+08, 0.35079E+08,
                                 0.36576E+08, 0.38123E+08, 0.39721E+08, 0.41371E+08, 0.43075E+08,
                                 0.44833E+08, 0.46647E+08, 0.48518E+08, 0.50448E+08, 0.52438E+08,
                                 0.54489E+08, 0.56603E+08, 0.58780E+08, 0.61023E+08, 0.63332E+08,
                                 0.65710E+08, 0.68157E+08, 0.70676E+08, 0.73266E+08, 0.75931E+08,
                                 0.78672E+08, 0.81490E+08, 0.84386E+08, 0.87363E+08, 0.90422E+08,
                                 0.93564E+08, 0.96791E+08, 0.10011E+09, 0.10351E+09, 0.10700E+09,
                                 0.11059E+09])


#  --------------- O3 777: M = 3, I = 18 ---------------------
M = 3
I = 18
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.71750E+05, 0.12094E+06, 0.17807E+06,
                                 0.24230E+06, 0.31324E+06, 0.39088E+06, 0.47550E+06, 0.56764E+06,
                                 0.66800E+06, 0.77740E+06, 0.89677E+06, 0.10271E+07, 0.11694E+07,
                                 0.13249E+07, 0.14945E+07, 0.16796E+07, 0.18813E+07, 0.21009E+07,
                                 0.23396E+07, 0.25989E+07, 0.28801E+07, 0.31847E+07, 0.35140E+07,
                                 0.38698E+07, 0.42535E+07, 0.46669E+07, 0.51115E+07, 0.55893E+07,
                                 0.61019E+07, 0.66513E+07, 0.72393E+07, 0.78680E+07, 0.85395E+07,
                                 0.92558E+07, 0.10019E+08, 0.10832E+08, 0.11696E+08, 0.12614E+08,
                                 0.13588E+08, 0.14621E+08, 0.15716E+08, 0.16875E+08, 0.18100E+08,
                                 0.19395E+08, 0.20762E+08, 0.22205E+08, 0.23726E+08, 0.25328E+08,
                                 0.27015E+08, 0.28789E+08, 0.30654E+08, 0.32614E+08, 0.34671E+08,
                                 0.36830E+08, 0.39093E+08, 0.41465E+08, 0.43949E+08, 0.46549E+08,
                                 0.49269E+08, 0.52112E+08, 0.55084E+08, 0.58188E+08, 0.61428E+08,
                                 0.64809E+08, 0.68335E+08, 0.72010E+08, 0.75840E+08, 0.79828E+08,
                                 0.83979E+08, 0.88299E+08, 0.92792E+08, 0.97463E+08, 0.10232E+09,
                                 0.10736E+09, 0.11260E+09, 0.11803E+09, 0.12367E+09, 0.12952E+09,
                                 0.13559E+09, 0.14187E+09, 0.14839E+09, 0.15513E+09, 0.16212E+09,
                                 0.16935E+09, 0.17683E+09, 0.18457E+09, 0.19257E+09, 0.20085E+09,
                                 0.20940E+09, 0.21824E+09, 0.22736E+09, 0.23678E+09, 0.24651E+09,
                                 0.25655E+09, 0.26691E+09, 0.27759E+09, 0.28861E+09, 0.29997E+09,
                                 0.31167E+09, 0.32374E+09, 0.33616E+09, 0.34896E+09, 0.36214E+09,
                                 0.37571E+09, 0.38967E+09, 0.40404E+09, 0.41882E+09, 0.43403E+09,
                                 0.44966E+09, 0.46573E+09, 0.48226E+09, 0.49923E+09, 0.51668E+09,
                                 0.53460E+09, 0.55301E+09, 0.57191E+09, 0.59131E+09, 0.61123E+09,
                                 0.63167E+09])


#  --------------- N2O 446: M = 4, I = 1 ---------------------
M = 4
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(9.)
TIPS_ISO_HASH[(M, I)] = float32([0.89943E+03, 0.12734E+04, 0.16489E+04,
                                 0.20293E+04, 0.24205E+04, 0.28289E+04, 0.32609E+04, 0.37222E+04,
                                 0.42180E+04, 0.47529E+04, 0.53312E+04, 0.59572E+04, 0.66348E+04,
                                 0.73683E+04, 0.81616E+04, 0.90190E+04, 0.99450E+04, 0.10944E+05,
                                 0.12021E+05, 0.13180E+05, 0.14426E+05, 0.15766E+05, 0.17203E+05,
                                 0.18745E+05, 0.20396E+05, 0.22162E+05, 0.24051E+05, 0.26069E+05,
                                 0.28222E+05, 0.30517E+05, 0.32962E+05, 0.35564E+05, 0.38331E+05,
                                 0.41271E+05, 0.44393E+05, 0.47704E+05, 0.51214E+05, 0.54932E+05,
                                 0.58868E+05, 0.63030E+05, 0.67429E+05, 0.72075E+05, 0.76979E+05,
                                 0.82151E+05, 0.87604E+05, 0.93348E+05, 0.99395E+05, 0.10576E+06,
                                 0.11245E+06, 0.11948E+06, 0.12686E+06, 0.13461E+06, 0.14275E+06,
                                 0.15128E+06, 0.16021E+06, 0.16958E+06, 0.17938E+06, 0.18964E+06,
                                 0.20037E+06, 0.21159E+06, 0.22331E+06, 0.23556E+06, 0.24834E+06,
                                 0.26169E+06, 0.27561E+06, 0.29012E+06, 0.30525E+06, 0.32101E+06,
                                 0.33743E+06, 0.35452E+06, 0.37230E+06, 0.39080E+06, 0.41004E+06,
                                 0.43004E+06, 0.45082E+06, 0.47241E+06, 0.49483E+06, 0.51810E+06,
                                 0.54225E+06, 0.56730E+06, 0.59329E+06, 0.62022E+06, 0.64814E+06,
                                 0.67707E+06, 0.70703E+06, 0.73806E+06, 0.77018E+06, 0.80342E+06,
                                 0.83781E+06, 0.87338E+06, 0.91016E+06, 0.94818E+06, 0.98748E+06,
                                 0.10281E+07, 0.10700E+07, 0.11133E+07, 0.11581E+07, 0.12042E+07,
                                 0.12519E+07, 0.13010E+07, 0.13517E+07, 0.14040E+07, 0.14579E+07,
                                 0.15134E+07, 0.15707E+07, 0.16297E+07, 0.16905E+07, 0.17530E+07,
                                 0.18175E+07, 0.18838E+07, 0.19521E+07, 0.20224E+07, 0.20947E+07,
                                 0.21690E+07, 0.22455E+07, 0.23242E+07, 0.24050E+07, 0.24881E+07,
                                 0.25735E+07])


#  --------------- N2O 456: M = 4, I = 2 ---------------------
M = 4
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.59966E+03, 0.84903E+03, 0.10995E+04,
                                 0.13538E+04, 0.16158E+04, 0.18903E+04, 0.21815E+04, 0.24934E+04,
                                 0.28295E+04, 0.31927E+04, 0.35862E+04, 0.40128E+04, 0.44752E+04,
                                 0.49763E+04, 0.55189E+04, 0.61059E+04, 0.67404E+04, 0.74256E+04,
                                 0.81646E+04, 0.89609E+04, 0.98180E+04, 0.10740E+05, 0.11729E+05,
                                 0.12791E+05, 0.13930E+05, 0.15149E+05, 0.16453E+05, 0.17847E+05,
                                 0.19335E+05, 0.20922E+05, 0.22614E+05, 0.24416E+05, 0.26333E+05,
                                 0.28371E+05, 0.30535E+05, 0.32833E+05, 0.35269E+05, 0.37851E+05,
                                 0.40585E+05, 0.43478E+05, 0.46537E+05, 0.49769E+05, 0.53182E+05,
                                 0.56783E+05, 0.60580E+05, 0.64582E+05, 0.68796E+05, 0.73232E+05,
                                 0.77898E+05, 0.82803E+05, 0.87957E+05, 0.93369E+05, 0.99048E+05,
                                 0.10501E+06, 0.11125E+06, 0.11780E+06, 0.12465E+06, 0.13182E+06,
                                 0.13933E+06, 0.14718E+06, 0.15539E+06, 0.16396E+06, 0.17291E+06,
                                 0.18226E+06, 0.19201E+06, 0.20218E+06, 0.21278E+06, 0.22383E+06,
                                 0.23534E+06, 0.24733E+06, 0.25980E+06, 0.27278E+06, 0.28628E+06,
                                 0.30032E+06, 0.31491E+06, 0.33007E+06, 0.34581E+06, 0.36216E+06,
                                 0.37912E+06, 0.39673E+06, 0.41499E+06, 0.43392E+06, 0.45355E+06,
                                 0.47389E+06, 0.49496E+06, 0.51678E+06, 0.53937E+06, 0.56276E+06,
                                 0.58695E+06, 0.61199E+06, 0.63788E+06, 0.66464E+06, 0.69231E+06,
                                 0.72090E+06, 0.75044E+06, 0.78094E+06, 0.81244E+06, 0.84496E+06,
                                 0.87853E+06, 0.91316E+06, 0.94889E+06, 0.98573E+06, 0.10237E+07,
                                 0.10629E+07, 0.11033E+07, 0.11449E+07, 0.11877E+07, 0.12319E+07,
                                 0.12773E+07, 0.13241E+07, 0.13723E+07, 0.14219E+07, 0.14729E+07,
                                 0.15254E+07, 0.15793E+07, 0.16349E+07, 0.16919E+07, 0.17506E+07,
                                 0.18109E+07])


#  --------------- N2O 546: M = 4, I = 3 ---------------------
M = 4
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.62051E+03, 0.87856E+03, 0.11377E+04,
                                 0.14003E+04, 0.16705E+04, 0.19529E+04, 0.22518E+04, 0.25713E+04,
                                 0.29149E+04, 0.32859E+04, 0.36873E+04, 0.41220E+04, 0.45929E+04,
                                 0.51028E+04, 0.56547E+04, 0.62515E+04, 0.68963E+04, 0.75923E+04,
                                 0.83428E+04, 0.91511E+04, 0.10021E+05, 0.10956E+05, 0.11960E+05,
                                 0.13036E+05, 0.14190E+05, 0.15425E+05, 0.16746E+05, 0.18158E+05,
                                 0.19664E+05, 0.21271E+05, 0.22984E+05, 0.24806E+05, 0.26745E+05,
                                 0.28806E+05, 0.30995E+05, 0.33317E+05, 0.35780E+05, 0.38389E+05,
                                 0.41151E+05, 0.44073E+05, 0.47162E+05, 0.50425E+05, 0.53871E+05,
                                 0.57505E+05, 0.61338E+05, 0.65375E+05, 0.69628E+05, 0.74102E+05,
                                 0.78808E+05, 0.83755E+05, 0.88951E+05, 0.94407E+05, 0.10013E+06,
                                 0.10614E+06, 0.11243E+06, 0.11902E+06, 0.12593E+06, 0.13316E+06,
                                 0.14072E+06, 0.14862E+06, 0.15689E+06, 0.16552E+06, 0.17453E+06,
                                 0.18394E+06, 0.19376E+06, 0.20399E+06, 0.21466E+06, 0.22578E+06,
                                 0.23737E+06, 0.24942E+06, 0.26198E+06, 0.27503E+06, 0.28861E+06,
                                 0.30273E+06, 0.31741E+06, 0.33265E+06, 0.34848E+06, 0.36492E+06,
                                 0.38197E+06, 0.39967E+06, 0.41803E+06, 0.43706E+06, 0.45679E+06,
                                 0.47723E+06, 0.49840E+06, 0.52033E+06, 0.54303E+06, 0.56653E+06,
                                 0.59084E+06, 0.61599E+06, 0.64200E+06, 0.66888E+06, 0.69667E+06,
                                 0.72539E+06, 0.75506E+06, 0.78569E+06, 0.81733E+06, 0.84998E+06,
                                 0.88369E+06, 0.91846E+06, 0.95433E+06, 0.99132E+06, 0.10295E+07,
                                 0.10688E+07, 0.11093E+07, 0.11511E+07, 0.11941E+07, 0.12384E+07,
                                 0.12840E+07, 0.13310E+07, 0.13793E+07, 0.14291E+07, 0.14803E+07,
                                 0.15329E+07, 0.15871E+07, 0.16428E+07, 0.17000E+07, 0.17589E+07,
                                 0.18194E+07])


#  --------------- N2O 448: M = 4, I = 4 ---------------------
M = 4
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(9.)
TIPS_ISO_HASH[(M, I)] = float32([0.95253E+03, 0.13487E+04, 0.17465E+04,
                                 0.21498E+04, 0.25648E+04, 0.29986E+04, 0.34580E+04, 0.39493E+04,
                                 0.44779E+04, 0.50488E+04, 0.56669E+04, 0.63366E+04, 0.70625E+04,
                                 0.78488E+04, 0.87003E+04, 0.96216E+04, 0.10617E+05, 0.11692E+05,
                                 0.12852E+05, 0.14102E+05, 0.15447E+05, 0.16893E+05, 0.18446E+05,
                                 0.20112E+05, 0.21898E+05, 0.23811E+05, 0.25856E+05, 0.28042E+05,
                                 0.30377E+05, 0.32866E+05, 0.35520E+05, 0.38345E+05, 0.41351E+05,
                                 0.44545E+05, 0.47939E+05, 0.51540E+05, 0.55359E+05, 0.59405E+05,
                                 0.63689E+05, 0.68222E+05, 0.73015E+05, 0.78078E+05, 0.83424E+05,
                                 0.89064E+05, 0.95012E+05, 0.10128E+06, 0.10788E+06, 0.11482E+06,
                                 0.12213E+06, 0.12981E+06, 0.13788E+06, 0.14635E+06, 0.15524E+06,
                                 0.16456E+06, 0.17433E+06, 0.18457E+06, 0.19530E+06, 0.20652E+06,
                                 0.21827E+06, 0.23055E+06, 0.24338E+06, 0.25679E+06, 0.27079E+06,
                                 0.28541E+06, 0.30066E+06, 0.31656E+06, 0.33314E+06, 0.35042E+06,
                                 0.36841E+06, 0.38715E+06, 0.40666E+06, 0.42695E+06, 0.44805E+06,
                                 0.46999E+06, 0.49279E+06, 0.51649E+06, 0.54109E+06, 0.56664E+06,
                                 0.59315E+06, 0.62066E+06, 0.64919E+06, 0.67877E+06, 0.70943E+06,
                                 0.74121E+06, 0.77413E+06, 0.80822E+06, 0.84351E+06, 0.88004E+06,
                                 0.91783E+06, 0.95693E+06, 0.99737E+06, 0.10392E+07, 0.10824E+07,
                                 0.11270E+07, 0.11732E+07, 0.12208E+07, 0.12700E+07, 0.13208E+07,
                                 0.13732E+07, 0.14272E+07, 0.14830E+07, 0.15405E+07, 0.15999E+07,
                                 0.16610E+07, 0.17240E+07, 0.17890E+07, 0.18559E+07, 0.19248E+07,
                                 0.19957E+07, 0.20687E+07, 0.21439E+07, 0.22213E+07, 0.23009E+07,
                                 0.23828E+07, 0.24671E+07, 0.25537E+07, 0.26428E+07, 0.27343E+07,
                                 0.28284E+07])


#  --------------- N2O 447: M = 4, I = 5 ---------------------
M = 4
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(54.)
TIPS_ISO_HASH[(M, I)] = float32([0.55598E+04, 0.78718E+04, 0.10193E+05,
                                 0.12546E+05, 0.14966E+05, 0.17495E+05, 0.20171E+05, 0.23031E+05,
                                 0.26106E+05, 0.29426E+05, 0.33018E+05, 0.36908E+05, 0.41121E+05,
                                 0.45684E+05, 0.50622E+05, 0.55962E+05, 0.61731E+05, 0.67958E+05,
                                 0.74671E+05, 0.81902E+05, 0.89681E+05, 0.98043E+05, 0.10702E+06,
                                 0.11665E+06, 0.12697E+06, 0.13801E+06, 0.14983E+06, 0.16244E+06,
                                 0.17591E+06, 0.19028E+06, 0.20558E+06, 0.22188E+06, 0.23920E+06,
                                 0.25762E+06, 0.27718E+06, 0.29793E+06, 0.31993E+06, 0.34323E+06,
                                 0.36791E+06, 0.39401E+06, 0.42160E+06, 0.45074E+06, 0.48151E+06,
                                 0.51397E+06, 0.54819E+06, 0.58424E+06, 0.62221E+06, 0.66215E+06,
                                 0.70416E+06, 0.74832E+06, 0.79470E+06, 0.84340E+06, 0.89450E+06,
                                 0.94808E+06, 0.10042E+07, 0.10631E+07, 0.11247E+07, 0.11892E+07,
                                 0.12567E+07, 0.13272E+07, 0.14009E+07, 0.14779E+07, 0.15583E+07,
                                 0.16422E+07, 0.17298E+07, 0.18211E+07, 0.19163E+07, 0.20154E+07,
                                 0.21187E+07, 0.22263E+07, 0.23382E+07, 0.24546E+07, 0.25757E+07,
                                 0.27016E+07, 0.28324E+07, 0.29683E+07, 0.31095E+07, 0.32560E+07,
                                 0.34081E+07, 0.35659E+07, 0.37295E+07, 0.38991E+07, 0.40750E+07,
                                 0.42572E+07, 0.44459E+07, 0.46414E+07, 0.48437E+07, 0.50531E+07,
                                 0.52698E+07, 0.54939E+07, 0.57257E+07, 0.59653E+07, 0.62129E+07,
                                 0.64688E+07, 0.67331E+07, 0.70061E+07, 0.72880E+07, 0.75790E+07,
                                 0.78792E+07, 0.81891E+07, 0.85086E+07, 0.88382E+07, 0.91780E+07,
                                 0.95283E+07, 0.98893E+07, 0.10261E+08, 0.10644E+08, 0.11039E+08,
                                 0.11445E+08, 0.11864E+08, 0.12294E+08, 0.12738E+08, 0.13194E+08,
                                 0.13663E+08, 0.14145E+08, 0.14641E+08, 0.15151E+08, 0.15675E+08,
                                 0.16214E+08])


#  --------------- CO 26: M = 5, I = 1 ---------------------
M = 5
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.21948E+02, 0.30961E+02, 0.39980E+02,
                                 0.49004E+02, 0.58035E+02, 0.67071E+02, 0.76112E+02, 0.85160E+02,
                                 0.94213E+02, 0.10327E+03, 0.11234E+03, 0.12142E+03, 0.13050E+03,
                                 0.13960E+03, 0.14872E+03, 0.15787E+03, 0.16704E+03, 0.17624E+03,
                                 0.18548E+03, 0.19477E+03, 0.20411E+03, 0.21350E+03, 0.22295E+03,
                                 0.23248E+03, 0.24207E+03, 0.25175E+03, 0.26151E+03, 0.27136E+03,
                                 0.28130E+03, 0.29134E+03, 0.30148E+03, 0.31172E+03, 0.32207E+03,
                                 0.33253E+03, 0.34312E+03, 0.35381E+03, 0.36463E+03, 0.37557E+03,
                                 0.38663E+03, 0.39782E+03, 0.40914E+03, 0.42060E+03, 0.43218E+03,
                                 0.44389E+03, 0.45575E+03, 0.46774E+03, 0.47987E+03, 0.49213E+03,
                                 0.50454E+03, 0.51708E+03, 0.52978E+03, 0.54261E+03, 0.55559E+03,
                                 0.56871E+03, 0.58198E+03, 0.59540E+03, 0.60896E+03, 0.62267E+03,
                                 0.63653E+03, 0.65055E+03, 0.66470E+03, 0.67901E+03, 0.69347E+03,
                                 0.70808E+03, 0.72284E+03, 0.73776E+03, 0.75283E+03, 0.76805E+03,
                                 0.78342E+03, 0.79895E+03, 0.81463E+03, 0.83047E+03, 0.84646E+03,
                                 0.86260E+03, 0.87891E+03, 0.89536E+03, 0.91197E+03, 0.92874E+03,
                                 0.94566E+03, 0.96275E+03, 0.97998E+03, 0.99738E+03, 0.10149E+04,
                                 0.10326E+04, 0.10505E+04, 0.10685E+04, 0.10867E+04, 0.11051E+04,
                                 0.11236E+04, 0.11422E+04, 0.11611E+04, 0.11800E+04, 0.11992E+04,
                                 0.12185E+04, 0.12380E+04, 0.12576E+04, 0.12774E+04, 0.12973E+04,
                                 0.13174E+04, 0.13377E+04, 0.13581E+04, 0.13787E+04, 0.13994E+04,
                                 0.14203E+04, 0.14414E+04, 0.14627E+04, 0.14841E+04, 0.15056E+04,
                                 0.15273E+04, 0.15492E+04, 0.15713E+04, 0.15935E+04, 0.16159E+04,
                                 0.16384E+04, 0.16611E+04, 0.16840E+04, 0.17070E+04, 0.17302E+04,
                                 0.17536E+04])


#  --------------- CO 36: M = 5, I = 2 ---------------------
M = 5
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.45888E+02, 0.64745E+02, 0.83615E+02,
                                 0.10250E+03, 0.12139E+03, 0.14030E+03, 0.15921E+03, 0.17814E+03,
                                 0.19708E+03, 0.21604E+03, 0.23501E+03, 0.25400E+03, 0.27302E+03,
                                 0.29207E+03, 0.31117E+03, 0.33031E+03, 0.34952E+03, 0.36880E+03,
                                 0.38817E+03, 0.40764E+03, 0.42723E+03, 0.44694E+03, 0.46679E+03,
                                 0.48679E+03, 0.50696E+03, 0.52730E+03, 0.54783E+03, 0.56855E+03,
                                 0.58948E+03, 0.61061E+03, 0.63198E+03, 0.65357E+03, 0.67539E+03,
                                 0.69747E+03, 0.71979E+03, 0.74237E+03, 0.76521E+03, 0.78832E+03,
                                 0.81169E+03, 0.83534E+03, 0.85927E+03, 0.88348E+03, 0.90798E+03,
                                 0.93277E+03, 0.95784E+03, 0.98322E+03, 0.10089E+04, 0.10349E+04,
                                 0.10611E+04, 0.10877E+04, 0.11146E+04, 0.11418E+04, 0.11693E+04,
                                 0.11971E+04, 0.12253E+04, 0.12537E+04, 0.12825E+04, 0.13115E+04,
                                 0.13409E+04, 0.13707E+04, 0.14007E+04, 0.14311E+04, 0.14617E+04,
                                 0.14928E+04, 0.15241E+04, 0.15558E+04, 0.15877E+04, 0.16200E+04,
                                 0.16527E+04, 0.16857E+04, 0.17190E+04, 0.17526E+04, 0.17866E+04,
                                 0.18209E+04, 0.18555E+04, 0.18905E+04, 0.19258E+04, 0.19614E+04,
                                 0.19974E+04, 0.20337E+04, 0.20703E+04, 0.21073E+04, 0.21446E+04,
                                 0.21823E+04, 0.22203E+04, 0.22586E+04, 0.22973E+04, 0.23363E+04,
                                 0.23756E+04, 0.24153E+04, 0.24553E+04, 0.24957E+04, 0.25364E+04,
                                 0.25775E+04, 0.26189E+04, 0.26606E+04, 0.27027E+04, 0.27451E+04,
                                 0.27879E+04, 0.28310E+04, 0.28745E+04, 0.29183E+04, 0.29625E+04,
                                 0.30070E+04, 0.30518E+04, 0.30970E+04, 0.31425E+04, 0.31885E+04,
                                 0.32347E+04, 0.32813E+04, 0.33282E+04, 0.33755E+04, 0.34231E+04,
                                 0.34711E+04, 0.35194E+04, 0.35681E+04, 0.36172E+04, 0.36666E+04,
                                 0.37163E+04])


#  --------------- CO 28: M = 5, I = 3 ---------------------
M = 5
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.23030E+02, 0.32495E+02, 0.41966E+02,
                                 0.51443E+02, 0.60926E+02, 0.70415E+02, 0.79910E+02, 0.89410E+02,
                                 0.98918E+02, 0.10843E+03, 0.11795E+03, 0.12749E+03, 0.13703E+03,
                                 0.14659E+03, 0.15618E+03, 0.16579E+03, 0.17543E+03, 0.18511E+03,
                                 0.19483E+03, 0.20461E+03, 0.21444E+03, 0.22434E+03, 0.23430E+03,
                                 0.24435E+03, 0.25447E+03, 0.26468E+03, 0.27499E+03, 0.28540E+03,
                                 0.29591E+03, 0.30652E+03, 0.31725E+03, 0.32810E+03, 0.33906E+03,
                                 0.35014E+03, 0.36136E+03, 0.37270E+03, 0.38417E+03, 0.39577E+03,
                                 0.40752E+03, 0.41940E+03, 0.43142E+03, 0.44358E+03, 0.45589E+03,
                                 0.46834E+03, 0.48094E+03, 0.49369E+03, 0.50659E+03, 0.51964E+03,
                                 0.53284E+03, 0.54619E+03, 0.55971E+03, 0.57337E+03, 0.58719E+03,
                                 0.60117E+03, 0.61530E+03, 0.62959E+03, 0.64405E+03, 0.65866E+03,
                                 0.67343E+03, 0.68837E+03, 0.70346E+03, 0.71872E+03, 0.73414E+03,
                                 0.74972E+03, 0.76547E+03, 0.78138E+03, 0.79745E+03, 0.81369E+03,
                                 0.83010E+03, 0.84667E+03, 0.86341E+03, 0.88031E+03, 0.89738E+03,
                                 0.91462E+03, 0.93202E+03, 0.94960E+03, 0.96734E+03, 0.98524E+03,
                                 0.10033E+04, 0.10216E+04, 0.10400E+04, 0.10586E+04, 0.10773E+04,
                                 0.10962E+04, 0.11153E+04, 0.11346E+04, 0.11540E+04, 0.11737E+04,
                                 0.11934E+04, 0.12134E+04, 0.12335E+04, 0.12538E+04, 0.12743E+04,
                                 0.12949E+04, 0.13157E+04, 0.13367E+04, 0.13578E+04, 0.13792E+04,
                                 0.14007E+04, 0.14223E+04, 0.14442E+04, 0.14662E+04, 0.14884E+04,
                                 0.15108E+04, 0.15333E+04, 0.15560E+04, 0.15789E+04, 0.16020E+04,
                                 0.16252E+04, 0.16486E+04, 0.16722E+04, 0.16960E+04, 0.17199E+04,
                                 0.17441E+04, 0.17684E+04, 0.17928E+04, 0.18175E+04, 0.18423E+04,
                                 0.18673E+04])


#  --------------- CO 27: M = 5, I = 4 ---------------------
M = 5
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.13505E+03, 0.19054E+03, 0.24606E+03,
                                 0.30161E+03, 0.35720E+03, 0.41283E+03, 0.46848E+03, 0.52418E+03,
                                 0.57991E+03, 0.63568E+03, 0.69149E+03, 0.74737E+03, 0.80332E+03,
                                 0.85937E+03, 0.91553E+03, 0.97183E+03, 0.10283E+04, 0.10850E+04,
                                 0.11420E+04, 0.11992E+04, 0.12568E+04, 0.13147E+04, 0.13730E+04,
                                 0.14318E+04, 0.14910E+04, 0.15507E+04, 0.16110E+04, 0.16718E+04,
                                 0.17332E+04, 0.17952E+04, 0.18579E+04, 0.19212E+04, 0.19852E+04,
                                 0.20499E+04, 0.21153E+04, 0.21815E+04, 0.22484E+04, 0.23161E+04,
                                 0.23846E+04, 0.24539E+04, 0.25240E+04, 0.25949E+04, 0.26666E+04,
                                 0.27392E+04, 0.28127E+04, 0.28869E+04, 0.29621E+04, 0.30381E+04,
                                 0.31150E+04, 0.31928E+04, 0.32715E+04, 0.33511E+04, 0.34316E+04,
                                 0.35129E+04, 0.35952E+04, 0.36785E+04, 0.37626E+04, 0.38477E+04,
                                 0.39336E+04, 0.40206E+04, 0.41084E+04, 0.41972E+04, 0.42869E+04,
                                 0.43776E+04, 0.44692E+04, 0.45618E+04, 0.46553E+04, 0.47498E+04,
                                 0.48452E+04, 0.49416E+04, 0.50390E+04, 0.51373E+04, 0.52366E+04,
                                 0.53368E+04, 0.54381E+04, 0.55403E+04, 0.56435E+04, 0.57476E+04,
                                 0.58527E+04, 0.59588E+04, 0.60659E+04, 0.61739E+04, 0.62829E+04,
                                 0.63930E+04, 0.65040E+04, 0.66160E+04, 0.67290E+04, 0.68429E+04,
                                 0.69579E+04, 0.70739E+04, 0.71908E+04, 0.73088E+04, 0.74277E+04,
                                 0.75477E+04, 0.76686E+04, 0.77905E+04, 0.79135E+04, 0.80374E+04,
                                 0.81624E+04, 0.82883E+04, 0.84153E+04, 0.85432E+04, 0.86722E+04,
                                 0.88022E+04, 0.89331E+04, 0.90651E+04, 0.91982E+04, 0.93322E+04,
                                 0.94672E+04, 0.96033E+04, 0.97404E+04, 0.98785E+04, 0.10018E+05,
                                 0.10158E+05, 0.10299E+05, 0.10441E+05, 0.10584E+05, 0.10728E+05,
                                 0.10874E+05])


#  --------------- CO 38: M = 5, I = 5 ---------------------
M = 5
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.48264E+02, 0.68112E+02, 0.87974E+02,
                                 0.10785E+03, 0.12773E+03, 0.14763E+03, 0.16754E+03, 0.18747E+03,
                                 0.20741E+03, 0.22736E+03, 0.24733E+03, 0.26732E+03, 0.28735E+03,
                                 0.30741E+03, 0.32752E+03, 0.34770E+03, 0.36794E+03, 0.38828E+03,
                                 0.40871E+03, 0.42926E+03, 0.44994E+03, 0.47077E+03, 0.49175E+03,
                                 0.51290E+03, 0.53424E+03, 0.55578E+03, 0.57752E+03, 0.59948E+03,
                                 0.62166E+03, 0.64409E+03, 0.66676E+03, 0.68969E+03, 0.71287E+03,
                                 0.73633E+03, 0.76006E+03, 0.78407E+03, 0.80836E+03, 0.83295E+03,
                                 0.85784E+03, 0.88302E+03, 0.90851E+03, 0.93431E+03, 0.96042E+03,
                                 0.98686E+03, 0.10136E+04, 0.10407E+04, 0.10681E+04, 0.10958E+04,
                                 0.11238E+04, 0.11522E+04, 0.11809E+04, 0.12100E+04, 0.12393E+04,
                                 0.12691E+04, 0.12991E+04, 0.13295E+04, 0.13603E+04, 0.13914E+04,
                                 0.14228E+04, 0.14546E+04, 0.14867E+04, 0.15192E+04, 0.15520E+04,
                                 0.15852E+04, 0.16187E+04, 0.16526E+04, 0.16869E+04, 0.17215E+04,
                                 0.17564E+04, 0.17917E+04, 0.18274E+04, 0.18634E+04, 0.18998E+04,
                                 0.19365E+04, 0.19736E+04, 0.20111E+04, 0.20489E+04, 0.20871E+04,
                                 0.21256E+04, 0.21645E+04, 0.22038E+04, 0.22434E+04, 0.22834E+04,
                                 0.23238E+04, 0.23645E+04, 0.24056E+04, 0.24471E+04, 0.24889E+04,
                                 0.25311E+04, 0.25736E+04, 0.26166E+04, 0.26599E+04, 0.27035E+04,
                                 0.27476E+04, 0.27920E+04, 0.28368E+04, 0.28819E+04, 0.29275E+04,
                                 0.29733E+04, 0.30196E+04, 0.30662E+04, 0.31133E+04, 0.31606E+04,
                                 0.32084E+04, 0.32565E+04, 0.33050E+04, 0.33539E+04, 0.34032E+04,
                                 0.34528E+04, 0.35028E+04, 0.35532E+04, 0.36040E+04, 0.36551E+04,
                                 0.37067E+04, 0.37586E+04, 0.38108E+04, 0.38635E+04, 0.39165E+04,
                                 0.39699E+04])


#  --------------- CO 37: M = 5, I = 6 ---------------------
M = 5
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.28271E+03, 0.39894E+03, 0.51524E+03,
                                 0.63162E+03, 0.74807E+03, 0.86459E+03, 0.98119E+03, 0.10979E+04,
                                 0.12146E+04, 0.13314E+04, 0.14484E+04, 0.15654E+04, 0.16826E+04,
                                 0.18000E+04, 0.19176E+04, 0.20355E+04, 0.21538E+04, 0.22725E+04,
                                 0.23916E+04, 0.25114E+04, 0.26318E+04, 0.27529E+04, 0.28749E+04,
                                 0.29977E+04, 0.31215E+04, 0.32463E+04, 0.33721E+04, 0.34991E+04,
                                 0.36274E+04, 0.37568E+04, 0.38876E+04, 0.40197E+04, 0.41533E+04,
                                 0.42882E+04, 0.44247E+04, 0.45626E+04, 0.47022E+04, 0.48433E+04,
                                 0.49860E+04, 0.51304E+04, 0.52763E+04, 0.54240E+04, 0.55735E+04,
                                 0.57246E+04, 0.58775E+04, 0.60321E+04, 0.61886E+04, 0.63468E+04,
                                 0.65068E+04, 0.66687E+04, 0.68324E+04, 0.69980E+04, 0.71654E+04,
                                 0.73347E+04, 0.75058E+04, 0.76789E+04, 0.78539E+04, 0.80307E+04,
                                 0.82096E+04, 0.83903E+04, 0.85729E+04, 0.87576E+04, 0.89441E+04,
                                 0.91326E+04, 0.93230E+04, 0.95154E+04, 0.97098E+04, 0.99061E+04,
                                 0.10104E+05, 0.10305E+05, 0.10507E+05, 0.10711E+05, 0.10918E+05,
                                 0.11126E+05, 0.11336E+05, 0.11549E+05, 0.11763E+05, 0.11979E+05,
                                 0.12198E+05, 0.12418E+05, 0.12640E+05, 0.12865E+05, 0.13091E+05,
                                 0.13320E+05, 0.13550E+05, 0.13783E+05, 0.14018E+05, 0.14254E+05,
                                 0.14493E+05, 0.14734E+05, 0.14977E+05, 0.15221E+05, 0.15468E+05,
                                 0.15718E+05, 0.15969E+05, 0.16222E+05, 0.16477E+05, 0.16734E+05,
                                 0.16994E+05, 0.17255E+05, 0.17519E+05, 0.17784E+05, 0.18052E+05,
                                 0.18322E+05, 0.18594E+05, 0.18868E+05, 0.19144E+05, 0.19422E+05,
                                 0.19703E+05, 0.19985E+05, 0.20270E+05, 0.20556E+05, 0.20845E+05,
                                 0.21136E+05, 0.21429E+05, 0.21724E+05, 0.22021E+05, 0.22320E+05,
                                 0.22622E+05])


#  --------------- CH4 211: M = 6, I = 1 ---------------------
M = 6
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.54800E+02, 0.91500E+02, 0.13410E+03,
                                 0.18180E+03, 0.23410E+03, 0.29070E+03, 0.35140E+03, 0.41600E+03,
                                 0.48450E+03, 0.55720E+03, 0.63420E+03, 0.71600E+03, 0.80310E+03,
                                 0.89590E+03, 0.99520E+03, 0.11017E+04, 0.12161E+04, 0.13393E+04,
                                 0.14721E+04, 0.16155E+04, 0.17706E+04, 0.19384E+04, 0.21202E+04,
                                 0.23172E+04, 0.25307E+04, 0.27624E+04, 0.30137E+04, 0.32864E+04,
                                 0.35823E+04, 0.39034E+04, 0.42519E+04, 0.46300E+04, 0.50402E+04,
                                 0.54853E+04, 0.59679E+04, 0.64913E+04, 0.70588E+04, 0.76739E+04,
                                 0.83404E+04, 0.90625E+04, 0.98446E+04, 0.10691E+05, 0.11608E+05,
                                 0.12600E+05, 0.13674E+05, 0.14835E+05, 0.16090E+05, 0.17447E+05,
                                 0.18914E+05, 0.20500E+05, 0.22212E+05, 0.24063E+05, 0.26061E+05,
                                 0.28218E+05, 0.30548E+05, 0.33063E+05, 0.35778E+05, 0.38708E+05,
                                 0.41871E+05, 0.45284E+05, 0.48970E+05, 0.52940E+05, 0.57230E+05,
                                 0.61860E+05, 0.66860E+05, 0.72250E+05, 0.78070E+05, 0.84350E+05,
                                 0.91130E+05, 0.98450E+05, 0.10635E+06, 0.11488E+06, 0.12408E+06,
                                 0.13403E+06, 0.14480E+06, 0.15640E+06, 0.16890E+06, 0.18240E+06,
                                 0.19700E+06, 0.21280E+06, 0.22980E+06, 0.24830E+06, 0.26820E+06,
                                 0.28970E+06, 0.31290E+06, 0.33800E+06, 0.36520E+06, 0.39450E+06,
                                 0.42600E+06, 0.46000E+06, 0.49700E+06, 0.53700E+06, 0.58100E+06,
                                 0.62700E+06, 0.67800E+06, 0.73300E+06, 0.79200E+06, 0.85600E+06,
                                 0.92500E+06, 0.10000E+07, 0.10800E+07, 0.11670E+07, 0.12610E+07,
                                 0.13620E+07, 0.14720E+07, 0.15910E+07, 0.17190E+07, 0.18600E+07,
                                 0.20100E+07, 0.21700E+07, 0.23400E+07, 0.25300E+07, 0.27300E+07,
                                 0.29500E+07, 0.31800E+07, 0.34300E+07, 0.37000E+07, 0.39900E+07,
                                 0.42856E+07])


#  --------------- CH4 311: M = 6, I = 2 ---------------------
M = 6
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.10958E+03, 0.18304E+03, 0.26818E+03,
                                 0.36356E+03, 0.46820E+03, 0.58141E+03, 0.70270E+03, 0.83186E+03,
                                 0.96893E+03, 0.11142E+04, 0.12682E+04, 0.14316E+04, 0.16055E+04,
                                 0.17909E+04, 0.19891E+04, 0.22016E+04, 0.24297E+04, 0.26752E+04,
                                 0.29399E+04, 0.32255E+04, 0.35342E+04, 0.38680E+04, 0.42294E+04,
                                 0.46208E+04, 0.50449E+04, 0.55046E+04, 0.60030E+04, 0.65434E+04,
                                 0.71293E+04, 0.77646E+04, 0.84535E+04, 0.92004E+04, 0.10010E+05,
                                 0.10888E+05, 0.11838E+05, 0.12869E+05, 0.13984E+05, 0.15193E+05,
                                 0.16501E+05, 0.17916E+05, 0.19448E+05, 0.21104E+05, 0.22895E+05,
                                 0.24830E+05, 0.26921E+05, 0.29180E+05, 0.31618E+05, 0.34250E+05,
                                 0.37090E+05, 0.40152E+05, 0.43454E+05, 0.47012E+05, 0.50845E+05,
                                 0.54973E+05, 0.59416E+05, 0.64197E+05, 0.69340E+05, 0.74870E+05,
                                 0.80813E+05, 0.87198E+05, 0.94055E+05, 0.10142E+06, 0.10932E+06,
                                 0.11779E+06, 0.12688E+06, 0.13662E+06, 0.14706E+06, 0.15824E+06,
                                 0.17021E+06, 0.18302E+06, 0.19673E+06, 0.21139E+06, 0.22706E+06,
                                 0.24381E+06, 0.26171E+06, 0.28082E+06, 0.30122E+06, 0.32299E+06,
                                 0.34621E+06, 0.37097E+06, 0.39737E+06, 0.42551E+06, 0.45548E+06,
                                 0.48739E+06, 0.52136E+06, 0.55752E+06, 0.59598E+06, 0.63688E+06,
                                 0.68036E+06, 0.72657E+06, 0.77566E+06, 0.82780E+06, 0.88316E+06,
                                 0.94191E+06, 0.10043E+07, 0.10704E+07, 0.11405E+07, 0.12148E+07,
                                 0.12936E+07, 0.13770E+07, 0.14654E+07, 0.15589E+07, 0.16579E+07,
                                 0.17627E+07, 0.18736E+07, 0.19908E+07, 0.21147E+07, 0.22456E+07,
                                 0.23840E+07, 0.25301E+07, 0.26844E+07, 0.28474E+07, 0.30193E+07,
                                 0.32007E+07, 0.33921E+07, 0.35939E+07, 0.38067E+07, 0.40310E+07,
                                 0.42673E+07])


#  --------------- CH4 212: M = 6, I = 3 ---------------------
M = 6
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.)
TIPS_ISO_HASH[(M, I)] = float32([0.44079E+03, 0.73786E+03, 0.10822E+04,
                                 0.14679E+04, 0.18913E+04, 0.23497E+04, 0.28415E+04, 0.33665E+04,
                                 0.39257E+04, 0.45211E+04, 0.51562E+04, 0.58349E+04, 0.65624E+04,
                                 0.73445E+04, 0.81872E+04, 0.90978E+04, 0.10084E+05, 0.11153E+05,
                                 0.12315E+05, 0.13579E+05, 0.14955E+05, 0.16455E+05, 0.18089E+05,
                                 0.19871E+05, 0.21816E+05, 0.23937E+05, 0.26251E+05, 0.28776E+05,
                                 0.31531E+05, 0.34535E+05, 0.37811E+05, 0.41384E+05, 0.45278E+05,
                                 0.49521E+05, 0.54144E+05, 0.59178E+05, 0.64657E+05, 0.70621E+05,
                                 0.77108E+05, 0.84161E+05, 0.91828E+05, 0.10016E+06, 0.10921E+06,
                                 0.11903E+06, 0.12968E+06, 0.14124E+06, 0.15378E+06, 0.16736E+06,
                                 0.18207E+06, 0.19800E+06, 0.21524E+06, 0.23389E+06, 0.25405E+06,
                                 0.27585E+06, 0.29939E+06, 0.32482E+06, 0.35226E+06, 0.38186E+06,
                                 0.41379E+06, 0.44821E+06, 0.48529E+06, 0.52522E+06, 0.56821E+06,
                                 0.61447E+06, 0.66422E+06, 0.71771E+06, 0.77519E+06, 0.83693E+06,
                                 0.90323E+06, 0.97438E+06, 0.10507E+07, 0.11326E+07, 0.12203E+07,
                                 0.13143E+07, 0.14150E+07, 0.15228E+07, 0.16382E+07, 0.17616E+07,
                                 0.18935E+07, 0.20346E+07, 0.21853E+07, 0.23463E+07, 0.25181E+07,
                                 0.27016E+07, 0.28973E+07, 0.31060E+07, 0.33284E+07, 0.35655E+07,
                                 0.38181E+07, 0.40870E+07, 0.43733E+07, 0.46780E+07, 0.50020E+07,
                                 0.53467E+07, 0.57130E+07, 0.61023E+07, 0.65158E+07, 0.69549E+07,
                                 0.74211E+07, 0.79158E+07, 0.84407E+07, 0.89973E+07, 0.95874E+07,
                                 0.10213E+08, 0.10875E+08, 0.11577E+08, 0.12320E+08, 0.13107E+08,
                                 0.13940E+08, 0.14820E+08, 0.15752E+08, 0.16736E+08, 0.17777E+08,
                                 0.18877E+08, 0.20038E+08, 0.21265E+08, 0.22560E+08, 0.23927E+08,
                                 0.25369E+08])


#  --------------- CH4 312: M = 6, I = 4 ---------------------
M = 6
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.88231E+03, 0.14770E+04, 0.21661E+04,
                                 0.29384E+04, 0.37859E+04, 0.47034E+04, 0.56879E+04, 0.67388E+04,
                                 0.78581E+04, 0.90501E+04, 0.10321E+05, 0.11680E+05, 0.13136E+05,
                                 0.14702E+05, 0.16389E+05, 0.18212E+05, 0.20186E+05, 0.22328E+05,
                                 0.24654E+05, 0.27185E+05, 0.29941E+05, 0.32943E+05, 0.36216E+05,
                                 0.39786E+05, 0.43681E+05, 0.47930E+05, 0.52567E+05, 0.57625E+05,
                                 0.63144E+05, 0.69164E+05, 0.75730E+05, 0.82890E+05, 0.90693E+05,
                                 0.99198E+05, 0.10846E+06, 0.11855E+06, 0.12954E+06, 0.14149E+06,
                                 0.15450E+06, 0.16864E+06, 0.18402E+06, 0.20072E+06, 0.21886E+06,
                                 0.23856E+06, 0.25993E+06, 0.28312E+06, 0.30825E+06, 0.33550E+06,
                                 0.36501E+06, 0.39696E+06, 0.43155E+06, 0.46896E+06, 0.50942E+06,
                                 0.55315E+06, 0.60039E+06, 0.65141E+06, 0.70648E+06, 0.76589E+06,
                                 0.82997E+06, 0.89904E+06, 0.97346E+06, 0.10536E+07, 0.11399E+07,
                                 0.12327E+07, 0.13326E+07, 0.14400E+07, 0.15554E+07, 0.16793E+07,
                                 0.18124E+07, 0.19553E+07, 0.21085E+07, 0.22729E+07, 0.24490E+07,
                                 0.26378E+07, 0.28400E+07, 0.30565E+07, 0.32881E+07, 0.35360E+07,
                                 0.38010E+07, 0.40843E+07, 0.43870E+07, 0.47103E+07, 0.50555E+07,
                                 0.54239E+07, 0.58169E+07, 0.62361E+07, 0.66830E+07, 0.71592E+07,
                                 0.76666E+07, 0.82069E+07, 0.87820E+07, 0.93940E+07, 0.10045E+08,
                                 0.10737E+08, 0.11473E+08, 0.12256E+08, 0.13086E+08, 0.13969E+08,
                                 0.14905E+08, 0.15899E+08, 0.16954E+08, 0.18072E+08, 0.19258E+08,
                                 0.20515E+08, 0.21847E+08, 0.23257E+08, 0.24750E+08, 0.26331E+08,
                                 0.28004E+08, 0.29774E+08, 0.31646E+08, 0.33625E+08, 0.35716E+08,
                                 0.37926E+08, 0.40261E+08, 0.42726E+08, 0.45329E+08, 0.48077E+08,
                                 0.50975E+08])


#  --------------- O2 66: M = 7, I = 1 ---------------------
M = 7
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.44334E+02, 0.62460E+02, 0.80596E+02,
                                 0.98738E+02, 0.11688E+03, 0.13503E+03, 0.15319E+03, 0.17136E+03,
                                 0.18954E+03, 0.20775E+03, 0.22600E+03, 0.24431E+03, 0.26270E+03,
                                 0.28119E+03, 0.29981E+03, 0.31857E+03, 0.33750E+03, 0.35662E+03,
                                 0.37594E+03, 0.39550E+03, 0.41529E+03, 0.43535E+03, 0.45568E+03,
                                 0.47630E+03, 0.49722E+03, 0.51844E+03, 0.53998E+03, 0.56185E+03,
                                 0.58406E+03, 0.60660E+03, 0.62949E+03, 0.65274E+03, 0.67635E+03,
                                 0.70031E+03, 0.72465E+03, 0.74936E+03, 0.77444E+03, 0.79990E+03,
                                 0.82574E+03, 0.85197E+03, 0.87858E+03, 0.90558E+03, 0.93297E+03,
                                 0.96076E+03, 0.98895E+03, 0.10175E+04, 0.10465E+04, 0.10759E+04,
                                 0.11057E+04, 0.11359E+04, 0.11665E+04, 0.11976E+04, 0.12290E+04,
                                 0.12609E+04, 0.12931E+04, 0.13258E+04, 0.13590E+04, 0.13925E+04,
                                 0.14265E+04, 0.14609E+04, 0.14958E+04, 0.15311E+04, 0.15669E+04,
                                 0.16031E+04, 0.16397E+04, 0.16768E+04, 0.17144E+04, 0.17524E+04,
                                 0.17909E+04, 0.18298E+04, 0.18692E+04, 0.19091E+04, 0.19495E+04,
                                 0.19904E+04, 0.20318E+04, 0.20736E+04, 0.21160E+04, 0.21588E+04,
                                 0.22022E+04, 0.22461E+04, 0.22905E+04, 0.23354E+04, 0.23809E+04,
                                 0.24268E+04, 0.24734E+04, 0.25204E+04, 0.25680E+04, 0.26162E+04,
                                 0.26649E+04, 0.27142E+04, 0.27641E+04, 0.28145E+04, 0.28655E+04,
                                 0.29171E+04, 0.29693E+04, 0.30221E+04, 0.30755E+04, 0.31295E+04,
                                 0.31841E+04, 0.32393E+04, 0.32951E+04, 0.33516E+04, 0.34087E+04,
                                 0.34665E+04, 0.35249E+04, 0.35839E+04, 0.36436E+04, 0.37040E+04,
                                 0.37650E+04, 0.38267E+04, 0.38891E+04, 0.39522E+04, 0.40159E+04,
                                 0.40804E+04, 0.41455E+04, 0.42114E+04, 0.42780E+04, 0.43452E+04,
                                 0.44132E+04])


#  --------------- O2 68: M = 7, I = 2 ---------------------
M = 7
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.89206E+02, 0.12759E+03, 0.16600E+03,
                                 0.20442E+03, 0.24285E+03, 0.28128E+03, 0.31973E+03, 0.35821E+03,
                                 0.39672E+03, 0.43530E+03, 0.47398E+03, 0.51281E+03, 0.55183E+03,
                                 0.59108E+03, 0.63062E+03, 0.67051E+03, 0.71078E+03, 0.75148E+03,
                                 0.79265E+03, 0.83435E+03, 0.87659E+03, 0.91941E+03, 0.96285E+03,
                                 0.10069E+04, 0.10517E+04, 0.10971E+04, 0.11432E+04, 0.11901E+04,
                                 0.12377E+04, 0.12861E+04, 0.13352E+04, 0.13851E+04, 0.14358E+04,
                                 0.14872E+04, 0.15395E+04, 0.15926E+04, 0.16466E+04, 0.17013E+04,
                                 0.17569E+04, 0.18134E+04, 0.18706E+04, 0.19288E+04, 0.19877E+04,
                                 0.20476E+04, 0.21083E+04, 0.21698E+04, 0.22323E+04, 0.22956E+04,
                                 0.23598E+04, 0.24248E+04, 0.24908E+04, 0.25576E+04, 0.26253E+04,
                                 0.26940E+04, 0.27635E+04, 0.28339E+04, 0.29052E+04, 0.29775E+04,
                                 0.30506E+04, 0.31247E+04, 0.31997E+04, 0.32756E+04, 0.33524E+04,
                                 0.34302E+04, 0.35089E+04, 0.35885E+04, 0.36691E+04, 0.37506E+04,
                                 0.38331E+04, 0.39166E+04, 0.40010E+04, 0.40864E+04, 0.41727E+04,
                                 0.42601E+04, 0.43484E+04, 0.44377E+04, 0.45280E+04, 0.46193E+04,
                                 0.47116E+04, 0.48049E+04, 0.48992E+04, 0.49946E+04, 0.50909E+04,
                                 0.51883E+04, 0.52868E+04, 0.53863E+04, 0.54868E+04, 0.55884E+04,
                                 0.56911E+04, 0.57949E+04, 0.58997E+04, 0.60056E+04, 0.61126E+04,
                                 0.62207E+04, 0.63298E+04, 0.64401E+04, 0.65516E+04, 0.66641E+04,
                                 0.67778E+04, 0.68926E+04, 0.70085E+04, 0.71256E+04, 0.72439E+04,
                                 0.73633E+04, 0.74839E+04, 0.76056E+04, 0.77286E+04, 0.78527E+04,
                                 0.79781E+04, 0.81046E+04, 0.82324E+04, 0.83613E+04, 0.84915E+04,
                                 0.86229E+04, 0.87556E+04, 0.88895E+04, 0.90247E+04, 0.91611E+04,
                                 0.92988E+04])


#  --------------- O2 67: M = 7, I = 3 ---------------------
M = 7
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.52071E+03, 0.74484E+03, 0.96908E+03,
                                 0.11934E+04, 0.14177E+04, 0.16422E+04, 0.18667E+04, 0.20913E+04,
                                 0.23161E+04, 0.25413E+04, 0.27671E+04, 0.29936E+04, 0.32212E+04,
                                 0.34501E+04, 0.36806E+04, 0.39130E+04, 0.41476E+04, 0.43846E+04,
                                 0.46242E+04, 0.48668E+04, 0.51125E+04, 0.53615E+04, 0.56140E+04,
                                 0.58701E+04, 0.61300E+04, 0.63938E+04, 0.66617E+04, 0.69337E+04,
                                 0.72099E+04, 0.74904E+04, 0.77754E+04, 0.80647E+04, 0.83586E+04,
                                 0.86571E+04, 0.89602E+04, 0.92680E+04, 0.95805E+04, 0.98977E+04,
                                 0.10220E+05, 0.10547E+05, 0.10878E+05, 0.11215E+05, 0.11556E+05,
                                 0.11903E+05, 0.12254E+05, 0.12611E+05, 0.12972E+05, 0.13338E+05,
                                 0.13710E+05, 0.14086E+05, 0.14468E+05, 0.14855E+05, 0.15247E+05,
                                 0.15644E+05, 0.16046E+05, 0.16453E+05, 0.16866E+05, 0.17283E+05,
                                 0.17706E+05, 0.18135E+05, 0.18568E+05, 0.19007E+05, 0.19452E+05,
                                 0.19901E+05, 0.20356E+05, 0.20817E+05, 0.21283E+05, 0.21754E+05,
                                 0.22231E+05, 0.22713E+05, 0.23201E+05, 0.23695E+05, 0.24194E+05,
                                 0.24699E+05, 0.25209E+05, 0.25725E+05, 0.26247E+05, 0.26775E+05,
                                 0.27308E+05, 0.27847E+05, 0.28393E+05, 0.28944E+05, 0.29500E+05,
                                 0.30063E+05, 0.30632E+05, 0.31207E+05, 0.31788E+05, 0.32375E+05,
                                 0.32968E+05, 0.33568E+05, 0.34173E+05, 0.34785E+05, 0.35403E+05,
                                 0.36028E+05, 0.36659E+05, 0.37296E+05, 0.37939E+05, 0.38590E+05,
                                 0.39246E+05, 0.39909E+05, 0.40579E+05, 0.41256E+05, 0.41939E+05,
                                 0.42629E+05, 0.43325E+05, 0.44029E+05, 0.44739E+05, 0.45456E+05,
                                 0.46180E+05, 0.46911E+05, 0.47649E+05, 0.48394E+05, 0.49146E+05,
                                 0.49905E+05, 0.50671E+05, 0.51445E+05, 0.52226E+05, 0.53014E+05,
                                 0.53809E+05])


#  --------------- NO 46: M = 8, I = 1 ---------------------
M = 8
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.)
TIPS_ISO_HASH[(M, I)] = float32([0.15840E+03, 0.23971E+03, 0.33080E+03,
                                 0.42907E+03, 0.53251E+03, 0.63972E+03, 0.74975E+03, 0.86195E+03,
                                 0.97582E+03, 0.10911E+04, 0.12074E+04, 0.13248E+04, 0.14430E+04,
                                 0.15621E+04, 0.16820E+04, 0.18027E+04, 0.19243E+04, 0.20468E+04,
                                 0.21703E+04, 0.22948E+04, 0.24204E+04, 0.25472E+04, 0.26753E+04,
                                 0.28046E+04, 0.29354E+04, 0.30676E+04, 0.32013E+04, 0.33365E+04,
                                 0.34734E+04, 0.36120E+04, 0.37522E+04, 0.38942E+04, 0.40379E+04,
                                 0.41835E+04, 0.43310E+04, 0.44803E+04, 0.46316E+04, 0.47849E+04,
                                 0.49400E+04, 0.50972E+04, 0.52564E+04, 0.54176E+04, 0.55809E+04,
                                 0.57462E+04, 0.59137E+04, 0.60832E+04, 0.62548E+04, 0.64286E+04,
                                 0.66045E+04, 0.67825E+04, 0.69628E+04, 0.71451E+04, 0.73297E+04,
                                 0.75164E+04, 0.77053E+04, 0.78964E+04, 0.80897E+04, 0.82853E+04,
                                 0.84830E+04, 0.86830E+04, 0.88852E+04, 0.90896E+04, 0.92963E+04,
                                 0.95052E+04, 0.97164E+04, 0.99297E+04, 0.10145E+05, 0.10363E+05,
                                 0.10583E+05, 0.10806E+05, 0.11031E+05, 0.11258E+05, 0.11487E+05,
                                 0.11718E+05, 0.11952E+05, 0.12188E+05, 0.12426E+05, 0.12667E+05,
                                 0.12910E+05, 0.13155E+05, 0.13403E+05, 0.13652E+05, 0.13905E+05,
                                 0.14159E+05, 0.14416E+05, 0.14675E+05, 0.14936E+05, 0.15199E+05,
                                 0.15465E+05, 0.15733E+05, 0.16004E+05, 0.16277E+05, 0.16552E+05,
                                 0.16829E+05, 0.17109E+05, 0.17391E+05, 0.17675E+05, 0.17962E+05,
                                 0.18251E+05, 0.18542E+05, 0.18836E+05, 0.19131E+05, 0.19430E+05,
                                 0.19730E+05, 0.20033E+05, 0.20338E+05, 0.20646E+05, 0.20955E+05,
                                 0.21268E+05, 0.21582E+05, 0.21899E+05, 0.22218E+05, 0.22539E+05,
                                 0.22863E+05, 0.23189E+05, 0.23518E+05, 0.23848E+05, 0.24181E+05,
                                 0.24517E+05])


#  --------------- NO 56: M = 8, I = 2 ---------------------
M = 8
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.10942E+03, 0.16560E+03, 0.22856E+03,
                                 0.29647E+03, 0.36795E+03, 0.44204E+03, 0.51808E+03, 0.59561E+03,
                                 0.67432E+03, 0.75396E+03, 0.83439E+03, 0.91551E+03, 0.99725E+03,
                                 0.10796E+04, 0.11625E+04, 0.12460E+04, 0.13302E+04, 0.14150E+04,
                                 0.15005E+04, 0.15868E+04, 0.16739E+04, 0.17618E+04, 0.18506E+04,
                                 0.19404E+04, 0.20311E+04, 0.21229E+04, 0.22158E+04, 0.23098E+04,
                                 0.24050E+04, 0.25013E+04, 0.25989E+04, 0.26976E+04, 0.27977E+04,
                                 0.28991E+04, 0.30018E+04, 0.31058E+04, 0.32112E+04, 0.33180E+04,
                                 0.34262E+04, 0.35358E+04, 0.36468E+04, 0.37593E+04, 0.38732E+04,
                                 0.39885E+04, 0.41054E+04, 0.42237E+04, 0.43436E+04, 0.44649E+04,
                                 0.45877E+04, 0.47121E+04, 0.48379E+04, 0.49654E+04, 0.50943E+04,
                                 0.52248E+04, 0.53568E+04, 0.54904E+04, 0.56255E+04, 0.57622E+04,
                                 0.59004E+04, 0.60403E+04, 0.61816E+04, 0.63246E+04, 0.64692E+04,
                                 0.66152E+04, 0.67630E+04, 0.69123E+04, 0.70631E+04, 0.72156E+04,
                                 0.73696E+04, 0.75253E+04, 0.76825E+04, 0.78414E+04, 0.80018E+04,
                                 0.81638E+04, 0.83275E+04, 0.84927E+04, 0.86596E+04, 0.88280E+04,
                                 0.89981E+04, 0.91698E+04, 0.93430E+04, 0.95180E+04, 0.96945E+04,
                                 0.98726E+04, 0.10052E+05, 0.10234E+05, 0.10417E+05, 0.10601E+05,
                                 0.10788E+05, 0.10975E+05, 0.11165E+05, 0.11356E+05, 0.11549E+05,
                                 0.11743E+05, 0.11939E+05, 0.12137E+05, 0.12336E+05, 0.12537E+05,
                                 0.12739E+05, 0.12943E+05, 0.13149E+05, 0.13356E+05, 0.13565E+05,
                                 0.13776E+05, 0.13988E+05, 0.14202E+05, 0.14418E+05, 0.14635E+05,
                                 0.14853E+05, 0.15074E+05, 0.15296E+05, 0.15520E+05, 0.15745E+05,
                                 0.15972E+05, 0.16200E+05, 0.16431E+05, 0.16663E+05, 0.16896E+05,
                                 0.17131E+05])


#  --------------- NO 48: M = 8, I = 3 ---------------------
M = 8
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.)
TIPS_ISO_HASH[(M, I)] = float32([0.16695E+03, 0.25269E+03, 0.34876E+03,
                                 0.45239E+03, 0.56148E+03, 0.67455E+03, 0.79059E+03, 0.90891E+03,
                                 0.10290E+04, 0.11506E+04, 0.12733E+04, 0.13971E+04, 0.15219E+04,
                                 0.16476E+04, 0.17742E+04, 0.19017E+04, 0.20302E+04, 0.21598E+04,
                                 0.22904E+04, 0.24223E+04, 0.25553E+04, 0.26897E+04, 0.28255E+04,
                                 0.29628E+04, 0.31016E+04, 0.32420E+04, 0.33842E+04, 0.35280E+04,
                                 0.36736E+04, 0.38211E+04, 0.39704E+04, 0.41217E+04, 0.42750E+04,
                                 0.44302E+04, 0.45876E+04, 0.47469E+04, 0.49084E+04, 0.50720E+04,
                                 0.52378E+04, 0.54058E+04, 0.55759E+04, 0.57483E+04, 0.59230E+04,
                                 0.60999E+04, 0.62791E+04, 0.64605E+04, 0.66443E+04, 0.68304E+04,
                                 0.70187E+04, 0.72095E+04, 0.74026E+04, 0.75980E+04, 0.77958E+04,
                                 0.79960E+04, 0.81986E+04, 0.84036E+04, 0.86109E+04, 0.88207E+04,
                                 0.90328E+04, 0.92474E+04, 0.94644E+04, 0.96839E+04, 0.99057E+04,
                                 0.10130E+05, 0.10357E+05, 0.10586E+05, 0.10817E+05, 0.11052E+05,
                                 0.11288E+05, 0.11527E+05, 0.11768E+05, 0.12012E+05, 0.12259E+05,
                                 0.12507E+05, 0.12759E+05, 0.13012E+05, 0.13269E+05, 0.13527E+05,
                                 0.13788E+05, 0.14052E+05, 0.14318E+05, 0.14587E+05, 0.14858E+05,
                                 0.15131E+05, 0.15408E+05, 0.15686E+05, 0.15967E+05, 0.16251E+05,
                                 0.16537E+05, 0.16825E+05, 0.17116E+05, 0.17410E+05, 0.17706E+05,
                                 0.18004E+05, 0.18305E+05, 0.18609E+05, 0.18915E+05, 0.19224E+05,
                                 0.19535E+05, 0.19848E+05, 0.20164E+05, 0.20483E+05, 0.20804E+05,
                                 0.21127E+05, 0.21453E+05, 0.21782E+05, 0.22113E+05, 0.22447E+05,
                                 0.22783E+05, 0.23122E+05, 0.23463E+05, 0.23807E+05, 0.24153E+05,
                                 0.24502E+05, 0.24853E+05, 0.25207E+05, 0.25563E+05, 0.25922E+05,
                                 0.26283E+05])


#  --------------- SO2 626: M = 9, I = 1 ---------------------
M = 9
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.52899E+03, 0.89171E+03, 0.13139E+04,
                                 0.17915E+04, 0.23246E+04, 0.29155E+04, 0.35675E+04, 0.42848E+04,
                                 0.50723E+04, 0.59352E+04, 0.68794E+04, 0.79109E+04, 0.90366E+04,
                                 0.10264E+05, 0.11599E+05, 0.13052E+05, 0.14629E+05, 0.16340E+05,
                                 0.18193E+05, 0.20199E+05, 0.22366E+05, 0.24704E+05, 0.27225E+05,
                                 0.29938E+05, 0.32855E+05, 0.35987E+05, 0.39346E+05, 0.42944E+05,
                                 0.46794E+05, 0.50909E+05, 0.55302E+05, 0.59986E+05, 0.64977E+05,
                                 0.70288E+05, 0.75934E+05, 0.81931E+05, 0.88294E+05, 0.95040E+05,
                                 0.10219E+06, 0.10975E+06, 0.11774E+06, 0.12619E+06, 0.13511E+06,
                                 0.14452E+06, 0.15443E+06, 0.16487E+06, 0.17586E+06, 0.18742E+06,
                                 0.19957E+06, 0.21234E+06, 0.22573E+06, 0.23978E+06, 0.25451E+06,
                                 0.26995E+06, 0.28611E+06, 0.30302E+06, 0.32071E+06, 0.33920E+06,
                                 0.35852E+06, 0.37869E+06, 0.39974E+06, 0.42171E+06, 0.44461E+06,
                                 0.46848E+06, 0.49334E+06, 0.51922E+06, 0.54617E+06, 0.57419E+06,
                                 0.60334E+06, 0.63363E+06, 0.66511E+06, 0.69780E+06, 0.73174E+06,
                                 0.76696E+06, 0.80349E+06, 0.84138E+06, 0.88066E+06, 0.92136E+06,
                                 0.96352E+06, 0.10072E+07, 0.10524E+07, 0.10992E+07, 0.11475E+07,
                                 0.11976E+07, 0.12493E+07, 0.13028E+07, 0.13580E+07, 0.14151E+07,
                                 0.14741E+07, 0.15349E+07, 0.15977E+07, 0.16625E+07, 0.17293E+07,
                                 0.17982E+07, 0.18693E+07, 0.19425E+07, 0.20180E+07, 0.20958E+07,
                                 0.21758E+07, 0.22583E+07, 0.23432E+07, 0.24305E+07, 0.25204E+07,
                                 0.26129E+07, 0.27080E+07, 0.28058E+07, 0.29064E+07, 0.30097E+07,
                                 0.31159E+07, 0.32250E+07, 0.33371E+07, 0.34522E+07, 0.35705E+07,
                                 0.36918E+07, 0.38164E+07, 0.39442E+07, 0.40754E+07, 0.42099E+07,
                                 0.43479E+07])


#  --------------- SO2 646: M = 9, I = 2 ---------------------
M = 9
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.53140E+03, 0.89578E+03, 0.13199E+04,
                                 0.17997E+04, 0.23353E+04, 0.29288E+04, 0.35837E+04, 0.43043E+04,
                                 0.50953E+04, 0.59621E+04, 0.69104E+04, 0.79465E+04, 0.90772E+04,
                                 0.10310E+05, 0.11651E+05, 0.13110E+05, 0.14694E+05, 0.16413E+05,
                                 0.18274E+05, 0.20289E+05, 0.22465E+05, 0.24814E+05, 0.27345E+05,
                                 0.30070E+05, 0.33000E+05, 0.36145E+05, 0.39519E+05, 0.43133E+05,
                                 0.46999E+05, 0.51132E+05, 0.55544E+05, 0.60248E+05, 0.65260E+05,
                                 0.70594E+05, 0.76264E+05, 0.82287E+05, 0.88678E+05, 0.95453E+05,
                                 0.10263E+06, 0.11022E+06, 0.11825E+06, 0.12674E+06, 0.13569E+06,
                                 0.14514E+06, 0.15510E+06, 0.16558E+06, 0.17662E+06, 0.18823E+06,
                                 0.20043E+06, 0.21325E+06, 0.22670E+06, 0.24081E+06, 0.25561E+06,
                                 0.27111E+06, 0.28733E+06, 0.30432E+06, 0.32208E+06, 0.34065E+06,
                                 0.36005E+06, 0.38031E+06, 0.40145E+06, 0.42351E+06, 0.44651E+06,
                                 0.47047E+06, 0.49544E+06, 0.52144E+06, 0.54849E+06, 0.57664E+06,
                                 0.60591E+06, 0.63633E+06, 0.66794E+06, 0.70077E+06, 0.73485E+06,
                                 0.77022E+06, 0.80691E+06, 0.84496E+06, 0.88440E+06, 0.92527E+06,
                                 0.96761E+06, 0.10115E+07, 0.10568E+07, 0.11038E+07, 0.11524E+07,
                                 0.12027E+07, 0.12546E+07, 0.13083E+07, 0.13638E+07, 0.14211E+07,
                                 0.14803E+07, 0.15414E+07, 0.16045E+07, 0.16695E+07, 0.17366E+07,
                                 0.18059E+07, 0.18772E+07, 0.19507E+07, 0.20265E+07, 0.21046E+07,
                                 0.21850E+07, 0.22678E+07, 0.23531E+07, 0.24408E+07, 0.25310E+07,
                                 0.26239E+07, 0.27194E+07, 0.28176E+07, 0.29186E+07, 0.30224E+07,
                                 0.31290E+07, 0.32386E+07, 0.33512E+07, 0.34668E+07, 0.35855E+07,
                                 0.37074E+07, 0.38324E+07, 0.39608E+07, 0.40925E+07, 0.42276E+07,
                                 0.43662E+07])

#  --------------- NO2 646: M = 10, I = 1 ---------------------
M = 10
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.)
TIPS_ISO_HASH[(M, I)] = float32([0.12046E+04, 0.20297E+04, 0.29875E+04,
                                 0.40626E+04, 0.52463E+04, 0.65350E+04, 0.79286E+04, 0.94298E+04,
                                 0.11043E+05, 0.12776E+05, 0.14634E+05, 0.16627E+05, 0.18765E+05,
                                 0.21056E+05, 0.23511E+05, 0.26143E+05, 0.28961E+05, 0.31979E+05,
                                 0.35209E+05, 0.38663E+05, 0.42355E+05, 0.46300E+05, 0.50510E+05,
                                 0.55001E+05, 0.59787E+05, 0.64884E+05, 0.70308E+05, 0.76075E+05,
                                 0.82201E+05, 0.88704E+05, 0.95602E+05, 0.10291E+06, 0.11065E+06,
                                 0.11884E+06, 0.12750E+06, 0.13665E+06, 0.14631E+06, 0.15650E+06,
                                 0.16724E+06, 0.17856E+06, 0.19047E+06, 0.20301E+06, 0.21618E+06,
                                 0.23002E+06, 0.24456E+06, 0.25981E+06, 0.27580E+06, 0.29256E+06,
                                 0.31012E+06, 0.32850E+06, 0.34773E+06, 0.36784E+06, 0.38886E+06,
                                 0.41082E+06, 0.43374E+06, 0.45766E+06, 0.48262E+06, 0.50863E+06,
                                 0.53574E+06, 0.56398E+06, 0.59339E+06, 0.62398E+06, 0.65581E+06,
                                 0.68891E+06, 0.72331E+06, 0.75905E+06, 0.79617E+06, 0.83470E+06,
                                 0.87469E+06, 0.91617E+06, 0.95919E+06, 0.10038E+07, 0.10500E+07,
                                 0.10979E+07, 0.11474E+07, 0.11988E+07, 0.12519E+07, 0.13068E+07,
                                 0.13636E+07, 0.14224E+07, 0.14831E+07, 0.15459E+07, 0.16107E+07,
                                 0.16776E+07, 0.17467E+07, 0.18180E+07, 0.18916E+07, 0.19675E+07,
                                 0.20458E+07, 0.21265E+07, 0.22097E+07, 0.22954E+07, 0.23837E+07,
                                 0.24747E+07, 0.25684E+07, 0.26648E+07, 0.27641E+07, 0.28662E+07,
                                 0.29713E+07, 0.30794E+07, 0.31905E+07, 0.33048E+07, 0.34223E+07,
                                 0.35430E+07, 0.36670E+07, 0.37944E+07, 0.39253E+07, 0.40597E+07,
                                 0.41976E+07, 0.43393E+07, 0.44846E+07, 0.46337E+07, 0.47867E+07,
                                 0.49437E+07, 0.51046E+07, 0.52696E+07, 0.54388E+07, 0.56122E+07,
                                 0.57900E+07])


#  --------------- NH3 4111: M = 11, I = 1 ---------------------
M = 11
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.)
TIPS_ISO_HASH[(M, I)] = float32([0.16013E+03, 0.26692E+03, 0.39067E+03,
                                 0.52933E+03, 0.68153E+03, 0.84641E+03, 0.10234E+04, 0.12125E+04,
                                 0.14136E+04, 0.16272E+04, 0.18537E+04, 0.20937E+04, 0.23481E+04,
                                 0.26177E+04, 0.29035E+04, 0.32065E+04, 0.35279E+04, 0.38688E+04,
                                 0.42304E+04, 0.46141E+04, 0.50212E+04, 0.54531E+04, 0.59114E+04,
                                 0.63976E+04, 0.69133E+04, 0.74602E+04, 0.80401E+04, 0.86549E+04,
                                 0.93066E+04, 0.99971E+04, 0.10729E+05, 0.11504E+05, 0.12324E+05,
                                 0.13193E+05, 0.14112E+05, 0.15085E+05, 0.16114E+05, 0.17201E+05,
                                 0.18352E+05, 0.19567E+05, 0.20851E+05, 0.22208E+05, 0.23640E+05,
                                 0.25152E+05, 0.26747E+05, 0.28430E+05, 0.30205E+05, 0.32077E+05,
                                 0.34050E+05, 0.36128E+05, 0.38317E+05, 0.40623E+05, 0.43050E+05,
                                 0.45605E+05, 0.48292E+05, 0.51119E+05, 0.54091E+05, 0.57215E+05,
                                 0.60498E+05, 0.63947E+05, 0.67569E+05, 0.71372E+05, 0.75364E+05,
                                 0.79552E+05, 0.83946E+05, 0.88553E+05, 0.93384E+05, 0.98447E+05,
                                 0.10375E+06, 0.10931E+06, 0.11513E+06, 0.12122E+06, 0.12760E+06,
                                 0.13427E+06, 0.14125E+06, 0.14855E+06, 0.15619E+06, 0.16417E+06,
                                 0.17250E+06, 0.18121E+06, 0.19031E+06, 0.19981E+06, 0.20973E+06,
                                 0.22008E+06, 0.23088E+06, 0.24215E+06, 0.25390E+06, 0.26615E+06,
                                 0.27892E+06, 0.29223E+06, 0.30610E+06, 0.32055E+06, 0.33559E+06,
                                 0.35125E+06, 0.36756E+06, 0.38453E+06, 0.40219E+06, 0.42056E+06,
                                 0.43967E+06, 0.45953E+06, 0.48019E+06, 0.50165E+06, 0.52396E+06,
                                 0.54714E+06, 0.57122E+06, 0.59622E+06, 0.62218E+06, 0.64913E+06,
                                 0.67710E+06, 0.70613E+06, 0.73624E+06, 0.76748E+06, 0.79988E+06,
                                 0.83347E+06, 0.86829E+06, 0.90439E+06, 0.94180E+06, 0.98056E+06,
                                 0.10207E+07])


#  --------------- NH3 5111: M = 11, I = 2 ---------------------
M = 11
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.10697E+03, 0.17832E+03, 0.26100E+03,
                                 0.35364E+03, 0.45533E+03, 0.56549E+03, 0.68377E+03, 0.81007E+03,
                                 0.94447E+03, 0.10872E+04, 0.12385E+04, 0.13988E+04, 0.15688E+04,
                                 0.17490E+04, 0.19399E+04, 0.21424E+04, 0.23571E+04, 0.25848E+04,
                                 0.28264E+04, 0.30828E+04, 0.33548E+04, 0.36434E+04, 0.39496E+04,
                                 0.42745E+04, 0.46190E+04, 0.49845E+04, 0.53720E+04, 0.57828E+04,
                                 0.62182E+04, 0.66796E+04, 0.71684E+04, 0.76862E+04, 0.82344E+04,
                                 0.88149E+04, 0.94292E+04, 0.10079E+05, 0.10767E+05, 0.11494E+05,
                                 0.12262E+05, 0.13074E+05, 0.13932E+05, 0.14839E+05, 0.15796E+05,
                                 0.16806E+05, 0.17872E+05, 0.18997E+05, 0.20183E+05, 0.21434E+05,
                                 0.22752E+05, 0.24141E+05, 0.25604E+05, 0.27145E+05, 0.28767E+05,
                                 0.30475E+05, 0.32271E+05, 0.34160E+05, 0.36146E+05, 0.38234E+05,
                                 0.40428E+05, 0.42733E+05, 0.45154E+05, 0.47696E+05, 0.50364E+05,
                                 0.53163E+05, 0.56100E+05, 0.59180E+05, 0.62408E+05, 0.65792E+05,
                                 0.69339E+05, 0.73053E+05, 0.76943E+05, 0.81016E+05, 0.85279E+05,
                                 0.89740E+05, 0.94406E+05, 0.99287E+05, 0.10439E+06, 0.10972E+06,
                                 0.11530E+06, 0.12112E+06, 0.12720E+06, 0.13355E+06, 0.14018E+06,
                                 0.14711E+06, 0.15433E+06, 0.16186E+06, 0.16971E+06, 0.17791E+06,
                                 0.18645E+06, 0.19534E+06, 0.20462E+06, 0.21428E+06, 0.22434E+06,
                                 0.23481E+06, 0.24572E+06, 0.25706E+06, 0.26887E+06, 0.28116E+06,
                                 0.29393E+06, 0.30722E+06, 0.32103E+06, 0.33539E+06, 0.35031E+06,
                                 0.36581E+06, 0.38191E+06, 0.39864E+06, 0.41600E+06, 0.43403E+06,
                                 0.45274E+06, 0.47215E+06, 0.49230E+06, 0.51319E+06, 0.53487E+06,
                                 0.55734E+06, 0.58064E+06, 0.60478E+06, 0.62981E+06, 0.65574E+06,
                                 0.68260E+06])


#  --------------- HNO3 146: M = 12, I = 1 ---------------------
M = 12
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.15010E+05, 0.25316E+05, 0.37374E+05,
                                 0.51216E+05, 0.67105E+05, 0.85473E+05, 0.10688E+06, 0.13201E+06,
                                 0.16165E+06, 0.19671E+06, 0.23825E+06, 0.28749E+06, 0.34583E+06,
                                 0.41490E+06, 0.49657E+06, 0.59302E+06, 0.70673E+06, 0.84054E+06,
                                 0.99775E+06, 0.11821E+07, 0.13978E+07, 0.16498E+07, 0.19436E+07,
                                 0.22855E+07, 0.26825E+07, 0.31428E+07, 0.36753E+07, 0.42903E+07,
                                 0.49993E+07, 0.58151E+07, 0.67523E+07, 0.78269E+07, 0.90572E+07,
                                 0.10463E+08, 0.12067E+08, 0.13895E+08, 0.15973E+08, 0.18333E+08,
                                 0.21009E+08, 0.24039E+08, 0.27464E+08, 0.31331E+08, 0.35690E+08,
                                 0.40597E+08, 0.46115E+08, 0.52310E+08, 0.59257E+08, 0.67037E+08,
                                 0.75739E+08, 0.85461E+08, 0.96310E+08, 0.10840E+09, 0.12186E+09,
                                 0.13683E+09, 0.15346E+09, 0.17191E+09, 0.19236E+09, 0.21501E+09,
                                 0.24006E+09, 0.26774E+09, 0.29830E+09, 0.33200E+09, 0.36914E+09,
                                 0.41002E+09, 0.45498E+09, 0.50438E+09, 0.55862E+09, 0.61812E+09,
                                 0.68332E+09, 0.75473E+09, 0.83286E+09, 0.91828E+09, 0.10116E+10,
                                 0.11134E+10, 0.12245E+10, 0.13456E+10, 0.14775E+10, 0.16210E+10,
                                 0.17771E+10, 0.19467E+10, 0.21309E+10, 0.23309E+10, 0.25477E+10,
                                 0.27827E+10, 0.30372E+10, 0.33127E+10, 0.36107E+10, 0.39329E+10,
                                 0.42809E+10, 0.46567E+10, 0.50623E+10, 0.54997E+10, 0.59711E+10,
                                 0.64789E+10, 0.70257E+10, 0.76140E+10, 0.82468E+10, 0.89269E+10,
                                 0.96575E+10, 0.10442E+11, 0.11284E+11, 0.12187E+11, 0.13155E+11,
                                 0.14193E+11, 0.15304E+11, 0.16494E+11, 0.17767E+11, 0.19129E+11,
                                 0.20585E+11, 0.22140E+11, 0.23802E+11, 0.25576E+11, 0.27469E+11,
                                 0.29489E+11, 0.31642E+11, 0.33937E+11, 0.36382E+11, 0.38985E+11,
                                 0.41757E+11])


#  --------------- HNO3 156: M = 12, I = 2 --------------------- NOT IN TIPS-2011
M = 12
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- OH 61: M = 13, I = 1 ---------------------
M = 13
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.20066E+02, 0.24774E+02, 0.30309E+02,
                                 0.36357E+02, 0.42745E+02, 0.49371E+02, 0.56168E+02, 0.63093E+02,
                                 0.70116E+02, 0.77217E+02, 0.84380E+02, 0.91594E+02, 0.98850E+02,
                                 0.10614E+03, 0.11346E+03, 0.12081E+03, 0.12818E+03, 0.13557E+03,
                                 0.14298E+03, 0.15041E+03, 0.15785E+03, 0.16531E+03, 0.17278E+03,
                                 0.18027E+03, 0.18778E+03, 0.19530E+03, 0.20284E+03, 0.21040E+03,
                                 0.21797E+03, 0.22556E+03, 0.23318E+03, 0.24082E+03, 0.24848E+03,
                                 0.25617E+03, 0.26389E+03, 0.27163E+03, 0.27941E+03, 0.28721E+03,
                                 0.29505E+03, 0.30292E+03, 0.31084E+03, 0.31878E+03, 0.32677E+03,
                                 0.33480E+03, 0.34287E+03, 0.35099E+03, 0.35915E+03, 0.36736E+03,
                                 0.37561E+03, 0.38391E+03, 0.39227E+03, 0.40067E+03, 0.40913E+03,
                                 0.41764E+03, 0.42620E+03, 0.43482E+03, 0.44350E+03, 0.45223E+03,
                                 0.46102E+03, 0.46987E+03, 0.47878E+03, 0.48775E+03, 0.49679E+03,
                                 0.50588E+03, 0.51503E+03, 0.52425E+03, 0.53354E+03, 0.54288E+03,
                                 0.55229E+03, 0.56177E+03, 0.57132E+03, 0.58092E+03, 0.59060E+03,
                                 0.60035E+03, 0.61016E+03, 0.62004E+03, 0.62999E+03, 0.64001E+03,
                                 0.65010E+03, 0.66025E+03, 0.67049E+03, 0.68078E+03, 0.69115E+03,
                                 0.70160E+03, 0.71211E+03, 0.72269E+03, 0.73335E+03, 0.74408E+03,
                                 0.75488E+03, 0.76576E+03, 0.77671E+03, 0.78773E+03, 0.79883E+03,
                                 0.81000E+03, 0.82124E+03, 0.83256E+03, 0.84396E+03, 0.85542E+03,
                                 0.86696E+03, 0.87858E+03, 0.89027E+03, 0.90204E+03, 0.91389E+03,
                                 0.92580E+03, 0.93781E+03, 0.94988E+03, 0.96203E+03, 0.97425E+03,
                                 0.98656E+03, 0.99893E+03, 0.10114E+04, 0.10239E+04, 0.10365E+04,
                                 0.10492E+04, 0.10620E+04, 0.10748E+04, 0.10878E+04, 0.11007E+04,
                                 0.11138E+04])

#  --------------- OH 81: M = 13, I = 2 ---------------------
M = 13
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.20124E+02, 0.24876E+02, 0.30457E+02,
                                 0.36553E+02, 0.42991E+02, 0.49666E+02, 0.56513E+02, 0.63489E+02,
                                 0.70563E+02, 0.77715E+02, 0.84929E+02, 0.92195E+02, 0.99504E+02,
                                 0.10685E+03, 0.11423E+03, 0.12164E+03, 0.12907E+03, 0.13654E+03,
                                 0.14403E+03, 0.15154E+03, 0.15909E+03, 0.16666E+03, 0.17427E+03,
                                 0.18191E+03, 0.18959E+03, 0.19731E+03, 0.20507E+03, 0.21287E+03,
                                 0.22073E+03, 0.22863E+03, 0.23658E+03, 0.24459E+03, 0.25266E+03,
                                 0.26078E+03, 0.26897E+03, 0.27722E+03, 0.28554E+03, 0.29393E+03,
                                 0.30238E+03, 0.31091E+03, 0.31952E+03, 0.32820E+03, 0.33696E+03,
                                 0.34579E+03, 0.35471E+03, 0.36371E+03, 0.37279E+03, 0.38196E+03,
                                 0.39121E+03, 0.40055E+03, 0.40998E+03, 0.41949E+03, 0.42910E+03,
                                 0.43879E+03, 0.44858E+03, 0.45845E+03, 0.46843E+03, 0.47849E+03,
                                 0.48865E+03, 0.49890E+03, 0.50924E+03, 0.51969E+03, 0.53022E+03,
                                 0.54086E+03, 0.55159E+03, 0.56242E+03, 0.57335E+03, 0.58437E+03,
                                 0.59550E+03, 0.60673E+03, 0.61805E+03, 0.62947E+03, 0.64100E+03,
                                 0.65263E+03, 0.66435E+03, 0.67618E+03, 0.68811E+03, 0.70014E+03,
                                 0.71228E+03, 0.72451E+03, 0.73685E+03, 0.74929E+03, 0.76184E+03,
                                 0.77449E+03, 0.78724E+03, 0.80009E+03, 0.81306E+03, 0.82612E+03,
                                 0.83929E+03, 0.85256E+03, 0.86594E+03, 0.87942E+03, 0.89301E+03,
                                 0.90670E+03, 0.92050E+03, 0.93440E+03, 0.94841E+03, 0.96253E+03,
                                 0.97675E+03, 0.99108E+03, 0.10055E+04, 0.10201E+04, 0.10347E+04,
                                 0.10495E+04, 0.10643E+04, 0.10793E+04, 0.10944E+04, 0.11096E+04,
                                 0.11248E+04, 0.11402E+04, 0.11558E+04, 0.11714E+04, 0.11871E+04,
                                 0.12029E+04, 0.12189E+04, 0.12349E+04, 0.12511E+04, 0.12673E+04,
                                 0.12837E+04])

#  --------------- OH 62: M = 13, I = 3 ---------------------
M = 13
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.)
TIPS_ISO_HASH[(M, I)] = float32([0.41032E+02, 0.54704E+02, 0.70201E+02,
                                 0.86985E+02, 0.10469E+03, 0.12306E+03, 0.14194E+03, 0.16119E+03,
                                 0.18075E+03, 0.20054E+03, 0.22053E+03, 0.24068E+03, 0.26096E+03,
                                 0.28135E+03, 0.30183E+03, 0.32241E+03, 0.34305E+03, 0.36376E+03,
                                 0.38453E+03, 0.40535E+03, 0.42622E+03, 0.44714E+03, 0.46811E+03,
                                 0.48913E+03, 0.51019E+03, 0.53131E+03, 0.55246E+03, 0.57368E+03,
                                 0.59495E+03, 0.61627E+03, 0.63766E+03, 0.65912E+03, 0.68064E+03,
                                 0.70223E+03, 0.72390E+03, 0.74565E+03, 0.76749E+03, 0.78941E+03,
                                 0.81143E+03, 0.83355E+03, 0.85578E+03, 0.87810E+03, 0.90054E+03,
                                 0.92310E+03, 0.94577E+03, 0.96857E+03, 0.99149E+03, 0.10145E+04,
                                 0.10377E+04, 0.10611E+04, 0.10845E+04, 0.11081E+04, 0.11319E+04,
                                 0.11558E+04, 0.11798E+04, 0.12040E+04, 0.12284E+04, 0.12529E+04,
                                 0.12776E+04, 0.13025E+04, 0.13275E+04, 0.13527E+04, 0.13781E+04,
                                 0.14036E+04, 0.14293E+04, 0.14552E+04, 0.14813E+04, 0.15076E+04,
                                 0.15340E+04, 0.15606E+04, 0.15874E+04, 0.16144E+04, 0.16416E+04,
                                 0.16690E+04, 0.16965E+04, 0.17243E+04, 0.17522E+04, 0.17804E+04,
                                 0.18087E+04, 0.18373E+04, 0.18660E+04, 0.18949E+04, 0.19241E+04,
                                 0.19534E+04, 0.19829E+04, 0.20127E+04, 0.20426E+04, 0.20727E+04,
                                 0.21031E+04, 0.21336E+04, 0.21644E+04, 0.21954E+04, 0.22266E+04,
                                 0.22579E+04, 0.22895E+04, 0.23213E+04, 0.23534E+04, 0.23856E+04,
                                 0.24180E+04, 0.24506E+04, 0.24835E+04, 0.25166E+04, 0.25499E+04,
                                 0.25834E+04, 0.26171E+04, 0.26510E+04, 0.26852E+04, 0.27195E+04,
                                 0.27541E+04, 0.27889E+04, 0.28239E+04, 0.28592E+04, 0.28946E+04,
                                 0.29303E+04, 0.29661E+04, 0.30023E+04, 0.30386E+04, 0.30751E+04,
                                 0.31119E+04])


#  --------------- HF 19: M = 14, I = 1 ---------------------
M = 14
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.95958E+01, 0.12933E+02, 0.16295E+02,
                                 0.19666E+02, 0.23043E+02, 0.26425E+02, 0.29809E+02, 0.33195E+02,
                                 0.36584E+02, 0.39974E+02, 0.43366E+02, 0.46759E+02, 0.50154E+02,
                                 0.53550E+02, 0.56947E+02, 0.60346E+02, 0.63746E+02, 0.67148E+02,
                                 0.70550E+02, 0.73955E+02, 0.77361E+02, 0.80769E+02, 0.84179E+02,
                                 0.87591E+02, 0.91006E+02, 0.94424E+02, 0.97846E+02, 0.10127E+03,
                                 0.10470E+03, 0.10813E+03, 0.11157E+03, 0.11502E+03, 0.11847E+03,
                                 0.12193E+03, 0.12540E+03, 0.12888E+03, 0.13236E+03, 0.13586E+03,
                                 0.13936E+03, 0.14288E+03, 0.14641E+03, 0.14995E+03, 0.15351E+03,
                                 0.15708E+03, 0.16066E+03, 0.16426E+03, 0.16788E+03, 0.17151E+03,
                                 0.17516E+03, 0.17882E+03, 0.18251E+03, 0.18621E+03, 0.18994E+03,
                                 0.19368E+03, 0.19745E+03, 0.20123E+03, 0.20504E+03, 0.20887E+03,
                                 0.21272E+03, 0.21659E+03, 0.22049E+03, 0.22441E+03, 0.22836E+03,
                                 0.23233E+03, 0.23632E+03, 0.24034E+03, 0.24439E+03, 0.24846E+03,
                                 0.25255E+03, 0.25668E+03, 0.26083E+03, 0.26501E+03, 0.26921E+03,
                                 0.27344E+03, 0.27770E+03, 0.28199E+03, 0.28631E+03, 0.29066E+03,
                                 0.29503E+03, 0.29944E+03, 0.30387E+03, 0.30833E+03, 0.31282E+03,
                                 0.31735E+03, 0.32190E+03, 0.32648E+03, 0.33110E+03, 0.33574E+03,
                                 0.34042E+03, 0.34512E+03, 0.34986E+03, 0.35463E+03, 0.35943E+03,
                                 0.36426E+03, 0.36913E+03, 0.37402E+03, 0.37895E+03, 0.38391E+03,
                                 0.38891E+03, 0.39393E+03, 0.39899E+03, 0.40408E+03, 0.40921E+03,
                                 0.41436E+03, 0.41955E+03, 0.42478E+03, 0.43004E+03, 0.43533E+03,
                                 0.44065E+03, 0.44601E+03, 0.45140E+03, 0.45683E+03, 0.46229E+03,
                                 0.46779E+03, 0.47332E+03, 0.47888E+03, 0.48448E+03, 0.49011E+03,
                                 0.49578E+03])

#  --------------- HF 29: M = 14, I = 2 --------------------- not in TIPS-2011
M = 14
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- HСl 15: M = 15, I = 1 --------------------
M = 15
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.)
TIPS_ISO_HASH[(M, I)] = float32([0.34775E+02, 0.48060E+02, 0.61370E+02,
                                 0.74692E+02, 0.88024E+02, 0.10136E+03, 0.11471E+03, 0.12806E+03,
                                 0.14141E+03, 0.15478E+03, 0.16814E+03, 0.18151E+03, 0.19489E+03,
                                 0.20827E+03, 0.22166E+03, 0.23506E+03, 0.24847E+03, 0.26189E+03,
                                 0.27533E+03, 0.28878E+03, 0.30225E+03, 0.31575E+03, 0.32928E+03,
                                 0.34284E+03, 0.35645E+03, 0.37009E+03, 0.38378E+03, 0.39753E+03,
                                 0.41134E+03, 0.42521E+03, 0.43914E+03, 0.45316E+03, 0.46725E+03,
                                 0.48142E+03, 0.49568E+03, 0.51003E+03, 0.52448E+03, 0.53902E+03,
                                 0.55368E+03, 0.56843E+03, 0.58330E+03, 0.59829E+03, 0.61339E+03,
                                 0.62862E+03, 0.64396E+03, 0.65944E+03, 0.67504E+03, 0.69078E+03,
                                 0.70665E+03, 0.72265E+03, 0.73880E+03, 0.75508E+03, 0.77151E+03,
                                 0.78809E+03, 0.80481E+03, 0.82168E+03, 0.83870E+03, 0.85587E+03,
                                 0.87320E+03, 0.89068E+03, 0.90832E+03, 0.92611E+03, 0.94407E+03,
                                 0.96218E+03, 0.98046E+03, 0.99889E+03, 0.10175E+04, 0.10363E+04,
                                 0.10552E+04, 0.10743E+04, 0.10936E+04, 0.11130E+04, 0.11326E+04,
                                 0.11524E+04, 0.11723E+04, 0.11924E+04, 0.12127E+04, 0.12332E+04,
                                 0.12538E+04, 0.12746E+04, 0.12956E+04, 0.13168E+04, 0.13381E+04,
                                 0.13597E+04, 0.13814E+04, 0.14032E+04, 0.14253E+04, 0.14475E+04,
                                 0.14700E+04, 0.14926E+04, 0.15153E+04, 0.15383E+04, 0.15615E+04,
                                 0.15848E+04, 0.16083E+04, 0.16320E+04, 0.16559E+04, 0.16800E+04,
                                 0.17043E+04, 0.17287E+04, 0.17533E+04, 0.17782E+04, 0.18032E+04,
                                 0.18284E+04, 0.18538E+04, 0.18794E+04, 0.19051E+04, 0.19311E+04,
                                 0.19573E+04, 0.19836E+04, 0.20102E+04, 0.20369E+04, 0.20638E+04,
                                 0.20910E+04, 0.21183E+04, 0.21458E+04, 0.21735E+04, 0.22014E+04,
                                 0.22295E+04])


#  --------------- HСl 17: M = 15, I = 2 ---------------------
M = 15
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.)
TIPS_ISO_HASH[(M, I)] = float32([0.34823E+02, 0.48128E+02, 0.61458E+02,
                                 0.74801E+02, 0.88152E+02, 0.10151E+03, 0.11488E+03, 0.12825E+03,
                                 0.14162E+03, 0.15500E+03, 0.16839E+03, 0.18178E+03, 0.19518E+03,
                                 0.20858E+03, 0.22199E+03, 0.23541E+03, 0.24884E+03, 0.26228E+03,
                                 0.27574E+03, 0.28921E+03, 0.30270E+03, 0.31622E+03, 0.32977E+03,
                                 0.34336E+03, 0.35698E+03, 0.37065E+03, 0.38436E+03, 0.39813E+03,
                                 0.41196E+03, 0.42585E+03, 0.43981E+03, 0.45384E+03, 0.46796E+03,
                                 0.48215E+03, 0.49644E+03, 0.51081E+03, 0.52528E+03, 0.53986E+03,
                                 0.55453E+03, 0.56932E+03, 0.58421E+03, 0.59922E+03, 0.61435E+03,
                                 0.62960E+03, 0.64498E+03, 0.66048E+03, 0.67611E+03, 0.69187E+03,
                                 0.70777E+03, 0.72381E+03, 0.73998E+03, 0.75630E+03, 0.77276E+03,
                                 0.78936E+03, 0.80612E+03, 0.82302E+03, 0.84007E+03, 0.85727E+03,
                                 0.87463E+03, 0.89215E+03, 0.90982E+03, 0.92765E+03, 0.94563E+03,
                                 0.96378E+03, 0.98209E+03, 0.10006E+04, 0.10192E+04, 0.10380E+04,
                                 0.10570E+04, 0.10761E+04, 0.10954E+04, 0.11149E+04, 0.11345E+04,
                                 0.11543E+04, 0.11743E+04, 0.11945E+04, 0.12148E+04, 0.12353E+04,
                                 0.12560E+04, 0.12768E+04, 0.12979E+04, 0.13191E+04, 0.13405E+04,
                                 0.13620E+04, 0.13838E+04, 0.14057E+04, 0.14278E+04, 0.14501E+04,
                                 0.14726E+04, 0.14952E+04, 0.15180E+04, 0.15410E+04, 0.15642E+04,
                                 0.15876E+04, 0.16112E+04, 0.16349E+04, 0.16589E+04, 0.16830E+04,
                                 0.17073E+04, 0.17318E+04, 0.17565E+04, 0.17814E+04, 0.18064E+04,
                                 0.18317E+04, 0.18572E+04, 0.18828E+04, 0.19086E+04, 0.19346E+04,
                                 0.19609E+04, 0.19873E+04, 0.20139E+04, 0.20406E+04, 0.20676E+04,
                                 0.20948E+04, 0.21222E+04, 0.21498E+04, 0.21775E+04, 0.22055E+04,
                                 0.22337E+04])


#  --------------- HСl 25: M = 15, I = 3 --------------------- not in TIPS-2011
M = 15
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- HСl 27: M = 15, I = 4 --------------------- not in TIPS-2011
M = 15
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- HBr 19: M = 16, I = 1 ---------------------
M = 16
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.)
TIPS_ISO_HASH[(M, I)] = float32([0.42744E+02, 0.59373E+02, 0.76023E+02,
                                 0.92685E+02, 0.10936E+03, 0.12604E+03, 0.14272E+03, 0.15942E+03,
                                 0.17612E+03, 0.19282E+03, 0.20954E+03, 0.22626E+03, 0.24299E+03,
                                 0.25973E+03, 0.27648E+03, 0.29325E+03, 0.31004E+03, 0.32686E+03,
                                 0.34371E+03, 0.36060E+03, 0.37753E+03, 0.39451E+03, 0.41156E+03,
                                 0.42868E+03, 0.44587E+03, 0.46314E+03, 0.48051E+03, 0.49798E+03,
                                 0.51556E+03, 0.53325E+03, 0.55106E+03, 0.56900E+03, 0.58708E+03,
                                 0.60530E+03, 0.62367E+03, 0.64219E+03, 0.66088E+03, 0.67972E+03,
                                 0.69874E+03, 0.71793E+03, 0.73730E+03, 0.75685E+03, 0.77659E+03,
                                 0.79652E+03, 0.81664E+03, 0.83696E+03, 0.85748E+03, 0.87820E+03,
                                 0.89914E+03, 0.92028E+03, 0.94163E+03, 0.96319E+03, 0.98498E+03,
                                 0.10070E+04, 0.10292E+04, 0.10516E+04, 0.10743E+04, 0.10972E+04,
                                 0.11203E+04, 0.11437E+04, 0.11673E+04, 0.11911E+04, 0.12151E+04,
                                 0.12394E+04, 0.12640E+04, 0.12887E+04, 0.13137E+04, 0.13390E+04,
                                 0.13645E+04, 0.13902E+04, 0.14162E+04, 0.14424E+04, 0.14689E+04,
                                 0.14956E+04, 0.15226E+04, 0.15498E+04, 0.15773E+04, 0.16050E+04,
                                 0.16330E+04, 0.16612E+04, 0.16897E+04, 0.17185E+04, 0.17475E+04,
                                 0.17767E+04, 0.18062E+04, 0.18360E+04, 0.18660E+04, 0.18963E+04,
                                 0.19269E+04, 0.19577E+04, 0.19888E+04, 0.20202E+04, 0.20518E+04,
                                 0.20837E+04, 0.21158E+04, 0.21482E+04, 0.21809E+04, 0.22139E+04,
                                 0.22471E+04, 0.22806E+04, 0.23143E+04, 0.23484E+04, 0.23827E+04,
                                 0.24173E+04, 0.24521E+04, 0.24873E+04, 0.25227E+04, 0.25584E+04,
                                 0.25943E+04, 0.26306E+04, 0.26671E+04, 0.27039E+04, 0.27409E+04,
                                 0.27783E+04, 0.28159E+04, 0.28538E+04, 0.28920E+04, 0.29305E+04,
                                 0.29693E+04])


#  --------------- HBr 11: M = 16, I = 2 ---------------------
M = 16
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.)
TIPS_ISO_HASH[(M, I)] = float32([0.42756E+02, 0.59390E+02, 0.76045E+02,
                                 0.92713E+02, 0.10939E+03, 0.12607E+03, 0.14277E+03, 0.15947E+03,
                                 0.17617E+03, 0.19288E+03, 0.20960E+03, 0.22633E+03, 0.24306E+03,
                                 0.25981E+03, 0.27656E+03, 0.29334E+03, 0.31014E+03, 0.32696E+03,
                                 0.34381E+03, 0.36071E+03, 0.37764E+03, 0.39464E+03, 0.41169E+03,
                                 0.42881E+03, 0.44601E+03, 0.46329E+03, 0.48066E+03, 0.49813E+03,
                                 0.51572E+03, 0.53341E+03, 0.55123E+03, 0.56918E+03, 0.58727E+03,
                                 0.60549E+03, 0.62387E+03, 0.64240E+03, 0.66109E+03, 0.67994E+03,
                                 0.69896E+03, 0.71816E+03, 0.73754E+03, 0.75710E+03, 0.77684E+03,
                                 0.79678E+03, 0.81691E+03, 0.83724E+03, 0.85776E+03, 0.87850E+03,
                                 0.89943E+03, 0.92058E+03, 0.94194E+03, 0.96352E+03, 0.98531E+03,
                                 0.10073E+04, 0.10295E+04, 0.10520E+04, 0.10747E+04, 0.10976E+04,
                                 0.11207E+04, 0.11441E+04, 0.11677E+04, 0.11915E+04, 0.12156E+04,
                                 0.12399E+04, 0.12644E+04, 0.12892E+04, 0.13142E+04, 0.13395E+04,
                                 0.13650E+04, 0.13907E+04, 0.14167E+04, 0.14429E+04, 0.14694E+04,
                                 0.14961E+04, 0.15231E+04, 0.15504E+04, 0.15778E+04, 0.16056E+04,
                                 0.16336E+04, 0.16618E+04, 0.16903E+04, 0.17191E+04, 0.17481E+04,
                                 0.17773E+04, 0.18069E+04, 0.18367E+04, 0.18667E+04, 0.18970E+04,
                                 0.19276E+04, 0.19584E+04, 0.19895E+04, 0.20209E+04, 0.20525E+04,
                                 0.20844E+04, 0.21166E+04, 0.21490E+04, 0.21817E+04, 0.22147E+04,
                                 0.22479E+04, 0.22814E+04, 0.23152E+04, 0.23492E+04, 0.23835E+04,
                                 0.24181E+04, 0.24530E+04, 0.24882E+04, 0.25236E+04, 0.25593E+04,
                                 0.25952E+04, 0.26315E+04, 0.26680E+04, 0.27048E+04, 0.27419E+04,
                                 0.27793E+04, 0.28169E+04, 0.28549E+04, 0.28931E+04, 0.29316E+04,
                                 0.29703E+04])


#  --------------- HBr 29: M = 16, I = 3 --------------------- not in TIPS-2011
M = 16
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- HBr 21: M = 16, I = 4 --------------------- not in TIPS-2011
M = 16
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- HI 17: M = 17, I = 1 ---------------------
M = 17
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.82031E+02, 0.11447E+03, 0.14694E+03,
                                 0.17943E+03, 0.21194E+03, 0.24445E+03, 0.27699E+03, 0.30953E+03,
                                 0.34209E+03, 0.37466E+03, 0.40725E+03, 0.43986E+03, 0.47249E+03,
                                 0.50517E+03, 0.53789E+03, 0.57068E+03, 0.60354E+03, 0.63650E+03,
                                 0.66957E+03, 0.70278E+03, 0.73614E+03, 0.76967E+03, 0.80340E+03,
                                 0.83735E+03, 0.87153E+03, 0.90596E+03, 0.94067E+03, 0.97566E+03,
                                 0.10110E+04, 0.10466E+04, 0.10826E+04, 0.11189E+04, 0.11555E+04,
                                 0.11926E+04, 0.12300E+04, 0.12679E+04, 0.13061E+04, 0.13448E+04,
                                 0.13839E+04, 0.14235E+04, 0.14635E+04, 0.15039E+04, 0.15448E+04,
                                 0.15862E+04, 0.16280E+04, 0.16704E+04, 0.17132E+04, 0.17565E+04,
                                 0.18003E+04, 0.18446E+04, 0.18894E+04, 0.19347E+04, 0.19806E+04,
                                 0.20269E+04, 0.20738E+04, 0.21212E+04, 0.21691E+04, 0.22176E+04,
                                 0.22666E+04, 0.23162E+04, 0.23662E+04, 0.24169E+04, 0.24680E+04,
                                 0.25198E+04, 0.25720E+04, 0.26249E+04, 0.26783E+04, 0.27322E+04,
                                 0.27867E+04, 0.28418E+04, 0.28975E+04, 0.29537E+04, 0.30105E+04,
                                 0.30678E+04, 0.31258E+04, 0.31843E+04, 0.32434E+04, 0.33031E+04,
                                 0.33633E+04, 0.34242E+04, 0.34856E+04, 0.35477E+04, 0.36103E+04,
                                 0.36735E+04, 0.37373E+04, 0.38018E+04, 0.38668E+04, 0.39324E+04,
                                 0.39986E+04, 0.40654E+04, 0.41329E+04, 0.42009E+04, 0.42696E+04,
                                 0.43388E+04, 0.44087E+04, 0.44792E+04, 0.45503E+04, 0.46221E+04,
                                 0.46944E+04, 0.47674E+04, 0.48410E+04, 0.49152E+04, 0.49901E+04,
                                 0.50656E+04, 0.51417E+04, 0.52185E+04, 0.52959E+04, 0.53739E+04,
                                 0.54526E+04, 0.55319E+04, 0.56118E+04, 0.56924E+04, 0.57736E+04,
                                 0.58555E+04, 0.59380E+04, 0.60212E+04, 0.61050E+04, 0.61895E+04,
                                 0.62746E+04])


#  --------------- HI 27: M = 17, I = 2 --------------------- not in TIPS-2011
M = 17
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- ClO 56: M = 18, I = 1 ---------------------
M = 18
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.53847E+03, 0.76580E+03, 0.10017E+04,
                                 0.12511E+04, 0.15168E+04, 0.18001E+04, 0.21014E+04, 0.24206E+04,
                                 0.27577E+04, 0.31127E+04, 0.34857E+04, 0.38765E+04, 0.42854E+04,
                                 0.47124E+04, 0.51575E+04, 0.56208E+04, 0.61025E+04, 0.66026E+04,
                                 0.71211E+04, 0.76582E+04, 0.82138E+04, 0.87882E+04, 0.93813E+04,
                                 0.99932E+04, 0.10624E+05, 0.11273E+05, 0.11942E+05, 0.12629E+05,
                                 0.13336E+05, 0.14061E+05, 0.14806E+05, 0.15570E+05, 0.16353E+05,
                                 0.17155E+05, 0.17976E+05, 0.18816E+05, 0.19676E+05, 0.20555E+05,
                                 0.21453E+05, 0.22371E+05, 0.23308E+05, 0.24264E+05, 0.25240E+05,
                                 0.26236E+05, 0.27250E+05, 0.28284E+05, 0.29338E+05, 0.30412E+05,
                                 0.31505E+05, 0.32617E+05, 0.33749E+05, 0.34901E+05, 0.36072E+05,
                                 0.37263E+05, 0.38474E+05, 0.39705E+05, 0.40955E+05, 0.42225E+05,
                                 0.43515E+05, 0.44825E+05, 0.46154E+05, 0.47504E+05, 0.48873E+05,
                                 0.50262E+05, 0.51672E+05, 0.53101E+05, 0.54549E+05, 0.56019E+05,
                                 0.57508E+05, 0.59017E+05, 0.60546E+05, 0.62095E+05, 0.63665E+05,
                                 0.65254E+05, 0.66864E+05, 0.68494E+05, 0.70144E+05, 0.71814E+05,
                                 0.73504E+05, 0.75215E+05, 0.76946E+05, 0.78698E+05, 0.80470E+05,
                                 0.82261E+05, 0.84074E+05, 0.85907E+05, 0.87760E+05, 0.89633E+05,
                                 0.91527E+05, 0.93442E+05, 0.95377E+05, 0.97333E+05, 0.99309E+05,
                                 0.10131E+06, 0.10332E+06, 0.10536E+06, 0.10742E+06, 0.10950E+06,
                                 0.11160E+06, 0.11372E+06, 0.11586E+06, 0.11802E+06, 0.12020E+06,
                                 0.12241E+06, 0.12463E+06, 0.12688E+06, 0.12914E+06, 0.13143E+06,
                                 0.13374E+06, 0.13607E+06, 0.13842E+06, 0.14079E+06, 0.14318E+06,
                                 0.14559E+06, 0.14802E+06, 0.15048E+06, 0.15295E+06, 0.15545E+06,
                                 0.15797E+06])


#  --------------- ClO 76: M = 18, I = 2 ---------------------
M = 18
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.54775E+03, 0.77899E+03, 0.10189E+04,
                                 0.12726E+04, 0.15430E+04, 0.18313E+04, 0.21378E+04, 0.24627E+04,
                                 0.28059E+04, 0.31674E+04, 0.35472E+04, 0.39454E+04, 0.43621E+04,
                                 0.47972E+04, 0.52508E+04, 0.57232E+04, 0.62143E+04, 0.67242E+04,
                                 0.72531E+04, 0.78010E+04, 0.83678E+04, 0.89537E+04, 0.95589E+04,
                                 0.10183E+05, 0.10827E+05, 0.11490E+05, 0.12172E+05, 0.12874E+05,
                                 0.13595E+05, 0.14335E+05, 0.15095E+05, 0.15875E+05, 0.16674E+05,
                                 0.17493E+05, 0.18332E+05, 0.19190E+05, 0.20068E+05, 0.20965E+05,
                                 0.21882E+05, 0.22820E+05, 0.23776E+05, 0.24753E+05, 0.25750E+05,
                                 0.26766E+05, 0.27803E+05, 0.28859E+05, 0.29935E+05, 0.31032E+05,
                                 0.32148E+05, 0.33284E+05, 0.34441E+05, 0.35617E+05, 0.36814E+05,
                                 0.38031E+05, 0.39267E+05, 0.40524E+05, 0.41802E+05, 0.43099E+05,
                                 0.44417E+05, 0.45755E+05, 0.47113E+05, 0.48492E+05, 0.49891E+05,
                                 0.51310E+05, 0.52750E+05, 0.54210E+05, 0.55690E+05, 0.57191E+05,
                                 0.58713E+05, 0.60255E+05, 0.61817E+05, 0.63400E+05, 0.65004E+05,
                                 0.66628E+05, 0.68272E+05, 0.69938E+05, 0.71624E+05, 0.73331E+05,
                                 0.75058E+05, 0.76806E+05, 0.78575E+05, 0.80364E+05, 0.82175E+05,
                                 0.84006E+05, 0.85858E+05, 0.87731E+05, 0.89625E+05, 0.91539E+05,
                                 0.93475E+05, 0.95431E+05, 0.97409E+05, 0.99407E+05, 0.10143E+06,
                                 0.10347E+06, 0.10553E+06, 0.10761E+06, 0.10972E+06, 0.11184E+06,
                                 0.11399E+06, 0.11615E+06, 0.11834E+06, 0.12055E+06, 0.12278E+06,
                                 0.12503E+06, 0.12731E+06, 0.12960E+06, 0.13192E+06, 0.13425E+06,
                                 0.13661E+06, 0.13899E+06, 0.14139E+06, 0.14382E+06, 0.14626E+06,
                                 0.14873E+06, 0.15121E+06, 0.15372E+06, 0.15625E+06, 0.15880E+06,
                                 0.16138E+06])


#  --------------- OCS 622: M = 19, I = 1 ---------------------
M = 19
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.20609E+03, 0.29199E+03, 0.37861E+03,
                                 0.46737E+03, 0.56024E+03, 0.65929E+03, 0.76649E+03, 0.88361E+03,
                                 0.10123E+04, 0.11541E+04, 0.13105E+04, 0.14829E+04, 0.16728E+04,
                                 0.18818E+04, 0.21113E+04, 0.23629E+04, 0.26383E+04, 0.29391E+04,
                                 0.32672E+04, 0.36245E+04, 0.40128E+04, 0.44343E+04, 0.48911E+04,
                                 0.53853E+04, 0.59193E+04, 0.64956E+04, 0.71166E+04, 0.77849E+04,
                                 0.85033E+04, 0.92746E+04, 0.10102E+05, 0.10988E+05, 0.11936E+05,
                                 0.12949E+05, 0.14032E+05, 0.15186E+05, 0.16416E+05, 0.17726E+05,
                                 0.19120E+05, 0.20601E+05, 0.22173E+05, 0.23842E+05, 0.25611E+05,
                                 0.27484E+05, 0.29468E+05, 0.31566E+05, 0.33783E+05, 0.36124E+05,
                                 0.38595E+05, 0.41202E+05, 0.43949E+05, 0.46842E+05, 0.49888E+05,
                                 0.53092E+05, 0.56460E+05, 0.59999E+05, 0.63716E+05, 0.67616E+05,
                                 0.71708E+05, 0.75997E+05, 0.80491E+05, 0.85197E+05, 0.90124E+05,
                                 0.95278E+05, 0.10067E+06, 0.10630E+06, 0.11219E+06, 0.11833E+06,
                                 0.12475E+06, 0.13144E+06, 0.13842E+06, 0.14570E+06, 0.15328E+06,
                                 0.16117E+06, 0.16940E+06, 0.17795E+06, 0.18686E+06, 0.19611E+06,
                                 0.20574E+06, 0.21574E+06, 0.22613E+06, 0.23692E+06, 0.24813E+06,
                                 0.25975E+06, 0.27182E+06, 0.28433E+06, 0.29730E+06, 0.31074E+06,
                                 0.32467E+06, 0.33909E+06, 0.35403E+06, 0.36950E+06, 0.38551E+06,
                                 0.40207E+06, 0.41920E+06, 0.43691E+06, 0.45522E+06, 0.47415E+06,
                                 0.49370E+06, 0.51390E+06, 0.53476E+06, 0.55629E+06, 0.57852E+06,
                                 0.60146E+06, 0.62513E+06, 0.64954E+06, 0.67471E+06, 0.70067E+06,
                                 0.72742E+06, 0.75499E+06, 0.78339E+06, 0.81265E+06, 0.84279E+06,
                                 0.87381E+06, 0.90576E+06, 0.93863E+06, 0.97246E+06, 0.10073E+07,
                                 0.10431E+07])


#  --------------- OCS 624: M = 19, I = 2 ---------------------
M = 19
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.21125E+03, 0.29930E+03, 0.38809E+03,
                                 0.47911E+03, 0.57437E+03, 0.67603E+03, 0.78610E+03, 0.90643E+03,
                                 0.10387E+04, 0.11846E+04, 0.13456E+04, 0.15231E+04, 0.17188E+04,
                                 0.19342E+04, 0.21709E+04, 0.24304E+04, 0.27145E+04, 0.30250E+04,
                                 0.33638E+04, 0.37328E+04, 0.41339E+04, 0.45694E+04, 0.50415E+04,
                                 0.55524E+04, 0.61045E+04, 0.67004E+04, 0.73427E+04, 0.80340E+04,
                                 0.87773E+04, 0.95755E+04, 0.10432E+05, 0.11349E+05, 0.12330E+05,
                                 0.13380E+05, 0.14500E+05, 0.15696E+05, 0.16970E+05, 0.18327E+05,
                                 0.19770E+05, 0.21305E+05, 0.22934E+05, 0.24663E+05, 0.26497E+05,
                                 0.28439E+05, 0.30495E+05, 0.32669E+05, 0.34968E+05, 0.37396E+05,
                                 0.39958E+05, 0.42661E+05, 0.45510E+05, 0.48511E+05, 0.51669E+05,
                                 0.54993E+05, 0.58487E+05, 0.62159E+05, 0.66014E+05, 0.70061E+05,
                                 0.74306E+05, 0.78757E+05, 0.83421E+05, 0.88305E+05, 0.93418E+05,
                                 0.98767E+05, 0.10436E+06, 0.11021E+06, 0.11632E+06, 0.12270E+06,
                                 0.12936E+06, 0.13631E+06, 0.14355E+06, 0.15111E+06, 0.15898E+06,
                                 0.16718E+06, 0.17572E+06, 0.18460E+06, 0.19385E+06, 0.20346E+06,
                                 0.21346E+06, 0.22385E+06, 0.23464E+06, 0.24585E+06, 0.25748E+06,
                                 0.26956E+06, 0.28209E+06, 0.29509E+06, 0.30856E+06, 0.32252E+06,
                                 0.33699E+06, 0.35198E+06, 0.36750E+06, 0.38357E+06, 0.40020E+06,
                                 0.41741E+06, 0.43521E+06, 0.45362E+06, 0.47264E+06, 0.49231E+06,
                                 0.51263E+06, 0.53362E+06, 0.55529E+06, 0.57768E+06, 0.60078E+06,
                                 0.62462E+06, 0.64922E+06, 0.67459E+06, 0.70075E+06, 0.72773E+06,
                                 0.75554E+06, 0.78419E+06, 0.81372E+06, 0.84413E+06, 0.87546E+06,
                                 0.90771E+06, 0.94092E+06, 0.97509E+06, 0.10103E+07, 0.10464E+07,
                                 0.10837E+07])


#  --------------- OCS 632: M = 19, I = 3 ---------------------
M = 19
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.41351E+03, 0.58591E+03, 0.76004E+03,
                                 0.93907E+03, 0.11273E+04, 0.13289E+04, 0.15481E+04, 0.17884E+04,
                                 0.20533E+04, 0.23459E+04, 0.26692E+04, 0.30264E+04, 0.34205E+04,
                                 0.38547E+04, 0.43323E+04, 0.48565E+04, 0.54309E+04, 0.60592E+04,
                                 0.67451E+04, 0.74928E+04, 0.83064E+04, 0.91903E+04, 0.10149E+05,
                                 0.11187E+05, 0.12310E+05, 0.13523E+05, 0.14831E+05, 0.16240E+05,
                                 0.17756E+05, 0.19384E+05, 0.21132E+05, 0.23005E+05, 0.25011E+05,
                                 0.27157E+05, 0.29449E+05, 0.31896E+05, 0.34506E+05, 0.37286E+05,
                                 0.40245E+05, 0.43392E+05, 0.46735E+05, 0.50284E+05, 0.54048E+05,
                                 0.58038E+05, 0.62263E+05, 0.66733E+05, 0.71460E+05, 0.76455E+05,
                                 0.81728E+05, 0.87292E+05, 0.93159E+05, 0.99341E+05, 0.10585E+06,
                                 0.11270E+06, 0.11991E+06, 0.12748E+06, 0.13543E+06, 0.14378E+06,
                                 0.15255E+06, 0.16174E+06, 0.17137E+06, 0.18146E+06, 0.19202E+06,
                                 0.20308E+06, 0.21465E+06, 0.22674E+06, 0.23937E+06, 0.25257E+06,
                                 0.26635E+06, 0.28073E+06, 0.29573E+06, 0.31137E+06, 0.32767E+06,
                                 0.34466E+06, 0.36235E+06, 0.38076E+06, 0.39992E+06, 0.41985E+06,
                                 0.44057E+06, 0.46211E+06, 0.48450E+06, 0.50775E+06, 0.53189E+06,
                                 0.55695E+06, 0.58295E+06, 0.60992E+06, 0.63789E+06, 0.66688E+06,
                                 0.69693E+06, 0.72806E+06, 0.76030E+06, 0.79368E+06, 0.82823E+06,
                                 0.86399E+06, 0.90097E+06, 0.93923E+06, 0.97878E+06, 0.10197E+07,
                                 0.10619E+07, 0.11056E+07, 0.11506E+07, 0.11972E+07, 0.12453E+07,
                                 0.12949E+07, 0.13460E+07, 0.13988E+07, 0.14533E+07, 0.15094E+07,
                                 0.15673E+07, 0.16270E+07, 0.16884E+07, 0.17518E+07, 0.18170E+07,
                                 0.18842E+07, 0.19533E+07, 0.20245E+07, 0.20978E+07, 0.21732E+07,
                                 0.22507E+07])


#  --------------- OCS 623: M = 19, I = 4 ---------------------
M = 19
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.83485E+03, 0.11828E+04, 0.15337E+04,
                                 0.18934E+04, 0.22697E+04, 0.26712E+04, 0.31059E+04, 0.35809E+04,
                                 0.41030E+04, 0.46785E+04, 0.53133E+04, 0.60135E+04, 0.67850E+04,
                                 0.76338E+04, 0.85663E+04, 0.95888E+04, 0.10708E+05, 0.11931E+05,
                                 0.13265E+05, 0.14718E+05, 0.16298E+05, 0.18012E+05, 0.19870E+05,
                                 0.21881E+05, 0.24054E+05, 0.26399E+05, 0.28926E+05, 0.31646E+05,
                                 0.34570E+05, 0.37710E+05, 0.41077E+05, 0.44685E+05, 0.48545E+05,
                                 0.52672E+05, 0.57078E+05, 0.61780E+05, 0.66790E+05, 0.72125E+05,
                                 0.77801E+05, 0.83833E+05, 0.90239E+05, 0.97036E+05, 0.10424E+06,
                                 0.11188E+06, 0.11996E+06, 0.12850E+06, 0.13754E+06, 0.14708E+06,
                                 0.15715E+06, 0.16777E+06, 0.17896E+06, 0.19076E+06, 0.20317E+06,
                                 0.21623E+06, 0.22996E+06, 0.24438E+06, 0.25953E+06, 0.27543E+06,
                                 0.29211E+06, 0.30959E+06, 0.32791E+06, 0.34710E+06, 0.36718E+06,
                                 0.38820E+06, 0.41017E+06, 0.43314E+06, 0.45713E+06, 0.48219E+06,
                                 0.50835E+06, 0.53564E+06, 0.56409E+06, 0.59376E+06, 0.62468E+06,
                                 0.65688E+06, 0.69041E+06, 0.72530E+06, 0.76161E+06, 0.79937E+06,
                                 0.83862E+06, 0.87941E+06, 0.92179E+06, 0.96581E+06, 0.10115E+07,
                                 0.10589E+07, 0.11081E+07, 0.11591E+07, 0.12120E+07, 0.12669E+07,
                                 0.13237E+07, 0.13825E+07, 0.14435E+07, 0.15066E+07, 0.15718E+07,
                                 0.16394E+07, 0.17093E+07, 0.17815E+07, 0.18562E+07, 0.19334E+07,
                                 0.20132E+07, 0.20956E+07, 0.21807E+07, 0.22685E+07, 0.23592E+07,
                                 0.24528E+07, 0.25494E+07, 0.26490E+07, 0.27517E+07, 0.28576E+07,
                                 0.29667E+07, 0.30792E+07, 0.31951E+07, 0.33145E+07, 0.34374E+07,
                                 0.35640E+07, 0.36943E+07, 0.38285E+07, 0.39665E+07, 0.41085E+07,
                                 0.42546E+07])


#  --------------- OCS 822: M = 19, I = 5 ---------------------
M = 19
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.21967E+03, 0.31126E+03, 0.40370E+03,
                                 0.49862E+03, 0.59823E+03, 0.70481E+03, 0.82050E+03, 0.94724E+03,
                                 0.10868E+04, 0.12409E+04, 0.14112E+04, 0.15993E+04, 0.18067E+04,
                                 0.20353E+04, 0.22866E+04, 0.25624E+04, 0.28645E+04, 0.31950E+04,
                                 0.35558E+04, 0.39490E+04, 0.43767E+04, 0.48413E+04, 0.53452E+04,
                                 0.58909E+04, 0.64810E+04, 0.71182E+04, 0.78053E+04, 0.85454E+04,
                                 0.93413E+04, 0.10196E+05, 0.11114E+05, 0.12098E+05, 0.13151E+05,
                                 0.14277E+05, 0.15480E+05, 0.16764E+05, 0.18133E+05, 0.19592E+05,
                                 0.21144E+05, 0.22794E+05, 0.24548E+05, 0.26409E+05, 0.28383E+05,
                                 0.30475E+05, 0.32689E+05, 0.35033E+05, 0.37511E+05, 0.40128E+05,
                                 0.42892E+05, 0.45808E+05, 0.48882E+05, 0.52121E+05, 0.55532E+05,
                                 0.59121E+05, 0.62895E+05, 0.66861E+05, 0.71028E+05, 0.75402E+05,
                                 0.79991E+05, 0.84803E+05, 0.89847E+05, 0.95130E+05, 0.10066E+06,
                                 0.10645E+06, 0.11251E+06, 0.11883E+06, 0.12545E+06, 0.13236E+06,
                                 0.13957E+06, 0.14710E+06, 0.15495E+06, 0.16313E+06, 0.17166E+06,
                                 0.18055E+06, 0.18980E+06, 0.19944E+06, 0.20946E+06, 0.21989E+06,
                                 0.23073E+06, 0.24200E+06, 0.25371E+06, 0.26587E+06, 0.27850E+06,
                                 0.29161E+06, 0.30521E+06, 0.31931E+06, 0.33394E+06, 0.34910E+06,
                                 0.36482E+06, 0.38109E+06, 0.39795E+06, 0.41541E+06, 0.43348E+06,
                                 0.45217E+06, 0.47151E+06, 0.49151E+06, 0.51219E+06, 0.53356E+06,
                                 0.55565E+06, 0.57847E+06, 0.60204E+06, 0.62637E+06, 0.65149E+06,
                                 0.67742E+06, 0.70417E+06, 0.73176E+06, 0.76023E+06, 0.78957E+06,
                                 0.81982E+06, 0.85100E+06, 0.88313E+06, 0.91622E+06, 0.95031E+06,
                                 0.98541E+06, 0.10216E+07, 0.10587E+07, 0.10970E+07, 0.11364E+07,
                                 0.11769E+07])


#  --------------- H2CO 126: M = 20, I = 2 ---------------------
M = 20
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.25934E+03, 0.43623E+03, 0.64143E+03,
                                 0.87152E+03, 0.11241E+04, 0.13975E+04, 0.16906E+04, 0.20029E+04,
                                 0.23344E+04, 0.26857E+04, 0.30577E+04, 0.34518E+04, 0.38698E+04,
                                 0.43138E+04, 0.47860E+04, 0.52890E+04, 0.58256E+04, 0.63985E+04,
                                 0.70109E+04, 0.76660E+04, 0.83673E+04, 0.91184E+04, 0.99230E+04,
                                 0.10785E+05, 0.11710E+05, 0.12700E+05, 0.13762E+05, 0.14900E+05,
                                 0.16119E+05, 0.17425E+05, 0.18823E+05, 0.20320E+05, 0.21923E+05,
                                 0.23637E+05, 0.25471E+05, 0.27432E+05, 0.29527E+05, 0.31765E+05,
                                 0.34155E+05, 0.36706E+05, 0.39428E+05, 0.42330E+05, 0.45424E+05,
                                 0.48720E+05, 0.52231E+05, 0.55968E+05, 0.59945E+05, 0.64175E+05,
                                 0.68672E+05, 0.73450E+05, 0.78526E+05, 0.83915E+05, 0.89634E+05,
                                 0.95701E+05, 0.10213E+06, 0.10895E+06, 0.11618E+06, 0.12383E+06,
                                 0.13193E+06, 0.14049E+06, 0.14956E+06, 0.15914E+06, 0.16927E+06,
                                 0.17997E+06, 0.19127E+06, 0.20320E+06, 0.21578E+06, 0.22906E+06,
                                 0.24306E+06, 0.25782E+06, 0.27336E+06, 0.28974E+06, 0.30698E+06,
                                 0.32513E+06, 0.34422E+06, 0.36430E+06, 0.38542E+06, 0.40761E+06,
                                 0.43093E+06, 0.45542E+06, 0.48114E+06, 0.50813E+06, 0.53646E+06,
                                 0.56617E+06, 0.59733E+06, 0.63000E+06, 0.66423E+06, 0.70010E+06,
                                 0.73767E+06, 0.77701E+06, 0.81818E+06, 0.86127E+06, 0.90635E+06,
                                 0.95349E+06, 0.10028E+07, 0.10543E+07, 0.11082E+07, 0.11644E+07,
                                 0.12232E+07, 0.12845E+07, 0.13485E+07, 0.14154E+07, 0.14851E+07,
                                 0.15578E+07, 0.16337E+07, 0.17127E+07, 0.17952E+07, 0.18810E+07,
                                 0.19705E+07, 0.20637E+07, 0.21607E+07, 0.22617E+07, 0.23669E+07,
                                 0.24763E+07, 0.25901E+07, 0.27085E+07, 0.28316E+07, 0.29596E+07,
                                 0.30926E+07])


#  --------------- H2CO 136: M = 20, I = 2 ---------------------
M = 20
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.53173E+03, 0.89447E+03, 0.13153E+04,
                                 0.17871E+04, 0.23051E+04, 0.28658E+04, 0.34669E+04, 0.41073E+04,
                                 0.47872E+04, 0.55074E+04, 0.62702E+04, 0.70785E+04, 0.79357E+04,
                                 0.88462E+04, 0.98147E+04, 0.10846E+05, 0.11946E+05, 0.13121E+05,
                                 0.14377E+05, 0.15721E+05, 0.17159E+05, 0.18699E+05, 0.20349E+05,
                                 0.22118E+05, 0.24013E+05, 0.26045E+05, 0.28222E+05, 0.30555E+05,
                                 0.33055E+05, 0.35733E+05, 0.38601E+05, 0.41671E+05, 0.44958E+05,
                                 0.48474E+05, 0.52235E+05, 0.56255E+05, 0.60552E+05, 0.65142E+05,
                                 0.70043E+05, 0.75275E+05, 0.80856E+05, 0.86808E+05, 0.93152E+05,
                                 0.99913E+05, 0.10711E+06, 0.11478E+06, 0.12293E+06, 0.13161E+06,
                                 0.14083E+06, 0.15063E+06, 0.16104E+06, 0.17209E+06, 0.18382E+06,
                                 0.19626E+06, 0.20945E+06, 0.22343E+06, 0.23825E+06, 0.25394E+06,
                                 0.27054E+06, 0.28812E+06, 0.30671E+06, 0.32636E+06, 0.34713E+06,
                                 0.36907E+06, 0.39224E+06, 0.41671E+06, 0.44252E+06, 0.46975E+06,
                                 0.49845E+06, 0.52872E+06, 0.56060E+06, 0.59418E+06, 0.62954E+06,
                                 0.66676E+06, 0.70591E+06, 0.74710E+06, 0.79040E+06, 0.83591E+06,
                                 0.88373E+06, 0.93395E+06, 0.98669E+06, 0.10421E+07, 0.11001E+07,
                                 0.11611E+07, 0.12250E+07, 0.12920E+07, 0.13622E+07, 0.14357E+07,
                                 0.15128E+07, 0.15934E+07, 0.16779E+07, 0.17662E+07, 0.18587E+07,
                                 0.19554E+07, 0.20565E+07, 0.21621E+07, 0.22725E+07, 0.23879E+07,
                                 0.25084E+07, 0.26342E+07, 0.27655E+07, 0.29026E+07, 0.30456E+07,
                                 0.31947E+07, 0.33502E+07, 0.35124E+07, 0.36814E+07, 0.38575E+07,
                                 0.40410E+07, 0.42321E+07, 0.44311E+07, 0.46382E+07, 0.48538E+07,
                                 0.50782E+07, 0.53116E+07, 0.55544E+07, 0.58068E+07, 0.60693E+07,
                                 0.63421E+07])


#  --------------- H2CO 128: M = 20, I = 3 ---------------------
M = 20
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.27198E+03, 0.45755E+03, 0.67282E+03,
                                 0.91421E+03, 0.11792E+04, 0.14660E+04, 0.17735E+04, 0.21012E+04,
                                 0.24490E+04, 0.28175E+04, 0.32077E+04, 0.36212E+04, 0.40598E+04,
                                 0.45256E+04, 0.50211E+04, 0.55488E+04, 0.61116E+04, 0.67127E+04,
                                 0.73552E+04, 0.80426E+04, 0.87783E+04, 0.95663E+04, 0.10410E+05,
                                 0.11315E+05, 0.12285E+05, 0.13324E+05, 0.14438E+05, 0.15632E+05,
                                 0.16911E+05, 0.18281E+05, 0.19748E+05, 0.21319E+05, 0.23000E+05,
                                 0.24799E+05, 0.26723E+05, 0.28780E+05, 0.30978E+05, 0.33326E+05,
                                 0.35834E+05, 0.38510E+05, 0.41365E+05, 0.44410E+05, 0.47656E+05,
                                 0.51115E+05, 0.54798E+05, 0.58719E+05, 0.62891E+05, 0.67329E+05,
                                 0.72047E+05, 0.77060E+05, 0.82385E+05, 0.88039E+05, 0.94039E+05,
                                 0.10040E+06, 0.10715E+06, 0.11431E+06, 0.12189E+06, 0.12991E+06,
                                 0.13841E+06, 0.14740E+06, 0.15691E+06, 0.16696E+06, 0.17759E+06,
                                 0.18882E+06, 0.20067E+06, 0.21318E+06, 0.22639E+06, 0.24032E+06,
                                 0.25501E+06, 0.27049E+06, 0.28680E+06, 0.30398E+06, 0.32207E+06,
                                 0.34111E+06, 0.36114E+06, 0.38221E+06, 0.40436E+06, 0.42765E+06,
                                 0.45211E+06, 0.47781E+06, 0.50479E+06, 0.53311E+06, 0.56283E+06,
                                 0.59400E+06, 0.62669E+06, 0.66097E+06, 0.69688E+06, 0.73451E+06,
                                 0.77393E+06, 0.81520E+06, 0.85840E+06, 0.90360E+06, 0.95090E+06,
                                 0.10004E+07, 0.10521E+07, 0.11061E+07, 0.11626E+07, 0.12216E+07,
                                 0.12833E+07, 0.13476E+07, 0.14148E+07, 0.14849E+07, 0.15581E+07,
                                 0.16344E+07, 0.17140E+07, 0.17969E+07, 0.18834E+07, 0.19735E+07,
                                 0.20674E+07, 0.21651E+07, 0.22669E+07, 0.23729E+07, 0.24832E+07,
                                 0.25980E+07, 0.27174E+07, 0.28416E+07, 0.29708E+07, 0.31050E+07,
                                 0.32446E+07])


#  --------------- HOCl 165: M = 21, I = 1 ---------------------
M = 21
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.)
TIPS_ISO_HASH[(M, I)] = float32([0.17041E+04, 0.28708E+04, 0.42250E+04,
                                 0.57456E+04, 0.74211E+04, 0.92470E+04, 0.11225E+05, 0.13359E+05,
                                 0.15657E+05, 0.18129E+05, 0.20785E+05, 0.23637E+05, 0.26696E+05,
                                 0.29974E+05, 0.33484E+05, 0.37239E+05, 0.41252E+05, 0.45536E+05,
                                 0.50105E+05, 0.54973E+05, 0.60152E+05, 0.65659E+05, 0.71507E+05,
                                 0.77711E+05, 0.84286E+05, 0.91249E+05, 0.98614E+05, 0.10640E+06,
                                 0.11462E+06, 0.12330E+06, 0.13244E+06, 0.14208E+06, 0.15222E+06,
                                 0.16289E+06, 0.17411E+06, 0.18589E+06, 0.19825E+06, 0.21123E+06,
                                 0.22483E+06, 0.23908E+06, 0.25400E+06, 0.26962E+06, 0.28596E+06,
                                 0.30303E+06, 0.32087E+06, 0.33950E+06, 0.35895E+06, 0.37923E+06,
                                 0.40038E+06, 0.42243E+06, 0.44539E+06, 0.46930E+06, 0.49419E+06,
                                 0.52008E+06, 0.54700E+06, 0.57498E+06, 0.60406E+06, 0.63426E+06,
                                 0.66562E+06, 0.69816E+06, 0.73192E+06, 0.76692E+06, 0.80322E+06,
                                 0.84083E+06, 0.87979E+06, 0.92014E+06, 0.96192E+06, 0.10052E+07,
                                 0.10499E+07, 0.10961E+07, 0.11440E+07, 0.11934E+07, 0.12445E+07,
                                 0.12973E+07, 0.13518E+07, 0.14081E+07, 0.14661E+07, 0.15261E+07,
                                 0.15879E+07, 0.16516E+07, 0.17174E+07, 0.17851E+07, 0.18550E+07,
                                 0.19269E+07, 0.20010E+07, 0.20773E+07, 0.21559E+07, 0.22367E+07,
                                 0.23200E+07, 0.24056E+07, 0.24936E+07, 0.25842E+07, 0.26773E+07,
                                 0.27730E+07, 0.28714E+07, 0.29724E+07, 0.30763E+07, 0.31829E+07,
                                 0.32924E+07, 0.34049E+07, 0.35203E+07, 0.36387E+07, 0.37603E+07,
                                 0.38850E+07, 0.40129E+07, 0.41441E+07, 0.42786E+07, 0.44165E+07,
                                 0.45579E+07, 0.47028E+07, 0.48512E+07, 0.50033E+07, 0.51592E+07,
                                 0.53187E+07, 0.54822E+07, 0.56495E+07, 0.58208E+07, 0.59961E+07,
                                 0.61755E+07])


#  --------------- HOCl 167: M = 21, I = 2 ---------------------
M = 21
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.)
TIPS_ISO_HASH[(M, I)] = float32([0.17342E+04, 0.29215E+04, 0.42998E+04,
                                 0.58473E+04, 0.75524E+04, 0.94107E+04, 0.11423E+05, 0.13595E+05,
                                 0.15935E+05, 0.18450E+05, 0.21154E+05, 0.24056E+05, 0.27168E+05,
                                 0.30505E+05, 0.34077E+05, 0.37899E+05, 0.41983E+05, 0.46343E+05,
                                 0.50993E+05, 0.55947E+05, 0.61218E+05, 0.66822E+05, 0.72774E+05,
                                 0.79088E+05, 0.85780E+05, 0.92866E+05, 0.10036E+06, 0.10829E+06,
                                 0.11665E+06, 0.12548E+06, 0.13479E+06, 0.14460E+06, 0.15492E+06,
                                 0.16578E+06, 0.17719E+06, 0.18918E+06, 0.20177E+06, 0.21497E+06,
                                 0.22881E+06, 0.24332E+06, 0.25851E+06, 0.27440E+06, 0.29102E+06,
                                 0.30840E+06, 0.32656E+06, 0.34552E+06, 0.36531E+06, 0.38595E+06,
                                 0.40748E+06, 0.42991E+06, 0.45328E+06, 0.47762E+06, 0.50295E+06,
                                 0.52929E+06, 0.55669E+06, 0.58517E+06, 0.61477E+06, 0.64550E+06,
                                 0.67741E+06, 0.71053E+06, 0.74489E+06, 0.78052E+06, 0.81745E+06,
                                 0.85573E+06, 0.89539E+06, 0.93645E+06, 0.97897E+06, 0.10230E+07,
                                 0.10685E+07, 0.11156E+07, 0.11643E+07, 0.12146E+07, 0.12666E+07,
                                 0.13203E+07, 0.13757E+07, 0.14330E+07, 0.14921E+07, 0.15531E+07,
                                 0.16160E+07, 0.16809E+07, 0.17478E+07, 0.18168E+07, 0.18878E+07,
                                 0.19611E+07, 0.20365E+07, 0.21141E+07, 0.21941E+07, 0.22764E+07,
                                 0.23611E+07, 0.24482E+07, 0.25378E+07, 0.26300E+07, 0.27248E+07,
                                 0.28222E+07, 0.29223E+07, 0.30251E+07, 0.31308E+07, 0.32393E+07,
                                 0.33508E+07, 0.34652E+07, 0.35827E+07, 0.37032E+07, 0.38269E+07,
                                 0.39539E+07, 0.40840E+07, 0.42176E+07, 0.43545E+07, 0.44948E+07,
                                 0.46387E+07, 0.47861E+07, 0.49372E+07, 0.50920E+07, 0.52506E+07,
                                 0.54130E+07, 0.55793E+07, 0.57496E+07, 0.59239E+07, 0.61024E+07,
                                 0.62850E+07])


#  --------------- N2 44: M = 22, I = 1 ---------------------
M = 22
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.95487E+02, 0.13466E+03, 0.17386E+03,
                                 0.21307E+03, 0.25230E+03, 0.29154E+03, 0.33080E+03, 0.37008E+03,
                                 0.40937E+03, 0.44868E+03, 0.48800E+03, 0.52736E+03, 0.56674E+03,
                                 0.60616E+03, 0.64562E+03, 0.68515E+03, 0.72475E+03, 0.76445E+03,
                                 0.80426E+03, 0.84420E+03, 0.88430E+03, 0.92457E+03, 0.96505E+03,
                                 0.10057E+04, 0.10467E+04, 0.10879E+04, 0.11293E+04, 0.11711E+04,
                                 0.12132E+04, 0.12556E+04, 0.12984E+04, 0.13416E+04, 0.13851E+04,
                                 0.14291E+04, 0.14734E+04, 0.15182E+04, 0.15635E+04, 0.16091E+04,
                                 0.16553E+04, 0.17019E+04, 0.17490E+04, 0.17965E+04, 0.18446E+04,
                                 0.18932E+04, 0.19422E+04, 0.19918E+04, 0.20419E+04, 0.20926E+04,
                                 0.21437E+04, 0.21954E+04, 0.22477E+04, 0.23004E+04, 0.23538E+04,
                                 0.24077E+04, 0.24621E+04, 0.25171E+04, 0.25727E+04, 0.26288E+04,
                                 0.26856E+04, 0.27428E+04, 0.28007E+04, 0.28591E+04, 0.29181E+04,
                                 0.29777E+04, 0.30379E+04, 0.30986E+04, 0.31600E+04, 0.32219E+04,
                                 0.32844E+04, 0.33475E+04, 0.34112E+04, 0.34755E+04, 0.35404E+04,
                                 0.36059E+04, 0.36720E+04, 0.37387E+04, 0.38060E+04, 0.38739E+04,
                                 0.39424E+04, 0.40115E+04, 0.40812E+04, 0.41515E+04, 0.42224E+04,
                                 0.42939E+04, 0.43661E+04, 0.44388E+04, 0.45122E+04, 0.45861E+04,
                                 0.46607E+04, 0.47359E+04, 0.48117E+04, 0.48882E+04, 0.49652E+04,
                                 0.50428E+04, 0.51211E+04, 0.52000E+04, 0.52795E+04, 0.53596E+04,
                                 0.54404E+04, 0.55217E+04, 0.56037E+04, 0.56863E+04, 0.57695E+04,
                                 0.58533E+04, 0.59378E+04, 0.60229E+04, 0.61086E+04, 0.61950E+04,
                                 0.62819E+04, 0.63695E+04, 0.64577E+04, 0.65465E+04, 0.66360E+04,
                                 0.67261E+04, 0.68168E+04, 0.69081E+04, 0.70001E+04, 0.70927E+04,
                                 0.71859E+04])


#  --------------- N2 45: M = 22, I = 2 --------------------- not in TIPS-2011
M = 22
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- HCN 124: M = 23, I = 1 ---------------------
M = 23
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.17143E+03, 0.24209E+03, 0.31285E+03,
                                 0.38392E+03, 0.45582E+03, 0.52929E+03, 0.60515E+03, 0.68424E+03,
                                 0.76731E+03, 0.85505E+03, 0.94805E+03, 0.10468E+04, 0.11519E+04,
                                 0.12637E+04, 0.13826E+04, 0.15090E+04, 0.16435E+04, 0.17863E+04,
                                 0.19378E+04, 0.20985E+04, 0.22689E+04, 0.24492E+04, 0.26401E+04,
                                 0.28418E+04, 0.30550E+04, 0.32801E+04, 0.35176E+04, 0.37680E+04,
                                 0.40318E+04, 0.43097E+04, 0.46021E+04, 0.49097E+04, 0.52330E+04,
                                 0.55727E+04, 0.59294E+04, 0.63038E+04, 0.66964E+04, 0.71081E+04,
                                 0.75396E+04, 0.79915E+04, 0.84646E+04, 0.89596E+04, 0.94774E+04,
                                 0.10019E+05, 0.10585E+05, 0.11176E+05, 0.11793E+05, 0.12437E+05,
                                 0.13108E+05, 0.13809E+05, 0.14540E+05, 0.15301E+05, 0.16094E+05,
                                 0.16919E+05, 0.17779E+05, 0.18673E+05, 0.19603E+05, 0.20570E+05,
                                 0.21575E+05, 0.22619E+05, 0.23704E+05, 0.24831E+05, 0.26000E+05,
                                 0.27213E+05, 0.28472E+05, 0.29778E+05, 0.31131E+05, 0.32534E+05,
                                 0.33987E+05, 0.35493E+05, 0.37052E+05, 0.38666E+05, 0.40336E+05,
                                 0.42064E+05, 0.43852E+05, 0.45701E+05, 0.47612E+05, 0.49587E+05,
                                 0.51629E+05, 0.53738E+05, 0.55916E+05, 0.58165E+05, 0.60486E+05,
                                 0.62883E+05, 0.65355E+05, 0.67905E+05, 0.70536E+05, 0.73249E+05,
                                 0.76045E+05, 0.78927E+05, 0.81897E+05, 0.84957E+05, 0.88108E+05,
                                 0.91354E+05, 0.94696E+05, 0.98136E+05, 0.10168E+06, 0.10532E+06,
                                 0.10907E+06, 0.11292E+06, 0.11689E+06, 0.12096E+06, 0.12516E+06,
                                 0.12946E+06, 0.13389E+06, 0.13844E+06, 0.14311E+06, 0.14791E+06,
                                 0.15284E+06, 0.15790E+06, 0.16310E+06, 0.16843E+06, 0.17391E+06,
                                 0.17953E+06, 0.18529E+06, 0.19120E+06, 0.19726E+06, 0.20348E+06,
                                 0.20986E+06])


#  --------------- HCN 134: M = 23, I = 2 ---------------------
M = 23
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.35186E+03, 0.49693E+03, 0.64221E+03,
                                 0.78815E+03, 0.93585E+03, 0.10868E+04, 0.12428E+04, 0.14056E+04,
                                 0.15766E+04, 0.17574E+04, 0.19491E+04, 0.21528E+04, 0.23695E+04,
                                 0.26002E+04, 0.28457E+04, 0.31068E+04, 0.33845E+04, 0.36795E+04,
                                 0.39926E+04, 0.43249E+04, 0.46770E+04, 0.50500E+04, 0.54447E+04,
                                 0.58621E+04, 0.63032E+04, 0.67690E+04, 0.72606E+04, 0.77789E+04,
                                 0.83252E+04, 0.89005E+04, 0.95062E+04, 0.10143E+05, 0.10813E+05,
                                 0.11517E+05, 0.12256E+05, 0.13032E+05, 0.13846E+05, 0.14699E+05,
                                 0.15593E+05, 0.16530E+05, 0.17511E+05, 0.18538E+05, 0.19612E+05,
                                 0.20734E+05, 0.21908E+05, 0.23134E+05, 0.24414E+05, 0.25750E+05,
                                 0.27145E+05, 0.28599E+05, 0.30115E+05, 0.31694E+05, 0.33340E+05,
                                 0.35054E+05, 0.36838E+05, 0.38694E+05, 0.40625E+05, 0.42633E+05,
                                 0.44720E+05, 0.46889E+05, 0.49142E+05, 0.51481E+05, 0.53910E+05,
                                 0.56430E+05, 0.59045E+05, 0.61757E+05, 0.64568E+05, 0.67482E+05,
                                 0.70502E+05, 0.73630E+05, 0.76869E+05, 0.80223E+05, 0.83694E+05,
                                 0.87285E+05, 0.91000E+05, 0.94843E+05, 0.98815E+05, 0.10292E+06,
                                 0.10716E+06, 0.11155E+06, 0.11608E+06, 0.12075E+06, 0.12558E+06,
                                 0.13056E+06, 0.13570E+06, 0.14100E+06, 0.14647E+06, 0.15211E+06,
                                 0.15793E+06, 0.16392E+06, 0.17009E+06, 0.17646E+06, 0.18301E+06,
                                 0.18976E+06, 0.19671E+06, 0.20387E+06, 0.21123E+06, 0.21881E+06,
                                 0.22660E+06, 0.23462E+06, 0.24287E+06, 0.25135E+06, 0.26007E+06,
                                 0.26903E+06, 0.27824E+06, 0.28771E+06, 0.29743E+06, 0.30742E+06,
                                 0.31767E+06, 0.32820E+06, 0.33901E+06, 0.35011E+06, 0.36150E+06,
                                 0.37319E+06, 0.38518E+06, 0.39749E+06, 0.41010E+06, 0.42304E+06,
                                 0.43631E+06])


#  --------------- HCN 135: M = 23, I = 3 ---------------------
M = 23
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.11863E+03, 0.16755E+03, 0.21653E+03,
                                 0.26576E+03, 0.31559E+03, 0.36656E+03, 0.41926E+03, 0.47428E+03,
                                 0.53214E+03, 0.59333E+03, 0.65824E+03, 0.72727E+03, 0.80074E+03,
                                 0.87898E+03, 0.96227E+03, 0.10509E+04, 0.11452E+04, 0.12454E+04,
                                 0.13518E+04, 0.14647E+04, 0.15844E+04, 0.17112E+04, 0.18455E+04,
                                 0.19875E+04, 0.21377E+04, 0.22962E+04, 0.24636E+04, 0.26402E+04,
                                 0.28263E+04, 0.30224E+04, 0.32289E+04, 0.34461E+04, 0.36745E+04,
                                 0.39145E+04, 0.41667E+04, 0.44314E+04, 0.47092E+04, 0.50005E+04,
                                 0.53059E+04, 0.56259E+04, 0.59609E+04, 0.63116E+04, 0.66785E+04,
                                 0.70622E+04, 0.74633E+04, 0.78823E+04, 0.83200E+04, 0.87769E+04,
                                 0.92536E+04, 0.97509E+04, 0.10269E+05, 0.10810E+05, 0.11373E+05,
                                 0.11959E+05, 0.12570E+05, 0.13205E+05, 0.13866E+05, 0.14554E+05,
                                 0.15268E+05, 0.16011E+05, 0.16782E+05, 0.17583E+05, 0.18415E+05,
                                 0.19279E+05, 0.20174E+05, 0.21103E+05, 0.22067E+05, 0.23065E+05,
                                 0.24100E+05, 0.25172E+05, 0.26282E+05, 0.27432E+05, 0.28622E+05,
                                 0.29853E+05, 0.31127E+05, 0.32445E+05, 0.33807E+05, 0.35215E+05,
                                 0.36670E+05, 0.38174E+05, 0.39727E+05, 0.41330E+05, 0.42986E+05,
                                 0.44695E+05, 0.46459E+05, 0.48278E+05, 0.50155E+05, 0.52091E+05,
                                 0.54086E+05, 0.56143E+05, 0.58263E+05, 0.60447E+05, 0.62696E+05,
                                 0.65013E+05, 0.67399E+05, 0.69856E+05, 0.72384E+05, 0.74986E+05,
                                 0.77663E+05, 0.80416E+05, 0.83249E+05, 0.86161E+05, 0.89156E+05,
                                 0.92233E+05, 0.95397E+05, 0.98648E+05, 0.10199E+06, 0.10542E+06,
                                 0.10894E+06, 0.11256E+06, 0.11627E+06, 0.12009E+06, 0.12400E+06,
                                 0.12802E+06, 0.13214E+06, 0.13636E+06, 0.14070E+06, 0.14515E+06,
                                 0.14971E+06])


#  --------------- CH3Cl 215: M = 24, I = 1 ---------------------
M = 24
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.50529E+04, 0.85123E+04, 0.12528E+05,
                                 0.17036E+05, 0.22005E+05, 0.27429E+05, 0.33325E+05, 0.39734E+05,
                                 0.46713E+05, 0.54336E+05, 0.62690E+05, 0.71876E+05, 0.82006E+05,
                                 0.93204E+05, 0.10560E+06, 0.11936E+06, 0.13463E+06, 0.15158E+06,
                                 0.17043E+06, 0.19137E+06, 0.21464E+06, 0.24049E+06, 0.26920E+06,
                                 0.30107E+06, 0.33642E+06, 0.37563E+06, 0.41907E+06, 0.46719E+06,
                                 0.52045E+06, 0.57936E+06, 0.64448E+06, 0.71641E+06, 0.79582E+06,
                                 0.88341E+06, 0.97997E+06, 0.10863E+07, 0.12034E+07, 0.13323E+07,
                                 0.14739E+07, 0.16295E+07, 0.18003E+07, 0.19877E+07, 0.21932E+07,
                                 0.24183E+07, 0.26649E+07, 0.29346E+07, 0.32296E+07, 0.35519E+07,
                                 0.39039E+07, 0.42881E+07, 0.47072E+07, 0.51639E+07, 0.56615E+07,
                                 0.62032E+07, 0.67926E+07, 0.74335E+07, 0.81299E+07, 0.88862E+07,
                                 0.97071E+07, 0.10598E+08, 0.11563E+08, 0.12609E+08, 0.13742E+08,
                                 0.14968E+08, 0.16294E+08, 0.17728E+08, 0.19277E+08, 0.20950E+08,
                                 0.22756E+08, 0.24704E+08, 0.26805E+08, 0.29069E+08, 0.31507E+08,
                                 0.34132E+08, 0.36957E+08, 0.39995E+08, 0.43260E+08, 0.46769E+08,
                                 0.50538E+08, 0.54583E+08, 0.58923E+08, 0.63578E+08, 0.68568E+08,
                                 0.73914E+08, 0.79640E+08, 0.85770E+08, 0.92329E+08, 0.99345E+08,
                                 0.10685E+09, 0.11486E+09, 0.12342E+09, 0.13257E+09, 0.14233E+09,
                                 0.15274E+09, 0.16384E+09, 0.17568E+09, 0.18829E+09, 0.20173E+09,
                                 0.21604E+09, 0.23127E+09, 0.24748E+09, 0.26471E+09, 0.28304E+09,
                                 0.30252E+09, 0.32322E+09, 0.34520E+09, 0.36853E+09, 0.39330E+09,
                                 0.41958E+09, 0.44745E+09, 0.47701E+09, 0.50833E+09, 0.54151E+09,
                                 0.57667E+09, 0.61389E+09, 0.65329E+09, 0.69498E+09, 0.73909E+09,
                                 0.78573E+09])


#  --------------- CH3Cl 217: M = 24, I = 2 ---------------------
M = 24
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.51327E+04, 0.86469E+04, 0.12726E+05,
                                 0.17306E+05, 0.22354E+05, 0.27863E+05, 0.33853E+05, 0.40364E+05,
                                 0.47453E+05, 0.55197E+05, 0.63684E+05, 0.73016E+05, 0.83306E+05,
                                 0.94681E+05, 0.10728E+06, 0.12125E+06, 0.13676E+06, 0.15399E+06,
                                 0.17313E+06, 0.19441E+06, 0.21804E+06, 0.24430E+06, 0.27347E+06,
                                 0.30584E+06, 0.34176E+06, 0.38158E+06, 0.42572E+06, 0.47460E+06,
                                 0.52871E+06, 0.58855E+06, 0.65471E+06, 0.72778E+06, 0.80844E+06,
                                 0.89743E+06, 0.99552E+06, 0.11036E+07, 0.12225E+07, 0.13534E+07,
                                 0.14973E+07, 0.16553E+07, 0.18289E+07, 0.20193E+07, 0.22280E+07,
                                 0.24567E+07, 0.27072E+07, 0.29812E+07, 0.32808E+07, 0.36083E+07,
                                 0.39659E+07, 0.43562E+07, 0.47819E+07, 0.52459E+07, 0.57514E+07,
                                 0.63017E+07, 0.69005E+07, 0.75515E+07, 0.82590E+07, 0.90273E+07,
                                 0.98613E+07, 0.10766E+08, 0.11747E+08, 0.12809E+08, 0.13960E+08,
                                 0.15206E+08, 0.16553E+08, 0.18010E+08, 0.19584E+08, 0.21283E+08,
                                 0.23118E+08, 0.25097E+08, 0.27231E+08, 0.29531E+08, 0.32008E+08,
                                 0.34674E+08, 0.37544E+08, 0.40630E+08, 0.43948E+08, 0.47513E+08,
                                 0.51341E+08, 0.55451E+08, 0.59860E+08, 0.64589E+08, 0.69658E+08,
                                 0.75089E+08, 0.80906E+08, 0.87134E+08, 0.93797E+08, 0.10092E+09,
                                 0.10854E+09, 0.11669E+09, 0.12539E+09, 0.13467E+09, 0.14459E+09,
                                 0.15517E+09, 0.16645E+09, 0.17847E+09, 0.19129E+09, 0.20494E+09,
                                 0.21948E+09, 0.23495E+09, 0.25141E+09, 0.26893E+09, 0.28754E+09,
                                 0.30733E+09, 0.32836E+09, 0.35069E+09, 0.37440E+09, 0.39956E+09,
                                 0.42626E+09, 0.45457E+09, 0.48460E+09, 0.51642E+09, 0.55013E+09,
                                 0.58585E+09, 0.62366E+09, 0.66369E+09, 0.70605E+09, 0.75085E+09,
                                 0.79824E+09])


#  --------------- H2O2 1661: M = 25, I = 1 ---------------------
M = 25
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.62392E+03, 0.10958E+04, 0.16692E+04,
                                 0.23492E+04, 0.31427E+04, 0.40574E+04, 0.51014E+04, 0.62840E+04,
                                 0.76157E+04, 0.91085E+04, 0.10776E+05, 0.12633E+05, 0.14696E+05,
                                 0.16983E+05, 0.19515E+05, 0.22312E+05, 0.25396E+05, 0.28792E+05,
                                 0.32526E+05, 0.36625E+05, 0.41118E+05, 0.46036E+05, 0.51410E+05,
                                 0.57275E+05, 0.63667E+05, 0.70623E+05, 0.78185E+05, 0.86394E+05,
                                 0.95295E+05, 0.10493E+06, 0.11536E+06, 0.12662E+06, 0.13878E+06,
                                 0.15188E+06, 0.16600E+06, 0.18118E+06, 0.19750E+06, 0.21503E+06,
                                 0.23383E+06, 0.25398E+06, 0.27556E+06, 0.29864E+06, 0.32333E+06,
                                 0.34970E+06, 0.37784E+06, 0.40786E+06, 0.43985E+06, 0.47392E+06,
                                 0.51018E+06, 0.54874E+06, 0.58972E+06, 0.63324E+06, 0.67943E+06,
                                 0.72843E+06, 0.78037E+06, 0.83540E+06, 0.89366E+06, 0.95530E+06,
                                 0.10205E+07, 0.10894E+07, 0.11622E+07, 0.12391E+07, 0.13202E+07,
                                 0.14057E+07, 0.14959E+07, 0.15909E+07, 0.16910E+07, 0.17963E+07,
                                 0.19072E+07, 0.20237E+07, 0.21463E+07, 0.22750E+07, 0.24102E+07,
                                 0.25522E+07, 0.27012E+07, 0.28575E+07, 0.30213E+07, 0.31931E+07,
                                 0.33730E+07, 0.35615E+07, 0.37588E+07, 0.39653E+07, 0.41813E+07,
                                 0.44072E+07, 0.46433E+07, 0.48901E+07, 0.51479E+07, 0.54171E+07,
                                 0.56982E+07, 0.59915E+07, 0.62976E+07, 0.66167E+07, 0.69495E+07,
                                 0.72963E+07, 0.76577E+07, 0.80342E+07, 0.84262E+07, 0.88343E+07,
                                 0.92591E+07, 0.97011E+07, 0.10161E+08, 0.10639E+08, 0.11136E+08,
                                 0.11652E+08, 0.12189E+08, 0.12746E+08, 0.13325E+08, 0.13926E+08,
                                 0.14550E+08, 0.15198E+08, 0.15870E+08, 0.16566E+08, 0.17289E+08,
                                 0.18038E+08, 0.18814E+08, 0.19619E+08, 0.20452E+08, 0.21315E+08,
                                 0.22209E+08])


#  --------------- C2H2 1221: M = 26, I = 1 ---------------------
M = 26
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.71617E+02, 0.10121E+03, 0.13092E+03,
                                 0.16104E+03, 0.19218E+03, 0.22509E+03, 0.26062E+03, 0.29959E+03,
                                 0.34281E+03, 0.39103E+03, 0.44503E+03, 0.50558E+03, 0.57346E+03,
                                 0.64950E+03, 0.73457E+03, 0.82960E+03, 0.93557E+03, 0.10535E+04,
                                 0.11846E+04, 0.13301E+04, 0.14911E+04, 0.16692E+04, 0.18658E+04,
                                 0.20825E+04, 0.23211E+04, 0.25833E+04, 0.28711E+04, 0.31867E+04,
                                 0.35323E+04, 0.39102E+04, 0.43230E+04, 0.47735E+04, 0.52645E+04,
                                 0.57991E+04, 0.63807E+04, 0.70127E+04, 0.76988E+04, 0.84430E+04,
                                 0.92495E+04, 0.10123E+05, 0.11067E+05, 0.12088E+05, 0.13191E+05,
                                 0.14381E+05, 0.15664E+05, 0.17047E+05, 0.18536E+05, 0.20137E+05,
                                 0.21859E+05, 0.23710E+05, 0.25696E+05, 0.27827E+05, 0.30112E+05,
                                 0.32561E+05, 0.35183E+05, 0.37990E+05, 0.40991E+05, 0.44199E+05,
                                 0.47626E+05, 0.51285E+05, 0.55189E+05, 0.59353E+05, 0.63791E+05,
                                 0.68518E+05, 0.73551E+05, 0.78908E+05, 0.84604E+05, 0.90661E+05,
                                 0.97095E+05, 0.10393E+06, 0.11118E+06, 0.11888E+06, 0.12704E+06,
                                 0.13569E+06, 0.14486E+06, 0.15457E+06, 0.16485E+06, 0.17572E+06,
                                 0.18722E+06, 0.19938E+06, 0.21223E+06, 0.22581E+06, 0.24014E+06,
                                 0.25527E+06, 0.27123E+06, 0.28807E+06, 0.30582E+06, 0.32452E+06,
                                 0.34423E+06, 0.36498E+06, 0.38683E+06, 0.40982E+06, 0.43401E+06,
                                 0.45944E+06, 0.48618E+06, 0.51428E+06, 0.54380E+06, 0.57480E+06,
                                 0.60735E+06, 0.64151E+06, 0.67735E+06, 0.71494E+06, 0.75436E+06,
                                 0.79568E+06, 0.83898E+06, 0.88434E+06, 0.93184E+06, 0.98158E+06,
                                 0.10336E+07, 0.10881E+07, 0.11451E+07, 0.12047E+07, 0.12670E+07,
                                 0.13321E+07, 0.14002E+07, 0.14713E+07, 0.15455E+07, 0.16231E+07,
                                 0.17040E+07])


#  --------------- C2H2 1231: M = 26, I = 2 ---------------------
M = 26
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.)
TIPS_ISO_HASH[(M, I)] = float32([0.28647E+03, 0.40486E+03, 0.52369E+03,
                                 0.64419E+03, 0.76874E+03, 0.90040E+03, 0.10425E+04, 0.11984E+04,
                                 0.13713E+04, 0.15642E+04, 0.17802E+04, 0.20223E+04, 0.22939E+04,
                                 0.25981E+04, 0.29384E+04, 0.33185E+04, 0.37424E+04, 0.42142E+04,
                                 0.47386E+04, 0.53203E+04, 0.59646E+04, 0.66769E+04, 0.74633E+04,
                                 0.83302E+04, 0.92845E+04, 0.10333E+05, 0.11485E+05, 0.12747E+05,
                                 0.14129E+05, 0.15641E+05, 0.17292E+05, 0.19094E+05, 0.21058E+05,
                                 0.23197E+05, 0.25523E+05, 0.28051E+05, 0.30796E+05, 0.33773E+05,
                                 0.36999E+05, 0.40492E+05, 0.44270E+05, 0.48354E+05, 0.52765E+05,
                                 0.57525E+05, 0.62658E+05, 0.68189E+05, 0.74144E+05, 0.80551E+05,
                                 0.87439E+05, 0.94840E+05, 0.10279E+06, 0.11131E+06, 0.12045E+06,
                                 0.13025E+06, 0.14074E+06, 0.15196E+06, 0.16397E+06, 0.17680E+06,
                                 0.19051E+06, 0.20514E+06, 0.22076E+06, 0.23742E+06, 0.25517E+06,
                                 0.27408E+06, 0.29421E+06, 0.31564E+06, 0.33842E+06, 0.36265E+06,
                                 0.38839E+06, 0.41572E+06, 0.44474E+06, 0.47553E+06, 0.50818E+06,
                                 0.54278E+06, 0.57945E+06, 0.61829E+06, 0.65940E+06, 0.70289E+06,
                                 0.74890E+06, 0.79754E+06, 0.84894E+06, 0.90324E+06, 0.96057E+06,
                                 0.10211E+07, 0.10849E+07, 0.11523E+07, 0.12233E+07, 0.12981E+07,
                                 0.13769E+07, 0.14599E+07, 0.15473E+07, 0.16393E+07, 0.17361E+07,
                                 0.18378E+07, 0.19447E+07, 0.20571E+07, 0.21752E+07, 0.22992E+07,
                                 0.24294E+07, 0.25661E+07, 0.27094E+07, 0.28598E+07, 0.30175E+07,
                                 0.31828E+07, 0.33560E+07, 0.35374E+07, 0.37274E+07, 0.39264E+07,
                                 0.41346E+07, 0.43525E+07, 0.45805E+07, 0.48188E+07, 0.50681E+07,
                                 0.53286E+07, 0.56008E+07, 0.58852E+07, 0.61823E+07, 0.64924E+07,
                                 0.68162E+07])


#  --------------- C2H2 1222: M = 26, I = 3 ---------------------
M = 26
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.24843E+03, 0.35373E+03, 0.45997E+03,
                                 0.56930E+03, 0.68497E+03, 0.81065E+03, 0.94999E+03, 0.11065E+04,
                                 0.12837E+04, 0.14848E+04, 0.17135E+04, 0.19731E+04, 0.22675E+04,
                                 0.26205E+04, 0.29999E+04, 0.34276E+04, 0.39086E+04, 0.44486E+04,
                                 0.50533E+04, 0.57294E+04, 0.64837E+04, 0.73237E+04, 0.82576E+04,
                                 0.92941E+04, 0.10443E+05, 0.11714E+05, 0.13117E+05, 0.14666E+05,
                                 0.16373E+05, 0.18250E+05, 0.20313E+05, 0.22578E+05, 0.25060E+05,
                                 0.27777E+05, 0.30750E+05, 0.33997E+05, 0.37541E+05, 0.41405E+05,
                                 0.45614E+05, 0.50192E+05, 0.55170E+05, 0.60576E+05, 0.66441E+05,
                                 0.72799E+05, 0.79686E+05, 0.87140E+05, 0.95199E+05, 0.10391E+06,
                                 0.11331E+06, 0.12345E+06, 0.13438E+06, 0.14615E+06, 0.15882E+06,
                                 0.17245E+06, 0.18710E+06, 0.20283E+06, 0.21972E+06, 0.23783E+06,
                                 0.25724E+06, 0.27804E+06, 0.30030E+06, 0.32411E+06, 0.34958E+06,
                                 0.37679E+06, 0.40585E+06, 0.43686E+06, 0.46994E+06, 0.50521E+06,
                                 0.54280E+06, 0.58282E+06, 0.62542E+06, 0.67074E+06, 0.71892E+06,
                                 0.77013E+06, 0.82453E+06, 0.88228E+06, 0.94356E+06, 0.10086E+07,
                                 0.10775E+07, 0.11505E+07, 0.12279E+07, 0.13098E+07, 0.13964E+07,
                                 0.14881E+07, 0.15850E+07, 0.16875E+07, 0.17957E+07, 0.19100E+07,
                                 0.20307E+07, 0.21580E+07, 0.22923E+07, 0.24339E+07, 0.25831E+07,
                                 0.27404E+07, 0.29060E+07, 0.30803E+07, 0.32638E+07, 0.34568E+07,
                                 0.36598E+07, 0.38733E+07, 0.40976E+07, 0.43332E+07, 0.45807E+07,
                                 0.48406E+07, 0.51133E+07, 0.53995E+07, 0.56997E+07, 0.60144E+07,
                                 0.63444E+07, 0.66901E+07, 0.70524E+07, 0.74317E+07, 0.78289E+07,
                                 0.82447E+07, 0.86797E+07, 0.91348E+07, 0.96108E+07, 0.10108E+08,
                                 0.10629E+08])


#  --------------- C2H6 1221: M = 27, I = 1 ---------------------
M = 27
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.47267E+04, 0.80011E+04, 0.11928E+05,
                                 0.16564E+05, 0.21985E+05, 0.28287E+05, 0.35590E+05, 0.44049E+05,
                                 0.53862E+05, 0.65277E+05, 0.78597E+05, 0.94191E+05, 0.11250E+06,
                                 0.13407E+06, 0.15952E+06, 0.18962E+06, 0.22526E+06, 0.26751E+06,
                                 0.31763E+06, 0.37714E+06, 0.44780E+06, 0.53174E+06, 0.63145E+06,
                                 0.74989E+06, 0.89056E+06, 0.10576E+07, 0.12559E+07, 0.14912E+07,
                                 0.17704E+07, 0.21013E+07, 0.24936E+07, 0.29582E+07, 0.35083E+07,
                                 0.41591E+07, 0.49286E+07, 0.58379E+07, 0.69116E+07, 0.81787E+07,
                                 0.96728E+07, 0.11433E+08, 0.13506E+08, 0.15945E+08, 0.18812E+08,
                                 0.22180E+08, 0.26134E+08, 0.30770E+08, 0.36204E+08, 0.42565E+08,
                                 0.50008E+08, 0.58708E+08, 0.68868E+08, 0.80725E+08, 0.94548E+08,
                                 0.11065E+09, 0.12940E+09, 0.15119E+09, 0.17652E+09, 0.20593E+09,
                                 0.24003E+09, 0.27956E+09, 0.32533E+09, 0.37829E+09, 0.43951E+09,
                                 0.51021E+09, 0.59180E+09, 0.68588E+09, 0.79427E+09, 0.91904E+09,
                                 0.10625E+10, 0.12275E+10, 0.14168E+10, 0.16341E+10, 0.18831E+10,
                                 0.21684E+10, 0.24949E+10, 0.28684E+10, 0.32951E+10, 0.37823E+10,
                                 0.43382E+10, 0.49719E+10, 0.56938E+10, 0.65156E+10, 0.74502E+10,
                                 0.85125E+10, 0.97190E+10, 0.11088E+11, 0.12641E+11, 0.14401E+11,
                                 0.16393E+11, 0.18648E+11, 0.21198E+11, 0.24079E+11, 0.27332E+11,
                                 0.31003E+11, 0.35142E+11, 0.39807E+11, 0.45060E+11, 0.50972E+11,
                                 0.57620E+11, 0.65091E+11, 0.73483E+11, 0.82902E+11, 0.93467E+11,
                                 0.10531E+12, 0.11858E+12, 0.13343E+12, 0.15005E+12, 0.16864E+12,
                                 0.18941E+12, 0.21260E+12, 0.23849E+12, 0.26737E+12, 0.29957E+12,
                                 0.33545E+12, 0.37541E+12, 0.41987E+12, 0.46934E+12, 0.52432E+12,
                                 0.58542E+12])


#  --------------- C2H6 1231: M = 27, I = 2 ---------------------
M = 27
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.24128E+04, 0.40845E+04, 0.60896E+04,
                                 0.84564E+04, 0.11224E+05, 0.14442E+05, 0.18170E+05, 0.22490E+05,
                                 0.27501E+05, 0.33329E+05, 0.40131E+05, 0.48094E+05, 0.57446E+05,
                                 0.68459E+05, 0.81458E+05, 0.96828E+05, 0.11503E+06, 0.13661E+06,
                                 0.16221E+06, 0.19260E+06, 0.22869E+06, 0.27156E+06, 0.32249E+06,
                                 0.38298E+06, 0.45483E+06, 0.54015E+06, 0.64144E+06, 0.76164E+06,
                                 0.90423E+06, 0.10733E+07, 0.12737E+07, 0.15110E+07, 0.17920E+07,
                                 0.21245E+07, 0.25176E+07, 0.29821E+07, 0.35307E+07, 0.41780E+07,
                                 0.49414E+07, 0.58408E+07, 0.68999E+07, 0.81461E+07, 0.96110E+07,
                                 0.11332E+08, 0.13352E+08, 0.15721E+08, 0.18497E+08, 0.21748E+08,
                                 0.25551E+08, 0.29997E+08, 0.35189E+08, 0.41248E+08, 0.48313E+08,
                                 0.56542E+08, 0.66122E+08, 0.77262E+08, 0.90206E+08, 0.10523E+09,
                                 0.12267E+09, 0.14287E+09, 0.16626E+09, 0.19333E+09, 0.22462E+09,
                                 0.26076E+09, 0.30247E+09, 0.35056E+09, 0.40596E+09, 0.46974E+09,
                                 0.54310E+09, 0.62740E+09, 0.72420E+09, 0.83527E+09, 0.96260E+09,
                                 0.11084E+10, 0.12754E+10, 0.14663E+10, 0.16845E+10, 0.19336E+10,
                                 0.22178E+10, 0.25418E+10, 0.29109E+10, 0.33311E+10, 0.38090E+10,
                                 0.43522E+10, 0.49691E+10, 0.56693E+10, 0.64633E+10, 0.73631E+10,
                                 0.83821E+10, 0.95352E+10, 0.10839E+11, 0.12312E+11, 0.13976E+11,
                                 0.15854E+11, 0.17971E+11, 0.20357E+11, 0.23043E+11, 0.26067E+11,
                                 0.29467E+11, 0.33289E+11, 0.37581E+11, 0.42399E+11, 0.47804E+11,
                                 0.53862E+11, 0.60649E+11, 0.68247E+11, 0.76750E+11, 0.86257E+11,
                                 0.96882E+11, 0.10875E+12, 0.12199E+12, 0.13677E+12, 0.15325E+12,
                                 0.17160E+12, 0.19204E+12, 0.21480E+12, 0.24010E+12, 0.26824E+12,
                                 0.29950E+12])


#  --------------- PH3 1111: M = 28, I = 1 ---------------------
M = 28
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.29652E+03, 0.49643E+03, 0.72810E+03,
                                 0.98777E+03, 0.12729E+04, 0.15820E+04, 0.19145E+04, 0.22708E+04,
                                 0.26520E+04, 0.30600E+04, 0.34971E+04, 0.39662E+04, 0.44702E+04,
                                 0.50126E+04, 0.55970E+04, 0.62273E+04, 0.69075E+04, 0.76421E+04,
                                 0.84357E+04, 0.92933E+04, 0.10220E+05, 0.11222E+05, 0.12304E+05,
                                 0.13473E+05, 0.14736E+05, 0.16099E+05, 0.17571E+05, 0.19160E+05,
                                 0.20873E+05, 0.22720E+05, 0.24710E+05, 0.26854E+05, 0.29162E+05,
                                 0.31646E+05, 0.34317E+05, 0.37188E+05, 0.40273E+05, 0.43585E+05,
                                 0.47140E+05, 0.50953E+05, 0.55040E+05, 0.59419E+05, 0.64108E+05,
                                 0.69127E+05, 0.74496E+05, 0.80236E+05, 0.86369E+05, 0.92918E+05,
                                 0.99909E+05, 0.10737E+06, 0.11532E+06, 0.12380E+06, 0.13282E+06,
                                 0.14244E+06, 0.15266E+06, 0.16354E+06, 0.17511E+06, 0.18739E+06,
                                 0.20044E+06, 0.21430E+06, 0.22900E+06, 0.24459E+06, 0.26111E+06,
                                 0.27862E+06, 0.29716E+06, 0.31680E+06, 0.33757E+06, 0.35954E+06,
                                 0.38277E+06, 0.40733E+06, 0.43326E+06, 0.46065E+06, 0.48955E+06,
                                 0.52005E+06, 0.55222E+06, 0.58614E+06, 0.62188E+06, 0.65953E+06,
                                 0.69917E+06, 0.74091E+06, 0.78483E+06, 0.83103E+06, 0.87960E+06,
                                 0.93067E+06, 0.98432E+06, 0.10407E+07, 0.10999E+07, 0.11620E+07,
                                 0.12272E+07, 0.12956E+07, 0.13673E+07, 0.14425E+07, 0.15212E+07,
                                 0.16038E+07, 0.16902E+07, 0.17808E+07, 0.18755E+07, 0.19746E+07,
                                 0.20784E+07, 0.21868E+07, 0.23002E+07, 0.24187E+07, 0.25425E+07,
                                 0.26719E+07, 0.28070E+07, 0.29480E+07, 0.30952E+07, 0.32488E+07,
                                 0.34091E+07, 0.35762E+07, 0.37504E+07, 0.39320E+07, 0.41213E+07,
                                 0.43185E+07, 0.45239E+07, 0.47378E+07, 0.49605E+07, 0.51923E+07,
                                 0.54335E+07])


#  --------------- COF2 269: M = 29, I = 1 ---------------------
M = 29
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.54999E+04, 0.92749E+04, 0.13668E+05,
                                 0.18643E+05, 0.24224E+05, 0.30487E+05, 0.37547E+05, 0.45543E+05,
                                 0.54639E+05, 0.65019E+05, 0.76886E+05, 0.90462E+05, 0.10600E+06,
                                 0.12377E+06, 0.14407E+06, 0.16723E+06, 0.19363E+06, 0.22367E+06,
                                 0.25780E+06, 0.29650E+06, 0.34031E+06, 0.38982E+06, 0.44568E+06,
                                 0.50859E+06, 0.57932E+06, 0.65872E+06, 0.74770E+06, 0.84724E+06,
                                 0.95844E+06, 0.10825E+07, 0.12205E+07, 0.13741E+07, 0.15446E+07,
                                 0.17336E+07, 0.19428E+07, 0.21742E+07, 0.24296E+07, 0.27113E+07,
                                 0.30214E+07, 0.33626E+07, 0.37373E+07, 0.41484E+07, 0.45989E+07,
                                 0.50921E+07, 0.56313E+07, 0.62202E+07, 0.68626E+07, 0.75628E+07,
                                 0.83251E+07, 0.91542E+07, 0.10055E+08, 0.11033E+08, 0.12093E+08,
                                 0.13242E+08, 0.14486E+08, 0.15831E+08, 0.17284E+08, 0.18853E+08,
                                 0.20546E+08, 0.22371E+08, 0.24335E+08, 0.26450E+08, 0.28724E+08,
                                 0.31167E+08, 0.33790E+08, 0.36605E+08, 0.39623E+08, 0.42856E+08,
                                 0.46318E+08, 0.50022E+08, 0.53983E+08, 0.58215E+08, 0.62735E+08,
                                 0.67558E+08, 0.72702E+08, 0.78186E+08, 0.84028E+08, 0.90247E+08,
                                 0.96865E+08, 0.10390E+09, 0.11138E+09, 0.11933E+09, 0.12777E+09,
                                 0.13672E+09, 0.14622E+09, 0.15629E+09, 0.16695E+09, 0.17825E+09,
                                 0.19021E+09, 0.20287E+09, 0.21625E+09, 0.23039E+09, 0.24534E+09,
                                 0.26113E+09, 0.27779E+09, 0.29538E+09, 0.31392E+09, 0.33348E+09,
                                 0.35409E+09, 0.37580E+09, 0.39867E+09, 0.42274E+09, 0.44806E+09,
                                 0.47470E+09, 0.50271E+09, 0.53215E+09, 0.56308E+09, 0.59557E+09,
                                 0.62968E+09, 0.66548E+09, 0.70304E+09, 0.74243E+09, 0.78374E+09,
                                 0.82703E+09, 0.87240E+09, 0.91992E+09, 0.96967E+09, 0.10218E+10,
                                 0.10763E+10])


#  --------------- COF2 369: M = 29, I = 2 --------------------- not in TIPS-2011
M = 29
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- SF6 29: M = 30, I = 1 ---------------------
M = 30
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.46373E+05, 0.78844E+05, 0.11939E+06,
                                 0.17183E+06, 0.24247E+06, 0.34059E+06, 0.47963E+06, 0.67906E+06,
                                 0.96713E+06, 0.13848E+07, 0.19911E+07, 0.28714E+07, 0.41481E+07,
                                 0.59956E+07, 0.86617E+07, 0.12496E+08, 0.17991E+08, 0.25832E+08,
                                 0.36971E+08, 0.52724E+08, 0.74895E+08, 0.10595E+09, 0.14923E+09,
                                 0.20925E+09, 0.29208E+09, 0.40582E+09, 0.56124E+09, 0.77259E+09,
                                 0.10586E+10, 0.14439E+10, 0.19605E+10, 0.26500E+10, 0.35662E+10,
                                 0.47781E+10, 0.63747E+10, 0.84689E+10, 0.11205E+11, 0.14765E+11,
                                 0.19378E+11, 0.25336E+11, 0.32998E+11, 0.42819E+11, 0.55361E+11,
                                 0.71323E+11, 0.91569E+11, 0.11716E+12, 0.14941E+12, 0.18992E+12,
                                 0.24065E+12, 0.30398E+12, 0.38283E+12, 0.48069E+12, 0.60182E+12,
                                 0.75136E+12, 0.93546E+12, 0.11615E+13, 0.14384E+13, 0.17767E+13,
                                 0.21890E+13, 0.26903E+13, 0.32984E+13, 0.40344E+13, 0.49232E+13,
                                 0.59942E+13, 0.72819E+13, 0.88272E+13, 0.10678E+14, 0.12889E+14,
                                 0.15527E+14, 0.18666E+14, 0.22397E+14, 0.26823E+14, 0.32062E+14,
                                 0.38253E+14, 0.45558E+14, 0.54161E+14, 0.64277E+14, 0.76153E+14,
                                 0.90072E+14, 0.10636E+15, 0.12539E+15, 0.14759E+15, 0.17345E+15,
                                 0.20354E+15, 0.23848E+15, 0.27902E+15, 0.32597E+15, 0.38028E+15,
                                 0.44303E+15, 0.51542E+15, 0.59883E+15, 0.69482E+15, 0.80516E+15,
                                 0.93182E+15, 0.10770E+16, 0.12434E+16, 0.14336E+16, 0.16511E+16,
                                 0.18992E+16, 0.21821E+16, 0.25043E+16, 0.28709E+16, 0.32875E+16,
                                 0.37604E+16, 0.42968E+16, 0.49046E+16, 0.55925E+16, 0.63704E+16,
                                 0.72492E+16, 0.82411E+16, 0.93596E+16, 0.10620E+17, 0.12038E+17,
                                 0.13633E+17, 0.15425E+17, 0.17438E+17, 0.19694E+17, 0.22224E+17,
                                 0.25057E+17])


#  --------------- H2S 121: M = 31, I = 1 ---------------------
M = 31
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.47192E+02, 0.78671E+02, 0.11510E+03,
                                 0.15589E+03, 0.20061E+03, 0.24896E+03, 0.30070E+03, 0.35571E+03,
                                 0.41386E+03, 0.47513E+03, 0.53951E+03, 0.60703E+03, 0.67772E+03,
                                 0.75167E+03, 0.82896E+03, 0.90969E+03, 0.99396E+03, 0.10819E+04,
                                 0.11736E+04, 0.12692E+04, 0.13689E+04, 0.14727E+04, 0.15809E+04,
                                 0.16937E+04, 0.18111E+04, 0.19333E+04, 0.20606E+04, 0.21931E+04,
                                 0.23309E+04, 0.24744E+04, 0.26236E+04, 0.27788E+04, 0.29403E+04,
                                 0.31081E+04, 0.32825E+04, 0.34638E+04, 0.36522E+04, 0.38478E+04,
                                 0.40510E+04, 0.42619E+04, 0.44808E+04, 0.47080E+04, 0.49437E+04,
                                 0.51881E+04, 0.54415E+04, 0.57042E+04, 0.59764E+04, 0.62584E+04,
                                 0.65505E+04, 0.68529E+04, 0.71660E+04, 0.74899E+04, 0.78251E+04,
                                 0.81718E+04, 0.85303E+04, 0.89008E+04, 0.92838E+04, 0.96795E+04,
                                 0.10088E+05, 0.10510E+05, 0.10946E+05, 0.11396E+05, 0.11860E+05,
                                 0.12339E+05, 0.12833E+05, 0.13342E+05, 0.13867E+05, 0.14408E+05,
                                 0.14966E+05, 0.15540E+05, 0.16132E+05, 0.16741E+05, 0.17368E+05,
                                 0.18013E+05, 0.18677E+05, 0.19361E+05, 0.20064E+05, 0.20786E+05,
                                 0.21529E+05, 0.22293E+05, 0.23078E+05, 0.23885E+05, 0.24714E+05,
                                 0.25565E+05, 0.26439E+05, 0.27337E+05, 0.28258E+05, 0.29204E+05,
                                 0.30174E+05, 0.31170E+05, 0.32191E+05, 0.33239E+05, 0.34313E+05,
                                 0.35414E+05, 0.36543E+05, 0.37700E+05, 0.38886E+05, 0.40101E+05,
                                 0.41346E+05, 0.42621E+05, 0.43926E+05, 0.45263E+05, 0.46631E+05,
                                 0.48033E+05, 0.49466E+05, 0.50934E+05, 0.52435E+05, 0.53971E+05,
                                 0.55542E+05, 0.57149E+05, 0.58792E+05, 0.60472E+05, 0.62190E+05,
                                 0.63946E+05, 0.65740E+05, 0.67574E+05, 0.69448E+05, 0.71362E+05,
                                 0.73318E+05])


#  --------------- H2S 141: M = 31, I = 2 ---------------------
M = 31
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.47310E+02, 0.78869E+02, 0.11539E+03,
                                 0.15628E+03, 0.20112E+03, 0.24959E+03, 0.30147E+03, 0.35661E+03,
                                 0.41491E+03, 0.47634E+03, 0.54088E+03, 0.60857E+03, 0.67945E+03,
                                 0.75359E+03, 0.83107E+03, 0.91201E+03, 0.99649E+03, 0.10846E+04,
                                 0.11766E+04, 0.12724E+04, 0.13724E+04, 0.14765E+04, 0.15850E+04,
                                 0.16980E+04, 0.18157E+04, 0.19382E+04, 0.20658E+04, 0.21987E+04,
                                 0.23369E+04, 0.24807E+04, 0.26303E+04, 0.27859E+04, 0.29478E+04,
                                 0.31160E+04, 0.32909E+04, 0.34727E+04, 0.36615E+04, 0.38576E+04,
                                 0.40613E+04, 0.42728E+04, 0.44923E+04, 0.47200E+04, 0.49563E+04,
                                 0.52013E+04, 0.54554E+04, 0.57188E+04, 0.59917E+04, 0.62744E+04,
                                 0.65672E+04, 0.68704E+04, 0.71843E+04, 0.75090E+04, 0.78451E+04,
                                 0.81926E+04, 0.85520E+04, 0.89236E+04, 0.93075E+04, 0.97042E+04,
                                 0.10114E+05, 0.10537E+05, 0.10974E+05, 0.11425E+05, 0.11890E+05,
                                 0.12370E+05, 0.12866E+05, 0.13376E+05, 0.13903E+05, 0.14445E+05,
                                 0.15004E+05, 0.15580E+05, 0.16173E+05, 0.16784E+05, 0.17412E+05,
                                 0.18059E+05, 0.18725E+05, 0.19410E+05, 0.20115E+05, 0.20839E+05,
                                 0.21584E+05, 0.22350E+05, 0.23137E+05, 0.23946E+05, 0.24777E+05,
                                 0.25630E+05, 0.26507E+05, 0.27407E+05, 0.28330E+05, 0.29278E+05,
                                 0.30251E+05, 0.31249E+05, 0.32273E+05, 0.33324E+05, 0.34401E+05,
                                 0.35505E+05, 0.36637E+05, 0.37797E+05, 0.38985E+05, 0.40204E+05,
                                 0.41451E+05, 0.42729E+05, 0.44038E+05, 0.45379E+05, 0.46751E+05,
                                 0.48155E+05, 0.49593E+05, 0.51064E+05, 0.52569E+05, 0.54109E+05,
                                 0.55684E+05, 0.57295E+05, 0.58943E+05, 0.60627E+05, 0.62349E+05,
                                 0.64109E+05, 0.65908E+05, 0.67747E+05, 0.69625E+05, 0.71544E+05,
                                 0.73505E+05])


#  --------------- H2S 131: M = 30, I = 3 ---------------------
M = 31
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.18901E+03, 0.31509E+03, 0.46102E+03,
                                 0.62437E+03, 0.80349E+03, 0.99713E+03, 0.12044E+04, 0.14247E+04,
                                 0.16576E+04, 0.19030E+04, 0.21609E+04, 0.24313E+04, 0.27145E+04,
                                 0.30106E+04, 0.33202E+04, 0.36436E+04, 0.39811E+04, 0.43332E+04,
                                 0.47005E+04, 0.50835E+04, 0.54827E+04, 0.58987E+04, 0.63321E+04,
                                 0.67836E+04, 0.72538E+04, 0.77434E+04, 0.82532E+04, 0.87838E+04,
                                 0.93360E+04, 0.99106E+04, 0.10508E+05, 0.11130E+05, 0.11777E+05,
                                 0.12449E+05, 0.13147E+05, 0.13874E+05, 0.14628E+05, 0.15412E+05,
                                 0.16225E+05, 0.17070E+05, 0.17947E+05, 0.18857E+05, 0.19801E+05,
                                 0.20780E+05, 0.21795E+05, 0.22847E+05, 0.23937E+05, 0.25067E+05,
                                 0.26236E+05, 0.27448E+05, 0.28702E+05, 0.29999E+05, 0.31342E+05,
                                 0.32730E+05, 0.34166E+05, 0.35650E+05, 0.37184E+05, 0.38769E+05,
                                 0.40406E+05, 0.42097E+05, 0.43842E+05, 0.45644E+05, 0.47503E+05,
                                 0.49421E+05, 0.51399E+05, 0.53439E+05, 0.55542E+05, 0.57709E+05,
                                 0.59942E+05, 0.62242E+05, 0.64611E+05, 0.67051E+05, 0.69563E+05,
                                 0.72148E+05, 0.74808E+05, 0.77545E+05, 0.80360E+05, 0.83255E+05,
                                 0.86232E+05, 0.89291E+05, 0.92435E+05, 0.95667E+05, 0.98986E+05,
                                 0.10240E+06, 0.10590E+06, 0.10949E+06, 0.11318E+06, 0.11697E+06,
                                 0.12086E+06, 0.12484E+06, 0.12893E+06, 0.13313E+06, 0.13743E+06,
                                 0.14184E+06, 0.14637E+06, 0.15100E+06, 0.15575E+06, 0.16062E+06,
                                 0.16560E+06, 0.17071E+06, 0.17594E+06, 0.18129E+06, 0.18677E+06,
                                 0.19238E+06, 0.19813E+06, 0.20400E+06, 0.21002E+06, 0.21617E+06,
                                 0.22246E+06, 0.22890E+06, 0.23548E+06, 0.24221E+06, 0.24909E+06,
                                 0.25612E+06, 0.26331E+06, 0.27065E+06, 0.27816E+06, 0.28583E+06,
                                 0.29366E+06])


#  --------------- HCOOH 126: M = 32, I = 1 ---------------------
M = 32
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.31899E+04, 0.53773E+04, 0.79205E+04,
                                 0.10792E+05, 0.13993E+05, 0.17550E+05, 0.21509E+05, 0.25930E+05,
                                 0.30885E+05, 0.36460E+05, 0.42750E+05, 0.49864E+05, 0.57926E+05,
                                 0.67071E+05, 0.77453E+05, 0.89243E+05, 0.10263E+06, 0.11783E+06,
                                 0.13507E+06, 0.15462E+06, 0.17676E+06, 0.20183E+06, 0.23018E+06,
                                 0.26221E+06, 0.29836E+06, 0.33911E+06, 0.38501E+06, 0.43664E+06,
                                 0.49467E+06, 0.55981E+06, 0.63286E+06, 0.71470E+06, 0.80628E+06,
                                 0.90865E+06, 0.10230E+07, 0.11505E+07, 0.12927E+07, 0.14509E+07,
                                 0.16269E+07, 0.18225E+07, 0.20396E+07, 0.22804E+07, 0.25472E+07,
                                 0.28425E+07, 0.31692E+07, 0.35301E+07, 0.39285E+07, 0.43681E+07,
                                 0.48525E+07, 0.53858E+07, 0.59727E+07, 0.66178E+07, 0.73265E+07,
                                 0.81042E+07, 0.89571E+07, 0.98918E+07, 0.10915E+08, 0.12035E+08,
                                 0.13259E+08, 0.14597E+08, 0.16057E+08, 0.17650E+08, 0.19387E+08,
                                 0.21279E+08, 0.23339E+08, 0.25579E+08, 0.28016E+08, 0.30663E+08,
                                 0.33536E+08, 0.36655E+08, 0.40037E+08, 0.43701E+08, 0.47671E+08,
                                 0.51967E+08, 0.56614E+08, 0.61639E+08, 0.67068E+08, 0.72930E+08,
                                 0.79257E+08, 0.86082E+08, 0.93439E+08, 0.10137E+09, 0.10990E+09,
                                 0.11909E+09, 0.12898E+09, 0.13960E+09, 0.15102E+09, 0.16329E+09,
                                 0.17646E+09, 0.19059E+09, 0.20575E+09, 0.22200E+09, 0.23941E+09,
                                 0.25806E+09, 0.27802E+09, 0.29938E+09, 0.32223E+09, 0.34666E+09,
                                 0.37276E+09, 0.40064E+09, 0.43041E+09, 0.46218E+09, 0.49607E+09,
                                 0.53221E+09, 0.57074E+09, 0.61179E+09, 0.65551E+09, 0.70206E+09,
                                 0.75159E+09, 0.80430E+09, 0.86034E+09, 0.91992E+09, 0.98324E+09,
                                 0.10505E+10, 0.11219E+10, 0.11977E+10, 0.12782E+10, 0.13635E+10,
                                 0.14540E+10])


#  --------------- HO2 166: M = 33, I = 1 ---------------------
M = 33
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.39277E+03, 0.66062E+03, 0.97123E+03,
                                 0.13194E+04, 0.17014E+04, 0.21148E+04, 0.25578E+04, 0.30296E+04,
                                 0.35297E+04, 0.40585E+04, 0.46167E+04, 0.52055E+04, 0.58264E+04,
                                 0.64809E+04, 0.71707E+04, 0.78978E+04, 0.86641E+04, 0.94715E+04,
                                 0.10322E+05, 0.11218E+05, 0.12161E+05, 0.13154E+05, 0.14198E+05,
                                 0.15296E+05, 0.16449E+05, 0.17661E+05, 0.18933E+05, 0.20267E+05,
                                 0.21666E+05, 0.23133E+05, 0.24669E+05, 0.26277E+05, 0.27960E+05,
                                 0.29720E+05, 0.31560E+05, 0.33482E+05, 0.35489E+05, 0.37584E+05,
                                 0.39769E+05, 0.42048E+05, 0.44423E+05, 0.46898E+05, 0.49475E+05,
                                 0.52157E+05, 0.54948E+05, 0.57850E+05, 0.60868E+05, 0.64003E+05,
                                 0.67261E+05, 0.70643E+05, 0.74154E+05, 0.77797E+05, 0.81575E+05,
                                 0.85492E+05, 0.89553E+05, 0.93760E+05, 0.98118E+05, 0.10263E+06,
                                 0.10730E+06, 0.11213E+06, 0.11713E+06, 0.12230E+06, 0.12765E+06,
                                 0.13317E+06, 0.13888E+06, 0.14478E+06, 0.15086E+06, 0.15715E+06,
                                 0.16363E+06, 0.17032E+06, 0.17723E+06, 0.18434E+06, 0.19168E+06,
                                 0.19924E+06, 0.20704E+06, 0.21506E+06, 0.22333E+06, 0.23185E+06,
                                 0.24061E+06, 0.24963E+06, 0.25891E+06, 0.26846E+06, 0.27828E+06,
                                 0.28838E+06, 0.29876E+06, 0.30943E+06, 0.32039E+06, 0.33166E+06,
                                 0.34323E+06, 0.35512E+06, 0.36732E+06, 0.37985E+06, 0.39271E+06,
                                 0.40590E+06, 0.41944E+06, 0.43333E+06, 0.44758E+06, 0.46219E+06,
                                 0.47717E+06, 0.49252E+06, 0.50826E+06, 0.52439E+06, 0.54091E+06,
                                 0.55784E+06, 0.57518E+06, 0.59293E+06, 0.61112E+06, 0.62973E+06,
                                 0.64878E+06, 0.66828E+06, 0.68824E+06, 0.70866E+06, 0.72955E+06,
                                 0.75091E+06, 0.77276E+06, 0.79511E+06, 0.81795E+06, 0.84131E+06,
                                 0.86518E+06])


#  --------------- O 6: M = 34, I = 1 --------------------- not in TIPS-2011
M = 34
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- ClONO2 5646: M = 35, I = 1 ---------------------
M = 35
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.11444E+06, 0.21121E+06, 0.34858E+06,
                                 0.53934E+06, 0.80041E+06, 0.11539E+07, 0.16286E+07, 0.22614E+07,
                                 0.30992E+07, 0.42015E+07, 0.56426E+07, 0.75152E+07, 0.99344E+07,
                                 0.13042E+08, 0.17012E+08, 0.22058E+08, 0.28437E+08, 0.36463E+08,
                                 0.46514E+08, 0.59042E+08, 0.74589E+08, 0.93801E+08, 0.11744E+09,
                                 0.14643E+09, 0.18181E+09, 0.22486E+09, 0.27705E+09, 0.34009E+09,
                                 0.41598E+09, 0.50705E+09, 0.61599E+09, 0.74590E+09, 0.90037E+09,
                                 0.10835E+10, 0.13001E+10, 0.15554E+10, 0.18556E+10, 0.22079E+10,
                                 0.26200E+10, 0.31012E+10, 0.36615E+10, 0.43126E+10, 0.50675E+10,
                                 0.59409E+10, 0.69492E+10, 0.81110E+10, 0.94469E+10, 0.10980E+11,
                                 0.12736E+11, 0.14745E+11, 0.17037E+11, 0.19649E+11, 0.22620E+11,
                                 0.25994E+11, 0.29819E+11, 0.34150E+11, 0.39044E+11, 0.44568E+11,
                                 0.50794E+11, 0.57799E+11, 0.65672E+11, 0.74506E+11, 0.84408E+11,
                                 0.95490E+11, 0.10788E+12, 0.12171E+12, 0.13713E+12, 0.15431E+12,
                                 0.17342E+12, 0.19465E+12, 0.21822E+12, 0.24435E+12, 0.27329E+12,
                                 0.30530E+12, 0.34069E+12, 0.37976E+12, 0.42286E+12, 0.47034E+12,
                                 0.52262E+12, 0.58012E+12, 0.64330E+12, 0.71267E+12, 0.78875E+12,
                                 0.87214E+12, 0.96344E+12, 0.10633E+13, 0.11725E+13, 0.12918E+13,
                                 0.14220E+13, 0.15640E+13, 0.17188E+13, 0.18873E+13, 0.20706E+13,
                                 0.22700E+13, 0.24866E+13, 0.27218E+13, 0.29771E+13, 0.32538E+13,
                                 0.35537E+13, 0.38784E+13, 0.42299E+13, 0.46100E+13, 0.50208E+13,
                                 0.54645E+13, 0.59435E+13, 0.64603E+13, 0.70175E+13, 0.76180E+13,
                                 0.82647E+13, 0.89608E+13, 0.97097E+13, 0.10515E+14, 0.11380E+14,
                                 0.12310E+14, 0.13307E+14, 0.14378E+14, 0.15526E+14, 0.16756E+14,
                                 0.18075E+14])


#  --------------- ClONO2 7646: M = 35, I = 2 ---------------------
M = 35
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.11735E+06, 0.21659E+06, 0.35745E+06,
                                 0.55307E+06, 0.82078E+06, 0.11833E+07, 0.16700E+07, 0.23189E+07,
                                 0.31781E+07, 0.43084E+07, 0.57862E+07, 0.77065E+07, 0.10187E+08,
                                 0.13374E+08, 0.17445E+08, 0.22619E+08, 0.29161E+08, 0.37391E+08,
                                 0.47698E+08, 0.60545E+08, 0.76487E+08, 0.96188E+08, 0.12043E+09,
                                 0.15015E+09, 0.18644E+09, 0.23059E+09, 0.28410E+09, 0.34874E+09,
                                 0.42657E+09, 0.51995E+09, 0.63167E+09, 0.76489E+09, 0.92329E+09,
                                 0.11111E+10, 0.13331E+10, 0.15950E+10, 0.19029E+10, 0.22641E+10,
                                 0.26867E+10, 0.31801E+10, 0.37547E+10, 0.44224E+10, 0.51965E+10,
                                 0.60921E+10, 0.71261E+10, 0.83174E+10, 0.96873E+10, 0.11260E+11,
                                 0.13061E+11, 0.15120E+11, 0.17471E+11, 0.20149E+11, 0.23196E+11,
                                 0.26656E+11, 0.30578E+11, 0.35019E+11, 0.40038E+11, 0.45703E+11,
                                 0.52087E+11, 0.59270E+11, 0.67343E+11, 0.76403E+11, 0.86556E+11,
                                 0.97921E+11, 0.11062E+12, 0.12481E+12, 0.14062E+12, 0.15824E+12,
                                 0.17783E+12, 0.19961E+12, 0.22377E+12, 0.25057E+12, 0.28024E+12,
                                 0.31308E+12, 0.34936E+12, 0.38943E+12, 0.43362E+12, 0.48232E+12,
                                 0.53593E+12, 0.59489E+12, 0.65968E+12, 0.73081E+12, 0.80883E+12,
                                 0.89434E+12, 0.98797E+12, 0.10904E+13, 0.12024E+13, 0.13247E+13,
                                 0.14582E+13, 0.16038E+13, 0.17625E+13, 0.19353E+13, 0.21233E+13,
                                 0.23278E+13, 0.25499E+13, 0.27911E+13, 0.30528E+13, 0.33366E+13,
                                 0.36442E+13, 0.39772E+13, 0.43376E+13, 0.47273E+13, 0.51486E+13,
                                 0.56036E+13, 0.60948E+13, 0.66248E+13, 0.71962E+13, 0.78119E+13,
                                 0.84751E+13, 0.91889E+13, 0.99569E+13, 0.10783E+14, 0.11670E+14,
                                 0.12623E+14, 0.13646E+14, 0.14744E+14, 0.15921E+14, 0.17183E+14,
                                 0.18535E+14])


#  --------------- NOp 46: M = 36, I = 1 ---------------------
M = 36
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.)
TIPS_ISO_HASH[(M, I)] = float32([0.63956E+02, 0.90185E+02, 0.11642E+03,
                                 0.14265E+03, 0.16889E+03, 0.19513E+03, 0.22138E+03, 0.24763E+03,
                                 0.27388E+03, 0.30013E+03, 0.32639E+03, 0.35266E+03, 0.37894E+03,
                                 0.40523E+03, 0.43155E+03, 0.45790E+03, 0.48429E+03, 0.51074E+03,
                                 0.53725E+03, 0.56383E+03, 0.59052E+03, 0.61731E+03, 0.64422E+03,
                                 0.67127E+03, 0.69846E+03, 0.72582E+03, 0.75335E+03, 0.78108E+03,
                                 0.80901E+03, 0.83715E+03, 0.86552E+03, 0.89413E+03, 0.92298E+03,
                                 0.95208E+03, 0.98144E+03, 0.10111E+04, 0.10410E+04, 0.10712E+04,
                                 0.11017E+04, 0.11325E+04, 0.11636E+04, 0.11950E+04, 0.12268E+04,
                                 0.12588E+04, 0.12912E+04, 0.13239E+04, 0.13570E+04, 0.13903E+04,
                                 0.14241E+04, 0.14581E+04, 0.14926E+04, 0.15273E+04, 0.15624E+04,
                                 0.15979E+04, 0.16337E+04, 0.16699E+04, 0.17065E+04, 0.17434E+04,
                                 0.17806E+04, 0.18183E+04, 0.18563E+04, 0.18947E+04, 0.19334E+04,
                                 0.19725E+04, 0.20120E+04, 0.20519E+04, 0.20921E+04, 0.21327E+04,
                                 0.21737E+04, 0.22151E+04, 0.22568E+04, 0.22990E+04, 0.23415E+04,
                                 0.23844E+04, 0.24276E+04, 0.24713E+04, 0.25153E+04, 0.25598E+04,
                                 0.26046E+04, 0.26497E+04, 0.26953E+04, 0.27413E+04, 0.27876E+04,
                                 0.28343E+04, 0.28815E+04, 0.29290E+04, 0.29769E+04, 0.30251E+04,
                                 0.30738E+04, 0.31229E+04, 0.31723E+04, 0.32222E+04, 0.32724E+04,
                                 0.33230E+04, 0.33740E+04, 0.34254E+04, 0.34772E+04, 0.35294E+04,
                                 0.35819E+04, 0.36349E+04, 0.36883E+04, 0.37420E+04, 0.37961E+04,
                                 0.38507E+04, 0.39056E+04, 0.39609E+04, 0.40166E+04, 0.40727E+04,
                                 0.41292E+04, 0.41861E+04, 0.42434E+04, 0.43010E+04, 0.43591E+04,
                                 0.44176E+04, 0.44764E+04, 0.45357E+04, 0.45953E+04, 0.46554E+04,
                                 0.47158E+04])


#  --------------- HOBr 169: M = 37, I = 1 ---------------------
M = 37
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.)
TIPS_ISO_HASH[(M, I)] = float32([0.24445E+04, 0.41206E+04, 0.60683E+04,
                                 0.82610E+04, 0.10689E+05, 0.13352E+05, 0.16261E+05, 0.19427E+05,
                                 0.22867E+05, 0.26600E+05, 0.30643E+05, 0.35018E+05, 0.39745E+05,
                                 0.44844E+05, 0.50338E+05, 0.56249E+05, 0.62599E+05, 0.69410E+05,
                                 0.76706E+05, 0.84509E+05, 0.92845E+05, 0.10174E+06, 0.11121E+06,
                                 0.12128E+06, 0.13199E+06, 0.14335E+06, 0.15540E+06, 0.16815E+06,
                                 0.18165E+06, 0.19591E+06, 0.21096E+06, 0.22684E+06, 0.24358E+06,
                                 0.26120E+06, 0.27974E+06, 0.29922E+06, 0.31969E+06, 0.34118E+06,
                                 0.36372E+06, 0.38735E+06, 0.41210E+06, 0.43800E+06, 0.46511E+06,
                                 0.49345E+06, 0.52307E+06, 0.55400E+06, 0.58628E+06, 0.61997E+06,
                                 0.65509E+06, 0.69170E+06, 0.72984E+06, 0.76954E+06, 0.81087E+06,
                                 0.85386E+06, 0.89856E+06, 0.94502E+06, 0.99329E+06, 0.10434E+07,
                                 0.10955E+07, 0.11495E+07, 0.12055E+07, 0.12636E+07, 0.13238E+07,
                                 0.13862E+07, 0.14508E+07, 0.15177E+07, 0.15870E+07, 0.16587E+07,
                                 0.17328E+07, 0.18095E+07, 0.18888E+07, 0.19707E+07, 0.20554E+07,
                                 0.21428E+07, 0.22331E+07, 0.23263E+07, 0.24225E+07, 0.25217E+07,
                                 0.26241E+07, 0.27296E+07, 0.28385E+07, 0.29506E+07, 0.30662E+07,
                                 0.31853E+07, 0.33079E+07, 0.34341E+07, 0.35641E+07, 0.36979E+07,
                                 0.38355E+07, 0.39771E+07, 0.41228E+07, 0.42725E+07, 0.44265E+07,
                                 0.45848E+07, 0.47474E+07, 0.49145E+07, 0.50862E+07, 0.52624E+07,
                                 0.54435E+07, 0.56293E+07, 0.58201E+07, 0.60159E+07, 0.62168E+07,
                                 0.64229E+07, 0.66343E+07, 0.68511E+07, 0.70734E+07, 0.73013E+07,
                                 0.75349E+07, 0.77742E+07, 0.80196E+07, 0.82709E+07, 0.85283E+07,
                                 0.87920E+07, 0.90620E+07, 0.93385E+07, 0.96215E+07, 0.99112E+07,
                                 0.10208E+08])


#  --------------- HOBr 161: M = 37, I = 2 ---------------------
M = 37
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(8.)
TIPS_ISO_HASH[(M, I)] = float32([0.24350E+04, 0.41047E+04, 0.60448E+04,
                                 0.82291E+04, 0.10648E+05, 0.13301E+05, 0.16200E+05, 0.19355E+05,
                                 0.22784E+05, 0.26504E+05, 0.30534E+05, 0.34895E+05, 0.39607E+05,
                                 0.44691E+05, 0.50169E+05, 0.56063E+05, 0.62394E+05, 0.69186E+05,
                                 0.76461E+05, 0.84243E+05, 0.92555E+05, 0.10142E+06, 0.11087E+06,
                                 0.12091E+06, 0.13159E+06, 0.14292E+06, 0.15494E+06, 0.16766E+06,
                                 0.18112E+06, 0.19534E+06, 0.21036E+06, 0.22620E+06, 0.24289E+06,
                                 0.26047E+06, 0.27896E+06, 0.29840E+06, 0.31882E+06, 0.34025E+06,
                                 0.36274E+06, 0.38630E+06, 0.41099E+06, 0.43683E+06, 0.46387E+06,
                                 0.49215E+06, 0.52169E+06, 0.55255E+06, 0.58475E+06, 0.61836E+06,
                                 0.65340E+06, 0.68992E+06, 0.72796E+06, 0.76757E+06, 0.80880E+06,
                                 0.85169E+06, 0.89628E+06, 0.94263E+06, 0.99079E+06, 0.10408E+07,
                                 0.10927E+07, 0.11466E+07, 0.12025E+07, 0.12605E+07, 0.13205E+07,
                                 0.13828E+07, 0.14472E+07, 0.15140E+07, 0.15831E+07, 0.16546E+07,
                                 0.17286E+07, 0.18051E+07, 0.18842E+07, 0.19660E+07, 0.20504E+07,
                                 0.21377E+07, 0.22277E+07, 0.23207E+07, 0.24167E+07, 0.25157E+07,
                                 0.26178E+07, 0.27231E+07, 0.28317E+07, 0.29436E+07, 0.30589E+07,
                                 0.31777E+07, 0.33001E+07, 0.34260E+07, 0.35557E+07, 0.36892E+07,
                                 0.38265E+07, 0.39678E+07, 0.41131E+07, 0.42626E+07, 0.44162E+07,
                                 0.45741E+07, 0.47364E+07, 0.49031E+07, 0.50744E+07, 0.52503E+07,
                                 0.54309E+07, 0.56164E+07, 0.58067E+07, 0.60021E+07, 0.62025E+07,
                                 0.64081E+07, 0.66191E+07, 0.68354E+07, 0.70572E+07, 0.72846E+07,
                                 0.75177E+07, 0.77565E+07, 0.80013E+07, 0.82521E+07, 0.85090E+07,
                                 0.87721E+07, 0.90415E+07, 0.93173E+07, 0.95997E+07, 0.98888E+07,
                                 0.10185E+08])


#  --------------- C2H4 221: M = 38, I = 1 ---------------------
M = 38
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.95843E+03, 0.16137E+04, 0.23744E+04,
                                 0.32285E+04, 0.41694E+04, 0.51963E+04, 0.63143E+04, 0.75337E+04,
                                 0.88702E+04, 0.10344E+05, 0.11978E+05, 0.13802E+05, 0.15846E+05,
                                 0.18145E+05, 0.20740E+05, 0.23675E+05, 0.27000E+05, 0.30770E+05,
                                 0.35048E+05, 0.39905E+05, 0.45420E+05, 0.51680E+05, 0.58786E+05,
                                 0.66850E+05, 0.75997E+05, 0.86369E+05, 0.98123E+05, 0.11144E+06,
                                 0.12651E+06, 0.14356E+06, 0.16284E+06, 0.18463E+06, 0.20923E+06,
                                 0.23699E+06, 0.26831E+06, 0.30360E+06, 0.34334E+06, 0.38808E+06,
                                 0.43840E+06, 0.49495E+06, 0.55847E+06, 0.62976E+06, 0.70973E+06,
                                 0.79935E+06, 0.89973E+06, 0.10121E+07, 0.11378E+07, 0.12782E+07,
                                 0.14351E+07, 0.16102E+07, 0.18055E+07, 0.20231E+07, 0.22656E+07,
                                 0.25354E+07, 0.28356E+07, 0.31692E+07, 0.35398E+07, 0.39511E+07,
                                 0.44074E+07, 0.49132E+07, 0.54736E+07, 0.60940E+07, 0.67803E+07,
                                 0.75392E+07, 0.83776E+07, 0.93035E+07, 0.10325E+08, 0.11452E+08,
                                 0.12694E+08, 0.14062E+08, 0.15567E+08, 0.17224E+08, 0.19045E+08,
                                 0.21046E+08, 0.23243E+08, 0.25655E+08, 0.28300E+08, 0.31200E+08,
                                 0.34377E+08, 0.37856E+08, 0.41662E+08, 0.45826E+08, 0.50378E+08,
                                 0.55351E+08, 0.60781E+08, 0.66707E+08, 0.73172E+08, 0.80219E+08,
                                 0.87899E+08, 0.96262E+08, 0.10537E+09, 0.11527E+09, 0.12604E+09,
                                 0.13775E+09, 0.15047E+09, 0.16428E+09, 0.17927E+09, 0.19553E+09,
                                 0.21316E+09, 0.23226E+09, 0.25296E+09, 0.27537E+09, 0.29963E+09,
                                 0.32587E+09, 0.35425E+09, 0.38492E+09, 0.41805E+09, 0.45383E+09,
                                 0.49246E+09, 0.53413E+09, 0.57908E+09, 0.62754E+09, 0.67977E+09,
                                 0.73602E+09, 0.79660E+09, 0.86179E+09, 0.93194E+09, 0.10074E+10,
                                 0.10885E+10])


#  --------------- C2H4 231: M = 38, I = 2 ---------------------
M = 38
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.39228E+04, 0.66051E+04, 0.97190E+04,
                                 0.13215E+05, 0.17066E+05, 0.21270E+05, 0.25846E+05, 0.30838E+05,
                                 0.36309E+05, 0.42341E+05, 0.49032E+05, 0.56496E+05, 0.64862E+05,
                                 0.74275E+05, 0.84897E+05, 0.96912E+05, 0.11052E+06, 0.12595E+06,
                                 0.14347E+06, 0.16335E+06, 0.18592E+06, 0.21155E+06, 0.24064E+06,
                                 0.27365E+06, 0.31109E+06, 0.35354E+06, 0.40166E+06, 0.45615E+06,
                                 0.51785E+06, 0.58765E+06, 0.66657E+06, 0.75575E+06, 0.85646E+06,
                                 0.97011E+06, 0.10983E+07, 0.12428E+07, 0.14055E+07, 0.15886E+07,
                                 0.17945E+07, 0.20260E+07, 0.22861E+07, 0.25779E+07, 0.29052E+07,
                                 0.32721E+07, 0.36830E+07, 0.41429E+07, 0.46573E+07, 0.52323E+07,
                                 0.58744E+07, 0.65912E+07, 0.73906E+07, 0.82816E+07, 0.92740E+07,
                                 0.10379E+08, 0.11607E+08, 0.12973E+08, 0.14490E+08, 0.16174E+08,
                                 0.18042E+08, 0.20112E+08, 0.22406E+08, 0.24945E+08, 0.27755E+08,
                                 0.30861E+08, 0.34293E+08, 0.38083E+08, 0.42266E+08, 0.46878E+08,
                                 0.51961E+08, 0.57560E+08, 0.63724E+08, 0.70504E+08, 0.77959E+08,
                                 0.86150E+08, 0.95145E+08, 0.10502E+09, 0.11585E+09, 0.12772E+09,
                                 0.14072E+09, 0.15496E+09, 0.17054E+09, 0.18759E+09, 0.20622E+09,
                                 0.22658E+09, 0.24880E+09, 0.27306E+09, 0.29952E+09, 0.32837E+09,
                                 0.35981E+09, 0.39404E+09, 0.43131E+09, 0.47186E+09, 0.51595E+09,
                                 0.56387E+09, 0.61594E+09, 0.67247E+09, 0.73382E+09, 0.80038E+09,
                                 0.87255E+09, 0.95076E+09, 0.10355E+10, 0.11272E+10, 0.12265E+10,
                                 0.13339E+10, 0.14501E+10, 0.15756E+10, 0.17113E+10, 0.18577E+10,
                                 0.20159E+10, 0.21865E+10, 0.23705E+10, 0.25688E+10, 0.27826E+10,
                                 0.30129E+10, 0.32608E+10, 0.35277E+10, 0.38149E+10, 0.41237E+10,
                                 0.44557E+10])


#  --------------- CH3OH 2161: M = 39, I = 1 --------------------- not in TIPS-2011
M = 39
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


#  --------------- CH3Br 219: M = 40, I = 1 ---------------------
M = 40
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.70299E+04, 0.11847E+05, 0.17442E+05,
                                 0.23741E+05, 0.30723E+05, 0.38408E+05, 0.46851E+05, 0.56138E+05,
                                 0.66375E+05, 0.77692E+05, 0.90239E+05, 0.10418E+06, 0.11972E+06,
                                 0.13704E+06, 0.15639E+06, 0.17801E+06, 0.20218E+06, 0.22920E+06,
                                 0.25940E+06, 0.29316E+06, 0.33087E+06, 0.37296E+06, 0.41992E+06,
                                 0.47229E+06, 0.53062E+06, 0.59557E+06, 0.66781E+06, 0.74812E+06,
                                 0.83731E+06, 0.93629E+06, 0.10461E+07, 0.11677E+07, 0.13023E+07,
                                 0.14513E+07, 0.16159E+07, 0.17978E+07, 0.19985E+07, 0.22199E+07,
                                 0.24638E+07, 0.27324E+07, 0.30280E+07, 0.33529E+07, 0.37099E+07,
                                 0.41019E+07, 0.45319E+07, 0.50034E+07, 0.55199E+07, 0.60853E+07,
                                 0.67039E+07, 0.73801E+07, 0.81189E+07, 0.89255E+07, 0.98056E+07,
                                 0.10765E+08, 0.11811E+08, 0.12949E+08, 0.14188E+08, 0.15535E+08,
                                 0.17000E+08, 0.18590E+08, 0.20317E+08, 0.22190E+08, 0.24220E+08,
                                 0.26421E+08, 0.28804E+08, 0.31383E+08, 0.34173E+08, 0.37189E+08,
                                 0.40448E+08, 0.43967E+08, 0.47765E+08, 0.51862E+08, 0.56280E+08,
                                 0.61040E+08, 0.66167E+08, 0.71686E+08, 0.77624E+08, 0.84009E+08,
                                 0.90873E+08, 0.98247E+08, 0.10616E+09, 0.11466E+09, 0.12378E+09,
                                 0.13356E+09, 0.14403E+09, 0.15526E+09, 0.16728E+09, 0.18014E+09,
                                 0.19391E+09, 0.20863E+09, 0.22436E+09, 0.24117E+09, 0.25913E+09,
                                 0.27830E+09, 0.29875E+09, 0.32057E+09, 0.34384E+09, 0.36864E+09,
                                 0.39506E+09, 0.42320E+09, 0.45316E+09, 0.48504E+09, 0.51896E+09,
                                 0.55502E+09, 0.59336E+09, 0.63410E+09, 0.67738E+09, 0.72334E+09,
                                 0.77212E+09, 0.82388E+09, 0.87879E+09, 0.93701E+09, 0.99873E+09,
                                 0.10641E+10, 0.11334E+10, 0.12068E+10, 0.12845E+10, 0.13667E+10,
                                 0.14536E+10])


#  --------------- CH3Br 211: M = 40, I = 2 ---------------------
M = 40
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.70566E+04, 0.11892E+05, 0.17508E+05,
                                 0.23832E+05, 0.30841E+05, 0.38557E+05, 0.47036E+05, 0.56362E+05,
                                 0.66644E+05, 0.78011E+05, 0.90615E+05, 0.10462E+06, 0.12023E+06,
                                 0.13763E+06, 0.15707E+06, 0.17880E+06, 0.20308E+06, 0.23023E+06,
                                 0.26059E+06, 0.29451E+06, 0.33240E+06, 0.37471E+06, 0.42191E+06,
                                 0.47453E+06, 0.53316E+06, 0.59843E+06, 0.67104E+06, 0.75176E+06,
                                 0.84141E+06, 0.94090E+06, 0.10512E+07, 0.11735E+07, 0.13088E+07,
                                 0.14585E+07, 0.16241E+07, 0.18069E+07, 0.20086E+07, 0.22312E+07,
                                 0.24764E+07, 0.27464E+07, 0.30435E+07, 0.33702E+07, 0.37291E+07,
                                 0.41231E+07, 0.45554E+07, 0.50294E+07, 0.55486E+07, 0.61171E+07,
                                 0.67389E+07, 0.74188E+07, 0.81616E+07, 0.89725E+07, 0.98573E+07,
                                 0.10822E+08, 0.11873E+08, 0.13018E+08, 0.14263E+08, 0.15618E+08,
                                 0.17090E+08, 0.18689E+08, 0.20425E+08, 0.22308E+08, 0.24350E+08,
                                 0.26563E+08, 0.28959E+08, 0.31552E+08, 0.34357E+08, 0.37389E+08,
                                 0.40666E+08, 0.44204E+08, 0.48023E+08, 0.52143E+08, 0.56585E+08,
                                 0.61371E+08, 0.66526E+08, 0.72076E+08, 0.78046E+08, 0.84467E+08,
                                 0.91369E+08, 0.98783E+08, 0.10674E+09, 0.11529E+09, 0.12446E+09,
                                 0.13429E+09, 0.14482E+09, 0.15611E+09, 0.16820E+09, 0.18113E+09,
                                 0.19497E+09, 0.20978E+09, 0.22560E+09, 0.24250E+09, 0.26056E+09,
                                 0.27983E+09, 0.30040E+09, 0.32234E+09, 0.34574E+09, 0.37068E+09,
                                 0.39725E+09, 0.42555E+09, 0.45567E+09, 0.48773E+09, 0.52184E+09,
                                 0.55811E+09, 0.59666E+09, 0.63763E+09, 0.68115E+09, 0.72736E+09,
                                 0.77642E+09, 0.82847E+09, 0.88368E+09, 0.94223E+09, 0.10043E+10,
                                 0.10701E+10, 0.11397E+10, 0.12135E+10, 0.12916E+10, 0.13743E+10,
                                 0.14618E+10])


#  --------------- CH3CN 2124: M = 41, I = 1 ---------------------
M = 41
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(3.)
TIPS_ISO_HASH[(M, I)] = float32([0.54361E+04, 0.91953E+04, 0.13708E+05,
                                 0.19097E+05, 0.25531E+05, 0.33206E+05, 0.42337E+05, 0.53173E+05,
                                 0.66002E+05, 0.81163E+05, 0.99053E+05, 0.12014E+06, 0.14496E+06,
                                 0.17414E+06, 0.20843E+06, 0.24866E+06, 0.29580E+06, 0.35099E+06,
                                 0.41551E+06, 0.49085E+06, 0.57871E+06, 0.68104E+06, 0.80008E+06,
                                 0.93836E+06, 0.10988E+07, 0.12848E+07, 0.14999E+07, 0.17487E+07,
                                 0.20359E+07, 0.23670E+07, 0.27484E+07, 0.31871E+07, 0.36912E+07,
                                 0.42697E+07, 0.49328E+07, 0.56921E+07, 0.65605E+07, 0.75526E+07,
                                 0.86847E+07, 0.99753E+07, 0.11445E+08, 0.13116E+08, 0.15016E+08,
                                 0.17172E+08, 0.19617E+08, 0.22386E+08, 0.25520E+08, 0.29063E+08,
                                 0.33064E+08, 0.37578E+08, 0.42667E+08, 0.48397E+08, 0.54844E+08,
                                 0.62090E+08, 0.70228E+08, 0.79358E+08, 0.89592E+08, 0.10105E+09,
                                 0.11388E+09, 0.12822E+09, 0.14424E+09, 0.16212E+09, 0.18205E+09,
                                 0.20427E+09, 0.22900E+09, 0.25652E+09, 0.28710E+09, 0.32107E+09,
                                 0.35877E+09, 0.40059E+09, 0.44692E+09, 0.49822E+09, 0.55500E+09,
                                 0.61777E+09, 0.68712E+09, 0.76370E+09, 0.84819E+09, 0.94135E+09,
                                 0.10440E+10, 0.11570E+10, 0.12814E+10, 0.14181E+10, 0.15684E+10,
                                 0.17334E+10, 0.19145E+10, 0.21131E+10, 0.23308E+10, 0.25693E+10,
                                 0.28304E+10, 0.31161E+10, 0.34285E+10, 0.37698E+10, 0.41426E+10,
                                 0.45496E+10, 0.49935E+10, 0.54776E+10, 0.60051E+10, 0.65796E+10,
                                 0.72049E+10, 0.78853E+10, 0.86251E+10, 0.94291E+10, 0.10303E+11,
                                 0.11251E+11, 0.12280E+11, 0.13396E+11, 0.14606E+11, 0.15916E+11,
                                 0.17336E+11, 0.18873E+11, 0.20536E+11, 0.22334E+11, 0.24278E+11,
                                 0.26379E+11, 0.28647E+11, 0.31096E+11, 0.33739E+11, 0.36589E+11,
                                 0.39661E+11])


#  --------------- CH3CN 2134: M = 41, I = 2 --------------------- not in HITRAN-2012
M = 41
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.10906E+05, 0.18458E+05, 0.27552E+05,
                                 0.38455E+05, 0.51523E+05, 0.67161E+05, 0.85818E+05, 0.10801E+06,
                                 0.13434E+06, 0.16550E+06, 0.20234E+06, 0.24581E+06, 0.29705E+06,
                                 0.35737E+06, 0.42831E+06, 0.51162E+06, 0.60936E+06, 0.72387E+06,
                                 0.85786E+06, 0.10145E+07, 0.11972E+07, 0.14102E+07, 0.16582E+07,
                                 0.19465E+07, 0.22813E+07, 0.26695E+07, 0.31190E+07, 0.36390E+07,
                                 0.42397E+07, 0.49328E+07, 0.57314E+07, 0.66507E+07, 0.77076E+07,
                                 0.89211E+07, 0.10313E+08, 0.11907E+08, 0.13732E+08, 0.15817E+08,
                                 0.18198E+08, 0.20914E+08, 0.24007E+08, 0.27527E+08, 0.31529E+08,
                                 0.36073E+08, 0.41228E+08, 0.47070E+08, 0.53683E+08, 0.61162E+08,
                                 0.69612E+08, 0.79149E+08, 0.89903E+08, 0.10202E+09, 0.11565E+09,
                                 0.13098E+09, 0.14820E+09, 0.16753E+09, 0.18921E+09, 0.21349E+09,
                                 0.24066E+09, 0.27106E+09, 0.30502E+09, 0.34293E+09, 0.38523E+09,
                                 0.43237E+09, 0.48486E+09, 0.54328E+09, 0.60823E+09, 0.68039E+09,
                                 0.76049E+09, 0.84935E+09, 0.94784E+09, 0.10569E+10, 0.11777E+10,
                                 0.13112E+10, 0.14588E+10, 0.16217E+10, 0.18016E+10, 0.19999E+10,
                                 0.22185E+10, 0.24592E+10, 0.27241E+10, 0.30155E+10, 0.33357E+10,
                                 0.36875E+10, 0.40736E+10, 0.44971E+10, 0.49615E+10, 0.54702E+10,
                                 0.60273E+10, 0.66369E+10, 0.73035E+10, 0.80322E+10, 0.88282E+10,
                                 0.96972E+10, 0.10645E+11, 0.11679E+11, 0.12806E+11, 0.14034E+11,
                                 0.15370E+11, 0.16824E+11, 0.18406E+11, 0.20125E+11, 0.21992E+11,
                                 0.24020E+11, 0.26221E+11, 0.28608E+11, 0.31197E+11, 0.34002E+11,
                                 0.37040E+11, 0.40330E+11, 0.43889E+11, 0.47739E+11, 0.51902E+11,
                                 0.56400E+11, 0.61259E+11, 0.66504E+11, 0.72165E+11, 0.78272E+11,
                                 0.84856E+11])


#  --------------- CH3CN 3124: M = 41, I = 3 --------------------- not in HITRAN-2012
M = 41
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.11223E+05, 0.18985E+05, 0.28307E+05,
                                 0.39441E+05, 0.52744E+05, 0.68620E+05, 0.87523E+05, 0.10997E+06,
                                 0.13658E+06, 0.16806E+06, 0.20524E+06, 0.24910E+06, 0.30080E+06,
                                 0.36165E+06, 0.43319E+06, 0.51722E+06, 0.61579E+06, 0.73127E+06,
                                 0.86640E+06, 0.10243E+07, 0.12086E+07, 0.14234E+07, 0.16735E+07,
                                 0.19642E+07, 0.23017E+07, 0.26931E+07, 0.31464E+07, 0.36706E+07,
                                 0.42762E+07, 0.49749E+07, 0.57801E+07, 0.67069E+07, 0.77722E+07,
                                 0.89955E+07, 0.10398E+08, 0.12006E+08, 0.13845E+08, 0.15947E+08,
                                 0.18346E+08, 0.21083E+08, 0.24201E+08, 0.27748E+08, 0.31781E+08,
                                 0.36361E+08, 0.41556E+08, 0.47442E+08, 0.54106E+08, 0.61643E+08,
                                 0.70157E+08, 0.79767E+08, 0.90604E+08, 0.10281E+09, 0.11655E+09,
                                 0.13199E+09, 0.14935E+09, 0.16882E+09, 0.19065E+09, 0.21512E+09,
                                 0.24250E+09, 0.27312E+09, 0.30733E+09, 0.34553E+09, 0.38814E+09,
                                 0.43562E+09, 0.48851E+09, 0.54736E+09, 0.61279E+09, 0.68548E+09,
                                 0.76617E+09, 0.85568E+09, 0.95489E+09, 0.10648E+10, 0.11864E+10,
                                 0.13209E+10, 0.14695E+10, 0.16337E+10, 0.18148E+10, 0.20146E+10,
                                 0.22348E+10, 0.24772E+10, 0.27441E+10, 0.30375E+10, 0.33601E+10,
                                 0.37143E+10, 0.41032E+10, 0.45298E+10, 0.49975E+10, 0.55099E+10,
                                 0.60709E+10, 0.66849E+10, 0.73563E+10, 0.80902E+10, 0.88918E+10,
                                 0.97670E+10, 0.10722E+11, 0.11763E+11, 0.12898E+11, 0.14134E+11,
                                 0.15480E+11, 0.16945E+11, 0.18537E+11, 0.20269E+11, 0.22149E+11,
                                 0.24191E+11, 0.26408E+11, 0.28812E+11, 0.31419E+11, 0.34244E+11,
                                 0.37303E+11, 0.40616E+11, 0.44201E+11, 0.48078E+11, 0.52269E+11,
                                 0.56799E+11, 0.61692E+11, 0.66974E+11, 0.72675E+11, 0.78824E+11,
                                 0.85454E+11])


#  --------------- CH3CN 3134: M = 41, I = 4 --------------------- not in HITRAN-2012
M = 41
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.22522E+05, 0.38117E+05, 0.56899E+05,
                                 0.79412E+05, 0.10640E+06, 0.13870E+06, 0.17726E+06, 0.22314E+06,
                                 0.27761E+06, 0.34214E+06, 0.41847E+06, 0.50862E+06, 0.61497E+06,
                                 0.74028E+06, 0.88774E+06, 0.10611E+07, 0.12646E+07, 0.15031E+07,
                                 0.17825E+07, 0.21092E+07, 0.24908E+07, 0.29358E+07, 0.34541E+07,
                                 0.40571E+07, 0.47576E+07, 0.55703E+07, 0.65120E+07, 0.76018E+07,
                                 0.88614E+07, 0.10315E+08, 0.11992E+08, 0.13922E+08, 0.16142E+08,
                                 0.18693E+08, 0.21619E+08, 0.24973E+08, 0.28812E+08, 0.33202E+08,
                                 0.38216E+08, 0.43936E+08, 0.50455E+08, 0.57876E+08, 0.66315E+08,
                                 0.75901E+08, 0.86779E+08, 0.99110E+08, 0.11307E+09, 0.12887E+09,
                                 0.14672E+09, 0.16688E+09, 0.18961E+09, 0.21523E+09, 0.24407E+09,
                                 0.27651E+09, 0.31295E+09, 0.35387E+09, 0.39975E+09, 0.45118E+09,
                                 0.50875E+09, 0.57315E+09, 0.64512E+09, 0.72549E+09, 0.81517E+09,
                                 0.91514E+09, 0.10265E+10, 0.11504E+10, 0.12883E+10, 0.14414E+10,
                                 0.16115E+10, 0.18001E+10, 0.20093E+10, 0.22410E+10, 0.24975E+10,
                                 0.27812E+10, 0.30948E+10, 0.34412E+10, 0.38235E+10, 0.42452E+10,
                                 0.47101E+10, 0.52220E+10, 0.57856E+10, 0.64055E+10, 0.70869E+10,
                                 0.78355E+10, 0.86574E+10, 0.95591E+10, 0.10548E+11, 0.11631E+11,
                                 0.12817E+11, 0.14116E+11, 0.15536E+11, 0.17088E+11, 0.18785E+11,
                                 0.20636E+11, 0.22657E+11, 0.24861E+11, 0.27264E+11, 0.29881E+11,
                                 0.32730E+11, 0.35832E+11, 0.39205E+11, 0.42871E+11, 0.46855E+11,
                                 0.51182E+11, 0.55878E+11, 0.60973E+11, 0.66497E+11, 0.72484E+11,
                                 0.78970E+11, 0.85992E+11, 0.93592E+11, 0.10181E+12, 0.11070E+12,
                                 0.12031E+12, 0.13069E+12, 0.14189E+12, 0.15398E+12, 0.16703E+12,
                                 0.18110E+12])


#  --------------- CF4 29: M = 42, I = 1 ---------------------
M = 42
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.76233E+04, 0.12867E+05, 0.19059E+05,
                                 0.26316E+05, 0.34895E+05, 0.45145E+05, 0.57461E+05, 0.72259E+05,
                                 0.89950E+05, 0.11092E+06, 0.13550E+06, 0.16399E+06, 0.19658E+06,
                                 0.23341E+06, 0.27457E+06, 0.32004E+06, 0.36978E+06, 0.42369E+06,
                                 0.48161E+06, 0.54338E+06, 0.60880E+06, 0.67764E+06, 0.55684E+07,
                                 0.71250E+07, 0.90615E+07, 0.11458E+08, 0.14407E+08, 0.18021E+08,
                                 0.22428E+08, 0.27778E+08, 0.34247E+08, 0.42038E+08, 0.51386E+08,
                                 0.62559E+08, 0.75869E+08, 0.91670E+08, 0.11037E+09, 0.13242E+09,
                                 0.15836E+09, 0.18878E+09, 0.22436E+09, 0.26584E+09, 0.31410E+09,
                                 0.37008E+09, 0.43488E+09, 0.50970E+09, 0.59589E+09, 0.69496E+09,
                                 0.80858E+09, 0.93863E+09, 0.10872E+10, 0.12565E+10, 0.14491E+10,
                                 0.16679E+10, 0.19159E+10, 0.21966E+10, 0.25136E+10, 0.28711E+10,
                                 0.32740E+10, 0.37260E+10, 0.42340E+10, 0.48030E+10, 0.54400E+10,
                                 0.61520E+10, 0.69470E+10, 0.78320E+10, 0.88170E+10, 0.99120E+10,
                                 0.11130E+11, 0.12470E+11, 0.13970E+11, 0.15620E+11, 0.17440E+11,
                                 0.19450E+11, 0.21670E+11, 0.24100E+11, 0.26790E+11, 0.29730E+11,
                                 0.33000E+11, 0.36500E+11, 0.40400E+11, 0.44600E+11, 0.49300E+11,
                                 0.54300E+11, 0.59800E+11, 0.65800E+11, 0.72400E+11, 0.79500E+11,
                                 0.87200E+11, 0.95500E+11, 0.10500E+12, 0.11400E+12, 0.12500E+12,
                                 0.13600E+12, 0.14900E+12, 0.16200E+12, 0.17700E+12, 0.19200E+12,
                                 0.21000E+12, 0.23000E+12, 0.25000E+12, 0.27000E+12, 0.29000E+12,
                                 0.31000E+12, 0.34000E+12, 0.36000E+12, 0.39000E+12, 0.42000E+12,
                                 0.46000E+12, 0.49000E+12, 0.53000E+12, 0.57000E+12, 0.61000E+12,
                                 0.66000E+12, 0.70000E+12, 0.75000E+12, 0.81000E+12, 0.86000E+12,
                                 0.93000E+12])


#  --------------- C4H2 1221: M = 43, I = 1 ---------------------
M = 43
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.57628E+03, 0.84874E+03, 0.11789E+04,
                                 0.15952E+04, 0.21317E+04, 0.28324E+04, 0.37543E+04, 0.49705E+04,
                                 0.65754E+04, 0.86894E+04, 0.11466E+05, 0.15099E+05, 0.19834E+05,
                                 0.25980E+05, 0.33920E+05, 0.44132E+05, 0.57210E+05, 0.73884E+05,
                                 0.95049E+05, 0.12180E+06, 0.15548E+06, 0.19771E+06, 0.25045E+06,
                                 0.31606E+06, 0.39739E+06, 0.49786E+06, 0.62152E+06, 0.77324E+06,
                                 0.95878E+06, 0.11850E+07, 0.14599E+07, 0.17930E+07, 0.21956E+07,
                                 0.26807E+07, 0.32637E+07, 0.39626E+07, 0.47983E+07, 0.57951E+07,
                                 0.69813E+07, 0.83896E+07, 0.10058E+08, 0.12030E+08, 0.14356E+08,
                                 0.17093E+08, 0.20309E+08, 0.24079E+08, 0.28491E+08, 0.33644E+08,
                                 0.39651E+08, 0.46642E+08, 0.54764E+08, 0.64184E+08, 0.75091E+08,
                                 0.87699E+08, 0.10225E+09, 0.11902E+09, 0.13832E+09, 0.16049E+09,
                                 0.18593E+09, 0.21507E+09, 0.24841E+09, 0.28650E+09, 0.32996E+09,
                                 0.37949E+09, 0.43586E+09, 0.49993E+09, 0.57266E+09, 0.65513E+09,
                                 0.74852E+09, 0.85418E+09, 0.97356E+09, 0.11083E+10, 0.12602E+10,
                                 0.14313E+10, 0.16238E+10, 0.18401E+10, 0.20829E+10, 0.23553E+10,
                                 0.26605E+10, 0.30021E+10, 0.33841E+10, 0.38109E+10, 0.42874E+10,
                                 0.48187E+10, 0.54107E+10, 0.60698E+10, 0.68029E+10, 0.76176E+10,
                                 0.85223E+10, 0.95260E+10, 0.10639E+11, 0.11871E+11, 0.13236E+11,
                                 0.14744E+11, 0.16412E+11, 0.18253E+11, 0.20285E+11, 0.22526E+11,
                                 0.24995E+11, 0.27714E+11, 0.30705E+11, 0.33995E+11, 0.37609E+11,
                                 0.41579E+11, 0.45934E+11, 0.50711E+11, 0.55947E+11, 0.61681E+11,
                                 0.67957E+11, 0.74824E+11, 0.82330E+11, 0.90532E+11, 0.99487E+11,
                                 0.10926E+12, 0.11992E+12, 0.13154E+12, 0.14420E+12, 0.15799E+12,
                                 0.17299E+12])


#  --------------- HC3N 12224: M = 44, I = 1 --------------------- 1224 in HITRAN, 12224 in TIPS
M = 44
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.16683E+04, 0.24538E+04, 0.33995E+04,
                                 0.45769E+04, 0.60637E+04, 0.79533E+04, 0.10360E+05, 0.13422E+05,
                                 0.17311E+05, 0.22232E+05, 0.28434E+05, 0.36215E+05, 0.45932E+05,
                                 0.58011E+05, 0.72958E+05, 0.91370E+05, 0.11395E+06, 0.14153E+06,
                                 0.17507E+06, 0.21570E+06, 0.26475E+06, 0.32372E+06, 0.39440E+06,
                                 0.47881E+06, 0.57930E+06, 0.69856E+06, 0.83968E+06, 0.10062E+07,
                                 0.12021E+07, 0.14320E+07, 0.17011E+07, 0.20153E+07, 0.23812E+07,
                                 0.28065E+07, 0.32996E+07, 0.38701E+07, 0.45287E+07, 0.52876E+07,
                                 0.61602E+07, 0.71616E+07, 0.83088E+07, 0.96206E+07, 0.11118E+08,
                                 0.12824E+08, 0.14765E+08, 0.16969E+08, 0.19469E+08, 0.22299E+08,
                                 0.25498E+08, 0.29110E+08, 0.33181E+08, 0.37763E+08, 0.42914E+08,
                                 0.48697E+08, 0.55180E+08, 0.62440E+08, 0.70558E+08, 0.79627E+08,
                                 0.89743E+08, 0.10102E+09, 0.11356E+09, 0.12752E+09, 0.14301E+09,
                                 0.16020E+09, 0.17925E+09, 0.20035E+09, 0.22367E+09, 0.24945E+09,
                                 0.27790E+09, 0.30928E+09, 0.34385E+09, 0.38191E+09, 0.42376E+09,
                                 0.46975E+09, 0.52023E+09, 0.57562E+09, 0.63632E+09, 0.70279E+09,
                                 0.77553E+09, 0.85506E+09, 0.94195E+09, 0.10368E+10, 0.11403E+10,
                                 0.12531E+10, 0.13759E+10, 0.15097E+10, 0.16552E+10, 0.18133E+10,
                                 0.19851E+10, 0.21715E+10, 0.23738E+10, 0.25931E+10, 0.28307E+10,
                                 0.30879E+10, 0.33662E+10, 0.36672E+10, 0.39926E+10, 0.43439E+10,
                                 0.47233E+10, 0.51325E+10, 0.55738E+10, 0.60493E+10, 0.65615E+10,
                                 0.71129E+10, 0.77061E+10, 0.83441E+10, 0.90298E+10, 0.97664E+10,
                                 0.10557E+11, 0.11406E+11, 0.12317E+11, 0.13293E+11, 0.14339E+11,
                                 0.15459E+11, 0.16659E+11, 0.17942E+11, 0.19316E+11, 0.20784E+11,
                                 0.22353E+11])


#  --------------- HC3N 12234: M = 44, I = 2 --------------------- see above
M = 44
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.33507E+04, 0.49290E+04, 0.68293E+04,
                                 0.91959E+04, 0.12185E+05, 0.15986E+05, 0.20828E+05, 0.26993E+05,
                                 0.34824E+05, 0.44739E+05, 0.57239E+05, 0.72931E+05, 0.92539E+05,
                                 0.11693E+06, 0.14713E+06, 0.18435E+06, 0.23004E+06, 0.28588E+06,
                                 0.35384E+06, 0.43625E+06, 0.53580E+06, 0.65562E+06, 0.79933E+06,
                                 0.97115E+06, 0.11759E+07, 0.14191E+07, 0.17073E+07, 0.20476E+07,
                                 0.24486E+07, 0.29196E+07, 0.34716E+07, 0.41169E+07, 0.48696E+07,
                                 0.57453E+07, 0.67621E+07, 0.79402E+07, 0.93022E+07, 0.10874E+08,
                                 0.12684E+08, 0.14764E+08, 0.17150E+08, 0.19884E+08, 0.23009E+08,
                                 0.26576E+08, 0.30641E+08, 0.35265E+08, 0.40518E+08, 0.46477E+08,
                                 0.53225E+08, 0.60856E+08, 0.69475E+08, 0.79195E+08, 0.90143E+08,
                                 0.10246E+09, 0.11629E+09, 0.13182E+09, 0.14921E+09, 0.16868E+09,
                                 0.19045E+09, 0.21477E+09, 0.24189E+09, 0.27211E+09, 0.30575E+09,
                                 0.34316E+09, 0.38471E+09, 0.43083E+09, 0.48196E+09, 0.53858E+09,
                                 0.60125E+09, 0.67052E+09, 0.74704E+09, 0.83148E+09, 0.92459E+09,
                                 0.10272E+10, 0.11401E+10, 0.12643E+10, 0.14007E+10, 0.15506E+10,
                                 0.17150E+10, 0.18953E+10, 0.20928E+10, 0.23090E+10, 0.25456E+10,
                                 0.28042E+10, 0.30867E+10, 0.33951E+10, 0.37316E+10, 0.40984E+10,
                                 0.44981E+10, 0.49332E+10, 0.54067E+10, 0.59216E+10, 0.64812E+10,
                                 0.70890E+10, 0.77488E+10, 0.84645E+10, 0.92405E+10, 0.10081E+11,
                                 0.10992E+11, 0.11978E+11, 0.13044E+11, 0.14197E+11, 0.15443E+11,
                                 0.16789E+11, 0.18243E+11, 0.19810E+11, 0.21501E+11, 0.23324E+11,
                                 0.25288E+11, 0.27403E+11, 0.29680E+11, 0.32130E+11, 0.34764E+11,
                                 0.37596E+11, 0.40639E+11, 0.43907E+11, 0.47416E+11, 0.51181E+11,
                                 0.55220E+11])


#  --------------- HC3N 12324: M = 44, I = 3 --------------------- see above
M = 44
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.33506E+04, 0.49280E+04, 0.68267E+04,
                                 0.91901E+04, 0.12174E+05, 0.15966E+05, 0.20793E+05, 0.26936E+05,
                                 0.34734E+05, 0.44598E+05, 0.57026E+05, 0.72612E+05, 0.92071E+05,
                                 0.11625E+06, 0.14616E+06, 0.18298E+06, 0.22813E+06, 0.28323E+06,
                                 0.35022E+06, 0.43133E+06, 0.52918E+06, 0.64677E+06, 0.78761E+06,
                                 0.95571E+06, 0.11557E+07, 0.13929E+07, 0.16734E+07, 0.20041E+07,
                                 0.23929E+07, 0.28488E+07, 0.33820E+07, 0.40040E+07, 0.47280E+07,
                                 0.55686E+07, 0.65423E+07, 0.76678E+07, 0.89661E+07, 0.10460E+08,
                                 0.12177E+08, 0.14145E+08, 0.16397E+08, 0.18970E+08, 0.21903E+08,
                                 0.25242E+08, 0.29036E+08, 0.33339E+08, 0.38214E+08, 0.43726E+08,
                                 0.49949E+08, 0.56965E+08, 0.64864E+08, 0.73743E+08, 0.83711E+08,
                                 0.94886E+08, 0.10740E+09, 0.12139E+09, 0.13701E+09, 0.15443E+09,
                                 0.17384E+09, 0.19543E+09, 0.21943E+09, 0.24607E+09, 0.27561E+09,
                                 0.30832E+09, 0.34452E+09, 0.38453E+09, 0.42870E+09, 0.47742E+09,
                                 0.53110E+09, 0.59020E+09, 0.65518E+09, 0.72659E+09, 0.80496E+09,
                                 0.89092E+09, 0.98510E+09, 0.10882E+10, 0.12010E+10, 0.13242E+10,
                                 0.14588E+10, 0.16056E+10, 0.17657E+10, 0.19401E+10, 0.21299E+10,
                                 0.23363E+10, 0.25606E+10, 0.28043E+10, 0.30687E+10, 0.33553E+10,
                                 0.36660E+10, 0.40024E+10, 0.43665E+10, 0.47601E+10, 0.51856E+10,
                                 0.56450E+10, 0.61408E+10, 0.66756E+10, 0.72520E+10, 0.78729E+10,
                                 0.85413E+10, 0.92604E+10, 0.10034E+11, 0.10864E+11, 0.11757E+11,
                                 0.12714E+11, 0.13742E+11, 0.14843E+11, 0.16023E+11, 0.17287E+11,
                                 0.18640E+11, 0.20087E+11, 0.21634E+11, 0.23288E+11, 0.25054E+11,
                                 0.26939E+11, 0.28950E+11, 0.31096E+11, 0.33382E+11, 0.35819E+11,
                                 0.38413E+11])


#  --------------- HC3N 13224: M = 44, I = 4 --------------------- see above
M = 44
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(12.)
TIPS_ISO_HASH[(M, I)] = float32([0.34439E+04, 0.50672E+04, 0.70230E+04,
                                 0.94603E+04, 0.12542E+05, 0.16462E+05, 0.21461E+05, 0.27833E+05,
                                 0.35935E+05, 0.46204E+05, 0.59168E+05, 0.75463E+05, 0.95854E+05,
                                 0.12126E+06, 0.15276E+06, 0.19165E+06, 0.23947E+06, 0.29802E+06,
                                 0.36943E+06, 0.45619E+06, 0.56121E+06, 0.68789E+06, 0.84018E+06,
                                 0.10227E+07, 0.12407E+07, 0.15003E+07, 0.18086E+07, 0.21738E+07,
                                 0.26052E+07, 0.31134E+07, 0.37106E+07, 0.44109E+07, 0.52300E+07,
                                 0.61861E+07, 0.72996E+07, 0.85939E+07, 0.10095E+08, 0.11833E+08,
                                 0.13841E+08, 0.16158E+08, 0.18825E+08, 0.21890E+08, 0.25407E+08,
                                 0.29436E+08, 0.34045E+08, 0.39308E+08, 0.45309E+08, 0.52143E+08,
                                 0.59912E+08, 0.68734E+08, 0.78737E+08, 0.90065E+08, 0.10288E+09,
                                 0.11735E+09, 0.13367E+09, 0.15206E+09, 0.17277E+09, 0.19604E+09,
                                 0.22217E+09, 0.25148E+09, 0.28432E+09, 0.32108E+09, 0.36218E+09,
                                 0.40809E+09, 0.45932E+09, 0.51644E+09, 0.58004E+09, 0.65082E+09,
                                 0.72950E+09, 0.81690E+09, 0.91388E+09, 0.10214E+10, 0.11405E+10,
                                 0.12724E+10, 0.14182E+10, 0.15794E+10, 0.17573E+10, 0.19536E+10,
                                 0.21701E+10, 0.24086E+10, 0.26711E+10, 0.29599E+10, 0.32774E+10,
                                 0.36262E+10, 0.40090E+10, 0.44290E+10, 0.48895E+10, 0.53939E+10,
                                 0.59462E+10, 0.65504E+10, 0.72111E+10, 0.79332E+10, 0.87217E+10,
                                 0.95823E+10, 0.10521E+11, 0.11544E+11, 0.12659E+11, 0.13874E+11,
                                 0.15195E+11, 0.16632E+11, 0.18194E+11, 0.19892E+11, 0.21735E+11,
                                 0.23736E+11, 0.25907E+11, 0.28260E+11, 0.30810E+11, 0.33572E+11,
                                 0.36563E+11, 0.39799E+11, 0.43299E+11, 0.47083E+11, 0.51172E+11,
                                 0.55588E+11, 0.60355E+11, 0.65500E+11, 0.71049E+11, 0.77031E+11,
                                 0.83478E+11])


#  --------------- HC3N 12225: M = 44, I = 5 --------------------- see above
M = 44
I = 5
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.11455E+04, 0.16850E+04, 0.23345E+04,
                                 0.31432E+04, 0.41647E+04, 0.54630E+04, 0.71168E+04, 0.92219E+04,
                                 0.11895E+05, 0.15279E+05, 0.19545E+05, 0.24897E+05, 0.31584E+05,
                                 0.39899E+05, 0.50190E+05, 0.62871E+05, 0.78428E+05, 0.97434E+05,
                                 0.12056E+06, 0.14859E+06, 0.18243E+06, 0.22314E+06, 0.27194E+06,
                                 0.33026E+06, 0.39972E+06, 0.48219E+06, 0.57983E+06, 0.69509E+06,
                                 0.83077E+06, 0.99009E+06, 0.11767E+07, 0.13946E+07, 0.16487E+07,
                                 0.19441E+07, 0.22868E+07, 0.26836E+07, 0.31420E+07, 0.36704E+07,
                                 0.42786E+07, 0.49770E+07, 0.57776E+07, 0.66938E+07, 0.77404E+07,
                                 0.89339E+07, 0.10293E+08, 0.11837E+08, 0.13590E+08, 0.15576E+08,
                                 0.17823E+08, 0.20362E+08, 0.23227E+08, 0.26454E+08, 0.30085E+08,
                                 0.34166E+08, 0.38745E+08, 0.43877E+08, 0.49622E+08, 0.56046E+08,
                                 0.63219E+08, 0.71222E+08, 0.80138E+08, 0.90062E+08, 0.10110E+09,
                                 0.11335E+09, 0.12695E+09, 0.14202E+09, 0.15870E+09, 0.17716E+09,
                                 0.19756E+09, 0.22009E+09, 0.24493E+09, 0.27232E+09, 0.30247E+09,
                                 0.33565E+09, 0.37211E+09, 0.41217E+09, 0.45613E+09, 0.50433E+09,
                                 0.55714E+09, 0.61497E+09, 0.67823E+09, 0.74739E+09, 0.82293E+09,
                                 0.90540E+09, 0.99536E+09, 0.10934E+10, 0.12002E+10, 0.13165E+10,
                                 0.14430E+10, 0.15805E+10, 0.17299E+10, 0.18922E+10, 0.20682E+10,
                                 0.22591E+10, 0.24660E+10, 0.26901E+10, 0.29326E+10, 0.31951E+10,
                                 0.34788E+10, 0.37854E+10, 0.41166E+10, 0.44741E+10, 0.48598E+10,
                                 0.52758E+10, 0.57240E+10, 0.62069E+10, 0.67269E+10, 0.72864E+10,
                                 0.78882E+10, 0.85352E+10, 0.92305E+10, 0.99773E+10, 0.10779E+11,
                                 0.11639E+11, 0.12562E+11, 0.13552E+11, 0.14612E+11, 0.15748E+11,
                                 0.16964E+11])


#  --------------- HC3N 22224: M = 44, I = 6 --------------------- see above
M = 44
I = 6
TIPS_GSI_HASH[(M, I)] = __FloatType__(9.)
TIPS_ISO_HASH[(M, I)] = float32([0.27029E+04, 0.39999E+04, 0.55894E+04,
                                 0.76092E+04, 0.10219E+05, 0.13616E+05, 0.18042E+05, 0.23798E+05,
                                 0.31255E+05, 0.40867E+05, 0.53189E+05, 0.68897E+05, 0.88807E+05,
                                 0.11390E+06, 0.14537E+06, 0.18461E+06, 0.23330E+06, 0.29342E+06,
                                 0.36733E+06, 0.45779E+06, 0.56802E+06, 0.70182E+06, 0.86361E+06,
                                 0.10585E+07, 0.12925E+07, 0.15725E+07, 0.19064E+07, 0.23034E+07,
                                 0.27739E+07, 0.33302E+07, 0.39858E+07, 0.47566E+07, 0.56604E+07,
                                 0.67176E+07, 0.79511E+07, 0.93872E+07, 0.11055E+08, 0.12989E+08,
                                 0.15225E+08, 0.17806E+08, 0.20779E+08, 0.24197E+08, 0.28119E+08,
                                 0.32612E+08, 0.37749E+08, 0.43612E+08, 0.50294E+08, 0.57895E+08,
                                 0.66528E+08, 0.76318E+08, 0.87403E+08, 0.99937E+08, 0.11409E+09,
                                 0.13004E+09, 0.14800E+09, 0.16819E+09, 0.19086E+09, 0.21629E+09,
                                 0.24476E+09, 0.27661E+09, 0.31219E+09, 0.35189E+09, 0.39615E+09,
                                 0.44542E+09, 0.50021E+09, 0.56108E+09, 0.62862E+09, 0.70350E+09,
                                 0.78641E+09, 0.87814E+09, 0.97952E+09, 0.10915E+10, 0.12149E+10,
                                 0.13510E+10, 0.15008E+10, 0.16656E+10, 0.18468E+10, 0.20457E+10,
                                 0.22640E+10, 0.25032E+10, 0.27653E+10, 0.30522E+10, 0.33659E+10,
                                 0.37088E+10, 0.40832E+10, 0.44917E+10, 0.49371E+10, 0.54224E+10,
                                 0.59508E+10, 0.65256E+10, 0.71507E+10, 0.78298E+10, 0.85671E+10,
                                 0.93672E+10, 0.10235E+11, 0.11175E+11, 0.12193E+11, 0.13295E+11,
                                 0.14487E+11, 0.15776E+11, 0.17168E+11, 0.18671E+11, 0.20293E+11,
                                 0.22043E+11, 0.23929E+11, 0.25960E+11, 0.28148E+11, 0.30502E+11,
                                 0.33034E+11, 0.35756E+11, 0.38681E+11, 0.41823E+11, 0.45195E+11,
                                 0.48812E+11, 0.52692E+11, 0.56850E+11, 0.61306E+11, 0.66076E+11,
                                 0.71183E+11])


#  --------------- H2 11: M = 45, I = 1 ---------------------
M = 45
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.15265E+01, 0.22243E+01, 0.29619E+01,
                                 0.36724E+01, 0.43456E+01, 0.49880E+01, 0.56090E+01, 0.62165E+01,
                                 0.68161E+01, 0.74113E+01, 0.80044E+01, 0.85966E+01, 0.91887E+01,
                                 0.97810E+01, 0.10374E+02, 0.10967E+02, 0.11561E+02, 0.12156E+02,
                                 0.12751E+02, 0.13347E+02, 0.13944E+02, 0.14541E+02, 0.15139E+02,
                                 0.15738E+02, 0.16337E+02, 0.16937E+02, 0.17538E+02, 0.18140E+02,
                                 0.18743E+02, 0.19346E+02, 0.19951E+02, 0.20556E+02, 0.21163E+02,
                                 0.21771E+02, 0.22379E+02, 0.22990E+02, 0.23601E+02, 0.24214E+02,
                                 0.24829E+02, 0.25445E+02, 0.26063E+02, 0.26683E+02, 0.27304E+02,
                                 0.27928E+02, 0.28553E+02, 0.29181E+02, 0.29811E+02, 0.30443E+02,
                                 0.31078E+02, 0.31715E+02, 0.32355E+02, 0.32997E+02, 0.33643E+02,
                                 0.34291E+02, 0.34942E+02, 0.35596E+02, 0.36253E+02, 0.36914E+02,
                                 0.37578E+02, 0.38245E+02, 0.38916E+02, 0.39590E+02, 0.40268E+02,
                                 0.40949E+02, 0.41635E+02, 0.42324E+02, 0.43017E+02, 0.43715E+02,
                                 0.44416E+02, 0.45122E+02, 0.45831E+02, 0.46546E+02, 0.47264E+02,
                                 0.47987E+02, 0.48714E+02, 0.49446E+02, 0.50183E+02, 0.50925E+02,
                                 0.51671E+02, 0.52422E+02, 0.53178E+02, 0.53939E+02, 0.54705E+02,
                                 0.55476E+02, 0.56252E+02, 0.57033E+02, 0.57820E+02, 0.58612E+02,
                                 0.59409E+02, 0.60212E+02, 0.61020E+02, 0.61833E+02, 0.62652E+02,
                                 0.63477E+02, 0.64308E+02, 0.65144E+02, 0.65986E+02, 0.66833E+02,
                                 0.67687E+02, 0.68546E+02, 0.69411E+02, 0.70283E+02, 0.71160E+02,
                                 0.72043E+02, 0.72933E+02, 0.73829E+02, 0.74730E+02, 0.75638E+02,
                                 0.76553E+02, 0.77473E+02, 0.78400E+02, 0.79333E+02, 0.80273E+02,
                                 0.81219E+02, 0.82172E+02, 0.83131E+02, 0.84097E+02, 0.85069E+02,
                                 0.86048E+02])


#  --------------- H2 12: M = 45, I = 2 ---------------------
M = 45
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(6.)
TIPS_ISO_HASH[(M, I)] = float32([0.81692E+01, 0.10308E+02, 0.12557E+02,
                                 0.14848E+02, 0.17159E+02, 0.19482E+02, 0.21815E+02, 0.24153E+02,
                                 0.26497E+02, 0.28845E+02, 0.31197E+02, 0.33552E+02, 0.35910E+02,
                                 0.38272E+02, 0.40636E+02, 0.43002E+02, 0.45372E+02, 0.47744E+02,
                                 0.50119E+02, 0.52496E+02, 0.54877E+02, 0.57261E+02, 0.59649E+02,
                                 0.62040E+02, 0.64435E+02, 0.66835E+02, 0.69240E+02, 0.71650E+02,
                                 0.74066E+02, 0.76489E+02, 0.78918E+02, 0.81354E+02, 0.83799E+02,
                                 0.86252E+02, 0.88715E+02, 0.91187E+02, 0.93669E+02, 0.96163E+02,
                                 0.98668E+02, 0.10118E+03, 0.10371E+03, 0.10626E+03, 0.10881E+03,
                                 0.11138E+03, 0.11397E+03, 0.11657E+03, 0.11919E+03, 0.12182E+03,
                                 0.12447E+03, 0.12714E+03, 0.12982E+03, 0.13252E+03, 0.13524E+03,
                                 0.13798E+03, 0.14074E+03, 0.14352E+03, 0.14632E+03, 0.14914E+03,
                                 0.15198E+03, 0.15484E+03, 0.15772E+03, 0.16062E+03, 0.16355E+03,
                                 0.16649E+03, 0.16946E+03, 0.17246E+03, 0.17547E+03, 0.17851E+03,
                                 0.18157E+03, 0.18466E+03, 0.18777E+03, 0.19090E+03, 0.19406E+03,
                                 0.19725E+03, 0.20045E+03, 0.20369E+03, 0.20695E+03, 0.21023E+03,
                                 0.21354E+03, 0.21687E+03, 0.22024E+03, 0.22362E+03, 0.22704E+03,
                                 0.23048E+03, 0.23394E+03, 0.23744E+03, 0.24096E+03, 0.24451E+03,
                                 0.24808E+03, 0.25169E+03, 0.25532E+03, 0.25897E+03, 0.26266E+03,
                                 0.26638E+03, 0.27012E+03, 0.27389E+03, 0.27769E+03, 0.28152E+03,
                                 0.28537E+03, 0.28926E+03, 0.29317E+03, 0.29712E+03, 0.30109E+03,
                                 0.30509E+03, 0.30913E+03, 0.31319E+03, 0.31728E+03, 0.32140E+03,
                                 0.32555E+03, 0.32974E+03, 0.33395E+03, 0.33819E+03, 0.34246E+03,
                                 0.34677E+03, 0.35110E+03, 0.35547E+03, 0.35987E+03, 0.36429E+03,
                                 0.36875E+03])


#  --------------- CS 22: M = 46, I = 1 ---------------------
M = 46
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.51416E+02, 0.72723E+02, 0.94044E+02,
                                 0.11538E+03, 0.13673E+03, 0.15810E+03, 0.17949E+03, 0.20093E+03,
                                 0.22245E+03, 0.24407E+03, 0.26582E+03, 0.28776E+03, 0.30992E+03,
                                 0.33233E+03, 0.35504E+03, 0.37807E+03, 0.40147E+03, 0.42525E+03,
                                 0.44944E+03, 0.47406E+03, 0.49914E+03, 0.52468E+03, 0.55071E+03,
                                 0.57723E+03, 0.60427E+03, 0.63183E+03, 0.65991E+03, 0.68854E+03,
                                 0.71771E+03, 0.74743E+03, 0.77771E+03, 0.80855E+03, 0.83996E+03,
                                 0.87193E+03, 0.90449E+03, 0.93762E+03, 0.97134E+03, 0.10056E+04,
                                 0.10405E+04, 0.10760E+04, 0.11121E+04, 0.11487E+04, 0.11860E+04,
                                 0.12239E+04, 0.12623E+04, 0.13014E+04, 0.13410E+04, 0.13813E+04,
                                 0.14222E+04, 0.14637E+04, 0.15057E+04, 0.15484E+04, 0.15917E+04,
                                 0.16357E+04, 0.16802E+04, 0.17253E+04, 0.17711E+04, 0.18175E+04,
                                 0.18645E+04, 0.19121E+04, 0.19603E+04, 0.20091E+04, 0.20586E+04,
                                 0.21087E+04, 0.21594E+04, 0.22107E+04, 0.22626E+04, 0.23152E+04,
                                 0.23684E+04, 0.24222E+04, 0.24767E+04, 0.25317E+04, 0.25874E+04,
                                 0.26438E+04, 0.27007E+04, 0.27583E+04, 0.28165E+04, 0.28754E+04,
                                 0.29348E+04, 0.29949E+04, 0.30557E+04, 0.31170E+04, 0.31790E+04,
                                 0.32417E+04, 0.33049E+04, 0.33688E+04, 0.34334E+04, 0.34986E+04,
                                 0.35644E+04, 0.36308E+04, 0.36979E+04, 0.37656E+04, 0.38340E+04,
                                 0.39030E+04, 0.39727E+04, 0.40430E+04, 0.41139E+04, 0.41855E+04,
                                 0.42577E+04, 0.43306E+04, 0.44041E+04, 0.44782E+04, 0.45530E+04,
                                 0.46284E+04, 0.47045E+04, 0.47813E+04, 0.48587E+04, 0.49367E+04,
                                 0.50154E+04, 0.50947E+04, 0.51747E+04, 0.52553E+04, 0.53366E+04,
                                 0.54185E+04, 0.55011E+04, 0.55844E+04, 0.56683E+04, 0.57528E+04,
                                 0.58380E+04])


#  --------------- CS 24: M = 46, I = 2 ---------------------
M = 46
I = 2
TIPS_GSI_HASH[(M, I)] = __FloatType__(1.)
TIPS_ISO_HASH[(M, I)] = float32([0.52247E+02, 0.73900E+02, 0.95568E+02,
                                 0.11725E+03, 0.13895E+03, 0.16066E+03, 0.18241E+03, 0.20420E+03,
                                 0.22607E+03, 0.24805E+03, 0.27018E+03, 0.29249E+03, 0.31503E+03,
                                 0.33784E+03, 0.36096E+03, 0.38442E+03, 0.40824E+03, 0.43247E+03,
                                 0.45712E+03, 0.48221E+03, 0.50778E+03, 0.53382E+03, 0.56037E+03,
                                 0.58743E+03, 0.61501E+03, 0.64312E+03, 0.67179E+03, 0.70100E+03,
                                 0.73077E+03, 0.76111E+03, 0.79202E+03, 0.82351E+03, 0.85559E+03,
                                 0.88824E+03, 0.92149E+03, 0.95533E+03, 0.98977E+03, 0.10248E+04,
                                 0.10605E+04, 0.10967E+04, 0.11336E+04, 0.11710E+04, 0.12091E+04,
                                 0.12478E+04, 0.12871E+04, 0.13270E+04, 0.13675E+04, 0.14087E+04,
                                 0.14505E+04, 0.14929E+04, 0.15359E+04, 0.15795E+04, 0.16238E+04,
                                 0.16687E+04, 0.17142E+04, 0.17604E+04, 0.18071E+04, 0.18546E+04,
                                 0.19026E+04, 0.19513E+04, 0.20006E+04, 0.20505E+04, 0.21011E+04,
                                 0.21523E+04, 0.22042E+04, 0.22566E+04, 0.23098E+04, 0.23635E+04,
                                 0.24179E+04, 0.24730E+04, 0.25286E+04, 0.25850E+04, 0.26419E+04,
                                 0.26995E+04, 0.27578E+04, 0.28167E+04, 0.28762E+04, 0.29364E+04,
                                 0.29972E+04, 0.30587E+04, 0.31208E+04, 0.31836E+04, 0.32470E+04,
                                 0.33111E+04, 0.33758E+04, 0.34412E+04, 0.35072E+04, 0.35739E+04,
                                 0.36412E+04, 0.37092E+04, 0.37778E+04, 0.38471E+04, 0.39171E+04,
                                 0.39877E+04, 0.40589E+04, 0.41309E+04, 0.42034E+04, 0.42767E+04,
                                 0.43505E+04, 0.44251E+04, 0.45003E+04, 0.45762E+04, 0.46527E+04,
                                 0.47299E+04, 0.48077E+04, 0.48863E+04, 0.49654E+04, 0.50453E+04,
                                 0.51258E+04, 0.52070E+04, 0.52888E+04, 0.53713E+04, 0.54545E+04,
                                 0.55383E+04, 0.56229E+04, 0.57080E+04, 0.57939E+04, 0.58804E+04,
                                 0.59676E+04])


#  --------------- CS 32: M = 46, I = 3 ---------------------
M = 46
I = 3
TIPS_GSI_HASH[(M, I)] = __FloatType__(2.)
TIPS_ISO_HASH[(M, I)] = float32([0.10889E+03, 0.15403E+03, 0.19920E+03,
                                 0.24440E+03, 0.28964E+03, 0.33491E+03, 0.38026E+03, 0.42571E+03,
                                 0.47134E+03, 0.51722E+03, 0.56342E+03, 0.61005E+03, 0.65719E+03,
                                 0.70493E+03, 0.75334E+03, 0.80249E+03, 0.85245E+03, 0.90329E+03,
                                 0.95504E+03, 0.10078E+04, 0.10615E+04, 0.11163E+04, 0.11721E+04,
                                 0.12291E+04, 0.12872E+04, 0.13464E+04, 0.14068E+04, 0.14684E+04,
                                 0.15311E+04, 0.15951E+04, 0.16604E+04, 0.17268E+04, 0.17945E+04,
                                 0.18635E+04, 0.19337E+04, 0.20051E+04, 0.20779E+04, 0.21519E+04,
                                 0.22272E+04, 0.23038E+04, 0.23817E+04, 0.24609E+04, 0.25414E+04,
                                 0.26232E+04, 0.27064E+04, 0.27908E+04, 0.28765E+04, 0.29636E+04,
                                 0.30520E+04, 0.31417E+04, 0.32327E+04, 0.33251E+04, 0.34188E+04,
                                 0.35138E+04, 0.36102E+04, 0.37079E+04, 0.38070E+04, 0.39074E+04,
                                 0.40091E+04, 0.41122E+04, 0.42166E+04, 0.43224E+04, 0.44295E+04,
                                 0.45380E+04, 0.46478E+04, 0.47590E+04, 0.48715E+04, 0.49854E+04,
                                 0.51007E+04, 0.52173E+04, 0.53353E+04, 0.54547E+04, 0.55754E+04,
                                 0.56975E+04, 0.58210E+04, 0.59458E+04, 0.60720E+04, 0.61996E+04,
                                 0.63285E+04, 0.64589E+04, 0.65906E+04, 0.67236E+04, 0.68581E+04,
                                 0.69940E+04, 0.71312E+04, 0.72698E+04, 0.74098E+04, 0.75512E+04,
                                 0.76940E+04, 0.78381E+04, 0.79837E+04, 0.81307E+04, 0.82790E+04,
                                 0.84287E+04, 0.85799E+04, 0.87324E+04, 0.88864E+04, 0.90417E+04,
                                 0.91984E+04, 0.93566E+04, 0.95161E+04, 0.96771E+04, 0.98394E+04,
                                 0.10003E+05, 0.10168E+05, 0.10335E+05, 0.10503E+05, 0.10672E+05,
                                 0.10843E+05, 0.11015E+05, 0.11189E+05, 0.11364E+05, 0.11541E+05,
                                 0.11719E+05, 0.11898E+05, 0.12079E+05, 0.12261E+05, 0.12444E+05,
                                 0.12630E+05])


#  --------------- CS 23: M = 46, I = 4 ---------------------
M = 46
I = 4
TIPS_GSI_HASH[(M, I)] = __FloatType__(4.)
TIPS_ISO_HASH[(M, I)] = float32([0.20737E+03, 0.29330E+03, 0.37930E+03,
                                 0.46535E+03, 0.55145E+03, 0.63764E+03, 0.72394E+03, 0.81043E+03,
                                 0.89722E+03, 0.98443E+03, 0.10722E+04, 0.11607E+04, 0.12501E+04,
                                 0.13406E+04, 0.14323E+04, 0.15253E+04, 0.16197E+04, 0.17158E+04,
                                 0.18135E+04, 0.19129E+04, 0.20142E+04, 0.21174E+04, 0.22226E+04,
                                 0.23298E+04, 0.24391E+04, 0.25504E+04, 0.26639E+04, 0.27796E+04,
                                 0.28976E+04, 0.30177E+04, 0.31401E+04, 0.32648E+04, 0.33918E+04,
                                 0.35211E+04, 0.36527E+04, 0.37867E+04, 0.39231E+04, 0.40618E+04,
                                 0.42029E+04, 0.43463E+04, 0.44922E+04, 0.46405E+04, 0.47912E+04,
                                 0.49443E+04, 0.50999E+04, 0.52579E+04, 0.54183E+04, 0.55812E+04,
                                 0.57465E+04, 0.59143E+04, 0.60846E+04, 0.62573E+04, 0.64325E+04,
                                 0.66102E+04, 0.67903E+04, 0.69729E+04, 0.71581E+04, 0.73457E+04,
                                 0.75358E+04, 0.77284E+04, 0.79235E+04, 0.81211E+04, 0.83212E+04,
                                 0.85239E+04, 0.87290E+04, 0.89367E+04, 0.91469E+04, 0.93596E+04,
                                 0.95748E+04, 0.97926E+04, 0.10013E+05, 0.10236E+05, 0.10461E+05,
                                 0.10689E+05, 0.10920E+05, 0.11153E+05, 0.11388E+05, 0.11626E+05,
                                 0.11867E+05, 0.12110E+05, 0.12356E+05, 0.12604E+05, 0.12855E+05,
                                 0.13109E+05, 0.13365E+05, 0.13623E+05, 0.13884E+05, 0.14148E+05,
                                 0.14415E+05, 0.14683E+05, 0.14955E+05, 0.15229E+05, 0.15506E+05,
                                 0.15785E+05, 0.16067E+05, 0.16351E+05, 0.16638E+05, 0.16928E+05,
                                 0.17220E+05, 0.17515E+05, 0.17813E+05, 0.18113E+05, 0.18416E+05,
                                 0.18721E+05, 0.19029E+05, 0.19340E+05, 0.19653E+05, 0.19969E+05,
                                 0.20287E+05, 0.20608E+05, 0.20932E+05, 0.21258E+05, 0.21587E+05,
                                 0.21919E+05, 0.22253E+05, 0.22590E+05, 0.22930E+05, 0.23272E+05,
                                 0.23617E+05])


#  --------------- SO3 26: M = 46, I = 1 --------------------- not in TIPS-2011
M = 47
I = 1
TIPS_GSI_HASH[(M, I)] = __FloatType__(0.)
TIPS_ISO_HASH[(M, I)] = float32([0.])


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
    if T < 70. or T > 3000.:
        #Qt = -1.
        #gi = 0.
        # return gi,Qt
        raise Exception('TIPS: T must be between 70K and 3000K.')

    try:
        # get statistical weight for specified isotopologue
        gi = TIPS_GSI_HASH[(M, I)]
        # interpolate partition sum for specified isotopologue
        Qt = AtoB(T, Tdat, TIPS_ISO_HASH[(M, I)], TIPS_NPT)
    except KeyError:
        raise Exception('TIPS: no data for M,I = %d,%d.' % (M, I))

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
        #n = (T[1]-T[0])/step
        #TT = linspace(T[0],T[1],n)
        TT = arange(T[0], T[1], step)
        return TT, array([BD_TIPS_2011_PYTHON(M, I, temp)[1] for temp in TT])

# ------------------ partition sum --------------------------------------


# ------------------ LINESHAPES -----------------------------------------

# ------------------ complex probability function -----------------------
# define static data
zone = __ComplexType__(1.0e0 + 0.0e0j)
zi = __ComplexType__(0.0e0 + 1.0e0j)
tt = __FloatType__([0.5e0, 1.5e0, 2.5e0, 3.5e0, 4.5e0, 5.5e0, 6.5e0,
                    7.5e0, 8.5e0, 9.5e0, 10.5e0, 11.5e0, 12.5e0, 13.5e0, 14.5e0])
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

    zm1 = zone/__ComplexType__(X + zi*Y)  # maybe redundant
    zm2 = zm1**2
    zsum = zone
    zterm = zone

    for tt_i in tt:
        zterm *= zm2*tt_i
        zsum += zterm

    zsum *= zi*zm1*pipwoeronehalf

    return zsum.real, zsum.imag


T = __FloatType__([0.314240376e0, 0.947788391e0, 1.59768264e0,
                   2.27950708e0, 3.02063703e0, 3.8897249e0])
U = __FloatType__([1.01172805e0, -0.75197147e0, 1.2557727e-2,
                   1.00220082e-2, -2.42068135e-4, 5.00848061e-7])
S = __FloatType__([1.393237e0, 0.231152406e0, -0.155351466e0,
                   6.21836624e-3, 9.19082986e-5, -6.27525958e-7])

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
    index_REGION3 = where(sqrt(X**2 + Y**2) > __FloatType__(8.0e0))
    X_REGION3 = X[index_REGION3]
    Y_REGION3 = Y[index_REGION3]
    zm1 = zone/__ComplexType__(X_REGION3 + zi*Y_REGION3)
    zm2 = zm1**2
    zsum_REGION3 = zone
    zterm = zone
    for tt_i in tt:
        zterm *= zm2*tt_i
        zsum_REGION3 += zterm
    zsum_REGION3 *= zi*zm1*pipwoeronehalf

    index_REGION12 = setdiff1d(array(arange(len(X))), array(index_REGION3))
    X_REGION12 = X[index_REGION12]
    Y_REGION12 = Y[index_REGION12]

    WR = __FloatType__(0.0e0)
    WI = __FloatType__(0.0e0)

    # REGION12
    Y1_REGION12 = Y_REGION12 + __FloatType__(1.5e0)
    Y2_REGION12 = Y1_REGION12**2

    # REGION2
    subindex_REGION2 = where((Y_REGION12 <= 0.85e0) &
                             (abs(X_REGION12) >= (18.1e0*Y_REGION12 + 1.65e0)))

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
    WR_REGION2[ii] = exp(-X_REGION2[ii]**2)
    WR_REGION2[~ii] = WR

    for I in range(6):
        R_REGION2 = X_REGION2 - T[I]
        R2_REGION2 = R_REGION2**2
        D_REGION2 = __FloatType__(1.0e0) / (R2_REGION2 + Y2_REGION2)
        D1_REGION2 = Y1_REGION2 * D_REGION2
        D2_REGION2 = R_REGION2 * D_REGION2
        WR_REGION2 = WR_REGION2 + Y_REGION2 * (U[I]*(R_REGION2*D2_REGION2 - 1.5e0*D1_REGION2) +
                                               S[I]*Y3_REGION2*D2_REGION2)/(R2_REGION2 + 2.25e0)
        R_REGION2 = X_REGION2 + T[I]
        R2_REGION2 = R_REGION2**2
        D_REGION2 = __FloatType__(1.0e0) / (R2_REGION2 + Y2_REGION2)
        D3_REGION2 = Y1_REGION2 * D_REGION2
        D4_REGION2 = R_REGION2 * D_REGION2
        WR_REGION2 = WR_REGION2 + Y_REGION2 * (U[I]*(R_REGION2*D4_REGION2 - 1.5e0*D3_REGION2) -
                                               S[I]*Y3_REGION2*D4_REGION2)/(R2_REGION2 + 2.25e0)
        WI_REGION2 = WI_REGION2 + \
            U[I]*(D2_REGION2 + D4_REGION2) + S[I]*(D1_REGION2 - D3_REGION2)

    # REGION3
    index_REGION1 = setdiff1d(array(index_REGION12), array(index_REGION2))
    X_REGION1 = X[index_REGION1]
    Y_REGION1 = X[index_REGION1]

    subindex_REGION1 = setdiff1d(
        array(arange(len(index_REGION12))), array(subindex_REGION2))
    Y1_REGION1 = Y1_REGION12[subindex_REGION1]
    Y2_REGION1 = Y2_REGION12[subindex_REGION1]

    WR_REGION1 = WR
    WI_REGION1 = WI

    for I in range(6):
        R_REGION1 = X_REGION1 - T[I]
        D_REGION1 = __FloatType__(1.0e0) / (R_REGION1**2 + Y2_REGION1)
        D1_REGION1 = Y1_REGION1 * D_REGION1
        D2_REGION1 = R_REGION1 * D_REGION1
        R_REGION1 = X_REGION1 + T[I]
        D_REGION1 = __FloatType__(1.0e0) / (R_REGION1**2 + Y2_REGION1)
        D3_REGION1 = Y1_REGION1 * D_REGION1
        D4_REGION1 = R_REGION1 * D_REGION1

        WR_REGION1 = WR_REGION1 + \
            U[I]*(D1_REGION1 + D3_REGION1) - S[I]*(D2_REGION1 - D4_REGION1)
        WI_REGION1 = WI_REGION1 + \
            U[I]*(D2_REGION1 + D4_REGION1) + S[I]*(D1_REGION1 - D3_REGION1)

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
    z = x + 1.0j*y
    M = 2*N
    M2 = 2*M
    k = arange(-M+1, M)  # '; # M2 = no. of sampling points.
    L = sqrt(N/sqrt(2))  # Optimal choice of L.
    theta = k*pi/M
    t = L*tan(theta/2)  # Variables theta and t.
    # f = exp(-t.A2)*(LA2+t.A2); f = [0; f]; # Function to be transformed.
    f = zeros(len(t)+1)
    f[0] = 0
    f[1:] = exp(-t**2)*(L**2+t**2)
    #f = insert(exp(-t**2)*(L**2+t**2),0,0)
    a = real(fft(fftshift(f)))/M2  # Coefficients of transform.
    a = flipud(a[1:N+1])  # Reorder coefficients.
    Z = (L+1.0j*z)/(L-1.0j*z)
    p = polyval(a, Z)  # Polynomial evaluation.
    w = 2*p/(L-1.0j*z)**2+(1/sqrt(pi))/(L-1.0j*z)  # Evaluate w(z).
    return w


# weideman24 by default
#weideman24 = lambda x,y: cef(x,y,24)
def weideman(x, y, n): return cef(x, y, n)


def hum1_wei(x, y, n=24):
    t = y-1.0j*x
    cerf = 1/sqrt(pi)*t/(0.5+t**2)
    """
    z = x+1j*y
    cerf = 1j*z/sqrt(pi)/(z**2-0.5)
    """
    mask = abs(x)+y < 15.0
    if any(mask):
        w24 = weideman(x[mask], y[mask], n)
        place(cerf, mask, w24)
    return cerf.real, cerf.imag


VARIABLES['CPF'] = hum1_wei
#VARIABLES['CPF'] = cpf

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

    cte = sqrt(log(2.0e0))/GamD
    rpi = sqrt(pi)
    iz = __ComplexType__(0.0e0 + 1.0e0j)

    c0 = __ComplexType__(Gam0 + 1.0e0j*Shift0)
    c2 = __ComplexType__(Gam2 + 1.0e0j*Shift2)
    c0t = __ComplexType__((1.0e0 - eta) * (c0 - 1.5e0 * c2) + anuVC)
    c2t = __ComplexType__((1.0e0 - eta) * c2)

    # PART1
    if abs(c2t) == 0.0e0:
        # print('PART1') # DEBUG
        Z1 = (iz*(sg0 - sg) + c0t) * cte
        xZ1 = -Z1.imag
        yZ1 = Z1.real
        #WR1,WI1 = cpf(xZ1,yZ1)
        WR1, WI1 = VARIABLES['CPF'](xZ1, yZ1)
        Aterm_GLOBAL = rpi*cte*__ComplexType__(WR1 + 1.0e0j*WI1)
        index_Z1 = abs(Z1) <= 4.0e3
        index_NOT_Z1 = ~index_Z1
        if any(index_Z1):
            # print('PART1/Z1') # DEBUG
            Bterm_GLOBAL = rpi*cte * \
                ((1.0e0 - Z1**2)*__ComplexType__(WR1 + 1.0e0j*WI1) + Z1/rpi)
        if any(index_NOT_Z1):
            # print('PART1/~Z1') # DEBUG
            Bterm_GLOBAL = cte * \
                (rpi*__ComplexType__(WR1 + 1.0e0j*WI1) + 0.5e0/Z1 - 0.75e0/(Z1**3))
    else:
        # PART2, PART3 AND PART4   (PART4 IS A MAIN PART)

        # X - vector, Y - scalar
        X = (iz * (sg0 - sg) + c0t) / c2t
        Y = __ComplexType__(1.0e0 / ((2.0e0*cte*c2t))**2)
        csqrtY = (Gam2 - iz*Shift2) / (2.0e0*cte *
                                       (1.0e0-eta) * (Gam2**2 + Shift2**2))

        index_PART2 = abs(X) <= 3.0e-8 * abs(Y)
        index_PART3 = (abs(Y) <= 1.0e-15 * abs(X)) & ~index_PART2
        index_PART4 = ~ (index_PART2 | index_PART3)

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
            SZ1 = sqrt(xZ1**2 + yZ1**2)
            SZ2 = sqrt(xZ2**2 + yZ2**2)
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
                #WR1,WI1 = cpf(xZ1[index_CPF],yZ1[index_CPF])
                #WR2,WI2 = cpf(xZ2[index_CPF],yZ2[index_CPF])
                WR1, WI1 = VARIABLES['CPF'](xZ1[index_CPF], yZ1[index_CPF])
                WR2, WI2 = VARIABLES['CPF'](xZ2[index_CPF], yZ2[index_CPF])
                WR1_PART4[index_CPF] = WR1
                WI1_PART4[index_CPF] = WI1
                WR2_PART4[index_CPF] = WR2
                WI2_PART4[index_CPF] = WI2

            Aterm = rpi*cte*(__ComplexType__(WR1_PART4 + 1.0e0j *
                                             WI1_PART4) - __ComplexType__(WR2_PART4+1.0e0j*WI2_PART4))
            Bterm = (-1.0e0 +
                     rpi/(2.0e0*csqrtY)*(1.0e0 - Z1**2)*__ComplexType__(WR1_PART4 + 1.0e0j*WI1_PART4) -
                     rpi/(2.0e0*csqrtY)*(1.0e0 - Z2**2)*__ComplexType__(WR2_PART4 + 1.0e0j*WI2_PART4)) / c2t
            Aterm_GLOBAL[index_PART4] = Aterm
            Bterm_GLOBAL[index_PART4] = Bterm

        # PART2
        if any(index_PART2):
            # print('PART2') # DEBUG
            X_TMP = X[index_PART2]
            Z1 = (iz*(sg0 - sg[index_PART2]) + c0t) * cte
            Z2 = sqrt(X_TMP + Y) + csqrtY
            xZ1 = -Z1.imag
            yZ1 = Z1.real
            xZ2 = -Z2.imag
            yZ2 = Z2.real
            #WR1_PART2,WI1_PART2 = cpf(xZ1,yZ1)
            #WR2_PART2,WI2_PART2 = cpf(xZ2,yZ2)
            WR1_PART2, WI1_PART2 = VARIABLES['CPF'](xZ1, yZ1)
            WR2_PART2, WI2_PART2 = VARIABLES['CPF'](xZ2, yZ2)
            Aterm = rpi*cte*(__ComplexType__(WR1_PART2 + 1.0e0j*WI1_PART2) -
                             __ComplexType__(WR2_PART2 + 1.0e0j*WI2_PART2))
            Bterm = (-1.0e0 +
                     rpi/(2.0e0*csqrtY)*(1.0e0 - Z1**2)*__ComplexType__(WR1_PART2 + 1.0e0j*WI1_PART2) -
                     rpi/(2.0e0*csqrtY)*(1.0e0 - Z2**2)*__ComplexType__(WR2_PART2 + 1.0e0j*WI2_PART2)) / c2t
            Aterm_GLOBAL[index_PART2] = Aterm
            Bterm_GLOBAL[index_PART2] = Bterm

        # PART3
        if any(index_PART3):
            # print('PART3') # DEBUG
            X_TMP = X[index_PART3]
            xZ1 = -sqrt(X_TMP + Y).imag
            yZ1 = sqrt(X_TMP + Y).real
            #WR1_PART3,WI1_PART3 =  cpf(xZ1,yZ1)
            WR1_PART3, WI1_PART3 = VARIABLES['CPF'](xZ1, yZ1)
            index_ABS = abs(sqrt(X_TMP)) <= 4.0e3
            index_NOT_ABS = ~index_ABS
            Aterm = zeros(len(index_PART3), dtype=__ComplexType__)
            Bterm = zeros(len(index_PART3), dtype=__ComplexType__)
            if any(index_ABS):
                xXb = -sqrt(X).imag
                yXb = sqrt(X).real
                #WRb,WIb = cpf(xXb,yXb)
                WRb, WIb = VARIABLES['CPF'](xXb, yXb)
                Aterm[index_ABS] = (
                    2.0e0*rpi/c2t)*(1.0e0/rpi - sqrt(X_TMP[index_ABS])*__ComplexType__(WRb + 1.0e0j*WIb))
                Bterm[index_ABS] = (1.0e0/c2t)*(-1.0e0 +
                                                2.0e0*rpi*(1.0e0 - X_TMP[index_ABS]-2.0e0*Y)*(1.0e0/rpi-sqrt(X_TMP[index_ABS])*__ComplexType__(WRb + 1.0e0j*WIb)) +
                                                2.0e0*rpi*sqrt(X_TMP[index_ABS] + Y)*__ComplexType__(WR1_PART3 + 1.0e0j*WI1_PART3))
            if any(index_NOT_ABS):
                Aterm[index_NOT_ABS] = (
                    1.0e0/c2t)*(1.0e0/X_TMP[index_NOT_ABS] - 1.5e0/(X_TMP[index_NOT_ABS]**2))
                Bterm[index_NOT_ABS] = (1.0e0/c2t)*(-1.0e0 + (1.0e0 - X_TMP[index_NOT_ABS] - 2.0e0*Y) *
                                                    (1.0e0/X_TMP[index_NOT_ABS] - 1.5e0/(X_TMP[index_NOT_ABS]**2)) +
                                                    2.0e0*rpi*sqrt(X_TMP[index_NOT_ABS] + Y)*__ComplexType__(WR1 + 1.0e0j*WI1))
            Aterm_GLOBAL[index_PART3] = Aterm
            Bterm_GLOBAL[index_PART3] = Bterm

    # common part
    LS_pCqSDHC = (1.0e0/pi) * (Aterm_GLOBAL / (1.0e0 -
                                               (anuVC-eta*(c0-1.5e0*c2))*Aterm_GLOBAL + eta*c2*Bterm_GLOBAL))
    # print('pcqsdc_end',sum(LS_pCqSDHC.real),sum(LS_pCqSDHC.imag))
    return LS_pCqSDHC.real, LS_pCqSDHC.imag


# ------------------  CROSS-SECTIONS, XSECT.PY --------------------------------

# set interfaces for TIPS(M,I,T)
def PYTIPS(M, I, T): return BD_TIPS_2011_PYTHON(M, I, T)[1]

# set interfaces for profiles
#PYHTP = pcqsdhc
#PROFILE_HTP = PYHTP
#PROFILE_VOIGT = lambda sg0,GamD,Gam0,sg: PROFILE_HTP(sg0,GamD,Gam0,cZero,cZero,cZero,cZero,cZero,sg)
#PROFILE_LORENTZ = lambda sg0,Gam0,sg: Gam0/(pi*(Gam0**2+(sg-sg0)**2))
#PROFILE_DOPPLER = lambda sg0,GamD,sg: cSqrtLn2divSqrtPi*exp(-cLn2*((sg-sg0)/GamD)**2)/GamD


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
    return Gam0/(pi*(Gam0**2+(sg-sg0)**2))


def PROFILE_DOPPLER(sg0, GamD, sg):
    """
    # Doppler profile.
    # Input parameters:
    #   sg0: Unperturbed line position in cm-1 (Input).
    #   GamD: Doppler HWHM in cm-1 (Input)
    #   sg: Current WaveNumber of the Computation in cm-1 (Input).
    """
    return cSqrtLn2divSqrtPi*exp(-cLn2*((sg-sg0)/GamD)**2)/GamD

# Volume concentration of all gas molecules at the pressure p and temperature T


def volumeConcentration(p, T):
    return (p/9.869233e-7)/(cBolts*T)  # CGS

# ------------------------------- PARAMETER DEPENDENCIES --------------------------------

# temperature dependence for intencities (HITRAN)


def EnvironmentDependency_Intensity(LineIntensityRef, T, Tref, SigmaT, SigmaTref,
                                    LowerStateEnergy, LineCenter):
    const = __FloatType__(1.4388028496642257)
    ch = exp(-const*LowerStateEnergy/T)*(1-exp(-const*LineCenter/T))
    zn = exp(-const*LowerStateEnergy/Tref)*(1-exp(-const*LineCenter/Tref))
    LineIntensity = LineIntensityRef*SigmaTref/SigmaT*ch/zn
    return LineIntensity

# environmental dependence for GammaD (HTP, Voigt)    # Tref/T ????


def EnvironmentDependency_GammaD(GammaD_ref, T, Tref):
    # Doppler parameters do not depend on pressure!
    return GammaD_ref*sqrt(T/Tref)

# environmental dependence for Gamma0 (HTP, Voigt)


def EnvironmentDependency_Gamma0(Gamma0_ref, T, Tref, p, pref, TempRatioPower):
    return Gamma0_ref*p/pref*(Tref/T)**TempRatioPower

# environmental dependence for Gamma2 (HTP)


def EnvironmentDependency_Gamma2(Gamma2_ref, T, Tref, p, pref, TempRatioPower):
    return Gamma2_ref*p/pref*(Tref/T)**TempRatioPower

# environmental dependence for Delta0 (HTP)


def EnvironmentDependency_Delta0(Delta0_ref, p, pref):
    return Delta0_ref*p/pref

# environmental dependence for Delta2 (HTP)


def EnvironmentDependency_Delta2(Delta2_ref, p, pref):
    return Delta2_ref*p/pref

# environmental dependence for anuVC (HTP)


def EnvironmentDependency_anuVC(anuVC_ref, T, Tref, p, pref):
    return anuVC_ref*Tref/T*p/pref

# ------------------------------- /PARAMETER DEPENDENCIES --------------------------------

# ------------------------------- BINGINGS --------------------------------


# default parameter bindings
DefaultParameterBindings = {}

# default temperature dependencies
DefaultEnvironmentDependencyBindings = {}

# ------------------------------- /BINGINGS --------------------------------

# default values for intensity threshold
DefaultIntensityThreshold = 0.  # cm*molec

# default value for omega wing in halfwidths (from center)
DefaultOmegaWingHW = 50.  # cm-1    HOTW default


# check and argument for being a tuple or list
# this is connected with a "bug" that in Python
# (val) is not a tuple, but (val,) is a tuple
def listOfTuples(a):
    if type(a) not in set([list, tuple]):
        a = [a]
    return a


# determine default parameters from those which are passed to absorptionCoefficient_...
def getDefaultValuesForXsect(Components, SourceTables, Environment, OmegaRange,
                             OmegaStep, OmegaWing, IntensityThreshold, Format):
    if SourceTables[0] == None:
        SourceTables = ['__BUFFER__', ]
    if Environment == None:
        Environment = {'T': 296., 'p': 1.}
    if Components == [None]:
        CompDict = {}
        for TableName in SourceTables:
            # check table existance
            if TableName not in list(LOCAL_TABLE_CACHE.keys()):
                raise Exception(
                    '%s: no such table. Check tableList() for more info.' % TableName)
            mol_ids = LOCAL_TABLE_CACHE[TableName]['data']['molec_id']
            iso_ids = LOCAL_TABLE_CACHE[TableName]['data']['local_iso_id']
            if len(mol_ids) != len(iso_ids):
                raise Exception('Lengths if mol_ids and iso_ids differ!')
            MI_zip = list(zip(mol_ids, iso_ids))
            MI_zip = set(MI_zip)
            for mol_id, iso_id in MI_zip:
                CompDict[(mol_id, iso_id)] = None
        Components = list(CompDict.keys())
    if OmegaRange == None:
        omega_min = float('inf')
        omega_max = float('-inf')
        for TableName in SourceTables:
            nu = LOCAL_TABLE_CACHE[TableName]['data']['nu']
            numin = min(nu)
            numax = max(nu)
            if omega_min > numin:
                omega_min = numin
            if omega_max < numax:
                omega_max = numax
        OmegaRange = (omega_min, omega_max)
    OmegaDelta = OmegaRange[1]-OmegaRange[0]
    if OmegaStep == None:
        #OmegaStep = OmegaDelta/100.
        OmegaStep = 0.01  # cm-1
    if OmegaWing == None:
        #OmegaWing = OmegaDelta/10.
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
        Format = '%.12f %e'
    return Components, SourceTables, Environment, OmegaRange,\
        OmegaStep, OmegaWing, IntensityThreshold, Format


# save numpy arrays to file
# arrays must have same dimensions
def save_to_file(fname, fformat, *arg):
    f = open(fname, 'w')
    for i in range(len(arg[0])):
        argline = []
        for j in range(len(arg)):
            argline.append(arg[j][i])
        f.write((fformat+'\n') % tuple(argline))
    f.close()

# ==========================================================================================
# =========================== NEW ABSORPTION COEFFICIENT ===================================
# ==========================================================================================


def absorptionCoefficient_HT(Components=None, SourceTables=None, partitionFunction=PYTIPS,
                             Environment=None, OmegaRange=None, OmegaStep=None, OmegaWing=None,
                             IntensityThreshold=DefaultIntensityThreshold,
                             OmegaWingHW=DefaultOmegaWingHW,
                             GammaL='gamma_air', HITRAN_units=True, LineShift=True,
                             File=None, Format=None, OmegaGrid=None,
                             WavenumberRange=None, WavenumberStep=None, WavenumberWing=None,
                             WavenumberWingHW=None, WavenumberGrid=None,
                             Diluent={}, EnvDependences=None):
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

    warn('To get the most up-to-date version please check http://hitran.org/hapi')

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
    Components, SourceTables, Environment, OmegaRange, OmegaStep, OmegaWing,\
        IntensityThreshold, Format = \
        getDefaultValuesForXsect(Components, SourceTables, Environment, OmegaRange,
                                 OmegaStep, OmegaWing, IntensityThreshold, Format)

    # warn user about too large omega step
    if OmegaStep > 0.1:
        warn('Big wavenumber step: possible accuracy decline')

    # get uniform linespace for cross-section
    #number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    #Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        #Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.)  # K
    pref = __FloatType__(1.)  # atm

    # actual temperature and pressure
    T = Environment['T']  # K
    p = Environment['p']  # atm

    # Find reference temperature
    TRanges = [(0, 100), (100, 200), (200, 400), (400, float('inf'))]
    Trefs = [50., 150., 296., 700.]
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
                ni = ISO[(M, I)][ISO_INDEX['abundance']]
            except KeyError:
                raise Exception('cannot find component M,I = %d,%d.' % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX['abundance']]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # setup the default empty environment dependence function
    if not EnvDependences:
        def EnvDependences(ENV, LINE): return {}
    Env = Environment.copy()
    Env['Tref'] = Tref
    Env['pref'] = pref

    # setup the Diluent variable
    GammaL = GammaL.lower()
    if not Diluent:
        if GammaL == 'gamma_air':
            Diluent = {'air': 1.}
        elif GammaL == 'gamma_self':
            Diluent = {'self': 1.}
        else:
            raise Exception('Unknown GammaL value: %s' % GammaL)

    # Simple check
    for key in Diluent:
        val = Diluent[key]
        if val < 0 and val > 1:
            raise Exception('Diluent fraction must be in [0,1]')

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get the number of rows
        nline = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']

        # get parameter names for each table
        parnames = list(LOCAL_TABLE_CACHE[TableName]['data'].keys())

        # loop through line centers (single stream)
        for RowID in range(nline):

            # Get the custom environment dependences
            Line = {parname: LOCAL_TABLE_CACHE[TableName]
                    ['data'][parname][RowID] for parname in parnames}
            CustomEnvDependences = EnvDependences(Env, Line)

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]['data']['nu'][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]['data']['sw'][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]['data']['elower'][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['molec_id'][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['local_iso_id'][RowID]

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES:
                continue

            # partition functions for T and Tref
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            if 'sw' in CustomEnvDependences:
                LineIntensity = CustomEnvDependences['sw']
            else:
                LineIntensity = EnvironmentDependency_Intensity(LineIntensityDB, T, Tref, SigmaT, SigmaTref,
                                                                LowerStateEnergyDB, LineCenterDB)

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            if LineIntensity < IntensityThreshold:
                continue

            #   doppler broadening coefficient (GammaD)
            cMassMol = 1.66053873e-27  # hapi
            m = molecularMass(MoleculeNumberDB, IsoNumberDB) * cMassMol * 1000
            GammaD = sqrt(2*cBolts*T*log(2)/m/cc**2)*LineCenterDB

            #   pressure broadening coefficient
            Gamma0 = 0.
            Shift0 = 0.
            Gamma2 = 0.
            Shift2 = 0.
            Delta2 = 0.
            NuVC = 0.
            EtaNumer = 0.
            for species in Diluent:
                species_lower = species.lower()

                abun = Diluent[species]

                # Search for broadening HWHM.
                try:
                    # search for HT-style name
                    Gamma0DB = LOCAL_TABLE_CACHE[TableName]['data']['gamma_HT_0_%s_%d' % (
                        species_lower, TrefHT)][RowID]
                    if Gamma0DB == 0.:
                        raise Exception
                except:
                    try:
                        # search for Voigt-style name
                        Gamma0DB = LOCAL_TABLE_CACHE[TableName]['data']['gamma_%s' %
                                                                        species_lower][RowID]
                    except:
                        Gamma0DB = 0.0

                # Search for temperature exponent for broadening HWHM.
                try:
                    # search for HT-style name
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_HT_%s_%d' % (
                        species_lower, TrefHT)][RowID]
                    if TempRatioPowerDB == 0.:
                        raise Exception
                    Tref = TrefHT
                except:
                    Tref = 296.
                    try:
                        # search for Voigt-style name
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_%s' %
                                                                                species_lower][RowID]
                        if species_lower == 'self' and TempRatioPowerDB == 0.:
                            # same for self as for air
                            TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_air'][RowID]
                    except:
                        #print('TempRatioPowerDB is set to zero')
                        #TempRatioPowerDB = 0
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_air'][RowID]

                # Add to the final Gamma0
                Gamma0T = CustomEnvDependences.get('gamma_HT_0_%s_%d' % (species_lower, TrefHT),
                                                   CustomEnvDependences.get('gamma_%s' % species_lower,
                                                                            EnvironmentDependency_Gamma0(Gamma0DB, T, Tref, p, pref, TempRatioPowerDB)))
                Gamma0 += abun*Gamma0T

                # Search for shift.
                try:
                    # search for HT-style name
                    Shift0DB = LOCAL_TABLE_CACHE[TableName]['data']['delta_HT_0_%s_%d' % (
                        species_lower, TrefHT)][RowID]
                    if Shift0DB == 0.:
                        raise Exception
                except:
                    try:
                        # search for Voigt-style name
                        Shift0DB = LOCAL_TABLE_CACHE[TableName]['data']['delta_%s' %
                                                                        species_lower][RowID]
                    except:
                        Shift0DB = 0.0

                # Search for temperature dependence for shift.
                try:
                    # search for HT-style name
                    deltap = LOCAL_TABLE_CACHE[TableName]['data']['deltap_HT_%s_%d' % (
                        species_lower, TrefHT)][RowID]
                    if deltap == 0.:
                        raise Exception
                    Tref = TrefHT
                except:
                    Tref = 296.
                    try:
                        # search for Voigt-style name
                        deltap = LOCAL_TABLE_CACHE[TableName]['data']['deltap_%s' %
                                                                      species_lower][RowID]
                    except:
                        deltap = 0.0

                Shift0T = CustomEnvDependences.get('deltap_HT_%s_%d' % (species_lower, TrefHT),
                                                   CustomEnvDependences.get('deltap_%s' % species_lower,
                                                                            ((Shift0DB + deltap*(T-Tref))*p/pref)))
                Shift0 += abun*Shift0T

                # Search for speed dependence for HWHM.
                try:
                    Gamma2DB = LOCAL_TABLE_CACHE[TableName]['data']['gamma_HT_2_%s_%d' % (
                        species_lower, TrefHT)][RowID]
                    if Gamma2DB == 0.:
                        raise Exception
                except:
                    try:
                        SDDB = LOCAL_TABLE_CACHE[TableName]['data']['SD_%s' %
                                                                    species_lower][RowID]
                        Gamma2DB = SDDB*Gamma0DB
                    except:
                        Gamma2DB = 0.0

                Gamma2 += abun * \
                    CustomEnvDependences.get('gamma_HT_2_%s_%d' % (
                        species_lower, TrefHT), Gamma2DB*(p/pref))

                # Search for speed dependence for shift.
                try:
                    Delta2DB = LOCAL_TABLE_CACHE[TableName]['data']['delta_HT_2_%s_%d' % (
                        species_lower, TrefHT)][RowID]
                except:
                    Delta2DB = 0.

                Delta2 += abun*CustomEnvDependences.get('delta_HT_2_%s_%d' % (species_lower, TrefHT),
                                                        Delta2DB*p/pref)

                # Search for frequency of VC
                try:
                    NuVCDB = LOCAL_TABLE_CACHE[TableName]['data']['nu_HT_%s' %
                                                                  species_lower][RowID]
                except:
                    NuVCDB = 0.

                # Search for temperature exponent for frequency of VC
                try:
                    KappaDB = LOCAL_TABLE_CACHE[TableName]['data']['kappa_HT_%s' %
                                                                   species_lower][RowID]
                except:
                    KappaDB = 0.

                NuVC += abun*CustomEnvDependences.get('nu_HT_%s' % species_lower,
                                                      NuVCDB*(Tref/T)**KappaDB*p)

                # Setup correlation parameter
                try:
                    EtaDB = LOCAL_TABLE_CACHE[TableName]['data']['eta_HT_%s' %
                                                                 species_lower][RowID]
                except:
                    EtaDB = 0.

                EtaNumer += EtaDB*abun*(Gamma0T+1j*Shift0T)

            Eta = EtaNumer/(Gamma0 + 1j*Shift0)

            #   get final wing of the line according to Gamma0, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW*Gamma0, OmegaWingHW*GammaD)

            #   shift coefficient
            #Shift0 = Shift0DB*p/pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB-OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB+OmegaWingF)
            lineshape_vals = PROFILE_HT(LineCenterDB, GammaD, Gamma0, Gamma2, Shift0,
                                        Shift2, NuVC, Eta, Omegas[BoundIndexLower:BoundIndexUpper])[0]
            Xsect[BoundIndexLower:BoundIndexUpper] += factor / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                LineIntensity * lineshape_vals
            # print(LineCenterDB,GammaD,Gamma0,Gamma2,Shift0,sum(lineshape_vals),sum(Xsect))
            #raise Exception

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def absorptionCoefficient_SDVoigt(Components=None, SourceTables=None, partitionFunction=PYTIPS,
                                  Environment=None, OmegaRange=None, OmegaStep=None, OmegaWing=None,
                                  IntensityThreshold=DefaultIntensityThreshold,
                                  OmegaWingHW=DefaultOmegaWingHW,
                                  GammaL='gamma_air', HITRAN_units=True, LineShift=True,
                                  File=None, Format=None, OmegaGrid=None,
                                  WavenumberRange=None, WavenumberStep=None, WavenumberWing=None,
                                  WavenumberWingHW=None, WavenumberGrid=None,
                                  Diluent={}, EnvDependences=None):
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

    warn('To get the most up-to-date version please check http://hitran.org/hapi')

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
    Components, SourceTables, Environment, OmegaRange, OmegaStep, OmegaWing,\
        IntensityThreshold, Format = \
        getDefaultValuesForXsect(Components, SourceTables, Environment, OmegaRange,
                                 OmegaStep, OmegaWing, IntensityThreshold, Format)

    # warn user about too large omega step
    if OmegaStep > 0.1:
        warn('Big wavenumber step: possible accuracy decline')

    # get uniform linespace for cross-section
    #number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    #Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        #Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.)  # K
    pref = __FloatType__(1.)  # atm

    # actual temperature and pressure
    T = Environment['T']  # K
    p = Environment['p']  # atm

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
                ni = ISO[(M, I)][ISO_INDEX['abundance']]
            except KeyError:
                raise Exception('cannot find component M,I = %d,%d.' % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX['abundance']]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # setup the default empty environment dependence function
    if not EnvDependences:
        def EnvDependences(ENV, LINE): return {}
    Env = Environment.copy()
    Env['Tref'] = Tref
    Env['pref'] = pref

    # setup the Diluent variable
    GammaL = GammaL.lower()
    if not Diluent:
        if GammaL == 'gamma_air':
            Diluent = {'air': 1.}
        elif GammaL == 'gamma_self':
            Diluent = {'self': 1.}
        else:
            raise Exception('Unknown GammaL value: %s' % GammaL)

    # Simple check
    for key in Diluent:
        val = Diluent[key]
        if val < 0 and val > 1:
            raise Exception('Diluent fraction must be in [0,1]')

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get the number of rows
        nline = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']

        # get parameter names for each table
        parnames = list(LOCAL_TABLE_CACHE[TableName]['data'].keys())

        # loop through line centers (single stream)
        for RowID in range(nline):

            # Get the custom environment dependences
            Line = {parname: LOCAL_TABLE_CACHE[TableName]
                    ['data'][parname][RowID] for parname in parnames}
            CustomEnvDependences = EnvDependences(Env, Line)

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]['data']['nu'][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]['data']['sw'][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]['data']['elower'][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['molec_id'][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['local_iso_id'][RowID]

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES:
                continue

            # partition functions for T and Tref
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            if 'sw' in CustomEnvDependences:
                LineIntensity = CustomEnvDependences['sw']
            else:
                LineIntensity = EnvironmentDependency_Intensity(LineIntensityDB, T, Tref, SigmaT, SigmaTref,
                                                                LowerStateEnergyDB, LineCenterDB)

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            if LineIntensity < IntensityThreshold:
                continue

            #   doppler broadening coefficient (GammaD)
            cMassMol = 1.66053873e-27  # hapi
            m = molecularMass(MoleculeNumberDB, IsoNumberDB) * cMassMol * 1000
            GammaD = sqrt(2*cBolts*T*log(2)/m/cc**2)*LineCenterDB

            #   pressure broadening coefficient
            Gamma0 = 0.
            Shift0 = 0.
            Gamma2 = 0.
            Shift2 = 0.
            for species in Diluent:
                species_lower = species.lower()

                abun = Diluent[species]

                gamma_name = 'gamma_' + species_lower
                try:
                    Gamma0DB = LOCAL_TABLE_CACHE[TableName]['data'][gamma_name][RowID]
                except:
                    Gamma0DB = 0.0

                n_name = 'n_' + species_lower
                try:
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data'][n_name][RowID]
                    if species_lower == 'self' and TempRatioPowerDB == 0.:
                        # same for self as for air
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_air'][RowID]
                except:
                    #TempRatioPowerDB = 0
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_air'][RowID]

                # Add to the final Gamma0
                Gamma0 += abun*CustomEnvDependences.get(gamma_name,  # default ->
                                                        EnvironmentDependency_Gamma0(Gamma0DB, T, Tref, p, pref, TempRatioPowerDB))

                delta_name = 'delta_' + species_lower
                try:
                    Shift0DB = LOCAL_TABLE_CACHE[TableName]['data'][delta_name][RowID]
                except:
                    Shift0DB = 0.0

                deltap_name = 'deltap_' + species_lower
                try:
                    deltap = LOCAL_TABLE_CACHE[TableName]['data'][deltap_name][RowID]
                except:
                    deltap = 0.0

                Shift0 += abun*CustomEnvDependences.get(delta_name,  # default ->
                                                        ((Shift0DB + deltap*(T-Tref))*p/pref))

                SD_name = 'SD_' + species_lower
                try:
                    SDDB = LOCAL_TABLE_CACHE[TableName]['data'][SD_name][RowID]
                except:
                    SDDB = 0.0

                Gamma2 += abun*CustomEnvDependences.get(SD_name,  # default ->
                                                        SDDB*p/pref) * Gamma0DB

                #   get final wing of the line according to Gamma0, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW*Gamma0, OmegaWingHW*GammaD)

            #   shift coefficient
            #Shift0 = Shift0DB*p/pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB-OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB+OmegaWingF)
            lineshape_vals = PROFILE_SDVOIGT(
                LineCenterDB, GammaD, Gamma0, Gamma2, Shift0, Shift2, Omegas[BoundIndexLower:BoundIndexUpper])[0]
            Xsect[BoundIndexLower:BoundIndexUpper] += factor / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                LineIntensity * lineshape_vals
            # print(LineCenterDB,GammaD,Gamma0,Gamma2,Shift0,sum(lineshape_vals),sum(Xsect))
            #raise Exception

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def absorptionCoefficient_Voigt(Components=None, SourceTables=None, partitionFunction=PYTIPS,
                                Environment=None, OmegaRange=None, OmegaStep=None, OmegaWing=None,
                                IntensityThreshold=DefaultIntensityThreshold,
                                OmegaWingHW=DefaultOmegaWingHW,
                                GammaL='gamma_air', HITRAN_units=True, LineShift=True,
                                File=None, Format=None, OmegaGrid=None,
                                WavenumberRange=None, WavenumberStep=None, WavenumberWing=None,
                                WavenumberWingHW=None, WavenumberGrid=None,
                                Diluent={}, EnvDependences=None):
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
    Components, SourceTables, Environment, OmegaRange, OmegaStep, OmegaWing,\
        IntensityThreshold, Format = \
        getDefaultValuesForXsect(Components, SourceTables, Environment, OmegaRange,
                                 OmegaStep, OmegaWing, IntensityThreshold, Format)

    # warn user about too large omega step
    if OmegaStep > 0.1:
        warn('Big wavenumber step: possible accuracy decline')

    # get uniform linespace for cross-section
    #number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    #Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        #Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.)  # K
    pref = __FloatType__(1.)  # atm

    # actual temperature and pressure
    T = Environment['T']  # K
    p = Environment['p']  # atm

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
                ni = ISO[(M, I)][ISO_INDEX['abundance']]
            except KeyError:
                raise Exception('cannot find component M,I = %d,%d.' % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX['abundance']]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # setup the default empty environment dependence function
    if not EnvDependences:
        def EnvDependences(ENV, LINE): return {}
    Env = Environment.copy()
    Env['Tref'] = Tref
    Env['pref'] = pref

    # setup the Diluent variable
    GammaL = GammaL.lower()
    if not Diluent:
        if GammaL == 'gamma_air':
            Diluent = {'air': 1.}
        elif GammaL == 'gamma_self':
            Diluent = {'self': 1.}
        else:
            raise Exception('Unknown GammaL value: %s' % GammaL)

    # Simple check
    for key in Diluent:
        val = Diluent[key]
        if val < 0 and val > 1:
            raise Exception('Diluent fraction must be in [0,1]')

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get the number of rows
        nline = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']

        # get parameter names for each table
        parnames = list(LOCAL_TABLE_CACHE[TableName]['data'].keys())

        # loop through line centers (single stream)
        for RowID in range(nline):

            # Get the custom environment dependences
            Line = {parname: LOCAL_TABLE_CACHE[TableName]
                    ['data'][parname][RowID] for parname in parnames}
            CustomEnvDependences = EnvDependences(Env, Line)

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]['data']['nu'][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]['data']['sw'][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]['data']['elower'][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['molec_id'][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['local_iso_id'][RowID]

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES:
                continue

            # partition functions for T and Tref
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            if 'sw' in CustomEnvDependences:
                LineIntensity = CustomEnvDependences['sw']
            else:
                LineIntensity = EnvironmentDependency_Intensity(LineIntensityDB, T, Tref, SigmaT, SigmaTref,
                                                                LowerStateEnergyDB, LineCenterDB)

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            if LineIntensity < IntensityThreshold:
                continue

            #   doppler broadening coefficient (GammaD)
            cMassMol = 1.66053873e-27  # hapi
            m = molecularMass(MoleculeNumberDB, IsoNumberDB) * cMassMol * 1000
            GammaD = sqrt(2*cBolts*T*log(2)/m/cc**2)*LineCenterDB

            #   pressure broadening coefficient
            Gamma0 = 0.
            Shift0 = 0.
            for species in Diluent:
                species_lower = species.lower()

                abun = Diluent[species]

                gamma_name = 'gamma_' + species_lower
                try:
                    Gamma0DB = LOCAL_TABLE_CACHE[TableName]['data'][gamma_name][RowID]
                except:
                    Gamma0DB = 0.0

                n_name = 'n_' + species_lower
                try:
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data'][n_name][RowID]
                    if species_lower == 'self' and TempRatioPowerDB == 0.:
                        # same for self as for air
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_air'][RowID]
                except:
                    #TempRatioPowerDB = 0
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_air'][RowID]

                # Add to the final Gamma0
                Gamma0 += abun*CustomEnvDependences.get(gamma_name,  # default ->
                                                        EnvironmentDependency_Gamma0(Gamma0DB, T, Tref, p, pref, TempRatioPowerDB))

                delta_name = 'delta_' + species_lower
                try:
                    Shift0DB = LOCAL_TABLE_CACHE[TableName]['data'][delta_name][RowID]
                except:
                    Shift0DB = 0.0

                deltap_name = 'deltap_' + species_lower
                try:
                    deltap = LOCAL_TABLE_CACHE[TableName]['data'][deltap_name][RowID]
                except:
                    deltap = 0.0

                Shift0 += abun*CustomEnvDependences.get(delta_name,  # default ->
                                                        ((Shift0DB + deltap*(T-Tref))*p/pref))

            #   get final wing of the line according to Gamma0, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW*Gamma0, OmegaWingHW*GammaD)

            #   shift coefficient
            #Shift0 = Shift0DB*p/pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB-OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB+OmegaWingF)
            lineshape_vals = PROFILE_VOIGT(
                LineCenterDB+Shift0, GammaD, Gamma0, Omegas[BoundIndexLower:BoundIndexUpper])[0]
            Xsect[BoundIndexLower:BoundIndexUpper] += factor / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                LineIntensity * lineshape_vals

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def absorptionCoefficient_Lorentz(Components=None, SourceTables=None, partitionFunction=PYTIPS,
                                  Environment=None, OmegaRange=None, OmegaStep=None, OmegaWing=None,
                                  IntensityThreshold=DefaultIntensityThreshold,
                                  OmegaWingHW=DefaultOmegaWingHW,
                                  GammaL='gamma_air', HITRAN_units=True, LineShift=True,
                                  File=None, Format=None, OmegaGrid=None,
                                  WavenumberRange=None, WavenumberStep=None, WavenumberWing=None,
                                  WavenumberWingHW=None, WavenumberGrid=None,
                                  Diluent={}, EnvDependences=None):
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
    Components, SourceTables, Environment, OmegaRange, OmegaStep, OmegaWing,\
        IntensityThreshold, Format = \
        getDefaultValuesForXsect(Components, SourceTables, Environment, OmegaRange,
                                 OmegaStep, OmegaWing, IntensityThreshold, Format)

    # warn user about too large omega step
    if OmegaStep > 0.1:
        warn('Big wavenumber step: possible accuracy decline')

    # get uniform linespace for cross-section
    #number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    #Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        #Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.)  # K
    pref = __FloatType__(1.)  # atm

    # actual temperature and pressure
    T = Environment['T']  # K
    p = Environment['p']  # atm

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
                ni = ISO[(M, I)][ISO_INDEX['abundance']]
            except KeyError:
                raise Exception('cannot find component M,I = %d,%d.' % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX['abundance']]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # setup the default empty environment dependence function
    if not EnvDependences:
        def EnvDependences(ENV, LINE): return {}
    Env = Environment.copy()
    Env['Tref'] = Tref
    Env['pref'] = pref

    # setup the Diluent variable
    GammaL = GammaL.lower()
    if not Diluent:
        if GammaL == 'gamma_air':
            Diluent = {'air': 1.}
        elif GammaL == 'gamma_self':
            Diluent = {'self': 1.}
        else:
            raise Exception('Unknown GammaL value: %s' % GammaL)

    # Simple check
    for key in Diluent:
        val = Diluent[key]
        if val < 0 and val > 1:
            raise Exception('Diluent fraction must be in [0,1]')

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get the number of rows
        nline = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']

        # get parameter names for each table
        parnames = list(LOCAL_TABLE_CACHE[TableName]['data'].keys())

        # loop through line centers (single stream)
        for RowID in range(nline):

            # Get the custom environment dependences
            Line = {parname: LOCAL_TABLE_CACHE[TableName]
                    ['data'][parname][RowID] for parname in parnames}
            CustomEnvDependences = EnvDependences(Env, Line)

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]['data']['nu'][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]['data']['sw'][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]['data']['elower'][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['molec_id'][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['local_iso_id'][RowID]

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES:
                continue

            # partition functions for T and Tref
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            if 'sw' in CustomEnvDependences:
                LineIntensity = CustomEnvDependences['sw']
            else:
                LineIntensity = EnvironmentDependency_Intensity(LineIntensityDB, T, Tref, SigmaT, SigmaTref,
                                                                LowerStateEnergyDB, LineCenterDB)

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            if LineIntensity < IntensityThreshold:
                continue

            #   pressure broadening coefficient
            Gamma0 = 0.
            Shift0 = 0.
            for species in Diluent:
                species_lower = species.lower()

                abun = Diluent[species]

                gamma_name = 'gamma_' + species_lower
                try:
                    Gamma0DB = LOCAL_TABLE_CACHE[TableName]['data'][gamma_name][RowID]
                except:
                    Gamma0DB = 0.0

                n_name = 'n_' + species_lower
                try:
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data'][n_name][RowID]
                    if species_lower == 'self' and TempRatioPowerDB == 0.:
                        # same for self as for air
                        TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_air'][RowID]
                except:
                    #TempRatioPowerDB = 0
                    TempRatioPowerDB = LOCAL_TABLE_CACHE[TableName]['data']['n_air'][RowID]

                # Add to the final Gamma0
                Gamma0 += abun*CustomEnvDependences.get(gamma_name,  # default ->
                                                        EnvironmentDependency_Gamma0(Gamma0DB, T, Tref, p, pref, TempRatioPowerDB))

                delta_name = 'delta_' + species_lower
                try:
                    Shift0DB = LOCAL_TABLE_CACHE[TableName]['data'][delta_name][RowID]
                except:
                    Shift0DB = 0.0

                deltap_name = 'deltap_' + species_lower
                try:
                    deltap = LOCAL_TABLE_CACHE[TableName]['data'][deltap_name][RowID]
                except:
                    deltap = 0.0

                Shift0 += abun*CustomEnvDependences.get(delta_name,  # default ->
                                                        ((Shift0DB + deltap*(T-Tref))*p/pref))

            #   get final wing of the line according to Gamma0, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW*Gamma0)

            #   shift coefficient
            #Shift0 = Shift0DB*p/pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB-OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB+OmegaWingF)
            lineshape_vals = PROFILE_LORENTZ(
                LineCenterDB+Shift0, Gamma0, Omegas[BoundIndexLower:BoundIndexUpper])
            Xsect[BoundIndexLower:BoundIndexUpper] += factor / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                LineIntensity * lineshape_vals

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


# Alias for a profile selector
absorptionCoefficient = absorptionCoefficient_HT

# ==========================================================================================
# =========================== /NEW ABSORPTION COEFFICIENT ===================================
# ==========================================================================================

# calculate apsorption for Doppler profile


def absorptionCoefficient_Doppler(Components=None, SourceTables=None, partitionFunction=PYTIPS,
                                  Environment=None, OmegaRange=None, OmegaStep=None, OmegaWing=None,
                                  IntensityThreshold=DefaultIntensityThreshold,
                                  OmegaWingHW=DefaultOmegaWingHW,
                                  ParameterBindings=DefaultParameterBindings,
                                  EnvironmentDependencyBindings=DefaultEnvironmentDependencyBindings,
                                  GammaL='dummy', HITRAN_units=True, LineShift=True,
                                  File=None, Format=None, OmegaGrid=None,
                                  WavenumberRange=None, WavenumberStep=None, WavenumberWing=None,
                                  WavenumberWingHW=None, WavenumberGrid=None):
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
    Components, SourceTables, Environment, OmegaRange, OmegaStep, OmegaWing,\
        IntensityThreshold, Format = \
        getDefaultValuesForXsect(Components, SourceTables, Environment, OmegaRange,
                                 OmegaStep, OmegaWing, IntensityThreshold, Format)
    # special for Doppler case: set OmegaStep to a smaller value
    if not OmegaStep:
        OmegaStep = 0.001

    # warn user about too large omega step
    if OmegaStep > 0.005:
        warn('Big wavenumber step: possible accuracy decline')

    # get uniform linespace for cross-section
    #number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    #Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        #Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = __FloatType__(296.)  # K
    pref = __FloatType__(1.)  # atm

    # actual temperature and pressure
    T = Environment['T']  # K
    p = Environment['p']  # atm

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
                ni = ISO[(M, I)][ISO_INDEX['abundance']]
            except KeyError:
                raise Exception('cannot find component M,I = %d,%d.' % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX['abundance']]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get line centers
        nline = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']

        # loop through line centers (single stream)
        for RowID in range(nline):

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]['data']['nu'][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]['data']['sw'][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]['data']['elower'][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['molec_id'][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['local_iso_id'][RowID]
            if LineShift:
                Shift0DB = LOCAL_TABLE_CACHE[TableName]['data']['delta_air'][RowID]
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
            LineIntensity = EnvironmentDependency_Intensity(LineIntensityDB, T, Tref, SigmaT, SigmaTref,
                                                            LowerStateEnergyDB, LineCenterDB)

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            # TODO: apply wing narrowing instead of filtering, this would be more appropriate
            if LineIntensity < IntensityThreshold:
                continue

            #   doppler broadening coefficient (GammaD)
            #GammaDDB = cSqrtLn2*LineCenterDB/cc*sqrt(2*cBolts*T/molecularMass(MoleculeNumberDB,IsoNumberDB))
            #GammaD = EnvironmentDependency_GammaD(GammaDDB,T,Tref)
            # print(GammaD)

            cMassMol = 1.66053873e-27
            #cSqrt2Ln2 = 1.1774100225
            fSqrtMass = sqrt(molecularMass(MoleculeNumberDB, IsoNumberDB))
            #fSqrtMass = sqrt(32831.2508809)
            cc_ = 2.99792458e8
            cBolts_ = 1.3806503e-23
            #cBolts_ = 1.3806488E-23
            GammaD = (cSqrt2Ln2/cc_)*sqrt(cBolts_/cMassMol) * \
                sqrt(T) * LineCenterDB/fSqrtMass

            #GammaD = 4.30140e-7*LineCenterDB*sqrt(T/molecularMass(MoleculeNumberDB,IsoNumberDB))

            # cc_ = 2.99792458e8 # 2.99792458e10 # 2.99792458e8
            # cBolts_ = 1.3806503e-23 #1.3806488E-16 # 1.380648813E-16 # 1.3806503e-23 # 1.3806488E-23
            #GammaD = sqrt(log(2))*LineCenterDB*sqrt(2*cBolts_*T/(cMassMol*molecularMass(MoleculeNumberDB,IsoNumberDB)*cc_**2))
            # print(GammaD)

            #   get final wing of the line according to GammaD, OmegaWingHW and OmegaWing
            # XXX min or max?
            OmegaWingF = max(OmegaWing, OmegaWingHW*GammaD)

            #   shift coefficient
            Shift0 = Shift0DB*p/pref

            # XXX other parameter (such as Delta0, Delta2, anuVC etc.) will be included in HTP version

            # PROFILE_VOIGT(sg0,GamD,Gam0,sg)
            #      sg0           : Unperturbed line position in cm-1 (Input).
            #      GamD          : Doppler HWHM in cm-1 (Input)
            #      Gam0          : Speed-averaged line-width in cm-1 (Input).
            #      sg            : Current WaveNumber of the Computation in cm-1 (Input).

            # XXX time?
            BoundIndexLower = bisect(Omegas, LineCenterDB-OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB+OmegaWingF)
            lineshape_vals = PROFILE_DOPPLER(
                LineCenterDB+Shift0, GammaD, Omegas[BoundIndexLower:BoundIndexUpper])
            #lineshape_vals = PROFILE_VOIGT(LineCenterDB,GammaD,cZero,Omegas[BoundIndexLower:BoundIndexUpper])[0]
            # Xsect[BoundIndexLower:BoundIndexUpper] += lineshape_vals # DEBUG
            Xsect[BoundIndexLower:BoundIndexUpper] += factor / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] * \
                LineIntensity * lineshape_vals

    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


# ---------------------------------------------------------------------------
# SHORTCUTS AND ALIASES FOR ABSORPTION COEFFICIENTS
# ---------------------------------------------------------------------------

absorptionCoefficient_Gauss = absorptionCoefficient_Doppler


def abscoef_HT(table=None, step=None, grid=None, env={'T': 296., 'p': 1.}, file=None):
    return absorptionCoefficient_HT(SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file)


def abscoef_Voigt(table=None, step=None, grid=None, env={'T': 296., 'p': 1.}, file=None):
    return absorptionCoefficient_Voigt(SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file)


def abscoef_Lorentz(table=None, step=None, grid=None, env={'T': 296., 'p': 1.}, file=None):
    return absorptionCoefficient_Lorentz(SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file)


def abscoef_Doppler(table=None, step=None, grid=None, env={'T': 296., 'p': 1.}, file=None):
    return absorptionCoefficient_Doppler(SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file)


abscoef_Gauss = abscoef_Doppler


# default
def abscoef(table=None, step=None, grid=None, env={'T': 296., 'p': 1.}, file=None):
    return absorptionCoefficient_Lorentz(SourceTables=table, OmegaStep=step, OmegaGrid=grid, Environment=env, File=file)

# ---------------------------------------------------------------------------


def transmittanceSpectrum(Omegas, AbsorptionCoefficient, Environment={'l': 100.},
                          File=None, Format='%e %e', Wavenumber=None):
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
    l = Environment['l']
    Xsect = exp(-AbsorptionCoefficient*l)
    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def absorptionSpectrum(Omegas, AbsorptionCoefficient, Environment={'l': 100.},
                       File=None, Format='%e %e', Wavenumber=None):
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
    l = Environment['l']
    Xsect = 1-exp(-AbsorptionCoefficient*l)
    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


def radianceSpectrum(Omegas, AbsorptionCoefficient, Environment={'l': 100., 'T': 296.},
                     File=None, Format='%e %e', Wavenumber=None):
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
    l = Environment['l']
    T = Environment['T']
    Alw = 1-exp(-AbsorptionCoefficient*l)
    LBBTw = 2*hh*cc**2*Omegas**3 / (exp(hh*cc*Omegas/(cBolts*T)) - 1) * 1.0E-7
    Xsect = Alw*LBBTw  # W/sr/cm**2/cm**-1
    if File:
        save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect


# GET X,Y FOR FINE PLOTTING OF A STICK SPECTRUM
def getStickXY(TableName):
    """
    Get X and Y for fine plotting of a stick spectrum.
    Usage: X,Y = getStickXY(TableName).
    """
    cent, intens = getColumns(TableName, ('nu', 'sw'))
    n = len(cent)
    cent_ = zeros(n*3)
    intens_ = zeros(n*3)
    for i in range(n):
        intens_[3*i] = 0
        intens_[3*i+1] = intens[i]
        intens_[3*i+2] = 0
        cent_[(3*i):(3*i+3)] = cent[i]
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
    f = open(filename, 'r')
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
    index_inner = abs(x) <= g/2
    index_outer = ~index_inner
    y = zeros(len(x))
    y[index_inner] = 1/g
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
    y[index_inner] = 1/g * (1 - abs(x[index_inner])/g)
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
    return sqrt(log(2))/(sqrt(pi)*g)*exp(-log(2)*(x/g)**2)

# dispersion slit function


def SLIT_DISPERSION(x, g):
    """
    Instrumental (slit) function.
    B(x) = γ/pi/(x**2+γ**2),
    where γ/2 is a lorentzian half-width at half-maximum.
    """
    g /= 2
    return g/pi/(x**2+g**2)

# cosinus slit function


def SLIT_COSINUS(x, g):
    return (cos(pi/g*x)+1)/(2*g)

# diffraction slit function


def SLIT_DIFFRACTION(x, g):
    """
    Instrumental (slit) function.
    """
    y = zeros(len(x))
    index_zero = x == 0
    index_nonzero = ~index_zero
    dk_ = pi/g
    x_ = dk_*x[index_nonzero]
    w_ = sin(x_)
    r_ = w_**2/x_**2
    y[index_zero] = 1
    y[index_nonzero] = r_/g
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
    dk_ = 2*pi/g
    x_ = dk_*x[index_nonzero]
    y[index_zero] = 1
    y[index_nonzero] = 2/g*sin(x_)/x_
    return y

# spectral convolution with an apparatus (slit) function


def convolveSpectrum(Omega, CrossSection, Resolution=0.1, AF_wing=10.,
                     SlitFunction=SLIT_RECTANGULAR, Wavenumber=None):
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
    step = Omega[1]-Omega[0]
    if step >= Resolution:
        raise Exception('step must be less than resolution')
    #x = arange(-AF_wing,AF_wing+step,step)
    x = arange_(-AF_wing, AF_wing+step, step)  # fix
    slit = SlitFunction(x, Resolution)
    # FIXING THE BUG: normalize slit function
    slit /= sum(slit)*step  # simple normalization
    left_bnd = len(slit)/2
    right_bnd = len(Omega) - len(slit)/2
    #CrossSectionLowRes = convolve(CrossSection,slit,mode='valid')*step
    CrossSectionLowRes = convolve(CrossSection, slit, mode='same')*step
    # return Omega[left_bnd:right_bnd],CrossSectionLowRes,left_bnd,right_bnd,slit
    return Omega[left_bnd:right_bnd], CrossSectionLowRes[left_bnd:right_bnd], left_bnd, right_bnd, slit

# DEBUG
# spectral convolution with an apparatus (slit) function


def convolveSpectrumSame(Omega, CrossSection, Resolution=0.1, AF_wing=10.,
                         SlitFunction=SLIT_RECTANGULAR):
    """
    Convolves cross section with a slit function with given parameters.
    """
    step = Omega[1]-Omega[0]
    x = arange(-AF_wing, AF_wing+step, step)
    slit = SlitFunction(x, Resolution)
    print('step=')
    print(step)
    print('x=')
    print(x)
    print('slitfunc=')
    print(SlitFunction)
    CrossSectionLowRes = convolve(CrossSection, slit, mode='same')*step
    return Omega, CrossSectionLowRes, None, None, slit

# DEBUG


def convolveSpectrumFull(Omega, CrossSection, Resolution=0.1, AF_wing=10., SlitFunction=SLIT_RECTANGULAR):
    """
    Convolves cross section with a slit function with given parameters.
    """
    step = Omega[1]-Omega[0]
    x = arange(-AF_wing, AF_wing+step, step)
    slit = SlitFunction(x, Resolution)
    print('step=')
    print(step)
    print('x=')
    print(x)
    print('slitfunc=')
    print(SlitFunction)
    CrossSectionLowRes = convolve(CrossSection, slit, mode='full')*step
    return Omega, CrossSectionLowRes, None, None

# ------------------------------------------------------------------
