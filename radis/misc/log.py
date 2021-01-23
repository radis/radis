# -*- coding: utf-8 -*-
"""
@author: Erwan

Logging and warning functions
"""


import sys
from time import gmtime, strftime
from warnings import warn

from termcolor import colored


def printlog(msg, logfile="log.txt", stdout=True):
    """Write a message to the logfile, adding date and time.

    Also print it to screen if stdout is True
    """
    with open(logfile, "a+") as f:
        f.write("{0} {1}\n".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()), msg))
    if stdout:
        print(msg)
    return


def warnlog(msg, logfile="log.txt", stdout=True):
    """Write a WARNING message to the logfile, adding date and time.

    Also print it to screen if stdout is True
    """
    with open(logfile, "a+") as f:
        f.write(
            "{0} WARNING. {1}\n".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()), msg)
        )
    if stdout:
        warn(msg)
    return


def printwarn(msg, verbose=True, warnings=True):
    """Send msg in the console if verbose is True, and as warning if warnings
    is True.

    Warnings is added at the end of the process but you may miss it if
    interrupted.
    """

    if verbose:
        print(colored("\n" + msg + "\n", "red"))
        sys.stdout.flush()
    if warnings:
        warn(msg)
    return
