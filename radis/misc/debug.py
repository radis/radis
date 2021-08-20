# -*- coding: utf-8 -*-
"""Created on Sun Jan  3 17:52:04 2016.

@author: Erwan

Debug functions

-------------------------------------------------------------------------------
"""


# %%
# ==============================================================================
# Debug functions
# ==============================================================================


def printdbg(*args, **kwargs):
    """Function that prints only in debug mode. change this at runtime with.

    >>> radis.debug = True

    Examples
    --------

    Embed this print in a if __debug__ statement::

        if __debug__: printdbg(...)

    so that printdbg are removed by the Python preprocessor when running in
    optimize mode::

        python -O *.py

    See the ``"DEBUG_MODE"`` key in :py:attr:`radis.config`
    """

    import radis

    if radis.config["DEBUG_MODE"]:
        print("DEBUG:", *args, **kwargs)


def export(var=locals()):
    """Export local variables. Useful for debugging.

    Debugging inside a function may be tedious because you can't access the
    local variables. One of the option is to use the ipython magic::

        %debug

    Or the pdb equivalent::

        import pdb
        pdb.pm()

    Another option is to insert this export() call in the troubled function,
    before the exception occurs.


    Examples
    --------

        debug_export(locals())

    Note: you can also use  'globals().update(locals())' directly in your
    function to debug

    Note
    ----

    - seems not to work for functions nested in functions
    - 01/05 : doesn't seem to work at all.. @Erwan
    """

    globals().update(var)

    return


#
# def get_locals(func):
#    ''' Initially from :
#    http://stackoverflow.com/questions/20577806/keep-function-namespace-alive-for-debugging-in-ipython
#
#    '''
#
#    def wrap(*args, **kw):
#        sys.settrace(tracefunc)
#        try:
#            res = func(*args, **kw)
#        finally:
#            sys.settrace(None)
#        return res
#
#    def tracefunc(frame, event, arg):
#        if event == "return":
#            if frame.f_code is func.func_code:
#                wrap.last_res = frame.f_locals
#        return tracefunc
#
#    return wrap
