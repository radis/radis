"""
This file exists to ensure that the project root is added to sys.path during test collection.
This fixes 'ModuleNotFoundError: No module named radis' errors that can occur in some
environments (including local development and CI) even when the package is installed in editable mode.
"""
