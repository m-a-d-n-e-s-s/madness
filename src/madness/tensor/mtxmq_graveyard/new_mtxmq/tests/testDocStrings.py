import unittest
import doctest
import sys
import os

def grep(filename, string):
    """
    Input:
        A filename to read and a string to search for.

    Output:
        True if the string appears in the file
        False otherwise
    """
    for line in open(filename):
        if string in line:
            return True
    return False

def load_tests(loader, tests, ignore):
    """
    Import files with doctests into the unittesting framework.
    """
    # Find the base of the project 
    base = sys.path[0]

    # Walk all files in this project
    for dirpath, dirnames, filenames in os.walk(base):
        for fn in filenames:
            # Check python files for the doctest prompt
            if fn[-3:] == ".py" and grep(os.path.join(dirpath,fn), '>'*3):
                # remove the project name and .py
                module = os.path.join(dirpath[len(base):], fn[:-3])
                # use '.' to make it a string that represents a module
                module = module.replace(os.sep,'.').strip('.')
                # add the module as a test
                tests.addTests(doctest.DocTestSuite(module))

    return tests
