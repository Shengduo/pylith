#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/2d/slipdir/TestFaultY.py
##
## @brief Test suite for testing sense of slip in 2-D for fault
## aligned with x-axis.

import numpy

from pylith.tests import run_pylith

from TestTri3 import TestTri3
from solution import SolnFaultY
from genspatialdb import GenDBFaultY

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class FaultYApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="faulty")
    return


class TestFaultY(TestTri3):
  """
  Test suite for testing sense of slip.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestTri3.setUp(self)
    run_pylith(FaultYApp, GenDBFaultY, nprocs=7)
    self.mesh['nvertices'] += 17
    self.outputRoot = "faulty"

    self.soln = SolnFaultY()
    return


  def calcStateVar(self, name, vertices, cells):
    """
    Calculate state variable.
    """
    ncells = self.mesh['ncells']
    pts = numpy.zeros( (ncells, 3), dtype=numpy.float64)
    if name == "total_strain":
      stateVar = self.soln.strain(pts)
    elif name == "stress":
      stateVar = self.soln.stress(pts)
    else:
      raise ValueError("Unknown state variable '%s'." % name)

    return stateVar


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestFaultY import TestFaultY as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
