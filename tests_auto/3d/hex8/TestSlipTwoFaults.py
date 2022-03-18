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

## @file tests/3d/hex8/TestSlipTwoFaults.py
##
## @brief Test suite for testing pylith with shear slip.

import numpy

from pylith.tests import run_pylith

from TestHex8 import TestHex8
from sliptwofaults_soln import AnalyticalSoln

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class SlipTwoFaultsApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="sliptwofaults")
    return


class TestSlipTwoFaults(TestHex8):
  """
  Test suite for testing pylith with shear slip on two faults.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestHex8.setUp(self)
    self.nverticesO = self.mesh['nvertices']
    self.mesh['nvertices'] += 2*55
    self.faultMesh = {'nvertices': 55,
                      'spaceDim': 3,
                      'ncells': 40,
                      'ncorners': 4}
    run_pylith(SlipTwoFaultsApp)
    self.outputRoot = "sliptwofaults"

    self.soln = AnalyticalSoln()
    return


  def test_fault_info(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    from pylith.tests.Fault import check_vertex_fields
    fields = ["normal_dir", "final_slip", "slip_time"]

    self.fault = 1
    filename = "%s-fault1_info.h5" % self.outputRoot
    check_vertex_fields(self, filename, self.faultMesh, fields)

    self.fault = 2
    filename = "%s-fault2_info.h5" % self.outputRoot
    check_vertex_fields(self, filename, self.faultMesh, fields)

    return


  def test_fault_data(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    from pylith.tests.Fault import check_vertex_fields
    fields = ["slip"]

    filename = "%s-fault1.h5" % self.outputRoot
    self.fault = 1
    check_vertex_fields(self, filename, self.faultMesh, fields)

    filename = "%s-fault2.h5" % self.outputRoot
    self.fault = 2
    check_vertex_fields(self, filename, self.faultMesh, fields)

    return


  def calcDisplacements(self, vertices):
    """
    Calculate displacement field given coordinates of vertices.
    """
    return self.soln.displacement(vertices, self.nverticesO)


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


  def calcFaultField(self, name, vertices):
    """
    Calculate fault info.
    """

    normalDir = (-1.0, 0.0, 0.0)
    finalSlip = -2.0
    slipTime = 0.0
    dim = 3

    nvertices = self.faultMesh['nvertices']

    if name == "normal_dir":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = normalDir[0]
      field[0,:,1] = normalDir[1]
      field[0,:,2] = normalDir[2]

    elif name == "final_slip":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = finalSlip

      if self.fault == 2:
        field *= -1
      
    elif name == "slip_time":
      field = slipTime*numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)
      
    elif name == "slip":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = finalSlip

      if self.fault == 2:
        field *= -1

    elif name == "traction_change":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = 0.0
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestSlipTwoFaults import TestSlipTwoFaults as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
