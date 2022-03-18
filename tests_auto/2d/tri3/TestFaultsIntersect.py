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

## @file tests/2d/tri3/TestFaultsIntersect.py
##
## @brief Test suite for testing pylith with shear slip.

import numpy
from pylith.tests import run_pylith
from TestTri3 import TestTri3
from sliponefault_soln import AnalyticalSoln

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class FaultsIntersectApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="faultsintersect")
    return


class TestFaultsIntersect(TestTri3):
  """
  Test suite for testing pylith with shear faultsintersect for 2-D box.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestTri3.setUp(self)
    self.nverticesO = self.mesh['nvertices']

    # Fault x
    self.mesh['nvertices'] += 9
    self.faultMeshX = {'nvertices': 9,
                      'spaceDim': 2,
                      'ncells': 8,
                      'ncorners': 2}
    # Fault y
    self.mesh['nvertices'] += 1
    self.faultMeshY = {'nvertices': 3,
                      'spaceDim': 2,
                      'ncells': 2,
                      'ncorners': 2}

    run_pylith(FaultsIntersectApp)
    self.outputRoot = "faultsintersect"

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

    self.fault = "x"
    filename = "%s-faultx_info.h5" % self.outputRoot
    check_vertex_fields(self, filename, self.faultMeshX, fields)

    self.fault = "y"
    filename = "%s-faulty_info.h5" % self.outputRoot
    check_vertex_fields(self, filename, self.faultMeshY, fields)

    return


  def test_fault_data(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    from pylith.tests.Fault import check_vertex_fields
    fields = ["slip"]

    filename = "%s-faultx.h5" % self.outputRoot
    self.fault = "x"
    check_vertex_fields(self, filename, self.faultMeshX, fields)

    filename = "%s-faulty.h5" % self.outputRoot
    self.fault = "y"
    check_vertex_fields(self, filename, self.faultMeshY, fields)

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
    elif name == "stress" or name == "cauchy_stress":
      stateVar = self.soln.stress(pts)
    else:
      raise ValueError("Unknown state variable '%s'." % name)

    return stateVar


  def calcFaultField(self, name, vertices):
    """
    Calculate fault info.
    """

    if self.fault == "x":
      normalDir = (-1.0, 0.0)
      finalSlip = -2.0
      faultMesh = self.faultMeshX
    elif self.fault == "y":
      normalDir = (0.0, +1.0)
      finalSlip = 0.0
      faultMesh = self.faultMeshY
    slipTime = 0.0

    nvertices = faultMesh['nvertices']

    if name == "normal_dir":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = normalDir[0]
      field[0,:,1] = normalDir[1]

    elif name == "final_slip":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = finalSlip
      
    elif name == "slip_time":
      field = slipTime*numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)
      
    elif name == "slip":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = finalSlip

    elif name == "traction_change":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = 0.0
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestFaultsIntersect import TestFaultsIntersect as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
