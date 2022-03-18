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

## @file tests/3d/tet4/TestSlipOneFault.py
##
## @brief Test suite for testing pylith with shear slip.

import numpy

from pylith.tests import run_pylith

from TestTet4 import TestTet4
from sliponefault_soln import AnalyticalSoln

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class SlipOneFaultApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="sliponefault")
    return


class TestSlipOneFault(TestTet4):
  """
  Test suite for testing pylith with shear sliponefault for 2-D box.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestTet4.setUp(self)
    self.nverticesO = self.mesh['nvertices']
    self.mesh['nvertices'] += 50
    self.faultMesh = {'nvertices': 50,
                      'spaceDim': 3,
                      'ncells': 72,
                      'ncorners': 3}
    run_pylith(SlipOneFaultApp)
    self.outputRoot = "sliponefault"

    self.soln = AnalyticalSoln()
    return


  def test_fault_info(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault_info.h5" % self.outputRoot
    fields = ["normal_dir", "final_slip", "slip_time"]

    from pylith.tests.Fault import check_vertex_fields
    check_vertex_fields(self, filename, self.faultMesh, fields)

    return


  def test_fault_data(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault.h5" % self.outputRoot
    fields = ["slip"]

    from pylith.tests.Fault import check_vertex_fields
    check_vertex_fields(self, filename, self.faultMesh, fields)

    return


  def test_points_data(self):
    """
    Check points information.
    """
    if not self.checkResults:
      return

    filename = "%s-points.h5" % self.outputRoot
    fields = ["displacement"]
    stations = ["ZZ.AAAA", "ZZ.BBBB", "ZZ.CCCC", "ZZ.DDDD"]

    from pylith.tests.SolutionPoints import check_displacements
    check_displacements(self, filename, npoints=len(stations), spaceDim=2)

    from pylith.tests.SolutionPoints import check_stations
    check_stations(self, filename, stations)

    return


  def calcDisplacements(self, vertices):
    """
    Calculate displacement field given coordinates of vertices.
    """
    return self.soln.displacement(vertices, self.nverticesO)


  def calcDisplacementPoints(self, vertices):
    """
    Calculate displacement field given coordinates of vertices.
    """
    return self.soln.displacement(vertices)


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
    dim = 3

    normalDir = (+1.0, 0.0, 0.0)
    finalSlip = -2.0
    slipTime = 0.0

    nvertices = self.faultMesh['nvertices']

    if name == "normal_dir":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = normalDir[0]
      field[0,:,1] = normalDir[1]
      field[0,:,2] = normalDir[2]

    elif name == "final_slip":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = finalSlip
      
    elif name == "slip_time":
      field = slipTime*numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)
      
    elif name == "slip":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = finalSlip

    elif name == "traction_change":
      field = numpy.zeros( (1, nvertices, dim), dtype=numpy.float64)
      field[0,:,0] = 0.0
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestSlipOneFault import TestSlipOneFault as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
