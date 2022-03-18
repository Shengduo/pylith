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

## @file tests/2d/tri3/TestShearDispFriction.py
##
## @brief Test suite for testing pylith with 2-D shear motion with no
## fault slip.

import numpy

from pylith.tests import run_pylith

from TestTri3 import TestTri3
from sheardisp_soln import AnalyticalSoln
from sheardisp_gendb import GenerateDB

# Local version of PyLithApp
from pylith.apps.PyLithApp import PyLithApp
class ShearApp(PyLithApp):
  def __init__(self):
    PyLithApp.__init__(self, name="sheardispfriction")
    return


class TestShearDispFriction(TestTri3):
  """
  Test suite for testing pylith with 2-D shear extension.
  """

  def setUp(self):
    """
    Setup for test.
    """
    TestTri3.setUp(self)
    self.nverticesO = self.mesh['nvertices']
    self.mesh['nvertices'] += 1
    self.faultMesh = {'nvertices': 3,
                      'spaceDim': 2,
                      'ncells': 2,
                      'ncorners': 2}
    run_pylith(ShearApp, GenerateDB, nprocs=3)
    self.outputRoot = "sheardispfriction"

    self.soln = AnalyticalSoln()
    return


  def test_fault_info(self):
    """
    Check fault information.
    """
    if not self.checkResults:
      return

    filename = "%s-fault_info.h5" % self.outputRoot
    fields = ["strike_dir", "normal_dir", "traction_initial","friction_coefficient","cohesion"]

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


  def calcDisplacements(self, vertices):
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

    strikeDir = (-1.0, 0.0)
    normalDir = (0.0, 1.0)
    initialTraction = (0.0, -100.0e+6)
    frictionCoefficient = 1.0

    nvertices = self.faultMesh['nvertices']

    if name == "strike_dir":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = strikeDir[0]
      field[0,:,1] = strikeDir[1]

    elif name == "normal_dir":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = normalDir[0]
      field[0,:,1] = normalDir[1]

    elif name == "traction_initial":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = initialTraction[0]
      field[0,:,1] = initialTraction[1]

    elif name == "friction_coefficient":
      field = numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)
      field[0,:,0] = frictionCoefficient

    elif name == "cohesion":
      field = numpy.zeros( (1, nvertices, 1), dtype=numpy.float64)

    elif name == "slip":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)

    elif name == "traction":
      field = numpy.zeros( (1, nvertices, 2), dtype=numpy.float64)
      field[0,:,0] = 1.0e+6
      field[0,:,1] = -1.0e+7
      
    else:
      raise ValueError("Unknown fault field '%s'." % name)

    # Mask clamped vertices
    maskX = numpy.bitwise_and(vertices[:,0] >= -2.0001e+3, vertices[:,0] <= 1.0)
    maskY = numpy.fabs(vertices[:,1]) <= 0.01e+3
    mask = numpy.bitwise_and(maskX,maskY)
    field[:,~mask,:] = 0

    return field


# ----------------------------------------------------------------------
if __name__ == '__main__':
  import unittest
  from TestShearDispFriction import TestShearDispFriction as Tester

  suite = unittest.TestSuite()
  suite.addTest(unittest.makeSuite(Tester))
  unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 
