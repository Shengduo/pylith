// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application quadratureapp.

#include "QuadratureData1Din2DLinear.hh"

const int pylith::feassemble::QuadratureData1Din2DLinear::_numVertices = 2;

const int pylith::feassemble::QuadratureData1Din2DLinear::_spaceDim = 2;

const int pylith::feassemble::QuadratureData1Din2DLinear::_numCells = 1;

const int pylith::feassemble::QuadratureData1Din2DLinear::_cellDim = 1;

const int pylith::feassemble::QuadratureData1Din2DLinear::_numBasis = 2;

const int pylith::feassemble::QuadratureData1Din2DLinear::_numQuadPts = 1;

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_vertices[] = {
 -2.00000000e-01, -5.00000000e-01,
  7.00000000e-01,  3.00000000e-01,
};

const int pylith::feassemble::QuadratureData1Din2DLinear::_cells[] = {
       0,       1,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_verticesRef[] = {
 -1.00000000e+00,
  1.00000000e+00,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_quadPtsRef[] = {
  0.00000000e+00,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_quadWts[] = {
  2.00000000e+00,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_quadPts[] = {
  2.50000000e-01, -1.00000000e-01,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_basis[] = {
  5.00000000e-01,
  5.00000000e-01,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_basisDerivRef[] = {
 -5.00000000e-01,
  5.00000000e-01,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_basisDeriv[] = {
 -1.11111111e+00, -1.25000000e+00,
  1.11111111e+00,  1.25000000e+00,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_jacobian[] = {
  4.50000000e-01,
  4.00000000e-01,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_jacobianDet[] = {
  6.02079729e-01,
};

const PylithScalar pylith::feassemble::QuadratureData1Din2DLinear::_jacobianInv[] = {
  2.22222222e+00,  2.50000000e+00,
};

pylith::feassemble::QuadratureData1Din2DLinear::QuadratureData1Din2DLinear(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  quadPtsRef = const_cast<PylithScalar*>(_quadPtsRef);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  quadPts = const_cast<PylithScalar*>(_quadPts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);
  basisDeriv = const_cast<PylithScalar*>(_basisDeriv);
  jacobian = const_cast<PylithScalar*>(_jacobian);
  jacobianDet = const_cast<PylithScalar*>(_jacobianDet);
  jacobianInv = const_cast<PylithScalar*>(_jacobianInv);
} // constructor

pylith::feassemble::QuadratureData1Din2DLinear::~QuadratureData1Din2DLinear(void)
{}


// End of file
