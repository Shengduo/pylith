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

/* Original mesh
 *
 * Cells are 0-1, vertices are 2-6.
 *
 * 2   3,4,5  6
 *
 *     ^^^^^ Face in x-y plane
 *
 * After adding cohesive elements
 *
 * Cells are 0-1,2, vertices are 3-10.
 *
 * 3   8,9,10  4,5,6    7
 *
 *     ^^^^^^^^^^^^ Cohesive element in x-y plane.
 */

#include "CohesiveDataTet4f.hh"

const int pylith::faults::CohesiveDataTet4f::_numVertices = 8;

const int pylith::faults::CohesiveDataTet4f::_spaceDim = 3;

const int pylith::faults::CohesiveDataTet4f::_numCells = 3;

const int pylith::faults::CohesiveDataTet4f::_cellDim = 3;

const int pylith::faults::CohesiveDataTet4f::_numCorners[3] = {
  4,
  4,
  6
};

const int pylith::faults::CohesiveDataTet4f::_materialIds[3] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTet4f::_numGroups = 2;

const int pylith::faults::CohesiveDataTet4f::_groupSizes[2] = {
  5+4+1, 6+6+2 // vertices+edges+faces
};

const char* pylith::faults::CohesiveDataTet4f::_groupNames[2] = {
  "output", "fault"
};

const char* pylith::faults::CohesiveDataTet4f::_groupTypes[2] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTet4f::_filename = "data/tet4f.mesh";

pylith::faults::CohesiveDataTet4f::CohesiveDataTet4f(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  numCorners = const_cast<int*>(_numCorners);
  materialIds = const_cast<int*>(_materialIds);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
} // constructor

pylith::faults::CohesiveDataTet4f::~CohesiveDataTet4f(void)
{}


// End of file
