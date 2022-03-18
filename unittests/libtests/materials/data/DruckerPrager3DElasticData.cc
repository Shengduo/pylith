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
// This file was generated from python application druckerprager3delastic.

#include "DruckerPrager3DElasticData.hh"

const int pylith::materials::DruckerPrager3DElasticData::_dimension = 3;

const int pylith::materials::DruckerPrager3DElasticData::_numLocs = 2;

const int pylith::materials::DruckerPrager3DElasticData::_numProperties = 6;

const int pylith::materials::DruckerPrager3DElasticData::_numStateVars = 1;

const int pylith::materials::DruckerPrager3DElasticData::_numDBProperties = 6;

const int pylith::materials::DruckerPrager3DElasticData::_numDBStateVars = 6;

const int pylith::materials::DruckerPrager3DElasticData::_numPropsQuadPt = 6;

const int pylith::materials::DruckerPrager3DElasticData::_numVarsQuadPt = 6;

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_lengthScale =   1.00000000e+03;

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_timeScale =   1.00000000e+00;

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_pressureScale =   2.25000000e+10;

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_densityScale =   2.25000000e+04;

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_dtStableImplicit =   1.00000000e+10;

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_dtStableExplicit =   1.92450090e-01;

const int pylith::materials::DruckerPrager3DElasticData::_numPropertyValues[] = {
1,
1,
1,
1,
1,
1,
};

const int pylith::materials::DruckerPrager3DElasticData::_numStateVarValues[] = {
6,
};

const char* pylith::materials::DruckerPrager3DElasticData::_dbPropertyValues[] = {
"density",
"vs",
"vp",
"friction-angle",
"cohesion",
"dilatation-angle",
};

const char* pylith::materials::DruckerPrager3DElasticData::_dbStateVarValues[] = {
"plastic-strain-xx",
"plastic-strain-yy",
"plastic-strain-zz",
"plastic-strain-xy",
"plastic-strain-yz",
"plastic-strain-xz",
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_dbProperties[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  5.23598776e-01,
  3.00000000e+05,
  3.49065850e-01,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
  4.36332313e-01,
  1.00000000e+05,
  4.36332313e-01,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_dbStateVars[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_properties[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  2.30940108e-01,
  3.60000000e+05,
  1.48583084e-01,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
  1.89338478e-01,
  1.21811303e+05,
  1.89338478e-01,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_stateVars[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_propertiesNondim[] = {
  1.11111111e-01,
  1.00000000e+00,
  1.00000000e+00,
  2.30940108e-01,
  1.60000000e-05,
  1.48583084e-01,
  8.88888889e-02,
  1.28000000e-01,
  1.28000000e-01,
  1.89338478e-01,
  5.41383567e-06,
  1.89338478e-01,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_stateVarsNondim[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_strain[] = {
 -1.10000000e-04,
 -1.20000000e-04,
 -1.30000000e-04,
  1.40000000e-04,
  1.50000000e-04,
  1.60000000e-04,
  4.10000000e-04,
  4.20000000e-04,
  4.30000000e-04,
  4.40000000e-04,
  4.50000000e-04,
  4.60000000e-04,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_stress[] = {
 -4.85790000e+07,
 -4.94780000e+07,
 -5.03770000e+07,
 -8.97600000e+06,
 -8.97500000e+06,
 -8.97400000e+06,
 -2.82900000e+06,
 -2.82800000e+06,
 -2.82700000e+06,
 -1.09800000e+06,
 -1.09700000e+06,
 -1.09600000e+06,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_elasticConsts[] = {
  6.75000000e+10,
  2.25000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.25000000e+10,
  6.75000000e+10,
  2.25000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.25000000e+10,
  2.25000000e+10,
  6.75000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  4.50000000e+10,
  8.64000000e+09,
  2.88000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.88000000e+09,
  8.64000000e+09,
  2.88000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  2.88000000e+09,
  2.88000000e+09,
  8.64000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  5.76000000e+09,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_initialStress[] = {
  2.10000000e+04,
  2.20000000e+04,
  2.30000000e+04,
  2.40000000e+04,
  2.50000000e+04,
  2.60000000e+04,
  5.10000000e+04,
  5.20000000e+04,
  5.30000000e+04,
  5.40000000e+04,
  5.50000000e+04,
  5.60000000e+04,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_initialStrain[] = {
  3.10000000e-04,
  3.20000000e-04,
  3.30000000e-04,
  3.40000000e-04,
  3.50000000e-04,
  3.60000000e-04,
  6.10000000e-04,
  6.20000000e-04,
  6.30000000e-04,
  6.40000000e-04,
  6.50000000e-04,
  6.60000000e-04,
};

const PylithScalar pylith::materials::DruckerPrager3DElasticData::_stateVarsUpdated[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

pylith::materials::DruckerPrager3DElasticData::DruckerPrager3DElasticData(void)
{ // constructor
  dimension = _dimension;
  numLocs = _numLocs;
  numProperties = _numProperties;
  numStateVars = _numStateVars;
  numDBProperties = _numDBProperties;
  numDBStateVars = _numDBStateVars;
  numPropsQuadPt = _numPropsQuadPt;
  numVarsQuadPt = _numVarsQuadPt;
  lengthScale = _lengthScale;
  timeScale = _timeScale;
  pressureScale = _pressureScale;
  densityScale = _densityScale;
  dtStableImplicit = _dtStableImplicit;
  dtStableExplicit = _dtStableExplicit;
  numPropertyValues = const_cast<int*>(_numPropertyValues);
  numStateVarValues = const_cast<int*>(_numStateVarValues);
  dbPropertyValues = const_cast<char**>(_dbPropertyValues);
  dbStateVarValues = const_cast<char**>(_dbStateVarValues);
  dbProperties = const_cast<PylithScalar*>(_dbProperties);
  dbStateVars = const_cast<PylithScalar*>(_dbStateVars);
  properties = const_cast<PylithScalar*>(_properties);
  stateVars = const_cast<PylithScalar*>(_stateVars);
  propertiesNondim = const_cast<PylithScalar*>(_propertiesNondim);
  stateVarsNondim = const_cast<PylithScalar*>(_stateVarsNondim);
  density = const_cast<PylithScalar*>(_density);
  strain = const_cast<PylithScalar*>(_strain);
  stress = const_cast<PylithScalar*>(_stress);
  elasticConsts = const_cast<PylithScalar*>(_elasticConsts);
  initialStress = const_cast<PylithScalar*>(_initialStress);
  initialStrain = const_cast<PylithScalar*>(_initialStrain);
  stateVarsUpdated = const_cast<PylithScalar*>(_stateVarsUpdated);
} // constructor

pylith::materials::DruckerPrager3DElasticData::~DruckerPrager3DElasticData(void)
{}


// End of file
