// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesiveDyn.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveDyn.hh" // USES FaultCohesiveDyn
#include "pylith/faults/TractionPerturbation.hh" // USES TractionPerturbation

#include "data/CohesiveDynData.hh" // USES CohesiveDynData

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/VisitorSubMesh.hh" // USES SubMeshIS
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/friction/StaticFriction.hh" // USES StaticFriction

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/spatialdb/SimpleDB.hh" // USES SimpleDB
#include "spatialdata/spatialdb/SimpleIOAscii.hh" // USES SimpleIOAscii
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES runtime_error

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::faults::TestFaultCohesiveDyn );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveDyn::setUp(void)
{ // setUp
  PYLITH_METHOD_BEGIN;

  _data = 0;
  _quadrature = new feassemble::Quadrature();CPPUNIT_ASSERT(_quadrature);
  _tractionPerturbation = 0;
  _dbInitialTract = 0;
  _friction = 0;
  _dbFriction = 0;

  PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::faults::TestFaultCohesiveDyn::tearDown(void)
{ // tearDown
  PYLITH_METHOD_BEGIN;

  delete _data; _data = 0;
  delete _quadrature; _quadrature = 0;
  delete _tractionPerturbation; _tractionPerturbation = 0;
  delete _dbInitialTract; _dbInitialTract = 0;
  delete _friction; _friction = 0;
  delete _dbFriction; _dbFriction = 0;

  PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveDyn::testConstructor(void)
{ // testConstructor
  PYLITH_METHOD_BEGIN;

  FaultCohesiveDyn fault;

  PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test tractionPerturbation().
void
pylith::faults::TestFaultCohesiveDyn::testTractionPerturbation(void)
{ // testTractionPerturbation
  PYLITH_METHOD_BEGIN;

  FaultCohesiveDyn fault;

  const std::string& label = "test database";
  TractionPerturbation tract;
  tract.label(label.c_str());
  fault.tractionPerturbation(&tract);
  CPPUNIT_ASSERT(fault._tractionPerturbation);

  PYLITH_METHOD_END;
} // testTractionPerturbation

// ----------------------------------------------------------------------
// Test zeroTolerance().
void
pylith::faults::TestFaultCohesiveDyn::testZeroTolerance(void)
{ // testZeroTolerance
  PYLITH_METHOD_BEGIN;

  FaultCohesiveDyn fault;

  CPPUNIT_ASSERT_EQUAL(PylithScalar(1.0e-10), fault._zeroTolerance); // default

  const PylithScalar value = 1.0e-20;
  fault.zeroTolerance(value);
  CPPUNIT_ASSERT_EQUAL(value, fault._zeroTolerance);

  PYLITH_METHOD_END;
} // zeroTolerance

// ----------------------------------------------------------------------
// Test openFreeSurf().
void
pylith::faults::TestFaultCohesiveDyn::testOpenFreeSurf(void)
{ // testOpenFreeSurf
  PYLITH_METHOD_BEGIN;

  FaultCohesiveDyn fault;

  CPPUNIT_ASSERT_EQUAL(true, fault._openFreeSurf); // default

  const bool value = false;
  fault.openFreeSurf(value);
  CPPUNIT_ASSERT_EQUAL(value, fault._openFreeSurf);
 } // testOpenFreeSurf

// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveDyn::testInitialize(void)
{ // testInitialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(_data);
  PetscErrorCode err;

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);

  topology::SubMeshIS subpointIS(*fault._faultMesh);
  const PetscInt numPoints = subpointIS.size();
  const PetscInt* points = subpointIS.points();CPPUNIT_ASSERT(points);

  PetscDM dmMesh = fault._faultMesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  const PetscInt vStart = verticesStratum.begin();
  const PetscInt vEnd = verticesStratum.end();

  for(PetscInt v = vStart; v < vEnd; ++v) {
    PetscInt faultPoint;

    err = PetscFindInt(_data->negativeVertices[v-vStart], numPoints, points, &faultPoint);PYLITH_CHECK_ERROR(err);
    CPPUNIT_ASSERT(faultPoint >= 0);
    CPPUNIT_ASSERT_EQUAL(faultPoint, v);
  } // for
  CPPUNIT_ASSERT_EQUAL(_data->numConstraintEdges, vEnd-vStart);

  // Check orientation
  //fault._fields->get("orientation").view("ORIENTATION"); // DEBUGGING
  topology::VecVisitorMesh orientationVisitor(fault._fields->get("orientation"));
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  const int spaceDim = _data->spaceDim;
  const int orientationSize = spaceDim*spaceDim;
  for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
    const PetscInt off = orientationVisitor.sectionOffset(v);
    CPPUNIT_ASSERT_EQUAL(orientationSize, orientationVisitor.sectionDof(v));

    const PylithScalar tolerance = 1.0e-06;
    for(PetscInt d = 0; d < orientationSize; ++d) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->orientation[iVertex*orientationSize+d], orientationArray[off+d], tolerance);
    } // for
  } // for

  // Prescribed traction perturbation
  if (fault._tractionPerturbation) {
    // :KLUDGE: Only check initial value
    topology::VecVisitorMesh tractionVisitor(fault.vertexField("traction_initial_value"));
    const PetscScalar* tractionArray = tractionVisitor.localArray();CPPUNIT_ASSERT(tractionArray);
    const PylithScalar tractionScale = _data->pressureScale;

    for(PetscInt v = vStart, iVertex = 0; v < vEnd; ++v, ++iVertex) {
      const PetscInt off = tractionVisitor.sectionOffset(v);
      CPPUNIT_ASSERT_EQUAL(spaceDim, tractionVisitor.sectionDof(v));

      const PylithScalar tolerance = 1.0e-06;
      for(PetscInt d = 0; d < spaceDim; ++d) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(_data->initialTractions[iVertex * spaceDim + d], tractionArray[off+d]*_data->pressureScale, tolerance);
      } // for
    } // for
  } // if

  PYLITH_METHOD_END;
} // testInitialize

// ----------------------------------------------------------------------
// Test integrateResidual() for sticking case.
void
pylith::faults::TestFaultCohesiveDyn::testIntegrateResidualStick(void)
{ // testIntegrateResidualStick
  PYLITH_METHOD_BEGIN;

  assert(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution());
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrStick);

  const int spaceDim = _data->spaceDim;

  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01 / _data->timeScale;
  fault.timeStep(dt);

  const topology::Field& residual = fields.get("residual");
  fault.integrateResidual(residual, t, &fields);
  
  residual.view("RESIDUAL"); // DEBUGGING
  { // Check residual values
    PetscInt pStart=0, pEnd=0;
    PetscErrorCode err = PetscSectionGetChart(residual.localSection(), &pStart, &pEnd);CPPUNIT_ASSERT(!err);
    topology::VecVisitorMesh residualVisitor(residual);
    const PetscScalar* residualArray = residualVisitor.localArray();CPPUNIT_ASSERT(residualArray);

    const PylithScalar* residualE = _data->residualStickE;CPPUNIT_ASSERT(residualE);
    const int fiberDimE = spaceDim; // number of values per point
    const PylithScalar tolerance = 1.0e-06;

    for (PetscInt p=pStart, iPoint=0; p < pEnd; ++p) {
      if (residualVisitor.sectionDof(p) > 0) {
	const PetscInt off = residualVisitor.sectionOffset(p);
	CPPUNIT_ASSERT_EQUAL(fiberDimE, residualVisitor.sectionDof(p));
	for (PetscInt iDim=0; iDim < fiberDimE; ++iDim) {
	  const PylithScalar valE = residualE[iPoint*fiberDimE+iDim];
	  if (fabs(valE) > tolerance) {
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+iDim]/valE, tolerance);
	  } else {
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, residualArray[off+iDim], tolerance);
	  } // if/else
	} // for
	++iPoint;
      } // if
    } // for
  } // Check residual values

  PYLITH_METHOD_END;
} // testIntegrateResidualStick

// ----------------------------------------------------------------------
// Test integrateResidual() for slipping case.
void
pylith::faults::TestFaultCohesiveDyn::testIntegrateResidualSlip(void)
{ // testIntegrateResidualSlip
  PYLITH_METHOD_BEGIN;

  assert(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution());
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrSlip);

  const int spaceDim = _data->spaceDim;

  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01 / _data->timeScale;
  fault.timeStep(dt);

  const topology::Field& residual = fields.get("residual");
  fault.integrateResidual(residual, t, &fields);

  residual.view("RESIDUAL"); // DEBUGGING
  { // Check residual values
    PetscInt pStart=0, pEnd=0;
    PetscErrorCode err = PetscSectionGetChart(residual.localSection(), &pStart, &pEnd);CPPUNIT_ASSERT(!err);
    topology::VecVisitorMesh residualVisitor(residual);
    const PetscScalar* residualArray = residualVisitor.localArray();CPPUNIT_ASSERT(residualArray);

    const PylithScalar* residualE = _data->residualSlipE;CPPUNIT_ASSERT(residualE);
    const int fiberDimE = spaceDim; // number of values per point
    const PylithScalar tolerance = 1.0e-06;

    for (PetscInt p=pStart, iPoint=0; p < pEnd; ++p) {
      if (residualVisitor.sectionDof(p) > 0) {
	const PetscInt off = residualVisitor.sectionOffset(p);
	CPPUNIT_ASSERT_EQUAL(fiberDimE, residualVisitor.sectionDof(p));
	for (PetscInt iDim=0; iDim < fiberDimE; ++iDim) {
	  const PylithScalar valE = residualE[iPoint*fiberDimE+iDim];
	  if (fabs(valE) > tolerance) {
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+iDim]/valE, tolerance);
	  } else {
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, residualArray[off+iDim], tolerance);
	  } // if/else
	} // for
	++iPoint;
      } // if
    } // for
  } // Check residual values

  PYLITH_METHOD_END;
} // testIntegrateResidualSlip

// ----------------------------------------------------------------------
// Test integrateResidual() for opening case.
void
pylith::faults::TestFaultCohesiveDyn::testIntegrateResidualOpen(void)
{ // testIntegrateResidualOpen
  PYLITH_METHOD_BEGIN;

  assert(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution());
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrOpen);

  const int spaceDim = _data->spaceDim;

  const PylithScalar t = 2.134 / _data->timeScale;
  const PylithScalar dt = 0.01 / _data->timeScale;
  fault.timeStep(dt);

  const topology::Field& residual = fields.get("residual");
  fault.integrateResidual(residual, t, &fields);

  topology::Field& solution = fields.solution();
  const topology::Field& dispIncrAdj = fields.get("dispIncr adjust");
  solution += dispIncrAdj;

  fault.updateStateVars(t, &fields);

  residual.view("RESIDUAL"); // DEBUGGING
  { // Check residual values
    PetscInt pStart=0, pEnd=0;
    PetscErrorCode err = PetscSectionGetChart(residual.localSection(), &pStart, &pEnd);CPPUNIT_ASSERT(!err);
    topology::VecVisitorMesh residualVisitor(residual);
    const PetscScalar* residualArray = residualVisitor.localArray();CPPUNIT_ASSERT(residualArray);

    const PylithScalar* residualE = _data->residualOpenE;CPPUNIT_ASSERT(residualE);
    const int fiberDimE = spaceDim; // number of values per point
    const PylithScalar tolerance = 1.0e-06;

    for (PetscInt p=pStart, iPoint=0; p < pEnd; ++p) {
      if (residualVisitor.sectionDof(p) > 0) {
	const PetscInt off = residualVisitor.sectionOffset(p);
	CPPUNIT_ASSERT_EQUAL(fiberDimE, residualVisitor.sectionDof(p));
	for (PetscInt iDim=0; iDim < fiberDimE; ++iDim) {
	  const PylithScalar valE = residualE[iPoint*fiberDimE+iDim];
	  if (fabs(valE) > tolerance) {
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, residualArray[off+iDim]/valE, tolerance);
	  } else {
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(valE, residualArray[off+iDim], tolerance);
	  } // if/else
	} // for
	++iPoint;
      } // if
    } // for
  } // Check residual values

  PYLITH_METHOD_END;
} // testIntegrateResidualOpen

// ----------------------------------------------------------------------
// Test updateStateVars().
void
pylith::faults::TestFaultCohesiveDyn::testUpdateStateVars(void)
{ // testUpdateStateVars
  PYLITH_METHOD_BEGIN;

  assert(_data);

  topology::Mesh mesh;
  FaultCohesiveDyn fault;
  topology::SolutionFields fields(mesh);
  _initialize(&mesh, &fault, &fields);
  topology::Jacobian jacobian(fields.solution());
  _setFieldsJacobian(&mesh, &fault, &fields, &jacobian, _data->fieldIncrSlip);

  const int spaceDim = _data->spaceDim;

  const PylithScalar t = 2.134;
  const PylithScalar dt = 0.01;
  fault.timeStep(dt);
  fault.updateStateVars(t, &fields);

  // :TODO: Need to verify that fault constitutive updateStateVars is called.
  // We don't have a way to verify state variables inside friction object.

  PYLITH_METHOD_END;
} // testUpdateStateVars

// ----------------------------------------------------------------------
// Initialize FaultCohesiveDyn interface condition.
void
pylith::faults::TestFaultCohesiveDyn::_initialize(topology::Mesh* const mesh,
						  FaultCohesiveDyn* const fault,
						  topology::SolutionFields* const fields)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(fault);
  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(_quadrature);

  meshio::MeshIOAscii iohandler;
  iohandler.filename(_data->meshFilename);
  iohandler.read(mesh);

  // Set coordinate system
  spatialdata::geocoords::CSCart cs;
  cs.setSpaceDim(mesh->dimension());
  cs.initialize();
  mesh->coordsys(&cs);

  // Set scales
  // Most test data is insensitive to the scales because we set the fields directly.
  spatialdata::units::Nondimensional normalizer;
  normalizer.lengthScale(_data->lengthScale);
  normalizer.pressureScale(_data->pressureScale);
  normalizer.densityScale(_data->densityScale);
  normalizer.timeScale(_data->timeScale);
  topology::MeshOps::nondimensionalize(mesh, normalizer);
  
  _quadrature->initialize(_data->basis, _data->numQuadPts, _data->numBasis,
			  _data->basisDeriv,
			  _data->numQuadPts, _data->numBasis, _data->cellDim,
			  _data->quadPts, _data->numQuadPts, _data->cellDim,
			  _data->quadWts, _data->numQuadPts,
			  _data->spaceDim);
  
  // Setup prescribed traction perturbation
  delete _tractionPerturbation; _tractionPerturbation = new TractionPerturbation();CPPUNIT_ASSERT(_tractionPerturbation);
  _tractionPerturbation->label("traction perturbation");
  spatialdata::spatialdb::SimpleDB* db = new spatialdata::spatialdb::SimpleDB("initial tractions");CPPUNIT_ASSERT(db);
  spatialdata::spatialdb::SimpleIOAscii ioInitialTract;
  ioInitialTract.filename(_data->initialTractFilename);
  db->ioHandler(&ioInitialTract);
  delete _dbInitialTract; _dbInitialTract = db;
  _tractionPerturbation->dbInitial(db);
  fault->tractionPerturbation(_tractionPerturbation);

  // Setup friction
  spatialdata::spatialdb::SimpleDB* dbFriction = new spatialdata::spatialdb::SimpleDB("static friction");CPPUNIT_ASSERT(dbFriction);
  spatialdata::spatialdb::SimpleIOAscii ioFriction;
  if (2 == _data->spaceDim)
    ioFriction.filename("data/static_friction_2d.spatialdb");
  else if (3 == _data->spaceDim)
    ioFriction.filename("data/static_friction_3d.spatialdb");
  dbFriction->ioHandler(&ioFriction);
  delete _dbFriction; _dbFriction = dbFriction;
  friction::StaticFriction* friction = new pylith::friction::StaticFriction();CPPUNIT_ASSERT(friction);
  friction->label("static friction");
  friction->dbProperties(dbFriction);
  friction->normalizer(normalizer);
  _friction = friction;
  fault->frictionModel(friction);

  PetscInt labelSize;
  PetscErrorCode err;
  err = DMPlexGetStratumSize(mesh->dmMesh(), _data->label, 1, &labelSize);PYLITH_CHECK_ERROR(err);

  PetscInt firstFaultVertex = 0;
  PetscInt firstLagrangeVertex = labelSize;
  PetscInt firstFaultCell = labelSize;
  if (fault->useLagrangeConstraints())
    firstFaultCell += labelSize;
  fault->id(_data->id);
  fault->label(_data->label);
  fault->quadrature(_quadrature);
  
  fault->adjustTopology(mesh, &firstFaultVertex, &firstLagrangeVertex, &firstFaultCell);
  
  const PylithScalar upDir[3] = { 0.0, 0.0, 1.0 };
  
  fault->normalizer(normalizer);
  fault->initialize(*mesh, upDir);
  
  // Setup fields
  fields->add("residual", "residual");
  fields->add("disp(t)", "displacement");
  fields->add("dispIncr(t->t+dt)", "displacement_increment");
  fields->add("velocity(t)", "velocity");
  fields->add("dispIncr adjust", "dispIncr_adjust");
  fields->solutionName("dispIncr(t->t+dt)");

  const int spaceDim = _data->spaceDim;
  topology::Field& residual = fields->get("residual");
  residual.subfieldAdd("displacement", spaceDim, topology::Field::VECTOR);
  residual.subfieldAdd("lagrange_multiplier", spaceDim, topology::Field::VECTOR);
  residual.subfieldsSetup();
  residual.setupSolnChart();
  residual.setupSolnDof(spaceDim);
  fault->setupSolnDof(&residual);
  residual.allocate();
  residual.zero();

  fields->copyLayout("residual");

  fault->verifyConfiguration(*mesh);

  PYLITH_METHOD_END;
} // _initialize

// ----------------------------------------------------------------------
// Set values for fields and Jacobian.
void
pylith::faults::TestFaultCohesiveDyn::_setFieldsJacobian(topology::Mesh* const mesh,
							 FaultCohesiveDyn* const fault,
							 topology::SolutionFields* const fields,
							 topology::Jacobian* const jacobian,
							 const PylithScalar* const fieldIncr)
{ // _initialize
  PYLITH_METHOD_BEGIN;

  CPPUNIT_ASSERT(mesh);
  CPPUNIT_ASSERT(fault);
  CPPUNIT_ASSERT(fields);
  CPPUNIT_ASSERT(jacobian);
  CPPUNIT_ASSERT(_data);
  CPPUNIT_ASSERT(fieldIncr);

  const int spaceDim = _data->spaceDim;
  const PylithScalar lengthScale = 1.0;

  // Get vertices in mesh
  PetscDM dmMesh = mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
  topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
  PetscErrorCode err;
  PetscInt pStart, pEnd;

  // Set displacement values
  topology::Field& disp = fields->get("disp(t)");
  topology::VecVisitorMesh dispVisitor(disp);
  err = PetscSectionGetChart(disp.localSection(), &pStart, &pEnd);CPPUNIT_ASSERT(!err);
  PetscScalar* dispArray = dispVisitor.localArray();CPPUNIT_ASSERT(dispArray);
  for (PetscInt p = pStart, iVertex = 0; p < pEnd; ++p) {
    if (dispVisitor.sectionDof(p) > 0) {
      const PetscInt off = dispVisitor.sectionOffset(p);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dispVisitor.sectionDof(p));
      for(PetscInt d = 0; d < spaceDim; ++d) {
	dispArray[off+d] = _data->fieldT[iVertex*spaceDim+d] / lengthScale;
      } // for
      ++iVertex;
    } // if
  } // for

  // Set increment values
  topology::Field& dispIncr = fields->get("dispIncr(t->t+dt)");
  topology::VecVisitorMesh dispIncrVisitor(dispIncr);
  err = PetscSectionGetChart(dispIncr.localSection(), &pStart, &pEnd);CPPUNIT_ASSERT(!err);
  PetscScalar* dispIncrArray = dispIncrVisitor.localArray();CPPUNIT_ASSERT(dispIncrArray);
  for (PetscInt p = pStart, iVertex = 0; p < pEnd; ++p) {
    if (dispIncrVisitor.sectionDof(p) > 0) {
      const PetscInt off = dispIncrVisitor.sectionOffset(p);
      CPPUNIT_ASSERT_EQUAL(spaceDim, dispIncrVisitor.sectionDof(p));
      for(PetscInt d = 0; d < spaceDim; ++d) {
	dispIncrArray[off+d] = fieldIncr[iVertex*spaceDim+d] / lengthScale;
      } // for
      ++iVertex;
    } // if
  } // for

  // Setup Jacobian matrix
  const PetscInt nrows = (verticesStratum.size()+_data->numConstraintEdges) * spaceDim;
  const PetscInt ncols = nrows;
  int nrowsM = 0;
  int ncolsM = 0;
  PetscMat jacobianMat = jacobian->matrix();
  err = MatGetSize(jacobianMat, &nrowsM, &ncolsM);PYLITH_CHECK_ERROR(err);
  CPPUNIT_ASSERT_EQUAL(nrows, nrowsM);
  CPPUNIT_ASSERT_EQUAL(ncols, ncolsM);
  // We ignore the sparsity patterns in our tests
  err = MatSetOption(jacobianMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);PYLITH_CHECK_ERROR(err);

  int_array rows(nrows);
  int_array cols(ncols);
  for (int iRow=0; iRow < nrows; ++iRow)
    rows[iRow] = iRow;
  for (int iCol=0; iCol < ncols; ++iCol)
    cols[iCol] = iCol;
  err = MatSetValues(jacobianMat, nrows, &rows[0], ncols, &cols[0], _data->jacobian, INSERT_VALUES);PYLITH_CHECK_ERROR(err);
  jacobian->assemble("final_assembly");

  PYLITH_METHOD_END;
} // _setFieldsJacobian

// ----------------------------------------------------------------------
// Determine if point is a Lagrange multiplier constraint point.
bool
pylith::faults::TestFaultCohesiveDyn::_isConstraintEdge(const int point) const
{ // _isConstraintEdge
  PYLITH_METHOD_BEGIN;

  assert(_data);

  const int numConstraintEdges = _data->numConstraintEdges;
  bool isFound = false;
  for (int i=0; i < _data->numConstraintEdges; ++i)
    if (_data->constraintEdges[i] == point) {
      isFound = true;
      break;
    } // if
  PYLITH_METHOD_RETURN(isFound);
} // _isConstraintEdge


// End of file 
