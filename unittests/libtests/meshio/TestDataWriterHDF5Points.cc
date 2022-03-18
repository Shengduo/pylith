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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestDataWriterHDF5Points.hh" // Implementation of class methods

#include "data/DataWriterDataPoints.hh" // USES DataWriterDataPoints

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/meshio/OutputSolnPoints.hh" // USES OutputSolnPoints
#include "pylith/meshio/DataWriterHDF5.hh" // USES DataWriterHDF5
#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::meshio::TestDataWriterHDF5Points );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::meshio::TestDataWriterHDF5Points::setUp(void)
{ // setUp
    PYLITH_METHOD_BEGIN;

    TestDataWriterPoints::setUp();

    PYLITH_METHOD_END;
} // setUp

// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::meshio::TestDataWriterHDF5Points::tearDown(void)
{ // tearDown
    PYLITH_METHOD_BEGIN;

    TestDataWriterPoints::tearDown();

    PYLITH_METHOD_END;
} // tearDown

// ----------------------------------------------------------------------
// Test constructor
void
pylith::meshio::TestDataWriterHDF5Points::testConstructor(void)
{ // testConstructor
    PYLITH_METHOD_BEGIN;

    DataWriterHDF5 writer;

    CPPUNIT_ASSERT(!writer._viewer);

    PYLITH_METHOD_END;
} // testConstructor

// ----------------------------------------------------------------------
// Test openTimeStep() and closeTimeStep()
void
pylith::meshio::TestDataWriterHDF5Points::testTimeStep(void)
{ // testTimeStep
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    OutputSolnPoints output;
    DataWriterHDF5 writer;
    spatialdata::units::Nondimensional normalizer;
    normalizer.lengthScale(10.0);

    writer.filename(_data->timestepFilename);
    output.writer(&writer);
    output.setupInterpolator(_mesh, _data->points, _data->numPoints, _data->spaceDim, _data->names, _data->numPoints, normalizer);

    const PylithScalar t = _data->time;
    const int numTimeSteps = 1;
    output.open(*_mesh, numTimeSteps);
    output.writePointNames();
    output.openTimeStep(t, *_mesh);

    output.closeTimeStep();
    output.close();

    checkFile(_data->timestepFilename);

    PYLITH_METHOD_END;
} // testTimeStep

// ----------------------------------------------------------------------
// Test writeVertexField.
void
pylith::meshio::TestDataWriterHDF5Points::testWriteVertexField(void)
{ // testWriteVertexField
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);

    OutputSolnPoints output;
    DataWriterHDF5 writer;
    spatialdata::units::Nondimensional normalizer;
    normalizer.lengthScale(10.0);

    topology::Fields vertexFields(*_mesh);
    _createVertexFields(&vertexFields);

    writer.filename(_data->vertexFilename);
    output.writer(&writer);
    output.setupInterpolator(_mesh, _data->points, _data->numPoints, _data->spaceDim, _data->names, _data->numPoints, normalizer);

    const int nfields = _data->numVertexFields;

    const PylithScalar t = _data->time;
    const int numTimeSteps = 1;
    output.open(*_mesh, numTimeSteps);
    output.writePointNames();
    output.openTimeStep(t, *_mesh);
    for (int i=0; i < nfields; ++i) {
        topology::Field& field = vertexFields.get(_data->vertexFieldsInfo[i].name);
        output.appendVertexField(t, field, *_mesh);
    } // for
    output.closeTimeStep();
    output.close();

    checkFile(_data->vertexFilename);

    PYLITH_METHOD_END;
} // testWriteVertexField


// End of file
