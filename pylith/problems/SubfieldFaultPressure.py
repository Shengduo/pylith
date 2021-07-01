# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/problems/SubfieldFaultPressure.py
#
# @brief Python object for fault pressure subfield.
#
# Factory: subfield.

from .SolutionSubfield import SolutionSubfield


class SubfieldFaultPressure(SolutionSubfield):
    """Python object for fault pressure subfield.

    FACTORY: soln_subfield
    """

    import pythia.pyre.inventory

    from .SolutionSubfield import validateAlias
    userAlias = pythia.pyre.inventory.str("alias", default="fault_pressure", validator=validateAlias)
    userAlias.meta['tip'] = "Name for subfield."

    fieldName = "fault_pressure"

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="subfieldfaultpressure"):
        """Constructor.
        """
        SolutionSubfield.__init__(self, name)
        return

    def initialize(self, normalizer, spaceDim):
        """Initialize subfield metadata.
        """
        from pylith.topology.Field import Field
        self.vectorFieldType = Field.SCALAR
        self.scale = normalizer.getPressureScale()
        self._setComponents(spaceDim)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """Set members based using inventory.
        """
        SolutionSubfield._configure(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def soln_subfield():
    """Factory associated with SubfieldFaultPressure.
    """
    return SubfieldFaultPressure()


# End of file
