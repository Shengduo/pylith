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
# @file pylith/faults/FaultPoroDiffusionCohesiveKin.py
#
# @brief Python object for a poroelastic fault surface with kinematic
# (prescribed) slip, fluid diffusion implemented with cohesive elements.
#
# Factory: fault

from .FaultCohesive import FaultCohesive
from .faults import FaultPoroDiffusionCohesiveKin as ModuleFaultPoroDiffusionCohesiveKin

# ITEM FACTORIES ///////////////////////////////////////////////////////


def eqsrcFactory(name):
    """Factory for earthquake source items.
    """
    from pythia.pyre.inventory import facility
    from .KinSrcStep import KinSrcStep
    return facility(name, family="eq_kinematic_src", factory=KinSrcStep)


class FaultPoroDiffusionCohesiveKin(FaultCohesive, ModuleFaultPoroDiffusionCohesiveKin):
    """Python object for a fault surface with kinematic (prescribed) slip and fluid diffusion
    implemented with cohesive elements.

    FACTORY: fault
    """

    import pythia.pyre.inventory

    from .SingleRupture import SingleRupture
    eqRuptures = pythia.pyre.inventory.facilityArray("eq_ruptures", itemFactory=eqsrcFactory, factory=SingleRupture)
    eqRuptures.meta['tip'] = "Kinematic earthquake sources information."

    from pylith.utils.NullComponent import NullComponent
    auxiliaryFieldDB = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=NullComponent)

    #from pylith.meshio.OutputFaultKin import OutputFaultKin
    #outputManager = pythia.pyre.inventory.facility("output", family="output_manager", factory=OutputFaultKin)
    #output.meta['tip'] = "Output manager associated with fault information."

    useBodyForce = pythia.pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in poroelastic fault equation."

    useSource = pythia.pyre.inventory.bool("use_source", default=False)
    useSource.meta['tip'] = "Include source_density term in poroelastic fault equation."

    useConstantPressureSource = pythia.pyre.inventory.bool("use_constant_pressure_source", default=False)
    useConstantPressureSource.meta['tip'] = "Include constant_pressure_source term in poroelastic fault equation."

    def __init__(self, name="faultporodiffusioncohesivekin"):
        """Initialize configuration.
        """
        FaultCohesive.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Pre-initializing fault '%s'." % self.label)

        FaultCohesive.preinitialize(self, problem)

        for eqsrc in self.eqRuptures.components():
            eqsrc.preinitialize()
        ModuleFaultPoroDiffusionCohesiveKin.setEqRuptures(
            self, self.eqRuptures.inventory.facilityNames(), self.eqRuptures.components())

        ModuleFaultPoroDiffusionCohesiveKin.useBodyForce(self, self.useBodyForce)
        ModuleFaultPoroDiffusionCohesiveKin.useSourceDensity(self, self.useSourceDensity)
        ModuleFaultPoroDiffusionCohesiveKin.useConstantPressureSource(self, self.useConstantPressureSource)
        return

    def verifyConfiguration(self):
        """Verify compatibility of configuration.
        """
        FaultCohesive.verifyConfiguration(self)
        ModuleFaultPoroDiffusionCohesiveKin.verifyConfiguration(self, self.mesh())

        for eqsrc in self.eqRuptures.components():
            eqsrc.verifyConfiguration()

        return

    def finalize(self):
        """Cleanup.
        """
        for eqsrc in self.eqRuptures.components():
            eqsrc.finalize()
        FaultCohesive.finalize(self)
        # self.output.close()
        # self.output.finalize()
        return

    def _configure(self):
        """Setup members using inventory.
        """
        FaultCohesive._configure(self)
        return

    def _createModuleObj(self):
        """Create handle to C++ FaultPoroDiffusionCohesiveKin.
        """
        ModuleFaultPoroDiffusionCohesiveKin.__init__(self)
        return


# Factories

def fault():
    """Factory associated with FaultPoroDiffusionCohesiveKin.
    """
    return FaultPoroDiffusionCohesiveKin()


# End of file
