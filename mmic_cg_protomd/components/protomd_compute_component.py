"""
Components in mmic_cg_protomd.
"""
import mmelemental
import numpy as np 
import os
import MDAnalysis as mda 
import proto_md.subsystems as ss
import tempfile
#from mmic_cg_protomd.models import ComputeProtomdInput, ComputeProtomdOutput
from mmic_cg.models.proc import InputCoarse, OutputCoarse
from typing import List, Tuple, Optional
from mmic.components.blueprints.generic_component import GenericComponent
from cmselemental.util.decorators import classproperty
from ..models import *


__all__ = ["Component"]


class CoarseProtoMDComponent(GenericComponent):
    """ A sample component that defines the 3 required methods. """

    @classproperty
    def input(cls):
        return InputCoarse

    @classproperty
    def output(cls):
        return OutputCoarse

    @classproperty
    def version(cls):
        return None

    def execute(
        self,
        inputs: InputCoarse,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, OutputCoarse]:

        # Convert input dictionary to model
        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        # Load the cg method
        supported_methods = {"spacewarping":"ss.SpaceWarpingSubsystemFactory", 
            "comtinuumfieldvariable":"ss.ContinuumSubsystemFactory",
            "centerofmass":"ss.RigidSubsystemFactory"}

        for key, val in supported_methods.items():
            if key == inputs.method:
                factory = val
                break
        else: 
            raise ValueError("Sorry, the method has not been supported by Proto_MD")

        # write the gro file
        mols = inputs.molecule
        mol_name, mol = list(mols.items()).pop()
        gro_file = tempfile.NamedTemporaryFile(suffix=".gro")
        clean_files = [gro_file] 
        mol.to_file(gro_file, translator="mmic_parmed")

        cg_compute = InputComputeProtomd(
            proc_input=inputs,
            molecule=gro_file,
            schema_name=inputs.schema_name,
            schema_version=inputs.schema_version,
        )

        universe = mda.Universe(cg_compute.molecule)
        vel = True # cg velocity is going to be computed by default

        nCG, SS = eval(factory)(**inputs.method_keywords)
        [sub.universe_changed(universe) for sub in SS]
        [sub.equilibrated() for sub in SS]
        cg_pos = [sub.ComputeCG(universe.atoms.positions) for sub in SS]
        # Need to write a 'if' to deal with velocities here

        try:
            cg_vel = [sub.ComputeCG_Vel(universe.atoms.velocities) for sub in SS]
        except Exception:
            cg_vel = cg_pos * 0
            vel = False
        """
        if inputs.cg_options["velocities"] == True:
            cg_vel = [sub.ComputeCG_Vel(universe.atoms.velocities) for sub in SS]
        """

        mols = {}
        j = 1# j means the jth sub system
        for i in cg_pos: # Go through every subsystem
            num_cg_atoms = len(i)
            symbols = np.array(["cg"]*num_cg_atoms)
            if vel == True:
                mol = mmelemental.models.Molecule(schema_name="mmschema_molecule", schema_version=1.0, symbols=symbols,name="cg_atoms"+str(j), geometry=i, velocities=cg_vel[j-1])
            else:
                mol = mmelemental.models.Molecule(schema_name="mmschema_molecule", schema_version=1.0, symbols=symbols,name="cg_atoms"+str(j), geometry=i)
            mols["cg_mol"+str(j)] = mol
            j = j+1


        self.cleanup(clean_files)
        
        return True, OutputCoarse(
            proc_input=inputs,
            schema_name=inputs.schema_name,
            schema_version=inputs.schema_version,
            molecule=mols,
            success=True,
            )


    @staticmethod
    def cleanup(remove: List[str]):
        for item in remove:
            if os.path.isdir(item):
                shutil.rmtree(item)
            elif os.path.isfile(item):
                os.remove(item)