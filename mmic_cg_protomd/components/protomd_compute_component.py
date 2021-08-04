"""
Components in mmic_cg_protomd.
"""
import mmelemental
import numpy as np 
import MDAnalysis as mda 
import proto_md.subsystems as ss
#from mmic_cg_protomd.models import ComputeProtomdInput, ComputeProtomdOutput
from mmelemental.util.files import random_file
from mmic_cg.models.proc import CoarseInput, CoarseOutput
from typing import List, Tuple, Optional
from mmic.components.blueprints.generic_component import GenericComponent
from ..models import *

__all__ = ["Component"]


class CoarseProtoMDComponent(GenericComponent):
    """ A sample component that defines the 3 required methods. """

    @classmethod
    def input(cls):
        return CoarseInput

    @classmethod
    def output(cls):
        return CoarseOutput

    def execute(
        self,
        inputs: CoarseInput,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, CoarseOutput]:

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
        gro_file = random_file(suffix=".gro") 
        mol.to_file(gro_file, translator="mmic_parmed")

        cg_compute = ComputeProtomdInput(
            proc_input=inputs,
            molecule=gro_file,
            schema_name=inputs.schema_name,
            schema_version=inputs.schema_version,
        )

        universe = mda.Universe(cg_compute.molecule)

        nCG, SS = eval(factory)(**inputs.keywords)
        [sub.universe_changed(universe) for sub in SS]
        [sub.equilibrated() for sub in SS]
        cg_pos = [sub.ComputeCG(universe.atoms.positions) for sub in SS]
        # Need to write a 'if' to deal with velocities here
        #cg_vel = [sub.ComputeCG_Vel(universe.atoms.)]

        mols = {}
        for i in cg_pos: # Go through every subsystem
            num_cg_atoms = len(i)
            symbols = np.array(["cg"]*num_cg_atoms)
            mol = mmelemental.models.Molecule(schema_name="mmschema_molecule", schema_version=1.0, symbols=symbols,name="cg_atoms"+str(i), geometry=i)
            mols["cg_mol"+str(i)] = mol


        return True, CoarseOutput(
            proc_input=inputs,
            schema_name=inputs.schema_name,
            schema_version=inputs.schema_version,
            molecule=mols,
            success=True,
            )
