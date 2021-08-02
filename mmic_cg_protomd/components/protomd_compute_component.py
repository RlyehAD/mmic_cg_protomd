"""
Components in mmic_cg_protomd.
"""
import numpy as np 
import MDAnalysis as mda 
import proto_md.subsystems as ss
#from mmic_cg_protomd.models import ComputeProtomdInput, ComputeProtomdOutput
from typing import List, Tuple, Optional
from mmic.components.blueprints.generic_component import GenericComponent
from ..models import *

__all__ = ["Component"]


class Component(GenericComponent):
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
    ) -> Tuple[bool, OutputModel]:

        # Convert input dictionary to model
        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        # Load the cg method
        supported_methods = ["spacewarping":"ss.SpaceWarpingSubsystemFactory", 
            "comtinuumfieldvariable":"ss.ContinuumSubsystemFactory",
            "centerofmass":"ss.RigidSubsystemFactory",]

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

        num_cg_atoms = np.size(cg_pos,1)
        #resindex = [0]*num_cg_atoms
        #u_cg = mda.Universe.empty(num_cg_atoms, 1, resindex, trajectory=True)
        #u.add_TopologyAttr('name',['cg']*num_cg_atoms)
        #u.add_TopologyAttr('resname', ['cg_atoms'])
        #u.add_TopologyAttr('resnum', CG) 
        #u_cg.load_new(cg_pos)
        #mol = Molecule.from_data(data=u,dtype="mdanalysis")

        # Populate kwargs from inputs
        return True, OutputModel(**kwargs)
