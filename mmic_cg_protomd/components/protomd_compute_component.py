"""
Components in mmic_cg_protomd.
"""
import mmelemental
import numpy as np 
import os
import MDAnalysis as mda 
import tempfile
import proto_md.subsystems as ss
#from mmic_cg_protomd.models import ComputeProtomdInput, ComputeProtomdOutput
from mmic_cg.models.proc import InputCoarse, OutputCoarse
from typing import List, Tuple, Optional
from mmic.components.blueprints.generic_component import GenericComponent
from cmselemental.util.decorators import classproperty
from ..models import *


__all__ = ["CoarseProtoMDComponent"]


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
            inputs = self.input(**inputs)

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
        gro_file = tempfile.NamedTemporaryFile(suffix=".gro").name
        clean_files = [gro_file] 
        mol.to_file(gro_file, translator="mmic_parmed")

        cg_compute = InputComputeProtomd(
            proc_input=inputs,
            molecule=gro_file,
            schema_name=inputs.schema_name,
            schema_version=inputs.schema_version,
        )

        universe = mda.Universe(cg_compute.molecule)

        nCG, SS = eval(factory)(**inputs.method_keywords)
 
        [sub.universe_changed(universe) for sub in SS]
        [sub.equilibrated() for sub in SS]
        cg_pos = [sub.computeCG_pos(universe.atoms.positions) for sub in SS]
        

        try:
            universe.atoms.velocities
            cg_vel = [sub.computeCG_vel(universe.atoms.velocities) for sub in SS]
        except Exception:
            cg_vel = [pos*0 for pos in cg_pos]

        mols = {}
        for i, j in zip(cg_pos, cg_vel):
            num_cg_atoms = len(i)
            symbols = np.array(["cg"]*num_cg_atoms)
            mol = mmelemental.models.Molecule(schema_name="mmschema_molecule", schema_version=1.0, symbols=symbols,name="cg_atoms"+str(j), geometry=i, velocities=j)
            mols["cg_mol"+mol.name] = mol



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