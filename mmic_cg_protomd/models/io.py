from cmselemental.models.procedures import InputProc
from mmelemental.models.base import ProtoModel
from mmic_cg.models import InputCoarse, OutputCoarse
from mmelemental.models import Molecule
from pydantic import Field
from typing import Optional


__all__ = ["InputComputeProtomd", "OutputComputeProtomd"]


class InputComputeProtomd(InputProc):
    """ An input model for mmic_cg_protomd. """
    proc_input: InputCoarse = Field(..., description="Procedure input schema.")

    molecule: str = Field(
        None,
        description="The file of the coordinates of the atoms in the system. Should be a .gro file.",
    )
    """
    velocities: bool = Field(
    	False,
    	description="If set as True, velocites are going to be coarse-grained. Please set it to False if no velocities is included in this step."
    	)
	"""

#molecule is None by default because ensemble might be added later

class OutputComputeProtomd(ProtoModel):
    """ An output model for mmic_cg_protomd. """

    proc_input: InputCoarse = Field(..., description="Procedure input schema.")

    molecule: str = Field(None, description="Molecule file string object that stores the coarse grained molecules")
