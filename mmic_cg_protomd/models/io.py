from cmselemental.models.procedures import ProcInput
from mmelemental.models.base import ProtoModel
from mmic_cg.models import CoarseInput, CoarseOutput
from mmelemental.models import Molecule
from pydantic import Field
from typing import Optional


__all__ = ["InputModel", "OutputModel"]


class ComputeProtomdInput(ProcInput):
    """ An input model for mmic_cg_protomd. """
    proc_input: CoarseInput = Field(..., description="Procedure input schema.")

    molecule: str = Field(
        None,
        description="The file of the coordinates of the atoms in the system. Should be a .gro file.",
    )


#molecule is None by default because ensemble might be added later

class ComputeProtomdOutput(ProtoModel):
    """ An output model for mmic_cg_protomd. """

    proc_input: CoarseOutput = Field(..., description="Procedure input schema.")

    molecule: str = Field(None, description="Molecule file string object that stores the coarse grained molecules")
