from cmselemental.models.procedures import ProcInput
from mmelemental.models.base import ProtoModel
from mmic_cg.models import CoarseInput, CoarseOutput
from mmelemental.models import Molecule
from pydantic import Field


__all__ = ["InputModel", "OutputModel"]


class ComputeProtomdInput(ProcInput):
    """ An input model for mmic_cg_protomd. """
    proc_input: CoarseInput = Field(..., description="Procedure input schema.")

    molecule: str = Field(
        ...,
        description="The file of the coordinates of the atoms in the system. Should be a .gro file.",
    )



class ComputeProtomdOutput(ProtoModel):
    """ An output model for mmic_cg_protomd. """

    proc_input: CoarseOutput = Field(..., description="Procedure input schema.")

    molecule: str = Field(..., description="Molecule file string object that stores the coarse grained molecules")
