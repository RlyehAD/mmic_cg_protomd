from cmselemental.models.procedures import ProcInput, ProcOutput
from mmelemental.models.base import ProtoModel
from mmelemental.models import Molecule, ForceField, TrajInput, ForcesInput
from pydantic import Field, validator
from typing import Optional, Dict, List, Tuple, Union

__all__ = ["InputModel", "OutputModel"]


class InputModel(ProtoModel):
    """ An input model for mmic_cg_protomd. """

    ...


class OutputModel(ProtoModel):
    """ An output model for mmic_cg_protomd. """

    ...
