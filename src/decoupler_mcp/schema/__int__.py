from __future__ import annotations

from .ccc import *
from .io import *
from .pl import *
from scmcp_shared.schema import AdataModel


ExpAdataModel = AdataModel


class ActivityAdataModel(AdataModel):
    """Input schema for the adata tool."""
    sampleid: str | None = Field(default=None, description="adata sampleid")
    adtype: str = Field(default="activity", description="the data type of stored adata.X, for analysis and plotting ")

    model_config = ConfigDict(
        extra="ignore"
    )
