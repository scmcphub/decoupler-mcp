from __future__ import annotations

from pydantic import BaseModel, Field, ConfigDict

class ActivityAdataInfo(BaseModel):
    """Input schema for the adata tool."""
    sampleid: str | None = Field(default=None, description="adata sampleid")
    adtype: str = Field(default="activity", description="the data type of stored adata.X, for analysis and plotting ")

    model_config = ConfigDict(
        extra="ignore"
    )
