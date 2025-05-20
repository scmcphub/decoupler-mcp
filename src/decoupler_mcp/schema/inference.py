from typing import Optional

from pydantic import Field
from scmcp_shared.schema import AdataModel

class ListInferenceMethodModel(AdataModel):
    """ListInferenceMethodModel"""
    pass


class ULMModel(AdataModel):
    """Input schema for decoupler's Univariate Linear Model (ULM) method.
    
    ULM fits a linear model for each sample and regulator, where the observed molecular 
    readouts are the response variable and the regulator weights are the explanatory one. 
    Target features with no associated weight are set to zero. The obtained t-value from 
    the fitted model is the activity of a given regulator.
    """
    
    source: str = Field(
        default="source",
        description="Column name in net with source nodes."
    )
    
    target: str = Field(
        default="target",
        description="Column name in net with target nodes."
    )
    
    weight: str = Field(
        default="weight",
        description="Column name in net with weights."
    )
    
    batch_size: int = Field(
        default=10000,
        description="Size of the samples to use for each batch. Increasing this will consume more memory but it will run faster."
    )
    
    min_n: int = Field(
        default=5,
        description="Minimum of targets per source. If less, sources are removed."
    )
    
    use_raw: bool = Field(
        default=True,
        description="Use raw attribute of mat if present."
    )
    
    key_added: Optional[str] = Field(
        default=None,
        description="Key under which the results will be stored in adata.obsm."
    )


class MLMModel(AdataModel):
    """Input schema for decoupler's Multivariate Linear Model (MLM) method.
    
    MLM fits a multivariate linear model for each sample, where the observed molecular readouts 
    are the response variable and the regulator weights are the covariates. Target features with 
    no associated weight are set to zero. The obtained t-values from the fitted model are the 
    activities of the regulators in the network.
    """
    
    source: str = Field(
        default="source",
        description="Column name in net with source nodes."
    )
    
    target: str = Field(
        default="target",
        description="Column name in net with target nodes."
    )
    
    weight: str = Field(
        default="weight",
        description="Column name in net with weights."
    )
    
    batch_size: int = Field(
        default=10000,
        description="Size of the samples to use for each batch. Increasing this will consume more memory but it will run faster."
    )
    
    min_n: int = Field(
        default=5,
        description="Minimum of targets per source. If less, sources are removed."
    )
    
    use_raw: bool = Field(
        default=True,
        description="Use raw attribute of mat if present."
    )
    
    key_added: Optional[str] = Field(
        default=None,
        description="Key under which the results will be stored in adata.obsm."
    )

class PathwayActivityModel(AdataModel):
    """Input schema for decoupler's Multivariate Linear Model (MLM) method.
    
    MLM fits a multivariate linear model for each sample, where the observed molecular readouts 
    are the response variable and the regulator weights are the covariates. Target features with 
    no associated weight are set to zero. The obtained t-values from the fitted model are the 
    activities of the regulators in the network.
    """
    source: str = Field(
        default="source",
        description="Column name in net with source nodes."
    )
    target: str = Field(
        default="target",
        description="Column name in net with target nodes."
    )
    weight: str = Field(
        default="weight",
        description="Column name in net with weights."
    )
    batch_size: int = Field(
        default=10000,
        description="Size of the samples to use for each batch. Increasing this will consume more memory but it will run faster."
    )
    min_n: int = Field(
        default=5,
        description="Minimum of targets per source. If less, sources are removed."
    )
    use_raw: bool = Field(
        default=True,
        description="Use raw attribute of mat if present."
    )    
    key_added: Optional[str] = Field(
        default=None,
        description="Key under which the results will be stored in adata.obsm."
    )
    organism: str = Field(
        default="human",
        description="The organism of interest. By default human."
    )
    
    top: int = Field(
        default=None,
        description="Number of genes per pathway to return. By default all of them."
    )

class TFActivityModel(AdataModel):
    """Input schema for decoupler's Univariate Linear Model (ULM) method.
    
    ULM fits a linear model for each sample and regulator, where the observed molecular 
    readouts are the response variable and the regulator weights are the explanatory one. 
    Target features with no associated weight are set to zero. The obtained t-value from 
    the fitted model is the activity of a given regulator.
    """
    source: str = Field(
        default="source",
        description="Column name in net with source nodes."
    )
    target: str = Field(
        default="target",
        description="Column name in net with target nodes."
    )
    weight: str = Field(
        default="weight",
        description="Column name in net with weights."
    )
    batch_size: int = Field(
        default=10000,
        description="Size of the samples to use for each batch. Increasing this will consume more memory but it will run faster."
    )
    min_n: int = Field(
        default=5,
        description="Minimum of targets per source. If less, sources are removed."
    )
    use_raw: bool = Field(
        default=True,
        description="Use raw attribute of mat if present."
    )
    key_added: Optional[str] = Field(
        default=None,
        description="Key under which the results will be stored in adata.obsm."
    )
    organism: str = Field(
        default="human",
        description="The organism of interest. By default human."
    )
