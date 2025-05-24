from typing import Optional

from pydantic import Field, BaseModel

class ListInferenceMethodModel(BaseModel):
    """ListInferenceMethodModel"""
    pass


class ULMModel(BaseModel):
    """Input schema for decoupler's Univariate Linear Model (ULM) method.
    
    ULM fits a linear model for each sample and regulator, where the observed molecular 
    readouts are the response variable and the regulator weights are the explanatory one. 
    Target features with no associated weight are set to zero. The obtained t-value from 
    the fitted model is the activity of a given regulator.
    """
    tmin: int = Field(
        default=5,
        description="Minimum number of targets per source. Sources with fewer targets will be removed."
    )
    empty: bool = Field(
        default=True,
        description="Whether to remove empty observations (rows) or features (columns)."
    )
    raw : bool = Field(
        default=True,
        description="Whether to use the .raw attribute of anndata.AnnData."
    )
    bsize: int = Field(
        default=250000,
        description="controls how many observations are processed at once. "
    )

class MLMModel(BaseModel):
    """Input schema for decoupler's Multivariate Linear Model (MLM) method.
    
    MLM fits a multivariate linear model for each sample, where the observed molecular readouts 
    are the response variable and the regulator weights are the covariates. Target features with 
    no associated weight are set to zero. The obtained t-values from the fitted model are the 
    activities of the regulators in the network.
    """
    tmin: int = Field(
        default=5,
        description="Minimum number of targets per source. Sources with fewer targets will be removed."
    )
    empty: bool = Field(
        default=True,
        description="Whether to remove empty observations (rows) or features (columns)."
    )
    raw : bool = Field(
        default=True,
        description="Whether to use the .raw attribute of anndata.AnnData."
    )
    bsize: int = Field(
        default=250000,
        description="controls how many observations are processed at once. "
    )


class PathwayActivityModel(MLMModel):
    """Input schema for decoupler's Multivariate Linear Model (MLM) method.
    
    MLM fits a multivariate linear model for each sample, where the observed molecular readouts 
    are the response variable and the regulator weights are the covariates. Target features with 
    no associated weight are set to zero. The obtained t-values from the fitted model are the 
    activities of the regulators in the network.
    """
    organism: str = Field(
        default="human",
        description="The organism of interest. By default human."
    )
    top: int | float = Field(
        default=float("inf"),
        description="Number of genes per pathway to return. By default all of them."
    )
    thr_padj: float = Field(
        default=0.05,
        description="Significance threshold to trim interactions."
    )
    license: str = Field(
        default="academic",
        description="Which license to use, available options are: academic, commercial, or nonprofit."
    )


class TFActivityModel(ULMModel):
    """Input schema for decoupler's Univariate Linear Model (ULM) method.
    
    ULM fits a linear model for each sample and regulator, where the observed molecular 
    readouts are the response variable and the regulator weights are the explanatory one. 
    Target features with no associated weight are set to zero. The obtained t-value from 
    the fitted model is the activity of a given regulator.
    """

    organism: str = Field(
        default="human",
        description="The organism of interest. By default human."
    )
    remove_complexes: bool = Field(
        default=False,
        description="Whether to remove complexes."
    )
    license: str = Field(
        default="academic",
        description="Which license to use, available options are: academic, commercial, or nonprofit."
    )