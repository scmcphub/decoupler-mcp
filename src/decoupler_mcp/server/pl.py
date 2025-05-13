import os
import inspect
from functools import partial
import scanpy as sc
from fastmcp import FastMCP, Context
from ..schema.pl import *
from pathlib import Path
from ..logging_config import setup_logger
from ..util import filter_args, set_fig_path, add_op_log,forward_request,obsm2adata
from ..logging_config import setup_logger


logger = setup_logger()



pl_mcp = FastMCP("DecouplerMCP-pl-Server")


@pl_mcp.tool()
async def violin(
    request: ViolinModel, 
    ctx: Context, 
    dtype: str = Field(default="activity", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for plotting")
):
    """Plot violin plot of one or more variables in adata.var and adata.obs or adata.obsm"""
    try:
        result = await forward_request("pl_violin", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        func_kwargs = filter_args(request, sc.pl.violin)
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        func_kwargs.pop("return_fig", True)
        func_kwargs["show"] = False
        func_kwargs["save"] = ".png"
        fig = sc.pl.violin(adata, **func_kwargs)
        fig_path = set_fig_path("violin", **func_kwargs)
        add_op_log(adata, sc.pl.violin, func_kwargs)
        return {"figpath": fig_path}
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e


@pl_mcp.tool()
async def matrixplot(
    request: MatrixplotModel, 
    ctx: Context, 
    dtype: str = Field(default="activity", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for plotting")
):
    """matrixplot, Create a heatmap of the mean expression values per group of each var_names."""
    try:
        result = await forward_request("pl_matrixplot", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        func_kwargs = filter_args(request, sc.pl.matrixplot)
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        func_kwargs.pop("return_fig", True)
        func_kwargs["show"] = False
        func_kwargs["save"] = ".png"
        if request.use_obsm is not None:
            adata = obsm2adata(adata, request.use_obsm)
        fig = sc.pl.matrixplot(adata, **func_kwargs)
        fig_path = set_fig_path("matrixplot", **func_kwargs)
        add_op_log(adata, sc.pl.matrixplot, func_kwargs)
        return {"figpath": fig_path}
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e 


@pl_mcp.tool()
async def embedding(
    request: EmbeddingModel, 
    ctx: Context, 
    dtype: str = Field(default="activity", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for plotting")
):
    """Scatter plot for user specified embedding basis (e.g. umap, tsne, etc)."""
    try:
        result = await forward_request("pl_embedding", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result       
        func_kwargs = filter_args(request, sc.pl.embedding)
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        func_kwargs.pop("return_fig", True)
        func_kwargs["show"] = False
        func_kwargs["save"] = ".png"
        try:      
            fig = sc.pl.embedding(adata, **func_kwargs)
        except KeyError as e:
            raise KeyError(f"Key '{e}' not found in adata.var and adata.obs. please check {e}, or {dtype}.")
        except Exception as e:
            raise e
        fig_path = set_fig_path("embedding", **func_kwargs)
        add_op_log(adata, sc.pl.embedding, func_kwargs)
        return {"figpath": fig_path}
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e
