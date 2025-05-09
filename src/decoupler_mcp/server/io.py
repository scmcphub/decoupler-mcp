import os
import inspect
from pathlib import Path
import scanpy as sc
from fastmcp import FastMCP , Context
from ..schema.io import *
from ..util import filter_args, forward_request
from ..logging_config import setup_logger
logger = setup_logger()

io_mcp = FastMCP("DECOUPLERNCP-IO-Server")


@io_mcp.tool()
async def read(
    request: ReadModel, 
    ctx: Context,
    sampleid: str = Field(default=None, description="adata sampleid"),
    sdtype: str = Field(default="exp", description="adata.X data type")
):
    """
    Read data from various file formats (h5ad, 10x, text files, etc.) or directory path.
    """
    kwargs = request.model_dump()

    result = await forward_request("io_read", kwargs)
    if result is not None:
        return result
    
    ads = ctx.request_context.lifespan_context
    if sampleid is not None:
        ads.active_id = sampleid
    else:
        ads.active_id = f"adata{len(ads.adata_dic[sdtype])}"

    file = Path(kwargs.get("filename", None))
    if file.is_dir():
        kwargs["path"] = kwargs["filename"]
        func_kwargs = filter_args(request, sc.read_10x_mtx)
        adata = sc.read_10x_mtx(kwargs["path"], **func_kwargs)
    elif file.is_file():
        func_kwargs = filter_args(request, sc.read)
        adata = sc.read(**func_kwargs)
        if not kwargs.get("first_column_obs", True):
            adata = adata.T
    else:
        raise ValueError("filename must be a file or a directory")
    adata.layers["counts"] = adata.X
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    ads.set_adata(adata, sampleid=sampleid, sdtype="exp")
    return {"sampleid": sampleid or ads.active_id, "adata": adata, "dtype": sdtype}



@io_mcp.tool()
async def write(
    request: WriteModel, 
    ctx: Context,
    sampleid: str = Field(default=None, description="adata sampleid"),
    dtype: str = Field(default="exp", description="adata.X data type")
):
    """save adata into a file.
    """
    result = await forward_request("io_write", request.model_dump())
    if result is not None:
        return result
    ads = ctx.request_context.lifespan_context
    adata = ads.get_adata(sampleid=sampleid, dtype=dtype)    
    kwargs = request.model_dump()
    sc.write(kwargs["filename"], adata)
    return {"filename": kwargs["filename"], "msg": "success to save file"}
