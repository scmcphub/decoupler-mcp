import os
import inspect
from pathlib import Path
import scanpy as sc
from fastmcp import FastMCP , Context
from ..schema.io import *
from ..util import filter_args, forward_request
from ..logging_config import setup_logger
logger = setup_logger()

io_mcp = FastMCP("decouplerMCP-IO-Server")

@io_mcp.tool()
async def read(
    request: ReadModel, 
    ctx: Context,
    sampleid: str = Field(default=None, description="adata sampleid"),
    dtype: str = Field(default="exp", description="adata.X data type")
):
    """
    Read data from various file formats (h5ad, 10x, text files, etc.) or directory path.
    """
    try:
        result = await forward_request("io_read", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result        
        kwargs = request.model_dump()
        ads = ctx.request_context.lifespan_context
        if sampleid is not None:
            ads.active_id = sampleid
        else:
            ads.active_id = f"adata{len(ads.adata_dic[dtype])}"

        file = Path(kwargs.get("filename", None))
        logger.error(f"filename")
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
            raise FileNotFoundError(f"{kwargs['filename']} does not exist")
        adata.layers["counts"] = adata.X
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        ads.set_adata(adata, sampleid=sampleid, sdtype=dtype)
        return {"sampleid": sampleid or ads.active_id, "dtype": dtype, "adata": adata}
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e


@io_mcp.tool()
async def write(
    request: WriteModel, 
    ctx: Context,
    sampleid: str = Field(default=None, description="adata sampleid"),
    dtype: str = Field(default="exp", description="adata.X data type")
):
    """save adata into a file.
    """
    try:
        result = await forward_request("io_write", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(sampleid=sampleid, dtype=dtype)
        kwargs = request.model_dump()
        sc.write(kwargs["filename"], adata)
        return {"filename": kwargs["filename"], "msg": "success to save file"}
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e