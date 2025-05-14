import os
import inspect
from pathlib import Path
import scanpy as sc
from fastmcp import FastMCP , Context
from ..schema.io import *
from scmcp_shared.util import filter_args, forward_request
from scmcp_shared.logging_config import setup_logger
logger = setup_logger()

util_mcp = FastMCP("DECOUPLERNCP-Util-Server")


@util_mcp.tool()
async def check_samples(ctx: Context):
    """check the stored samples    
    """
    ads = ctx.request_context.lifespan_context
    adata_dic = ads.adata_dic
    return {"sampleid": [list(adata_dic[dk].keys()) for dk in adata_dic.keys()]}
    
    