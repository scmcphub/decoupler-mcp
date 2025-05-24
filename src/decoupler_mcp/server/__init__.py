import asyncio
from fastmcp import FastMCP
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from typing import Any


import scmcp_shared.server as shs

from scmcp_shared.server import ul_mcp
from scmcp_shared.util import filter_tools
from .inference import if_mcp



ads = shs.AdataState()

@asynccontextmanager
async def adata_lifespan(server: FastMCP) -> AsyncIterator[Any]:
    yield ads


decoupler_mcp = FastMCP(lifespan=adata_lifespan)


async def setup(modules=None):
    ul_mcp = await filter_tools(shs.ul_mcp, include_tools=["query_op_log", "check_samples"])
    pl_mcp = await filter_tools(shs.pl_mcp, exclude_tools=["highly_variable_genes", "diffmap"])
    mcp_dic = {
        "io": shs.io_mcp, 
        "if": if_mcp, 
        "pl": pl_mcp,
        "ul": ul_mcp
        }
    if modules is None or modules == "all":
        modules = mcp_dic.keys()
    for module in modules:
        await decoupler_mcp.import_server(module, mcp_dic[module])
 