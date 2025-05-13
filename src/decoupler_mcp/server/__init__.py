import asyncio
from fastmcp import FastMCP
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from typing import Any


from .io import io_mcp
from .inference import if_mcp
from .pl import pl_mcp
from .util import util_mcp

class AdataState:
    def __init__(self):
        self.adata_dic = {"exp": {}, "activity": {}}
        self.active_id = None
        self.metadata = {}
        
    def get_adata(self, sampleid=None, dtype="exp"):
        try:
            if self.active_id is None:
                return None
            sampleid = sampleid or self.active_id
            return self.adata_dic[dtype][sampleid]
        except KeyError as e:
            raise KeyError(f"Key {e} not found in adata_dic[{dtype}].Please check the sampleid or dtype.")
        except Exception as e:
            raise Exception(f"Error: {e}")
    
    def set_adata(self, adata, sampleid=None, sdtype="exp"):
        sampleid = sampleid or self.active_id
        self.adata_dic[sdtype][sampleid] = adata



ads = AdataState()

@asynccontextmanager
async def adata_lifespan(server: FastMCP) -> AsyncIterator[Any]:
    yield ads


decoupler_mcp = FastMCP("decoupler-MCP-Server", lifespan=adata_lifespan)

async def setup(modules=None):
    mcp_dic = {
        "io": io_mcp, 
        "if": if_mcp, 
        "pl": pl_mcp, 
        "util": util_mcp
        }
    if modules is None or modules == "all":
        modules = ["io", "if", "pl", "util"]
    for module in modules:
        await decoupler_mcp.import_server(module, mcp_dic[module])
