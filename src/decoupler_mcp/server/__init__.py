from scmcp_shared.server.preset import ScanpyUtilMCP, ScanpyPlottingMCP
from scmcp_shared.mcp_base import BaseMCPManager
from scmcp_shared.server.preset import io_mcp

from ..schema import ActivityAdataInfo
from .inference import if_mcp
from scmcp_shared.server.code import nb_mcp
from .rag import rag_mcp
from scmcp_shared.server.auto import auto_mcp

ul_mcp = ScanpyUtilMCP(
    include_tools=["query_op_log", "check_samples"],
).mcp

pl_mcp = ScanpyPlottingMCP(
    exclude_tools=["highly_variable_genes", "diffmap"], AdataInfo=ActivityAdataInfo
).mcp


class DecouplerMCPManager(BaseMCPManager):
    """Manager class for Scanpy MCP modules."""

    def init_mcp(self):
        """Initialize available Decoupler MCP modules."""
        self.available_modules = {
            "io": io_mcp,
            "if": if_mcp,
            "ul": ul_mcp,
            "pl": pl_mcp,
            "auto": auto_mcp,
            "nb": nb_mcp,
            "rag": rag_mcp,
        }
