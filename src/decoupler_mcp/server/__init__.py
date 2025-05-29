from scmcp_shared.server import ScanpyUtilMCP, ScanpyPlottingMCP
from scmcp_shared.server import BaseMCPManager
from scmcp_shared.server import io_mcp

from ..schema import ActivityAdataInfo
from .inference import if_mcp

ul_mcp = ScanpyUtilMCP(
    include_tools=["query_op_log", "check_samples"],
).mcp

pl_mcp = ScanpyPlottingMCP(
    exclude_tools=["highly_variable_genes", "diffmap"],
    AdataInfo=ActivityAdataInfo
).mcp
 

class DecouplerMCPManager(BaseMCPManager):
    """Manager class for Scanpy MCP modules."""
    
    def _init_modules(self):
        """Initialize available Decoupler MCP modules."""
        self.available_modules = {
            "io": io_mcp, 
            "if": if_mcp, 
            "ul": ul_mcp,
            "pl": pl_mcp,
        }
