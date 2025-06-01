"""
Command-line interface for decoupler-mcp.

This module provides a CLI entry point for the decoupler-mcp package.
"""

from scmcp_shared.cli import MCPCLI
from .server import DecouplerMCPManager

cli = MCPCLI(
    name="decoupler-mcp", 
    help_text="Decoupler MCP Server CLI",
    manager=DecouplerMCPManager
)
