"""
Command-line interface for decoupler-mcp.

This module provides a CLI entry point for the decoupler-mcp package.
"""

import asyncio
import os
import sys
import typer
from enum import Enum
from typing import Optional




app = typer.Typer(
    name="decoupler-mcp",
    help="Decoupler MCP Server CLI",
    add_completion=False,
    no_args_is_help=True,  # Show help if no args provided    
    )

# Define enums for choices
class Module(str, Enum):
    ALL = "all"


class Transport(str, Enum):
    STDIO = "stdio"
    SSE = "sse"
    SHTTP = "shttp"


@app.command(name="run")
def run(
    log_file: Optional[str] = typer.Option(None, "--log-file", help="Log file path, use stdout if None"),
    module: Module = typer.Option(Module.ALL, "-m", "--module", help="Specify modules to load", 
                              case_sensitive=False),
    transport: Transport = typer.Option(Transport.STDIO, "-t", "--transport", help="Specify transport type", 
                                 case_sensitive=False),
    port: int = typer.Option(8000, "-p", "--port", help="transport port"),
    host: str = typer.Option("127.0.0.1", "--host", help="transport host"),
    forward: str = typer.Option(None, "-f", "--forward", help="forward request to another server"),
):
    """Start Decoupler MCP Server"""
    
    # Set environment variables
    if log_file is not None:
        os.environ['SCMCP_LOG_FILE'] = log_file
    if forward is not None:
        os.environ['SCMCP_FORWARD'] = forward
    os.environ['SCMCP_TRANSPORT'] = transport.value
    os.environ['SCMCP_HOST'] = host
    os.environ['SCMCP_PORT'] = str(port)
    os.environ['SCMCP_MODULE'] = module.value
    
    from .server import decoupler_mcp, setup
    asyncio.run(setup())
    from scmcp_shared.util import add_figure_route

    if transport == Transport.STDIO:
        decoupler_mcp.run()
    elif transport == Transport.SSE:
        add_figure_route(decoupler_mcp)
        decoupler_mcp.run(
                transport="sse",
                host=host, 
                port=port, 
                log_level="info"
            )
    elif transport == Transport.SHTTP:
        add_figure_route(decoupler_mcp)
        decoupler_mcp.run(
                transport="streamable-http",
                host=host, 
                port=port, 
                log_level="info"
            )

@app.callback()
def main():
    """decoupler MCP CLI root command."""
    pass