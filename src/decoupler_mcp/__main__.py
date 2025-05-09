
import asyncio
from fastmcp import FastMCP
from .server import io_mcp, if_mcp, pl_mcp

mcp = FastMCP("Decoupler-MCP-Server")


async def setup():
    await mcp.import_server("io", io_mcp)
    await mcp.import_server("if", if_mcp)
    await mcp.import_server("pl", pl_mcp)



if __name__ == "__main__":
    asyncio.run(setup())
    mcp.run()
    