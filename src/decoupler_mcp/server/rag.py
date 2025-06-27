from fastmcp import FastMCP
from scmcp_shared.agent import rag_agent
from pydantic import Field

rag_mcp = FastMCP("RAG-Server")


@rag_mcp.tool(tags={"rag"})
def retrieve_knowledge(
    task: str = Field(description="The tasks or questions that needs to be solved"),
):
    """search function and parameters that can be used to solve the user's tasks or questions"""
    return rag_agent(task, software="decoupler")
