import pytest


@pytest.fixture
def mcp():
    from decoupler_mcp.server import DecouplerMCPManager
    from scmcp_shared.backend import AdataManager

    mcp = DecouplerMCPManager("decoupler-mcp", backend=AdataManager).mcp
    return mcp
