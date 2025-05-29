
import pytest

@pytest.fixture
def mcp():
    from decoupler_mcp.server import DecouplerMCPManager
    mcp = DecouplerMCPManager("decoupler-mcp").mcp
    return mcp