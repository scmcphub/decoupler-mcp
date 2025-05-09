
import pytest

@pytest.fixture
def mcp_config():
    return {
        "mcpServers": {
            "decoupler-mcp": {
                "command": "decoupler-mcp",
                "args": ["run"]
            }
        }
    }