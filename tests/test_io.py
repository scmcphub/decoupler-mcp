import pytest
from fastmcp import Client
from pathlib import Path


@pytest.mark.asyncio 
async def test_read_and_write(mcp_config):
    # Pass the server directly to the Client constructor
    test_dir = Path(__file__).parent / "data/pbmc3k_processed.h5ad"
    outfile = Path(__file__).parent / "data/test.h5ad"
    async with Client(mcp_config) as client:
        result = await client.call_tool("io_read", {"request":{"filename": test_dir}, "sampleid": "pbmc3k"})
        assert "AnnData" in result[0].text

        result = await client.call_tool("io_write", {"request":{"filename": outfile}, "sampleid": "pbmc3k"})
        assert outfile.exists()
        outfile.unlink()
