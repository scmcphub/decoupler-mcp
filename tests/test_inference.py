import pytest
from fastmcp import Client
import anndata
from pathlib import Path


@pytest.mark.asyncio 
async def test_activity(mcp_config):
    # Pass the server directly to the Client constructor
    testfile = Path(__file__).parent / "data/pbmc3k_processed.h5ad"
    outfile = Path(__file__).parent / "data/test.h5ad"
    async with Client(mcp_config) as client:
        result = await client.call_tool("io_read", {"request":{"filename": testfile}})
        assert "AnnData" in result[0].text

        result = await client.call_tool("if_pathway_activity", {"request":{"top": 500}})
        assert "progeny_mlm_estimate" in result[0].text

        result = await client.call_tool("if_tf_activity", {"request":{}})
        assert "collectri_ulm_estimate" in result[0].text

