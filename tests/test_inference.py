import pytest
from fastmcp import Client
from pathlib import Path
from decoupler_mcp.server import decoupler_mcp,setup
import asyncio

asyncio.run(setup())

@pytest.mark.asyncio 
async def test_activity():
    # Pass the server directly to the Client constructor
    testfile = Path(__file__).parent / "data/pbmc3k_processed.h5ad"
    outfile = Path(__file__).parent / "data/test.h5ad"
    async with Client(decoupler_mcp) as client:
        result = await client.call_tool("io_read", {"request":{"filename": testfile}, "adinfo":{ "sampleid": "pbmc3k", "adtype": "exp"}})
        assert "AnnData" in result[0].text

        result = await client.call_tool(
            "if_pathway_activity", 
            {"request":{"top": 500} , "adinfo":{ "sampleid": "pbmc3k", "adtype": "exp"}}
        )
        assert "score_mlm" in result[0].text

        result = await client.call_tool(
            "if_tf_activity", 
            {"request":{}, "adinfo":{ "sampleid": "pbmc3k", "adtype": "exp"}}
        )
        assert "score_ulm" in result[0].text

