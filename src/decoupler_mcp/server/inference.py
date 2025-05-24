
from pydantic import Field
from fastmcp import FastMCP, Context
from fastmcp.exceptions import ToolError
import decoupler as dc
from ..schema.inference import *
from scmcp_shared.schema import AdataModel
from scmcp_shared.util import add_op_log, filter_args,forward_request, obsm2adata, get_ads
from scmcp_shared.logging_config import setup_logger

if_mcp = FastMCP("decoupler-mcp-inference-Server")

logger = setup_logger()


@if_mcp.tool()
async def pathway_activity(
    request: PathwayActivityModel,
    adinfo: AdataModel = AdataModel()
):
    """Pathway activity inference"""
    try:
        result = await forward_request("if_pathway_activity", request, adinfo)
        if result is not None:
            return result
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        kwargs = request.model_dump()
        progeny = dc.op.progeny(organism=kwargs["organism"], top=kwargs.get("top", None))
        func_kwargs = filter_args(request, dc.mt.mlm)
        dc.mt.mlm(data=adata, net=progeny, **func_kwargs)
        score = dc.pp.get_obsm(adata=adata, key='score_mlm')
        add_op_log(adata, dc.mt.mlm, func_kwargs, adinfo)
        sdtype = "activity"
        ads.set_adata(score, sampleid="score_mlm", sdtype=sdtype)
        return [
            {"sampleid": adinfo.sampleid or ads.active_id, "adtype": adinfo.adtype, "adata": adata},
            {"sampleid": "score_mlm", "adtype": sdtype, "adata": score},
        ]
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@if_mcp.tool()
async def tf_activity(
    request: TFActivityModel, 
    adinfo: AdataModel = AdataModel()
):
    """Transcription factor activity inference"""
    try:
        result = await forward_request("if_tf_activity", request, adinfo)
        if result is not None:
            return result
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        kwargs = request.model_dump()
        net = dc.op.collectri(organism=kwargs["organism"])
        func_kwargs = filter_args(request, dc.mt.ulm)
        dc.mt.ulm(data=adata, net=net, **func_kwargs)
        score = dc.pp.get_obsm(adata=adata, key='score_ulm')
        add_op_log(adata, dc.mt.ulm, func_kwargs, adinfo)
        sdtype = "activity"
        ads.set_adata(score, sampleid="score_ulm", sdtype=sdtype)
        return [    
            {"sampleid": adinfo.sampleid or ads.active_id, "adtype": adinfo.adtype, "adata": adata},
            {"sampleid": "score_ulm", "adtype": sdtype, "adata": score}
        ]
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

 