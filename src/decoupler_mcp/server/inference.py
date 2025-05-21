
from pydantic import Field
from fastmcp import FastMCP, Context
from fastmcp.exceptions import ToolError
import decoupler as dc
from ..schema.inference import *
from scmcp_shared.util import add_op_log, filter_args,forward_request, obsm2adata, get_ads
from scmcp_shared.logging_config import setup_logger

if_mcp = FastMCP("decoupler-mcp-inference-Server")

logger = setup_logger()


@if_mcp.tool()
async def pathway_activity(
    request: PathwayActivityModel
):
    """Pathway activity inference"""
    try:
        result = await forward_request("if_pathway_activity", request)
        if result is not None:
            return result
        ads = get_ads()
        adata = ads.get_adata(request=request)
        kwargs = request.model_dump()
        progeny = dc.get_progeny(organism=kwargs["organism"], top=kwargs.get("top", None))
        func_kwargs = filter_args(request, dc.run_mlm)
        dc.run_mlm(mat=adata, net=progeny, **func_kwargs)
        adata.obsm['progeny_mlm_estimate'] = adata.obsm['mlm_estimate'].copy()
        adata.obsm['progeny_mlm_pvals'] = adata.obsm['mlm_pvals'].copy()
        add_op_log(adata, dc.run_mlm, func_kwargs)
        estimate_adata = obsm2adata(adata, "progeny_mlm_estimate")
        pvals_adata = obsm2adata(adata, "progeny_mlm_pvals")
        sdtype = "activity"
        ads.set_adata(pvals_adata, sampleid="progeny_mlm_pvals", sdtype=sdtype)
        ads.set_adata(estimate_adata, sampleid="progeny_mlm_estimate", sdtype=sdtype)
        return [
            {"sampleid": sampleid or ads.active_id, "adtype": request.adtype, "adata": adata},
            {"sampleid": "progeny_mlm_pvals", "adtype": sdtype, "adata": pvals_adata},
            {"sampleid": "progeny_mlm_estimate", "adtype": sdtype, "adata": estimate_adata}
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
):
    """Transcription factor activity inference"""
    try:
        result = await forward_request("if_tf_activity", request)
        if result is not None:
            return result
        ads = get_ads()
        adata = ads.get_adata(request=request)
        kwargs = request.model_dump()
        net = dc.get_collectri(organism=kwargs["organism"], split_complexes=False) 
        func_kwargs = filter_args(request, dc.run_ulm)
        dc.run_ulm(mat=adata, net=net, **func_kwargs)
        adata.obsm['collectri_ulm_estimate'] = adata.obsm['ulm_estimate'].copy()
        adata.obsm['collectri_ulm_pvals'] = adata.obsm['ulm_pvals'].copy()
        add_op_log(adata, dc.run_ulm, func_kwargs)
        estimate_adata = obsm2adata(adata, "collectri_ulm_estimate")
        pvals_adata = obsm2adata(adata, "collectri_ulm_pvals")
        sdtype = "activity"
        ads.set_adata(pvals_adata, sampleid="collectri_ulm_pvals", sdtype=sdtype)
        ads.set_adata(estimate_adata, sampleid="collectri_ulm_estimate", sdtype=sdtype)
        return [    
            {"sampleid": request.sampleid or ads.active_id, "adtype": request.adtype, "adata": adata},
            {"sampleid": "collectri_ulm_pvals", "adtype": sdtype, "adata": pvals_adata},
            {"sampleid": "collectri_ulm_estimate", "adtype": sdtype, "adata": estimate_adata}
        ]
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)

