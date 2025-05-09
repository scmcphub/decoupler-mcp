import os
from pathlib import Path

from fastmcp import FastMCP, Context
import decoupler as dc
from ..schema.inference import *
from ..util import add_op_log, filter_args,forward_request
from ..logging_config import setup_logger

if_mcp = FastMCP("decoupler-mcp-inference-Server")

logger = setup_logger()


@if_mcp.tool()
async def pathway_activity(request: PathwayActivityModel, ctx: Context):
    """Pathway activity inference"""
    result = await forward_request("infer_pathway_activity", request.model_dump())
    if result is not None:
        return result
    
    kwargs = request.model_dump()
    progeny = dc.get_progeny(organism=kwargs["organism"], top=kwargs.get("top", None))
    ads = ctx.request_context.lifespan_context
    adata = ads.get_adata()
    func_kwargs = filter_args(request, dc.run_mlm)
    dc.run_mlm(mat=adata, net=progeny, **func_kwargs)
    adata.obsm['progeny_mlm_estimate'] = adata.obsm['mlm_estimate'].copy()
    adata.obsm['progeny_mlm_pvals'] = adata.obsm['mlm_pvals'].copy()
    add_op_log(adata, dc.run_mlm, func_kwargs)
    ads.set_adata(adata)
    return adata


@if_mcp.tool()
async def tf_activity(request: TFActivityModel, ctx: Context):
    """Transcription factor activity inference"""
    result = await forward_request("infer_tf_activity", request.model_dump())
    if result is not None:
        return result
    
    kwargs = request.model_dump()
    net = dc.get_collectri(organism=kwargs["organism"], split_complexes=False) 
    ads = ctx.request_context.lifespan_context
    adata = ads.get_adata()
    func_kwargs = filter_args(request, dc.run_ulm)
    dc.run_ulm(mat=adata, net=net, **func_kwargs)
    adata.obsm['collectri_ulm_estimate'] = adata.obsm['ulm_estimate'].copy()
    adata.obsm['collectri_ulm_pvals'] = adata.obsm['ulm_pvals'].copy()
    add_op_log(adata, dc.run_ulm, func_kwargs)
    ads.set_adata(adata)
    return adata
