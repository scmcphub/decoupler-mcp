import inspect
import os
from pathlib import Path
from starlette.responses import FileResponse, Response

PKG = __package__.split('.')[0].upper()

def filter_args(request, func):
    # sometime,it is a bit redundant, but I think it adds robustness in case the function parameters change
    kwargs = request.model_dump()
    args = request.model_fields_set
    parameters = inspect.signature(func).parameters
    func_kwargs = {k: kwargs.get(k) for k in args if k in parameters}
    return func_kwargs


def add_op_log(adata, func, kwargs):
    import hashlib
    import json
    
    if "operation" not in adata.uns:
        adata.uns["operation"] = {}
        adata.uns["operation"]["op"] = {}
        adata.uns["operation"]["opid"] = []
    # Handle different function types to get the function name
    if hasattr(func, "func") and hasattr(func.func, "__name__"):
        # For partial functions, use the original function name
        func_name = func.func.__name__
    elif hasattr(func, "__name__"):
        func_name = func.__name__
    elif hasattr(func, "__class__"):
        func_name = func.__class__.__name__
    else:
        func_name = str(func)
    new_kwargs = {}
    for k,v in kwargs.items():
        if isinstance(v, tuple):
            new_kwargs[k] = list(v)
        else:
            new_kwargs[k] = v
    try:
        kwargs_str = json.dumps(new_kwargs, sort_keys=True)
    except:
        kwargs_str = str(new_kwargs)
    hash_input = f"{func_name}:{kwargs_str}"
    hash_key = hashlib.md5(hash_input.encode()).hexdigest()
    adata.uns["operation"]["op"][hash_key] = {func_name: new_kwargs}
    adata.uns["operation"]["opid"].append(hash_key)
    from .logging_config import setup_logger
    logger = setup_logger()
    logger.info(f"{func}: {new_kwargs}")


def savefig(fig, file):
    try:
        file_path = Path(file)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        if hasattr(fig, 'figure'):  # if Axes
            fig.figure.savefig(file_path)
        elif hasattr(fig, 'save'):  # for plotnine.ggplot.ggplot
            fig.save(file_path)
        else:  # if Figure 
            fig.savefig(file_path)
        return file_path
    except Exception as e:
        raise e


def set_fig_path(func, fig=None, **kwargs):
    "maybe I need to save figure by myself, instead of using scanpy save function..."
    fig_dir = Path(os.getcwd()) / "figures"

    kwargs.pop("save", None)
    kwargs.pop("show", None)
    args = []
    for k,v in kwargs.items():
        if isinstance(v, (tuple, list, set)):
            args.append(f"{k}-{'-'.join([str(i) for i in v])}")
        else:
            args.append(f"{k}-{v}")
    args_str = "_".join(args)
    if func == "rank_genes_groups_dotplot":
        old_path = fig_dir / 'dotplot_.png'
        fig_path = fig_dir / f"{func}_{args_str}.png"
    elif func in ["scatter", "embedding"]:
        if "basis" in kwargs and kwargs['basis'] is not None:
            old_path = fig_dir / f"{kwargs['basis']}.png"
            fig_path = fig_dir / f"{func}_{args_str}.png"
        else:
            old_path = fig_dir / f"{func}.png"
            fig_path = fig_dir / f"{func}_{args_str}.png"
    elif func == "highly_variable_genes":        
        old_path = fig_dir / 'filter_genes_dispersion.png'
        fig_path = fig_dir / f"{func}_{args_str}.png"
    else:
        if (fig_dir / f"{func}_.png").is_file():
            old_path = fig_dir / f"{func}_.png"
        else:
            old_path = fig_dir / f"{func}.png"
        fig_path = fig_dir / f"{func}_{args_str}.png"
    try:
        if fig is not None:
            savefig(fig, fig_path)
        else:
            os.rename(old_path, fig_path)
        return fig_path
    except FileNotFoundError:
        print(f"The file {old_path} does not exist")
    except FileExistsError:
        print(f"The file {fig_path} already exists")
    except PermissionError:
        print("You don't have permission to rename this file")

    if os.environ.get(f"{PKG}_TRANSPORT") == "stdio":
        return fig_path
    else:
        host = os.environ.get(f"{PKG}_HOST")
        port = os.environ.get(f"{PKG}_PORT")
        fig_path = f"http://{host}:{port}/figures/{Path(fig_path).name}"
        return fig_path


async def get_figure(request):
    figure_name = request.path_params["figure_name"]
    figure_path = f"./figures/{figure_name}"
    
    # 检查文件是否存在
    if not os.path.isfile(figure_path):
        return Response(content={"error": "figure not found"}, media_type="application/json")
    
    return FileResponse(figure_path)


async def forward_request(func, request, **kwargs):
    from fastmcp import Client

    forward_url = os.environ.get(f"{PKG}_FORWARD")
    request_kwargs = request.model_dump()
    request_args = request.model_fields_set
    func_kwargs = {"request": {k: request_kwargs.get(k) for k in request_args}}
    func_kwargs.update({k:v for k,v in kwargs.items() if v is not None})
    if not forward_url:
        return None
        
    client = Client(forward_url)
    async with client:
        tools = await client.list_tools()
        func = [t.name for t in tools if t.name.endswith(func)][0]
        try:
            result = await client.call_tool(func, func_kwargs)
            return result
        except Exception as e:
            raise e

def obsm2adata(adata, obsm_key):
    from anndata import AnnData

    if obsm_key not in adata.obsm_keys():
        raise ValueError(f"key {obsm_key} not found in adata.obsm")
    else:
        return AnnData(adata.obsm[obsm_key], obs=adata.obs, obsm=adata.obsm)