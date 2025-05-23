# decoupler-MCP

Natural language interface for scRNA-Seq analysis with decoupler through MCP.

## 🪩 What can it do?

- IO module like read and write scRNA-Seq data
- Pathway activity/Transcription factor  inference
- Tool module, like clustering, differential expression etc.
- Plotting module, like violin, umap/tsne

## ❓ Who is this for?

- Anyone who wants to do scRNA-Seq analysis natural language!
- Agent developers who want to call decoupler's functions for their applications

## 🌐 Where to use it?

You can use decoupler-mcp in most AI clients, plugins, or agent frameworks that support the MCP:

- AI clients, like Cherry Studio
- Plugins, like Cline
- Agent frameworks, like Agno 

## 🎬 Demo

A demo showing scRNA-Seq cell cluster analysis in a AI client Cherry Studio using natural language based on decoupler-mcp


## 🏎️ Quickstart

### Install

Install from PyPI
```
pip install decoupler-mcp
```
you can test it by running
```
decoupler-mcp run
```

#### run scnapy-server locally
Refer to the following configuration in your MCP client:

```
"mcpServers": {
  "decoupler-mcp": {
    "command": "decoupler-mcp",
    "args": [
      "run"
    ]
  }
}
```

#### run scnapy-server remotely
Refer to the following configuration in your MCP client:

run it in your server
```
decoupler-mcp run --transport shttp --port 8000
```

Then configure your MCP client, like this:
```
http://localhost:8000/mcp
```

## 🤝 Contributing

If you have any questions, welcome to submit an issue, or contact me(hsh-me@outlook.com). Contributions to the code are also welcome!

## Citing

If you use decoupler-mcp in for your research, please consider citing  following work: 
> Badia-i-Mompel P., Vélez Santiago J., Braunger J., Geiss C., Dimitrov D., Müller-Dott S., Taus P., Dugourd A., Holland C.H., Ramirez Flores R.O. and Saez-Rodriguez J. 2022. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances. https://doi.org/10.1093/bioadv/vbac016

