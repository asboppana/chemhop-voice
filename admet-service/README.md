# ADMET AI MCP Service

A standalone Model Context Protocol (MCP) server for ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) property predictions using the ADMET AI model.

## Purpose

This service provides isolated ADMET prediction capabilities to avoid dependency conflicts with other chemistry tools. It runs as a separate MCP server with its own UV environment and dependencies.

## Features

- Predict ADMET properties for molecules from SMILES strings
- Isolated dependency environment (separate from main backend)
- MCP protocol for clean client-server communication
- Fast inference using pre-trained ADMET AI models

## Installation

### Prerequisites

- Python 3.11 or higher
- UV package manager

### Setup

1. Navigate to the admet-service directory:
```bash
cd admet-service
```

2. Create and activate the virtual environment:
```bash
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

3. Install dependencies:
```bash
uv pip install -e .
```

## Running the Service

### Standalone Mode (stdio)

Run the MCP server in stdio mode (default):

```bash
python main.py
```

### As a Service

The service can be integrated into the main backend application or run independently. See the main backend documentation for integration details.

## MCP Tools

### `predict_admet_properties`

Predict ADMET properties for a molecule.

**Parameters:**
- `smiles` (str): SMILES string representation of the molecule

**Returns:**
- Dictionary containing:
  - `smiles`: The input SMILES string
  - `predictions`: Dictionary of predicted ADMET properties

**Example:**
```python
{
    "smiles": "CCO",
    "predictions": {
        "Caco2_Wang": 0.234,
        "Solubility_AqSolDB": -0.234,
        "HIA_Hou": 0.95,
        ...
    }
}
```

## Dependencies

Core dependencies:
- `admet-ai==1.4.0` - ADMET property prediction model
- `chemprop==1.6.1` - Chemistry property prediction framework
- `rdkit>=2025.9.3` - Chemistry toolkit
- `torch==2.5.0` - PyTorch for model inference
- `fastmcp` - MCP server framework

## Architecture

This service is designed to be completely independent from the main backend:

```
admet-service/
├── .venv/              # Isolated UV environment
├── src/
│   ├── admet_tool.py   # ADMET AI wrapper
│   └── server.py       # MCP server implementation
├── main.py             # Entry point
└── pyproject.toml      # Dependencies
```

## Development

### Testing

Run tests:
```bash
python -m pytest tests/
```

### Adding New Tools

To add new MCP tools, edit `src/server.py` and add new `@mcp.tool()` decorated functions.

## License

Part of the devTir-us project.

