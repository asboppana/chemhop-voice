"""
ADMET AI MCP Server Entry Point.

Run this script to start the ADMET AI MCP server:
    python main.py

The server will run in stdio mode by default, communicating via standard input/output.
"""
from src.server import run

def main():
    run()

if __name__ == "__main__":
    main()

