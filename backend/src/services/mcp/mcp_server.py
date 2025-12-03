"""
Simple MCP Server using the official SDK.
Run with: python mcp_server.py
"""

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent
import asyncio
# server.py
from fastmcp import FastMCP

mcp_server = FastMCP("drugdiscoveryagent_mcp")

@mcp_server.tool
def add(a: int, b: int) -> int:
    """Add two numbers"""
    return a + b

if __name__ == "__main__":
    mcp_server.run()