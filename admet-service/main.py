from src.server import mcp


if __name__ == "__main__":
    # Use SSE transport for OpenAI MCP integration
    # OpenAI expects Server-Sent Events (SSE) protocol
    # Note: SSE creates its own routes, don't specify custom path
    mcp.run(transport="sse", host="0.0.0.0", port=8002)
    # mcp.run()

