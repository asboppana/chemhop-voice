#!/bin/bash
# Development startup script for Tir-us backend with MCP servers

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}🚀 Starting Tir-us Development Environment${NC}"
echo ""

# Check for virtual environments
if [ ! -d ".venv" ]; then
    echo -e "${RED}❌ Backend virtual environment not found at backend/.venv${NC}"
    echo -e "${YELLOW}   Please run: python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt${NC}"
    exit 1
fi

if [ ! -d "../admet-service/.venv" ]; then
    echo -e "${RED}❌ ADMET service virtual environment not found at admet-service/.venv${NC}"
    echo -e "${YELLOW}   Please run: cd ../admet-service && python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt${NC}"
    exit 1
fi

# Activate backend virtual environment for this script
source .venv/bin/activate

# Cleanup function
cleanup() {
    echo -e "\n${YELLOW}🛑 Shutting down services...${NC}"
    kill $(jobs -p) 2>/dev/null || true
    exit 0
}

trap cleanup SIGINT SIGTERM

# Function to wait for a service to be ready
wait_for_service() {
    local port=$1
    local name=$2
    local max_attempts=30
    local attempt=0
    
    echo -e "${YELLOW}⏳ Waiting for ${name} on port ${port}...${NC}"
    while ! nc -z localhost $port 2>/dev/null; do
        attempt=$((attempt + 1))
        if [ $attempt -ge $max_attempts ]; then
            echo -e "${RED}❌ Timeout waiting for ${name}${NC}"
            return 1
        fi
        sleep 1
    done
    echo -e "${GREEN}✅ ${name} is ready${NC}"
}

# Function to get ngrok URL from API
get_ngrok_url() {
    local api_port=$1
    local max_attempts=30
    local attempt=0
    
    while [ $attempt -lt $max_attempts ]; do
        local url=$(curl -s http://localhost:${api_port}/api/tunnels 2>/dev/null | python3 -c "import sys, json; data=json.load(sys.stdin); print(data['tunnels'][0]['public_url'] if data.get('tunnels') else '')" 2>/dev/null || echo "")
        
        if [ ! -z "$url" ]; then
            echo "$url"
            return 0
        fi
        
        attempt=$((attempt + 1))
        sleep 1
    done
    
    echo ""
    return 1
}

# 1. Start MCP server (phase1 - port 8001)
echo -e "${BLUE}1️⃣  Starting Phase 1 MCP Server (port 8001)...${NC}"
# Uses backend's activated virtual environment
python3 mcp_start.py > logs/mcp_phase1.log 2>&1 &
MCP1_PID=$!
wait_for_service 8001 "Phase 1 MCP Server"

# 2. Start ADMET service (phase2 - port 8002)
echo -e "${BLUE}2️⃣  Starting Phase 2 ADMET Service (port 8002)...${NC}"
# Use admet-service's own virtual environment
../admet-service/.venv/bin/python3 ../admet-service/main.py > logs/admet_service.log 2>&1 &
ADMET_PID=$!
wait_for_service 8002 "ADMET Service"

# 3. Start ngrok for port 8001
echo -e "${BLUE}3️⃣  Starting ngrok tunnel for port 8001...${NC}"
ngrok http 8001 --log=stdout > logs/ngrok_8001.log 2>&1 &
NGROK1_PID=$!
sleep 2

# 4. Start ngrok for port 8002 (uses different web interface port)
echo -e "${BLUE}4️⃣  Starting ngrok tunnel for port 8002...${NC}"
ngrok http 8002 --log=stdout --web-addr=localhost:4041 > logs/ngrok_8002.log 2>&1 &
NGROK2_PID=$!
sleep 2

# 5. Get ngrok URLs
echo -e "${BLUE}5️⃣  Retrieving ngrok URLs...${NC}"
PHASE1_URL=$(get_ngrok_url 4040)
PHASE2_URL=$(get_ngrok_url 4041)

if [ -z "$PHASE1_URL" ] || [ -z "$PHASE2_URL" ]; then
    echo -e "${RED}❌ Failed to get ngrok URLs${NC}"
    cleanup
fi

# Convert http to https if needed
PHASE1_URL=${PHASE1_URL/http:/https:}
PHASE2_URL=${PHASE2_URL/http:/https:}

echo -e "${GREEN}✅ Phase 1 MCP Server URL: ${PHASE1_URL}/sse${NC}"
echo -e "${GREEN}✅ Phase 2 ADMET Service URL: ${PHASE2_URL}/sse${NC}"

# 6. Export environment variables and start main backend
echo -e "${BLUE}6️⃣  Starting main backend with ngrok URLs...${NC}"
export PHASE1_MCP_URL="${PHASE1_URL}/sse"
export PHASE2_MCP_URL="${PHASE2_URL}/sse"

echo ""
echo -e "${GREEN}════════════════════════════════════════════════${NC}"
echo -e "${GREEN}🎉 All services started successfully!${NC}"
echo -e "${GREEN}════════════════════════════════════════════════${NC}"
echo -e "Phase 1 MCP:    ${BLUE}${PHASE1_URL}/sse${NC}"
echo -e "Phase 2 ADMET:  ${BLUE}${PHASE2_URL}/sse${NC}"
echo -e "Main Backend:   ${BLUE}http://localhost:8000${NC}"
echo -e "${GREEN}════════════════════════════════════════════════${NC}"
echo ""
echo -e "${YELLOW}Press Ctrl+C to stop all services${NC}"
echo ""

# Start uvicorn in foreground (uses backend's activated virtual environment)
# Ctrl+C will trigger cleanup and stop all background services
uvicorn main:app --reload --host 0.0.0.0 --port 8000

