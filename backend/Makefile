# Makefile for Development
# Quick access to common commands

.PHONY: help dev-local dev-docker frontend

# Default target
help:
	@echo "Development Commands"
	@echo ""
	@echo "Backend:"
	@echo "  make dev-local        🚀 Start backend locally (no Docker)"
	@echo "  make dev-docker       🐳 Start backend with Docker"
	@echo "  make stop-docker      Stop Docker services"
	@echo ""
	@echo "Frontend:"
	@echo "  make frontend         🎨 Start frontend in dev mode"
	@echo ""

# ============================================================================
# Backend Development
# ============================================================================

# Local development without Docker
dev-local:
	@echo "🚀 Starting backend locally..."
	@echo "Starting FastAPI server..."
	@uvicorn main:app --reload --host 0.0.0.0 --port 8000

# Docker development
dev-docker:
	@echo "🐳 Starting backend with Docker..."
	@docker compose up -d

stop-docker:
	@echo "🛑 Stopping Docker services..."
	@docker compose down

# ============================================================================
# Frontend Development
# ============================================================================

# Start frontend development server
frontend:
	@echo "🎨 Starting frontend development server..."
	@cd ../frontend && pnpm run dev