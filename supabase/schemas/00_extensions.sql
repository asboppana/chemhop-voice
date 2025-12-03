-- PostgreSQL Extensions

-- UUID generation
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Crypto functions for gen_random_uuid()
CREATE EXTENSION IF NOT EXISTS "pgcrypto";

-- Full text search
CREATE EXTENSION IF NOT EXISTS "pg_trgm";
