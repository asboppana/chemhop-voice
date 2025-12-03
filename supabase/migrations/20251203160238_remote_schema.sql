-- Initial Schema Migration
-- This migration combines all schema files in order:
-- 1. Extensions (00_extensions.sql)
-- 2. Schema permissions (00_schemas.sql)
-- 3. Helper functions (01_functions.sql)
-- 4. Public schema tables (10_public/*.sql)

-- ============================================================================
-- EXTENSIONS
-- ============================================================================

-- UUID generation
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Crypto functions for gen_random_uuid()
CREATE EXTENSION IF NOT EXISTS "pgcrypto";

-- Full text search
CREATE EXTENSION IF NOT EXISTS "pg_trgm";

-- ============================================================================
-- SCHEMA PERMISSIONS
-- ============================================================================

-- Grant permissions on public schema
GRANT ALL ON SCHEMA public TO postgres;
GRANT USAGE ON SCHEMA public TO authenticated;
GRANT USAGE ON SCHEMA public TO anon;
GRANT USAGE ON SCHEMA public TO service_role;

-- ============================================================================
-- HELPER FUNCTIONS
-- ============================================================================

-- Function to automatically update updated_at timestamp
CREATE OR REPLACE FUNCTION fn_set_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- ============================================================================
-- PUBLIC SCHEMA TABLES
-- ============================================================================

-- User Information Table
CREATE TABLE IF NOT EXISTS public.user_information (
  id uuid DEFAULT gen_random_uuid() NOT NULL PRIMARY KEY,
  person_id uuid NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
  first_name text,
  last_name text,
  email character varying,
  admin_label character varying,
  created_at timestamp with time zone DEFAULT now() NOT NULL,
  updated_at timestamp with time zone DEFAULT now() NOT NULL
);

ALTER TABLE public.user_information OWNER TO "postgres";

-- Table comments
COMMENT ON COLUMN public.user_information.first_name IS 'First name of the user';
COMMENT ON COLUMN public.user_information.last_name IS 'Last name of the user';
COMMENT ON COLUMN public.user_information.email IS 'Email of the user';
COMMENT ON COLUMN public.user_information.admin_label IS 'Admin label for the user';

-- Constraints
ALTER TABLE ONLY "public"."user_information"
    ADD CONSTRAINT user_information_person_id_key UNIQUE (person_id);

-- Indexes
CREATE INDEX ix_user_information_id ON public.user_information USING btree (id);

-- Triggers
CREATE TRIGGER trg_user_information_updated_at 
  BEFORE UPDATE ON public.user_information 
  FOR EACH ROW 
  EXECUTE FUNCTION fn_set_updated_at();

-- Permissions
GRANT ALL ON TABLE public.user_information TO "anon";
GRANT ALL ON TABLE public.user_information TO "authenticated";
GRANT ALL ON TABLE public.user_information TO "service_role";