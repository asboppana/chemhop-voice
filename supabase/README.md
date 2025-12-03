# Supabase Local Development

This directory contains all configuration and scripts for running Supabase locally.

## 🚀 Quick Start

```bash
# First time setup (resets database, runs all migrations, verifies schemas)
./supabase/verify.sh
```

## 📁 Directory Structure

```
supabase/
├── config.toml              # Supabase configuration (ports, schemas, auth, etc.)
├── seed.sql                 # Database seed data
├── migrations/              # SQL migrations (auto-applied on db reset)
│   ├── 20251027211017_initial_setup.sql
│   ├── 20251027211330_add_comms_schema.sql
│   ├── 20251027211745_add_nutrition_schema.sql
│   └── ... (more migrations)
├── schemas/                 # Schema definitions (alternative to migrations)
│   ├── 00_extensions.sql
│   ├── 00_schemas.sql
│   ├── 01_enums.sql
│   └── ... (schema files)
├── seed/                    # Seed scripts and data
│   ├── scripts/
│   │   ├── seed_local.py    # Main seed script
│   │   ├── create_user.py   # Create test users
│   │   └── upsert_survey.py # Load survey data
│   └── data/
│       ├── clinical_survey.json
│       └── onboarding_survey.json
└── README.md               # This file
```
