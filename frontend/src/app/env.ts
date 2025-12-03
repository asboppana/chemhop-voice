import { createEnv } from "@t3-oss/env-core";
import { z } from "zod";

export const env = createEnv({
  clientPrefix: "VITE_",

  client: {
    // API
    VITE_API_BASE_URL: z.string().url(),
    VITE_API_PREFIX: z.string(),
		VITE_SUPABASE_URL: z.string().url(),
		VITE_SUPABASE_ANON_KEY: z.string(),
	},
	runtimeEnv: import.meta.env,
	emptyStringAsUndefined: true,
});
