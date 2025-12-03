CREATE TABLE IF NOT EXISTS public.user_information (
  id uuid DEFAULT uuid_generate_v4() NOT NULL PRIMARY KEY,
  person_id uuid NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
  first_name text,
  last_name text,
  email character varying,
  admin_label character varying,
  created_at timestamp with time zone DEFAULT now() NOT NULL,
  updated_at timestamp with time zone DEFAULT now() NOT NULL
);

ALTER TABLE public.user_information OWNER TO "postgres";

COMMENT ON COLUMN public.user_information.first_name IS 'First name of the user';
COMMENT ON COLUMN public.user_information.last_name IS 'Last name of the user';
COMMENT ON COLUMN public.user_information.email IS 'Email of the user';
COMMENT ON COLUMN public.user_information.admin_label IS 'Admin label for the user';

ALTER TABLE ONLY "public"."user_information"
    ADD CONSTRAINT user_information_person_id_key UNIQUE (person_id);

CREATE INDEX ix_user_information_id ON public.user_information USING btree (id);


CREATE TRIGGER trg_user_information_updated_at BEFORE UPDATE ON public.user_information FOR EACH ROW EXECUTE FUNCTION fn_set_updated_at();
GRANT ALL ON TABLE public.user_information TO "anon";
GRANT ALL ON TABLE public.user_information TO "authenticated";
GRANT ALL ON TABLE public.user_information TO "service_role";
