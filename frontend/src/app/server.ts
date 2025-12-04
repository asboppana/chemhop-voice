import axios from 'axios';
import { env } from './env';

const BASE_URL = env.VITE_API_BASE_URL;
// Configurable API prefix (defaults to /api/v1)
const API_PREFIX = '/api/v1';

// Final API URL with prefix
const API_BASE_URL = `${BASE_URL}${API_PREFIX}`;

const server = axios.create({
  baseURL: API_BASE_URL,
});


server.interceptors.request.use((config) => {
  const token = localStorage.getItem('authToken');
  if (token) {
    config.headers.Authorization = `Bearer ${token}`;
  }
  
  // Add impersonation header if active, but NEVER for admin messaging endpoints
  // Admin messaging should always run as the admin, not as an impersonated user
  const urlStr = (config.url || '').toLowerCase();
  const isAdminMessagesEndpoint = urlStr.includes('/admin/messages');
  if (!isAdminMessagesEndpoint) {
    const impersonationStr = sessionStorage.getItem('impersonation');
    if (impersonationStr) {
      try {
        const impersonation = JSON.parse(impersonationStr);
        if (impersonation.isImpersonating && impersonation.targetUserId) {
          config.headers['X-Target-User-ID'] = impersonation.targetUserId;

          // Optional: Log for debugging (can be removed in production)
          if (config.url && !config.url.includes('/bulk-connection-status')) {
            console.log(`🎭 API call with impersonation: ${config.method?.toUpperCase()} ${config.url} → ${impersonation.targetUserId}`);
          }
        }
      } catch (error) {
        console.error('Failed to parse impersonation state:', error);
        sessionStorage.removeItem('impersonation');
      }
    }
  }
  
  return config;
});

// -----------------------------------------------------------
// Global auth invalidation handling (JWT/sub-claim issues)
// -----------------------------------------------------------
function isAuthInvalidError(error: any): boolean {
  const status = error?.response?.status;
  const data = error?.response?.data;
  const msg = typeof data === 'string' ? data : (data?.message || data?.error || '');

  if (status === 401 || status === 403) return true;

  // Be specific about 500s we treat as auth invalid
  if (status === 500 && /sub claim.*does not exist|jwt|token|authorization/i.test(String(msg))) {
    return true;
  }

  return false;
}

function isPublicPath(pathname: string): boolean {
  const publics = [
    /^\/$/,
    /^\/login(\/.*)?$/,
    /^\/forgot-password$/,
    /^\/reset-password$/,
    /^\/auth\/confirm(\/.*)?$/,
    /^\/invite\/[^/]+(\/(payment|wearables))?$/,
    /^\/unlock$/,
    /^\/intro-scroll$/,
    /^\/welcomemerged$/,
    /^\/connect-devices$/,
    /^\/connect-ehr$/,
    /^\/admin-invite\/[^/]+$/,
    /^\/twilio-verify$/,
    /^\/surveys\/magic\/[^/]+$/,
    /^\/raffle$/,
  ];
  return publics.some((rx) => rx.test(pathname));
}

function clearAuthSideEffects() {
  try { sessionStorage.removeItem('impersonation'); } catch {}
  try { localStorage.removeItem('authToken'); } catch {}
  try { localStorage.removeItem('refreshToken'); } catch {}
  try { localStorage.removeItem('userInfo'); } catch {}
}

function makeAuthInvalidHandler(client: any) {
  return async (error: any) => {
    if (!isAuthInvalidError(error)) return Promise.reject(error);

    const original = error.config || {};
    
    // If this is a 401 and we haven't tried refreshing yet, attempt token refresh
    if (error.response?.status === 401 && !original._retriedRefresh) {
      try {
        // No refresh needed for API tokens
        return Promise.reject(error);
      } catch (refreshErr) {
        console.error('Token refresh failed:', refreshErr);
        // Fall through to logout logic
      }
    }
    
    // If refresh failed or this wasn't a 401, proceed with logout
    if (original._retriedAuthClear) return Promise.reject(error);

    clearAuthSideEffects();

    const onPublic = isPublicPath(window.location.pathname);
    if (onPublic) {
      // Retry once without token
      original._retriedAuthClear = true;
      if (original.headers) delete original.headers.Authorization;
      return client.request(original);
    } else {
      // Sign out from your application (if applicable)
      clearAuthSideEffects();
      const redirect = encodeURIComponent(window.location.pathname + window.location.search);
      window.location.href = `/login?redirect=${redirect}`;
      // Prevent further promise chains; navigation will take over
      return new Promise(() => {});
    }
  };
}

server.interceptors.response.use(
  (response: any) => response,
  makeAuthInvalidHandler(server)
);


export const authAPI = {
  signIn: (data: { email: string; password: string }) => 
    server.post('/signin', data),  // ← Updated: no /users prefix
  
  signUp: (user_data: { 
    email: string; 
    password: string; 
    first_name?: string; 
    last_name?: string; 
    invite_token?: string;
    phone_number?: string;
    street_line1?: string;
    street_line2?: string;
    city?: string;
    state?: string;
    zip_code?: string;
    country?: string;
    date_of_birth?: string;
    gender?: string;
  }) => server.post('/signup', user_data),  // ← Updated
}

const llmAPI = {  
  stream: (data: { messages: { role: 'user' | 'assistant', content: string }[], model: string }) => server.post('/chat/general', data).then(res => res.data),
};

const voiceAPI = {
  callChatAgent: (data: { messages: { role: 'user' | 'assistant', content: string }[], model?: string }) => server.post('/chat', data).then(res => res.data),
};

export default {
  ...server,
  authAPI,
  llmAPI,
  voiceAPI,
};
