/**
 * ElevenLabs Voice Agent Service
 * Simple service to fetch agent ID from backend
 */

import server from "@/app/server";

/**
 * Fetches the voice agent ID from the backend
 */
export async function fetchVoiceAgentId(): Promise<string> {
  const response = await server.voiceAPI.getVoiceAgentId();
  console.log('Voice Agent ID:', response.data.agent_id);
  return response.data.agent_id;
}

