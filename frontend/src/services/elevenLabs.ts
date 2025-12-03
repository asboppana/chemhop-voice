/**
 * ElevenLabs Voice Agent Service
 * Handles initialization and management of ElevenLabs conversational AI
 */

import server from "@/app/server";
/**
 * Fetches the voice agent ID from the backend
 */
export async function fetchVoiceAgentId(): Promise<string> {
  const response = await server.voiceAPI.getVoiceAgentId();
  console.log(response)
  return response.data.agent_id;
}

/**
 * Initializes the ElevenLabs Conversational AI widget
 * @param agentId - The ElevenLabs agent ID
 * @param containerId - The ID of the container element to mount the widget
 */
export function initializeElevenLabsWidget(
  agentId: string,
  containerId: string
): void {
  // Check if the ElevenLabs Convai object is available
  if (typeof window !== "undefined" && (window as any).Convai) {
    const Convai = (window as any).Convai;

    // Initialize the widget
    Convai.mount({
      agentId: agentId,
      element: document.getElementById(containerId),
      // Optional configuration
      // You can customize the widget appearance and behavior here
      // See: https://elevenlabs.io/docs/agents-platform/guides/quickstarts/java-script
    });
  } else {
    console.error("ElevenLabs Convai library not loaded");
  }
}

/**
 * Loads the ElevenLabs Convai script
 * @returns Promise that resolves when the script is loaded
 */
export function loadElevenLabsScript(): Promise<void> {
  return new Promise((resolve, reject) => {
    // Check if script is already loaded
    if ((window as any).Convai) {
      resolve();
      return;
    }

    // Create script element
    const script = document.createElement("script");
    script.src = "https://elevenlabs.io/convai-widget/index.js";
    script.async = true;
    script.type = "module";

    script.onload = () => resolve();
    script.onerror = () =>
      reject(new Error("Failed to load ElevenLabs Convai script"));

    document.head.appendChild(script);
  });
}

