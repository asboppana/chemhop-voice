import React, { useEffect, useState } from 'react';
import {
  fetchVoiceAgentId,
  loadElevenLabsScript,
  initializeElevenLabsWidget,
} from '@/services/elevenLabs';
import RotatingSedonaLogo from '@/components/animations/RotatingSedonaLogo';

export const ChemistryMainPage: React.FC = () => {
  const [agentId, setAgentId] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let mounted = true;

    const initializeVoiceAgent = async () => {
      try {
        // Load the ElevenLabs script
        await loadElevenLabsScript();

        // Fetch the agent ID from backend
        const id = await fetchVoiceAgentId();

        if (!mounted) return;

        setAgentId(id);
        setLoading(false);

        // Initialize the widget after a short delay to ensure DOM is ready
        setTimeout(() => {
          if (mounted && id) {
            initializeElevenLabsWidget(id, 'elevenlabs-widget-container');
          }
        }, 100);
      } catch (err) {
        if (!mounted) return;
        console.error('Failed to initialize voice agent:', err);
        setError(
          err instanceof Error ? err.message : 'Failed to load voice agent'
        );
        setLoading(false);
      }
    };

    initializeVoiceAgent();

    return () => {
      mounted = false;
    };
  }, []);

  return (
    <div className="min-h-screen bg-white flex items-center justify-center p-8">
      <div className="w-full max-w-4xl">
        {loading && (
          <div className="text-center">
            <RotatingSedonaLogo size={60} />
            <p className="mt-4 text-gray-1000">Loading voice agent...</p>
          </div>
        )}

        {error && (
          <div className="bg-red-50 border border-red-200 rounded-lg p-4 text-center">
            <p className="text-red-800">Error: {error}</p>
            <button
              onClick={() => window.location.reload()}
              className="mt-4 px-4 py-2 bg-red-600 text-white rounded hover:bg-red-700 transition"
            >
              Retry
            </button>
          </div>
        )}

        {!loading && !error && agentId && (
          <div className="text-center">
            <h1 className="text-3xl font-bold mb-8 text-gray-900">
              Drug Discovery Assistant
            </h1>
            <div
              id="elevenlabs-widget-container"
              className="w-full min-h-[500px] flex items-center justify-center bg-gray-50 rounded-lg border border-gray-200"
            >
              {/* ElevenLabs widget will be mounted here */}
            </div>
            <p className="mt-4 text-gray-600 text-sm">
              Start a conversation with our AI-powered drug discovery assistant
            </p>
          </div>
        )}
      </div>
    </div>
  );
};

export default ChemistryMainPage;
