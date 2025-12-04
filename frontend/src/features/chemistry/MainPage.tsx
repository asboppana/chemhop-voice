import React, { useEffect, useState } from 'react';
import FloatingNullState from '@/components/animations/FloatingNullState';
import { useChatContext } from '@/contexts/ChatContext';

export const ChemistryMainPage: React.FC = () => {
  const { 
    state: { voiceActive, voiceConnectionStatus, voiceAgentStatus, messages, isOpen },
    startVoice,
    openChat
  } = useChatContext();
  
  const [error, setError] = useState<string | null>(null);
  const [initializing, setInitializing] = useState(false);

  useEffect(() => {
    // Auto-start voice conversation on mount if not already active
    // Note: Chat panel will auto-open when first message is received
    const initVoice = async () => {
      if (!voiceActive && !initializing) {
        setInitializing(true);
        const started = await startVoice();
        if (!started) {
          setError('Failed to start voice conversation. Please try again.');
        }
        setInitializing(false);
      }
    };

    // Small delay to ensure chat context is ready
    const timer = setTimeout(initVoice, 500);
    
    // DON'T cleanup voice on unmount - voice should persist across navigation
    return () => {
      clearTimeout(timer);
      // Voice conversation continues even if MainPage unmounts
    };
  }, []); // Only run once on mount

  const isConnected = voiceConnectionStatus === 'connected';
  const isSpeaking = voiceAgentStatus === 'speaking';
  
  // Check if there's actual conversation beyond the initial system greeting
  // Hide animation if there are user messages or multiple AI messages
  const hasUserMessages = messages.some(m => m.type === 'user');
  const hasConversation = hasUserMessages || messages.length > 1;
  const shouldShowAnimation = !error && !hasConversation;

  // Auto-open chat panel when conversation starts
  useEffect(() => {
    if (hasConversation && !isOpen) {
      console.log('🎯 Conversation detected, auto-opening chat panel');
      openChat();
    }
  }, [hasConversation, isOpen, openChat]);

  return (
    <div className="max-h-[90vh] bg-white flex flex-col items-center justify-center p-8">
      <div className="w-full max-w-4xl flex flex-col items-center">
        {/* Show animation only until conversation starts */}
        {shouldShowAnimation && (
          <div className="relative flex justify-center">
            <div className="flex justify-center">
              <FloatingNullState 
                hideIcons={isSpeaking}
                timing={
                  isSpeaking
                    ? {
                        pushDurationMs: 800,
                        relaxDurationMs: 500,
                        itemHoldMs: 200,
                        idleMinMs: 600,
                        idleJitterMs: 300,
                      }
                    : undefined
                }
              />
            </div>
            {/* Pulsing green dot in center when connected */}
            {isConnected && (
              <div className="absolute inset-0 flex items-center justify-center pointer-events-none">
                {/* Outer glow - larger, subtle */}
                <div 
                  className="absolute w-5 h-5 rounded-full bg-biomarker-green/20"
                  style={{
                    animation: 'pulse-glow 2s ease-in-out infinite',
                  }}
                />
                {/* Inner core - bright green */}
                <div 
                  className="absolute w-2 h-2 rounded-full bg-biomarker-green shadow-lg"
                  style={{
                    animation: 'pulse-core 2s ease-in-out infinite',
                    boxShadow: '0 0 8px rgba(34, 197, 94, 0.6)',
                  }}
                />
              </div>
            )}
          </div>
        )}
      </div>

      {/* Pulse animations */}
      <style>{`
        @keyframes pulse-glow {
          0%, 100% {
            transform: scale(1);
            opacity: 0.4;
          }
          50% {
            transform: scale(1.4);
            opacity: 0.2;
          }
        }
        
        @keyframes pulse-core {
          0%, 100% {
            transform: scale(1);
          }
          50% {
            transform: scale(1.2);
          }
        }
      `}</style>
    </div>
  );
};

export default ChemistryMainPage;
