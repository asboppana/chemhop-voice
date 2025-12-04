
import React, { useEffect, useRef, useState } from 'react';
import { XClose } from '@/components/ui/XClose';
import { AnimatePresence, motion } from 'framer-motion';

import RotatingSedonaLogo from '@/components/animations/RotatingSedonaLogo';
import { useChatContext, type Message } from '@/contexts/ChatContext';


// Configuration for chat panel positioning and behavior
const CHAT_CONFIG = {
  // Floating mode positioning
  FLOATING_TOP_OFFSET: 120,
  FLOATING_BOTTOM_OFFSET: 4,
  FLOATING_RIGHT_OFFSET: 6,
  FLOATING_WIDTH: 384, // w-96 = 24rem = 384px
  
  // Resize constraints
  MIN_WIDTH: 300,
  MAX_WIDTH: 600,
  MIN_HEIGHT: 400,
  
  // Inline mode - should align with BottomDrawer states
  DRAWER_EXPANDED_Y: 100,      // Matches DRAWER_CONFIG.EXPANDED_Y
  DRAWER_HALF_Y_PERCENT: 0.55, // Matches DRAWER_CONFIG.HALF_Y_PERCENT  
  DRAWER_COLLAPSED_HEIGHT: 100, // Matches DRAWER_CONFIG.COLLAPSED_VISIBLE_HEIGHT
  DRAWER_BAR_HEIGHT: 32,      // Height of the drawer's drag bar
};

interface ChatPanelProps {
  isInlineMode?: boolean;
  chatMode?: 'floating-only' | 'inline-allowed';
  // When true, parent controls slide/position; suppress internal inline animations
  parentControlsInlineAnimation?: boolean;
}

export const ChatPanel: React.FC<ChatPanelProps> = ({ 
  isInlineMode = false, 
  chatMode = 'inline-allowed',
  parentControlsInlineAnimation = false,
}) => {
  const messagesEndRef = useRef<HTMLDivElement>(null);
  const [inputValue, setInputValue] = useState('');
  const [isResizing, setIsResizing] = useState(false);
  const [panelSize, setPanelSize] = useState({
    width: CHAT_CONFIG.FLOATING_WIDTH,
    height: window.innerHeight - CHAT_CONFIG.FLOATING_TOP_OFFSET - CHAT_CONFIG.FLOATING_BOTTOM_OFFSET
  });
  
  const {
    state: { messages, isLoading, mode, voiceActive, voiceConnectionStatus, voiceAgentStatus },
    shouldShowInlineChat,
    shouldShowFloatingChat,
    isMobile,
    closeChat,
    toggleMode,
    addMessage,
    toggleVoice,
    sendTextToVoice,
    startVoice
  } = useChatContext();

  // Don't render if this mode isn't active
  const shouldRender = isInlineMode ? shouldShowInlineChat : shouldShowFloatingChat;

  // Auto-scroll to bottom when new messages arrive
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  if (!shouldRender) return null;

  const handleResizeStart = (e: React.MouseEvent) => {
    e.preventDefault();
    setIsResizing(true);
    
    const startX = e.clientX;
    const startY = e.clientY;
    const startWidth = panelSize.width;
    const startHeight = panelSize.height;
    
    const handleMouseMove = (e: MouseEvent) => {
      const deltaX = startX - e.clientX; // Inverted because we're dragging from top-left
      const deltaY = startY - e.clientY; // Inverted because we want top to move, not bottom
      
      const newWidth = Math.max(
        CHAT_CONFIG.MIN_WIDTH,
        Math.min(CHAT_CONFIG.MAX_WIDTH, startWidth + deltaX)
      );
      
      const maxHeight = window.innerHeight - CHAT_CONFIG.FLOATING_TOP_OFFSET - CHAT_CONFIG.FLOATING_BOTTOM_OFFSET;
      const newHeight = Math.max(
        CHAT_CONFIG.MIN_HEIGHT,
        Math.min(maxHeight, startHeight + deltaY)
      );
      
      setPanelSize({ width: newWidth, height: newHeight });
    };
    
    const handleMouseUp = () => {
      setIsResizing(false);
      document.removeEventListener('mousemove', handleMouseMove);
      document.removeEventListener('mouseup', handleMouseUp);
    };
    
    document.addEventListener('mousemove', handleMouseMove);
    document.addEventListener('mouseup', handleMouseUp);
  };

  const sendMessage = async (content: string) => {
    if (!content.trim() || isLoading) return;

    // Auto-start voice if not active
    if (!voiceActive) {
      console.log('Voice not active, starting voice connection...');
      await startVoice();
      // Give it a moment to connect before sending
      await new Promise(resolve => setTimeout(resolve, 100));
    }

    // Always use ElevenLabs agent for all messages
    // sendTextToVoice will handle adding the user message to chat
    const sent = sendTextToVoice(content.trim());
    
    if (sent) {
      setInputValue('');
      // ElevenLabs will handle the response through voice callbacks
    } else {
      console.error('Failed to send text to voice agent');
      // Add error message
      const errorMessage: Message = {
        id: Date.now().toString(),
        type: 'assistant',
        content: 'Sorry, I encountered an error connecting to the voice agent. Please try again.',
        timestamp: new Date()
      };
      addMessage(errorMessage);
    }
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    sendMessage(inputValue);
  };

  const handleKeyPress = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      sendMessage(inputValue);
    }
  };

  const renderMessages = () => (
    <div className={`flex-1 overflow-y-auto p-4 space-y-4 bg-[#1C1C1C]`}>
      {messages.length === 0 ? (
        <div className="text-center text-white py-8">
          <div className="w-12 h-12 mx-auto mb-3 text-white flex items-center justify-center">
            <img src="/icons/CoolChat.svg" alt="Sedona Health" className="w-6 h-6" style={{ filter: 'brightness(0) saturate(100%) invert(100%)' }} />
          </div>
          <p className="text-sm font-medium text-white/70">Ask a question about your health data.</p>
        </div>
      ) : (
        messages.map((message) => (
          <div
            key={message.id}
            className={`flex gap-3 ${message.type === 'user' ? 'justify-end' : 'justify-start'}`}
          >
            {message.type === 'user' ? (
              // User message - light blue bubble
              <div className="max-w-[75%] px-4 py-2 rounded-2xl text-sm" style={{ backgroundColor: 'rgb(232, 243, 254)' }}>
                {message.isVoice && (
                  <div className="flex items-center gap-1 mb-1 text-xs" style={{ color: 'rgb(59, 130, 246)' }}>
                    <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
                      <path fillRule="evenodd" d="M7 4a3 3 0 016 0v4a3 3 0 11-6 0V4zm4 10.93A7.001 7.001 0 0017 8a1 1 0 10-2 0A5 5 0 015 8a1 1 0 00-2 0 7.001 7.001 0 006 6.93V17H6a1 1 0 100 2h8a1 1 0 100-2h-3v-2.07z" clipRule="evenodd" />
                    </svg>
                    <span>Voice</span>
                  </div>
                )}
                <p className="whitespace-pre-wrap text-gray-1500 font-medium">{message.content}</p>
                <p className="text-xs mt-1 text-gray-1000">
                  {message.timestamp.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}
                </p>
              </div>
            ) : (
              // Assistant message - no bubble, just text
              <div className="max-w-[85%] text-sm">
                {message.isVoice && (
                  <div className="flex items-center gap-1 mb-1 text-xs text-gray-500">
                    <svg className="w-3 h-3" fill="currentColor" viewBox="0 0 20 20">
                      <path fillRule="evenodd" d="M7 4a3 3 0 016 0v4a3 3 0 11-6 0V4zm4 10.93A7.001 7.001 0 0017 8a1 1 0 10-2 0A5 5 0 015 8a1 1 0 00-2 0 7.001 7.001 0 006 6.93V17H6a1 1 0 100 2h8a1 1 0 100-2h-3v-2.07z" clipRule="evenodd" />
                    </svg>
                    <span>Voice</span>
                  </div>
                )}
                <p className="whitespace-pre-wrap text-white font-medium leading-relaxed">{message.content}</p>
                <p className="text-xs mt-1 text-gray-500">
                  {message.timestamp.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}
                </p>
              </div>
            )}
          </div>
        ))
      )}
      
      {/* Loading indicator */}
      {isLoading && (
        <div className="flex gap-3">
          <div className="px-3 py-2 rounded-lg flex items-center gap-2">
            <RotatingSedonaLogo size={16} color="white" />
            <span className="text-white text-sm">Thinking...</span>
          </div>
        </div>
      )}
      
      <div ref={messagesEndRef} />
    </div>
  );

  const handleVoiceToggle = async () => {
    await toggleVoice();
  };

  const renderInput = () => (
    <div className="bg-[#1C1C1C] p-4">
      {/* Voice status indicator with stop button */}
      {voiceActive && (
        <div className="mb-2 flex items-center justify-between text-xs text-gray-750">
          <div className="flex items-center gap-2">
            <div className={`w-2 h-2 rounded-full ${
              voiceConnectionStatus === 'connected' ? 'bg-biomarker-green' : 
              voiceConnectionStatus === 'connecting' ? 'bg-yellow-500' : 
              'bg-red-500'
            }`} 
            style={voiceConnectionStatus === 'connected' ? {
              animation: 'pulse 2s ease-in-out infinite'
            } : undefined}
            />
            <span>
              {voiceConnectionStatus === 'connected' 
                ? (voiceAgentStatus === 'speaking' ? 'AI Speaking...' : 'Listening...')
                : voiceConnectionStatus === 'connecting' 
                  ? 'Connecting...' 
                  : 'Disconnected'}
            </span>
          </div>
          
          {/* Stop voice button */}
          <button
            onClick={handleVoiceToggle}
            className="px-2 py-1 text-xs bg-red-500/20 hover:bg-red-500/30 text-red-400 rounded-full transition-colors"
            title="Stop voice conversation"
          >
            Stop Voice
          </button>
        </div>
      )}
      
      <form onSubmit={handleSubmit} className="flex gap-2">
        {/* Voice toggle button - only show when voice is NOT active */}
        {!voiceActive && (
          <button
            type="button"
            onClick={handleVoiceToggle}
            className="flex-shrink-0 px-3 py-2 rounded-full transition-all bg-[#1B1B1B] text-gray-750 hover:text-white border border-gray-1500"
            title="Start voice conversation"
          >
            <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
              <path fillRule="evenodd" d="M7 4a3 3 0 016 0v4a3 3 0 11-6 0V4zm4 10.93A7.001 7.001 0 0017 8a1 1 0 10-2 0A5 5 0 015 8a1 1 0 00-2 0 7.001 7.001 0 006 6.93V17H6a1 1 0 100 2h8a1 1 0 100-2h-3v-2.07z" clipRule="evenodd" />
            </svg>
          </button>
        )}

        {/* Text input - works alongside voice mode */}
        <input
          type="text"
          value={inputValue}
          onChange={(e) => setInputValue(e.target.value)}
          onKeyPress={handleKeyPress}
          placeholder={voiceActive ? "Type or speak your message..." : "Ask another question"}
          className="flex-1 px-3 py-2 bg-[#1B1B1B] text-white placeholder-gray-750 rounded-full border border-gray-1500 focus:outline-none focus:ring-1 focus:ring-gray-1000 text-sm disabled:opacity-50"
          disabled={isLoading}
        />
        <button
          type="submit"
          disabled={!inputValue.trim() || isLoading}
          className="px-3 py-2 text-white flex items-center justify-center disabled:opacity-50"
        >
          <img src="/icons/Arrow.svg" alt="Send" className="w-6 h-6 transform rotate-180" style={{ filter: 'brightness(0) saturate(100%) invert(100%)' }} />
        </button>
      </form>
    </div>
  );

  const renderHeader = () => (
    <div className={`flex items-center justify-between p-4 bg-[#1C1C1C] ${isInlineMode ? 'rounded-tl-3xl' : 'rounded-t-xl'} relative`}>
      {/* Resize handle - only show in floating mode */}
      {!isInlineMode && (
        <div
          className={`absolute top-0 left-0 w-4 h-4 cursor-nw-resize ${
            isResizing ? 'bg-gray-600' : 'hover:bg-gray-700'
          } transition-colors rounded-tl-xl`}
          onMouseDown={handleResizeStart}
          title="Drag to resize"
        />
      )}

      <div className="flex items-center justify-between gap-2 flex-1">
        {/* Detach/Attach button - only show on desktop and when inline mode is allowed */}
        {!isMobile && chatMode === 'inline-allowed' && (
          <button
            onClick={toggleMode}
            className="text-gray-750 hover:text-white transition-colors p-1"
            title={mode === 'attached' ? 'Detach to floating window' : 'Attach to sidebar'}
          >
            {mode === 'attached' ? (
              <img src="/icons/ArrowBack.svg" alt="Detach" className="mr-2 w-4 h-4 invert" />
            ) : (
              <img src="/icons/LeftPanel.svg" alt="Attach" className="mr-2 w-5 h-5 invert" />
            )}
          </button>
        )}
        <XClose 
          onClick={closeChat}
          className="ml-auto" 
          invertIcon={false}
          size="md"
        />
      </div>
    </div>
  );

  // Inline mode (attached) - overlay on the right side
  if (isInlineMode) {
    // When parent is orchestrating the inline animation, render without internal motion
    if (parentControlsInlineAnimation) {
      return (
        <div className="h-full bg-black flex flex-col overflow-hidden rounded-tl-[24px] rounded-bl-[24px]">
          {renderHeader()}
          {renderMessages()}
          {renderInput()}
        </div>
      );
    }

    return (
      <motion.div 
        className="h-full bg-black flex flex-col overflow-hidden rounded-tl-[24px] rounded-bl-[24px]"
        initial={{ x: '100%', opacity: 0 }}
        animate={{ x: 0, opacity: 1 }}
        exit={{ x: '100%', opacity: 0 }}
        transition={{
          duration: 0.3,
          ease: 'linear'
        }}
      >
        {renderHeader()}
        {renderMessages()}
        {renderInput()}
      </motion.div>
    );
  }

  // Floating mode (detached or mobile)
  return (
    <AnimatePresence mode="wait">
      <>
        {/* Backdrop */}
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
          onClick={closeChat}
          className="fixed inset-0 bg-black bg-opacity-20 z-50"
          transition={{
            duration: 0.3,
            ease: 'linear'
          }}
        />
        
        {/* Panel */}
        <motion.div
          initial={{ 
            x: '100%',
            opacity: 0
          }}
          animate={{ 
            x: 0,
            opacity: 1
          }}
          exit={{ 
            x: '100%',
            opacity: 0
          }}
          transition={{
            duration: 0.3,
            ease: 'linear'
          }}
          className="fixed bg-[#1C1C1C] z-50 flex flex-col rounded-[32px] overflow-hidden"
          style={{ 
            bottom: `${CHAT_CONFIG.FLOATING_BOTTOM_OFFSET}px`,
            right: `${CHAT_CONFIG.FLOATING_RIGHT_OFFSET}px`,
            width: `${panelSize.width}px`,
            height: `${panelSize.height}px`,
            maxWidth: 'calc(100vw - 2rem)',
            userSelect: isResizing ? 'none' : 'auto'
          }}
        >
          {renderHeader()}
          {renderMessages()}
          {renderInput()}
        </motion.div>
      </>
    </AnimatePresence>
  );
};

// For backward compatibility, export the ChatMessage type
export type { Message };
