import React, { useEffect, useRef, useState } from 'react';
import { XClose } from '@/components/ui/XClose';
import { AnimatePresence, motion } from 'framer-motion';
import { useChatContext, type Message } from '@/contexts/ChatContext';
import RotatingSedonaLogo from '@/components/animations/RotatingSedonaLogo';

// Configuration for chat panel positioning and behavior
const CHAT_CONFIG = {
  // Floating mode positioning
  FLOATING_TOP_OFFSET: 120,
  FLOATING_BOTTOM_OFFSET: 4,
  FLOATING_RIGHT_OFFSET: 8,
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
  
  // Progressive loading states
  const [loadingMessage, setLoadingMessage] = useState('Thinking...');
  
  const {
    state: { messages, mode, voiceActive, isLoading },
    shouldShowInlineChat,
    shouldShowFloatingChat,
    isMobile,
    closeChat,
    toggleMode,
    toggleVoice,
    sendMessage
  } = useChatContext();

  // Don't render if this mode isn't active
  const shouldRender = isInlineMode ? shouldShowInlineChat : shouldShowFloatingChat;

  // Progressive loading messages based on elapsed time
  useEffect(() => {
    if (!isLoading) {
      setLoadingMessage('Thinking...');
      return;
    }
    
    // Reset to initial message
    setLoadingMessage('Thinking...');
    
    // After 5 seconds, change to "Searching Literature..."
    const timer1 = setTimeout(() => {
      setLoadingMessage('Searching Literature...');
    }, 5000);
    
    // After 15 seconds, change to "Cross Checking Sources..."
    const timer2 = setTimeout(() => {
      setLoadingMessage('Cross Checking Sources...');
    }, 15000);
    
    // Cleanup timers when loading stops or component unmounts
    return () => {
      clearTimeout(timer1);
      clearTimeout(timer2);
    };
  }, [isLoading]);

  // Auto-scroll to bottom when new messages arrive or loading state changes
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages, isLoading, loadingMessage]);

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

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (!inputValue.trim() || isLoading) return;
    
    sendMessage(inputValue);
    setInputValue('');
  };

  const renderMessages = () => {
    return (
      <div className="flex-1 overflow-y-auto p-4 space-y-4 bg-[#1C1C1C]">
        {messages.length === 0 ? (
          <div className="text-center text-white py-8">
            <div className="w-12 h-12 mx-auto mb-3 text-white flex items-center justify-center">
              <img src="/icons/CoolChat.svg" alt="Sedona Health" className="w-6 h-6" style={{ filter: 'brightness(0) saturate(100%) invert(100%)' }} />
            </div>
            <p className="text-sm font-medium text-white/70">Let's start designing new molecules.</p>
          </div>
        ) : (
          messages.map((message, index) => {
            // Only show cursor for the last message when it's actively streaming
            const isLastMessage = index === messages.length - 1;
            const shouldShowCursor = isLastMessage && message.isStreaming;
            
            return (
              <div
                key={message.id}
                className={`flex gap-3 ${message.type === 'user' ? 'justify-end' : 'justify-start'}`}
              >
                {message.type === 'user' ? (
                  // User message - light blue bubble
                  <div className="max-w-[75%] px-4 py-2 rounded-2xl text-sm bg-gray-500">
                    <p className="whitespace-pre-wrap text-gray-1500 font-medium">
                      {message.content}
                      {shouldShowCursor && (
                        <span className="inline-block w-px h-4 bg-black ml-1 animate-pulse translate-y-px" />
                      )}
                    </p>
                    <p className="text-xs mt-1 text-gray-1000">
                      {message.timestamp.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}
                    </p>
                  </div>
                ) : (
                  // Assistant message - no bubble, just text
                  <div className="max-w-[85%] text-sm">
                    <p className="whitespace-pre-wrap text-white font-medium leading-relaxed">
                      {message.content}
                      {shouldShowCursor && (
                        <span className="inline-block w-px h-4 bg-white ml-1 animate-pulse translate-y-px" />
                      )}
                    </p>
                    <p className="text-xs mt-1 text-gray-500">
                      {message.timestamp.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}
                    </p>
                  </div>
                )}
              </div>
            );
          })
        )}
        
        {/* Thinking indicator - shows during backend processing */}
        {isLoading && (
          <div className="flex gap-3 justify-start">
            <div className="flex items-center gap-2 px-4 py-3 text-sm text-white">
              <RotatingSedonaLogo size={16} color="white" easing="exponential" />
              <span className="text-white/70 font-medium">{loadingMessage}</span>
            </div>
          </div>
        )}
        
        <div ref={messagesEndRef} />
      </div>
    );
  };

  const handleVoiceToggle = () => {
    toggleVoice();
  };

  const renderInput = () => (
    <div className="bg-[#1C1C1C] p-4">
      {/* Voice status indicator with stop button */}
      {voiceActive && (
        <div className="mb-2 flex items-center justify-between text-xs text-gray-750">
          <div className="flex items-center gap-2">
            <div className={`w-2 h-2 rounded-full ${
              voiceActive ? 'bg-biomarker-green' : 
              'bg-red-500'
            }`} 
            style={voiceActive ? {
              animation: 'pulse 2s ease-in-out infinite'
            } : undefined}
            />
            <span>
              {voiceActive 
                ? 'Listening...'
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
        {/* Voice toggle button */}
        <button
          type="button"
          onClick={handleVoiceToggle}
          className={`flex-shrink-0 px-3 py-2 rounded-full transition-all border border-gray-1500 ${
            voiceActive 
              ? 'bg-biomarker-green/20 text-biomarker-green hover:bg-biomarker-green/30' 
              : 'bg-[#1B1B1B] text-gray-750 hover:text-white'
          }`}
          title={voiceActive ? "Stop voice conversation" : "Start voice conversation"}
        >
          <svg className="w-4 h-4" fill="currentColor" viewBox="0 0 20 20">
            <path fillRule="evenodd" d="M7 4a3 3 0 016 0v4a3 3 0 11-6 0V4zm4 10.93A7.001 7.001 0 0017 8a1 1 0 10-2 0A5 5 0 015 8a1 1 0 00-2 0 7.001 7.001 0 006 6.93V17H6a1 1 0 100 2h8a1 1 0 100-2h-3v-2.07z" clipRule="evenodd" />
          </svg>
        </button>

        {/* Text input - always enabled */}
        <input
          type="text"
          value={inputValue}
          onChange={(e) => setInputValue(e.target.value)}
          placeholder={voiceActive ? "Type or speak your message..." : "Type your message..."}
          className="flex-1 px-3 py-2 bg-[#1B1B1B] text-white placeholder-gray-750 rounded-full border border-gray-1500 text-sm focus:outline-none focus:border-gray-1000"
          disabled={isLoading}
        />

        {/* Send button - only show when there's text to send */}
        {inputValue.trim() && (
          <button
            type="submit"
            disabled={isLoading}
            className="flex-shrink-0 px-3 py-2 rounded-full transition-all bg-white text-black hover:bg-gray-200 disabled:opacity-50 disabled:cursor-not-allowed"
            title="Send message"
          >
            <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 19l9 2-9-18-9 18 9-2zm0 0v-8" />
            </svg>
          </button>
        )}
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
              <img src="/icons/LeftPanel.svg" alt="Attach" className="mr-2 w-4 h-4 invert" />
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
        <div className="h-full bg-black flex flex-col overflow-hidden rounded-tl-[24px] rounded-tr-[24px] rounded-br-[24px]">
          {renderHeader()}
          {renderMessages()}
          {renderInput()}
        </div>
      );
    }

    return (
      <motion.div 
        className="h-full bg-black flex flex-col overflow-hidden rounded-tl-[24px] rounded-tr-[24px] rounded-br-[24px]"
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
          className="fixed inset-0 bg-black bg-opacity-20 z-40"
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
