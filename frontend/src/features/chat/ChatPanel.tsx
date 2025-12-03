
import React, { useEffect, useRef, useState } from 'react';
import { XClose } from '@/components/ui/XClose';
import { AnimatePresence, motion } from 'framer-motion';

import server  from '@/app/server';

import RotatingSedonaLogo from '@/components/animations/RotatingSedonaLogo';
import { useChatContext, type Message } from '@/contexts/ChatContext';
import { CustomIcon } from '@/components/ui/CustomIcon';


// Configuration for chat panel positioning and behavior
const CHAT_CONFIG = {
  // Floating mode positioning
  FLOATING_TOP_OFFSET: 120,
  FLOATING_BOTTOM_OFFSET: 4,
  FLOATING_RIGHT_OFFSET: 4,
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
  const [dataContext, setDataContext] = useState<string>('');
  const [inputValue, setInputValue] = useState('');
  const [isResizing, setIsResizing] = useState(false);
  const [panelSize, setPanelSize] = useState({
    width: CHAT_CONFIG.FLOATING_WIDTH,
    height: window.innerHeight - CHAT_CONFIG.FLOATING_TOP_OFFSET - CHAT_CONFIG.FLOATING_BOTTOM_OFFSET
  });
  
  const {
    state: { messages, isLoading, mode },
    shouldShowInlineChat,
    shouldShowFloatingChat,
    isMobile,
    closeChat,
    toggleMode,
    addMessage,
    setLoading, 
    setMessages
  } = useChatContext();

  const messagesRef = useRef(messages);
  useEffect(() => { messagesRef.current = messages; }, [messages]);

  // Auto-scroll to bottom when new messages arrive
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  // Don't render if this mode isn't active
  const shouldRender = isInlineMode ? shouldShowInlineChat : shouldShowFloatingChat;
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

    const currentMessages: {role: 'user' | 'assistant', content: string}[] = messages.map((m) => ({
      role: m.type,
      content: m.content
    }));
    currentMessages.push({role: 'user', content: content.trim()});
    
    // Add user message
    const userMessage: Message = {
      id: Date.now().toString(),
      type: 'user',
      content: content.trim(),
      timestamp: new Date()
    };
    addMessage(userMessage);
    setInputValue('');
    setLoading(true);

    let healthData: string = dataContext;

    if (!dataContext) {
      const healthData = 'User\'s health data:';
      setDataContext(healthData);
    }

    const assistantId = (Date.now() + 1).toString();
    addMessage({ id: assistantId, type: 'assistant', content: '', timestamp: new Date() });

    for await (const chunk of await server.llmChat.stream({
      messages: [{ role: "assistant", content: `User's health data: ${healthData}` }, ...currentMessages],
      model: "gpt-4o"
    })) {
      const list = messagesRef.current;
      const idx = list.findIndex(m => m.id === assistantId);
      if (idx >= 0) {
        const updated = [...list];
        updated[idx] = { ...updated[idx], content: updated[idx].content + chunk };
        setMessages(updated);
      }
    }
    setLoading(false);
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
            <div
              className={`max-w-[75%] px-3 py-2 rounded-lg text-sm ${
                message.type === 'user'
                  ? 'bg-[#1B1B1B] text-white rounded-br-sm'
                  : 'bg-[#1B1B1B] text-white rounded-bl-sm'
              }`}
            >
              <p className="whitespace-pre-wrap text-white">{message.content}</p>
              <p className="text-xs mt-1 text-gray-750">
                {message.timestamp.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}
              </p>
            </div>
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

  const renderInput = () => (
    <div className="bg-[#1C1C1C] p-4">
      <form onSubmit={handleSubmit} className="flex gap-2">
        <input
          type="text"
          value={inputValue}
          onChange={(e) => setInputValue(e.target.value)}
          onKeyPress={handleKeyPress}
          placeholder="Ask another question"
          className="flex-1 px-3 py-2 bg-[#1B1B1B] text-white placeholder-gray-750 rounded-full border border-gray-1500 focus:outline-none focus:ring-1 focus:ring-gray-1000 text-sm"
          disabled={isLoading}
        />
        <button
          type="submit"
          disabled={!inputValue.trim() || isLoading}
          className="px-3 py-2 text-white flex items-center justify-center"
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
        <button
          onClick={closeChat}
          className="text-gray-1000 transition-colors p-1 ml-auto"
        >
          <XClose className="w-5 h-5" invertIcon={false} />
        </button>
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
