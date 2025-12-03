import React from 'react';
import { motion } from 'framer-motion';
import { useChatContext } from '@/contexts/ChatContext';

export const FloatingChatBubble: React.FC = () => {
  const { 
    toggleChat, 
    state: { isOpen } 
  } = useChatContext();

  // Don't show bubble if chat panel is already open
  if (isOpen) return null;

  return (
    <motion.button
      whileHover={{ scale: 1.03 }}
      onClick={toggleChat}
      className="fixed bottom-6 right-6 z-50 w-12 h-12 bg-black text-white rounded-full shadow-xl flex items-center justify-center"
    >
      <img 
        src="/icons/CoolChat.svg" 
        alt="Chat" 
        className="w-6 h-6 filter brightness-0 invert"
      />
    </motion.button>
  );
}; 