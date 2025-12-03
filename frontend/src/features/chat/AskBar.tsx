import React, { useState, useRef } from 'react';
import HoverRotatingSedonaLogo from '@/components/animations/HoverRotatingSedonaLogo';
import { useChatContext } from '@/contexts/ChatContext';
import type { Message } from '@/contexts/ChatContext';

interface AskBarProps {
  className?: string;
}

export const AskBar: React.FC<AskBarProps> = ({ className = "" }) => {
  const [question, setQuestion] = useState('');
  const [isFocused, setIsFocused] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);
  
  const { 
    openChat, 
    setLoading, 
    addMessage
  } = useChatContext();

  const handleQuestionSubmit = async (question: string) => {
    // Add user message to chat
    const userMessage: Message = {
      id: Date.now().toString(),
      type: 'user',
      content: question,
      timestamp: new Date()
    };
    
    addMessage(userMessage);
    openChat(); // Use new action instead of setIsChatPanelOpen(true)
    setLoading(true);
    
    try {
      // TODO: Replace with actual API call to your backend
      console.log('Submitting question:', question);
      
      // Simulate API call
      await new Promise(resolve => setTimeout(resolve, 1500));
      
      // Simulate assistant response
      const assistantMessage: Message = {
        id: (Date.now() + 1).toString(),
        type: 'assistant',
        content: `Feature coming soon!`,
        timestamp: new Date()
      };
      
      addMessage(assistantMessage);
    } catch (error) {
      console.error('Error submitting question:', error);
      
      // Add error message
      const errorMessage: Message = {
        id: (Date.now() + 1).toString(),
        type: 'assistant',
        content: 'Sorry, I encountered an error while processing your question. Please try again.',
        timestamp: new Date()
      };
      
      addMessage(errorMessage);
    } finally {
      setLoading(false);
    }
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (!question.trim()) return;

    handleQuestionSubmit(question.trim());
    setQuestion('');
  };

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSubmit(e);
    }
  };

  return (
    <div className={`fixed bottom-6 left-1/2 transform -translate-x-1/2 z-40 ${className}`}>
      <div className={`bg-white rounded-full drop-shadow-lg pl-4 pr-0 py-0 max-w-[calc(100vw-2rem)] transition-all duration-300 ease-out ${
        isFocused ? 'w-[400px]' : 'w-[320px]'
      }`}>
        <form onSubmit={handleSubmit} className="flex items-center gap-4">
          <HoverRotatingSedonaLogo size={24} color="black" className="flex-shrink-0" />
          <input
            ref={inputRef}
            type="text"
            value={question}
            onChange={(e) => setQuestion(e.target.value)}
            onKeyDown={handleKeyDown}
            onFocus={() => setIsFocused(true)}
            onBlur={() => setIsFocused(false)}
            placeholder="Start designing a molecule"
            className="flex-1 border-none outline-none text-gray-1500 placeholder-gray-1000 text-sm font-medium bg-transparent py-3 caret-gray-1500"
          />
          
          <button
            type="submit"
            disabled={!question.trim()}
            className="flex-shrink-0 transition-all duration-200 pr-2"
          >
            {question.trim() ?      
            <img
              src="/icons/Arrow.svg"
              alt="Send"
              className={`w-8 h-8 rotate-180 filter invert`}
            /> : <img
              src="/icons/GrayArrow.svg"
              alt="Ask"
              className={`w-8 h-8`}
            />}
          </button>
        </form>
      </div>
    </div>
  );
};