import React, { useState, useRef, useEffect } from 'react';
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
  const dynamicWords = ['small molecule', 'protein', 'antibody', 'metabolite'];
  const [wordIndex, setWordIndex] = useState(0);
  const [charIndex, setCharIndex] = useState(0);
  const [isDeleting, setIsDeleting] = useState(false);
  
  const { 
    state: { voiceActive },
    openChat,
    addMessage,
    sendTextToVoice,
    startVoice
  } = useChatContext();

  const getBaseTextForWord = (word: string) => {
    if (!word) return "Let's design a ";
    const first = word[0]?.toLowerCase();
    const useAn = ['a', 'e', 'i', 'o', 'u'].includes(first);
    return `Let's design ${useAn ? 'an ' : 'a '}`;
  };

  useEffect(() => {
    // Don't animate while the user is actively typing something
    if (question.trim().length > 0) return;

    const currentWord = dynamicWords[wordIndex % dynamicWords.length];

    let timeout = isDeleting ? 70 : 110;

    // Pause with full word shown (cursor blinking) before deleting
    if (!isDeleting && charIndex === currentWord.length) {
      timeout = 2000; // ~2s pause
    }

    // Pause at the end of deletion (empty line + blinking caret)
    if (isDeleting && charIndex === 0) {
      timeout = 900;
    }

    const timer = setTimeout(() => {
      if (!isDeleting) {
        if (charIndex < currentWord.length) {
          setCharIndex((prev) => prev + 1);
        } else {
          setIsDeleting(true);
        }
      } else {
        if (charIndex > 0) {
          setCharIndex((prev) => prev - 1);
        } else {
          setIsDeleting(false);
          setWordIndex((prev) => (prev + 1) % dynamicWords.length);
        }
      }
    }, timeout);

    return () => clearTimeout(timer);
  }, [charIndex, isDeleting, wordIndex, question, dynamicWords]);

  const currentWord = dynamicWords[wordIndex % dynamicWords.length];
  const baseText = getBaseTextForWord(currentWord);
  const animatedPlaceholder = `${baseText}${currentWord.slice(0, charIndex)}`;

  const handleQuestionSubmit = async (question: string) => {
    openChat();
    
    // Auto-start voice if not active
    if (!voiceActive) {
      console.log('Voice not active, starting voice connection...');
      await startVoice();
      // Give it a moment to connect before sending
      await new Promise(resolve => setTimeout(resolve, 100));
    }
    
    // Always use ElevenLabs agent for all messages
    // sendTextToVoice will handle adding the user message to chat
    const sent = sendTextToVoice(question);
    
    if (!sent) {
      console.error('Failed to send to ElevenLabs agent');
      // Add error message
      const errorMessage: Message = {
        id: `error-${Date.now()}`,
        type: 'assistant',
        content: 'Sorry, I encountered an error connecting to the voice agent. Please try again.',
        timestamp: new Date(),
        isFinal: true
      };
      addMessage(errorMessage);
    }
    // ElevenLabs will handle the response through voice callbacks
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
          <div className="relative flex-1">
            <input
              ref={inputRef}
              type="text"
              value={question}
              onChange={(e) => setQuestion(e.target.value)}
              onKeyDown={handleKeyDown}
              onFocus={() => setIsFocused(true)}
              onBlur={() => setIsFocused(false)}
              placeholder=""
              className="w-full border-none outline-none text-gray-1500 placeholder-gray-1000 text-sm font-medium bg-transparent py-3 caret-gray-1500"
            />
            {question.trim().length === 0 && (
              <div className="pointer-events-none absolute inset-0 flex items-center text-gray-1000 text-sm font-medium">
                <span>{animatedPlaceholder}</span>
                <span className="typing-caret" />
              </div>
            )}
          </div>
          
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