import React, { useState, useEffect, useRef } from 'react';
import HoverRotatingSedonaLogo from '@/components/animations/HoverRotatingSedonaLogo';
import { useChatContext } from '@/contexts/ChatContext';

interface AskBarProps {
  className?: string;
}

export const AskBar: React.FC<AskBarProps> = ({ className = "" }) => {
  const dynamicWords = ['small molecule', 'protein', 'antibody', 'metabolite'];
  const [wordIndex, setWordIndex] = useState(0);
  const [charIndex, setCharIndex] = useState(0);
  const [isDeleting, setIsDeleting] = useState(false);
  const [inputValue, setInputValue] = useState('');
  const [isFocused, setIsFocused] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);
  
  const { 
    state: { voiceActive },
    openChat,
    startVoice,
    sendMessage
  } = useChatContext();

  const getBaseTextForWord = (word: string) => {
    if (!word) return "Let's design a ";
    const first = word[0]?.toLowerCase();
    const useAn = ['a', 'e', 'i', 'o', 'u'].includes(first);
    return `Let's design ${useAn ? 'an ' : 'a '}`;
  };

  // Only run animation when not focused and no input
  useEffect(() => {
    if (isFocused || inputValue) return;

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
  }, [charIndex, isDeleting, wordIndex, dynamicWords, isFocused, inputValue]);

  const currentWord = dynamicWords[wordIndex % dynamicWords.length];
  const baseText = getBaseTextForWord(currentWord);
  const animatedPlaceholder = `${baseText}${currentWord.slice(0, charIndex)}`;

  const handleVoiceClick = (e: React.MouseEvent) => {
    e.stopPropagation();
    
    // Open chat first
    openChat();
    
    // Then start voice listening
    const started = startVoice();
    if (!started) {
      console.error('Failed to start voice recognition. Check microphone permissions.');
      // You could show a toast notification here
    } else {
      console.log('✅ Voice recognition started from AskBar');
    }
  };

  const handleContainerClick = () => {
    // Focus the input when clicking anywhere on the bar
    inputRef.current?.focus();
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (!inputValue.trim()) return;

    // Open chat and send message
    openChat();
    sendMessage(inputValue);
    setInputValue('');
    setIsFocused(false);
    inputRef.current?.blur();
  };

  return (
    <div className={`fixed bottom-6 left-1/2 transform -translate-x-1/2 z-40 ${className}`}>
      <form
        onSubmit={handleSubmit}
        onClick={handleContainerClick}
        className="bg-white rounded-full drop-shadow-lg px-6 py-3 hover:shadow-xl transition-all duration-300 ease-out flex items-center gap-4 cursor-text"
        style={{
          width: isFocused || inputValue ? '450px' : '320px',
          transition: 'width 0.3s ease-out'
        }}
      >
        <HoverRotatingSedonaLogo size={24} color="black" className="flex-shrink-0" />
        
        <div className="flex-1 relative text-gray-1000 text-sm font-medium">
          {/* Animated placeholder - only show when not focused and no input */}
          {!isFocused && !inputValue && (
            <div className="absolute inset-0 flex items-center gap-2 pointer-events-none">
              <span>{animatedPlaceholder}</span>
              <span className="typing-caret" />
            </div>
          )}
          
          {/* Actual input field */}
          <input
            ref={inputRef}
            type="text"
            value={inputValue}
            onChange={(e) => setInputValue(e.target.value)}
            onFocus={() => setIsFocused(true)}
            onBlur={() => setIsFocused(false)}
            className="w-full bg-transparent outline-none text-gray-1000 placeholder-gray-500"
            placeholder={isFocused ? "Type your message..." : ""}
          />
        </div>

        <button
          type="button"
          onClick={handleVoiceClick}
          className={`flex-shrink-0 relative transition-all ${
            voiceActive 
              ? 'text-green-500 scale-110' 
              : 'text-gray-1500 hover:text-gray-1000'
          }`}
          title={voiceActive ? 'Voice active - listening...' : 'Start voice conversation'}
        >
          <svg className="w-5 h-5" fill="currentColor" viewBox="0 0 20 20">
            <path fillRule="evenodd" d="M7 4a3 3 0 016 0v4a3 3 0 11-6 0V4zm4 10.93A7.001 7.001 0 0017 8a1 1 0 10-2 0A5 5 0 015 8a1 1 0 00-2 0 7.001 7.001 0 006 6.93V17H6a1 1 0 100 2h8a1 1 0 100-2h-3v-2.07z" clipRule="evenodd" />
          </svg>
          {voiceActive && (
            <div className="absolute -top-1 -right-1 w-2 h-2 bg-green-500 rounded-full animate-pulse" />
          )}
        </button>
      </form>
    </div>
  );
};