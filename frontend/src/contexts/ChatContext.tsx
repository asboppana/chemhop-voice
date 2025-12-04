import React, { createContext, useContext, useReducer, useEffect, useRef, useCallback, type ReactNode } from 'react';
import { useIsMobile, useIsTablet } from '@/hooks/use-mobile';
import { voiceASR, type VoiceASRCallbacks } from '@/services/elevenLabs';
import server from '@/app/server';

export interface Message {
  id: string;
  type: 'user' | 'assistant';
  content: string;
  timestamp: Date;
  isStreaming?: boolean; // Only for user voice transcripts
}

export type ChatAttachmentMode = 'attached' | 'detached';

interface ChatState {
  messages: Message[];
  isOpen: boolean;
  mode: ChatAttachmentMode;
  voiceActive: boolean;
  isLoading: boolean;
}

type ChatAction =
  | { type: 'OPEN_CHAT'; mode?: ChatAttachmentMode }
  | { type: 'CLOSE_CHAT' }
  | { type: 'TOGGLE_CHAT' }
  | { type: 'SET_MODE'; mode: ChatAttachmentMode }
  | { type: 'TOGGLE_MODE' }
  | { type: 'ADD_MESSAGE'; message: Message }
  | { type: 'UPDATE_MESSAGE'; id: string; content: string }
  | { type: 'CLEAR_MESSAGES' }
  | { type: 'SET_VOICE_ACTIVE'; active: boolean }
  | { type: 'SET_LOADING'; loading: boolean };

// Load initial mode from localStorage
const getInitialMode = (): ChatAttachmentMode => {
  if (typeof window === 'undefined') return 'attached';
  try {
    const saved = localStorage.getItem('sedona-chat-mode');
    return saved === 'detached' ? 'detached' : 'attached';
  } catch {
    return 'attached';
  }
};

const initialState: ChatState = {
  messages: [],
  isOpen: false,
  mode: getInitialMode(),
  voiceActive: false,
  isLoading: false
};

function chatReducer(state: ChatState, action: ChatAction): ChatState {
  switch (action.type) {
    case 'OPEN_CHAT':
      return {
        ...state,
        isOpen: true,
        mode: action.mode || state.mode
      };
    
    case 'CLOSE_CHAT':
      return {
        ...state,
        isOpen: false
      };
    
    case 'TOGGLE_CHAT':
      return {
        ...state,
        isOpen: !state.isOpen
      };
    
    case 'SET_MODE':
      return {
        ...state,
        mode: action.mode
      };
    
    case 'TOGGLE_MODE':
      return {
        ...state,
        mode: state.mode === 'attached' ? 'detached' : 'attached'
      };
    
    case 'ADD_MESSAGE':
      return {
        ...state,
        messages: [...state.messages, action.message]
      };
    
    case 'UPDATE_MESSAGE':
      return {
        ...state,
        messages: state.messages.map(msg =>
          msg.id === action.id ? { ...msg, content: action.content } : msg
        )
      };
    
    case 'CLEAR_MESSAGES':
      return {
        ...state,
        messages: []
      };
    
    case 'SET_VOICE_ACTIVE':
      return {
        ...state,
        voiceActive: action.active
      };
    
    case 'SET_LOADING':
      return {
        ...state,
        isLoading: action.loading
      };
    
    default:
      return state;
  }
}

interface ChatContextType {
  // State
  state: ChatState;
  
  // Responsive states
  isMobile: boolean;
  isTablet: boolean;
  isDesktop: boolean;
  
  // Computed states
  shouldShowInlineChat: boolean;
  shouldShowFloatingChat: boolean;
  
  // Actions
  openChat: (mode?: ChatAttachmentMode) => void;
  closeChat: () => void;
  toggleChat: () => void;
  setMode: (mode: ChatAttachmentMode) => void;
  toggleMode: () => void;
  clearMessages: () => void;
  sendMessage: (content: string) => Promise<void>;
  
  // Voice Actions
  startVoice: () => boolean;
  stopVoice: () => void;
  toggleVoice: () => void;
}

const ChatContext = createContext<ChatContextType | undefined>(undefined);

interface ChatProviderProps {
  children: ReactNode;
}

// Message ID generator
let messageIdCounter = 0;
const generateMessageId = (): string => {
  return `msg-${Date.now()}-${messageIdCounter++}`;
};

export const ChatProvider: React.FC<ChatProviderProps> = ({ children }) => {
  const [state, dispatch] = useReducer(chatReducer, initialState);
  const currentUserMessageRef = useRef<string | null>(null);
  const chatAutoOpenedRef = useRef(false);
  
  // Responsive detection
  const isMobile = useIsMobile();
  const isTablet = useIsTablet();
  const isDesktop = !isMobile && !isTablet;
  
  // Persist mode to localStorage
  useEffect(() => {
    try {
      localStorage.setItem('sedona-chat-mode', state.mode);
    } catch (error) {
      console.warn('Failed to persist chat mode:', error);
    }
  }, [state.mode]);
  
  // Computed states
  const shouldShowInlineChat = state.isOpen && state.mode === 'attached' && !isMobile;
  const shouldShowFloatingChat = state.isOpen && (state.mode === 'detached' || isMobile);
  
  // Actions - Memoized to prevent unnecessary re-renders
  const openChat = useCallback((mode?: ChatAttachmentMode) => {
    dispatch({ type: 'OPEN_CHAT', mode });
  }, []);
  
  const closeChat = useCallback(() => {
    dispatch({ type: 'CLOSE_CHAT' });
  }, []);
  
  const toggleChat = useCallback(() => {
    dispatch({ type: 'TOGGLE_CHAT' });
  }, []);
  
  const setMode = useCallback((mode: ChatAttachmentMode) => {
    if (isMobile) return;
    dispatch({ type: 'SET_MODE', mode });
  }, [isMobile]);
  
  const toggleMode = useCallback(() => {
    if (isMobile) return;
    dispatch({ type: 'TOGGLE_MODE' });
  }, [isMobile]);
  
  const clearMessages = useCallback(() => {
    dispatch({ type: 'CLEAR_MESSAGES' });
  }, []);
  
  // Call backend API with user message
  const callBackendAPI = useCallback(async (userMessage: string) => {
    dispatch({ type: 'SET_LOADING', loading: true });
    
    try {
      // Build conversation history for context - only last 3 messages
      const conversationHistory = state.messages
        .filter(m => !m.isStreaming) // Only include finalized messages
        .slice(-3) // Keep only last 3 messages for context
        .map(m => ({
          role: m.type === 'user' ? 'user' as const : 'assistant' as const,
          content: m.content
        }));
      
      // Add current user message
      conversationHistory.push({
        role: 'user',
        content: userMessage
      });
      
      // Call backend
      const response = await server.voiceAPI.callChatAgent({
        messages: conversationHistory
      });
      
      // Add assistant response to chat
      dispatch({
        type: 'ADD_MESSAGE',
        message: {
          id: generateMessageId(),
          type: 'assistant',
          content: response.response || response.text || 'No response',
          timestamp: new Date(),
          isStreaming: false
        }
      });
      
      console.log('📦 Backend response:', response);
      
      // TODO: Handle structured JSON for main page display
      // You can emit an event or use a callback to pass structured data to MainPage
      if (response.structured) {
        console.log('📊 Structured data:', response.structured);
        // window.dispatchEvent(new CustomEvent('structured-data', { detail: response.structured }));
      }
      
    } catch (error) {
      console.error('Failed to call backend:', error);
      dispatch({
        type: 'ADD_MESSAGE',
        message: {
          id: generateMessageId(),
          type: 'assistant',
          content: 'Sorry, I encountered an error. Please try again.',
          timestamp: new Date(),
          isStreaming: false
        }
      });
    } finally {
      dispatch({ type: 'SET_LOADING', loading: false });
    }
  }, [state.messages]);
  
  // Send a text message (from typed input)
  const sendMessage = useCallback(async (content: string) => {
    if (!content.trim()) return;
    
    // Add user message to chat immediately
    dispatch({
      type: 'ADD_MESSAGE',
      message: {
        id: generateMessageId(),
        type: 'user',
        content: content.trim(),
        timestamp: new Date(),
        isStreaming: false
      }
    });
    
    // Call backend API
    await callBackendAPI(content.trim());
  }, [callBackendAPI]);
  
  // Voice Actions
  const startVoice = useCallback((): boolean => {
    if (voiceASR.isActive()) {
      console.log('Voice already active');
      return true;
    }
    
    // Setup callbacks
    const callbacks: VoiceASRCallbacks = {
      onTranscript: (text, isFinal) => {
        console.log('🎤 Transcript:', { text, isFinal });
        
        // Auto-open chat on first user speech
        if (!chatAutoOpenedRef.current && text.trim()) {
          dispatch({ type: 'OPEN_CHAT' });
          chatAutoOpenedRef.current = true;
        }
        
        if (isFinal) {
          // Finalize the message and call backend
          if (currentUserMessageRef.current) {
            // Update message to be final
            dispatch({
              type: 'UPDATE_MESSAGE',
              id: currentUserMessageRef.current,
              content: text
            });
            currentUserMessageRef.current = null;
            
            // Call backend API
            callBackendAPI(text);
          }
        } else {
          // Streaming transcript - show in real-time
          if (!currentUserMessageRef.current) {
            const id = generateMessageId();
            currentUserMessageRef.current = id;
            dispatch({
              type: 'ADD_MESSAGE',
              message: {
                id,
                type: 'user',
                content: text,
                timestamp: new Date(),
                isStreaming: true
              }
            });
          } else {
            dispatch({
              type: 'UPDATE_MESSAGE',
              id: currentUserMessageRef.current,
              content: text
            });
          }
        }
      },
      onError: (error) => {
        console.error('ASR error:', error);
        // Don't stop voice on error, just log it
      }
    };
    
    voiceASR.setCallbacks(callbacks);
    
    // Start listening
    const started = voiceASR.startListening();
    if (started) {
      dispatch({ type: 'SET_VOICE_ACTIVE', active: true });
    }
    return started;
  }, [callBackendAPI]);
  
  const stopVoice = useCallback((): void => {
    voiceASR.stopListening();
    dispatch({ type: 'SET_VOICE_ACTIVE', active: false });
    currentUserMessageRef.current = null;
    chatAutoOpenedRef.current = false;
  }, []);
  
  const toggleVoice = useCallback((): void => {
    if (state.voiceActive) {
      stopVoice();
    } else {
      startVoice();
    }
  }, [state.voiceActive, startVoice, stopVoice]);
  
  const value: ChatContextType = {
    state,
    isMobile,
    isTablet,
    isDesktop,
    shouldShowInlineChat,
    shouldShowFloatingChat,
    openChat,
    closeChat,
    toggleChat,
    setMode,
    toggleMode,
    clearMessages,
    sendMessage,
    startVoice,
    stopVoice,
    toggleVoice
  };
  
  return (
    <ChatContext.Provider value={value}>
      {children}
    </ChatContext.Provider>
  );
};

export const useChatContext = (): ChatContextType => {
  const context = useContext(ChatContext);
  if (!context) {
    throw new Error('useChatContext must be used within a ChatProvider');
  }
  return context;
}; 