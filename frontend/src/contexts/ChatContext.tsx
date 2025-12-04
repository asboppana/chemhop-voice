import React, { createContext, useContext, useReducer, useEffect, useRef, useCallback, type ReactNode } from 'react';
import { useIsMobile, useIsTablet } from '@/hooks/use-mobile';
import { voiceASR, type VoiceASRCallbacks } from '@/services/elevenLabs';
import server from '@/app/server';
import { parseAgentResponse, debugLogParsedResponse } from '@/utils/chatResponseHandler';
import { dispatchStructuredData, getMoleculeContextString } from '@/hooks/useStructuredData';

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
      const conversationHistory: Array<{ role: 'user' | 'assistant', content: string }> = [];
      
      // Add recent message history
      state.messages
        .filter(m => !m.isStreaming) // Only include finalized messages
        .slice(-3) // Keep only last 3 messages for context
        .forEach(m => {
          conversationHistory.push({
            role: m.type === 'user' ? 'user' as const : 'assistant' as const,
            content: m.content
          });
        });
      
      // =========================================================================
      // APPEND MOLECULE CONTEXT TO USER MESSAGE IF AVAILABLE
      // =========================================================================
      // This allows the backend to understand references like "structure A",
      // "permute pattern B", etc. by providing the current molecule state
      // We append it to the user message so the LLM has context about the molecule
      const moleculeContext = getMoleculeContextString();
      let finalUserMessage = userMessage;
      if (moleculeContext) {
        finalUserMessage = `${userMessage}\n\n${moleculeContext}`;
        console.log('🧬 Including molecule context in user message');
      }
      
      // Add current user message (with optional molecule context appended)
      conversationHistory.push({
        role: 'user',
        content: finalUserMessage
      });
      
      // Call backend
      const response = await server.voiceAPI.callChatAgent({
        messages: conversationHistory
      });
      
      // =========================================================================
      // PARSE RESPONSE USING chatResponseHandler
      // =========================================================================
      // This maps backend response types to frontend render actions:
      //   - text_only       → RENDER_TEXT
      //   - annotation      → RENDER_MOLECULE
      //   - query_set       → RENDER_QUERY_SET
      //   - bioisosteres    → RENDER_BIOISOSTERES
      //   - admet           → RENDER_ADMET
      //   - multi_bioisostere → RENDER_MULTI_BIO
      //   - multi_admet     → RENDER_MULTI_ADMET
      // =========================================================================
      const parsed = parseAgentResponse(response);
      
      // Debug log in development
      if (import.meta.env.DEV) {
        debugLogParsedResponse(parsed);
      }
      
      // Add assistant response text to chat
      dispatch({
        type: 'ADD_MESSAGE',
        message: {
          id: generateMessageId(),
          type: 'assistant',
          content: parsed.textContent || 'No response',
          timestamp: new Date(),
          isStreaming: false
        }
      });
      
      // Dispatch structured data to MainPage if present
      if (parsed.hasStructuredData) {
        dispatchStructuredData(parsed);
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
    
    // Count words to filter out short messages
    const wordCount = content.trim().split(/\s+/).filter(w => w.length > 0).length;
    
    if (wordCount < 5) {
      console.log(`🚫 Message too short (${wordCount} words), minimum 5 words required`);
      return;
    }
    
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
            // Count words to filter out short/noise transcripts
            const wordCount = text.trim().split(/\s+/).filter(w => w.length > 0).length;
            
            if (wordCount < 5) {
              console.log(`🚫 Ignoring short transcript (${wordCount} words): "${text}"`);
              // Remove the streaming message
              dispatch({
                type: 'UPDATE_MESSAGE',
                id: currentUserMessageRef.current,
                content: text
              });
              currentUserMessageRef.current = null;
              return;
            }
            
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