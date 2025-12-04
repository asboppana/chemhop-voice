import React, { createContext, useContext, useReducer, useEffect, useRef, useCallback, type ReactNode } from 'react';
import { useIsMobile, useIsTablet } from '@/hooks/use-mobile';
import { voiceConversationService, type VoiceConnectionStatus, type VoiceAgentStatus } from '@/services/voiceConversation';
import { fetchVoiceAgentId } from '@/services/elevenLabs';

export interface Message {
  id: string;
  type: 'user' | 'assistant';
  content: string;
  timestamp: Date;
  isVoice?: boolean; // Flag to indicate if message came from voice
}

export type ChatAttachmentMode = 'attached' | 'detached';
export type ChatInputMode = 'text' | 'voice' | 'both';

interface ChatState {
  messages: Message[];
  isOpen: boolean;
  isLoading: boolean;
  mode: ChatAttachmentMode;
  // Voice-related state
  voiceActive: boolean;
  voiceConnectionStatus: VoiceConnectionStatus;
  voiceAgentStatus: VoiceAgentStatus;
  inputMode: ChatInputMode;
}

type ChatAction =
  | { type: 'OPEN_CHAT'; mode?: ChatAttachmentMode }
  | { type: 'CLOSE_CHAT' }
  | { type: 'TOGGLE_CHAT' }
  | { type: 'SET_MODE'; mode: ChatAttachmentMode }
  | { type: 'TOGGLE_MODE' }
  | { type: 'SET_LOADING'; loading: boolean }
  | { type: 'ADD_MESSAGE'; message: Message }
  | { type: 'SET_MESSAGES'; messages: Message[] }
  | { type: 'CLEAR_MESSAGES' }
  | { type: 'SET_VOICE_ACTIVE'; active: boolean }
  | { type: 'SET_VOICE_CONNECTION_STATUS'; status: VoiceConnectionStatus }
  | { type: 'SET_VOICE_AGENT_STATUS'; status: VoiceAgentStatus }
  | { type: 'SET_INPUT_MODE'; mode: ChatInputMode };

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
  isLoading: false,
  mode: getInitialMode(),
  voiceActive: false,
  voiceConnectionStatus: 'disconnected',
  voiceAgentStatus: 'idle',
  inputMode: 'both'
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
        isOpen: false,
        isLoading: false
      };
    
    case 'TOGGLE_CHAT':
      return {
        ...state,
        isOpen: !state.isOpen,
        isLoading: state.isOpen ? false : state.isLoading
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
    
    case 'SET_LOADING':
      return {
        ...state,
        isLoading: action.loading
      };
    
    case 'ADD_MESSAGE':
      return {
        ...state,
        messages: [...state.messages, action.message]
      };
    
    case 'SET_MESSAGES':
      return {
        ...state,
        messages: action.messages
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
    
    case 'SET_VOICE_CONNECTION_STATUS':
      return {
        ...state,
        voiceConnectionStatus: action.status
      };
    
    case 'SET_VOICE_AGENT_STATUS':
      return {
        ...state,
        voiceAgentStatus: action.status
      };
    
    case 'SET_INPUT_MODE':
      return {
        ...state,
        inputMode: action.mode
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
  setLoading: (loading: boolean) => void;
  addMessage: (message: Message) => void;
  setMessages: (messages: Message[]) => void;
  clearMessages: () => void;
  
  // Voice Actions
  startVoice: () => Promise<boolean>;
  stopVoice: () => Promise<void>;
  toggleVoice: () => Promise<void>;
  setInputMode: (mode: ChatInputMode) => void;
  sendTextToVoice: (text: string) => boolean;
  
  // Legacy compatibility (for easier migration)
  chatMessages: Message[];
  isChatPanelOpen: boolean;
  isChatLoading: boolean;
  attachmentMode: ChatAttachmentMode;
  setIsChatPanelOpen: (open: boolean) => void;
  setIsChatLoading: (loading: boolean) => void;
  toggleChatPanel: () => void;
  closeChatPanel: () => void;
  toggleAttachmentMode: () => void;
}

const ChatContext = createContext<ChatContextType | undefined>(undefined);

interface ChatProviderProps {
  children: ReactNode;
}

export const ChatProvider: React.FC<ChatProviderProps> = ({ children }) => {
  const [state, dispatch] = useReducer(chatReducer, initialState);
  const voiceInitializedRef = useRef(false);
  const currentUserMessageRef = useRef<{ id: string; content: string } | null>(null);
  const currentAssistantMessageRef = useRef<{ id: string; content: string } | null>(null);
  const messagesRef = useRef<Message[]>(state.messages);
  const chatAutoOpenedForVoiceRef = useRef(false);
  const voiceActiveRef = useRef(state.voiceActive);
  
  // Keep refs in sync with state
  useEffect(() => {
    messagesRef.current = state.messages;
  }, [state.messages]);
  
  useEffect(() => {
    voiceActiveRef.current = state.voiceActive;
  }, [state.voiceActive]);
  
  // Reset auto-open flag when voice becomes inactive
  useEffect(() => {
    if (!state.voiceActive) {
      chatAutoOpenedForVoiceRef.current = false;
    }
  }, [state.voiceActive]);
  
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
    if (isMobile) return; // Don't allow mode switching on mobile
    dispatch({ type: 'SET_MODE', mode });
  }, [isMobile]);
  
  const toggleMode = useCallback(() => {
    if (isMobile) return; // Don't allow mode switching on mobile
    dispatch({ type: 'TOGGLE_MODE' });
  }, [isMobile]);
  
  const setLoading = useCallback((loading: boolean) => {
    dispatch({ type: 'SET_LOADING', loading });
  }, []);
  
  const addMessage = useCallback((message: Message) => {
    dispatch({ type: 'ADD_MESSAGE', message });
  }, []);
  
  const setMessages = useCallback((messages: Message[]) => {
    dispatch({ type: 'SET_MESSAGES', messages });
  }, []);
  
  const clearMessages = useCallback(() => {
    dispatch({ type: 'CLEAR_MESSAGES' });
  }, []);
  
  // Voice Actions - Memoized to prevent infinite re-renders
  const startVoice = useCallback(async (): Promise<boolean> => {
    // If already active, don't start again
    if (voiceConversationService.isConversationActive()) {
      console.log('Voice conversation already active, skipping initialization');
      return true;
    }
    
    // Initialize agent ID if not done yet
    if (!voiceInitializedRef.current) {
      try {
        const agentId = await fetchVoiceAgentId();
        voiceConversationService.setAgentId(agentId);
        voiceInitializedRef.current = true;
      } catch (error) {
        console.error('Failed to fetch agent ID:', error);
        return false;
      }
    }
    
    // Setup callbacks
    voiceConversationService.setCallbacks({
      onStatusChange: (status) => {
        dispatch({ type: 'SET_VOICE_CONNECTION_STATUS', status });
      },
      onAgentStatusChange: (status) => {
        dispatch({ type: 'SET_VOICE_AGENT_STATUS', status });
      },
      onMessage: (voiceMessage) => {
        console.log('Voice message received:', voiceMessage);
        
        // Auto-open chat ONLY on user's first message, NOT on system greeting
        // This ensures chat opens when user speaks, not when agent greets
        if (!chatAutoOpenedForVoiceRef.current && 
            voiceMessage.message.trim() && 
            voiceMessage.source === 'user') {
          dispatch({ type: 'OPEN_CHAT' });
          chatAutoOpenedForVoiceRef.current = true;
        }
        
        // Handle user transcripts
        // NOTE: When text is sent via sendUserMessage, ElevenLabs does NOT echo it back
        // This callback only fires for actual voice input transcriptions
        if (voiceMessage.source === 'user') {
          console.log('User transcript received:', voiceMessage);
          
          if (!currentUserMessageRef.current || voiceMessage.isFinal) {
            // Create new user message or finalize existing one
            const messageId = currentUserMessageRef.current?.id || Date.now().toString();
            const message: Message = {
              id: messageId,
              type: 'user',
              content: voiceMessage.message,
              timestamp: new Date(),
              isVoice: true // This is truly from voice
            };
            
            if (currentUserMessageRef.current) {
              // Update existing message
              dispatch({ 
                type: 'SET_MESSAGES', 
                messages: messagesRef.current.map(m => 
                  m.id === messageId ? message : m
                )
              });
              if (voiceMessage.isFinal) {
                currentUserMessageRef.current = null;
              }
            } else {
              // Add new message
              dispatch({ type: 'ADD_MESSAGE', message });
              if (!voiceMessage.isFinal) {
                currentUserMessageRef.current = { id: messageId, content: voiceMessage.message };
              }
            }
          } else {
            // Update in-progress user message
            currentUserMessageRef.current.content = voiceMessage.message;
            dispatch({ 
              type: 'SET_MESSAGES', 
              messages: messagesRef.current.map(m => 
                m.id === currentUserMessageRef.current!.id 
                  ? { ...m, content: voiceMessage.message }
                  : m
              )
            });
          }
        }
        
        // Handle AI responses
        if (voiceMessage.source === 'ai') {
          console.log('AI response received:', {
            isFinal: voiceMessage.isFinal,
            messageLength: voiceMessage.message.length,
            preview: voiceMessage.message.substring(0, 50)
          });
          
          if (!currentAssistantMessageRef.current || voiceMessage.isFinal) {
            // Create new assistant message or finalize existing one
            const messageId = currentAssistantMessageRef.current?.id || (Date.now() + 1).toString();
            const message: Message = {
              id: messageId,
              type: 'assistant',
              content: voiceMessage.message,
              timestamp: new Date(),
              isVoice: true
            };
            
            if (currentAssistantMessageRef.current) {
              // Update existing message
              dispatch({ 
                type: 'SET_MESSAGES', 
                messages: messagesRef.current.map(m => 
                  m.id === messageId ? message : m
                )
              });
              if (voiceMessage.isFinal) {
                currentAssistantMessageRef.current = null;
              }
            } else {
              // Add new message
              dispatch({ type: 'ADD_MESSAGE', message });
              if (!voiceMessage.isFinal) {
                currentAssistantMessageRef.current = { id: messageId, content: voiceMessage.message };
              }
            }
          } else {
            // Update in-progress assistant message
            currentAssistantMessageRef.current.content = voiceMessage.message;
            dispatch({ 
              type: 'SET_MESSAGES', 
              messages: messagesRef.current.map(m => 
                m.id === currentAssistantMessageRef.current!.id 
                  ? { ...m, content: voiceMessage.message }
                  : m
              )
            });
          }
        }
      },
      onError: (error) => {
        console.error('Voice error:', error);
        dispatch({ type: 'SET_VOICE_ACTIVE', active: false });
        // Optionally add error message to chat
        const errorMessage: Message = {
          id: Date.now().toString(),
          type: 'assistant',
          content: `Voice error: ${error}`,
          timestamp: new Date()
        };
        dispatch({ type: 'ADD_MESSAGE', message: errorMessage });
      }
    });
    
    // Start the conversation
    const started = await voiceConversationService.startConversation();
    if (started) {
      dispatch({ type: 'SET_VOICE_ACTIVE', active: true });
      dispatch({ type: 'SET_INPUT_MODE', mode: 'voice' });
    }
    return started;
  }, []); // Empty deps - uses refs and dispatch which are stable
  
  const stopVoice = useCallback(async (): Promise<void> => {
    await voiceConversationService.stopConversation();
    dispatch({ type: 'SET_VOICE_ACTIVE', active: false });
    dispatch({ type: 'SET_INPUT_MODE', mode: 'text' });
    currentUserMessageRef.current = null;
    currentAssistantMessageRef.current = null;
  }, []);
  
  const toggleVoice = useCallback(async (): Promise<void> => {
    if (voiceActiveRef.current) {
      await stopVoice();
    } else {
      await startVoice();
    }
  }, [startVoice, stopVoice]);
  
  const setInputMode = useCallback((mode: ChatInputMode) => {
    dispatch({ type: 'SET_INPUT_MODE', mode });
  }, []);
  
  const sendTextToVoice = useCallback((text: string): boolean => {
    if (!voiceActiveRef.current) {
      console.warn('Cannot send text to voice: voice conversation not active');
      return false;
    }
    
    // Add user message immediately (marked as text, not voice)
    const userMessage: Message = {
      id: Date.now().toString(),
      type: 'user',
      content: text,
      timestamp: new Date(),
      isVoice: false // Typed message, not spoken
    };
    dispatch({ type: 'ADD_MESSAGE', message: userMessage });
    
    // Send to ElevenLabs - they will NOT echo this back as a user transcript
    // They will only send back the AI response
    const sent = voiceConversationService.sendTextMessage(text);
    
    if (!sent) {
      console.error('Failed to send text message to voice agent');
    }
    
    return sent;
  }, []);
  
  // DON'T cleanup voice on unmount - voice should persist
  // Voice will only stop when explicitly called via stopVoice() or toggleVoice()
  // This ensures the conversation stays active across navigation and component remounts
  
  const value: ChatContextType = {
    // State
    state,
    
    // Responsive states
    isMobile,
    isTablet,
    isDesktop,
    
    // Computed states
    shouldShowInlineChat,
    shouldShowFloatingChat,
    
    // Actions
    openChat,
    closeChat,
    toggleChat,
    setMode,
    toggleMode,
    setLoading,
    addMessage,
    setMessages,
    clearMessages,
    
    // Voice Actions
    startVoice,
    stopVoice,
    toggleVoice,
    setInputMode,
    sendTextToVoice,
    
    // Legacy compatibility
    chatMessages: state.messages,
    isChatPanelOpen: state.isOpen,
    isChatLoading: state.isLoading,
    attachmentMode: state.mode,
    setIsChatPanelOpen: (open: boolean) => {
      if (open) openChat();
      else closeChat();
    },
    setIsChatLoading: setLoading,
    toggleChatPanel: toggleChat,
    closeChatPanel: closeChat,
    toggleAttachmentMode: toggleMode
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