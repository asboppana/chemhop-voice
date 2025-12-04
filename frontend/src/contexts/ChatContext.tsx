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
  isFinal: boolean; // Whether this message is complete
}

export type ChatAttachmentMode = 'attached' | 'detached';
export type ChatInputMode = 'text' | 'voice' | 'both';

interface ChatState {
  messages: Message[]; // Only FINAL messages
  streamingMessage: Message | null; // Current message being streamed
  streamingChunks: string[]; // Progressive chunks for display sync with voice
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
  | { type: 'START_STREAMING_MESSAGE'; message: Message }
  | { type: 'ADD_STREAMING_CHUNK'; chunk: string }
  | { type: 'UPDATE_STREAMING_CONTENT'; content: string }
  | { type: 'FINALIZE_STREAMING_MESSAGE'; finalContent?: string }
  | { type: 'INTERRUPT_STREAMING_MESSAGE' }
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
  streamingMessage: null,
  streamingChunks: [],
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
        messages: [],
        streamingMessage: null,
        streamingChunks: []
      };
    
    case 'START_STREAMING_MESSAGE':
      return {
        ...state,
        streamingMessage: action.message,
        streamingChunks: []
      };
    
    case 'ADD_STREAMING_CHUNK':
      return {
        ...state,
        streamingChunks: [...state.streamingChunks, action.chunk]
      };
    
    case 'UPDATE_STREAMING_CONTENT':
      if (!state.streamingMessage) return state;
      return {
        ...state,
        streamingMessage: {
          ...state.streamingMessage,
          content: action.content
        }
      };
    
    case 'FINALIZE_STREAMING_MESSAGE':
      if (!state.streamingMessage) return state;
      const finalMessage: Message = {
        ...state.streamingMessage,
        content: action.finalContent || state.streamingMessage.content,
        isFinal: true
      };
      return {
        ...state,
        messages: [...state.messages, finalMessage],
        streamingMessage: null,
        streamingChunks: []
      };
    
    case 'INTERRUPT_STREAMING_MESSAGE':
      if (!state.streamingMessage) return state;
      // Save interrupted message with whatever content was visible
      const interruptedMessage: Message = {
        ...state.streamingMessage,
        content: state.streamingChunks.join(' ') || state.streamingMessage.content,
        isFinal: true
      };
      return {
        ...state,
        messages: [...state.messages, interruptedMessage],
        streamingMessage: null,
        streamingChunks: []
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

// Message ID generator with monotonic counter
let messageIdCounter = 0;
const generateMessageId = (): string => {
  return `msg-${Date.now()}-${messageIdCounter++}`;
};

export const ChatProvider: React.FC<ChatProviderProps> = ({ children }) => {
  const [state, dispatch] = useReducer(chatReducer, initialState);
  const voiceInitializedRef = useRef(false);
  const currentUserMessageIdRef = useRef<string | null>(null);
  const currentAssistantMessageIdRef = useRef<string | null>(null);
  const currentAssistantFullTextRef = useRef<string>(''); // Store full AI text for interruptions
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
      onChunk: (chunk, source) => {
        console.log('📝 Chunk received:', { chunk: chunk.substring(0, 30), source });
        // Add chunk to streaming display
        dispatch({ type: 'ADD_STREAMING_CHUNK', chunk });
      },
      onMessage: (voiceMessage) => {
        console.log('📨 Voice message received:', voiceMessage);
        
        // Auto-open chat ONLY on user's first message, NOT on system greeting
        if (!chatAutoOpenedForVoiceRef.current && 
            voiceMessage.message.trim() && 
            voiceMessage.source === 'user') {
          dispatch({ type: 'OPEN_CHAT' });
          chatAutoOpenedForVoiceRef.current = true;
        }
        
        // ============================================================
        // HANDLE USER TRANSCRIPTS (Voice Input)
        // ============================================================
        if (voiceMessage.source === 'user') {
          console.log('🎤 User transcript:', { isFinal: voiceMessage.isFinal, content: voiceMessage.message });
          
          // User interrupted AI - finalize AI's message first
          if (currentAssistantMessageIdRef.current) {
            console.log('🚨 User interrupting AI, finalizing AI message');
            dispatch({ type: 'INTERRUPT_STREAMING_MESSAGE' });
            currentAssistantMessageIdRef.current = null;
            currentAssistantFullTextRef.current = '';
          }
          
          if (voiceMessage.isFinal) {
            // User finished speaking - finalize message
            console.log('✅ User transcript FINAL');
            dispatch({ 
              type: 'FINALIZE_STREAMING_MESSAGE',
              finalContent: voiceMessage.message 
            });
            currentUserMessageIdRef.current = null;
          } else {
            // User still speaking - stream the transcript
            if (!currentUserMessageIdRef.current) {
              // Start new user message
              const messageId = generateMessageId();
              currentUserMessageIdRef.current = messageId;
              dispatch({ 
                type: 'START_STREAMING_MESSAGE',
                message: {
                  id: messageId,
                  type: 'user',
                  content: voiceMessage.message,
                  timestamp: new Date(),
                  isVoice: true,
                  isFinal: false
                }
              });
            } else {
              // Update streaming user message content
              dispatch({ 
                type: 'UPDATE_STREAMING_CONTENT',
                content: voiceMessage.message 
              });
            }
          }
        }
        
        // ============================================================
        // HANDLE AI RESPONSES (Full text arrives immediately)
        // ============================================================
        if (voiceMessage.source === 'ai') {
          console.log('🤖 AI response:', {
            isFinal: voiceMessage.isFinal,
            messageLength: voiceMessage.message.length,
            preview: voiceMessage.message.substring(0, 50)
          });
          
          if (voiceMessage.isFinal) {
            // Final AI message - finalize streaming
            console.log('✅ AI response FINAL');
            dispatch({ 
              type: 'FINALIZE_STREAMING_MESSAGE',
              finalContent: voiceMessage.message 
            });
            currentAssistantMessageIdRef.current = null;
            currentAssistantFullTextRef.current = '';
          } else {
            // New AI response starting - store FULL text for reference
            if (!currentAssistantMessageIdRef.current) {
              const messageId = generateMessageId();
              currentAssistantMessageIdRef.current = messageId;
              currentAssistantFullTextRef.current = voiceMessage.message;
              
              dispatch({ 
                type: 'START_STREAMING_MESSAGE',
                message: {
                  id: messageId,
                  type: 'assistant',
                  content: voiceMessage.message, // Full text stored
                  timestamp: new Date(),
                  isVoice: true,
                  isFinal: false
                }
              });
            } else {
              // Update full text (in case it changes during streaming)
              currentAssistantFullTextRef.current = voiceMessage.message;
              dispatch({ 
                type: 'UPDATE_STREAMING_CONTENT',
                content: voiceMessage.message 
              });
            }
          }
        }
      },
      onError: (error) => {
        console.error('Voice error:', error);
        dispatch({ type: 'SET_VOICE_ACTIVE', active: false });
        // Optionally add error message to chat
        const errorMessage: Message = {
          id: generateMessageId(),
          type: 'assistant',
          content: `Voice error: ${error}`,
          timestamp: new Date(),
          isFinal: true
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
    currentUserMessageIdRef.current = null;
    currentAssistantMessageIdRef.current = null;
    currentAssistantFullTextRef.current = '';
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
    
    // Interrupt any streaming AI message if user types while AI is responding
    if (currentAssistantMessageIdRef.current) {
      console.log('🚨 User typing interrupted AI, finalizing AI message');
      dispatch({ type: 'INTERRUPT_STREAMING_MESSAGE' });
      currentAssistantMessageIdRef.current = null;
      currentAssistantFullTextRef.current = '';
    }
    
    // Add user message immediately (marked as text, not voice)
    // Typed messages are always final (not streamed)
    const userMessage: Message = {
      id: generateMessageId(),
      type: 'user',
      content: text,
      timestamp: new Date(),
      isVoice: false, // Typed message, not spoken
      isFinal: true // Typed messages are always complete
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