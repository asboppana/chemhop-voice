import React, { createContext, useContext, useReducer, useEffect, type ReactNode } from 'react';
import { useIsMobile, useIsTablet } from '@/hooks/use-mobile';

export interface Message {
  id: string;
  type: 'user' | 'assistant';
  content: string;
  timestamp: Date;
}

export type ChatAttachmentMode = 'attached' | 'detached';

interface ChatState {
  messages: Message[];
  isOpen: boolean;
  isLoading: boolean;
  mode: ChatAttachmentMode;
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
  | { type: 'CLEAR_MESSAGES' };

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
  mode: getInitialMode()
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
  
  // Actions
  const openChat = (mode?: ChatAttachmentMode) => {
    dispatch({ type: 'OPEN_CHAT', mode });
  };
  
  const closeChat = () => {
    dispatch({ type: 'CLOSE_CHAT' });
  };
  
  const toggleChat = () => {
    dispatch({ type: 'TOGGLE_CHAT' });
  };
  
  const setMode = (mode: ChatAttachmentMode) => {
    if (isMobile) return; // Don't allow mode switching on mobile
    dispatch({ type: 'SET_MODE', mode });
  };
  
  const toggleMode = () => {
    if (isMobile) return; // Don't allow mode switching on mobile
    dispatch({ type: 'TOGGLE_MODE' });
  };
  
  const setLoading = (loading: boolean) => {
    dispatch({ type: 'SET_LOADING', loading });
  };
  
  const addMessage = (message: Message) => {
    dispatch({ type: 'ADD_MESSAGE', message });
  };
  
  const setMessages = (messages: Message[]) => {
    dispatch({ type: 'SET_MESSAGES', messages });
  };
  
  const clearMessages = () => {
    dispatch({ type: 'CLEAR_MESSAGES' });
  };
  
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