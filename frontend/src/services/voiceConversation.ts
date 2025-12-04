/**
 * ElevenLabs Voice Conversation Service
 * Centralized management for voice conversations with transcript handling
 */

import { Conversation } from '@elevenlabs/client';

export type VoiceConnectionStatus = 'disconnected' | 'connecting' | 'connected' | 'error';
export type VoiceAgentStatus = 'idle' | 'listening' | 'speaking';

export interface VoiceMessage {
  message: string;
  source: 'user' | 'ai';
  isFinal?: boolean;
}

export interface VoiceConversationCallbacks {
  onStatusChange?: (status: VoiceConnectionStatus) => void;
  onAgentStatusChange?: (status: VoiceAgentStatus) => void;
  onMessage?: (message: VoiceMessage) => void;
  onChunk?: (chunk: string, source: 'user' | 'ai') => void; // Progressive text as voice speaks
  onError?: (error: string) => void;
  onVadScore?: (score: number) => void;
}

class VoiceConversationService {
  private conversation: any = null;
  private agentId: string | null = null;
  private callbacks: VoiceConversationCallbacks = {};
  private isActive = false;
  private isConnecting = false; // Prevent race conditions during connection
  private shouldStayConnected = false; // Track if user wants to stay connected
  private reconnectAttempts = 0;
  private maxReconnectAttempts = 3;

  /**
   * Initialize the conversation service with an agent ID
   */
  setAgentId(agentId: string) {
    this.agentId = agentId;
  }

  /**
   * Register callbacks for conversation events
   */
  setCallbacks(callbacks: VoiceConversationCallbacks) {
    this.callbacks = { ...this.callbacks, ...callbacks };
  }

  /**
   * Check if a conversation is currently active
   */
  isConversationActive(): boolean {
    return this.isActive && this.conversation !== null;
  }

  /**
   * Request microphone permission
   */
  async requestMicrophonePermission(): Promise<boolean> {
    try {
      await navigator.mediaDevices.getUserMedia({ audio: true });
      return true;
    } catch (error) {
      console.error('Microphone permission denied:', error);
      this.callbacks.onError?.('Microphone permission denied. Please allow microphone access to use voice features.');
      return false;
    }
  }

  /**
   * Start a new voice conversation
   */
  async startConversation(): Promise<boolean> {
    if (!this.agentId) {
      this.callbacks.onError?.('No agent ID configured');
      return false;
    }

    if (this.isActive) {
      console.warn('Conversation already active');
      return true;
    }

    if (this.isConnecting) {
      console.warn('Connection already in progress, waiting...');
      return true;
    }

    try {
      this.isConnecting = true;
      
      // Request microphone permission first
      const hasPermission = await this.requestMicrophonePermission();
      if (!hasPermission) {
        this.isConnecting = false;
        return false;
      }

      this.callbacks.onStatusChange?.('connecting');
      this.shouldStayConnected = true; // Mark that user wants to stay connected

      // Start the ElevenLabs conversation
      this.conversation = await (Conversation as any).startSession({
        agentId: this.agentId,
        
        onConnect: () => {
          console.log('✅ Voice conversation connected', {
            timestamp: new Date().toISOString(),
            agentId: this.agentId
          });
          this.isActive = true;
          this.isConnecting = false;
          this.reconnectAttempts = 0; // Reset reconnect counter on successful connection
          this.callbacks.onStatusChange?.('connected');
          this.callbacks.onAgentStatusChange?.('listening');
        },
        
        onDisconnect: (event: any) => {
          console.warn('⚠️ Voice conversation disconnected', { 
            event, 
            shouldStayConnected: this.shouldStayConnected,
            reconnectAttempts: this.reconnectAttempts,
            timestamp: new Date().toISOString()
          });
          console.trace('Disconnect trace');
          this.isActive = false;
          this.isConnecting = false;
          
          // If user wants to stay connected and we haven't exceeded retry limit, attempt reconnect
          if (this.shouldStayConnected && this.reconnectAttempts < this.maxReconnectAttempts) {
            this.reconnectAttempts++;
            console.log(`🔄 Attempting auto-reconnect (${this.reconnectAttempts}/${this.maxReconnectAttempts})...`);
            this.callbacks.onStatusChange?.('connecting');
            setTimeout(() => {
              this.startConversation().catch(err => {
                console.error('Auto-reconnect failed:', err);
              });
            }, 1000); // Wait 1 second before reconnecting
          } else {
            this.callbacks.onStatusChange?.('disconnected');
            this.callbacks.onAgentStatusChange?.('idle');
            if (this.reconnectAttempts >= this.maxReconnectAttempts) {
              console.error('❌ Max reconnection attempts reached');
              this.callbacks.onError?.('Connection lost. Please try again.');
            }
          }
        },
        
        onError: (error: any) => {
          console.error('❌ Voice conversation error:', error);
          this.isActive = false;
          this.isConnecting = false;
          this.callbacks.onStatusChange?.('error');
          this.callbacks.onError?.(error?.message || 'Voice conversation error');
        },
        
        onModeChange: (mode: any) => {
          console.log('Agent mode changed:', mode);
          // mode.mode can be 'speaking', 'listening', etc.
          const status: VoiceAgentStatus = mode.mode === 'speaking' ? 'speaking' : 'listening';
          this.callbacks.onAgentStatusChange?.(status);
        },
        
        onMessage: (message: any) => {
          console.log('📨 Message received (RAW):', {
            fullMessage: message,
            messageText: message.message,
            source: message.source,
            role: message.role,
            isFinal: message.isFinal,
            isFinalType: typeof message.isFinal,
            allKeys: Object.keys(message)
          });
          
          // Handle different message types from ElevenLabs
          // ElevenLabs sends: { message: string, source: 'user' | 'ai', role?: string }
          // isFinal might not always be present, especially for streaming responses
          if (message.message && message.message.trim()) {
            // More conservative isFinal detection:
            // - Only mark as final if explicitly set to true
            // - For AI responses (role: 'agent'), default to false (streaming)
            // - For user transcripts, use the isFinal flag
            const isUserMessage = message.source === 'user';
            const isFinal = isUserMessage 
              ? (message.isFinal === true) // User transcripts: respect the flag
              : (message.isFinal === true); // AI responses: only if explicitly true
            
            console.log('📤 Forwarding message:', {
              preview: message.message.substring(0, 50),
              source: message.source || 'ai',
              isFinal,
              messageLength: message.message.length
            });
            
            this.callbacks.onMessage?.({
              message: message.message,
              source: message.source || 'ai',
              isFinal
            });
          }
        },
        
        onDebug: (debugEvent: any) => {
          // Handle tentative/streaming agent responses for real-time display
          console.log('🔍 Debug event:', {
            type: debugEvent.type,
            hasResponse: !!debugEvent.response,
            event: debugEvent
          });
          
          // ElevenLabs sends tentative agent responses through debug events
          // type: 'tentative_agent_response' contains the streaming text AS IT'S BEING SPOKEN
          if (debugEvent.type === 'tentative_agent_response' && debugEvent.response) {
            console.log('📝 AI speaking chunk:', {
              preview: debugEvent.response.substring(0, 50),
              length: debugEvent.response.length
            });
            // This is text as it's being spoken - use onChunk to sync with voice
            this.callbacks.onChunk?.(debugEvent.response, 'ai');
          }
        },
        
        onVadScore: (vadData: any) => {
          // Voice Activity Detection score (0-1)
          // Higher score = more likely user is speaking
          if (this.callbacks.onVadScore && vadData.score !== undefined) {
            this.callbacks.onVadScore(vadData.score);
          }
        }
      });

      return true;
    } catch (error) {
      console.error('Failed to start conversation:', error);
      this.isActive = false;
      this.isConnecting = false;
      this.callbacks.onStatusChange?.('error');
      this.callbacks.onError?.(
        error instanceof Error ? error.message : 'Failed to start voice conversation'
      );
      return false;
    }
  }

  /**
   * Send a text message to the agent
   * This will make the agent respond with voice + transcript
   */
  sendTextMessage(text: string): boolean {
    if (!this.conversation || !this.isActive) {
      console.warn('Cannot send message: conversation not active', {
        hasConversation: !!this.conversation,
        isActive: this.isActive
      });
      return false;
    }

    try {
      console.log('Sending text message to voice agent:', text);
      this.conversation.sendUserMessage(text);
      console.log('Text message sent successfully');
      return true;
    } catch (error) {
      console.error('Failed to send text message:', error);
      this.callbacks.onError?.('Failed to send message');
      return false;
    }
  }

  /**
   * Stop the current conversation
   */
  async stopConversation(): Promise<void> {
    console.log('🛑 Stopping voice conversation (user requested)');
    this.shouldStayConnected = false; // User explicitly stopped, don't auto-reconnect
    this.reconnectAttempts = 0;
    
    if (this.conversation) {
      try {
        await this.conversation.endSession();
      } catch (error) {
        console.error('Error ending conversation:', error);
      }
      this.conversation = null;
    }
    this.isActive = false;
    this.isConnecting = false;
    this.callbacks.onStatusChange?.('disconnected');
    this.callbacks.onAgentStatusChange?.('idle');
  }

  /**
   * Clean up resources
   * WARNING: This should rarely be called - voice should persist across navigation
   * Only call this when explicitly stopping the app or logging out
   */
  cleanup(): void {
    console.warn('⚠️ cleanup() called - this should rarely happen');
    console.trace('Cleanup trace');
    this.shouldStayConnected = false;
    this.reconnectAttempts = 0;
    this.stopConversation();
    this.callbacks = {};
  }
}

// Export singleton instance
export const voiceConversationService = new VoiceConversationService();

