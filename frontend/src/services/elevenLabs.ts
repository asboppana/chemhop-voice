// Simple Web Speech API wrapper for continuous voice recognition
export interface VoiceASRCallbacks {
  onTranscript: (text: string, isFinal: boolean) => void;
  onError: (error: string) => void;
}

class VoiceASRService {
  private recognition: any = null;
  private callbacks: VoiceASRCallbacks | null = null;
  private isListening: boolean = false;
  private finalTranscriptTimer: number | null = null;
  private currentTranscript: string = '';

  setCallbacks(callbacks: VoiceASRCallbacks) {
    this.callbacks = callbacks;
  }

  startListening(): boolean {
    // Check if browser supports Web Speech API
    const SpeechRecognition = (window as any).SpeechRecognition || (window as any).webkitSpeechRecognition;
    
    if (!SpeechRecognition) {
      this.callbacks?.onError('Speech recognition not supported in this browser');
      return false;
    }

    try {
      this.recognition = new SpeechRecognition();
      this.recognition.continuous = true; // Keep listening
      this.recognition.interimResults = true; // Get partial results
      this.recognition.lang = 'en-US';

      this.recognition.onresult = (event: any) => {
        let interimTranscript = '';
        let finalTranscript = '';

        // Process all results
        for (let i = event.resultIndex; i < event.results.length; i++) {
          const transcript = event.results[i][0].transcript;
          if (event.results[i].isFinal) {
            finalTranscript += transcript + ' ';
          } else {
            interimTranscript += transcript;
          }
        }

        // Update current transcript
        if (finalTranscript) {
          this.currentTranscript += finalTranscript;
          this.callbacks?.onTranscript(this.currentTranscript.trim(), false);
          
          // Reset timer - wait 5 seconds of silence before finalizing
          if (this.finalTranscriptTimer) {
            clearTimeout(this.finalTranscriptTimer);
          }
          this.finalTranscriptTimer = setTimeout(() => {
            if (this.currentTranscript.trim()) {
              console.log('🎤 Finalizing after 3s pause:', this.currentTranscript.trim());
              this.callbacks?.onTranscript(this.currentTranscript.trim(), true);
              this.currentTranscript = '';
            }
          }, 3000); // 3 second pause
        } else if (interimTranscript) {
          // Show streaming interim results
          const fullText = (this.currentTranscript + ' ' + interimTranscript).trim();
          this.callbacks?.onTranscript(fullText, false);
        }
      };

      this.recognition.onerror = (event: any) => {
        console.error('Speech recognition error:', event.error);
        this.callbacks?.onError(event.error);
        
        // Auto-restart on certain errors
        if (event.error === 'no-speech' || event.error === 'audio-capture') {
          setTimeout(() => {
            if (this.isListening) {
              this.recognition?.start();
            }
          }, 1000);
        }
      };

      this.recognition.onend = () => {
        console.log('Speech recognition ended');
        // Auto-restart if we're supposed to be listening
        if (this.isListening) {
          console.log('Auto-restarting recognition...');
          setTimeout(() => {
            if (this.isListening) {
              this.recognition?.start();
            }
          }, 100);
        }
      };

      this.recognition.start();
      this.isListening = true;
      console.log('✅ Voice recognition started');
      return true;
    } catch (error) {
      console.error('Failed to start recognition:', error);
      this.callbacks?.onError(String(error));
      return false;
    }
  }

  stopListening(): void {
    if (this.recognition) {
      this.isListening = false;
      this.recognition.stop();
      this.recognition = null;
      
      // Clear any pending timers
      if (this.finalTranscriptTimer) {
        clearTimeout(this.finalTranscriptTimer);
        this.finalTranscriptTimer = null;
      }
      
      // Finalize any remaining transcript
      if (this.currentTranscript.trim()) {
        this.callbacks?.onTranscript(this.currentTranscript.trim(), true);
        this.currentTranscript = '';
      }
      
      console.log('❌ Voice recognition stopped');
    }
  }

  isActive(): boolean {
    return this.isListening;
  }

  // Manual finalize for "send" button
  finalizeCurrentTranscript(): void {
    if (this.finalTranscriptTimer) {
      clearTimeout(this.finalTranscriptTimer);
      this.finalTranscriptTimer = null;
    }
    
    if (this.currentTranscript.trim()) {
      this.callbacks?.onTranscript(this.currentTranscript.trim(), true);
      this.currentTranscript = '';
    }
  }
}

export const voiceASR = new VoiceASRService();