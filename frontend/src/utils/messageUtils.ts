export type BaseMessage = {
  id: string;
  type: 'admin' | 'user';
  content: string;
  timestamp: Date;
  senderName?: string;
  deliveryStatus?: string;
  audioTranscript?: string;
  files?: { id: string; [key: string]: any }[];
  twilioSid?: string;
  deliveryUpdatedAt?: Date;
};

export const normalizeMessagesFromApi = (items: any[]): BaseMessage[] => {
  return (items || [])
    .map((m: any) => ({
      id: String(m.message_id),
      type: m.sender_type as 'admin' | 'user',
      content: m.message_body || m.audio_transcript || 'No content',
      timestamp: new Date(m.created_at),
      senderName: m.sender_name,
      deliveryStatus: m.delivery_status,
      audioTranscript: m.audio_transcript,
      files: m.files || [],
      twilioSid: m.twilio_message_sid,
      deliveryUpdatedAt: m.delivery_updated_at ? new Date(m.delivery_updated_at) : undefined,
    }))
    .sort((a: BaseMessage, b: BaseMessage) => a.timestamp.getTime() - b.timestamp.getTime());
};


