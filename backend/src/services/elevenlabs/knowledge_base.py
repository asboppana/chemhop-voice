from elevenlabs.client import ElevenLabs
from dotenv import load_dotenv
import os

def get_knowledge_base():
    """
    Get the knowledge base for the ElevenLabs voice agent.
    """

    load_dotenv()

    client = ElevenLabs(
        api_key=os.getenv("ELEVENLABS_API_KEY"),
    )

    # --------------------
    # Add documents to knowledge base

    doc1 = client.conversational_ai.knowledge_base.documents.create_from_url(
        url="https://en.wikipedia.org/wiki/Drug_discovery",
        name="Drug Discovery Wikipedia page",
    )

    # --------------------
    # Concat all documents into the knowledge base

    knowledge_base = [
        {
            "type": "url",
            "name": doc1.name,
            "id": doc1.id,
            "usage": "auto"
        }
    ]

    return knowledge_base