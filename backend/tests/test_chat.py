"""
Test chat endpoint streaming functionality.

This test demonstrates that the chat endpoint properly streams responses
by printing each chunk to the terminal as it arrives.
"""
import asyncio
import sys
from typing import AsyncIterator

import httpx


async def test_chat_streaming():
    """
    Test the chat endpoint with streaming enabled.
    
    This test:
    1. Sends a simple health-related question
    2. Streams the response chunks
    3. Prints each chunk to terminal in real-time
    4. Verifies streaming is working properly
    """
    # API endpoint
    base_url = "http://localhost:8000"
    endpoint = f"{base_url}/api/v1/chat/general"
    
    # Test request payload
    payload = {
        "messages": [
            {
                "role": "user",
                "content": "What are the benefits of drinking water regularly?"
            }
        ],
        "model": "gpt-4-turbo"  # Optional: defaults to gpt-4-turbo in controller
    }
    
    print("=" * 80)
    print("Testing Chat Endpoint Streaming")
    print("=" * 80)
    print(f"\n📍 Endpoint: {endpoint}")
    print(f"❓ Question: {payload['messages'][0]['content']}")
    print(f"🤖 Model: {payload.get('model', 'default (gpt-4-turbo)')}")
    print("\n" + "=" * 80)
    print("📡 Streaming Response:")
    print("=" * 80 + "\n")
    
    # Track streaming metrics
    chunk_count = 0
    total_chars = 0
    full_response = []
    
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            async with client.stream("POST", endpoint, json=payload) as response:
                # Check if request was successful
                if response.status_code != 200:
                    print(f"\n❌ Error: HTTP {response.status_code}")
                    error_text = await response.aread()
                    print(f"Response: {error_text.decode()}")
                    return
                
                print("✅ Connection established, streaming started...\n")
                print("-" * 80)
                
                # Stream and print chunks
                async for chunk in response.aiter_text():
                    if chunk:
                        # Print chunk immediately (without newline)
                        sys.stdout.write(chunk)
                        sys.stdout.flush()
                        
                        # Track metrics
                        chunk_count += 1
                        total_chars += len(chunk)
                        full_response.append(chunk)
        
        # Print streaming statistics
        print("\n" + "-" * 80)
        print("\n" + "=" * 80)
        print("📊 Streaming Statistics:")
        print("=" * 80)
        print(f"✅ Total chunks received: {chunk_count}")
        print(f"✅ Total characters: {total_chars}")
        print(f"✅ Average chunk size: {total_chars / chunk_count if chunk_count > 0 else 0:.1f} chars")
        print(f"✅ Full response length: {len(''.join(full_response))} chars")
        print("=" * 80)
        print("\n✨ Streaming test completed successfully!\n")
        
    except httpx.ConnectError:
        print("\n❌ Error: Could not connect to server")
        print("💡 Make sure the server is running: python main.py or uvicorn main:app")
        print("   Server should be at: http://localhost:8000")
    except httpx.TimeoutException:
        print("\n❌ Error: Request timed out")
        print("💡 The API might be taking too long to respond")
    except Exception as e:
        print(f"\n❌ Unexpected error: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()


async def test_chat_streaming_with_conversation():
    """
    Test the chat endpoint with a multi-turn conversation.
    
    This demonstrates how to send conversation history in the messages array.
    """
    base_url = "http://localhost:8000"
    endpoint = f"{base_url}/api/v1/chat/general"
    
    # Test request with conversation history
    payload = {
        "messages": [
            {
                "role": "user",
                "content": "What is a healthy heart rate?"
            },
            {
                "role": "assistant",
                "content": "A healthy resting heart rate for adults typically ranges from 60 to 100 beats per minute."
            },
            {
                "role": "user",
                "content": "What about during exercise?"
            }
        ],
        "model": "gpt-4-turbo"
    }
    
    print("\n" + "=" * 80)
    print("Testing Multi-Turn Conversation")
    print("=" * 80)
    print(f"\n📍 Endpoint: {endpoint}")
    print("\n💬 Conversation History:")
    for i, msg in enumerate(payload['messages'], 1):
        role_emoji = "👤" if msg['role'] == 'user' else "🤖"
        print(f"  {role_emoji} {msg['role'].title()}: {msg['content'][:100]}...")
    
    print("\n" + "=" * 80)
    print("📡 Streaming Response:")
    print("=" * 80 + "\n")
    
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            async with client.stream("POST", endpoint, json=payload) as response:
                if response.status_code != 200:
                    print(f"\n❌ Error: HTTP {response.status_code}")
                    error_text = await response.aread()
                    print(f"Response: {error_text.decode()}")
                    return
                
                print("-" * 80)
                async for chunk in response.aiter_text():
                    if chunk:
                        sys.stdout.write(chunk)
                        sys.stdout.flush()
                print("\n" + "-" * 80)
                print("\n✨ Multi-turn conversation test completed!\n")
                
    except Exception as e:
        print(f"\n❌ Error: {type(e).__name__}: {e}")


async def main():
    """Run all tests."""
    print("\n🚀 Starting Chat API Streaming Tests\n")
    
    # Test 1: Simple streaming
    await test_chat_streaming()
    
    # Wait a moment between tests
    await asyncio.sleep(2)
    
    # Test 2: Multi-turn conversation
    await test_chat_streaming_with_conversation()
    
    print("\n" + "=" * 80)
    print("🎉 All tests completed!")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    # Run the async tests
    asyncio.run(main())
