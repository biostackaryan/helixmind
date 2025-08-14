# helixmind/chatgpt_module.py
"""
HelixMind - Bioinformatics Toolkit
Author: Aryan Dutt (https://github.com/biostackaryan)
License: GNU GPL v3
Copyright (C) 2025 Aryan Dutt
"""

import os
import requests
from dotenv import load_dotenv, find_dotenv

# Load .env file automatically
load_dotenv(find_dotenv())

def ask_chatgpt(prompt, api_key=None, model="mistralai/Mistral-7B-Instruct-v0.2"):
    """
    Send a prompt to Together.ai API (default: Mistral model).
    If no API key is provided, it is read from TOGETHER_API_KEY env var.
    """
    
    # Get API key from param or env
    api_key = api_key or os.getenv("TOGETHER_API_KEY")
    if not api_key:
        return "Error: TOGETHER_API_KEY not set in environment or .env file"

    url = "https://api.together.xyz/v1/chat/completions"
    
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }

    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": "You are a helpful bioinformatics assistant."},
            {"role": "user", "content": prompt}
        ],
        "max_tokens": 512,
        "temperature": 0.7
    }

    try:
        response = requests.post(url, headers=headers, json=payload, timeout=30)
        response.raise_for_status()
        data = response.json()
        return data["choices"][0]["message"]["content"].strip()
    except requests.exceptions.RequestException as e:
        return f"Error: {str(e)}"
    except (KeyError, IndexError):
        return "Error: Unexpected API response format"
def ask_llm(prompt: str) -> str:
    """
    Dummy ChatGPT call placeholder.
    Replace with actual OpenAI API call if needed.
    """
    return f"LLM Response for: {prompt}"
