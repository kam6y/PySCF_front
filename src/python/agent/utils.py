"""
Common utility functions for agent modules.

This module provides shared functionality across all agent implementations:
- API key management
- Language detection
- Prompt loading
"""

import logging
from pathlib import Path
from typing import Optional
from langchain_google_genai import ChatGoogleGenerativeAI

# Set up logging
logger = logging.getLogger(__name__)


def get_gemini_api_key() -> Optional[str]:
    """
    Get Gemini API key from environment variable or settings file.

    Priority:
        1. Environment variable: GEMINI_API_KEY
        2. Settings file: settings.gemini_api_key
        3. None (fallback)

    Returns:
        Optional[str]: API key if found, None otherwise
    """
    import os
    from quantum_calc.settings_manager import get_current_settings

    # Priority 1: Environment variable
    api_key = os.environ.get("GEMINI_API_KEY")
    if api_key:
        logger.debug("Using Gemini API key from environment variable")
        return api_key

    # Priority 2: Settings file
    try:
        settings = get_current_settings()
        if settings.gemini_api_key:
            logger.debug("Using Gemini API key from settings file")
            return settings.gemini_api_key
    except Exception as e:
        logger.warning(f"Failed to load settings for API key: {e}")

    # Priority 3: None (fallback)
    return None


def detect_language(text: str, llm: ChatGoogleGenerativeAI) -> str:
    """
    Detect the language of the given text using LLM.

    Args:
        text: Text to analyze
        llm: LLM instance for detection

    Returns:
        Language code (e.g., 'ja', 'en', 'es', 'zh', 'fr')
    """
    try:
        prompt = f"""Detect the language of the following text and return ONLY the ISO 639-1 language code (2 letters).

Examples:
- Japanese text → ja
- English text → en
- Spanish text → es
- Chinese text → zh
- French text → fr

Text: {text}

Return ONLY the 2-letter language code, nothing else:"""

        response = llm.invoke(prompt)
        lang_code = response.content.strip().lower()

        # Validate it's a 2-letter code
        if len(lang_code) == 2 and lang_code.isalpha():
            logger.info(f"Detected language: {lang_code}")
            return lang_code
        else:
            logger.warning(f"Invalid language code detected: {lang_code}, defaulting to 'en'")
            return 'en'
    except Exception as e:
        logger.error(f"Error detecting language: {e}, defaulting to 'en'")
        return 'en'


def load_prompt(prompt_dir: Path, filename: str) -> str:
    """
    Load a prompt from a prompts directory.

    Args:
        prompt_dir: Path to the prompts directory
        filename: Name of the prompt file to load

    Returns:
        str: Prompt content, or empty string if file not found
    """
    try:
        prompt_path = prompt_dir / filename
        if prompt_path.exists():
            with open(prompt_path, 'r', encoding='utf-8') as f:
                content = f.read()
            logger.debug(f"Loaded prompt from {prompt_path}")
            return content
        else:
            logger.warning(f"Prompt file not found: {prompt_path}")
            return ""
    except Exception as e:
        logger.error(f"Failed to load prompt {filename}: {e}")
        return ""


# Language name mapping for clarity in prompts
LANGUAGE_NAMES = {
    'ja': 'Japanese (日本語)',
    'en': 'English',
    'es': 'Spanish (Español)',
    'zh': 'Chinese (中文)',
    'fr': 'French (Français)',
    'de': 'German (Deutsch)',
    'ko': 'Korean (한국어)',
    'ru': 'Russian (Русский)',
    'pt': 'Portuguese (Português)',
    'it': 'Italian (Italiano)'
}


def get_language_name(lang_code: str) -> str:
    """
    Get the full language name from a language code.

    Args:
        lang_code: ISO 639-1 language code (e.g., 'ja', 'en')

    Returns:
        str: Full language name with native script, or the code itself if not found
    """
    return LANGUAGE_NAMES.get(lang_code, lang_code)
