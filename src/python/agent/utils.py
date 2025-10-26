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


def get_tavily_api_key() -> Optional[str]:
    """
    Get Tavily API key from environment variable or settings file.

    Priority:
        1. Environment variable: TAVILY_API_KEY
        2. Settings file: settings.tavily_api_key
        3. None (fallback)

    Returns:
        Optional[str]: API key if found, None otherwise
    """
    import os
    from quantum_calc.settings_manager import get_current_settings

    # Priority 1: Environment variable
    api_key = os.environ.get("TAVILY_API_KEY")
    if api_key:
        logger.debug("Using Tavily API key from environment variable")
        return api_key

    # Priority 2: Settings file
    try:
        settings = get_current_settings()
        if settings.tavily_api_key:
            logger.debug("Using Tavily API key from settings file")
            return settings.tavily_api_key
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


# Language detection cache: {text_hash: (language_code, timestamp)}
_language_cache: dict[int, tuple[str, float]] = {}
_CACHE_DURATION_SECONDS = 3600  # 1 hour


def detect_language_cached(text: str, llm: ChatGoogleGenerativeAI, cache_duration: int = _CACHE_DURATION_SECONDS) -> str:
    """
    Detect the language of the given text using LLM, with memory caching.

    This function caches language detection results to avoid redundant LLM API calls
    within the same session. The cache is based on a hash of the input text and has
    a configurable expiration time.

    Args:
        text: Text to analyze
        llm: LLM instance for detection
        cache_duration: Cache validity duration in seconds (default: 3600 = 1 hour)

    Returns:
        Language code (e.g., 'ja', 'en', 'es', 'zh', 'fr')
    """
    import time

    # Generate cache key from text hash
    text_hash = hash(text.lower().strip())

    # Check cache
    if text_hash in _language_cache:
        cached_lang, cached_time = _language_cache[text_hash]
        time_elapsed = time.time() - cached_time

        # Return cached result if still valid
        if time_elapsed < cache_duration:
            logger.debug(f"Language detection cache hit: {cached_lang} (cached {time_elapsed:.1f}s ago)")
            return cached_lang
        else:
            logger.debug(f"Language detection cache expired (age: {time_elapsed:.1f}s)")

    # Cache miss or expired - perform actual detection
    logger.debug("Language detection cache miss, calling LLM")
    detected_lang = detect_language(text, llm)

    # Update cache
    _language_cache[text_hash] = (detected_lang, time.time())

    # Clean up old cache entries (keep only last 100)
    if len(_language_cache) > 100:
        # Remove oldest entries
        sorted_entries = sorted(_language_cache.items(), key=lambda x: x[1][1])
        for old_hash, _ in sorted_entries[:-50]:  # Keep newest 50
            del _language_cache[old_hash]
        logger.debug("Language detection cache cleaned up")

    return detected_lang


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


def extract_language_from_messages(messages: list) -> str:
    """
    Extract language code from messages inserted by Supervisor.

    The Supervisor inserts a SystemMessage with format "DETECTED_LANGUAGE:xx"
    to share the detected language with all workers, avoiding redundant
    language detection LLM calls.

    Args:
        messages: List of messages (MessagesState)

    Returns:
        str: Language code (e.g., 'ja', 'en'), defaults to 'en' if not found
    """
    from langchain_core.messages import SystemMessage

    # Search for DETECTED_LANGUAGE marker in SystemMessages
    for msg in messages:
        if isinstance(msg, SystemMessage):
            content = msg.content
            if isinstance(content, str) and content.startswith("DETECTED_LANGUAGE:"):
                lang_code = content.split(":", 1)[1].strip()
                logger.debug(f"Extracted language from Supervisor: {lang_code}")
                return lang_code

    # Default to English if no language marker found
    logger.debug("No language marker found in messages, defaulting to 'en'")
    return 'en'
