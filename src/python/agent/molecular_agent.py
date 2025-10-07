"""
Molecular Agent for AI-powered molecular analysis and assistance.
Integrates with Google Gemini API to provide intelligent responses for quantum chemistry queries.
"""

import os
import json
import logging
import time
from pathlib import Path
from typing import List, Dict, Any, Optional, Iterator
from google import genai
from google.genai import types
from google.genai.errors import APIError, ServerError
from quantum_calc.settings_manager import get_current_settings
from generated_models import HistoryItem

# Set up logging
logger = logging.getLogger(__name__)


def get_gemini_api_key() -> Optional[str]:
    """
    Get Gemini API key from environment variable or settings file.
    
    This is a utility function that can be used by any module that needs
    to access the Gemini API key.
    
    Priority:
        1. Environment variable: GEMINI_API_KEY
        2. Settings file: settings.gemini_api_key
        3. None (fallback)
    
    Returns:
        Optional[str]: API key if found, None otherwise
    """
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


class MolecularAgent:
    """AI agent for molecular analysis and quantum chemistry assistance."""

    def __init__(self):
        """Initialize the Molecular Agent with Gemini API and tools."""
        # Load agent configuration from config file
        self.config = self._load_agent_config()

        # Get API key with priority: environment variable -> settings file -> fallback
        self.api_key = self._get_api_key()
        self.client = None
        self.available_tools = []

        if not self.api_key:
            logger.warning("Gemini API key not found in environment variable or settings. Agent will use fallback responses.")
        else:
            try:
                # Initialize with latest google.genai SDK
                self.client = genai.Client(api_key=self.api_key)

                # Initialize tools
                self.available_tools = self._initialize_tools()

                logger.info("Gemini API client initialized successfully with tool integration")
                logger.debug(f"Available tools: {[tool.__name__ for tool in self.available_tools]}")
            except Exception as e:
                logger.error(f"Failed to initialize Gemini API client: {e}", exc_info=True)
                logger.debug(f"API key length: {len(self.api_key) if self.api_key else 0}")
                logger.debug(f"Error type: {type(e).__name__}")
                self.client = None
                self.available_tools = []
    
    def _get_api_key(self) -> Optional[str]:
        """Get Gemini API key from environment variable or settings file."""
        return get_gemini_api_key()

    def _load_agent_config(self) -> Dict[str, Any]:
        """Load AI agent configuration from server-config.json with fallback defaults."""
        # Default configuration
        default_config = {
            "model_name": "gemini-2.5-flash",
            "max_retries": 3,
            "retry_delays": [1, 2, 4],
            "retryable_status_codes": [429, 500, 502, 503, 504]
        }

        try:
            # Try to find config file (go up from agent directory to project root)
            config_path = Path(__file__).parent.parent.parent.parent / "config" / "server-config.json"

            if config_path.exists():
                with open(config_path, 'r') as f:
                    server_config = json.load(f)

                # Extract ai_agent section if it exists
                if "ai_agent" in server_config:
                    config = {**default_config, **server_config["ai_agent"]}
                    logger.debug(f"Loaded AI agent config from {config_path}")
                    return config
                else:
                    logger.info("No 'ai_agent' section in server-config.json, using defaults")
                    return default_config
            else:
                logger.warning(f"Config file not found at {config_path}, using default configuration")
                return default_config

        except Exception as e:
            logger.error(f"Failed to load agent configuration: {e}, using defaults")
            return default_config

    def _initialize_tools(self) -> List[Any]:
        """Initialize and return the list of available tools."""
        from . import tools

        return [
            # Calculation management
            tools.list_all_calculations,
            tools.get_calculation_details,
            tools.start_quantum_calculation,
            tools.delete_calculation,  # HIL-enabled destructive tool

            # Molecular structure tools
            tools.search_pubchem_by_name,
            tools.convert_smiles_to_xyz,
            tools.validate_xyz_format,

            # Molecular orbital analysis
            tools.get_molecular_orbitals,
            tools.generate_orbital_cube,
            tools.list_cube_files,
            tools.delete_cube_files,

            # Spectroscopy
            tools.generate_ir_spectrum,

            # System and settings
            tools.get_supported_parameters,
            tools.get_app_settings,
            tools.update_app_settings,
            tools.get_system_resources,
        ]

    def _convert_history_to_gemini_format(self, history: List[HistoryItem]) -> List[Dict[str, Any]]:
        """Convert frontend chat history to Gemini API format.
        
        Handles both Pydantic HistoryItem objects and plain dictionaries,
        ensuring compatibility with both direct API calls and LangGraph dispatcher.
        """
        gemini_history = []
        
        for item in history:
            # Handle both dict and Pydantic model objects
            if isinstance(item, dict):
                role = item.get("role", "user")
                parts = item.get("parts", [])
            else:
                # Pydantic model (HistoryItem)
                role = item.role.value if item.role else 'user'
                parts = item.parts or []
            
            # Extract text from all parts
            text_content = ""
            for part in parts:
                if isinstance(part, dict) and "text" in part:
                    text_content += part["text"]
                elif hasattr(part, "text") and part.text:
                    text_content += part.text
            
            # Only add non-empty messages
            if text_content.strip():
                gemini_history.append({
                    "role": role,
                    "parts": [{"text": text_content.strip()}]
                })
        
        return gemini_history
    
    def _get_system_prompt(self) -> str:
        """Get the system prompt for molecular analysis context with ReAct framework."""
        try:
            # Try to load system prompt from file
            prompt_path = Path(__file__).parent / "prompts" / "system_prompt.txt"

            if prompt_path.exists():
                with open(prompt_path, 'r', encoding='utf-8') as f:
                    prompt = f.read()
                logger.debug(f"Loaded system prompt from {prompt_path}")
                return prompt
            else:
                logger.warning(f"System prompt file not found at {prompt_path}, using fallback")
                return self._get_fallback_system_prompt()

        except Exception as e:
            logger.error(f"Failed to load system prompt from file: {e}, using fallback")
            return self._get_fallback_system_prompt()

    def _get_fallback_system_prompt(self) -> str:
        """Fallback system prompt in case file loading fails."""
        return "You are PySCF Agent, an AI assistant specialized in molecular analysis and quantum chemistry. You have access to tools for calculation management, molecular structure analysis, and more. Use them proactively to help users."
    
    def reload_api_key(self) -> bool:
        """
        Reload the API key from environment variable or settings file.
        Returns True if API key was found and model was initialized successfully.
        """
        new_api_key = self._get_api_key()

        # Only reinitialize if the API key has changed
        if new_api_key != self.api_key:
            self.api_key = new_api_key

            if not self.api_key:
                logger.warning("Gemini API key not found after reload. Agent will use fallback responses.")
                self.client = None
                self.available_tools = []
                return False
            else:
                try:
                    self.client = genai.Client(api_key=self.api_key)

                    # Reinitialize tools using centralized method
                    self.available_tools = self._initialize_tools()

                    logger.info("Gemini API client reinitialized successfully after settings update")
                    return True
                except Exception as e:
                    logger.error(f"Failed to reinitialize Gemini API client after reload: {e}", exc_info=True)
                    logger.debug(f"New API key length: {len(new_api_key) if new_api_key else 0}")
                    logger.debug(f"Reload error type: {type(e).__name__}")
                    self.client = None
                    self.available_tools = []
                    return False

        # No change needed, return current status
        return self.client is not None

    def chat(self, message: str, history: List[Dict[str, Any]]) -> Iterator[str]:
        """
        Chat with AI agent using streaming with Function Calling support.

        Args:
            message (str): User's message
            history (List[Dict[str, Any]]): Chat history in frontend format

        Yields:
            str: AI agent's response chunks
        """
        # Load retry configuration from config
        max_retries = self.config.get("max_retries", 3)
        retry_delays = self.config.get("retry_delays", [1, 2, 4])
        retryable_status_codes = self.config.get("retryable_status_codes", [429, 500, 502, 503, 504])

        for attempt in range(max_retries):
            try:
                logger.debug(f"Starting Function Calling chat (attempt {attempt + 1}/{max_retries}) with message length: {len(message)}, history entries: {len(history)}")
                
                # If Gemini API is not available, provide fallback response
                if not self.client or not self.available_tools:
                    logger.warning("Gemini client or tools not available, using fallback response")
                    yield self._get_fallback_response(message)
                    return
                
                # Convert history to Gemini format
                gemini_history = self._convert_history_to_gemini_format(history)
                logger.debug(f"Converted history to Gemini format: {len(gemini_history)} entries")
                
                # Add new user message in proper format
                full_conversation = gemini_history + [{"role": "user", "parts": [{"text": message}]}]
                
                # Use system_instruction parameter for system prompt
                system_instruction = self._get_system_prompt()
                
                logger.debug(f"Calling Gemini API with {len(full_conversation)} conversation entries and {len(self.available_tools)} tools")
                
                # Use new SDK with Function Calling
                model_name = self.config.get("model_name", "gemini-2.5-flash")
                response_stream = self.client.models.generate_content_stream(
                    model=model_name,
                    contents=full_conversation,
                    config=types.GenerateContentConfig(
                        system_instruction=system_instruction,
                        tools=self.available_tools,  # Pass Python functions directly
                    )
                )
                
                # Stream response chunks
                chunk_count = 0
                for chunk in response_stream:
                    if chunk.text:
                        chunk_count += 1
                        logger.debug(f"Yielding chunk #{chunk_count}, length: {len(chunk.text)}")
                        yield chunk.text
                        
                logger.debug(f"Completed streaming with {chunk_count} chunks")
                return  # Success - exit retry loop
                    
            except ServerError as e:
                # Extract status code from error
                status_code = getattr(e, 'code', None)
                is_retryable = status_code in retryable_status_codes
                is_last_attempt = attempt == max_retries - 1

                # Log appropriately based on error type
                if is_retryable and not is_last_attempt:
                    delay = retry_delays[attempt] if attempt < len(retry_delays) else retry_delays[-1]
                    logger.warning(
                        f"Gemini API returned temporary error {status_code} (attempt {attempt + 1}/{max_retries}). "
                        f"Retrying in {delay}s... Error: {e}"
                    )
                    time.sleep(delay)
                    continue  # Retry
                else:
                    # Last attempt or non-retryable error
                    logger.error(f"Gemini API error: {status_code} {e}")
                    yield self._get_error_response(e)
                    return
                    
            except APIError as e:
                # Handle other API errors
                logger.error(f"Gemini API error: {e}")
                yield self._get_error_response(e)
                return
                
            except Exception as e:
                logger.error(f"Unexpected error during Gemini API stream call: {e}", exc_info=True)
                
                # Log additional context for debugging
                logger.debug(f"Message that caused error: {message}")
                logger.debug(f"History length: {len(history)}")
                
                yield self._get_error_response(e)
                return
    
    def _get_fallback_response(self, message: str) -> str:
        """Provide a fallback response when Gemini API is not available."""
        return (
            f"I received your message: '{message}'. "
            "However, the AI agent is currently not fully configured. "
            "To enable full AI capabilities, please set the GEMINI_API_KEY environment variable. "
            "In the meantime, I can help you with basic information about the PySCF application features."
        )
    
    def _get_error_response(self, error: Exception) -> str:
        """Provide an appropriate error response based on the type of error."""
        error_str = str(error).lower()
        
        # Handle ServerError (503, 502, etc.) specifically
        if isinstance(error, ServerError):
            status_code = getattr(error, 'code', None)
            
            if status_code == 503:
                return (
                    "🔄 The AI service is currently experiencing high demand. "
                    "The system automatically retried your request, but the service is still overloaded. "
                    "Please wait a few moments and try again."
                )
            elif status_code == 502:
                return (
                    "🔌 The AI service is experiencing connectivity issues. "
                    "Please try again in a moment."
                )
            elif status_code == 500:
                return (
                    "⚠️ The AI service is experiencing internal issues. "
                    "Please try again later."
                )
            elif status_code == 429:
                return (
                    "⏱️ Too many requests have been sent to the AI service. "
                    "Please wait a moment before trying again."
                )
        
        # Define error patterns and their corresponding messages
        error_patterns = {
            ("api key", "authentication", "unauthorized"): (
                "🔑 There was an authentication issue with the AI service. "
                "Please check that your API key is correctly configured in Settings."
            ),
            ("quota", "limit"): (
                "📊 The AI service quota has been reached. "
                "Please try again later or check your API quota."
            ),
            ("network", "connection", "timeout"): (
                "🌐 There was a network connectivity issue. "
                "Please check your internet connection and try again."
            ),
        }
        
        # Check error patterns
        for patterns, message in error_patterns.items():
            if any(pattern in error_str for pattern in patterns):
                return message
        
        # Handle generic APIError
        if isinstance(error, APIError):
            status_code = getattr(error, 'code', None)
            if status_code:
                return f"⚠️ The AI service returned an error (code {status_code}). Please try again later."
        
        # Generic error response
        logger.error(f"Unexpected error in AI agent: {error}", exc_info=True)
        return (
            "❌ I'm experiencing a temporary issue and cannot process your request right now. "
            "Please try again later."
        )