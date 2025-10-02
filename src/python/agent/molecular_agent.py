"""
Molecular Agent for AI-powered molecular analysis and assistance.
Integrates with Google Gemini API to provide intelligent responses for quantum chemistry queries.
"""

import os
import logging
import time
from typing import List, Dict, Any, Optional, Iterator
from google import genai
from google.genai import types
from google.genai.errors import APIError, ServerError
from quantum_calc.settings_manager import get_current_settings
from generated_models import HistoryItem

# Set up logging
logger = logging.getLogger(__name__)


class MolecularAgent:
    """AI agent for molecular analysis and quantum chemistry assistance."""
    
    def __init__(self):
        """Initialize the Molecular Agent with Gemini API and tools."""
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
                
                # Import and set up available tools
                from . import tools
                self.available_tools = [
                    tools.list_all_calculations,
                    tools.get_calculation_details,
                    tools.get_supported_parameters,
                    tools.start_quantum_calculation,
                    tools.search_pubchem_by_name,
                    tools.delete_calculation,  # HIL-enabled destructive tool
                ]
                
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
    
    def _convert_history_to_gemini_format(self, history: List[HistoryItem]) -> List[Dict[str, Any]]:
        """Convert frontend chat history to Gemini API format."""
        gemini_history = []
        
        for item in history:
            # Extract role (default to 'user')
            role = item.role.value if item.role else 'user'
            
            # Extract text from first part if available
            text_content = ""
            if item.parts and len(item.parts) > 0 and item.parts[0].text:
                text_content = item.parts[0].text.strip()
            
            # Only add non-empty messages
            if text_content:
                gemini_history.append({
                    "role": role,
                    "parts": [{"text": text_content}]
                })
        
        return gemini_history
    
    def _get_system_prompt(self) -> str:
        """Get the system prompt for molecular analysis context with ReAct framework."""
        return """You are an intelligent AI assistant specialized in molecular analysis and quantum chemistry, called "PySCF Agent". 
You can help users with quantum chemistry calculations, molecular structure analysis, and related scientific tasks.

IMPORTANT: You have access to powerful tools that can interact with the PySCF application backend. Use these tools actively to help users achieve their goals.

**Available Tools:**
- `list_all_calculations`: Retrieve all saved quantum chemistry calculations
- `get_calculation_details`: Get detailed results for a specific calculation ID
- `get_supported_parameters`: Get available calculation methods, basis sets, functionals, and solvents
- `start_quantum_calculation`: Start a new quantum chemistry calculation with molecular XYZ data
- `search_pubchem_by_name`: Search PubChem database for molecular structures by compound name
- `delete_calculation`: **[DESTRUCTIVE]** Request deletion of a calculation (requires user confirmation)

**ReAct Framework Instructions:**
Before taking any action, you MUST follow this pattern:

1. **[Thought]**: Analyze the user's request and plan your approach. Think step-by-step about what tools you need to use and why.
2. **[Action]**: Use the appropriate tools to gather information or perform tasks.
3. **[Observation]**: Analyze the results from your tool usage.
4. **[Response]**: Provide a comprehensive answer to the user based on your observations.

**Guidelines:**
- Always respond in clear, professional English
- Use tools proactively when they can help answer the user's question
- If a tool returns an error, explain what happened and suggest alternatives
- For calculation-related questions, always check the calculation list first
- Provide scientific explanations when discussing quantum chemistry results
- Be helpful, accurate, and educational in your responses

**CRITICAL - Destructive Operations:**
- The `delete_calculation` tool returns a JSON confirmation request, NOT a deletion result
- When you receive a confirmation request JSON (with "requires_confirmation": true), you MUST:
  1. Include the COMPLETE JSON in your response as a ```json code block FIRST
  2. Then provide a natural language explanation to the user
  3. Example format:
     ```json
     {{
       "requires_confirmation": true,
       "action": "delete_calculation",
       "calculation_id": "calc_id_here",
       "calculation_name": "calc name here",
       "message": "confirmation message here"
     }}
     ```
     
     [Your natural language explanation about the deletion request]
- NEVER tell the user that you "deleted" something when you only requested confirmation
- The actual deletion happens only after user approval through the UI confirmation modal
- The frontend will automatically parse the JSON block and display the confirmation modal

Remember: You are not just a chatbot - you are an active assistant that can access real application data and functionality to help users with their molecular analysis workflows."""
    
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
                    
                    # Reinitialize tools
                    from . import tools
                    self.available_tools = [
                    tools.list_all_calculations,
                    tools.get_calculation_details,
                    tools.get_supported_parameters,
                    tools.start_quantum_calculation,
                    tools.search_pubchem_by_name,
                    tools.delete_calculation,  # HIL-enabled destructive tool
                ]
                    
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
        # Retry configuration
        MAX_RETRIES = 3
        RETRY_DELAYS = [1, 2, 4]  # Exponential backoff in seconds
        RETRYABLE_STATUS_CODES = [429, 500, 502, 503, 504]  # Temporary errors
        
        for attempt in range(MAX_RETRIES):
            try:
                logger.debug(f"Starting Function Calling chat (attempt {attempt + 1}/{MAX_RETRIES}) with message length: {len(message)}, history entries: {len(history)}")
                
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
                response_stream = self.client.models.generate_content_stream(
                    model='gemini-2.5-pro',
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
                is_retryable = status_code in RETRYABLE_STATUS_CODES
                is_last_attempt = attempt == MAX_RETRIES - 1
                
                # Log appropriately based on error type
                if is_retryable and not is_last_attempt:
                    logger.warning(
                        f"Gemini API returned temporary error {status_code} (attempt {attempt + 1}/{MAX_RETRIES}). "
                        f"Retrying in {RETRY_DELAYS[attempt]}s... Error: {e}"
                    )
                    time.sleep(RETRY_DELAYS[attempt])
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