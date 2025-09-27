"""
Molecular Agent for AI-powered molecular analysis and assistance.
Integrates with Google Gemini API to provide intelligent responses for quantum chemistry queries.
"""

import os
import logging
from typing import List, Dict, Any, Optional, Iterator
import google.generativeai as genai
from quantum_calc.settings_manager import get_current_settings

# Set up logging
logger = logging.getLogger(__name__)


class MolecularAgent:
    """AI agent for molecular analysis and quantum chemistry assistance."""
    
    def __init__(self):
        """Initialize the Molecular Agent with Gemini API."""
        # Get API key with priority: environment variable -> settings file -> fallback
        self.api_key = self._get_api_key()
        if not self.api_key:
            logger.warning("Gemini API key not found in environment variable or settings. Agent will use fallback responses.")
            self.model = None
        else:
            try:
                genai.configure(api_key=self.api_key)
                self.model = genai.GenerativeModel('gemini-2.5-pro')
                logger.info("Gemini API initialized successfully")
            except Exception as e:
                logger.error(f"Failed to initialize Gemini API: {e}")
                self.model = None
    
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
    
    def _convert_history_to_gemini_format(self, history: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Convert frontend chat history to Gemini API format."""
        gemini_history = []
        
        for item in history:
            role = item.get('role', 'user')
            parts = item.get('parts', [])
            
            # Convert role mapping: 'model' -> 'model', 'user' -> 'user'
            gemini_role = 'model' if role == 'model' else 'user'
            
            # Extract text from parts
            text_content = ""
            if parts and len(parts) > 0:
                text_content = parts[0].get('text', '')
            
            if text_content.strip():  # Only add non-empty messages
                gemini_history.append({
                    "role": gemini_role,
                    "parts": [{"text": text_content}]
                })
        
        return gemini_history
    
    def _get_system_prompt(self) -> str:
        """Get the system prompt for molecular analysis context."""
        return """You are a helpful AI assistant specialized in molecular analysis and quantum chemistry. 
You can help with:
- Molecular structure analysis and interpretation
- Quantum chemistry calculation guidance (DFT, HF, MP2, CCSD, TDDFT, CASCI, CASSCF)
- Basis set and functional recommendations
- Interpretation of calculation results
- Molecular orbital analysis
- IR spectrum interpretation
- General chemistry and physics questions related to molecules

Please provide clear, accurate, and helpful responses while being concise and focused on the scientific aspects."""
    
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
                self.model = None
                return False
            else:
                try:
                    genai.configure(api_key=self.api_key)
                    self.model = genai.GenerativeModel('gemini-2.5-pro')
                    logger.info("Gemini API reinitialized successfully after settings update")
                    return True
                except Exception as e:
                    logger.error(f"Failed to reinitialize Gemini API after reload: {e}")
                    self.model = None
                    return False
        
        # No change needed, return current status
        return self.model is not None

    def chat(self, message: str, history: List[Dict[str, Any]]) -> Iterator[str]:
        """
        AIエージェントとストリーミングでチャットを行うジェネレータ関数。
        
        Args:
            message (str): User's message
            history (List[Dict[str, Any]]): Chat history in frontend format
            
        Yields:
            str: AI agent's response chunks
        """
        try:
            # If Gemini API is not available, provide fallback response
            if not self.model:
                yield self._get_fallback_response(message)
                return
            
            # Convert history to Gemini format
            gemini_history = self._convert_history_to_gemini_format(history)
            
            # Add system context if this is the first message
            if not gemini_history:
                system_message = {
                    "role": "user",
                    "parts": [{"text": f"{self._get_system_prompt()}\n\nUser question: {message}"}]
                }
                
                # stream=True を指定してストリーミング応答を開始
                response_stream = self.model.generate_content([system_message], stream=True)
                
                # ストリームからチャンクを一つずつyieldする
                for chunk in response_stream:
                    if chunk.text:
                        yield chunk.text
            else:
                # Start chat with history and enable streaming
                chat_session = self.model.start_chat(history=gemini_history)
                response_stream = chat_session.send_message(message, stream=True)
                
                # ストリームからチャンクを一つずつyieldする
                for chunk in response_stream:
                    if chunk.text:
                        yield chunk.text
                
        except Exception as e:
            logger.error(f"Error during Gemini API stream call: {e}", exc_info=True)
            
            # Log additional context for debugging
            logger.debug(f"Message that caused error: {message}")
            logger.debug(f"History length: {len(history)}")
            
            yield self._get_error_response(e)
    
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
        
        if "api key" in error_str or "authentication" in error_str:
            return (
                "There was an authentication issue with the AI service. "
                "Please check that your API key is correctly configured."
            )
        elif "quota" in error_str or "limit" in error_str:
            return (
                "The AI service is currently experiencing high demand. "
                "Please try again in a few moments."
            )
        elif "network" in error_str or "connection" in error_str:
            return (
                "There was a network connectivity issue. "
                "Please check your internet connection and try again."
            )
        else:
            return (
                "I'm experiencing a temporary issue and cannot process your request right now. "
                "Please try again later."
            )