"""
Chat history service module.

This service encapsulates chat history management logic,
providing a unified interface for both API endpoints and potential future features.
"""

import logging
from typing import Dict, Any, List, Optional

from database.chat_history import ChatHistoryDatabase
from .exceptions import ServiceError, ValidationError

logger = logging.getLogger(__name__)


class ChatHistoryService:
    """Service for chat history management."""

    def __init__(self, db: ChatHistoryDatabase = None):
        """
        Initialize the chat history service.

        Args:
            db: Optional ChatHistoryDatabase instance. If None, creates a new instance.
        """
        self.db = db if db is not None else ChatHistoryDatabase()

    # Session operations

    def create_session(self, name: str) -> Dict[str, Any]:
        """
        Create a new chat session.

        Args:
            name: Session name/title (provided by Pydantic validation with default)

        Returns:
            Dict containing the created session data

        Raises:
            ValidationError: If name is invalid
            ServiceError: For other errors

        Note:
            Validation of name (length, non-empty) is handled by Pydantic model.
            Additional business logic validation can be added here if needed.
        """
        try:
            # Additional business logic validation (if needed)
            # Note: Basic validation (length, non-empty) is already done by Pydantic
            name_stripped = name.strip()

            logger.info(f"Creating new chat session: '{name_stripped}'")
            session = self.db.create_session(name=name_stripped)
            logger.info(f"Successfully created session: {session['id']}")

            return session

        except ValidationError:
            raise
        except Exception as e:
            logger.error(f"Failed to create chat session: {e}", exc_info=True)
            raise ServiceError(f'Failed to create chat session: {str(e)}')

    def get_session(self, session_id: str) -> Optional[Dict[str, Any]]:
        """
        Get a chat session by ID.

        Args:
            session_id: Session ID to retrieve

        Returns:
            Dict containing session data, or None if not found

        Raises:
            ServiceError: If retrieval fails
        """
        try:
            logger.debug(f"Getting chat session: {session_id}")
            return self.db.get_session(session_id)

        except Exception as e:
            logger.error(f"Failed to get chat session: {e}", exc_info=True)
            raise ServiceError(f'Failed to get chat session: {str(e)}')

    def list_sessions(self, limit: int = 100, offset: int = 0) -> Dict[str, Any]:
        """
        List all chat sessions.

        Args:
            limit: Maximum number of sessions to return
            offset: Number of sessions to skip

        Returns:
            Dict with 'sessions' list and 'total_count'

        Raises:
            ServiceError: If listing fails
        """
        try:
            logger.info(f"Listing chat sessions (limit={limit}, offset={offset})")
            sessions, total_count = self.db.list_sessions(limit=limit, offset=offset)

            logger.info(f"Retrieved {len(sessions)} sessions (total: {total_count})")
            return {
                'sessions': sessions,
                'total_count': total_count
            }

        except Exception as e:
            logger.error(f"Failed to list chat sessions: {e}", exc_info=True)
            raise ServiceError(f'Failed to list chat sessions: {str(e)}')

    def update_session(self, session_id: str, name: str) -> Optional[Dict[str, Any]]:
        """
        Update a chat session's name.

        Args:
            session_id: Session ID to update
            name: New session name

        Returns:
            Dict containing updated session data, or None if not found

        Raises:
            ValidationError: If name is invalid
            ServiceError: For other errors
        """
        try:
            # Validate name
            if not name or not name.strip():
                raise ValidationError("Session name cannot be empty")

            if len(name) > 200:
                raise ValidationError("Session name is too long (maximum 200 characters)")

            logger.info(f"Updating chat session: {session_id} - '{name}'")
            session = self.db.update_session(session_id, name.strip())

            if session:
                logger.info(f"Successfully updated session: {session_id}")
            else:
                logger.warning(f"Session not found: {session_id}")

            return session

        except ValidationError:
            raise
        except Exception as e:
            logger.error(f"Failed to update chat session: {e}", exc_info=True)
            raise ServiceError(f'Failed to update chat session: {str(e)}')

    def delete_session(self, session_id: str) -> bool:
        """
        Delete a chat session and all its messages.

        Args:
            session_id: Session ID to delete

        Returns:
            True if session was deleted, False if not found

        Raises:
            ServiceError: If deletion fails
        """
        try:
            logger.info(f"Deleting chat session: {session_id}")
            deleted = self.db.delete_session(session_id)

            if deleted:
                logger.info(f"Successfully deleted session: {session_id}")
            else:
                logger.warning(f"Session not found: {session_id}")

            return deleted

        except Exception as e:
            logger.error(f"Failed to delete chat session: {e}", exc_info=True)
            raise ServiceError(f'Failed to delete chat session: {str(e)}')

    # Message operations

    def add_message(
        self,
        session_id: str,
        role: str,
        content: str
    ) -> Dict[str, Any]:
        """
        Add a message to a chat session.

        Args:
            session_id: Session ID to add message to
            role: Message role ('user' or 'model')
            content: Message content

        Returns:
            Dict containing the created message data

        Raises:
            ValidationError: If parameters are invalid
            ServiceError: For other errors
        """
        try:
            # Validate role
            if role not in ('user', 'model'):
                raise ValidationError(f"Invalid role: {role}. Must be 'user' or 'model'")

            # Validate content
            if not content:
                raise ValidationError("Message content cannot be empty")

            # Check if session exists
            session = self.db.get_session(session_id)
            if session is None:
                raise ValidationError(f"Session not found: {session_id}")

            logger.debug(f"Adding message to session {session_id}: {role}")
            message = self.db.add_message(session_id, role, content)
            logger.debug(f"Successfully added message: {message['id']}")

            return message

        except ValidationError:
            raise
        except Exception as e:
            logger.error(f"Failed to add message: {e}", exc_info=True)
            raise ServiceError(f'Failed to add message: {str(e)}')

    def get_messages(self, session_id: str) -> List[Dict[str, Any]]:
        """
        Get all messages for a chat session.

        Args:
            session_id: Session ID to get messages for

        Returns:
            List of message dicts

        Raises:
            ServiceError: If retrieval fails
        """
        try:
            logger.debug(f"Getting messages for session: {session_id}")
            return self.db.get_messages(session_id)

        except Exception as e:
            logger.error(f"Failed to get messages: {e}", exc_info=True)
            raise ServiceError(f'Failed to get messages: {str(e)}')

    def get_session_with_messages(self, session_id: str) -> Optional[Dict[str, Any]]:
        """
        Get a session with all its messages.

        Args:
            session_id: Session ID to retrieve

        Returns:
            Dict with 'session' and 'messages' keys, or None if session not found

        Raises:
            ServiceError: If retrieval fails
        """
        try:
            logger.debug(f"Getting session with messages: {session_id}")
            return self.db.get_session_with_messages(session_id)

        except Exception as e:
            logger.error(f"Failed to get session with messages: {e}", exc_info=True)
            raise ServiceError(f'Failed to get session with messages: {str(e)}')

    # Utility methods

    def get_stats(self) -> Dict[str, int]:
        """
        Get database statistics.

        Returns:
            Dict with session_count and message_count

        Raises:
            ServiceError: If retrieval fails
        """
        try:
            return self.db.get_stats()

        except Exception as e:
            logger.error(f"Failed to get stats: {e}", exc_info=True)
            raise ServiceError(f'Failed to get stats: {str(e)}')


# Global service instance (lazy-loaded)
_chat_history_service: Optional[ChatHistoryService] = None


def get_chat_history_service() -> ChatHistoryService:
    """
    Get the global ChatHistoryService instance.

    Returns:
        ChatHistoryService: Global service instance
    """
    global _chat_history_service
    if _chat_history_service is None:
        _chat_history_service = ChatHistoryService()
    return _chat_history_service
