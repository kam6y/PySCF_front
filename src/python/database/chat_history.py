"""
Chat history database module.

This module provides SQLite-based persistent storage for AI chat history.
Uses context managers for safe database access and transaction management.
"""

import sqlite3
import logging
import uuid
import os
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from contextlib import contextmanager

logger = logging.getLogger(__name__)

# Restrictive file-system permission modes for sensitive data
_DIR_MODE = 0o700  # Owner-only: rwx------
_FILE_MODE = 0o600  # Owner-only: rw-------
_RESTRICTIVE_UMASK = 0o077  # Complement: files born 0o600, dirs born 0o700


def _secure_permissions(path: Path, mode: int) -> None:
    """Set restrictive file-system permissions on a path (POSIX only).

    Silently skips on non-POSIX platforms (e.g. Windows) where
    ``os.chmod`` has different semantics.

    Args:
        path: File or directory path to secure.
        mode: Octal permission mode (e.g. 0o700 for directories).
    """
    if os.name != "posix":
        return
    try:
        os.chmod(path, mode)
    except OSError as exc:
        logger.warning(
            "Failed to set permissions %s on %s: %s",
            oct(mode),
            path,
            exc,
        )


class ChatHistoryDatabase:
    """
    SQLite database manager for chat history.

    Manages chat sessions and messages with automatic schema initialization
    and migration support.
    """

    def __init__(self, db_path: str = None):
        """
        Initialize the chat history database.

        Args:
            db_path: Path to the SQLite database file.
                    If None, determines path based on environment:
                    1. Environment variable PYSCF_DATA_DIR if set
                    2. Environment variable PYSCF_ENV to detect dev/production
                    3. Project root 'data/' directory (development)
                    4. User home directory '~/.pyscf_app/data/' (packaged app)
        """
        if db_path is None:
            # Priority 1: Check environment variable for explicit data directory
            data_dir_env = os.environ.get("PYSCF_DATA_DIR")
            if data_dir_env:
                data_dir = Path(data_dir_env)
                logger.info(f"Using data directory from PYSCF_DATA_DIR: {data_dir}")
            else:
                # Priority 2: Check PYSCF_ENV for explicit environment detection
                pyscf_env = os.environ.get("PYSCF_ENV", "").lower()

                if pyscf_env == "development":
                    # Explicitly in development mode
                    project_root = Path(__file__).parent.parent.parent.parent
                    data_dir = project_root / "data"
                    logger.info(
                        f"Using project data directory (PYSCF_ENV=development): {data_dir}"
                    )
                elif pyscf_env == "production":
                    # Explicitly in production mode
                    data_dir = Path.home() / ".pyscf_app" / "data"
                    logger.info(
                        f"Using user home data directory (PYSCF_ENV=production): {data_dir}"
                    )
                else:
                    # Priority 3: Auto-detect based on project structure (fallback)
                    try:
                        project_root = Path(__file__).parent.parent.parent.parent
                        data_dir = project_root / "data"

                        # Verify this is a valid development environment by checking multiple indicators
                        is_dev_env = (
                            (project_root / "package.json").exists()
                            and (project_root / "src" / "python" / "app.py").exists()
                            and not (
                                project_root / ".packaged"
                            ).exists()  # Marker file for packaged apps
                        )

                        if is_dev_env:
                            logger.info(
                                f"Using project data directory (auto-detected development): {data_dir}"
                            )
                        else:
                            raise FileNotFoundError("Not in development environment")
                    except (FileNotFoundError, OSError):
                        # Priority 4: Use user home directory (packaged app or fallback)
                        data_dir = Path.home() / ".pyscf_app" / "data"
                        logger.info(
                            f"Using user home data directory (packaged/fallback): {data_dir}"
                        )

            # F3: Create directory with restrictive mode from the start.
            # os.makedirs mode= only sets the leaf directory's perms (parents
            # use umask), which is acceptable — only the leaf holds sensitive
            # files.  The post-create chmod is an idempotent safety net for
            # directories that already existed with lax permissions.
            try:
                os.makedirs(str(data_dir), mode=_DIR_MODE, exist_ok=True)
            except OSError as e:
                logger.error(f"Failed to create data directory {data_dir}: {e}")
                # Fallback to temp directory if we can't create the desired location
                import tempfile

                data_dir = Path(tempfile.gettempdir()) / "pyscf_app_data"
                os.makedirs(str(data_dir), mode=_DIR_MODE, exist_ok=True)
                logger.warning(f"Using temporary directory as fallback: {data_dir}")

            db_path = str(data_dir / "chat_history.db")

        self.db_path = db_path
        logger.info(f"Initializing chat history database at: {self.db_path}")

        # Pre-create the DB file with restrictive permissions (POSIX only)
        # to eliminate a world-readable window between sqlite3.connect()
        # creating the file with default umask and the post-init chmod.
        db_file = Path(self.db_path)
        if os.name == "posix" and not db_file.exists():
            try:
                fd = os.open(
                    str(db_file),
                    os.O_CREAT | os.O_WRONLY,
                    _FILE_MODE,
                )
                os.close(fd)
            except OSError as exc:
                logger.warning(
                    "Failed to pre-create DB file with restrictive "
                    "permissions %s on %s: %s — falling back to default "
                    "creation via sqlite3",
                    oct(_FILE_MODE),
                    db_file,
                    exc,
                )

        self._initialize_database()

        # Harden directory and file permissions (idempotent for existing
        # installs).  Derive paths from self.db_path so this works
        # regardless of which branch above determined the path.
        _secure_permissions(db_file.parent, _DIR_MODE)
        if db_file.exists():
            _secure_permissions(db_file, _FILE_MODE)

    @contextmanager
    def _get_connection(self):
        """
        Context manager for database connections.

        On POSIX, temporarily sets umask to 0o077 so that any SQLite sidecar
        files (e.g. ``<db>-journal``) are created with mode 0o600 instead of
        the process default (typically 0o644).  The previous umask is restored
        in a ``finally`` block.

        NOTE: ``os.umask`` is process-global, not per-thread.  In a
        multi-threaded server the modified umask may briefly affect other
        threads' file creation.  Under concurrent ``_get_connection`` calls,
        interleaving can also leave the process umask parked at the
        restrictive value after both threads restore (thread B captures the
        restrictive mask as its "previous").  The drift is always toward
        *more* restrictive (never wider), so it is security-safe.  The
        0o700 parent directory is the primary protection; the per-file
        umask is defense-in-depth.

        Yields:
            sqlite3.Connection: Database connection with row factory set
        """
        # F2: Restrict umask so SQLite sidecar files are born 0o600.
        _prev_umask: Optional[int] = None
        if os.name == "posix":
            try:
                _prev_umask = os.umask(_RESTRICTIVE_UMASK)
            except OSError as exc:
                logger.warning("Failed to set umask: %s", exc)

        conn: Optional[sqlite3.Connection] = None
        try:
            conn = sqlite3.connect(self.db_path)
            conn.row_factory = sqlite3.Row  # Enable column access by name

            # Enable foreign key constraints (required for CASCADE deletion)
            conn.execute("PRAGMA foreign_keys = ON")

            yield conn
            conn.commit()
        except Exception:
            if conn is not None:
                conn.rollback()
            raise
        finally:
            # Restore the previous umask BEFORE closing the connection so
            # that the process umask is never left restricted if close()
            # raises unexpectedly.  This also covers the case where
            # sqlite3.connect() itself raised before conn was assigned.
            if _prev_umask is not None:
                try:
                    os.umask(_prev_umask)
                except OSError as exc:
                    logger.warning("Failed to restore umask: %s", exc)
            if conn is not None:
                conn.close()

    def _initialize_database(self):
        """
        Initialize database schema if it doesn't exist.

        Creates tables for chat_sessions and chat_messages with appropriate
        indexes for performance optimization.

        IMPORTANT: FOREIGN KEY constraints are enabled per-connection in _get_connection().
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()

            # Create chat_sessions table
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS chat_sessions (
                    id TEXT PRIMARY KEY,
                    name TEXT NOT NULL,
                    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
                )
            """
            )

            # Create chat_messages table
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS chat_messages (
                    id TEXT PRIMARY KEY,
                    session_id TEXT NOT NULL,
                    role TEXT NOT NULL CHECK(role IN ('user', 'model')),
                    content TEXT NOT NULL,
                    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    FOREIGN KEY (session_id) REFERENCES chat_sessions(id) ON DELETE CASCADE
                )
            """
            )

            # Create indexes for performance
            cursor.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_messages_session_id
                ON chat_messages(session_id)
            """
            )

            cursor.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_messages_created_at
                ON chat_messages(created_at)
            """
            )

            cursor.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_sessions_updated_at
                ON chat_sessions(updated_at DESC)
            """
            )

            logger.info(
                "Database schema initialized successfully with FOREIGN KEY constraints enabled"
            )

    # Session operations

    def create_session(self, name: str, session_id: str = None) -> Dict[str, Any]:
        """
        Create a new chat session.

        Args:
            name: Session name/title
            session_id: Optional custom session ID (UUID will be generated if None)

        Returns:
            Dict containing the created session data
        """
        if session_id is None:
            session_id = f"session_{uuid.uuid4()}"

        now = datetime.utcnow().isoformat()

        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                """
                INSERT INTO chat_sessions (id, name, created_at, updated_at)
                VALUES (?, ?, ?, ?)
            """,
                (session_id, name, now, now),
            )

            logger.info(f"Created chat session: {session_id} - '{name}'")

            return {
                "id": session_id,
                "name": name,
                "created_at": now,
                "updated_at": now,
            }

    def get_session(self, session_id: str) -> Optional[Dict[str, Any]]:
        """
        Get a chat session by ID.

        Args:
            session_id: Session ID to retrieve

        Returns:
            Dict containing session data, or None if not found
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                """
                SELECT id, name, created_at, updated_at
                FROM chat_sessions
                WHERE id = ?
            """,
                (session_id,),
            )

            row = cursor.fetchone()
            if row is None:
                return None

            return dict(row)

    def list_sessions(
        self, limit: int = 100, offset: int = 0
    ) -> Tuple[List[Dict[str, Any]], int]:
        """
        List all chat sessions, ordered by most recent update.

        Args:
            limit: Maximum number of sessions to return
            offset: Number of sessions to skip

        Returns:
            Tuple of (list of session summaries with message counts, total count)
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()

            # Get sessions with message counts
            cursor.execute(
                """
                SELECT
                    s.id,
                    s.name,
                    s.created_at,
                    s.updated_at,
                    COUNT(m.id) as message_count
                FROM chat_sessions s
                LEFT JOIN chat_messages m ON s.id = m.session_id
                GROUP BY s.id
                ORDER BY s.updated_at DESC
                LIMIT ? OFFSET ?
            """,
                (limit, offset),
            )

            sessions = [dict(row) for row in cursor.fetchall()]

            # Get total count
            cursor.execute("SELECT COUNT(*) FROM chat_sessions")
            total_count = cursor.fetchone()[0]

            return sessions, total_count

    def update_session(self, session_id: str, name: str) -> Optional[Dict[str, Any]]:
        """
        Update a chat session's name.

        Args:
            session_id: Session ID to update
            name: New session name

        Returns:
            Dict containing updated session data, or None if not found
        """
        now = datetime.utcnow().isoformat()

        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                """
                UPDATE chat_sessions
                SET name = ?, updated_at = ?
                WHERE id = ?
            """,
                (name, now, session_id),
            )

            if cursor.rowcount == 0:
                return None

            logger.info(f"Updated chat session: {session_id} - '{name}'")
            cursor.execute(
                """
                SELECT id, name, created_at, updated_at
                FROM chat_sessions
                WHERE id = ?
            """,
                (session_id,),
            )
            row = cursor.fetchone()
            return dict(row) if row else None

    def delete_session(self, session_id: str) -> bool:
        """
        Delete a chat session and all its messages.

        Args:
            session_id: Session ID to delete

        Returns:
            True if session was deleted, False if not found
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("DELETE FROM chat_sessions WHERE id = ?", (session_id,))

            if cursor.rowcount == 0:
                return False

            logger.info(f"Deleted chat session: {session_id}")
            return True

    # Message operations

    def add_message(
        self, session_id: str, role: str, content: str, message_id: str = None
    ) -> Dict[str, Any]:
        """
        Add a message to a chat session.

        Args:
            session_id: Session ID to add message to
            role: Message role ('user' or 'model')
            content: Message content
            message_id: Optional custom message ID (UUID will be generated if None)

        Returns:
            Dict containing the created message data

        Raises:
            ValueError: If role is not 'user' or 'model'
        """
        if role not in ("user", "model"):
            raise ValueError(f"Invalid role: {role}. Must be 'user' or 'model'")

        if message_id is None:
            message_id = f"msg_{uuid.uuid4()}"

        now = datetime.utcnow().isoformat()

        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                """
                INSERT INTO chat_messages (id, session_id, role, content, created_at)
                VALUES (?, ?, ?, ?, ?)
            """,
                (message_id, session_id, role, content, now),
            )

            # Update session's updated_at timestamp (within same connection)
            cursor.execute(
                """
                UPDATE chat_sessions
                SET updated_at = ?
                WHERE id = ?
            """,
                (now, session_id),
            )

            logger.debug(
                f"Added message to session {session_id}: {role} - {len(content)} chars"
            )

            return {
                "id": message_id,
                "session_id": session_id,
                "role": role,
                "content": content,
                "created_at": now,
            }

    def get_messages(self, session_id: str) -> List[Dict[str, Any]]:
        """
        Get all messages for a chat session, ordered by creation time.

        Args:
            session_id: Session ID to get messages for

        Returns:
            List of message dicts
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                """
                SELECT id, session_id, role, content, created_at
                FROM chat_messages
                WHERE session_id = ?
                ORDER BY created_at ASC
            """,
                (session_id,),
            )

            return [dict(row) for row in cursor.fetchall()]

    def get_session_with_messages(self, session_id: str) -> Optional[Dict[str, Any]]:
        """
        Get a session with all its messages.

        Args:
            session_id: Session ID to retrieve

        Returns:
            Dict with 'session' and 'messages' keys, or None if session not found
        """
        session = self.get_session(session_id)
        if session is None:
            return None

        messages = self.get_messages(session_id)

        return {"session": session, "messages": messages}
