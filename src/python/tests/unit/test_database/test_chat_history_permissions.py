"""
Unit tests for ChatHistoryDatabase file-system permission hardening.

All permission assertions are POSIX-only and are skipped on non-POSIX
platforms (e.g. Windows) where os.chmod has different semantics.
"""

import os
import stat

import pytest

from database.chat_history import ChatHistoryDatabase

pytestmark = pytest.mark.skipif(
    os.name != "posix",
    reason="Permission tests require POSIX",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SIDECAR_SUFFIXES = ("-journal", "-wal", "-shm")


def _file_mode(path: str) -> int:
    """Return the permission bits for *path*."""
    return stat.S_IMODE(os.stat(path).st_mode)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestChatHistoryDatabasePermissions:
    """Verify that ChatHistoryDatabase creates files/dirs with restrictive modes."""

    def test_db_file_mode_is_0600(self, tmp_path):
        """
        GIVEN a ChatHistoryDatabase pointed at a temp directory
        WHEN the database is initialised and a row is written
        THEN the database file has mode 0o600 (owner read/write only)
        """
        # Arrange
        db_file = tmp_path / "chat_history.db"

        # Act
        db = ChatHistoryDatabase(db_path=str(db_file))
        session = db.create_session("perm-test-session")
        db.add_message(session["id"], "user", "hello")

        # Assert
        actual_mode = _file_mode(str(db_file))
        assert actual_mode == 0o600, f"DB file mode {oct(actual_mode)} != 0o600"

    def test_db_parent_directory_mode_is_0700(self, tmp_path):
        """
        GIVEN a ChatHistoryDatabase pointed at a temp directory
        WHEN the database is initialised
        THEN the parent directory has mode 0o700 (owner rwx only)

        NOTE: When db_path is passed explicitly the __init__ code does NOT
        call os.makedirs (that branch runs only when db_path is None).
        However it unconditionally calls _secure_permissions on the parent
        directory, so the permission is hardened regardless.
        """
        # Arrange -- parent must already exist when db_path is explicit
        db_dir = tmp_path / "data"
        db_dir.mkdir()
        db_file = db_dir / "chat_history.db"

        # Act
        ChatHistoryDatabase(db_path=str(db_file))

        # Assert
        actual_mode = _file_mode(str(db_dir))
        assert (
            actual_mode == 0o700
        ), f"DB parent directory mode {oct(actual_mode)} != 0o700"

    def test_sidecar_files_are_0600_or_absent(self, tmp_path):
        """
        GIVEN a ChatHistoryDatabase with a written row
        WHEN the transaction completes
        THEN any surviving sidecar files (-journal/-wal/-shm) have mode 0o600.
             If no sidecars remain (DELETE journal mode cleans up on commit)
             the test still passes -- it asserts the .db file itself is 0o600
             and documents that sidecars were cleaned up as expected.
        """
        # Arrange
        db_file = tmp_path / "chat_history.db"
        db = ChatHistoryDatabase(db_path=str(db_file))

        # Act -- create a session and write a message to trigger journal I/O
        session = db.create_session("sidecar-test")
        db.add_message(session["id"], "user", "trigger journal")

        # Collect any surviving sidecar files
        db_path_str = str(db_file)
        sidecars = [
            db_path_str + suffix
            for suffix in _SIDECAR_SUFFIXES
            if os.path.exists(db_path_str + suffix)
        ]

        if sidecars:
            # Assert -- every sidecar must be 0o600
            for sidecar in sidecars:
                mode = _file_mode(sidecar)
                assert mode == 0o600, (
                    f"Sidecar {os.path.basename(sidecar)} mode " f"{oct(mode)} != 0o600"
                )
        else:
            # Sidecars cleaned up by DELETE journal mode (expected).
            # Verify the .db file itself is still locked down.
            db_mode = _file_mode(db_path_str)
            assert db_mode == 0o600, (
                f"DB file mode {oct(db_mode)} != 0o600 "
                "(sidecars absent -- DELETE journal mode cleaned up)"
            )

    def test_permissions_explicit_not_inherited_from_umask(self, tmp_path):
        """
        GIVEN a permissive umask (0o000)
        WHEN ChatHistoryDatabase is initialised and used
        THEN the DB file is 0o600 and the parent dir is 0o700
              (proving permissions are set explicitly, not inherited)
        """
        # Arrange
        db_dir = tmp_path / "umask_test"
        db_dir.mkdir()
        db_file = db_dir / "chat_history.db"
        old_umask = os.umask(0)

        try:
            # Act
            db = ChatHistoryDatabase(db_path=str(db_file))
            session = db.create_session("umask-proof")
            db.add_message(session["id"], "user", "verify umask bypass")

            # Assert -- directory
            dir_mode = _file_mode(str(db_dir))
            assert (
                dir_mode == 0o700
            ), f"Directory mode {oct(dir_mode)} != 0o700 under umask(0)"

            # Assert -- file
            file_mode = _file_mode(str(db_file))
            assert (
                file_mode == 0o600
            ), f"DB file mode {oct(file_mode)} != 0o600 under umask(0)"
        finally:
            os.umask(old_umask)
