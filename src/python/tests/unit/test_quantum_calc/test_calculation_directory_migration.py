from pathlib import Path

import pytest
import quantum_calc._calculation_directory_migration as migration_module
from quantum_calc import CalculationDirectoryMigration


def test_move_calculations_directory_moves_all_calculations_when_no_conflicts(
    tmp_path: Path,
) -> None:
    """
    GIVEN source calculations and an empty destination
    WHEN moving the calculations directory
    THEN all calculations move and the migration succeeds
    """
    old_dir = tmp_path / "old_calculations"
    new_dir = tmp_path / "new_calculations"
    old_calc_a = old_dir / "calc_a"
    old_calc_b = old_dir / "calc_b"

    old_calc_a.mkdir(parents=True)
    old_calc_b.mkdir()
    (old_calc_a / "metadata.json").write_text("old a", encoding="utf-8")
    (old_calc_b / "metadata.json").write_text("old b", encoding="utf-8")

    migration = CalculationDirectoryMigration(base_dir=str(old_dir))

    result = migration.move_calculations_directory(str(new_dir))

    assert result["success"] is True
    assert result["moved_count"] == 2
    assert result["failed_count"] == 0
    calc_a_metadata = new_dir / "calc_a" / "metadata.json"
    calc_b_metadata = new_dir / "calc_b" / "metadata.json"
    assert calc_a_metadata.read_text(encoding="utf-8") == "old a"
    assert calc_b_metadata.read_text(encoding="utf-8") == "old b"
    assert migration.base_dir == new_dir.resolve()


def test_move_calculations_directory_with_no_calculations_switches_base_directory(
    tmp_path: Path,
) -> None:
    """
    GIVEN the current calculations directory is empty
    WHEN moving the calculations directory
    THEN no files move and the base directory still switches
    """
    old_dir = tmp_path / "old_calculations"
    new_dir = tmp_path / "new_calculations"
    old_dir.mkdir()

    migration = CalculationDirectoryMigration(base_dir=str(old_dir))

    result = migration.move_calculations_directory(str(new_dir))

    assert result["success"] is True
    assert result["moved_count"] == 0
    assert result["message"] == "No calculations to move"
    assert migration.base_dir == new_dir.resolve()


def test_move_calculations_directory_to_pyscf_subfolder_excludes_destination(
    tmp_path: Path,
) -> None:
    """
    GIVEN the destination is the allowed PySCF_calculations child folder
    WHEN moving the calculations directory
    THEN the newly created destination folder is not treated as a calculation
    """
    old_dir = tmp_path / "project_calculations"
    new_dir = old_dir / "PySCF_calculations"
    old_calc = old_dir / "calc_a"

    old_calc.mkdir(parents=True)
    (old_calc / "metadata.json").write_text("old a", encoding="utf-8")

    migration = CalculationDirectoryMigration(base_dir=str(old_dir))

    result = migration.move_calculations_directory(str(new_dir))

    assert result["success"] is True
    assert result["moved_count"] == 1
    assert (new_dir / "calc_a" / "metadata.json").read_text(encoding="utf-8") == "old a"
    assert migration.base_dir == new_dir.resolve()


def test_move_calculations_directory_with_destination_conflict_does_not_move_anything(
    tmp_path: Path,
) -> None:
    """
    GIVEN the destination already contains a calculation directory with the same name
    WHEN moving the calculations directory
    THEN the migration fails before moving any source calculations
    """
    old_dir = tmp_path / "old_calculations"
    new_dir = tmp_path / "new_calculations"
    old_calc_a = old_dir / "calc_a"
    old_calc_b = old_dir / "calc_b"
    existing_new_calc_b = new_dir / "calc_b"

    old_calc_a.mkdir(parents=True)
    old_calc_b.mkdir()
    existing_new_calc_b.mkdir(parents=True)
    (old_calc_a / "metadata.json").write_text("old a", encoding="utf-8")
    (old_calc_b / "metadata.json").write_text("old b", encoding="utf-8")
    (existing_new_calc_b / "metadata.json").write_text("new b", encoding="utf-8")

    migration = CalculationDirectoryMigration(base_dir=str(old_dir))

    result = migration.move_calculations_directory(str(new_dir))

    assert result["success"] is False
    assert result["moved_count"] == 0
    assert result["failed_count"] == 1
    assert result["failed_moves"] == [
        {"name": "calc_b", "reason": "Destination already exists"}
    ]
    assert (old_calc_a / "metadata.json").read_text(encoding="utf-8") == "old a"
    assert (old_calc_b / "metadata.json").read_text(encoding="utf-8") == "old b"
    assert (existing_new_calc_b / "metadata.json").read_text(encoding="utf-8") == "new b"


def test_move_calculations_directory_rolls_back_successful_moves_after_io_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    GIVEN one calculation move fails after another calculation has already moved
    WHEN moving the calculations directory
    THEN the successful move is rolled back to the original directory
    """
    old_dir = tmp_path / "old_calculations"
    new_dir = tmp_path / "new_calculations"
    old_calc_a = old_dir / "calc_a"
    old_calc_b = old_dir / "calc_b"

    old_calc_a.mkdir(parents=True)
    old_calc_b.mkdir()
    (old_calc_a / "metadata.json").write_text("old a", encoding="utf-8")
    (old_calc_b / "metadata.json").write_text("old b", encoding="utf-8")

    real_move = migration_module.shutil.move

    def fail_calc_b_move(src: str, dest: str) -> str:
        if src.endswith("calc_b"):
            raise OSError("disk full")
        return real_move(src, dest)

    monkeypatch.setattr(migration_module.shutil, "move", fail_calc_b_move)
    migration = CalculationDirectoryMigration(base_dir=str(old_dir))

    result = migration.move_calculations_directory(str(new_dir))

    assert result["success"] is False
    assert result["moved_count"] == 0
    assert result["failed_count"] == 1
    assert result["failed_moves"] == [{"name": "calc_b", "reason": "disk full"}]
    assert (old_calc_a / "metadata.json").read_text(encoding="utf-8") == "old a"
    assert (old_calc_b / "metadata.json").read_text(encoding="utf-8") == "old b"
    assert not (new_dir / "calc_a").exists()
