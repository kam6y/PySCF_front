"""Unit tests for CubeArtifactService — security field-removal regression."""

from quantum_calc._cube_artifact_service import CubeArtifactService


def test_get_cube_files_info_does_not_expose_file_path(tmp_path):
    """
    GIVEN a calculation directory containing a valid CUBE file
    WHEN get_cube_files_info() is called
    THEN each returned item contains 'filename' but NOT 'file_path'
         (absolute-path leak removed)
    """
    # Arrange: create orbital subdirectory with a properly-named cube file
    orbital_dir = tmp_path / "orbital"
    orbital_dir.mkdir()
    cube_file = orbital_dir / "orbital_0_grid40.cube"
    cube_file.write_text("dummy cube content", encoding="utf-8")

    service = CubeArtifactService(base_dir=str(tmp_path))

    # Act
    cube_files = service.get_cube_files_info(str(tmp_path))

    # Assert
    assert len(cube_files) == 1
    item = cube_files[0]
    assert item["filename"] == "orbital_0_grid40.cube"
    assert item["orbital_index"] == 0
    assert item["grid_size"] == 40
    assert "file_path" not in item
