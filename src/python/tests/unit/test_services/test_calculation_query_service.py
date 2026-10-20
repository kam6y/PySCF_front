"""Unit tests for CalculationQueryService — security field-removal regression."""

from unittest.mock import Mock

from services.calculation_query_service import CalculationQueryService


def test_list_calculations_does_not_expose_base_directory():
    """
    GIVEN a CalculationQueryService with one calculation available
    WHEN list_calculations() is called
    THEN the returned dict contains 'calculations' and 'count'
         but does NOT contain 'base_directory' (absolute-path leak removed)
    """
    context = Mock()
    context.repository.list_calculations.return_value = [
        {"id": "calc-001", "name": "Water", "status": "completed"},
    ]

    service = CalculationQueryService(context)
    result = service.list_calculations()

    assert "calculations" in result
    assert result["count"] == 1
    assert len(result["calculations"]) == 1
    assert "base_directory" not in result
