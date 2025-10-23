"""
Supervisor Tool Module

This module provides tools for the Supervisor agent to search and identify
quantum chemistry calculations across the application.

The Supervisor uses these tools to:
- List all available calculations
- Search for calculations by name, status, method, basis set, or date
- Provide specific calculation IDs to worker agents (especially science_analyst)
"""

import json
import logging
from typing import Optional
from datetime import datetime

from agent.quantum_calculator.tools import list_all_calculations as _base_list_all_calculations

# Logger configuration
logger = logging.getLogger(__name__)


def list_all_calculations() -> str:
    """
    Retrieves a list of all quantum chemistry calculations that have been executed and saved.

    This is a wrapper around the base implementation in quantum_calculator.tools.
    Returns a list of calculation summaries including ID, name, status, and creation date.

    Use this when the user asks "show me all calculations" or "what calculations are available".

    Returns:
        str: JSON string of calculation list. On success, contains an array of calculation summaries.
             On error, contains structured text with error message.
    """
    # Delegate to the base implementation
    return _base_list_all_calculations()


def find_calculations(
    name_query: Optional[str] = None,
    status: Optional[str] = None,
    calculation_method: Optional[str] = None,
    basis_function: Optional[str] = None,
    date_from: Optional[str] = None,
    date_to: Optional[str] = None
) -> str:
    """
    Search for calculations using flexible filtering criteria.

    This function searches through all saved calculations and filters them based on
    the provided criteria. All parameters are optional - you can combine them as needed.
    Performs case-insensitive partial matching for name_query.

    Args:
        name_query (str, optional): Partial match search in calculation name.
                                   Case-insensitive. Example: "水" matches "水分子のDFT計算"
        status (str, optional): Filter by status. Options: "completed", "running", "error"
        calculation_method (str, optional): Filter by calculation method.
                                           Options: "DFT", "HF", "MP2", "CCSD", "TDDFT", "CASCI", "CASSCF"
        basis_function (str, optional): Filter by basis set. Case-insensitive.
                                       Examples: "6-31G(d)", "cc-pVDZ", "def2-SVP"
        date_from (str, optional): Start date for date range filtering (ISO format: YYYY-MM-DD)
        date_to (str, optional): End date for date range filtering (ISO format: YYYY-MM-DD)

    Returns:
        str: JSON string containing search results.

             Single match:
             {
               "success": true,
               "count": 1,
               "calculation": {...}  // Detailed calculation info
             }

             Multiple matches:
             {
               "success": true,
               "multiple_matches": true,
               "count": N,
               "message": "N件の計算が見つかりました。いずれかを選択してください:",
               "calculations": [
                 {
                   "id": "calc_...",
                   "name": "...",
                   "status": "...",
                   "calculation_method": "...",
                   "basis_function": "...",
                   "created_at": "..."
                 },
                 ...
               ]
             }

             No matches:
             {
               "success": false,
               "error": "検索条件に一致する計算が見つかりませんでした。"
             }
    """
    try:
        # Get all calculations
        logger.debug("Fetching all calculations for filtering")
        all_calculations_json = _base_list_all_calculations()
        all_calculations_result = json.loads(all_calculations_json)

        if not all_calculations_result.get('success'):
            return all_calculations_json  # Return error from base function

        # Extract calculations list from the nested data structure
        # Structure: {"success": true, "data": {"calculations": [...], "count": N, "base_directory": "..."}}
        data = all_calculations_result.get('data', {})
        calculations = data.get('calculations', [])
        logger.info(f"Retrieved {len(calculations)} total calculations for filtering")

        # Apply filters
        filtered = calculations

        # Filter by name (case-insensitive partial match)
        if name_query:
            name_lower = name_query.lower()
            filtered = [
                calc for calc in filtered
                if name_lower in calc.get('name', '').lower()
            ]
            logger.debug(f"After name filter '{name_query}': {len(filtered)} calculations")

        # Filter by status (exact match)
        if status:
            filtered = [
                calc for calc in filtered
                if calc.get('status', '').lower() == status.lower()
            ]
            logger.debug(f"After status filter '{status}': {len(filtered)} calculations")

        # Filter by calculation_method (case-insensitive exact match)
        if calculation_method:
            filtered = [
                calc for calc in filtered
                if calc.get('calculation_method', '').lower() == calculation_method.lower()
            ]
            logger.debug(f"After method filter '{calculation_method}': {len(filtered)} calculations")

        # Filter by basis_function (case-insensitive exact match)
        if basis_function:
            filtered = [
                calc for calc in filtered
                if calc.get('basis_function', '').lower() == basis_function.lower()
            ]
            logger.debug(f"After basis filter '{basis_function}': {len(filtered)} calculations")

        # Filter by date range
        if date_from or date_to:
            def parse_date(date_str: str) -> Optional[datetime]:
                """Parse ISO date string to datetime object."""
                try:
                    return datetime.fromisoformat(date_str.replace('Z', '+00:00'))
                except (ValueError, AttributeError):
                    return None

            date_from_dt = parse_date(date_from) if date_from else None
            date_to_dt = parse_date(date_to) if date_to else None

            filtered_by_date = []
            for calc in filtered:
                created_at = calc.get('created_at')
                if not created_at:
                    continue

                calc_date = parse_date(created_at)
                if not calc_date:
                    continue

                # Check date range
                if date_from_dt and calc_date < date_from_dt:
                    continue
                if date_to_dt and calc_date > date_to_dt:
                    continue

                filtered_by_date.append(calc)

            filtered = filtered_by_date
            logger.debug(f"After date filter (from={date_from}, to={date_to}): {len(filtered)} calculations")

        # Process results
        count = len(filtered)

        if count == 0:
            # No matches
            logger.info("No calculations matched the search criteria")
            return json.dumps({
                "success": False,
                "error": "検索条件に一致する計算が見つかりませんでした。条件を変更して再度お試しください。"
            }, ensure_ascii=False, indent=2)

        elif count == 1:
            # Single match - return detailed info
            logger.info(f"Found single matching calculation: {filtered[0].get('id')}")
            return json.dumps({
                "success": True,
                "count": 1,
                "calculation": filtered[0]
            }, ensure_ascii=False, indent=2)

        else:
            # Multiple matches - return list for user selection
            logger.info(f"Found {count} matching calculations")

            # Sort by created_at (newest first)
            sorted_calcs = sorted(
                filtered,
                key=lambda x: x.get('created_at', ''),
                reverse=True
            )

            # Prepare simplified list for presentation
            calc_list = []
            for calc in sorted_calcs:
                calc_list.append({
                    "id": calc.get('id'),
                    "name": calc.get('name'),
                    "status": calc.get('status'),
                    "calculation_method": calc.get('calculation_method'),
                    "basis_function": calc.get('basis_function'),
                    "created_at": calc.get('created_at')
                })

            return json.dumps({
                "success": True,
                "multiple_matches": True,
                "count": count,
                "message": f"{count}件の計算が見つかりました。いずれかを選択してください:",
                "calculations": calc_list
            }, ensure_ascii=False, indent=2)

    except json.JSONDecodeError as e:
        logger.error(f"JSON parsing error in find_calculations: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"計算データの解析中にエラーが発生しました: {str(e)}"
        }, ensure_ascii=False, indent=2)

    except Exception as e:
        logger.error(f"Unexpected error in find_calculations: {e}", exc_info=True)
        return json.dumps({
            "success": False,
            "error": f"計算検索中に予期しないエラーが発生しました: {str(e)}"
        }, ensure_ascii=False, indent=2)


__all__ = [
    'list_all_calculations',
    'find_calculations',
]
