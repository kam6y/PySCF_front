"""
Unit tests for Research Agent tools.

Tests the search_arxiv tool with mocked arxiv library to verify
search functionality, error handling, and output formatting.
"""

import pytest
from unittest.mock import MagicMock, patch
from datetime import datetime

from agent.literature_surveyor.tools import search_arxiv


# ============================================================================
# Mock Data Classes
# ============================================================================

class MockArxivAuthor:
    """Mock arxiv.Result.Author object."""
    def __init__(self, name: str):
        self.name = name


class MockArxivResult:
    """Mock arxiv.Result object."""
    def __init__(
        self,
        title: str,
        authors: list,
        published: datetime,
        summary: str,
        pdf_url: str
    ):
        self.title = title
        self.authors = [MockArxivAuthor(name) for name in authors]
        self.published = published
        self.summary = summary
        self.pdf_url = pdf_url


# ============================================================================
# search_arxiv() Success Tests
# ============================================================================

@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_success_single_result(mock_search_class, mock_client_class):
    """
    GIVEN arxiv API returns a single paper result
    WHEN search_arxiv is called with valid query
    THEN it should return formatted Markdown string with paper details
    """
    # ARRANGE
    mock_result = MockArxivResult(
        title="Quantum Chemistry with PySCF",
        authors=["John Doe", "Jane Smith"],
        published=datetime(2024, 1, 15),
        summary="This paper discusses quantum chemistry methods.",
        pdf_url="https://arxiv.org/pdf/2401.12345.pdf"
    )

    mock_client = MagicMock()
    mock_client.results.return_value = [mock_result]
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    result = search_arxiv.func("quantum chemistry", max_results=1)

    # ASSERT
    assert isinstance(result, str)
    assert "**Title:** Quantum Chemistry with PySCF" in result
    assert "**Authors:** John Doe, Jane Smith" in result
    assert "**Published:** 2024-01-15" in result
    assert "**Summary:** This paper discusses quantum chemistry methods." in result
    assert "**PDF Link:** https://arxiv.org/pdf/2401.12345.pdf" in result


@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_success_multiple_results(mock_search_class, mock_client_class):
    """
    GIVEN arxiv API returns multiple paper results
    WHEN search_arxiv is called with max_results=3
    THEN it should return all results separated by Markdown dividers
    """
    # ARRANGE
    mock_results = [
        MockArxivResult(
            title="Paper 1",
            authors=["Author A"],
            published=datetime(2024, 1, 1),
            summary="Summary 1",
            pdf_url="https://arxiv.org/pdf/1.pdf"
        ),
        MockArxivResult(
            title="Paper 2",
            authors=["Author B", "Author C"],
            published=datetime(2024, 2, 1),
            summary="Summary 2",
            pdf_url="https://arxiv.org/pdf/2.pdf"
        ),
        MockArxivResult(
            title="Paper 3",
            authors=["Author D"],
            published=datetime(2024, 3, 1),
            summary="Summary 3",
            pdf_url="https://arxiv.org/pdf/3.pdf"
        )
    ]

    mock_client = MagicMock()
    mock_client.results.return_value = mock_results
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    result = search_arxiv.func("test query", max_results=3)

    # ASSERT
    assert "Paper 1" in result
    assert "Paper 2" in result
    assert "Paper 3" in result
    assert result.count("---") == 2  # Two dividers for 3 results
    assert "Author A" in result
    assert "Author B, Author C" in result


@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_handles_newlines_in_summary(mock_search_class, mock_client_class):
    """
    GIVEN arxiv result contains newlines in summary
    WHEN search_arxiv is called
    THEN it should replace newlines with spaces for readability
    """
    # ARRANGE
    mock_result = MockArxivResult(
        title="Test Paper",
        authors=["Test Author"],
        published=datetime(2024, 1, 1),
        summary="This is a summary\nwith multiple\nlines in it.",
        pdf_url="https://arxiv.org/pdf/test.pdf"
    )

    mock_client = MagicMock()
    mock_client.results.return_value = [mock_result]
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    result = search_arxiv.func("test")

    # ASSERT
    assert "This is a summary with multiple lines in it." in result
    # Summary field should have newlines replaced with spaces (not the whole result)
    summary_part = result.split("**Summary:**")[1].split("**PDF Link:**")[0].strip()
    assert "This is a summary with multiple lines in it." in summary_part


# ============================================================================
# search_arxiv() Empty Results Tests
# ============================================================================

@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_no_results(mock_search_class, mock_client_class):
    """
    GIVEN arxiv API returns empty results
    WHEN search_arxiv is called
    THEN it should return appropriate message
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client.results.return_value = []  # Empty results
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    result = search_arxiv.func("nonexistent topic xyzabc123")

    # ASSERT
    assert "No relevant papers found" in result
    assert "arxiv" in result.lower()


# ============================================================================
# search_arxiv() Error Handling Tests
# ============================================================================

@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_handles_api_error(mock_search_class, mock_client_class):
    """
    GIVEN arxiv API raises an exception
    WHEN search_arxiv is called
    THEN it should return error message without crashing
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client.results.side_effect = Exception("ArXiv API connection error")
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    result = search_arxiv.func("test query")

    # ASSERT
    assert "error occurred while searching arXiv" in result
    assert "ArXiv API connection error" in result


@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_handles_search_creation_error(mock_search_class, mock_client_class):
    """
    GIVEN arxiv.Search raises an exception during initialization
    WHEN search_arxiv is called
    THEN it should return error message
    """
    # ARRANGE
    mock_search_class.side_effect = ValueError("Invalid query parameter")

    # ACT
    result = search_arxiv.func("test query")

    # ASSERT
    assert "error occurred while searching arXiv" in result
    assert "Invalid query parameter" in result


# ============================================================================
# search_arxiv() Parameter Tests
# ============================================================================

@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_respects_max_results_parameter(mock_search_class, mock_client_class):
    """
    GIVEN max_results parameter is specified
    WHEN search_arxiv is called
    THEN it should pass correct max_results to arxiv.Search
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client.results.return_value = []
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    search_arxiv.func("test", max_results=10)

    # ASSERT
    mock_search_class.assert_called_once()
    call_kwargs = mock_search_class.call_args.kwargs
    assert call_kwargs['max_results'] == 10


@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_uses_default_max_results(mock_search_class, mock_client_class):
    """
    GIVEN max_results parameter is not specified
    WHEN search_arxiv is called
    THEN it should use default value of 5
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client.results.return_value = []
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    search_arxiv("test")

    # ASSERT
    mock_search_class.assert_called_once()
    call_kwargs = mock_search_class.call_args.kwargs
    assert call_kwargs['max_results'] == 5


@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_uses_relevance_sorting(mock_search_class, mock_client_class):
    """
    GIVEN search_arxiv is called
    WHEN arxiv.Search is created
    THEN it should use Relevance as sort criterion
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client.results.return_value = []
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    search_arxiv("test")

    # ASSERT
    from agent.literature_surveyor.tools import arxiv
    mock_search_class.assert_called_once()
    call_kwargs = mock_search_class.call_args.kwargs
    assert call_kwargs['sort_by'] == arxiv.SortCriterion.Relevance


# ============================================================================
# LangChain Tool Integration Tests
# ============================================================================

def test_search_arxiv_is_langchain_tool():
    """
    GIVEN search_arxiv function
    WHEN inspecting its attributes
    THEN it should have LangChain tool metadata
    """
    # ASSERT
    # LangChain @tool decorator adds specific attributes
    assert hasattr(search_arxiv, 'name')
    assert hasattr(search_arxiv, 'description')
    assert search_arxiv.name == 'search_arxiv'


def test_search_arxiv_tool_has_proper_docstring():
    """
    GIVEN search_arxiv function
    WHEN checking its description
    THEN it should have informative documentation
    """
    # ASSERT
    description = search_arxiv.description
    assert description is not None
    assert len(description) > 0
    assert 'arXiv' in description or 'arxiv' in description


# ============================================================================
# Output Format Tests
# ============================================================================

@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_output_is_markdown_formatted(mock_search_class, mock_client_class):
    """
    GIVEN arxiv API returns results
    WHEN search_arxiv is called
    THEN output should be properly formatted with Markdown
    """
    # ARRANGE
    mock_result = MockArxivResult(
        title="Test",
        authors=["Author"],
        published=datetime(2024, 1, 1),
        summary="Summary",
        pdf_url="https://arxiv.org/pdf/test.pdf"
    )

    mock_client = MagicMock()
    mock_client.results.return_value = [mock_result]
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    result = search_arxiv.func("test")

    # ASSERT
    # Check for Markdown bold formatting
    assert "**Title:**" in result
    assert "**Authors:**" in result
    assert "**Published:**" in result
    assert "**Summary:**" in result
    assert "**PDF Link:**" in result


@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_multiple_authors_formatted_correctly(mock_search_class, mock_client_class):
    """
    GIVEN arxiv result has multiple authors
    WHEN search_arxiv is called
    THEN authors should be comma-separated
    """
    # ARRANGE
    mock_result = MockArxivResult(
        title="Test",
        authors=["First Author", "Second Author", "Third Author"],
        published=datetime(2024, 1, 1),
        summary="Summary",
        pdf_url="https://arxiv.org/pdf/test.pdf"
    )

    mock_client = MagicMock()
    mock_client.results.return_value = [mock_result]
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    result = search_arxiv.func("test")

    # ASSERT
    assert "First Author, Second Author, Third Author" in result


# ============================================================================
# Edge Cases Tests
# ============================================================================

@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_handles_empty_query(mock_search_class, mock_client_class):
    """
    GIVEN empty query string
    WHEN search_arxiv is called
    THEN it should handle gracefully (arxiv library may raise error)
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client.results.return_value = []
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT - Should not raise exception
    result = search_arxiv.func("")

    # ASSERT
    assert isinstance(result, str)


@patch('agent.research.tools.arxiv.Client')
@patch('agent.research.tools.arxiv.Search')
def test_search_arxiv_handles_special_characters(mock_search_class, mock_client_class):
    """
    GIVEN query with special characters
    WHEN search_arxiv is called
    THEN it should handle without errors
    """
    # ARRANGE
    mock_client = MagicMock()
    mock_client.results.return_value = []
    mock_client_class.return_value = mock_client

    mock_search = MagicMock()
    mock_search_class.return_value = mock_search

    # ACT
    result = search_arxiv.func("quantum chemistry: DFT & TDDFT (test)")

    # ASSERT
    assert isinstance(result, str)
    mock_search_class.assert_called_once()
