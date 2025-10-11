"""
Research Tools Module

This module defines LangChain tools for academic research, primarily focused on
searching and retrieving papers from arXiv.

Tools:
- search_arxiv: Search arXiv database for academic papers
"""

import logging
import arxiv
from langchain_core.tools import tool

# Set up logging
logger = logging.getLogger(__name__)


@tool
def search_arxiv(query: str, max_results: int = 5) -> str:
    """
    Search the arXiv database for papers matching the query.

    Args:
        query: Search query string (e.g., "density functional theory")
        max_results: Maximum number of results to return (default: 5)

    Returns:
        Formatted string containing paper titles, authors, summaries, and PDF links.
        Returns an error message if the search fails.

    Example:
        >>> result = search_arxiv("quantum chemistry CASSCF", max_results=3)
    """
    try:
        logger.info(f"Searching arXiv with query: '{query}', max_results: {max_results}")

        search = arxiv.Search(
            query=query,
            max_results=max_results,
            sort_by=arxiv.SortCriterion.Relevance
        )
        results = arxiv.Client().results(search)

        formatted_results = []
        paper_count = 0
        for r in results:
            paper_count += 1
            # Log PDF URL to verify it's being retrieved correctly
            logger.debug(f"Paper {paper_count}: {r.title[:50]}... | PDF URL: {r.pdf_url}")

            formatted_results.append(
                f"**Title:** {r.title}\n"
                f"**Authors:** {', '.join(author.name for author in r.authors)}\n"
                f"**Published:** {r.published.date()}\n"
                f"**Summary:** {r.summary.replace(chr(10), ' ')}\n"
                f"**PDF Link:** {r.pdf_url}"
            )

        if not formatted_results:
            logger.warning(f"No papers found for query: '{query}'")
            return "No relevant papers found on arXiv for the given query."

        logger.info(f"Successfully retrieved {len(formatted_results)} papers from arXiv")
        return "\n\n---\n\n".join(formatted_results)
    except Exception as e:
        logger.error(f"Error searching arXiv with query '{query}': {str(e)}", exc_info=True)
        return f"An error occurred while searching arXiv: {str(e)}"
