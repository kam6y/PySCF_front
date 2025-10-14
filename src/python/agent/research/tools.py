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


@tool
def search_tavily(query: str, max_results: int = 5) -> str:
    """
    Search the web using Tavily API for up-to-date information.

    Args:
        query: Search query string (e.g., "latest advances in quantum computing")
        max_results: Maximum number of results to return (default: 5)

    Returns:
        Formatted string containing web search results with titles, URLs, and content.
        Returns an error message if the search fails.

    Example:
        >>> result = search_tavily("TDDFT applications 2024", max_results=3)
    """
    try:
        import os
        from langchain_tavily import TavilySearch
        from quantum_calc.settings_manager import get_current_settings

        logger.info(f"Searching Tavily with query: '{query}', max_results: {max_results}")

        # Get Tavily API key from settings or environment
        api_key = None
        try:
            settings = get_current_settings()
            api_key = settings.tavily_api_key
        except Exception as e:
            logger.debug(f"Could not load settings for Tavily API key: {e}")

        # Fallback to environment variable
        if not api_key:
            api_key = os.environ.get("TAVILY_API_KEY")

        if not api_key:
            logger.warning("Tavily API key not configured")
            return "Tavily search is not available. Please configure TAVILY_API_KEY in Settings."

        # Initialize Tavily search tool with API key
        tavily_tool = TavilySearch(tavily_api_key=api_key, max_results=max_results)

        # Perform search
        search_results = tavily_tool.invoke(query)
        
        # Format results
        formatted_results = []
        if "results" in search_results and search_results["results"]:
            for idx, result in enumerate(search_results["results"], 1):
                title = result.get("title", "No title")
                url = result.get("url", "")
                content = result.get("content", "No content available")
                
                logger.debug(f"Result {idx}: {title[:50]}... | URL: {url}")
                
                formatted_results.append(
                    f"**{idx}. {title}**\n"
                    f"**URL:** {url}\n"
                    f"**Content:** {content}"
                )
        
        if not formatted_results:
            logger.warning(f"No results found for query: '{query}'")
            return "No relevant web results found for the given query."
        
        logger.info(f"Successfully retrieved {len(formatted_results)} results from Tavily")
        return "\n\n---\n\n".join(formatted_results)
        
    except ImportError as e:
        logger.error(f"Tavily package not installed: {str(e)}")
        return "Tavily search is not available. Please install langchain-tavily package."
    except Exception as e:
        logger.error(f"Error searching Tavily with query '{query}': {str(e)}", exc_info=True)
        return f"An error occurred while searching Tavily: {str(e)}"


@tool
def search_pubmed(query: str, max_results: int = 5) -> str:
    """
    Search PubMed database for biomedical and life sciences papers.

    Args:
        query: Search query string (e.g., "protein folding molecular dynamics")
        max_results: Maximum number of results to return (default: 5)

    Returns:
        Formatted string containing paper titles, authors, abstracts, and PubMed links.
        Returns an error message if the search fails.

    Example:
        >>> result = search_pubmed("CRISPR gene editing", max_results=3)
    """
    try:
        from Bio import Entrez
        
        logger.info(f"Searching PubMed with query: '{query}', max_results: {max_results}")
        
        # Set email for Entrez (required by NCBI)
        Entrez.email = "pyscf-research-agent@example.com"
        
        # Search PubMed
        search_handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            sort="relevance"
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = search_results.get("IdList", [])
        
        if not id_list:
            logger.warning(f"No papers found for query: '{query}'")
            return "No relevant papers found on PubMed for the given query."
        
        # Fetch details for each paper
        fetch_handle = Entrez.efetch(
            db="pubmed",
            id=id_list,
            rettype="abstract",
            retmode="xml"
        )
        papers = Entrez.read(fetch_handle)
        fetch_handle.close()
        
        formatted_results = []
        
        for idx, paper in enumerate(papers["PubmedArticle"], 1):
            try:
                article = paper["MedlineCitation"]["Article"]
                pmid = paper["MedlineCitation"]["PMID"]
                
                # Extract title
                title = article.get("ArticleTitle", "No title")
                
                # Extract authors
                authors_list = article.get("AuthorList", [])
                authors = ", ".join([
                    f"{author.get('LastName', '')} {author.get('Initials', '')}".strip()
                    for author in authors_list[:5]  # First 5 authors
                ])
                if len(authors_list) > 5:
                    authors += " et al."
                
                # Extract abstract
                abstract_sections = article.get("Abstract", {}).get("AbstractText", [])
                if abstract_sections:
                    if isinstance(abstract_sections, list):
                        abstract = " ".join([str(section) for section in abstract_sections])
                    else:
                        abstract = str(abstract_sections)
                else:
                    abstract = "No abstract available"
                
                # Extract publication year
                pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
                year = pub_date.get("Year", "Unknown year")
                
                # PubMed URL
                pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                
                logger.debug(f"Paper {idx}: {title[:50]}... | PMID: {pmid}")
                
                formatted_results.append(
                    f"**Title:** {title}\n"
                    f"**Authors:** {authors}\n"
                    f"**Year:** {year}\n"
                    f"**Abstract:** {abstract[:500]}{'...' if len(abstract) > 500 else ''}\n"
                    f"**PubMed Link:** {pubmed_url}"
                )
            except Exception as paper_error:
                logger.warning(f"Error parsing paper {idx}: {str(paper_error)}")
                continue
        
        if not formatted_results:
            logger.warning(f"Could not parse any papers for query: '{query}'")
            return "Found papers but could not parse their details."
        
        logger.info(f"Successfully retrieved {len(formatted_results)} papers from PubMed")
        return "\n\n---\n\n".join(formatted_results)
        
    except ImportError as e:
        logger.error(f"Biopython package not installed: {str(e)}")
        return "PubMed search is not available. Please install biopython package."
    except Exception as e:
        logger.error(f"Error searching PubMed with query '{query}': {str(e)}", exc_info=True)
        return f"An error occurred while searching PubMed: {str(e)}"
