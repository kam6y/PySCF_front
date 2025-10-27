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
        from langchain_tavily import TavilySearch
        from agent.utils import get_tavily_api_key

        logger.info(f"Searching Tavily with query: '{query}', max_results: {max_results}")

        # Get Tavily API key using unified utility function
        api_key = get_tavily_api_key()

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
        from agent.utils import get_research_email

        logger.info(f"Searching PubMed with query: '{query}', max_results: {max_results}")

        # Set email for Entrez (required by NCBI)
        Entrez.email = get_research_email()
        
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


@tool
def search_chemrxiv(query: str, max_results: int = 5) -> str:
    """
    Search ChemRxiv database for chemistry preprint papers.

    Args:
        query: Search query string (e.g., "catalysis mechanisms")
        max_results: Maximum number of results to return (default: 5)

    Returns:
        Formatted string containing paper titles, authors, DOIs, summaries, and PDF/HTML links.
        Returns an error message if the search fails.

    Example:
        >>> result = search_chemrxiv("organocatalysis", max_results=3)
    """
    try:
        import chemrxiv

        logger.info(f"Searching ChemRxiv with query: '{query}', max_results: {max_results}")

        # Create ChemRxiv client
        client = chemrxiv.Client()

        # Perform search
        search = chemrxiv.Search(
            term=query,
            limit=max_results,
            sort=chemrxiv.SortCriterion.PUBLISHED_DATE_DESC
        )

        results = list(client.results(search))

        if not results:
            logger.warning(f"No papers found for query: '{query}'")
            return "No relevant papers found on ChemRxiv for the given query."

        formatted_results = []

        for idx, paper in enumerate(results, 1):
            try:
                # Extract paper information
                title = paper.title or "No title"
                authors = ", ".join([str(author) for author in paper.authors]) if paper.authors else "Unknown authors"
                doi = paper.doi or "No DOI"
                abstract = paper.abstract or "No abstract available"

                # Extract PDF/HTML URLs
                pdf_url = "Not available"
                html_url = "Not available"

                if hasattr(paper, 'item') and hasattr(paper.item, 'asset'):
                    if hasattr(paper.item.asset, 'original'):
                        pdf_url = paper.item.asset.original
                    if hasattr(paper.item.asset, 'html'):
                        html_url = paper.item.asset.html

                logger.debug(f"Paper {idx}: {title[:50]}... | DOI: {doi}")

                formatted_results.append(
                    f"**Title:** {title}\n"
                    f"**Authors:** {authors}\n"
                    f"**DOI:** {doi}\n"
                    f"**Abstract:** {abstract[:500]}{'...' if len(abstract) > 500 else ''}\n"
                    f"**PDF Link:** {pdf_url}\n"
                    f"**HTML Link:** {html_url}"
                )
            except Exception as paper_error:
                logger.warning(f"Error parsing ChemRxiv paper {idx}: {str(paper_error)}")
                continue

        if not formatted_results:
            logger.warning(f"Could not parse any papers for query: '{query}'")
            return "Found papers but could not parse their details."

        logger.info(f"Successfully retrieved {len(formatted_results)} papers from ChemRxiv")
        return "\n\n---\n\n".join(formatted_results)

    except ImportError as e:
        logger.error(f"ChemRxiv package not installed: {str(e)}")
        return "ChemRxiv search is not available. Please install chemrxiv package."
    except Exception as e:
        logger.error(f"Error searching ChemRxiv with query '{query}': {str(e)}", exc_info=True)
        return f"An error occurred while searching ChemRxiv: {str(e)}"


@tool
def search_openalex(query: str, max_results: int = 5) -> str:
    """
    Search OpenAlex database for comprehensive academic papers across all fields.

    Args:
        query: Search query string (e.g., "quantum computing algorithms")
        max_results: Maximum number of results to return (default: 5)

    Returns:
        Formatted string containing paper titles, authors, publication years, citation counts,
        open access status, and DOIs. Returns an error message if the search fails.

    Example:
        >>> result = search_openalex("machine learning chemistry", max_results=3)
    """
    try:
        from pyalex import Works
        import pyalex
        from agent.utils import get_research_email

        logger.info(f"Searching OpenAlex with query: '{query}', max_results: {max_results}")

        # Configure polite pool for better performance
        pyalex.config.email = get_research_email()

        # Perform search
        search_results = Works().search(query).get()

        if not search_results:
            logger.warning(f"No papers found for query: '{query}'")
            return "No relevant papers found on OpenAlex for the given query."

        # Limit to max_results
        results = search_results[:max_results]

        formatted_results = []

        for idx, work in enumerate(results, 1):
            try:
                # Extract work information
                title = work.get('title', 'No title')

                # Extract authors
                authorships = work.get('authorships', [])
                authors = ", ".join([
                    authorship.get('author', {}).get('display_name', 'Unknown')
                    for authorship in authorships[:5]
                ])
                if len(authorships) > 5:
                    authors += " et al."

                # Publication year
                pub_year = work.get('publication_year', 'Unknown year')

                # Citation count
                citation_count = work.get('cited_by_count', 0)

                # Open access info
                oa_info = work.get('open_access', {})
                oa_status = oa_info.get('oa_status', 'closed')
                oa_url = oa_info.get('oa_url', '')

                # DOI
                doi = work.get('doi', 'No DOI')
                if doi and not doi.startswith('http'):
                    doi_url = f"https://doi.org/{doi.replace('https://doi.org/', '')}"
                else:
                    doi_url = doi

                logger.debug(f"Paper {idx}: {title[:50]}... | Citations: {citation_count}")

                formatted_results.append(
                    f"**Title:** {title}\n"
                    f"**Authors:** {authors}\n"
                    f"**Year:** {pub_year}\n"
                    f"**Citations:** {citation_count}\n"
                    f"**Open Access:** {oa_status.upper()}" +
                    (f" - {oa_url}" if oa_url else "") + "\n"
                    f"**DOI:** {doi_url}"
                )
            except Exception as work_error:
                logger.warning(f"Error parsing OpenAlex work {idx}: {str(work_error)}")
                continue

        if not formatted_results:
            logger.warning(f"Could not parse any works for query: '{query}'")
            return "Found works but could not parse their details."

        logger.info(f"Successfully retrieved {len(formatted_results)} works from OpenAlex")
        return "\n\n---\n\n".join(formatted_results)

    except ImportError as e:
        logger.error(f"PyAlex package not installed: {str(e)}")
        return "OpenAlex search is not available. Please install pyalex package."
    except Exception as e:
        logger.error(f"Error searching OpenAlex with query '{query}': {str(e)}", exc_info=True)
        return f"An error occurred while searching OpenAlex: {str(e)}"


@tool
def search_semantic_scholar(query: str, max_results: int = 5) -> str:
    """
    Search Semantic Scholar database for AI/ML/science papers with rich citation data.

    Args:
        query: Search query string (e.g., "neural networks for molecular dynamics")
        max_results: Maximum number of results to return (default: 5)

    Returns:
        Formatted string containing paper titles, authors, publication years, citation counts,
        abstracts, and paper URLs. Returns an error message if the search fails.

    Example:
        >>> result = search_semantic_scholar("deep learning protein folding", max_results=3)
    """
    try:
        import time
        from semanticscholar import SemanticScholar
        from semanticscholar.SemanticScholarException import BadQueryParametersException

        logger.info(f"Searching Semantic Scholar with query: '{query}', max_results: {max_results}")

        # Add delay to respect rate limits (1 second between requests)
        time.sleep(1.0)

        # Initialize Semantic Scholar client with retry disabled to avoid long waits on rate limits
        sch = SemanticScholar(retry=False)

        # Perform search
        try:
            search_results = sch.search_paper(query)
        except Exception as search_error:
            # Handle various API errors gracefully
            error_str = str(search_error)

            # Check for rate limit errors (HTTP 429)
            if ("429" in error_str or "Too Many Requests" in error_str or
                "rate limit" in error_str.lower() or "RetryError" in error_str):
                logger.warning(f"Semantic Scholar rate limit reached for query: '{query}'")
                return "Semantic Scholar rate limit reached. Skipping this source to continue research with other databases."

            # Check for timeout errors (HTTP 504 Gateway Timeout)
            elif ("504" in error_str or "Gateway Timeout" in error_str or
                  "timed out" in error_str.lower() or "timeout" in error_str.lower()):
                logger.warning(f"Semantic Scholar request timed out for query: '{query}'")
                return "Semantic Scholar request timed out. Skipping this source to continue research with other databases."

            # For other errors, log and return graceful message
            else:
                logger.error(f"Error searching Semantic Scholar: {str(search_error)}")
                return f"Semantic Scholar search temporarily unavailable. Continuing with other sources."

        # Collect results up to max_results
        papers = []
        for paper in search_results:
            if len(papers) >= max_results:
                break
            papers.append(paper)

        if not papers:
            logger.warning(f"No papers found for query: '{query}'")
            return "No relevant papers found on Semantic Scholar for the given query."

        formatted_results = []

        for idx, paper in enumerate(papers, 1):
            try:
                # Extract paper information
                title = paper.title or "No title"

                # Extract authors
                authors_list = paper.authors or []
                authors = ", ".join([author.name for author in authors_list[:5]])
                if len(authors_list) > 5:
                    authors += " et al."

                # Publication year
                year = paper.year or "Unknown year"

                # Citation count
                citation_count = paper.citationCount or 0

                # Abstract
                abstract = paper.abstract or "No abstract available"

                # Paper URL
                paper_url = paper.url or "No URL available"

                # PDF URL if available
                pdf_url = "Not available"
                if hasattr(paper, 'openAccessPdf') and paper.openAccessPdf:
                    pdf_url = paper.openAccessPdf.get('url', 'Not available')

                logger.debug(f"Paper {idx}: {title[:50]}... | Citations: {citation_count}")

                formatted_results.append(
                    f"**Title:** {title}\n"
                    f"**Authors:** {authors}\n"
                    f"**Year:** {year}\n"
                    f"**Citations:** {citation_count}\n"
                    f"**Abstract:** {abstract[:500]}{'...' if len(abstract) > 500 else ''}\n"
                    f"**Paper URL:** {paper_url}\n"
                    f"**PDF URL:** {pdf_url}"
                )
            except Exception as paper_error:
                logger.warning(f"Error parsing Semantic Scholar paper {idx}: {str(paper_error)}")
                continue

        if not formatted_results:
            logger.warning(f"Could not parse any papers for query: '{query}'")
            return "Found papers but could not parse their details."

        logger.info(f"Successfully retrieved {len(formatted_results)} papers from Semantic Scholar")
        return "\n\n---\n\n".join(formatted_results)

    except ImportError as e:
        logger.error(f"Semantic Scholar package not installed: {str(e)}")
        return "Semantic Scholar search is not available. Please install semanticscholar package."
    except Exception as e:
        logger.error(f"Error searching Semantic Scholar with query '{query}': {str(e)}", exc_info=True)
        return f"An error occurred while searching Semantic Scholar: {str(e)}"
