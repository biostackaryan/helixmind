"""
HelixMind - Bioinformatics Toolkit
Author: Aryan Dutt (https://github.com/biostackaryan)
License: GNU GPL v3
Copyright (C) 2025 Aryan Dutt
"""

import requests

def fetch_pubmed_articles(query, max_results=5):
    """
    Fetches PubMed article IDs and basic info based on a query.

    Args:
        query (str): Search term (gene name, keyword, etc.)
        max_results (int): Maximum number of articles to retrieve

    Returns:
        dict: {
            "status": "ok" or "error",
            "results": list of dicts with keys: id, title, source, pubdate
            "error_message": str (if error)
        }
    """
    try:
        # Step 1: Search PubMed IDs matching query
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json"
        }
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        id_list = response.json().get("esearchresult", {}).get("idlist", [])

        # Step 2: Fetch summaries/details for each ID
        if not id_list:
            return {"status": "ok", "results": []}

        ids_str = ",".join(id_list)
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {
            "db": "pubmed",
            "id": ids_str,
            "retmode": "json"
        }
        summary_response = requests.get(summary_url, params=params)
        summary_response.raise_for_status()
        summaries = summary_response.json().get("result", {})

        results = []
        for pmid in id_list:
            item = summaries.get(pmid, {})
            title = item.get("title", f"PubMed Article {pmid}")
            source = item.get("source", "PubMed")
            pubdate = item.get("pubdate", "Unknown")
            results.append({
                "id": pmid,
                "title": title,
                "source": source,
                "pubdate": pubdate
            })

        return {"status": "ok", "results": results}

    except requests.RequestException as e:
        return {"status": "error", "error_message": f"Network error: {e}"}
    except Exception as e:
        return {"status": "error", "error_message": str(e)}
