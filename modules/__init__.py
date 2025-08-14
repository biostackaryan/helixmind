"""
Helix Mind Modules
------------------
Core imports for the Helix Mind Bioinformatics Toolkit.
Author: Aryan Dutt (https://github.com/biostackaryan)
License: GNU GPL v3
"""

# Core functions
from .fasta_parser import parse_fasta_file
from .kegg_module import search_kegg_pathway
from .pubmed_module import fetch_pubmed_articles
from .structure_viewer import render_structure
from .blast_module import run_blast

# Optional extended functions
try:
    from .kegg_module import fetch_kegg_details, fetch_kegg_enzyme_info
except ImportError:
    fetch_kegg_details = None
    fetch_kegg_enzyme_info = None

# Backward compatibility
def search_kegg(query):
    return search_kegg_pathway(query)

__all__ = [
    "parse_fasta_file",
    "search_kegg",
    "search_kegg_pathway",
    "fetch_kegg_details",
    "fetch_kegg_enzyme_info",
    "fetch_pubmed_articles",
    "render_structure",
    "run_blast",
]
