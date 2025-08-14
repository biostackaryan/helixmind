"""
HelixMind - Bioinformatics Toolkit
Author: Aryan Dutt (https://github.com/biostackaryan)
License: GNU GPL v3
Copyright (C) 2025 Aryan Dutt
"""

from Bio import SeqIO

def parse_fasta_file(file_handle_or_path, desired_length=None):
    """
    Parse a FASTA file given a path or file-like object.
    Returns:
        - records: list of SeqRecord objects
        - stats: dict with counts and length info
    """
    try:
        nucleotide_records = list(SeqIO.parse(file_handle_or_path, "fasta"))
    except Exception as e:
        return [], {"error": str(e)}

    lengths = [len(record.seq) for record in nucleotide_records]

    stats = {
        "total_sequences": len(nucleotide_records),
        "shortest": min(lengths) if lengths else 0,
        "longest": max(lengths) if lengths else 0,
        "average": sum(lengths)//len(lengths) if lengths else 0,
        "filtered_count": 0
    }

    if desired_length is not None:
        filtered = [record for record in nucleotide_records if len(record.seq) == desired_length]
        stats["filtered_count"] = len(filtered)
        return filtered, stats

    return nucleotide_records, stats
