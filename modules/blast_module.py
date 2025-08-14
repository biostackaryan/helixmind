"""
HelixMind - Bioinformatics Toolkit
Author: Aryan Dutt (https://github.com/biostackaryan)
License: GNU GPL v3
Copyright (C) 2025 Aryan Dutt
"""


import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

def run_blast(
    fasta_path,
    blast_type,
    database,
    evalue=0.001,
    output_path="blast_output.txt",
    threads=4,
    chunk_size=500,
    cleanup=True
):
    """
    Run local BLAST with chunking and parallel execution for large FASTA files.

    Args:
        fasta_path (str or Path): Path to input FASTA file.
        blast_type (str): BLAST program name (e.g., 'blastn', 'blastp').
        database (str): Path to BLAST database prefix.
        evalue (float): E-value threshold.
        output_path (str or Path): File to write combined BLAST output.
        threads (int): Number of threads per BLAST run.
        chunk_size (int): Number of sequences per chunk.
        cleanup (bool): Whether to delete chunk files after completion.

    Returns:
        str: Path to combined BLAST output file.

    Raises:
        FileNotFoundError: If required BLAST DB files are missing.
        RuntimeError: If any BLAST chunk fails.
    """
    fasta_path = Path(fasta_path)
    output_path = Path(output_path)
    temp_dir = Path("blast_chunks")
    temp_dir.mkdir(exist_ok=True)

    # Check BLAST DB files
    required_db_exts = [".nin", ".nhr", ".nsq"]  # for nucleotide DBs; adjust if protein DB
    missing_files = [database + ext for ext in required_db_exts if not Path(database + ext).exists()]
    if missing_files:
        raise FileNotFoundError(f"Missing BLAST DB files: {missing_files}")

    # Split FASTA into chunks
    def split_fasta(file_path, chunk_size):
        chunk_files = []
        chunk = []
        seq_count = 0
        chunk_index = 1

        with open(file_path) as f:
            for line in f:
                if line.startswith(">"):
                    seq_count += 1
                    if seq_count > chunk_size:
                        chunk_file = temp_dir / f"chunk_{chunk_index}.fasta"
                        chunk_file.write_text("".join(chunk))
                        chunk_files.append(chunk_file)
                        chunk_index += 1
                        chunk = []
                        seq_count = 1  # current sequence counts for next chunk
                chunk.append(line)

            # Write last chunk
            if chunk:
                chunk_file = temp_dir / f"chunk_{chunk_index}.fasta"
                chunk_file.write_text("".join(chunk))
                chunk_files.append(chunk_file)

        return chunk_files

    chunk_files = split_fasta(fasta_path, chunk_size)

    # Function to run BLAST on a chunk
    def blast_chunk(chunk_file):
        chunk_output = temp_dir / (chunk_file.stem + "_out.txt")
        cmd = [
            blast_type,
            "-query", str(chunk_file),
            "-db", database,
            "-evalue", str(evalue),
            "-outfmt", "6 qseqid sseqid stitle pident length evalue bitscore",
            "-num_threads", str(threads),
            "-out", str(chunk_output)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"BLAST failed on {chunk_file}:\n{result.stderr}")
        return chunk_output

    # Run all chunks in parallel
    with ThreadPoolExecutor(max_workers=min(len(chunk_files), threads)) as executor:
        outputs = list(executor.map(blast_chunk, chunk_files))

    # Merge chunk outputs
    with open(output_path, "w") as combined_out:
        for out_file in outputs:
            combined_out.write(out_file.read_text())

    # Optionally cleanup chunk files
    if cleanup:
        for f in chunk_files + outputs:
            try:
                f.unlink()
            except Exception:
                pass
        try:
            temp_dir.rmdir()
        except Exception:
            pass

    print(f"BLAST completed. Results saved to {output_path}")

    return str(output_path)
