
from csv import DictReader, DictWriter
from typing import TextIO

# Alignments with e-values greater than 1 are low-quality alignments and associated with
# a high rate of false-positives. These should be filtered at all alignment steps.
MAX_EVALUE_THRESHOLD = 1

# blastn output format 6 as documented in
# http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# it's also the format of our GSNAP and RAPSEARCH2 output
_BLASTN_OUTPUT_6_SCHEMA = [
    ("qseqid", str),
    ("sseqid", str),
    ("pident", float),
    ("length", int),
    ("mismatch", int),
    ("gapopen", int),
    ("qstart", int),
    ("qend", int),
    ("sstart", int),
    ("send", int),
    ("evalue", float),
    ("bitscore", float),
]

# Additional blastn output columns.
_BLASTN_OUTPUT_6_NT_SCHEMA = _BLASTN_OUTPUT_6_SCHEMA + [
    ("qlen", int),      # query sequence length, helpful for computing qcov
    ("slen", int),      # subject sequence length, so far unused in IDseq
]

# Re-ranked output of blastn.  One row per query.  Two additional columns.
_BLASTN_OUTPUT_6_NT_RERANKED_SCHEMA = _BLASTN_OUTPUT_6_NT_SCHEMA + [
    ("qcov", float),     # fraction of query covered by the optimal set of HSPs
    ("hsp_count", int),   # cardihnality of optimal fragment cover;  see BlastCandidate
]


class BlastnOutput6Reader(DictReader):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _BLASTN_OUTPUT_6_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")


class BlastnOutput6Writer(DictWriter):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _BLASTN_OUTPUT_6_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")


class BlastnOutput6NTReader(DictReader):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _BLASTN_OUTPUT_6_NT_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")


class BlastnOutput6NTWriter(DictWriter):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _BLASTN_OUTPUT_6_NT_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")


class BlastnOutput6NTRerankedReader(DictReader):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _BLASTN_OUTPUT_6_NT_RERANKED_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")


class BlastnOutput6NTRerankedWriter(DictWriter):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _BLASTN_OUTPUT_6_NT_RERANKED_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")


_HIT_SUMMARY_SCHEMA = [
    ("read_id", str),
    ("level", int),
    ("taxid", int),
    ("accession_id", str),
    ("species_taxid", int),
    ("genus_taxid", int),
    ("family_taxid", int),
]


_HIT_SUMMARY_MERGED_SCHEMA = _HIT_SUMMARY_SCHEMA + [
    ("contig_id", str),
    ("contig_accession_id", str),
    ("contig_species_taxid", int),
    ("contig_genus_taxid", int),
    ("contig_family_taxid", int),
    ("from_assembly", str),
    ("source_count_type", str),
]


class HitSummaryReader(DictReader):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _HIT_SUMMARY_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")


class HitSummaryWriter(DictWriter):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _HIT_SUMMARY_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")


class HitSummaryMergedReader(DictReader):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _HIT_SUMMARY_MERGED_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")


class HitSummaryMergedWriter(DictWriter):
    def __init__(self, f: TextIO) -> None:
        fieldnames = [name for name, _ in _HIT_SUMMARY_MERGED_SCHEMA]
        super().__init__(f, fieldnames, delimiter="\t")
