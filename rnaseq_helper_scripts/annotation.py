#!/usr/bin/env python3
"""Transcript-to-gene mapping shared by the expression and differential expression tables.

The reference biomart export is the same file kallisto/DESeq2/Sleuth already use, so gene-level
numbers in ``expression_gene.tsv`` are consistent with ``de_gene.tsv`` by construction.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from schemas import NA

ENSEMBL_GENE_RE = re.compile(r"^ENS[A-Z]*G\d+(\.\d+)?$")
ENSEMBL_TRANSCRIPT_RE = re.compile(r"^ENS[A-Z]*T\d+(\.\d+)?$")
# kallisto targets for viral contigs carry a per-transcript suffix (RPMS1_1) that the biomart
# export does not have.
TARGET_SUFFIX_RE = re.compile(r"_\d+$")


def strip_version(value: str) -> str:
    """``ENSG00000160310.12`` -> ``ENSG00000160310``. Leaves non-Ensembl ids alone."""
    if not isinstance(value, str):
        return value
    if "." in value and (ENSEMBL_GENE_RE.match(value) or ENSEMBL_TRANSCRIPT_RE.match(value)):
        return value.split(".", 1)[0]
    return value


class Annotation:
    """tx2gene mapping with the gene id column auto-detected.

    Most of the lab's biomart exports are two columns, ``target_id`` and ``gene``, where ``gene``
    is a *symbol* rather than an Ensembl accession -- that is what Sleuth's ``aggregation_column``
    expects. When a column of Ensembl gene accessions is present we prefer it for ``gene_id`` and
    report ``source == "ensembl"``; otherwise ``gene_id`` falls back to the symbol and ``source``
    says so, which the validator surfaces instead of silently accepting non-Ensembl ids.
    """

    def __init__(self, mart_file: Path):
        mart_file = Path(mart_file)
        if not mart_file.is_file():
            raise FileNotFoundError(f"biomart export not found: {mart_file}")

        frame = pd.read_csv(mart_file, sep="\t", dtype=str)
        frame = frame.fillna("")
        if frame.shape[1] < 2:
            raise ValueError(
                f"{mart_file} must have at least 2 columns (transcript, gene); "
                f"found {list(frame.columns)}"
            )

        transcript_col = frame.columns[0]
        symbol_col = frame.columns[1]
        gene_id_col = self._detect_gene_id_column(frame)
        self.source = "ensembl" if gene_id_col is not None else "gene_symbol"
        if gene_id_col is None:
            gene_id_col = symbol_col

        table = pd.DataFrame({
            "transcript_id": frame[transcript_col].str.strip(),
            "gene_id": frame[gene_id_col].str.strip().map(strip_version),
            "gene_symbol": frame[symbol_col].str.strip(),
        })
        table = table[table["transcript_id"] != ""]
        table = table.drop_duplicates(subset=["transcript_id"], keep="first")
        self.table = table.reset_index(drop=True)

        self._by_target = {
            row.transcript_id: (row.gene_id, row.gene_symbol)
            for row in self.table.itertuples(index=False)
        }
        self._by_target_stripped = {
            TARGET_SUFFIX_RE.sub("", target): value
            for target, value in self._by_target.items()
        }
        self._by_symbol = {}
        self._by_gene_id = {}
        for row in self.table.itertuples(index=False):
            self._by_symbol.setdefault(row.gene_symbol, row.gene_id)
            self._by_gene_id.setdefault(row.gene_id, row.gene_symbol)

    @staticmethod
    def _detect_gene_id_column(frame: pd.DataFrame) -> str | None:
        for column in frame.columns:
            sample = frame[column].head(500)
            sample = sample[sample != ""]
            if sample.empty:
                continue
            hits = sample.map(lambda v: bool(ENSEMBL_GENE_RE.match(v))).mean()
            if hits > 0.5:
                return column
        return None

    def lookup_target(self, target_id: str) -> tuple[str, str]:
        """Gene id and symbol for a kallisto/Sleuth target, with suffix fallback."""
        hit = self._by_target.get(target_id)
        if hit is None:
            hit = self._by_target_stripped.get(TARGET_SUFFIX_RE.sub("", str(target_id)))
        return hit if hit is not None else (NA, NA)

    def annotate_targets(self, target_ids: pd.Series) -> pd.DataFrame:
        """Vectorized :meth:`lookup_target` over a column of target ids."""
        direct = target_ids.map(self._by_target)
        unmatched = direct.isna()
        if unmatched.any():
            stripped = target_ids[unmatched].str.replace(TARGET_SUFFIX_RE, "", regex=True)
            direct.loc[unmatched] = stripped.map(self._by_target_stripped)
        pairs = direct.map(lambda v: v if isinstance(v, tuple) else (NA, NA))
        return pd.DataFrame({
            "gene_id": [p[0] for p in pairs],
            "gene_symbol": [p[1] for p in pairs],
        }, index=target_ids.index)

    def gene_id_for_symbol(self, symbol: str) -> str:
        """Ensembl gene id for a symbol, or the symbol itself when no mapping exists.

        DESeq2 and GSEA both key on the biomart ``gene`` column, so this is the only route from
        their output back to a stable gene identifier.
        """
        return self._by_symbol.get(symbol, symbol)

    def annotate_symbols(self, symbols: pd.Series) -> pd.Series:
        return symbols.map(lambda s: self._by_symbol.get(s, s))

    def symbol_for_gene_id(self, gene_id: str) -> str:
        return self._by_gene_id.get(gene_id, NA)

    def resolve_gene(self, label: str) -> tuple[str, str]:
        """``(gene_id, gene_symbol)`` for a DESeq2/GSEA row label of either kind."""
        label = str(label)
        if ENSEMBL_GENE_RE.match(label):
            gene_id = strip_version(label)
            return gene_id, self.symbol_for_gene_id(gene_id)
        return self._by_symbol.get(label, label), label
