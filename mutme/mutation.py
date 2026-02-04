from __future__ import annotations

import csv
import os

from pathlib import Path
from dataclasses import dataclass
from typing import Callable, Mapping, Iterable, Sequence


class InputFormatError(ValueError):
    """Raised when an input TSV is missing required columns or is otherwise malformed."""


def _split_csv_cell(cell: str | None) -> set[str]:
    """Split a comma-separated Nextclade cell into a set of tokens."""

    if cell is None:
        return set()
    cell = cell.strip()
    if not cell:
        return set()
    return {tok.strip() for tok in cell.split(",") if tok.strip()}


def _merge_field_value(existing: str, incoming: str) -> str:
    """
    Merge two non-empty string values for the same field.

    If they match (case-sensitive), keep one. If they differ, join with '; ' while
    avoiding duplicates.
    """

    if not existing:
        return incoming
    if not incoming or incoming == existing:
        return existing

    parts = [p.strip() for p in existing.split(";") if p.strip()]
    if incoming not in parts:
        parts.append(incoming)
    return "; ".join(parts)


@dataclass(frozen=True, slots=True)
class MutationAnnotations:
    """
    Annotations for a single mutation from the input annotation TSV.

    Attributes
    ----------
    mutation:
        The mutation key (string, should match Nextclade output format).
    fields:
        Mapping of arbitrary annotation column header -> string value.
        Example: {"phenotype": "resistant", "Ganciclovir": "8x", "Aciclovir": "2x"}
    """

    mutation: str
    fields: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class SequenceMutationHit:
    """
    A mutation detected in a sequence and the associated annotations.

    Attributes
    ----------
    mutation:
        Mutation string (e.g., 'N:R203K', 'S:214:EPE', 'N:E31-').
    annotations:
        Mapping of annotation fields for this mutation.
    """

    mutation: str
    annotations: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class SequenceAnnotationReport:
    """
    Per-sequence report of AA substitution/indel hits against an annotation table.
    """

    seq_name: str
    qc_status: str
    aa_substitutions: tuple[str, ...]
    aa_indels: tuple[str, ...]
    hits: tuple[SequenceMutationHit, ...]


def load_annotation_table(
    annotation_tsv: str | os.PathLike[str],
    *,
    mutation_col: str = "mutation",
    normalize: Callable[[str], str] | None = None,
    delimiter: str = ","
) -> dict[str, dict[str, str]]:
    """
    Load an annotation CSV with a required `mutation` column and any number of
    additional columns.

    Returns
    -------
    dict[str, dict[str, str]]
        Mapping: mutation -> {field_name -> field_value}

    Merge semantics
    ---------------
    If a mutation appears multiple times:
    - for each field, non-empty values are merged
    - if conflicting non-empty values appear, they are concatenated using '; '

    Raises
    ------
    FileNotFoundError:
        If the TSV does not exist.
    InputFormatError:
        If header is missing or `mutation` column is missing/empty.
    """

    path = Path(annotation_tsv)
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"Annotation CSV not found: {path}")

    out: dict[str, dict[str, str]] = {}

    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None:
            raise InputFormatError(f"Annotation CSV has no header row: {path}")

        fieldnames = list(reader.fieldnames)
        if mutation_col not in fieldnames:
            raise InputFormatError(f"Annotation CSV missing required column '{mutation_col}': {path}")

        annotation_fields = [c for c in fieldnames if c != mutation_col]
        if not annotation_fields:
            # 0 is allowed, but it's usually accidental. We'll allow it but be explicit.
            print("Warning: no annotation fields for provided in annotation file")
            pass

        for line_no, row in enumerate(reader, start=2):
            raw_mut = (row.get(mutation_col) or "").strip()
            if not raw_mut:
                raise InputFormatError(f"Empty '{mutation_col}' at {path}:{line_no}")

            mut = normalize(raw_mut) if normalize else raw_mut

            if mut not in out:
                out[mut] = {}

            for col in annotation_fields:
                val = (row.get(col) or "").strip()
                if not val:
                    continue
                existing = out[mut].get(col, "")
                out[mut][col] = _merge_field_value(existing, val)

    return out


def compare_nextclade_to_annotations(
    nextclade_tsv: str | os.PathLike[str],
    annotation_csv: str | os.PathLike[str],
    *,
    # Nextclade columns
    seq_name_col: str = "seqName",
    qc_status_col: str = "qc.OverallStatus",
    aa_sub_col: str = "aaSubstitutions",
    aa_del_col: str = "aaDeletions",
    aa_ins_col: str = "aaInsertions",
    # Annotation table column
    mutation_col: str = "mutation",
    # Optional normalization hook (applied to both sides)
    normalize: Callable[[str], str] | None = None,
    # Optional delimiter for annotation file
    delimiter: str = ","
) -> list[SequenceAnnotationReport]:
    """
    Compare Nextclade AA substitutions/indels to an annotation table CSV.

    The annotation CSV must contain a `mutation` column and may contain any number
    of additional columns with arbitrary headers (all treated as string annotations).

    Matching is performed by **exact string match** after optional normalization.

    Returns
    -------
    list[SequenceAnnotationReport]
        One report per sequence in Nextclade TSV, including matched mutation hits
        and their annotation fields.

    Raises
    ------
    FileNotFoundError, InputFormatError
    """
    next_path = Path(nextclade_tsv)
    if not next_path.exists() or not next_path.is_file():
        raise FileNotFoundError(f"Nextclade output not found: {next_path}")

    annotations = load_annotation_table(
        annotation_csv,
        mutation_col=mutation_col,
        normalize=normalize,
        delimiter=delimiter
    )
    annotation_keys = set(annotations.keys())

    reports: list[SequenceAnnotationReport] = []

    with next_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise InputFormatError(f"Nextclade output has no header row: {next_path}")

        required = {seq_name_col, aa_sub_col, aa_del_col, aa_ins_col}
        missing = required - set(reader.fieldnames)
        if missing:
            raise InputFormatError(
                f"Nextclade output missing required columns {sorted(missing)}: {next_path}"
            )

        for line_no, row in enumerate(reader, start=2):
            seq_name = (row.get(seq_name_col) or "").strip()
            if not seq_name:
                raise InputFormatError(f"Empty '{seq_name_col}' at {next_path}:{line_no}")
            
            qc_status = (row.get(qc_status_col) or "").strip()

            subs = _split_csv_cell(row.get(aa_sub_col))
            dels = _split_csv_cell(row.get(aa_del_col))
            ins = _split_csv_cell(row.get(aa_ins_col))

            if normalize is not None:
                subs = {normalize(x) for x in subs}
                dels = {normalize(x) for x in dels}
                ins = {normalize(x) for x in ins}

            aa_indels = dels | ins
            all_aa = subs | aa_indels

            hits: list[SequenceMutationHit] = []
            for mut in sorted(all_aa):
                if mut in annotation_keys:
                    hits.append(SequenceMutationHit(mutation=mut, annotations=annotations[mut]))

            reports.append(
                SequenceAnnotationReport(
                    seq_name=seq_name,
                    qc_status=qc_status,
                    aa_substitutions=tuple(sorted(subs)),
                    aa_indels=tuple(sorted(aa_indels)),
                    hits=tuple(hits),
                )
            )

    return reports



def write_long_format_hits_tsv(
    reports: Sequence[SequenceAnnotationReport],
    out_tsv: str | Path,
    *,
    include_sequences_with_no_hits: bool = False,
    seq_name_col: str = "seqName",
    mutation_col: str = "mutation",
    # Optional: add these for auditing/debugging
    include_detected_lists: bool = False,
    detected_subs_col: str = "detectedAaSubstitutions",
    detected_indels_col: str = "detectedAaIndels",
    # Ordering / stability
    annotation_field_order: Sequence[str] | None = None,
) -> Path:
    """
    Write a flat TSV in *long* format: one row per (sequence, hit mutation).

    Output columns
    --------------
    Always included:
    - seqName (configurable via `seq_name_col`)
    - mutation (configurable via `mutation_col`)
    - all annotation columns present across `reports[*].hits[*].annotations`

    Optionally included:
    - detectedAaSubstitutions: comma-joined list of all AA substitutions Nextclade reported
    - detectedAaIndels: comma-joined list of all AA indels (deletions + insertions)

    Parameters
    ----------
    reports:
        Output of `compare_nextclade_to_annotations()`.
    out_tsv:
        Output path for the TSV.
    include_sequences_with_no_hits:
        If True, emits one row per sequence even when there are no hits, with the
        mutation column empty and annotations empty.
    include_detected_lists:
        If True, includes the detected mutation lists for auditing/debugging.
    annotation_field_order:
        Optional explicit ordering for annotation columns (any unspecified fields
        will be appended in sorted order).

    Returns
    -------
    Path
        The written file path.

    Notes
    -----
    - This writer unions annotation keys across all hits to build the header.
    - For a stable, deterministic output, headers are sorted unless an explicit
      `annotation_field_order` is provided.
    """
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Union all annotation fields across all hits.
    all_fields: set[str] = set()
    for r in reports:
        for h in r.hits:
            all_fields.update(h.annotations.keys())

    # Determine annotation header order.
    if annotation_field_order is None:
        anno_cols = sorted(all_fields)
    else:
        specified = [c for c in annotation_field_order if c in all_fields]
        remaining = sorted(all_fields - set(specified))
        anno_cols = specified + remaining

    # Build full header.
    header: list[str] = [seq_name_col, mutation_col]
    if include_detected_lists:
        header.extend([detected_subs_col, detected_indels_col])
    header.extend(anno_cols)

    def _join(items: Iterable[str]) -> str:
        # deterministic, no spaces
        return ",".join(items)

    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for r in reports:
            detected_subs = _join(r.aa_substitutions)
            detected_indels = _join(r.aa_indels)

            if not r.hits:
                if include_sequences_with_no_hits:
                    row = {seq_name_col: r.seq_name, mutation_col: ""}
                    if include_detected_lists:
                        row[detected_subs_col] = detected_subs
                        row[detected_indels_col] = detected_indels
                    # annotation fields left blank
                    writer.writerow(row)
                continue

            for h in r.hits:
                row = {seq_name_col: r.seq_name, mutation_col: h.mutation}
                if include_detected_lists:
                    row[detected_subs_col] = detected_subs
                    row[detected_indels_col] = detected_indels

                # Fill annotation fields (missing -> empty)
                for col in anno_cols:
                    row[col] = h.annotations.get(col, "")

                writer.writerow(row)

    return out_path
