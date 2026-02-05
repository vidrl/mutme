from __future__ import annotations

import csv
import os
import re
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path


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
        Mutation string in Nextclade format (e.g., 'N:R203K', 'S:214:EPE', 'N:E31-').
    comments:
        Comments if provided in the annotations table for this mutation.
    annotations:
        Mapping of annotation fields for this mutation.
    """

    mutation: str
    comments: str
    annotations: Mapping[str, str]


@dataclass(frozen=True, slots=True)
class SequenceAnnotationReport:
    """
    Per-sequence report of AA substitution/indel/stop hits against an annotation table.
    """

    seq_name: str
    qc_status: str
    aa_substitutions: tuple[str, ...]
    aa_indels: tuple[str, ...]
    aa_stop_codons: tuple[str, ...]
    hits: tuple[SequenceMutationHit, ...]


def load_annotation_table(
    annotation_tsv: str | os.PathLike[str],
    *,
    mutation_col: str = "mutation",
    normalize: Callable[[str], str] | None = None,
    delimiter: str = ",",
) -> dict[str, dict[str, str]]:
    """
    Load an annotation CSV/TSV with a required mutation column and any number of
    additional columns.

    The mutation column is resolved case-insensitively, so either "mutation" or
    "Mutation" (or any other casing) will be accepted.

    Parameters
    ----------
    annotation_tsv:
        Path to the annotation file.
    mutation_col:
        Preferred name of the mutation column. Resolution is case-insensitive:
        e.g. "mutation" will match a header of "Mutation".
    normalize:
        Optional normalizer applied to each mutation string after stripping.
    delimiter:
        Delimiter used by the file ("," for CSV, "\\t" for TSV, etc).

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
        If the file does not exist.
    InputFormatError:
        If header is missing or the mutation column is missing/empty.
    """
    path = Path(annotation_tsv)
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"Annotation file not found: {path}")

    out: dict[str, dict[str, str]] = {}

    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None:
            raise InputFormatError(f"Annotation file has no header row: {path}")

        fieldnames = list(reader.fieldnames)

        # Resolve mutation column case-insensitively.
        # If there are multiple headers differing only by case, fail loudly.
        target = mutation_col.strip()
        if not target:
            raise ValueError("mutation_col must be a non-empty string")

        matches = [c for c in fieldnames if c.casefold() == target.casefold()]
        if not matches:
            raise InputFormatError(
                f"Annotation file missing required column '{mutation_col}' (case-insensitive): {path}"
            )
        if len(matches) > 1:
            raise InputFormatError(
                f"Ambiguous mutation column: multiple headers match '{mutation_col}' case-insensitively: "
                f"{matches} in {path}"
            )

        resolved_mutation_col = matches[0]

        comments_matches = [c for c in fieldnames if c.casefold() == "comment"]
        if len(comments_matches) > 1:
            raise InputFormatError(
                f"Ambiguous comment column: multiple headers match 'comment' case-insensitively: "
                f"{comments_matches} in {path}"
            )
        resolved_comments_col: str | None = comments_matches[0] if comments_matches else None

        annotation_fields = [c for c in fieldnames if c != resolved_mutation_col]
        if not annotation_fields:
            print("Warning: no annotation fields provided in annotation file")

        for line_no, row in enumerate(reader, start=2):
            raw_mut = (row.get(resolved_mutation_col) or "").strip()
            if not raw_mut:
                raise InputFormatError(f"Empty '{resolved_mutation_col}' at {path}:{line_no}")

            mut = normalize(raw_mut) if normalize else raw_mut

            bucket = out.setdefault(mut, {})

            for col in annotation_fields:
                val = (row.get(col) or "").strip()
                if not val:
                    continue

                key = "comment" if (resolved_comments_col and col == resolved_comments_col) else col
                existing = bucket.get(key, "")
                bucket[key] = _merge_field_value(existing, val)

    return out


def compare_nextclade_to_annotations(
    nextclade_tsv: str | os.PathLike[str],
    annotation_csv: str | os.PathLike[str],
    *,
    # Nextclade columns
    seq_name_col: str = "seqName",
    qc_status_col: str = "qc.overallStatus",
    aa_stop_col: str = "qc.stopCodons.stopCodons",
    aa_sub_col: str = "aaSubstitutions",
    aa_del_col: str = "aaDeletions",
    aa_ins_col: str = "aaInsertions",
    # Annotation table column
    mutation_col: str = "mutation",
    # Optional normalization hook (applied to both sides)
    normalize: Callable[[str], str] | None = None,
    # Optional delimiter for annotation file
    delimiter: str = ",",
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
        delimiter=delimiter,
    )
    annotation_keys = set(annotations.keys())

    # STOP CODONS: build lookup from annotation stop keys
    stop_lookup = _build_stop_lookup(annotation_keys)

    reports: list[SequenceAnnotationReport] = []

    with next_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise InputFormatError(f"Nextclade output has no header row: {next_path}")

        required = {seq_name_col, aa_sub_col, aa_del_col, aa_ins_col, aa_stop_col}
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

            stops = _split_csv_cell(row.get(aa_stop_col))

            if normalize is not None:
                subs = {normalize(x) for x in subs}
                dels = {normalize(x) for x in dels}
                ins = {normalize(x) for x in ins}

            aa_stops: set[str] = set()

            for raw in stops:
                parsed = _parse_nextclade_stop(raw)
                if parsed is None:
                    continue
                gene, pos = parsed

                mapped = stop_lookup.get((gene, pos))
                if mapped:
                    aa_stops.update(mapped)

            aa_indels = dels | ins
            all_aa = subs | aa_indels | aa_stops

            hits: list[SequenceMutationHit] = []

            for mut in sorted(all_aa):
                if mut not in annotation_keys:
                    continue

                anno = annotations[mut]
                mut_comment = (anno.get("comment") or "").strip()

                # exclude comments from generic annotations
                hit_annos = {k: v for k, v in anno.items() if k != "comment"}

                hits.append(
                    SequenceMutationHit(
                        mutation=mut,
                        comments=mut_comment,
                        annotations=hit_annos,
                    )
                )

            reports.append(
                SequenceAnnotationReport(
                    seq_name=seq_name,
                    qc_status=qc_status,
                    aa_substitutions=tuple(sorted(subs)),
                    aa_indels=tuple(sorted(aa_indels)),
                    aa_stop_codons=tuple(sorted(aa_stops)),
                    hits=tuple(hits),
                )
            )

    return reports


def write_long_format_table(
    reports: Sequence[SequenceAnnotationReport],
    output: str | Path,
    *,
    include_sequences_with_no_hits: bool = False,
    seq_name_col: str = "seq_name",
    seq_qc_col: str = "seq_quality",
    mutation_col: str = "mutation",
    # Optional: add these for auditing/debugging
    include_detected_lists: bool = False,
    detected_subs_col: str = "detected_aa_substitutions",
    detected_indels_col: str = "detected_aa_indels",
    include_mutation_comments: bool = False,
    comments_col: str = "comment",
    # Ordering / stability
    annotation_field_order: Sequence[str] | None = None,
    # Delimiter
    delimiter: str = ",",
) -> Path:
    """
    Write a flat CSV in *long* format: one row per (sequence, hit mutation).

    Output columns
    --------------
    Always included:
    - seq_name (configurable via `seq_name_col`)
    - seq_quality (configurable via `seq_quality_col`)
    - mutation (configurable via `mutation_col`)
    - all annotation columns present across `reports[*].hits[*].annotations`

    Optionally included:
    - detected_aa_substitutions: comma-joined list of all AA substitutions Nextclade reported
    - detected_aa_indels: comma-joined list of all AA indels (deletions + insertions)

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
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Union all annotation fields across hits (comments are now separate)
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
    header: list[str] = [seq_name_col, seq_qc_col, mutation_col]
    if include_detected_lists:
        header.extend([detected_subs_col, detected_indels_col])
    header.extend(anno_cols)
    if include_mutation_comments:
        header.append(comments_col)

    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()

        for r in reports:
            detected_subs = ",".join(r.aa_substitutions)
            detected_indels = ",".join(r.aa_indels)

            if not r.hits:
                if include_sequences_with_no_hits:
                    row = {
                        seq_name_col: r.seq_name,
                        seq_qc_col: r.qc_status,
                        mutation_col: "",
                    }
                    if include_detected_lists:
                        row[detected_subs_col] = detected_subs
                        row[detected_indels_col] = detected_indels
                    if include_mutation_comments:
                        row[comments_col] = ""
                    writer.writerow(row)
                continue

            for h in r.hits:
                row = {
                    seq_name_col: r.seq_name,
                    seq_qc_col: r.qc_status,
                    mutation_col: h.mutation,
                }

                if include_detected_lists:
                    row[detected_subs_col] = detected_subs
                    row[detected_indels_col] = detected_indels

                for col in anno_cols:
                    row[col] = h.annotations.get(col, "")

                if include_mutation_comments:
                    row[comments_col] = h.comments

                writer.writerow(row)

    return out_path


# Stop codons
# Because Nextclade doesn’t include the reference AA (it only gives ORF1a:4715), one cannot reconstruct N4715* reliably. So the clean approach is:
# - Pre-index annotation table’s stop-codon mutations into a lookup keyed by (gene, position)
# - Parse Nextclade’s stop list into (gene, position).
# - Use that lookup to convert Nextclade stop sites into the exact annotation mutation strings (so the existing exact-match logic works).

# Annotation: requires gene prefix; allows optional ref AA before position (ORF1a:N4715*)
_STOP_ANNO_RE = re.compile(r"^(?P<gene>[^:]+):(?P<aa>[A-Za-z\*])?(?P<pos>\d+)\*$")

# Nextclade: gene:pos (ORF1a:4715)
_STOP_NEXTCLADE_RE = re.compile(r"^(?P<gene>[^:]+):(?P<pos>\d+)$")


def _parse_annotation_stop(mut: str) -> tuple[str, int] | None:
    m = _STOP_ANNO_RE.match(mut.strip())
    if not m:
        return None
    return (m.group("gene"), int(m.group("pos")))


def _parse_nextclade_stop(item: str) -> tuple[str, int] | None:
    m = _STOP_NEXTCLADE_RE.match(item.strip())
    if not m:
        return None
    return (m.group("gene"), int(m.group("pos")))


def _build_stop_lookup(annotation_keys: set[str]) -> dict[tuple[str, int], list[str]]:
    """
    (gene, pos) -> [annotation stop mutation strings]
    """
    out: dict[tuple[str, int], list[str]] = defaultdict(list)
    for key in annotation_keys:
        if not key.endswith("*"):
            continue
        parsed = _parse_annotation_stop(key)
        if parsed is None:
            continue
        out[parsed].append(key)
    return dict(out)
