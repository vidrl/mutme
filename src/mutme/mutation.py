from __future__ import annotations

import csv
import os
import re
from collections import defaultdict
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
CONTROL_COLUMNS = {"comment", "linked_mutations"}


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
        Mutation string in Nextclade format (e.g. 'N:R203K', 'S:214:EPE', 'N:E31-').
        When wildcard matching is enabled, this is the *detected* mutation string
        (not the wildcard key from the annotation table).
    comments:
        Comments if provided in the annotations table for this mutation.
        When wildcard matching is enabled and multiple annotation keys match a
        detected mutation, comments are merged using the same merge semantics as
        other fields (concatenated with '; ' while avoiding duplicates).
    annotations:
        Mapping of annotation fields for this mutation (excluding comments).
        When wildcard matching is enabled and multiple annotation keys match a
        detected mutation, field values are merged.
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
    require_full_deletion_ranges: bool = True,
) -> dict[str, dict[str, str]]:
    """
    Load an annotation CSV/TSV with a required mutation column and any number of
    additional columns.

    The mutation column is resolved case-insensitively, so either "mutation" or
    "Mutation" (or any other casing) will be accepted.

    Additional supported input conveniences
    ---------------------------------------
    1) Deletion ranges:

       If the mutation cell matches "{gene}:del{start}-{stop}" (inclusive),
       the row is expanded into multiple canonical deletion keys:
           "{gene}:{start}-", "{gene}:{start+1}-", ..., "{gene}:{stop}-"

       Each expanded key receives the same annotation fields.

       This expansion is performed at load time so downstream matching logic is
       unchanged.

    2) Linked mutations (co-occurrence gating):

       If an optional column "linked_mutations" is present (case-insensitive),
       its value is stored under the special key "linked_mutations" and later
       used by the matcher to require co-occurrence. The value is merged using
       set-union semantics when a mutation key appears multiple times.

       For deletion range rows, this loader automatically injects an implied
       linked_mutations set consisting of all expanded deletion keys, giving the
       "all-or-nothing" semantics for the range without requiring the user to
       specify it manually.

    Parameters
    ----------
    annotation_tsv:
        Path to the annotation file.
    mutation_col:
        Preferred name of the mutation column. Resolution is case-insensitive.
    normalize:
        Optional normalizer applied to each non-range mutation string after stripping.
        (Range keys are expanded before normalization; expanded keys are not normalized.)
    delimiter:
        Delimiter used by the file ("," for CSV, "\\t" for TSV, etc).
    require_full_deletion_ranges:
        If True (default), deletion range inputs like "{gene}:del87-89" are expanded
        and also automatically treated as linked (all positions must co-occur).
        If False, ranges are expanded but not auto-linked; each position can match
        independently unless `linked_mutations` is explicitly provided.

    Returns
    -------
    dict[str, dict[str, str]]
        Mapping: mutation -> {field_name -> field_value}

    Merge semantics
    ---------------
    If a mutation appears multiple times:
    - for each annotation field, non-empty values are merged
    - if conflicting non-empty values appear, they are concatenated using '; '
    - for the control field "linked_mutations", values are merged by set union

    Raises
    ------
    FileNotFoundError:
        If the file does not exist.
    InputFormatError:
        If header is missing or the mutation column is missing/empty, or deletion
        range syntax is malformed.
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

        linked_matches = [c for c in fieldnames if c.casefold() == "linked_mutations"]
        if len(linked_matches) > 1:
            raise InputFormatError(
                f"Ambiguous linked_mutations column: multiple headers match 'linked_mutations' case-insensitively: "
                f"{linked_matches} in {path}"
            )
        resolved_linked_col: str | None = linked_matches[0] if linked_matches else None

        annotation_fields = [c for c in fieldnames if c != resolved_mutation_col]
        if not annotation_fields:
            print("Warning: no annotation fields provided in annotation file")

        for line_no, row in enumerate(reader, start=2):
            raw_mut = (row.get(resolved_mutation_col) or "").strip()
            if not raw_mut:
                raise InputFormatError(f"Empty '{resolved_mutation_col}' at {path}:{line_no}")

            # Expand deletion ranges (GENE:del{pos_start}-{pos_stop}) into per-position deletion keys (GENE:{pos}-).
            try:
                expanded = _expand_deletion_range_key(raw_mut)
            except InputFormatError as e:
                raise InputFormatError(f"{e} at {path}:{line_no}") from e

            if expanded:
                mut_keys = expanded
                auto_link_value = ";".join(mut_keys) if require_full_deletion_ranges else ""
            else:
                mut = normalize(raw_mut) if normalize else raw_mut
                mut_keys = [mut]
                auto_link_value = ""

            for mut_key in mut_keys:
                bucket = out.setdefault(mut_key, {})

                for col in annotation_fields:
                    val = (row.get(col) or "").strip()
                    if not val:
                        continue

                    if resolved_comments_col and col == resolved_comments_col:
                        key = "comment"
                    elif resolved_linked_col and col == resolved_linked_col:
                        key = "linked_mutations"
                    else:
                        key = col

                    if key == "linked_mutations":
                        bucket[key] = _merge_linked_mutations(bucket.get(key, ""), val)
                    else:
                        bucket[key] = _merge_field_value(bucket.get(key, ""), val)

                # Inject auto-linking for range-expanded rows even if the column wasn't present.
                if auto_link_value:
                    bucket["linked_mutations"] = _merge_linked_mutations(
                        bucket.get("linked_mutations", ""),
                        auto_link_value,
                    )

    return out


def compare_nextclade_to_annotations(
    nextclade_tsv: str | os.PathLike[str],
    annotation_csv: str | os.PathLike[str],
    # Nextclade columns
    seq_name_col: str = "seqName",
    qc_status_col: str = "qc.overallStatus",
    aa_stop_col: str = "qc.stopCodons.stopCodons",
    aa_sub_col: str = "aaSubstitutions",
    aa_del_col: str = "aaDeletions",
    aa_ins_col: str = "aaInsertions",
    # Annotation table column
    mutation_col: str = "mutation",
    # Optional delimiter for annotation file
    delimiter: str = ",",
    # Wildcard mutations
    allow_x_wildcards: bool = False,
    x_charset: str = AMINO_ACIDS,
    # Deletion range full or partial
    require_full_deletion_ranges: bool = True,
) -> list[SequenceAnnotationReport]:
    """
    Compare Nextclade AA substitutions/indels to an annotation table CSV.

    (Existing docstring omitted here for brevity; only behavior additions described.)

    Additions
    ---------
    - Supports co-occurrence gating via the annotation control field "linked_mutations":
      if a matched annotation dict contains "linked_mutations", the hit is emitted
      only if *all* required tokens are present in the sequence's detected mutation
      universe.

    Notes
    -----
    - "linked_mutations" is treated as a control field and is not emitted in
      SequenceMutationHit.annotations.
    """
    next_path = Path(nextclade_tsv)
    if not next_path.exists() or not next_path.is_file():
        raise FileNotFoundError(f"Nextclade output not found: {next_path}")

    annotations = load_annotation_table(
        annotation_csv,
        mutation_col=mutation_col,
        delimiter=delimiter,
        require_full_deletion_ranges=require_full_deletion_ranges,
    )
    annotation_keys = set(annotations.keys())

    if allow_x_wildcards:
        exact_keys, wildcard_rules = _compile_x_wildcard_rules(
            annotation_keys,
            x_charset=x_charset,
        )
    else:
        exact_keys = annotation_keys
        wildcard_rules = []

    # STOP CODONS: build lookup from annotation stop keys
    stop_lookup = _build_stop_lookup(annotation_keys)

    # DELETIONS: build lookup from annotation del keys
    del_lookup = _build_del_lookup(annotation_keys)

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

            aa_subs = _split_csv_cell(row.get(aa_sub_col))
            aa_ins = _split_csv_cell(row.get(aa_ins_col))

            aa_dels = _split_csv_cell(row.get(aa_del_col))
            aa_stops = _split_csv_cell(row.get(aa_stop_col))

            # Lookups because annotations accept different formats than Nextclade for DEL/STOP
            aa_dels_match: set[str] = set()
            for raw in aa_dels:
                parsed = _parse_nextclade_del(raw)
                if parsed is None:
                    continue
                mapped = del_lookup.get(parsed)
                if mapped:
                    aa_dels_match.update(mapped)

            aa_stops_match: set[str] = set()
            for raw in aa_stops:
                parsed = _parse_nextclade_stop(raw)
                if parsed is None:
                    continue
                mapped = stop_lookup.get(parsed)
                if mapped:
                    aa_stops_match.update(mapped)

            aa_indels = aa_dels | aa_ins  # keep raw Nextclade tokens for optional reporting

            # all_aa is the detected universe used both for matching and for linked gating.
            all_aa = aa_subs | aa_ins | aa_dels_match | aa_stops_match

            hits: list[SequenceMutationHit] = []

            for detected_mut in sorted(all_aa):
                matched_keys: list[str] = []

                # Exact match
                if detected_mut in exact_keys:
                    matched_keys.append(detected_mut)

                # Wildcard matches (annotation keys containing X)
                if wildcard_rules:
                    for rule in wildcard_rules:
                        if rule.pattern.match(detected_mut):
                            matched_keys.append(rule.anno_key)

                if not matched_keys:
                    continue

                # Merge all matching annotation dicts (exact + wildcard).
                merged = _merge_annotation_dicts_for_keys(annotations, matched_keys)

                # Linked gating: if linked_mutations is present, require all linked tokens
                # to be present in the sequence's detected mutation universe (all_aa).
                linked_raw = (merged.get("linked_mutations") or "").strip()
                if linked_raw:
                    required_links = _split_linked_mutations_cell(linked_raw)
                    if not required_links.issubset(all_aa):
                        continue

                merged_comment = (merged.get("comment") or "").strip()
                merged_annos = {
                    k: v for k, v in merged.items() if k not in ("comment", "linked_mutations")
                }

                hits.append(
                    SequenceMutationHit(
                        mutation=detected_mut,  # always the detected mutation
                        comments=merged_comment,
                        annotations=merged_annos,
                    )
                )

            reports.append(
                SequenceAnnotationReport(
                    seq_name=seq_name,
                    qc_status=qc_status,
                    aa_substitutions=tuple(sorted(aa_subs)),
                    aa_indels=tuple(sorted(aa_indels)),
                    aa_stop_codons=tuple(sorted(aa_stops)),
                    hits=tuple(hits),
                )
            )

    return reports


def write_long_format_table(
    reports: Sequence[SequenceAnnotationReport],
    output: str | Path,
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
    # Output columns
    all_annotation_fields: Sequence[str] | None = None,
) -> Path:
    """
    Write a flat CSV in *long* format: one row per (sequence, hit mutation).

    Output columns
    --------------
    Always included:
    - seq_name (configurable via `seq_name_col`)
    - seq_quality (configurable via `seq_qc_col`)
    - mutation (configurable via `mutation_col`)
    - one column per annotation field

    Annotation columns
    ------------------
    By default, annotation columns are the union of all annotation keys observed
    across `reports[*].hits[*].annotations`.

    If `all_annotation_fields` is provided, the output will *always* include those
    columns, even if some (or all) of them are not present in any hit. Missing
    values are emitted as empty strings.

    Column ordering is controlled by `annotation_field_order` if provided; any
    unspecified fields are appended in sorted order.

    The special annotation column "comment" is not included in `all_annotation_fields`
    and is emitted only when `include_mutation_comments=True`.

    Optionally included columns
    ----------------------------
    - detected_aa_substitutions: comma-joined list of all AA substitutions reported
      by Nextclade for the sequence
    - detected_aa_indels: comma-joined list of all AA indels (deletions + insertions)
    - comment: mutation-specific comments, if `include_mutation_comments=True`

    Parameters
    ----------
    reports:
        Output of `compare_nextclade_to_annotations()`.
    output:
        Output path for the CSV.
    include_sequences_with_no_hits:
        If True, emits one row per sequence even when there are no hits, with the
        mutation column empty and all annotation fields empty.
    include_detected_lists:
        If True, includes detected mutation lists for auditing/debugging.
    annotation_field_order:
        Optional explicit ordering for annotation columns. Fields not listed here
        will be appended in sorted order.
    all_annotation_fields:
        Optional explicit universe of annotation columns to emit. If provided,
        these columns are always included in the output header regardless of
        whether they appear in any hit.
    delimiter:
        Delimiter for the output file.

    Returns
    -------
    Path
        The written file path.

    Notes
    -----
    - This writer does not modify or normalize annotation values.
    - Providing `all_annotation_fields` is recommended when a stable schema is
      required across runs or datasets.
    """
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Determine the universe of annotation fields.
    if all_annotation_fields is None:
        # existing behavior: union across hits (comments are now separate)
        all_fields: set[str] = set()
        for r in reports:
            for h in r.hits:
                all_fields.update(h.annotations.keys())
    else:
        # new behavior: always output these columns
        all_fields = {f for f in all_annotation_fields if f}

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


def get_annotation_fields(
    annotations_csv: str | os.PathLike[str],
    mutation_col: str = "mutation",
    delimiter: str = ",",
) -> list[str]:
    """
    Extract annotation column names from an annotation CSV/TSV.

    This function reads only the header row and returns all column names except
    the mutation column. Resolution of the mutation column name is
    case-insensitive.

    Control columns
    ---------------
    The optional control columns "comment" and "linked_mutations" are excluded
    from the returned list. "comment" is emitted only when explicitly requested
    by writers. "linked_mutations" is used for co-occurrence gating and is not a
    user-facing annotation.

    Parameters
    ----------
    annotations_csv:
        Path to the annotation file.
    mutation_col:
        Name of the mutation column. Matching is case-insensitive.
    delimiter:
        Delimiter used by the file (\",\" for CSV, \"\\t\" for TSV, etc).

    Returns
    -------
    list[str]
        List of annotation column names in file order, excluding the mutation
        column and control columns.

    Raises
    ------
    FileNotFoundError:
        If the file does not exist.
    InputFormatError:
        If the file has no header row or the mutation column is missing or
        ambiguous.
    """
    path = Path(annotations_csv)
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None:
            raise InputFormatError(f"Annotation file has no header row: {path}")

    fieldnames = list(reader.fieldnames)

    matches = [c for c in fieldnames if c.casefold() == mutation_col.casefold()]
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

    return [
        c for c in fieldnames if c != resolved_mutation_col and c.casefold() not in CONTROL_COLUMNS
    ]


# Annotation accepts formats: {gene}:{pos} (Nextclade) and {gene}:{aa}{pos}* (general)
_STOP_NEXTCLADE_RE = re.compile(r"^(?P<gene>[^:]+):(?P<pos>\d+)$")
_STOP_ANNO_GENERAL_RE = re.compile(r"^(?P<gene>[^:]+):(?P<aa>[A-Za-z])(?P<pos>\d+)\*$")


def _parse_annotation_stop(mut: str) -> tuple[str, int] | None:
    s = mut.strip()

    m = _STOP_NEXTCLADE_RE.match(s)
    if m:
        return (m.group("gene"), int(m.group("pos")))

    m = _STOP_ANNO_GENERAL_RE.match(s)
    if m:
        return (m.group("gene"), int(m.group("pos")))

    return None


def _parse_nextclade_stop(item: str) -> tuple[str, int] | None:
    m = _STOP_NEXTCLADE_RE.match(item.strip())
    if not m:
        return None
    return (m.group("gene"), int(m.group("pos")))


def _build_stop_lookup(annotation_keys: set[str]) -> dict[tuple[str, int], list[str]]:
    out: dict[tuple[str, int], list[str]] = defaultdict(list)
    for key in annotation_keys:
        parsed = _parse_annotation_stop(key)
        if parsed is None:
            continue
        out[parsed].append(key)
    return dict(out)


# Annotation accepts: {gene}:{pos}-  OR  {gene}:{aa}{pos}-
_DEL_ANNO_RE = re.compile(r"^(?P<gene>[^:]+):(?P<aa>[A-Za-z])?(?P<pos>\d+)-$")

# Annotation convenience format: {gene}:del{start}-{stop} (inclusive)
# Permissive: allow whitespace after 'del' and around '-'
_DEL_RANGE_RE = re.compile(r"^(?P<gene>[^:]+):del\s*(?P<start>\d+)\s*-\s*(?P<stop>\d+)$")

# Nextclade accepts: {gene}:{aa}{pos}-
_DEL_NEXTCLADE_RE = re.compile(r"^(?P<gene>[^:]+):(?P<aa>[A-Za-z])(?P<pos>\d+)-$")


def _parse_annotation_del(mut: str) -> tuple[str, int] | None:
    m = _DEL_ANNO_RE.match(mut.strip())
    if not m:
        return None
    return (m.group("gene"), int(m.group("pos")))


def _parse_nextclade_del(item: str) -> tuple[str, int] | None:
    m = _DEL_NEXTCLADE_RE.match(item.strip())
    if not m:
        return None
    return (m.group("gene"), int(m.group("pos")))


def _build_del_lookup(annotation_keys: set[str]) -> dict[tuple[str, int], list[str]]:
    out: dict[tuple[str, int], list[str]] = defaultdict(list)
    for key in annotation_keys:
        parsed = _parse_annotation_del(key)
        if parsed is None:
            continue
        out[parsed].append(key)
    return dict(out)


def _split_linked_mutations_cell(cell: str | None) -> set[str]:
    """
    Parse the optional `linked_mutations` control column into a set of tokens.

    This column is intended to encode *co-occurrence requirements* for a hit:
    if a mutation's annotations include `linked_mutations`, then the current
    sequence must also contain every token listed in that field, otherwise the
    hit is suppressed.

    Input format
    ------------
    The cell may contain mutation tokens separated by commas and/or semicolons.
    Whitespace around tokens is ignored. Empty tokens are ignored.

    Examples
    --------
    - "S:E484K,S:K417N" -> {"S:E484K", "S:K417N"}
    - "S:E484K; S:K417N" -> {"S:E484K", "S:K417N"}
    - "" or None -> set()

    Notes
    -----
    This parser does not validate token formats; it treats them as opaque strings.
    Downstream matching uses exact membership in the detected mutation universe.
    """
    if cell is None:
        return set()
    s = cell.strip()
    if not s:
        return set()
    toks = re.split(r"[;,]", s)
    return {t.strip() for t in toks if t.strip()}


def _merge_linked_mutations(existing: str, incoming: str) -> str:
    """
    Merge two `linked_mutations` field values using set-union semantics.

    `linked_mutations` is a control field (not a human-facing annotation) and is
    best treated as a set of required tokens. When the same mutation appears in
    multiple rows, we union their required-token sets.

    Serialization
    -------------
    The merged set is serialized as a semicolon-separated string in sorted order to
    provide stable output.

    Examples
    --------
    existing="S:E484K", incoming="S:K417N" -> "S:E484K;S:K417N"
    existing="S:E484K;S:K417N", incoming="S:K417N" -> "S:E484K;S:K417N"
    """
    a = _split_linked_mutations_cell(existing)
    b = _split_linked_mutations_cell(incoming)
    merged = sorted(a | b)
    return ";".join(merged)


def _expand_deletion_range_key(raw: str) -> list[str] | None:
    """
    Expand an annotation-table deletion range key into canonical per-position deletions.

    This enables an input convenience syntax in the annotation table:

        "{gene}:del{start}-{stop}"

    where start and stop are 1-based amino-acid indices and the range is inclusive.

    The range is expanded ("rolled out") into annotation deletion keys in the
    existing canonical format accepted by downstream logic:

        "{gene}:{pos}-"

    Examples
    --------
    - "S:del87-89" -> ["S:87-", "S:88-", "S:89-"]
    - "ORF1a:del1-1" -> ["ORF1a:1-"]
    - "S:87-" -> None (not a range key; caller should treat as literal)

    Validation
    ----------
    - start and stop must be positive integers
    - stop must be >= start

    Returns
    -------
    list[str] | None
        A list of expanded mutation keys if `raw` is a range key, otherwise None.

    Raises
    ------
    InputFormatError
        If the key matches the range pattern but the indices are invalid.
    """
    s = raw.strip()
    m = _DEL_RANGE_RE.match(s)
    if not m:
        return None

    gene = m.group("gene")
    start = int(m.group("start"))
    stop = int(m.group("stop"))

    if start <= 0 or stop <= 0:
        raise InputFormatError(f"Deletion range indices must be positive: {raw!r}")
    if stop < start:
        raise InputFormatError(f"Deletion range stop < start: {raw!r}")

    return [f"{gene}:{i}-" for i in range(start, stop + 1)]


# Substitution (Nextclade style): {gene}:{ref}{pos}{alt}
# Example: S:N87Y
# Allow alt to be 'X' in annotation keys for wildcard matching.
_SUB_RE = re.compile(r"^(?P<gene>[^:]+):(?P<ref>[A-Za-z])(?P<pos>\d+)(?P<alt>[A-Za-z])$")

# Insertion (Nextclade style): {gene}:{pos}:{ins}
# Example: S:214:EPE
# Allow ins to contain one or more 'X' in annotation keys for wildcard matching.
_INS_RE = re.compile(r"^(?P<gene>[^:]+):(?P<pos>\d+):(?P<ins>[A-Za-z]+)$")


@dataclass(frozen=True, slots=True)
class _WildcardRule:
    """
    Compiled wildcard rule derived from an annotation-table mutation key containing 'X'.

    Attributes
    ----------
    pattern:
        Compiled regex matching Nextclade mutation tokens.
    anno_key:
        Original annotation-table mutation key this rule represents.
    """

    pattern: re.Pattern[str]
    anno_key: str


def _escape_except_x(s: str, *, x_charset: str) -> str:
    """
    Escape a literal string for regex, translating 'X' into a single-AA wildcard class.

    Parameters
    ----------
    s:
        Input string to translate.
    x_charset:
        Allowed characters for 'X'. Used to build a character class.

    Returns
    -------
    str
        Regex-safe string where literal content is escaped and 'X' is replaced with
        a single-character class.
    """
    # Escape everything, then un-escape the placeholder we insert for X.
    # We'll do this robustly by building the regex piecewise.
    aa_class = f"[{re.escape(x_charset)}]"
    parts: list[str] = []
    for ch in s:
        if ch == "X":
            parts.append(aa_class)
        else:
            parts.append(re.escape(ch))
    return "".join(parts)


def _compile_x_wildcard_rules(
    annotation_keys: set[str],
    *,
    x_charset: str,
) -> tuple[set[str], list[_WildcardRule]]:
    """
    Split annotation keys into exact keys and compiled 'X' wildcard rules.

    Only two mutation formats are eligible for wildcard compilation:
    - Substitution: {gene}:{ref}{pos}{alt} (alt may be 'X')
    - Insertion:   {gene}:{pos}:{ins}     (ins may contain 'X')

    Any annotation key containing 'X' that does not match either of these formats
    is treated as a literal exact key (i.e. no wildcard semantics applied).
    """
    exact: set[str] = set()
    rules: list[_WildcardRule] = []

    for key in annotation_keys:
        if "X" not in key:
            exact.add(key)
            continue

        s = key.strip()

        m = _SUB_RE.match(s)
        if m:
            # Only treat X as wildcard in the *alt* position (final AA).
            # If ref is X (weird), do not wildcard it.
            alt = m.group("alt")
            if alt != "X":
                exact.add(key)
                continue

            gene = m.group("gene")
            ref = m.group("ref")
            pos = m.group("pos")
            pat = re.compile(rf"^{re.escape(gene)}:{re.escape(ref)}{pos}[{re.escape(x_charset)}]$")
            rules.append(_WildcardRule(pattern=pat, anno_key=key))
            continue

        m = _INS_RE.match(s)
        if m:
            gene = m.group("gene")
            pos = m.group("pos")
            ins = m.group("ins")
            # Allow X anywhere in insertion string, each X == one AA.
            ins_pat = _escape_except_x(ins, x_charset=x_charset)
            pat = re.compile(rf"^{re.escape(gene)}:{pos}:{ins_pat}$")
            rules.append(_WildcardRule(pattern=pat, anno_key=key))
            continue

        # Unsupported 'X' placement/format: treat as literal exact key.
        exact.add(key)

    return exact, rules


def _merge_annotation_dicts_for_keys(
    annotations: Mapping[str, Mapping[str, str]],
    anno_keys: Sequence[str],
) -> dict[str, str]:
    """
    Merge annotation dictionaries for multiple matching annotation keys.

    Merge semantics:
    - Per field: merge non-empty values with `_merge_field_value`.
    - Missing keys are ignored.
    """
    merged: dict[str, str] = {}
    for k in anno_keys:
        d = annotations.get(k)
        if not d:
            continue
        for field, val in d.items():
            v = (val or "").strip()
            if not v:
                continue
            merged[field] = _merge_field_value(merged.get(field, ""), v)
    return merged
