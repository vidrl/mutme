"""
Example preset transformation functions for:
    - Stanford Database
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any

#: Canonical “full” names used for output columns.
MAB_ABBREVIATIONS: dict[str, str] = {
    "BAM": "Bamlanivimab",  # LY-CoV555 / LY3819253
    "ETE": "Etesevimab",  # LY-CoV016 / JS016 / CB6
    "CAS": "Casirivimab",  # REGN10933
    "IMD": "Imdevimab",  # REGN10987
    "CIL": "Cilgavimab",  # COV2-2130 / AZD1061
    "TIX": "Tixagevimab",  # COV2-2196 / AZD8895
    "SOT": "Sotrovimab",  # Vir-7831 / S309
    "BEB": "Bebtelovimab",  # LY-CoV1404 / LY3853113
    "REG": "Regdanvimab",  # CT-P59
    "AMU": "Amubarvimab",  # BRII-196 / P2C-1f11
    "ROM": "Romlusevimab",  # BRII-198 / P2B-1G5
    "ADI": "Adintrevimab",  # ADG20 / ADG-2
}

INH_ABBREVIATIONS: dict[str, str] = {
    "NTV": "Nirmatrelvir",  # PAXLovid / PF-07321332
    "ENS": "Ensitrelvir",  # Xocova / S-217622
}

RDRP_INH_ABBREVIATIONS: dict[str, str] = {
    "RDV": "Remdesivir",  # Veklury / GS-5734
}


@dataclass(frozen=True)
class MabColumns:
    """Resolved column names for one mAb in the source table."""

    abbrev: str
    full_name: str
    fold_col: str
    dms_col: str | None


_FOLD_RE = re.compile(r"^\s*([A-Z]{3})\s*:\s*fold\s*$")
_DMS_RE = re.compile(r"^\s*([A-Z]{3})\s*:\s*dms\s*$")


@dataclass(frozen=True)
class InhibitorColumns:
    """Resolved column names for one inhibitor in the source table."""

    abbrev: str
    full_name: str
    fold_col: str
    pocket_col: str | None


_INH_FOLD_RE = re.compile(r"^\s*([A-Z]{3})\s*:\s*fold\s*$")
_INH_POCKET_RE = re.compile(r"^\s*([A-Z]{3})\s*:\s*pocket\s*$")


@dataclass(frozen=True)
class RdRpInhibitorColumns:
    """Resolved column names for one RdRP inhibitor in the source table."""

    abbrev: str
    full_name: str
    fold_col: str


_RDRP_FOLD_RE = re.compile(r"^\s*([A-Z]{3})\s*:\s*fold\s*$")


def _parse_float(value: str | None) -> float | None:
    """
    Parse a numeric cell into a float.

    Returns None for empty/whitespace-only values. Raises ValueError if the
    value is non-empty but not parseable as float.
    """
    if value is None:
        return None
    s = value.strip()
    if not s:
        return None
    return float(s)


def classify_fold_susceptibility_reduction(fold: float) -> str:
    """
    Classify a fold susceptibility reduction value into a qualitative bucket.

    Parameters
    ----------
    fold:
        Fold susceptibility reduction (numeric).

    Returns
    -------
    str
        One of:
        - "high"   if fold >= 25
        - "medium" if 5 <= fold < 25
        - "low"    if fold < 5

    Notes
    -----
    This function assumes the input is a valid numeric fold value. If you want
    to treat missing values specially, handle that before calling.
    """
    if fold >= 25:
        return "high"
    if fold >= 5:
        return "medium"
    return "low"


def resolve_mab_columns(fieldnames: Iterable[str]) -> list[MabColumns]:
    """
    Resolve which mAbs are present in a SARS-CoV-2 spike mAb resistance CSV header.

    This inspects the header and identifies columns of the form:
      - "<ABBR>: fold"
      - "<ABBR>: dms"

    It then maps ABBR to a canonical full name (see MAB_ABBREVIATIONS) and returns
    a structured description for each detected mAb.

    Parameters
    ----------
    fieldnames:
        The CSV header fields (typically DictReader.fieldnames).

    Returns
    -------
    list[MabColumns]
        One entry per detected mAb that has a fold column. The DMS column may be
        absent (dms_col=None) if not present in the header.

    Raises
    ------
    KeyError
        If an abbreviation is encountered in the header but is not present in
        MAB_ABBREVIATIONS.

    Examples
    --------
    If the header contains "BAM: fold" and "BAM: dms", this returns an entry with:
      full_name="Bamlanivimab", fold_col="BAM: fold", dms_col="BAM: dms".
    """
    fields = list(fieldnames)
    fold_by_abbr: dict[str, str] = {}
    dms_by_abbr: dict[str, str] = {}

    for f in fields:
        m = _FOLD_RE.match(f)
        if m:
            fold_by_abbr[m.group(1)] = f
            continue
        m = _DMS_RE.match(f)
        if m:
            dms_by_abbr[m.group(1)] = f

    resolved: list[MabColumns] = []
    for abbr, fold_col in sorted(fold_by_abbr.items()):
        if abbr not in MAB_ABBREVIATIONS:
            raise KeyError(f"Unknown mAb abbreviation in header: {abbr!r}")
        resolved.append(
            MabColumns(
                abbrev=abbr,
                full_name=MAB_ABBREVIATIONS[abbr],
                fold_col=fold_col,
                dms_col=dms_by_abbr.get(abbr),
            )
        )
    return resolved


def transform_spike_mab_resistance_table(
    source: str | Path | IO[str],
    *,
    add_dms_plus: bool = False,
    dms_threshold: float = 0.1,
    mutation_prefix: str = "S:",
) -> list[dict[str, Any]]:
    """
    Transform a SARS-CoV-2 spike mAb resistance mutation table into a compact “per-drug” format.

    This is designed for CSV exports where the first column is "Mutation" and the remaining
    columns include per-mAb measurements using the naming convention:

      - "<ABBR>: fold"  (fold susceptibility reduction)
      - "<ABBR>: dms"   (DMS escape fraction; optional)

    The transformation produces one output row per input mutation, with:
      1) "Mutation" values prefixed (default "S:"), e.g. "R346K" -> "S:R346K"
      2) One output column per mAb (named by full canonical name, e.g. "Bamlanivimab")
         containing a qualitative summary based on the fold column:
            - "high"   if fold >= 25
            - "medium" if 5 <= fold < 25
            - "low"    if fold < 5
         If the fold cell is empty/missing, the output cell is "" (empty string).

      3) Optionally, if add_dms_plus=True, append "+" to the qualitative label when the
         corresponding DMS escape fraction is >= dms_threshold (default 0.1). If the DMS
         column is absent or empty, no "+" is appended.

    Parameters
    ----------
    source:
        Input CSV source. May be:
        - a filesystem path (str or Path)
        - an open text file handle (IO[str])
    add_dms_plus:
        If True, append "+" to "high"/"medium"/"low" when DMS escape fraction
        is >= dms_threshold for that mutation and mAb.
    dms_threshold:
        Threshold for DMS escape fraction that triggers "+" when add_dms_plus is True.
        Default is 0.1.
    mutation_prefix:
        Prefix prepended to each mutation value in the "Mutation" column.
        Default is "S:".

    Returns
    -------
    list[dict[str, Any]]
        A list of output rows (dicts). Keys are:
        - "Mutation"
        - One key per detected mAb using its full canonical name (e.g. "Etesevimab")

        Values are strings ("" or "low"/"medium"/"high" with optional "+").

    Raises
    ------
    ValueError
        If a non-empty fold or DMS cell cannot be parsed as float.
    KeyError
        If an abbreviation present in the header is not present in MAB_ABBREVIATIONS.
    RuntimeError
        If the CSV appears to be missing a "Mutation" header.

    Notes
    -----
    - This function is intentionally “lossy”: it bins numeric fold values into three buckets.
    - Any extra columns besides the recognized fold/dms columns are ignored.
    - If you need a stable output column order, you can post-process the dicts or
      wrap this to emit (fieldnames, rows).

    Examples
    --------
    Given an input row:
        Mutation="R346K", "BAM: fold"="1.3", "BAM: dms"="", "IMD: fold"="2.9", ...

    Output might include:
        {"Mutation": "S:R346K", "Bamlanivimab": "low", "Imdevimab": "low", ...}

    With add_dms_plus=True and a DMS value >= 0.1, you’d get e.g. "medium+".
    """
    # Open if a path was provided.
    must_close = False
    if isinstance(source, (str, Path)):
        fh = open(source, newline="", encoding="utf-8-sig")
        must_close = True
    else:
        fh = source

    try:
        reader = csv.DictReader(fh, delimiter=",")
        if not reader.fieldnames:
            raise RuntimeError("CSV appears to have no header row.")
        if "Mutation" not in reader.fieldnames:
            raise RuntimeError('CSV header is missing required column "Mutation".')

        mab_cols = resolve_mab_columns(reader.fieldnames)

        out_rows: list[dict[str, Any]] = []
        for row in reader:
            mut_raw = (row.get("Mutation") or "").strip()
            if not mut_raw:
                # Skip completely empty mutation rows.
                continue

            out: dict[str, Any] = {"Mutation": f"{mutation_prefix}{mut_raw}"}

            for mc in mab_cols:
                fold_val = _parse_float(row.get(mc.fold_col))
                if fold_val is None:
                    out[mc.full_name] = ""
                    continue

                label = classify_fold_susceptibility_reduction(fold_val)

                if add_dms_plus and mc.dms_col:
                    dms_val = _parse_float(row.get(mc.dms_col))
                    if dms_val is not None and dms_val >= dms_threshold:
                        label += "+"

                out[mc.full_name] = label

            out_rows.append(out)

        return out_rows

    finally:
        if must_close:
            fh.close()


def classify_3clpro_fold_susceptibility_reduction(fold: float) -> str:
    """
    Classify 3CLpro inhibitor fold susceptibility reduction according to the Stanford UI bins.

    Website bin definitions:
    - dark blue:   median >= 10-fold reduction
    - light blue:  median 5-10-fold reduction
    - very light:  median 2.5-5-fold reduction
    - white:       median < 2.5-fold reduction
    - gray:        absence of susceptibility data (handled upstream as missing)

    Parameters
    ----------
    fold:
        Fold susceptibility reduction (numeric).

    Returns
    -------
    str
        One of:
        - "high"       if fold >= 10
        - "medium"     if 5 <= fold < 10
        - "low"        if 2.5 <= fold < 5
        - "very_low"   if fold < 2.5
    """
    if fold >= 10:
        return "high"
    if fold >= 5:
        return "medium"
    if fold >= 2.5:
        return "low"
    return "very_low"


def _truthy_pocket(value: str | None) -> bool:
    """
    Interpret a 'pocket' cell as a boolean.

    The Stanford export often uses '1' to indicate pocket mutation, but we accept:
    - any non-empty string that parses to a non-zero float (e.g. '1', '1.0')
    - any non-empty non-numeric string (treated as True)

    Empty / whitespace -> False
    Numeric 0 / 0.0 -> False
    """
    if value is None:
        return False
    s = value.strip()
    if not s:
        return False
    try:
        return float(s) != 0.0
    except ValueError:
        return True


def resolve_inhibitor_columns(fieldnames: Iterable[str]) -> list[InhibitorColumns]:
    """
    Resolve which 3CLpro inhibitors are present in a resistance CSV header.

    Detects:
      - "<ABBR>: fold"
      - "<ABBR>: pocket" (optional)

    Uses INH_ABBREVIATIONS to map ABBR -> full output column name.

    Parameters
    ----------
    fieldnames:
        CSV header fields (typically DictReader.fieldnames).

    Returns
    -------
    list[InhibitorColumns]
        One entry per detected inhibitor that has a fold column. pocket_col may be None.

    Raises
    ------
    KeyError
        If an abbreviation is encountered in the header but is not present in INH_ABBREVIATIONS.
    """
    fields = list(fieldnames)
    fold_by_abbr: dict[str, str] = {}
    pocket_by_abbr: dict[str, str] = {}

    for f in fields:
        m = _INH_FOLD_RE.match(f)
        if m:
            fold_by_abbr[m.group(1)] = f
            continue
        m = _INH_POCKET_RE.match(f)
        if m:
            pocket_by_abbr[m.group(1)] = f

    resolved: list[InhibitorColumns] = []
    for abbr, fold_col in sorted(fold_by_abbr.items()):
        if abbr not in INH_ABBREVIATIONS:
            raise KeyError(f"Unknown inhibitor abbreviation in header: {abbr!r}")
        resolved.append(
            InhibitorColumns(
                abbrev=abbr,
                full_name=INH_ABBREVIATIONS[abbr],
                fold_col=fold_col,
                pocket_col=pocket_by_abbr.get(abbr),
            )
        )
    return resolved


def transform_3clpro_inhibitor_resistance_table(
    source: str | Path | IO[str],
    *,
    mutation_prefix: str = "nsp5:",
    add_pocket_suffix: bool = True,
    pocket_suffix: str = "-P",
) -> list[dict[str, Any]]:
    """
    Transform a 3CLpro inhibitor resistance table into a compact “per-drug” format.

    Input format assumptions
    ------------------------
    - First column is "Mutation"
    - Per-drug columns follow naming convention:
        * "<ABBR>: fold"
        * "<ABBR>: pocket" (optional)

    Output
    ------
    - "Mutation" value is prefixed (default "nsp5:")
    - One output column per inhibitor (full name from INH_ABBREVIATIONS)
      containing a qualitative label derived from the fold value:
        * "high"   if fold >= 10
        * "medium" if 5 <= fold < 10
        * "low"    if 2.5 <= fold < 5
        * "none"   if fold < 2.5
      If fold is missing/empty, output is "" (empty string).

    Pocket suffix
    -------------
    If add_pocket_suffix=True and the inhibitor has a "<ABBR>: pocket" column,
    then a trailing suffix (default "-P") is appended to the label when the pocket
    cell is truthy (e.g., "1").

    Examples:
      - fold=12.6, pocket empty -> "high"
      - fold=12.6, pocket=1     -> "high-P"
      - fold missing            -> ""

    Parameters
    ----------
    source:
        Input CSV source: path (str/Path) or open text file handle (IO[str]).
    mutation_prefix:
        Prefix prepended to mutation values, e.g. "nsp5:" or "3CLpro:".
    add_pocket_suffix:
        Whether to append the pocket suffix when pocket is truthy.
    pocket_suffix:
        Suffix appended when pocket is truthy (default "-P").

    Returns
    -------
    list[dict[str, Any]]
        Output rows as dictionaries.
    """

    must_close = False
    if isinstance(source, (str, Path)):
        fh = open(source, newline="", encoding="utf-8-sig")
        must_close = True
    else:
        fh = source

    try:
        reader = csv.DictReader(fh, delimiter=",")
        if not reader.fieldnames:
            raise RuntimeError("CSV appears to have no header row.")
        if "Mutation" not in reader.fieldnames:
            raise RuntimeError('CSV header is missing required column "Mutation".')
        else:
            mutation_col = "Mutation"

        inh_cols = resolve_inhibitor_columns(reader.fieldnames)

        out_rows: list[dict[str, Any]] = []
        for row in reader:
            mut_raw = (row.get(mutation_col) or "").strip()
            if not mut_raw:
                continue

            out: dict[str, Any] = {"Mutation": f"{mutation_prefix}{mut_raw}"}

            for ic in inh_cols:
                fold_cell = row.get(ic.fold_col)
                fold_val = _parse_float(fold_cell)
                if fold_val is None:
                    out[ic.full_name] = ""
                    continue

                label = classify_3clpro_fold_susceptibility_reduction(fold_val)

                if add_pocket_suffix and ic.pocket_col:
                    if _truthy_pocket(row.get(ic.pocket_col)):
                        label += pocket_suffix

                out[ic.full_name] = label

            out_rows.append(out)

        return out_rows

    finally:
        if must_close:
            fh.close()


def resolve_rdrp_inhibitor_columns(
    fieldnames: Iterable[str],
) -> list[RdRpInhibitorColumns]:
    """
    Resolve which RdRP inhibitors are present in a resistance CSV header.

    Detects columns named:
      - "<ABBR>: fold"

    Uses RDRP_INH_ABBREVIATIONS to map ABBR -> full output column name.

    Parameters
    ----------
    fieldnames:
        CSV header fields (typically DictReader.fieldnames).

    Returns
    -------
    list[RdRpInhibitorColumns]
        One entry per detected inhibitor that has a fold column.

    Raises
    ------
    KeyError
        If an abbreviation is encountered in the header but is not present in
        RDRP_INH_ABBREVIATIONS.
    """
    fold_by_abbr: dict[str, str] = {}

    for f in fieldnames:
        m = _RDRP_FOLD_RE.match(f)
        if m:
            fold_by_abbr[m.group(1)] = f

    resolved: list[RdRpInhibitorColumns] = []
    for abbr, fold_col in sorted(fold_by_abbr.items()):
        if abbr not in RDRP_INH_ABBREVIATIONS:
            raise KeyError(f"Unknown RdRP inhibitor abbreviation in header: {abbr!r}")
        resolved.append(
            RdRpInhibitorColumns(
                abbrev=abbr,
                full_name=RDRP_INH_ABBREVIATIONS[abbr],
                fold_col=fold_col,
            )
        )
    return resolved


def classify_rdrp_fold_susceptibility_reduction(fold: float) -> str:
    """
    Classify RdRP inhibitor fold susceptibility reduction using the Stanford UI bins.

    Bin definitions:
    - dark blue:   median >= 10-fold reduction  -> "high"
    - light blue:  median 5-10-fold reduction   -> "medium"
    - very light:  median 2.5-5-fold reduction  -> "low"
    - white:       median < 2.5-fold reduction  -> "very_low"
    - gray:        absence of susceptibility data (handled upstream as missing)

    Parameters
    ----------
    fold:
        Fold susceptibility reduction (numeric).

    Returns
    -------
    str
        One of: "high", "medium", "low", "none".
    """
    if fold >= 10:
        return "high"
    if fold >= 5:
        return "medium"
    if fold >= 2.5:
        return "low"
    return "very_low"


def transform_rdrp_inhibitor_resistance_table(
    source: str | Path | IO[str],
    *,
    mutation_prefix: str = "nsp12:",
) -> list[dict[str, Any]]:
    """
    Transform an RdRP (nsp12) inhibitor resistance table into a compact per-drug format.

    Input format assumptions
    ------------------------
    - First column is "Mutation" (case-insensitive: accepts "Mutation" or "mutation")
    - Per-drug columns follow naming convention:
        * "<ABBR>: fold"

    Output
    ------
    - "Mutation" value is prefixed (default "nsp12:")
    - One output column per inhibitor (full name from RDRP_INH_ABBREVIATIONS)
      containing a qualitative label derived from fold:
        * "high"   if fold >= 10
        * "medium" if 5 <= fold < 10
        * "low"    if 2.5 <= fold < 5
        * "none"   if fold < 2.5
      If fold is missing/empty, output is "".

    Parameters
    ----------
    source:
        Input CSV source: path (str/Path) or open text file handle (IO[str]).
    mutation_prefix:
        Prefix prepended to mutation values (default "nsp12:").

    Returns
    -------
    list[dict[str, Any]]
        Output rows as dictionaries.

    Raises
    ------
    RuntimeError
        If the CSV appears to have no header or no mutation column.
    ValueError
        If a non-empty fold cell cannot be parsed as float.
    KeyError
        If a header abbreviation is unknown.
    """
    must_close = False
    if isinstance(source, (str, Path)):
        fh = open(source, newline="", encoding="utf-8-sig")
        must_close = True
    else:
        fh = source

    try:
        reader = csv.DictReader(fh, delimiter=",")
        if not reader.fieldnames:
            raise RuntimeError("CSV appears to have no header row.")

        if "Mutation" in reader.fieldnames:
            mutation_col = "Mutation"
        else:
            raise RuntimeError('CSV header is missing required column "Mutation" (or "mutation").')

        drug_cols = resolve_rdrp_inhibitor_columns(reader.fieldnames)

        out_rows: list[dict[str, Any]] = []
        for row in reader:
            mut_raw = (row.get(mutation_col) or "").strip()
            if not mut_raw:
                continue

            out: dict[str, Any] = {"Mutation": f"{mutation_prefix}{mut_raw}"}

            for dc in drug_cols:
                fold_val = _parse_float(row.get(dc.fold_col))
                if fold_val is None:
                    out[dc.full_name] = ""
                    continue
                out[dc.full_name] = classify_rdrp_fold_susceptibility_reduction(fold_val)

            out_rows.append(out)

        return out_rows

    finally:
        if must_close:
            fh.close()
