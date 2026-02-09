from __future__ import annotations

import csv
import os
import shlex
import subprocess
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any, Literal

DEFAULT_MUTATION_COLUMNS = (
    "seqName",
    "qc.overallStatus",
    "qc.stopCodons.stopCodons",
    "substitutions",
    "deletions",
    "insertions",
    "aaSubstitutions",
    "aaDeletions",
    "aaInsertions",
)


DatabasePreset = Literal[
    "stanford-scov2-mab-resistance",
    "stanford-scov2-3clpro-inhibitor",
    "stanford-scov2-rdrp-inhibitor",
]
AlignmentPreset = Literal["default", "high-diversity", "short-sequences"]
QualityControlStatus = Literal["good", "mediocre", "bad"]


class CommandExecutionError(RuntimeError):
    """
    Raised when a system command returns a non-zero exit code (or fails to start).

    This is a *rich* error that carries the execution result so callers can log or
    inspect stderr/stdout without re-running anything.
    """

    def __init__(self, message: str, result: CommandResult | None = None) -> None:
        super().__init__(message)
        self.result = result


@dataclass(frozen=True, slots=True)
class CommandResult:
    """
    Structured outcome of a command execution.

    Attributes
    ----------
    args:
        The command executed, normalized as a list of strings (argv-style).
    returncode:
        Process exit status. Conventionally, 0 means success.
    stdout:
        Captured standard output as text (empty string if not captured / no output).
    stderr:
        Captured standard error as text (empty string if not captured / no output).
    cwd:
        Working directory used to run the command, if any.
    """

    args: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    cwd: str | None


def run_command(
    command: str | Sequence[str],
    *,
    cwd: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
    timeout: float | None = None,
    check: bool = True,
    text: bool = True,
    capture_output: bool = True,
    stdin: str | bytes | None = None,
    shell: bool = False,
) -> CommandResult:
    """
    Execute a system command in a production-friendly way.

    This function is a safer, more ergonomic wrapper around :func:`subprocess.run`.
    It returns a structured :class:`CommandResult` and (optionally) raises a
    :class:`CommandExecutionError` on failure.

    The default behavior is designed for typical backend/service use:
    - **No shell** execution (reduces injection risk)
    - **Captures stdout/stderr** for logging and debugging
    - **Text mode** decoding (UTF-8 using the current locale rules)
    - **Raises on non-zero** exit codes (configurable)

    Parameters
    ----------
    command:
        The command to execute. Use either:
        - a sequence of args: ``["git", "status", "--porcelain"]`` (preferred), or
        - a string: ``"git status --porcelain"``

        If a string is passed **and** ``shell=False`` (default), it is split using
        :func:`shlex.split` (POSIX-style parsing).

        If ``shell=True``, the command is passed through to the platform shell
        (e.g. ``/bin/sh`` or ``cmd.exe``). This is powerful but dangerous—avoid it
        unless you truly need shell features (pipes, redirects, globs).

    cwd:
        Working directory for the command. If provided, it is resolved to a string
        path and validated to exist.

    env:
        Environment variables for the command. If provided, these values overlay
        the current process environment.

    timeout:
        Maximum number of seconds to allow the process to run. If exceeded, a
        :class:`CommandExecutionError` is raised.

    check:
        If True (default), raise :class:`CommandExecutionError` when the command
        exits with a non-zero status.

    text:
        If True (default), operate in text mode and return stdout/stderr as strings.
        If False, stdout/stderr are treated as bytes. Note: this function always
        returns strings in :class:`CommandResult`; if you need bytes, set
        ``text=False`` and pass ``capture_output=True`` then decode yourself (or
        adapt the dataclass).

    capture_output:
        If True (default), capture stdout/stderr and return them in the result.
        If False, subprocess inherits parent stdout/stderr, and returned stdout/stderr
        will be empty strings.

    stdin:
        Data to send to stdin. Provide a string in text mode or bytes in binary mode.
        If None, stdin is not provided.

    shell:
        If True, execute via the system shell. **Use with extreme caution**.
        If False (default), execute directly without a shell.

    Returns
    -------
    CommandResult
        A structured result with stdout/stderr (when captured) and the exit code.

    Raises
    ------
    ValueError
        If ``cwd`` is provided but does not exist or is not a directory.
    TypeError
        If ``command`` is not a string or sequence of strings.
    CommandExecutionError
        If the command fails to start, times out, or returns non-zero with ``check=True``.

    Security Notes
    --------------
    - Prefer passing ``command`` as a list/tuple of args with ``shell=False``.
    - Avoid interpolating untrusted input into command strings.

    Examples
    --------
    Run a command (recommended argv style):

    >>> result = run_command(["python", "--version"])
    >>> result.returncode
    0

    Run a command with stdin:

    >>> run_command(["bash", "-lc", "cat"], stdin="hello").stdout
    'hello'

    Allow non-zero exit codes:

    >>> result = run_command(["bash", "-lc", "exit 2"], check=False)
    >>> result.returncode
    2
    """

    # --- Normalize and validate command ---
    if isinstance(command, str):
        if shell:
            args: list[str] | str = command
            args_tuple = (command,)
        else:
            parts = shlex.split(command)
            args = parts
            args_tuple = tuple(parts)
    elif isinstance(command, Sequence) and not isinstance(command, (bytes, bytearray)):
        if not all(isinstance(x, str) for x in command):
            raise TypeError("command must be a str or a sequence of str")
        args = list(command)
        args_tuple = tuple(args)
    else:
        raise TypeError("command must be a str or a sequence of str")

    # --- cwd validation ---
    cwd_str: str | None = None
    if cwd is not None:
        p = Path(cwd)
        if not p.exists():
            raise ValueError(f"cwd does not exist: {p}")
        if not p.is_dir():
            raise ValueError(f"cwd is not a directory: {p}")
        cwd_str = str(p)

    # --- env overlay ---
    env_final: dict[str, str] | None
    if env is None:
        env_final = None
    else:
        env_final = dict(os.environ)
        env_final.update({str(k): str(v) for k, v in env.items()})

    # --- stdin handling ---
    input_data: str | bytes | None = stdin
    if input_data is not None:
        if text and isinstance(input_data, bytes):
            raise TypeError("stdin is bytes but text=True; pass a str or set text=False")
        if (not text) and isinstance(input_data, str):
            raise TypeError("stdin is str but text=False; pass bytes or set text=True")

    try:
        completed = subprocess.run(
            args,  # type: ignore[arg-type]
            cwd=cwd_str,
            env=env_final,
            input=input_data,
            timeout=timeout,
            check=False,  # we implement our own check to raise a richer error
            text=text,
            capture_output=capture_output,
            shell=shell,
        )
    except FileNotFoundError as e:
        raise CommandExecutionError(
            f"Command not found: {args_tuple[0] if args_tuple else command!r}"
        ) from e
    except subprocess.TimeoutExpired as e:
        # TimeoutExpired may carry partial stdout/stderr depending on platform.
        stdout = ""
        stderr = ""
        if getattr(e, "stdout", None) is not None:
            stdout = (
                e.stdout.decode() if isinstance(e.stdout, (bytes, bytearray)) else str(e.stdout)
            )
        if getattr(e, "stderr", None) is not None:
            stderr = (
                e.stderr.decode() if isinstance(e.stderr, (bytes, bytearray)) else str(e.stderr)
            )

        partial = CommandResult(
            args=args_tuple,
            returncode=-1,
            stdout=stdout,
            stderr=stderr,
            cwd=cwd_str,
        )
        raise CommandExecutionError(
            f"Command timed out after {timeout} seconds: {args_tuple!r}", result=partial
        ) from e
    except OSError as e:
        raise CommandExecutionError(f"Failed to execute command: {args_tuple!r}") from e

    result = CommandResult(
        args=args_tuple,
        returncode=completed.returncode,
        stdout=completed.stdout if (capture_output and completed.stdout is not None) else "",
        stderr=completed.stderr if (capture_output and completed.stderr is not None) else "",
        cwd=cwd_str,
    )

    if check and result.returncode != 0:
        # Don't dump stdout/stderr into the message (logs can decide); keep message tight.
        raise CommandExecutionError(
            f"Command failed with exit code {result.returncode}: {result.args!r}",
            result=result,
        )

    return result


@dataclass(frozen=True, slots=True)
class NextcladeTabularOutput:
    """Paths to tabular outputs produced by Nextclade."""

    tsv: Path


def run_nextclade(
    sequences_fasta: str | os.PathLike[str],
    reference_fasta: str | os.PathLike[str],
    annotation_gff3: str | os.PathLike[str],
    out_tsv: str | os.PathLike[str],
    alignment_preset: AlignmentPreset = "default",
    nextclade_bin: str = "nextclade",
    # Optional but sometimes useful even without a full dataset:
    extra_args: Sequence[str] = (),
    env: Mapping[str, str] | None = None,
    timeout_s: float | None = None,
    threads: int | None = None,
) -> tuple[NextcladeTabularOutput, CommandResult]:
    """
    Run Nextclade (v3.x) *without a dataset*, using only reference FASTA + GFF3 annotation,
    and write a TSV suitable for mutation extraction.

    Parameters
    ----------
    sequences_fasta:
        Query sequences to analyze (FASTA).
    reference_fasta:
        Reference genome (FASTA) used for alignment and mutation calling.
    annotation_gff3:
        Genome annotation in GFF3 format. Required for AA mutation calling.
    out_tsv:
        Output TSV path written by Nextclade.
    alignment_preset:
        One of: "default", "high-diversity", "short-sequences".
    nextclade_bin:
        Executable name or absolute path to `nextclade`.
    extra_args:
        Additional CLI args to pass through (advanced; treated as trusted input).
    env, timeout_s:
        Process environment overlay and timeout.

    Returns
    -------
    (NextcladeTabularOutput, CommandResult)

    Raises
    ------
    FileNotFoundError:
        If required inputs are missing.
    ValueError:
        If `alignment_preset` is invalid.
    CommandExecutionError:
        If Nextclade fails or output is missing.
    """
    seq = Path(sequences_fasta)
    ref = Path(reference_fasta)
    gff = Path(annotation_gff3)
    out = Path(out_tsv)

    for p, label in [
        (seq, "sequences_fasta"),
        (ref, "reference_fasta"),
        (gff, "annotation_gff3"),
    ]:
        if not p.exists():
            raise FileNotFoundError(f"{label} not found: {p}")
        if not p.is_file():
            raise FileNotFoundError(f"{label} is not a file: {p}")

    if alignment_preset not in ("default", "high-diversity", "short-sequences"):
        raise ValueError(f"Invalid alignment_preset: {alignment_preset!r}")

    out.parent.mkdir(parents=True, exist_ok=True)

    args: list[str] = [
        nextclade_bin,
        "run",
        "--input-ref",
        str(ref),
        "--input-annotation",
        str(gff),
        "--alignment-preset",
        alignment_preset,
        "--output-tsv",
        str(out),
        str(seq),
    ]

    # Nextclade expects a comma-separated list.
    args.extend(["--output-columns-selection", ",".join(DEFAULT_MUTATION_COLUMNS)])

    # Jobs specification
    if threads is not None:
        args.extend(["--jobs", str(threads)])

    args.extend(list(extra_args))

    result = run_command(
        args,
        env=env,
        timeout=timeout_s,
        check=True,
        capture_output=True,
        text=True,
        shell=False,
    )

    if not out.exists():
        raise CommandExecutionError(
            f"Nextclade reported success but TSV output is missing: {out}",
            result=result,
        )

    return NextcladeTabularOutput(tsv=out), result


def write_rows_to_delimited_file(
    rows: Iterable[Mapping[str, Any]],
    destination: str | Path,
    *,
    delimiter: str,
    include_header: bool = True,
) -> None:
    """
    Write row-oriented data to a delimited text file (CSV/TSV).

    This is a generic helper that takes an iterable of mapping objects
    (e.g. dicts) and writes them to disk using the specified delimiter.
    Column order is inferred from the first row.

    Parameters
    ----------
    rows:
        Iterable of row mappings. Each mapping represents one output row,
        where keys are column names and values are cell contents.
        All rows are expected to share the same keys.
    destination:
        Output file path.
    delimiter:
        Field delimiter to use (e.g. "," for CSV, "\\t" for TSV).
    include_header:
        Whether to write the header row. Defaults to True.

    Raises
    ------
    ValueError
        If `rows` is empty.
    IOError
        If the file cannot be written.

    Notes
    -----
    - Values are written exactly as provided (converted to str by csv).
    - Newlines are handled according to csv module best practices.
    """
    rows = list(rows)
    if not rows:
        raise ValueError("No rows provided; refusing to write an empty file.")

    fieldnames = list(rows[0].keys())

    with open(destination, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter=delimiter,
            extrasaction="ignore",
        )
        if include_header:
            writer.writeheader()
        writer.writerows(rows)


def write_rows_to_csv(
    rows: Iterable[Mapping[str, Any]],
    destination: str | Path,
    *,
    include_header: bool = True,
) -> None:
    """
    Write transformed spike mAb resistance data to a CSV file.

    This is a thin wrapper around `write_rows_to_delimited_file`
    with delimiter set to comma.

    Parameters
    ----------
    rows:
        Iterable of row mappings produced by
        `transform_spike_mab_resistance_csv`.
    destination:
        Output CSV file path.
    include_header:
        Whether to write the header row. Defaults to True.
    """
    write_rows_to_delimited_file(
        rows,
        destination,
        delimiter=",",
        include_header=include_header,
    )


def write_rows_to_tsv(
    rows: Iterable[Mapping[str, Any]],
    destination: str | Path,
    *,
    include_header: bool = True,
) -> None:
    """
    Write transformed spike mAb resistance data to a TSV file.

    This is a thin wrapper around `write_rows_to_delimited_file`
    with delimiter set to tab.

    Parameters
    ----------
    rows:
        Iterable of row mappings produced by
        `transform_spike_mab_resistance_csv`.
    destination:
        Output TSV file path.
    include_header:
        Whether to write the header row. Defaults to True.
    """
    write_rows_to_delimited_file(
        rows,
        destination,
        delimiter="\t",
        include_header=include_header,
    )


def write_command_log(
    result: CommandResult,
    log_path: str | Path,
    *,
    include_stdout: bool = True,
    include_stderr: bool = True,
) -> None:
    """
    Write a detailed command execution log to disk.

    Parameters
    ----------
    result:
        CommandResult returned by run_command.
    log_path:
        Destination path for the log file.
    include_stdout:
        Whether to include stdout in the log.
    include_stderr:
        Whether to include stderr in the log.

    Notes
    -----
    - The log format is plain text and stable.
    - This function never raises due to command failure; it only raises
      on I/O errors.
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w", encoding="utf-8") as fh:
        fh.write("=== COMMAND EXECUTION LOG ===\n\n")
        fh.write(f"Command      : {' '.join(result.args)}\n")
        fh.write(f"Return code  : {result.returncode}\n")
        fh.write(f"Working dir : {result.cwd or '<inherit>'}\n\n")

        if include_stdout:
            fh.write("=== STDOUT ===\n")
            fh.write(result.stdout or "<empty>\n")
            fh.write("\n")

        if include_stderr:
            fh.write("=== STDERR ===\n")
            fh.write(result.stderr or "<empty>\n")
            fh.write("\n")


def log_command_result(
    result_or_error: CommandResult | CommandExecutionError,
    log_path: str | Path,
) -> None:
    """
    Write a command execution log from either a CommandResult or
    a CommandExecutionError.

    Parameters
    ----------
    result_or_error:
        Either:
        - CommandResult (successful or non-checked failure), or
        - CommandExecutionError carrying a result.
    log_path:
        Destination log file path.

    Notes
    -----
    - If the error does not contain a result (e.g. command not found),
      a minimal log is written.
    """
    if isinstance(result_or_error, CommandExecutionError):
        result = result_or_error.result
        if result is None:
            # Minimal log for failures before execution
            Path(log_path).write_text(
                f"Command failed before execution:\n{result_or_error}\n",
                encoding="utf-8",
            )
            return
    else:
        result = result_or_error

    write_command_log(result, log_path)


@dataclass(frozen=True, slots=True)
class CleanupReport:
    """
    Outcome of a cleanup attempt.

    Attributes
    ----------
    path:
        Path that was targeted for cleanup.
    removed:
        True if the file was deleted.
    reason:
        Human-readable explanation (e.g. 'deleted', 'missing', 'kept', 'permission denied').
    """

    path: Path
    removed: bool
    reason: str


def cleanup_file(
    path: str | os.PathLike[str],
    *,
    keep: bool = False,
    missing_ok: bool = True,
) -> CleanupReport:
    """
    Delete a file if it exists, optionally keeping it.

    Parameters
    ----------
    path:
        File path to remove.
    keep:
        If True, do not delete; return a report stating it was kept.
    missing_ok:
        If True, a missing file is not treated as an error.

    Returns
    -------
    CleanupReport
        Indicates whether deletion occurred and why.

    Notes
    -----
    - This function never raises for missing files when missing_ok=True.
    - Permission errors or other OS errors are captured into the report.
    """
    p = Path(path)

    if keep:
        return CleanupReport(path=p, removed=False, reason="kept")

    if not p.exists():
        if missing_ok:
            return CleanupReport(path=p, removed=False, reason="missing")
        return CleanupReport(path=p, removed=False, reason="missing (missing_ok=False)")

    if not p.is_file():
        return CleanupReport(path=p, removed=False, reason="not a file")

    try:
        p.unlink()
        return CleanupReport(path=p, removed=True, reason="deleted")
    except PermissionError as e:
        return CleanupReport(path=p, removed=False, reason=f"permission denied: {e}")
    except OSError as e:
        return CleanupReport(path=p, removed=False, reason=f"os error: {e}")


@dataclass(frozen=True)
class SubsetGffResult:
    """
    Result summary for `subset_gff3_by_mutation_prefixes`.
    """

    unique_prefixes: set[str]
    matched_gene_ids: set[str]
    genes_written: int
    cds_written: int
    rows_read: int
    warnings_emitted: int


def subset_gff3_by_mutation_prefixes(
    gff3_path: str | os.PathLike,
    table_path: str | os.PathLike,
    delimiter: str = ",",
    mutation_column: str = "mutation",
    prefix_separator: str = ":",
    strip_prefixes: bool = True,
    gff_gene_feature_type: str = "gene",
    gene_match_attr_key: str = "gene",
    include_cds: bool = False,
    cds_feature_type: str = "CDS",
    cds_gene_attr_key: str = "gene",
    output_path: str | os.PathLike | None = None,
    encoding: str = "utf-8",
) -> SubsetGffResult:
    """
    Subset a GFF3 file using unique gene-name prefixes derived from a table column.

    This helper reads a CSV/TSV-like table (configurable delimiter), locates a column
    named "mutation" case-insensitively (or a user-specified column name), and extracts
    *prefixes* from each cell by splitting on `prefix_separator` and taking the portion
    before the first separator. If a value does not contain the separator, a warning is
    emitted and the value is skipped. Prefixes are de-duplicated.

    The GFF3 is then streamed and filtered to keep:

      1) All `gff_gene_feature_type` features (default: "gene") whose attributes contain
         a value under any key in `gene_match_attr_keys` (default: ("Name", "gene"))
         matching one of the unique prefixes.

      2) Optionally, all `cds_feature_type` features (default: "CDS") whose attribute
         `cds_gene_attr_key` (default: "gene") matches one of the unique prefixes.

    Output is written as a valid GFF3 subset (preserving header/comment lines). If
    `output_path` is None, the function returns results without writing. If `output_path`
    is provided, the subset is written to that file.

    Parameters
    ----------
    gff3_path:
        Path to an input GFF3 file. Only plain text.
    table_path:
        Path to an input delimited table (CSV/TSV/etc.). Only plain text.
    delimiter:
        The table delimiter, e.g. "," for CSV or "\\t" for TSV.
    mutation_column:
        The name of the table column containing mutation strings. The column lookup is
        case-insensitive.
    prefix_separator:
        Separator used to split mutation strings, e.g. ":" for values like "GENE:K123N".
        The prefix is the substring before the first occurrence of this separator.
    strip_prefixes:
        If True, whitespace is stripped around extracted prefixes (and around whole
        values when no separator is present).
    gff_gene_feature_type:
        Feature type to treat as the "gene" record in column 3 of GFF3 (default "gene").
    gene_match_attr_keys:
        Attribute keys in the GFF3 attribute column (9th field) to try when matching gene
        features. Common choices include "Name", "gene", "ID", "locus_tag".
    include_cds:
        If True, also include CDS records whose attribute `cds_gene_attr_key` matches a
        unique prefix.
    cds_feature_type:
        Feature type to include for coding sequence features (default "CDS").
    cds_gene_attr_key:
        Attribute key to use when matching CDS features (default "gene").
    output_path:
        Path to write the subset GFF3. If None, no file is written and only summary
        statistics are returned.
    encoding:
        Text encoding for reading/writing.

    Returns
    -------
    SubsetGffResult
        Summary containing the set of unique prefixes, matched gene IDs, record counts,
        and warnings emitted.

    Raises
    ------
    FileNotFoundError
        If either input path does not exist.
    ValueError
        If the mutation column cannot be found (case-insensitive) in the table header,
        or if delimiter is invalid, or if GFF3 lines are malformed.
    OSError
        If an I/O error occurs reading or writing files.

    Notes
    -----
    - The function is streaming: it does not load the entire GFF3 into memory.

    - Attribute parsing follows the GFF3 convention of semicolon-separated key=value
      fields in the 9th column. Values are URL-decoded only minimally (no full percent
      decoding) to avoid surprising transformations; matching is performed on the raw
      value string after stripping surrounding whitespace.

    - Matching is exact, case-sensitive by default (consistent with typical gene IDs).
      If you need case-insensitive matching, normalize prefixes and attributes upstream
      (e.g. .upper()).

    Examples
    --------
    >>> res = subset_gff3_by_mutation_prefixes(
    ...     "annotations.gff3",
    ...     "mutations.tsv",
    ...     delimiter="\\t",
    ...     include_cds=True,
    ...     output_path="subset.gff3",
    ... )
    >>> sorted(list(res.unique_prefixes))[:5]
    """

    gff3_path = Path(gff3_path)
    table_path = Path(table_path)

    if not gff3_path.exists():
        raise FileNotFoundError(f"GFF3 not found: {gff3_path}")
    if not table_path.exists():
        raise FileNotFoundError(f"Table not found: {table_path}")
    if not delimiter or len(delimiter) != 1:
        raise ValueError(f"Delimiter must be a single character; got {delimiter!r}")

    warnings_emitted = 0

    def _open_text(path: Path, mode: str) -> IO[Any]:
        """
        Open plain text files transparently, returning a text handle.
        """
        if "b" in mode:
            raise ValueError("Binary mode not supported by _open_text; use text modes only.")
        return path.open(mode=mode, encoding=encoding, newline="")

    def _parse_gff3_attrs(attr_field: str) -> dict[str, str]:
        """
        Parse the 9th GFF3 column into a dict of key -> value.
        """
        attrs: dict[str, str] = {}
        s = attr_field.strip()
        if not s or s == ".":
            return attrs
        # GFF3: key=value;key2=value2
        parts = s.split(";")
        for part in parts:
            part = part.strip()
            if not part:
                continue
            if "=" not in part:
                # Non-conforming attribute token; ignore but keep predictable behavior.
                continue
            k, v = part.split("=", 1)
            k = k.strip()
            v = v.strip()
            if k:
                attrs[k] = v
        return attrs

    # Read table and extract unique prefixes
    unique_prefixes: set[str] = set()
    rows_read = 0

    with _open_text(table_path, "rt") as tf:
        reader = csv.reader(tf, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"Table appears empty: {table_path}")

        # Find mutation column case-insensitively
        target = mutation_column.strip().lower()
        col_idx = None
        for i, col in enumerate(header):
            if col.strip().lower() == target:
                col_idx = i
                break
        if col_idx is None:
            raise ValueError(
                f"Could not find a column named {mutation_column!r} (case-insensitive) "
                f"in table header: {header!r}"
            )

        for row in reader:
            rows_read += 1
            if col_idx >= len(row):
                # Ragged row: treat as missing
                continue
            raw = row[col_idx]
            if raw is None:
                continue
            val = raw.strip() if strip_prefixes else str(raw)
            if val == "":
                continue

            if prefix_separator in val:
                prefix = val.split(prefix_separator, 1)[0]
                prefix = prefix.strip() if strip_prefixes else prefix
            else:
                warnings_emitted += 1
                print(
                    f"Mutation value {val} did not contain separator {prefix_separator}; skipping entry."
                )
                continue

            if prefix:
                unique_prefixes.add(prefix)

    # Early exit: nothing to match
    if not unique_prefixes:
        if output_path is not None:
            # Still write header-only subset? Here we write only GFF header/comments if present.
            outp = Path(output_path)
            with _open_text(gff3_path, "rt") as gi, _open_text(outp, "wt") as go:
                for line in gi:
                    if line.startswith("#"):
                        go.write(line)
                    else:
                        break
        return SubsetGffResult(
            unique_prefixes=set(),
            matched_gene_ids=set(),
            genes_written=0,
            cds_written=0,
            rows_read=rows_read,
            warnings_emitted=warnings_emitted,
        )

    # Stream GFF3 and filter
    genes_written = 0
    cds_written = 0
    matched_gene_ids: set[str] = set()

    out_handle: IO[Any] | None = None
    if output_path is not None:
        outp = Path(output_path)
        out_handle = _open_text(outp, "wt")

    try:
        with _open_text(gff3_path, "rt") as gi:
            for line in gi:
                if line.startswith("#"):
                    if out_handle is not None:
                        out_handle.write(line)
                    continue

                stripped = line.rstrip("\n")
                if not stripped:
                    continue
                fields = stripped.split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"Malformed GFF3 line (expected 9 tab-separated fields): {stripped!r}"
                    )

                ftype = fields[2]
                attrs = _parse_gff3_attrs(fields[8])

                if ftype == gff_gene_feature_type:
                    keep, gid = _should_keep_gene(
                        attrs,
                        unique_prefixes,
                        gene_attr_key=gene_match_attr_key,
                    )
                    if keep:
                        genes_written += 1
                        if gid is not None:
                            matched_gene_ids.add(gid)
                        if out_handle is not None:
                            out_handle.write(line)
                elif include_cds and ftype == cds_feature_type:
                    keep, gid = _should_keep_cds(
                        attrs,
                        unique_prefixes,
                        cds_gene_attr_key=cds_gene_attr_key,
                    )
                    if keep:
                        # If both "Name" and "gene" attributes exist, force "Name" value to equal "gene" value
                        # otherwise Nextclade will use the "Name" value directly and it may not correspond to
                        # the gene prefixes in the annotation table
                        gene_val = attrs.get("gene")
                        if gene_val is not None and "Name" in attrs:
                            attrs["Name"] = gene_val.strip()

                        cds_written += 1
                        if gid is not None:
                            matched_gene_ids.add(gid)

                        if out_handle is not None:
                            # Write a modified line with updated attributes
                            fields[8] = _format_gff3_attrs(attrs)
                            out_handle.write("\t".join(fields) + "\n")
                else:
                    # Ignore other features
                    continue
    finally:
        if out_handle is not None:
            out_handle.close()

    return SubsetGffResult(
        unique_prefixes=unique_prefixes,
        matched_gene_ids=matched_gene_ids,
        genes_written=genes_written,
        cds_written=cds_written,
        rows_read=rows_read,
        warnings_emitted=warnings_emitted,
    )


def _should_keep_gene(
    attrs: Mapping[str, str],
    unique_prefixes: set[str],
    gene_attr_key: str,
) -> tuple[bool, str | None]:
    """
    Decide whether a gene feature should be kept based on a single attribute key.

    Returns (keep, matched_value).
    """
    v = attrs.get(gene_attr_key)
    if v is None:
        return False, None

    v = v.strip()
    if v in unique_prefixes:
        return True, v

    return False, None


def _should_keep_cds(
    attrs: Mapping[str, str],
    unique_prefixes: set[str],
    cds_gene_attr_key: str,
) -> tuple[bool, str | None]:
    """
    Decide whether a CDS feature should be kept based on a single attribute key.

    Returns (keep, matched_value).
    """
    v = attrs.get(cds_gene_attr_key)
    if v is None:
        return False, None

    v = v.strip()
    if v in unique_prefixes:
        return True, v

    return False, None


def _format_gff3_attrs(attrs: Mapping[str, str]) -> str:
    """
    Serialize a dict of GFF3 attributes into the 9th column format key=value;key2=value2.
    """
    if not attrs:
        return "."
    # Prefer a stable, conventional order: ID first if present, then the rest alphabetically.
    items = list(attrs.items())
    if "ID" in attrs:
        items = [("ID", attrs["ID"])] + [(k, v) for (k, v) in items if k != "ID"]
    # Sort remaining keys (excluding ID which is already first)
    head = items[:1] if items and items[0][0] == "ID" else []
    tail = items[1:] if head else items
    tail = sorted(tail, key=lambda kv: kv[0])
    items = head + tail
    return ";".join(f"{k}={v}" for k, v in items)
