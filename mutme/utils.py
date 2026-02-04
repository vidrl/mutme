from __future__ import annotations

import os
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence, Literal

DEFAULT_MUTATION_COLUMNS = (
    "seqName",
    "qc.overallStatus",
    "substitutions",
    "deletions",
    "insertions",
    "aaSubstitutions",
    "aaDeletions",
    "aaInsertions",
)



AlignmentPreset = Literal["default", "high-diversity", "short-sequences"]



QualityControlStatus = Literal["good", "mediocre", "bad"]

class CommandExecutionError(RuntimeError):
    """
    Raised when a system command returns a non-zero exit code (or fails to start).

    This is a *rich* error that carries the execution result so callers can log or
    inspect stderr/stdout without re-running anything.
    """

    def __init__(self, message: str, result: "CommandResult | None" = None) -> None:
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
            stdout = e.stdout.decode() if isinstance(e.stdout, (bytes, bytearray)) else str(e.stdout)
        if getattr(e, "stderr", None) is not None:
            stderr = e.stderr.decode() if isinstance(e.stderr, (bytes, bytearray)) else str(e.stderr)

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
    *,
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
) -> tuple[NextcladeTabularOutput, "CommandResult"]:
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

    for p, label in [(seq, "sequences_fasta"), (ref, "reference_fasta"), (gff, "annotation_gff3")]:
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
        "--input-ref", str(ref),
        "--input-annotation", str(gff),
        "--alignment-preset", alignment_preset,
        "--output-tsv", str(out),
        str(seq),
    ]

    # Nextclade expects a comma-separated list.
    args.extend(["--output-columns-selection", ",".join(DEFAULT_MUTATION_COLUMNS)])

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
