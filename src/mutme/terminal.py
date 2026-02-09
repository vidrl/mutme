import json
import shlex
from collections.abc import Sequence
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import typer

from .mutation import (
    compare_nextclade_to_annotations,
    get_annotation_fields,
    write_long_format_table,
)
from .presets import (
    transform_3clpro_inhibitor_resistance_table,
    transform_rdrp_inhibitor_resistance_table,
    transform_spike_mab_resistance_table,
)
from .utils import (
    AlignmentPreset,
    CommandExecutionError,
    DatabasePreset,
    cleanup_file,
    log_command_result,
    run_nextclade,
    subset_gff3_by_mutation_prefixes,
    write_rows_to_csv,
)

app = typer.Typer(add_completion=False, pretty_exceptions_enable=False)


def _parse_delimiter(delim: str) -> str:
    """
    Convert common escaped delimiter spellings (e.g. '\\t') into the actual character.
    """
    if delim == r"\t":
        return "\t"
    if len(delim) != 1:
        raise typer.BadParameter(f"Delimiter must be a single character (or '\\t'); got {delim!r}")
    return delim


def version_callback(value: bool):
    if value:
        typer.echo(f"{version('mutme')}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        "-v",
        help="Show version and exit",
        callback=version_callback,
        is_eager=True,
    ),
):
    pass


@dataclass(frozen=True)
class NextcladeInputs:
    tsv: Path | None
    reference: Path | None
    gff: Path | None


def resolve_nextclade_source(
    nextclade_tsv: Path | None,
    reference: Path | None,
    gff: Path | None,
) -> tuple[Path | None, Path | None, Path | None]:
    """
    Returns (precomputed_tsv, reference, gff) after validating combinations.

    If precomputed_tsv is not None => reference and gff must be None.
    If precomputed_tsv is None => reference and gff must both be provided.
    """
    if nextclade_tsv:
        if reference or gff:
            raise typer.BadParameter("Use --nextclade-tsv OR (--reference and --gff) not both.")
        return nextclade_tsv, None, None

    if not reference or not gff:
        missing = [opt for opt, val in (("--reference", reference), ("--gff", gff)) if not val]
        raise typer.BadParameter(
            f"Missing required option(s): {', '.join(missing)}. Provide them or use --nextclade-tsv."
        )
    return None, reference, gff


@app.command()
def run(
    sequences: Path = typer.Option(
        ..., "--sequences", "-s", exists=True, help="Input sequences FASTA"
    ),
    annotations: Path = typer.Option(
        ...,
        "--annotations",
        "-a",
        exists=True,
        help="Annotation table - must contain 'mutation' column that matches format of 'aaSubstitution', 'aaDeletion' or 'aaInsertion' from Nextclade",
    ),
    reference: Path | None = typer.Option(
        None, "--reference", "-r", exists=True, help="Reference FASTA"
    ),
    gff: Path | None = typer.Option(None, "--gff", "-g", exists=True, help="Annotation GFF3"),
    nextclade_tsv: Path | None = typer.Option(
        None,
        "--nextclade-tsv",
        "-n",
        exists=True,
        help="Precomputed Nextclade TSV. If provided, --reference/--gff are not used.",
    ),
    annotations_delimiter: str = typer.Option(
        ",",
        "--annotations-delimiter",
        help="Annotation delimiter - default is comma-delimited for CSV",
        show_default=True,
    ),
    output: Path = typer.Option(
        ..., "--output", "-o", help="Name of output file: {output}.csv", show_default=True
    ),
    output_delimiter: str = typer.Option(
        ",",
        "--output-delimiter",
        help="Annotation delimiter - default is comma-delimited for CSV",
        show_default=True,
    ),
    output_all_columns: bool = typer.Option(
        True,
        "--output-all-columns/--no-output-all-columns",
        help="Output all annotation columns for consistent data format",
        show_default=True,
    ),
    include_comments: bool = typer.Option(
        False,
        "--include-comments",
        "-c",
        help="Include comments in output table if present in annotations table",
        show_default=True,
    ),
    alignment_preset: AlignmentPreset = typer.Option(
        "default", "-p", "--alignment-preset", help="Nextclade alignment preset"
    ),
    nextclade_bin: str = typer.Option(
        "nextclade", "--nextclade-bin", help="Nextclade executable", show_default=True
    ),
    nextclade_extra_args: str = typer.Option(
        "", "--nextclade-extra-args", help="Extra args passed to Nextclade"
    ),
    nextclade_keep_tsv: bool = typer.Option(
        False,
        "--nextclade-keep-tsv",
        help="Keep the Nextclade TSV output in the current working directory on completion",
        show_default=True,
    ),
    nextclade_threads: int | None = typer.Option(
        None,
        "--nextclade-threads",
        "-t",
        min=1,
        help="Number of threads to use for Nextclade (default: all available)",
    ),
) -> None:
    """
    Run Nextclade (reference + GFF), match AA mutations against an annotation table,
    and write a long-format TSV (one row per sequence × mutation hit).
    """

    precomputed_tsv, reference, gff = resolve_nextclade_source(
        nextclade_tsv=nextclade_tsv,
        reference=reference,
        gff=gff,
    )

    produced_nextclade_tsv: Path | None = None

    if precomputed_tsv:
        nextclade_file = precomputed_tsv
    else:
        out_nextclade_tsv = output.with_suffix(".nextclade.tsv")
        produced_nextclade_tsv = out_nextclade_tsv

        typer.echo("Running Nextclade…")
        nextclade_args = shlex.split(nextclade_extra_args) if nextclade_extra_args else []

        assert (
            reference is not None and gff is not None
        )  # guaranteed non-None here but for type checking

        try:
            nextclade_output, _ = run_nextclade(
                sequences_fasta=sequences,
                reference_fasta=reference,
                annotation_gff3=gff,
                out_tsv=out_nextclade_tsv,
                alignment_preset=alignment_preset,
                nextclade_bin=nextclade_bin,
                extra_args=nextclade_args,
                threads=nextclade_threads,
            )
        except CommandExecutionError as e:
            log_command_result(e, "nextclade.error.log")
            raise

        nextclade_file = nextclade_output.tsv

    typer.echo("Comparing mutations to annotation table…")
    reports = compare_nextclade_to_annotations(
        nextclade_tsv=nextclade_file,
        annotation_csv=annotations,
        delimiter=_parse_delimiter(annotations_delimiter),
    )

    annotation_fields = (
        get_annotation_fields(
            annotations_csv=annotations, delimiter=_parse_delimiter(annotations_delimiter)
        )
        if output_all_columns
        else None
    )

    typer.echo("Writing long-format output table")
    write_long_format_table(
        reports,
        output=output,
        delimiter=_parse_delimiter(output_delimiter),
        include_mutation_comments=include_comments,
        all_annotation_fields=annotation_fields,
    )

    # Cleanup only if we actually produced a TSV ourselves
    if produced_nextclade_tsv is not None:
        cleanup = cleanup_file(produced_nextclade_tsv, keep=nextclade_keep_tsv, missing_ok=True)
        if cleanup.removed:
            typer.echo(f"Nextclade output removed: {cleanup.path}")
        else:
            typer.echo(f"Nextclade cleanup: {cleanup.reason} ({cleanup.path})")

    typer.echo(f"Done → {output}")


@app.command()
def preset_scov2(
    input: Path = typer.Option(
        ..., "--input", "-i", exists=True, help="Input database file to transform (CSV)"
    ),
    output: Path = typer.Option(..., "--output", "-o", help="Output annotations file (CSV)"),
    preset: DatabasePreset = typer.Option(
        ..., "--database", "-d", help="Database for which transform presets exist"
    ),
):
    """
    Build annotation presets from the Stanford SARS-CoV-2 resistance databases (2024)
    """

    if preset == "stanford-scov2-mab-resistance":
        write_rows_to_csv(transform_spike_mab_resistance_table(input, add_dms_plus=True), output)
    elif preset == "stanford-scov23clpro-inhibitor":
        write_rows_to_csv(
            transform_3clpro_inhibitor_resistance_table(input, add_pocket_suffix=True),
            output,
        )
    elif preset == "stanford-scov2-rdrp-inhibitor":
        write_rows_to_csv(transform_rdrp_inhibitor_resistance_table(input), output)
    else:
        typer.echo(f"Database preset '{preset}' not supported")


def _parse_attr_keys(keys: str) -> Sequence[str]:
    """
    Support a single comma-separated string (--gene-attr-key gene,Name).
    """
    if "," in keys:
        return tuple(k.strip() for k in keys.split(","))
    return keys


@app.command()
def subset_gff3(
    gff3_path: Path = typer.Argument(
        ..., exists=True, dir_okay=False, readable=True, help="Input GFF3 path (plain text)."
    ),
    table_path: Path = typer.Argument(
        ..., exists=True, dir_okay=False, readable=True, help="Input CSV/TSV path (plain text)."
    ),
    output_path: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        dir_okay=False,
        writable=True,
        help="Output GFF3 path. If omitted, nothing is written (stats still printed).",
    ),
    delimiter: str = typer.Option(
        ",",
        "--delimiter",
        "-d",
        help=r"Table delimiter. Use ',' for CSV or '\t' for TSV.",
        show_default=True,
    ),
    mutation_column: str = typer.Option(
        "mutation",
        "--mutation-column",
        help="Mutation column name (matched case-insensitively).",
        show_default=True,
    ),
    prefix_separator: str = typer.Option(
        ":",
        "--prefix-separator",
        help="Separator used to split mutation values; prefix is before the first separator.",
        show_default=True,
    ),
    strip_prefixes: bool = typer.Option(
        True,
        "--strip-prefixes/--no-strip-prefixes",
        help="Strip whitespace around extracted prefixes.",
        show_default=True,
    ),
    gff_gene_feature_type: str = typer.Option(
        "gene",
        "--gene-feature-type",
        help="GFF3 feature type to treat as genes (3rd column).",
        show_default=True,
    ),
    gene_match_attr_key: str = typer.Option(
        "gene",
        "--gene-attr-key",
        help="Attribute keys to match on gene features. Comma-separated for multiple.",
        show_default=True,
    ),
    include_cds: bool = typer.Option(
        False,
        "--include-cds",
        help="Also include CDS features whose gene attribute matches prefixes.",
        show_default=True,
    ),
    cds_feature_type: str = typer.Option(
        "CDS",
        "--cds-feature-type",
        help="GFF3 feature type to include as CDS (3rd column).",
        show_default=True,
    ),
    cds_gene_attr_key: str = typer.Option(
        "gene",
        "--cds-gene-attr-key",
        help="Attribute key to match on CDS features (9th column).",
        show_default=True,
    ),
    encoding: str = typer.Option(
        "utf-8",
        "--encoding",
        help="Text encoding for reading/writing.",
        show_default=True,
    ),
    json_out: bool = typer.Option(
        False,
        "--json",
        help="Print machine-readable JSON summary to stdout.",
        show_default=True,
    ),
) -> None:
    """
    Subset a GFF3 to only the mutations present in the annotations table
    """
    try:
        delim = _parse_delimiter(delimiter)

        res = subset_gff3_by_mutation_prefixes(
            gff3_path,
            table_path,
            delimiter=delim,
            mutation_column=mutation_column,
            prefix_separator=prefix_separator,
            strip_prefixes=strip_prefixes,
            gff_gene_feature_type=gff_gene_feature_type,
            gene_match_attr_key=gene_match_attr_key,
            include_cds=include_cds,
            cds_feature_type=cds_feature_type,
            cds_gene_attr_key=cds_gene_attr_key,
            output_path=output_path,
            encoding=encoding,
        )

        if json_out:
            payload = {
                "unique_prefixes_count": len(res.unique_prefixes),
                "unique_prefixes": sorted(res.unique_prefixes),
                "matched_gene_ids_count": len(res.matched_gene_ids),
                "matched_gene_ids": sorted(res.matched_gene_ids),
                "genes_written": res.genes_written,
                "cds_written": res.cds_written,
                "rows_read": res.rows_read,
                "warnings_emitted": res.warnings_emitted,
                "output_path": str(output_path) if output_path else None,
            }
            typer.echo(json.dumps(payload, indent=2, sort_keys=True))
        else:
            typer.echo(f"Unique prefixes: {len(res.unique_prefixes)}")
            typer.echo(f"Matched gene IDs: {len(res.matched_gene_ids)}")
            typer.echo(f"Genes written: {res.genes_written}")
            typer.echo(f"CDS written: {res.cds_written}")
            typer.echo(f"Rows read: {res.rows_read}")
            typer.echo(f"Warnings: {res.warnings_emitted}")
            if output_path:
                typer.echo(f"Wrote: {output_path}")

    except typer.BadParameter:
        raise
    except (FileNotFoundError, ValueError) as e:
        raise typer.BadParameter(str(e))
    except Exception as e:
        # Make unexpected errors clearly visible and return non-zero.
        raise typer.Exit(code=1) from e
