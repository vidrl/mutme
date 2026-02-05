import shlex
import typer

from pathlib import Path

from .utils import run_command, run_nextclade, write_rows_to_csv, log_command_result, cleanup_file, AlignmentPreset, DatabasePreset, CommandExecutionError
from .mutation import compare_nextclade_to_annotations, write_long_format_table
from .presets import transform_spike_mab_resistance_table, transform_3clpro_inhibitor_resistance_table, transform_rdrp_inhibitor_resistance_table



app = typer.Typer(add_completion=False, pretty_exceptions_enable=False)


@app.command()
def test():    
    """Test"""
    typer.echo("ok")


@app.command()
def nextclade_version():    
    """Print Nextclade version"""
    result = run_command(["nextclade", "--version"])
    print(result.stdout)



@app.command()
def run(
    sequences: Path = typer.Option(..., "--sequences", "-s", exists=True, help="Input sequences FASTA"),
    reference: Path = typer.Option(..., "--reference", "-r", exists=True, help="Reference FASTA"),
    gff: Path = typer.Option(..., "--gff", "-g", exists=True, help="Annotation GFF3"),
    annotations: Path = typer.Option(..., "--annotations", "-a", exists=True, help="Annotation table - must contain 'mutation' column that matches format of 'aaSubstitution', 'aaDeletion' or 'aaInsertion' from Nextclade"),
    annotations_delimiter: str = typer.Option(",", "--annotations-delimiter", help="Annotation delimiter - default is comma-delimited for CSV"),
    output: Path = typer.Option(..., "--output", "-o", help="Name of output file: {output}.csv"),
    output_delimiter: str = typer.Option(",", "--output-delimiter", help="Annotation delimiter - default is comma-delimited for CSV"),
    include_comments: bool = typer.Option(False, "--include-comments", "-c", help="Include comments in output table if present in annotations table"),
    alignment_preset: AlignmentPreset = typer.Option("default", "-p", "--alignment-preset", help="Nextclade alignment preset"),
    nextclade_bin: str = typer.Option("nextclade", "--nextclade-bin", help="Nextclade executable"),
    nextclade_extra_args: str = typer.Option("", "--nextclade-extra-args", help="Extra args passed to Nextclade"),
    nextclade_keep_tsv: bool = typer.Option(False, "--nextclade-keep-tsv", help="Keep the Nextclade TSV output in the current working directory on completion"),
) -> None:
    """
    Run Nextclade (reference + GFF), match AA mutations against an annotation table,
    and write a long-format TSV (one row per sequence × mutation hit).
    """

    out_nextclade_tsv = output.with_suffix(".nextclade.tsv")

    typer.echo("Running Nextclade…")
    nextclade_args = shlex.split(nextclade_extra_args) if nextclade_extra_args else []

    try:
        nextclade_output, _ = run_nextclade(
            sequences_fasta=sequences,
            reference_fasta=reference,
            annotation_gff3=gff,
            out_tsv=out_nextclade_tsv,
            alignment_preset=alignment_preset,
            nextclade_bin=nextclade_bin,
            extra_args=nextclade_args,
        )
    except CommandExecutionError as e:
        log_command_result(e, "nextclade.error.log")
        raise
    

    typer.echo("Comparing mutations to annotation table…")
    reports = compare_nextclade_to_annotations(
        nextclade_tsv=nextclade_output.tsv,
        annotation_csv=annotations,
        delimiter=annotations_delimiter
    )

    typer.echo("Writing long-format output table")
    write_long_format_table(
        reports,
        output=output,
        delimiter=output_delimiter,
        include_mutation_comments=include_comments
    )

    cleanup = cleanup_file(nextclade_output.tsv, keep=nextclade_keep_tsv, missing_ok=True)
    if cleanup.removed:
        typer.echo(f"Nextclade output removed: {cleanup.path}")
    else:
        # Failure of cleanup is not fatal
        typer.echo(f"Nextclade cleanup: {cleanup.reason} ({cleanup.path})")

    typer.echo(f"Done → {output}")



@app.command()
def prepare_database(
    input: Path = typer.Option(..., "--input", "-i", exists=True, help="Input database file to transform (CSV)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output annotations file (CSV)"),
    preset: DatabasePreset = typer.Option(..., "--database", "-d", help="Database for which transform presets exist"),
):
    if preset == "sars-cov-2-mab-resistance":
        write_rows_to_csv(transform_spike_mab_resistance_table(input, add_dms_plus=True), output)
    elif preset == "sars-cov-2-3clpro-inhibitor":
        write_rows_to_csv(transform_3clpro_inhibitor_resistance_table(input, add_pocket_suffix=True), output)
    elif preset == "sars-cov-2-rdrp-inhibitor":
        write_rows_to_csv(transform_rdrp_inhibitor_resistance_table(input), output)
    else:
        pass


app()