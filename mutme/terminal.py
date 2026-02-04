import shlex
import typer

from pathlib import Path

from .utils import run_command,run_nextclade,AlignmentPreset
from .mutation import compare_nextclade_to_annotations, write_long_format_hits_tsv

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
    annotations: Path = typer.Option(..., "--annotations", "-a", exists=True, help="Annotation CSV - must contain 'mutation' column that matches format of 'aaSubstitution', 'aaDeletion' or 'aaInsertion' from Nextclade"),
    delimiter: str = typer.Option(",", "--annotations-delimiter", "-d", exists=True, help="Annotation delimiter - default is comma-delimited for CSV"),
    output: Path = typer.Option(..., "--output", "-o", help="Name of output file: {output}.tsv"),
    alignment_preset: AlignmentPreset = typer.Option("default", "-p", "--alignment-preset", help="Nextclade alignment preset"),
    nextclade_bin: str = typer.Option("nextclade", "--nextclade-bin", help="Nextclade executable"),
    nextclade_extra_args: str = typer.Option("", "--nextclade-extra-args", help="Extra args passed to Nextclade"),
) -> None:
    """
    Run Nextclade (reference + GFF), match AA mutations against an annotation table,
    and write a long-format TSV (one row per sequence × mutation hit).
    """

    out_nextclade_tsv = output.with_suffix(".nextclade.tsv")

    typer.echo("Running Nextclade…")
    nextclade_args = shlex.split(nextclade_extra_args) if nextclade_extra_args else []

    nextclade_output, _ = run_nextclade(
        sequences_fasta=sequences,
        reference_fasta=reference,
        annotation_gff3=gff,
        out_tsv=out_nextclade_tsv,
        alignment_preset=alignment_preset,  # type: ignore[arg-type]
        nextclade_bin=nextclade_bin,
        extra_args=nextclade_args,
    )

    typer.echo("Comparing mutations to annotation table…")
    reports = compare_nextclade_to_annotations(
        nextclade_tsv=nextclade_output.tsv,
        annotation_csv=annotations,
        delimiter=delimiter
    )

    typer.echo("Writing long-format output table")
    write_long_format_hits_tsv(
        reports,
        out_tsv=output,
    )

    typer.echo(f"Done → {output}")

app()