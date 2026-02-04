import typer

from .utils import run_command

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

app()