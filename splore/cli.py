from tempfile import NamedTemporaryFile

import click
import rich
import uvicorn

from splore.utilities import set_env


@click.command()
@click.option(
    "--file",
    "file_path",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    help="The path to the file of molecules (.smi, .sdf, .sdf.gz) to display.",
)
@click.option(
    "--port",
    type=int,
    default=8000,
    show_default=True,
    required=True,
    help="The port to run the GUI on.",
)
def main(file_path, port):

    from splore.db import SploreDB
    from splore.io import molecules_from_file

    console = rich.get_console()
    console.rule("SPLORE")

    molecules = molecules_from_file(file_path)

    with NamedTemporaryFile(suffix=".sqlite") as db_file:

        db = SploreDB(db_file.name)

        with console.status(f"loading [file]{file_path}[/file]"):
            db.create(molecules)

        with set_env(
            SPLORE_DB_PATH=db_file.name,
            SPLORE_API_PORT=f"{port}",
        ):

            rich.print(
                f"The GUI will be available at http://localhost:{port} after a few "
                f"seconds."
            )

            uvicorn.run(
                "splore.app:app",
                host="0.0.0.0",
                port=port,
                log_level="error",
            )
