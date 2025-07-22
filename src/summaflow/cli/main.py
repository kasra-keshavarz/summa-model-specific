import click
from pathlib import Path
from ..core import SUMMAWorkflow


@click.command()
@click.option('--json', 'json_file', required=True, type=click.Path(exists=True, path_type=Path),
              help='Path to the JSON configuration file')
@click.option('--output-path', default=None, type=click.Path(path_type=Path),
              help='Output path for saving results (optional)')
def main(json_file, output_path):
    """Run SUMMAWorkflow from a JSON configuration file using MAF outputs."""
    try:
        # Load workflow from JSON file
        click.echo(f"Loading workflow from {json_file}")
        workflow = SUMMAWorkflow._from_maf_json_file(json_file)

        # Execute the workflow
        click.echo("Running workflow...")
        workflow.run(save=True, path=output_path)

        click.echo(f"Saving results to {output_path}")
        click.echo("Workflow completed successfully!")


    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


if __name__ == '__main__':
    main()