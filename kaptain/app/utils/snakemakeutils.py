from pathlib import Path
from typing import Dict, Any, List

import yaml

from kaptain.app.command import Command
from kaptain.app.utils.loggingutils import logger


def generate_config_file(config_data: Dict[str, Any], output_dir: Path, output_basename: str = 'config.yml') -> str:
    """
    Generates a configuration file for Snakemake in YAML file format.
    :param config_data: Configuration data
    :param output_dir: Output directory
    :param output_basename: Output basename
    :return: Path to config file
    """
    config_path = output_dir / output_basename
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    with config_path.open('w') as handle:
        yaml.dump(config_data, handle, sort_keys=False)
    logger.info(f"Configuration file created: {config_path}")
    return str(config_path)


def run_snakemake(snakefile: Path, config_path: Path, targets: List[Path], dir_: Path, threads: int = 8, dry_run: bool = False) -> Command:
    """
    Helper function to run snakemake workflows.
    :param snakefile: Workflow snakefile
    :param config_path: Path to configuration file
    :param targets: Target output files
    :param dir_: Working directory
    :param threads: Number of threads to use
    :param dry_run: Dry run snakemake pipeline
    :return: Snakemake command
    """
    # Create working directory
    if not dir_.exists():
        logger.debug(f'Creating working directory: {dir_}')
        dir_.mkdir(parents=True)

    # Create and run command
    cmd_str = ' '.join([
        'snakemake',
        *[str(x) for x in targets],
        '--snakefile', str(snakefile),
        '--configfile', str(config_path),
        '--cores', str(threads),
    ])

    if dry_run:
        cmd_str += " --dry-run"

    command = Command(cmd_str)
    command.run(dir_)

    # Check if command executed successfully
    if command.exit_code != 0:
        logger.error(command.stderr)
        raise RuntimeError(f'Error executing Snakemake')
    return command
