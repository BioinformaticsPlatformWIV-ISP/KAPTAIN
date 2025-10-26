#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Tuple, List

from kaptain import initialize_logging, logger, KAPTAIN


def parse_galaxy_args() -> Tuple[argparse.Namespace, List[str]]:
    """
    Parses the Galaxy arguments.
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--ont-in', nargs="+", required=True, type=Path)
    parser.add_argument('--ont-in-name', nargs="+", type=str)
    parser.add_argument('--dir-working', default=Path.cwd(), type=Path, help='Working directory')
    return parser.parse_known_args()


def _sanitize_input_name(name: str) -> str:
    """
    Sanitizes the input file name.
    :param name: Name
    :return: None
    """

    valid_extensions = [
        "fasta", "fa", "fastq", "fq",
        "fasta.gz", "fa.gz", "fastq.gz", "fq.gz"
    ]

    invalid_chars = '/!@#$\\'

    # Replace spaces by dashes
    name = name.replace(' ', '_')

    # Avoid double dot before the extension
    if name.endswith('.'):
        name = name[:-1]

    # Detect if the name already has one of the valid extensions
    lower_name = name.lower()
    matched_ext = next((ext for ext in valid_extensions if lower_name.endswith(f".{ext}")), None)

    # Sanitize invalid characters
    name = ''.join(c for c in name if c not in invalid_chars)

    # If it already has a valid extension, return as-is (cleaned)
    if matched_ext:
        return name

    # Otherwise, default to `.fasta`
    return f"{name}.fasta"


def main() -> None:
    """
    Runs the main script.
    :return: None
    """
    initialize_logging()
    logger.info(f'Running KAPTAIN through Galaxy')
    args, unparsed_args = parse_galaxy_args()

    # Create input directory
    dir_in = args.dir_working / 'input'
    dir_in.mkdir(exist_ok=True)

    # Symlink
    unparsed_args.append('--ont-in')
    for ont_in, ont_in_name in zip(args.ont_in, args.ont_in_name):
        sym_ont_in = dir_in / _sanitize_input_name(ont_in_name)
        if not sym_ont_in.exists():
            sym_ont_in.symlink_to(Path(ont_in))
        logger.info(f'({ont_in}, {ont_in_name}) --> {sym_ont_in}')
        unparsed_args.append(str(sym_ont_in.absolute()))

    # Use shm on Galaxy
    unparsed_args.append('--shm')

    # Run the main script
    workflow = KAPTAIN(unparsed_args)
    workflow.run()


if __name__ == '__main__':
    main()
