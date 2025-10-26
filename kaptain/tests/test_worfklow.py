from importlib.resources import files
from pathlib import Path
from typing import List

import pytest

from kaptain import KAPTAIN
from kaptain.app.command import Command


def build_args(queries: List[str], fdrs: List[int], subsamplings: List[str], tmp_path) -> List[str]:
    """
    Builds command-line arguments for the kaptain workflow.
    :param queries: Query file paths
    :param fdrs: FDR thresholds
    :param subsamplings: Subsampling levels; None values become "full"
    :param tmp_path: Temporary directory path for output
    :return: Command-line arguments
    """

    # Normalize inputs so they’re always lists (or None)
    if fdrs and not isinstance(fdrs, (list, tuple)):
        fdrs = [fdrs]

    if subsamplings and not isinstance(subsamplings, (list, tuple)):
        subsamplings = [subsamplings]

    path_testdata = Path(str(files('kaptain'))).joinpath('resources/testdata')

    fdr_name = '_'.join(map(str, fdrs)) if fdrs else 'default'
    sub_name = '_'.join(str(x or 'full') for x in subsamplings) if subsamplings else 'full'
    dir_name = f"output_fdr_{fdr_name}_sub_{sub_name}"
    dir_out = tmp_path / dir_name

    query_files = [str(path_testdata / query_file) for query_file in queries]

    args = [
        '--output', str(dir_out),
        '--dir-working', str(dir_out / 'working'),
        '--ont-in', *query_files,
        '--db', str(tmp_path / 'db' / 'kma_db'),
        '--db-lookup', str(path_testdata / 'database' / 'db_lookup.tsv'),
        '--threads', '4',
        '--log'
    ]

    if subsamplings:
        args += ['--subsampling', *map(str, subsamplings)]

    if fdrs:
        args += ['--fdr', *map(str, fdrs)]

    return args


def build_test_db(fasta_db: Path, output_path: Path) -> None:
    """
    Build the test database before running tests.
    :param fasta_db: Path to the fasta/q used for the test DB
    :param output_path: Output path for the test DB
    :return: None
    """
    # First check if KMA is available
    command = Command("kma -v")
    command.run(output_path, disable_logging=True)
    if command.exit_code == 127:
        raise RuntimeError(f"Dependency 'KMA' not available")

    # Build DB
    output_path_db = output_path / 'db'
    output_path_db.mkdir(parents=False, exist_ok=True)
    command = Command(f"kma index -i {fasta_db} -o {output_path_db / 'kma_db'}")
    command.run(output_path, disable_logging=False)
    if command.exit_code != 0:
        raise RuntimeError(f"Building test database failed")


VALID_COMBINATIONS = [
    (["query.fasta", "query2.fasta", "query2.fasta"], [1, 5, 10], [None, "200M", "2000M"]),
    (["query.fasta", "query2.fasta", "query2.fasta"], [15], None),
    (["query.fasta", "query2.fasta", "query2.fasta"], [15], "200M"),
    (["query.fasta"], [15], "200M"),
    (["query.fasta"], None, "200M"),
    (["query.fasta"], None, None),
]


@pytest.mark.parametrize("query, fdr, subsampling", VALID_COMBINATIONS)
def test_workflow_valid(tmp_path, query, fdr, subsampling):
    """Runs kaptain workflow with valid FDR/subsampling combinations."""
    args = build_args(query, fdr, subsampling, tmp_path)
    build_test_db(Path(str(files('kaptain'))).joinpath('resources/testdata/database/database.fasta'), tmp_path)
    workflow = KAPTAIN(args)
    workflow.run()


INVALID_COMBINATIONS = [
    (["query.fasta"], [1, 5], [None, "200M", "2000M"]),
    (["query.fasta"], [1, 5, 10], [None, "200M", "2000M"]),
    (["query.fasta"], [1, 5, 10], [None, "200M"]),
]


@pytest.mark.parametrize("query, fdr, subsampling", INVALID_COMBINATIONS)
def test_workflow_invalid(tmp_path, query, fdr, subsampling):
    """Runs kaptain workflow with invalid inputs; should raise RuntimeError."""
    args = build_args(query, fdr, subsampling, tmp_path)
    with pytest.raises(SystemExit):
        workflow = KAPTAIN(args)
        workflow.run()
