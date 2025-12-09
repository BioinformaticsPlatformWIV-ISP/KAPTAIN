import argparse
import shutil
from importlib.resources import files
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

from kaptain.app import __version__
from kaptain.app.command import Command
from kaptain.app.utils import snakemakeutils
from kaptain.app.utils.loggingutils import initialize_logging, logger


class KAPTAIN:
    """
    Main class to run the KAPTAIN pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        :return: None
        """
        self._parser = KAPTAIN._parser()
        self._args = self._parser.parse_args(args=args)
        initialize_logging(self._args.log)
        self._args.dir_working.mkdir(parents=True, exist_ok=True)
        self._path_html_out = self._args.output.absolute() / self._args.output_html
        if self._path_html_out.exists():
            self._path_html_out.unlink()

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        logger.info('Starting KAPTAIN workflow')
        self._check_dependencies()
        self._check_args()
        self._validate_lookup_file()
        self._check_database()
        path_snakefile = Path(str(files('kaptain').joinpath('resources/main_workflow.smk')))
        # Run Snakemake workflow
        path_config = self.__create_config_file()
        snakemakeutils.run_snakemake(path_snakefile,
                                     path_config,
                                     [self._path_html_out],
                                     self._args.dir_working,
                                     self._args.threads,
                                     self._args.dry_run
                                     )
        logger.info(f'Workflow finished successfully, output available in: {self._args.output}')

        # Copy the log file if it exists
        if (self._args.dir_working / 'kt.log').exists():
            shutil.copyfile(self._args.dir_working / 'kt.log', self._args.output / 'kt.log')

    @staticmethod
    def _parser() -> argparse.ArgumentParser:
        """
        Parses the command line arguments.
        :return: Parsed arguments.
        """

        class SmartFormatter(argparse.ArgumentDefaultsHelpFormatter,
                             argparse.RawDescriptionHelpFormatter):
            """Combine RawDescriptionHelpFormatter and ArgumentDefaultsHelpFormatter."""
            pass

        description = Path(str(files('kaptain') / 'resources' / 'usage.txt')).read_text()

        parser = argparse.ArgumentParser(prog="kaptain",
                                         add_help=True,
                                         description=description,
                                         formatter_class=SmartFormatter)

        # Input and output
        parser.add_argument('-i', '--ont-in', nargs="+", required=True,
                            help='ONT input FASTA/Q file')
        parser.add_argument('--db', required=True, type=Path,
                            help='Basename of KMA database. Database exists of four files named *.{comp.b, length.b, name, seq.b}')
        parser.add_argument('--db-lookup', required=True,
                            help='Lookup file of KMA database')
        parser.add_argument('--dir-working', type=Path, default=Path.cwd() / 'working',
                            help='Working directory')
        parser.add_argument('-o', '--output', required=True, type=Path,
                            help='Output directory')
        parser.add_argument('--output-html', type=Path, default='report.html',
                            help='Output report name')

        # Parameters
        parser.add_argument("--subsampling", nargs="+", type=lambda x: None if x == "None" else x, default=[None],
                            choices=["200M", "500M", "1000M", "1500M", "2000M", None],
                            help='Subsample input to number of bases before classification. Leave empty or use None for no downsampling.')
        parser.add_argument("--fdr", type=int, nargs="+", default=[5],
                            choices=[15, 10, 5, 1],
                            help='FDR setting.')

        # Advanced and hidden
        parser.add_argument('--thresholds', default=Path(str(files('kaptain').joinpath('resources/thresholds.csv'))),
                            help=argparse.SUPPRESS)
        parser.add_argument("--shm", action="store_true",
                            help=argparse.SUPPRESS)
        parser.add_argument("--dry-run", action="store_true",
                            help=argparse.SUPPRESS)
        parser.add_argument("--seed", type=int, default=0,
                            help=argparse.SUPPRESS)  # set seed for Rasusa

        # Supplements
        parser.add_argument('--threads', type=int, default=4,
                            help="Number of threads")
        parser.add_argument('--version', action='version', version=f'KAPTAIN {__version__}',
                            help='Print version and exit')
        parser.add_argument("--log", action="store_true",
                            help="Write out log information to file")

        return parser

    def __create_config_file(self) -> Path:
        """
        Creates the config file for the Snakemake workflow.
        :return: Path to the config file
        """
        queries = {}
        for ont_in, fdr, subsampling in zip(self._args.ont_in, self._args.fdr, self._args.subsampling):
            name = Path(ont_in).name[: -sum(len(s) for s in Path(ont_in).suffixes)]
            entry = queries.setdefault(name, {"file": ont_in, "combinations": []})
            entry["combinations"].append({
                "fdr": fdr,
                "subsampling": subsampling
            })

        self._config_data = {
            'output': {'dir': str(self._args.output.absolute()), 'html': str(self._path_html_out)},
            'wd': str(self._args.dir_working),
            'shm': self._args.shm,
            'db': str(self._args.db),
            'lookup': self._args.db_lookup,
            'subsampling_seed': self._args.seed,
            'thresholds': str(self._args.thresholds),
            "queries": queries
        }
        return Path(snakemakeutils.generate_config_file(self._config_data, self._args.dir_working.absolute()))

    def _check_args(self) -> None:
        """
        Checks conditional args.
        :return: None
        """

        def _normalize_arg_length(arg_list, reference_list, name):
            """Ensure arg_list has either one value or the same length as reference_list."""
            if len(arg_list) == 1:
                return arg_list * len(reference_list)
            if len(arg_list) != len(reference_list):
                raise self._parser.error(f"Either provide one {name} value or as many as input files.")
            return arg_list

        self._args.fdr = _normalize_arg_length(self._args.fdr, self._args.ont_in, "FDR")
        self._args.subsampling = _normalize_arg_length(self._args.subsampling, self._args.ont_in, "subsampling")

    def _check_dependencies(self) -> None:
        """
        Checks if the required dependencies are available.
        :return: None
        """
        commands = {
            'kma': 'kma -v',
            'rasusa': 'rasusa --version',
            'seqtk': 'seqtk'
        }
        if not self._args.subsampling:
            commands.pop('rasusa')
        logger.info('Checking dependencies')
        for tool, command in commands.items():
            command = Command(command)
            command.run(self._args.dir_working, disable_logging=True)
            # Use code 127. seqtk never returns 0 unless a full command is run successfully.
            if command.exit_code == 127:
                raise RuntimeError(f"Dependency '{tool}' not available")
            logger.info(f"{tool}: OK")

    def _validate_lookup_file(self) -> None:
        """
        Checks if the provided lookup file is valid.
        :return: None
        """
        df_lookup = pd.read_table(self._args.db_lookup,
                                  names=['full_name', 'length', 'accession', 'seq_type', 'species_name'],
                                  )

        if not df_lookup["length"].dtype == "int64":
            raise ValueError('The second column of --db_lookup file contains a non-integer')
        invalid_seq_type = df_lookup.loc[~df_lookup["seq_type"].str.lower().isin(["plasmid", "chromosome"]), "seq_type"].unique()
        if len(invalid_seq_type):
            raise ValueError(f"Invalid values in fourth column of --db-lookup: {invalid_seq_type}")

    def _check_database(self) -> None:
        """
       Checks if the database is available.
       :return: None
       """
        required_suffixes = ["comp.b", "length.b", "name", "seq.b"]

        missing_files = []
        for suffix in required_suffixes:
            file_path = self._args.db.with_suffix(f".{suffix}")
            if not file_path.exists():
                missing_files.append(file_path)

        if missing_files:
            missing_str = "\n  ".join(str(f) for f in missing_files)
            raise FileNotFoundError(
                f"The following KMA database files are missing:\n  {missing_str}"
            )
