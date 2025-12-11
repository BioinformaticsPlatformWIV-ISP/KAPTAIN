import datetime
from pathlib import Path
from typing import Any, Dict
import json
import pandas as pd

from kaptain.app import __version__
from kaptain.app.report.htmlcitation import HtmlCitation
from kaptain.app.report.htmlreportsection import HtmlReportSection
from kaptain.app.report.htmltablecell import HtmlTableCell


def create_analysis_info_section(config_: Dict[str, Any]) -> HtmlReportSection:
    """
    Creates the analysis info section.
    :param config_: Configuration data
    :return: Analysis info section
    """
    section = HtmlReportSection('Analysis info')

    # Analysis info
    section.add_table([
        ['Workflow version:', __version__],
        ['Analysis date:', datetime.datetime.now().strftime('%d/%m/%Y - %X')],
        ['Nb. of samples:', str(len(config_['queries'].keys()))],
    ], table_attributes=[('class', 'information')])

    # TODO once paper is available
    # Citation disclaimer
    # section.add_header('Disclaimer', 2)
    # section.add_paragraph('If you use this pipeline for your scientific work, please cite:')
    # section.add_html_object(HtmlCitation.parse_from_json('Bogaerts_2023-ont_outbreak'))

    return section


def create_thresholds_table(thresholds: pd.DataFrame) -> HtmlReportSection:
    """
    Creates the analysis info section.
    :param thresholds: Thresholds table
    :return: Threshold info section
    """
    section = HtmlReportSection('Reference thresholds for template ID')

    thresholds = thresholds.map('{:,.2f}'.format)

    # Analysis info
    section.add_table(thresholds.to_records().tolist(),
                      column_names=[""] + thresholds.columns.tolist(),
                      table_attributes=[('class', 'matrix')])
    section.add_paragraph(
        f"Based on ten defined mock communities, the expected precision is 100% - FDR X% if the above thresholds are applied.\n"
        f"E.g., the precision is 95% for FDR 5% (100% - 5%)."
        f"This means 95% of detected species above this threshold can be considered correct."
    )

    return section


def warning_message(query_info: dict) -> str:
    """
    Generate warning messages based on subsampling and yield thresholds.
    """
    icon = {"note": "&#8505;&#65039;",
            "warning": "&#10071;"}

    n_bases_m = query_info["n_bases"] / 1_000_000
    subsampling = query_info["subsampling"]
    subsample_value = int(subsampling.rstrip("M")) if subsampling else 0

    LOWER_YIELD = 200
    UPPER_YIELD = 2000
    ERROR_BOUND = 0.01  # 1%

    LOWER_YIELD_LIMIT = LOWER_YIELD * (1 - ERROR_BOUND)
    UPPER_YIELD_LIMIT = UPPER_YIELD * (1 + ERROR_BOUND)
    SUBSAMPLE_VALUE_LIMIT = subsample_value * (1 - ERROR_BOUND)

    rules = [
        (
            "note",
            n_bases_m < SUBSAMPLE_VALUE_LIMIT,
            f"The sample's yield is below the requested subsample setting ({subsample_value}M). "
        ),
        (
            "warning",
            n_bases_m < LOWER_YIELD_LIMIT,
            f"The sample's yield is below the lowest possible setting ({LOWER_YIELD}M). "
            "The set thresholds may not be appropriate."
        ),
        (
            "warning",
            n_bases_m > UPPER_YIELD_LIMIT,
            f"The sample's yield is above the highest possible setting ({UPPER_YIELD}M). "
            "The set thresholds may not be appropriate. Try downsampling if possible."
        ),
    ]

    warnings = [
        f"{icon[message_type]} {message}"
        for message_type, condition, message in rules if condition
    ]

    return "\n".join(warnings)


def create_parameter_section(query_information_files: list[Path]) -> HtmlReportSection:
    """
    Creates the parameter section.
    :param query_information_files: List of Paths pointing to information on the queries
    :return: Section
    """
    section = HtmlReportSection('Parameters')

    data_table = []
    for query_information_file in query_information_files:
        with open(query_information_file, "r") as f:
            query_information = json.load(f)
        data_table.append([
            query_information["query"],
            query_information["subsampling"] or "NA",
            query_information["n_bases"],
            query_information["closest_yield"],
            query_information["fdr"],
            query_information["closest_threshold"],
            warning_message(query_information)
        ])
    section.add_table(data_table,
                      column_names=["Input Sample", "Requested Subsampling", "Final Yield (bases)", "Closest Yield Setting",
                                    "Requested FDR", "Selected Template ID Threshold", "Note"],
                      table_attributes=[('class', 'data')]
                      )

    return section


def __get_colored_cell_template_id(value: float, threshold: float) -> HtmlTableCell:
    """
    Returns a colored cell.
    :param value: Template ID value
    :return: Table cell
    """
    value_str = f'{value:.2f}%'
    if value < threshold:
        color = 'red'
    else:
        color = 'green'

    return HtmlTableCell(value_str, color=color)


def create_mapping_section(query_information_files: list[Path]) -> HtmlReportSection:
    """
    Creates the read mapping section.
    :param query_information_files: List of Paths pointing to information on the queries
    :return: Section
    """
    section = HtmlReportSection("KMA MAPPING")

    for query_information_file in query_information_files:
        with open(query_information_file, "r") as f:
            query_information = json.load(f)
        mapping_result_file = query_information.get("query_file_path")
        sample_name = Path(mapping_result_file).parts[-4]
        kma_results = (pd.read_table(mapping_result_file,
                                     usecols=["genome", "species_name", "Read_Count_Aln",
                                              "Template_Identity", 'Template_Coverage', 'Template_Depth'],
                                     )
                       .sort_values(by=["Template_Identity"], ascending=False)
                       .round(3)
                       .rename(columns={"genome": "Accession",
                                        "species_name": "Species",
                                        "Read_Count_Aln": "Aligned reads",
                                        "Template_Identity": "Template ID (%)",
                                        'Template_Coverage': "Template Coverage (%)",
                                        'Template_Depth': "Template Depth"
                                        }
                               )
                       )
        header = kma_results.columns
        section.add_header(f"{sample_name.upper()} - {query_information['subsampling'] or 'FULL'} - FDR {query_information['fdr']}%",
                           level=3)

        section.add_table([[row["Accession"],
                            row["Species"],
                            row["Aligned reads"],
                            __get_colored_cell_template_id(row["Template ID (%)"], query_information['closest_threshold']),
                            row["Template Coverage (%)"],
                            row["Template Depth"],
                            ] for row in kma_results.to_dict('records')], header, [('class', 'data')])
        section.add_paragraph(
            f"Species with a template ID below or above the threshold {query_information['closest_threshold']:.2f}% are displayed in red "
            f"or green, respectively."
        )

    return section


def create_citations_section() -> HtmlReportSection:
    """
    Creates the report section with the citations.
    :return: Section
    """
    section = HtmlReportSection('Citations')
    keys = [
        'Clausen_2018-kma',
        'Hall_2022-rasusa',
    ]
    for citation_key in keys:
        section.add_html_object(HtmlCitation.parse_from_json(citation_key))

    return section
