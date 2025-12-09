from pathlib import Path

import numpy as np
import pandas as pd

from kaptain.app.report.htmlreport import HtmlReport

WORKING_DIR = Path(config['wd']) / '{sample_name}' / '{sample_yield}'

rule all:
    """
    Entry point for the SnakeMake workflow
    """
    input:
        HTML= Path(config['output']['html'])

rule subsample_input_fastq:
    """
    Subsamples the input fastq file to a lower yield (bases).
    """
    input:
        fq = lambda wildcards: config['queries'][wildcards.sample_name]["file"]
    output:
        subsampled_fq =WORKING_DIR / 'filtered_fastq' / 'subsampled.fq'
    log:
        WORKING_DIR / 'filtered_fastq' / 'subsampled.log'
    params:
        sample_yield = lambda wildcards: wildcards.sample_yield,
        seed = config['subsampling_seed']
    shell:
        """
        rasusa reads \
            -v \
            -s {params.seed} \
            -b {params.sample_yield} \
            -o {output} \
            {input} 2> {log}
        """

def get_fastq_input(wildcards):
    """
    Determines whether the downsampled of full input is used.
    """
    if wildcards.sample_yield != "full":
        return rules.subsample_input_fastq.output.subsampled_fq
    else:
        return config['queries'][wildcards.sample_name]["file"]

rule get_nb_of_bases:
    """
    Collect the number of bases on the sample to be classified (downsampled or not)
    """
    input:
        get_fastq_input
    output:
        N_BASES =WORKING_DIR / 'n_bases.txt'
    shell:
        """
        seqtk size {input} | cut -f 2 > {output.N_BASES}
        """

rule kma_mapping:
    """
    Maps the input to the reference database with KMA.
    """
    input:
        get_fastq_input
    output:
        mapstat =WORKING_DIR / 'classification' / 'kma.mapstat',
        tsv =WORKING_DIR / 'classification' / 'kma.tsv'
    params:
        DB = lambda wildcards: config["db"],
        prefix = lambda wildcards, output: Path(output.mapstat).parent / 'kma',
        MP = 20,
        MRS = 0.0,
        BC = 0.7,
        SHM = lambda wildcards: "-shm" if config['shm'] else ""
    threads: 4
    shell:
        """
        kma \
        -i {input} \
        -o {params.prefix} \
        -t_db {params.DB} \
        -t {threads} \
        -mrs {params.MRS} \
        -bcNano \
        -bc {params.BC} \
        -ef \
        -a \
        -mem_mode \
        -1t1 \
        -matrix \
        -ID 0.01 \
        -tsv \
        {params.SHM} 
        """

rule merge_kma_tsv_and_mapstat:
    """
    Merges the *.tsv file with the *.mapstat file. 
    """
    input:
        tsv = rules.kma_mapping.output.tsv,
        mapstat = rules.kma_mapping.output.mapstat
    output:
        tsv_mapstat =WORKING_DIR / 'classification' / 'kma_tsv_mapstat.tsv'
    run:
        tsv = pd.read_table(input.tsv, usecols=['Template_Name', 'Template_Length', 'Template_Identity',
                                                'Template_Coverage', 'Template_Depth', 'Query_Identity',
                                                'Query_Coverage', 'Query_Depth', 'Read_Count_Map', 'Read_Count_Aln',
                                                'Score']
        )

        mapstat = pd.read_table(input.mapstat, skiprows=6)
        mapstat.rename(columns={'# refSequence': 'Template_Name'}, inplace=True)

        tsv_mapstat = tsv.merge(
            mapstat.loc[:, ["Template_Name", "refConsensusSum", "refCoveredPositions", "bpTotal"]],
            on="Template_Name"
        )

        tsv_mapstat.to_csv(output.tsv_mapstat, sep='\t', header=True, index=False)

rule merge_kma_tsv_mapstat_and_lookup:
    """
    Add additional info on found templates.
    Also add sequences of genomes if another part (sequence) of the genome was found. 
    """
    input:
        tsv_mapstat=rules.merge_kma_tsv_and_mapstat.output.tsv_mapstat
    output:
        tsv_mapstat_lookup =WORKING_DIR / 'classification' / 'kma_tsv_mapstat_extended.tsv'
    params:
        lookup = config['lookup']
    run:
        tsv_mapstat = pd.read_table(input.tsv_mapstat)

        lookup = pd.read_table(params.lookup, names=['full_name', 'length', 'accession', 'seq_type', 'species_name'])

        tsv_mapstat['Template_Name_Key'] = tsv_mapstat['Template_Name'].str.split().str[0]
        lookup['full_name_key'] = lookup['full_name'].str.split().str[0]


        missing = set(tsv_mapstat['Template_Name_Key']) - set(lookup['full_name_key'])
        if missing:
            raise AssertionError(f"Missing Template_Names in lookup: {missing}")

        present_accessions = tsv_mapstat['Template_Name_Key'].map(lookup.set_index('full_name_key')['accession'])
        lookup_present_accessions = lookup.loc[lookup['accession'].isin(present_accessions)]
        tsv_mapstat_lookup = tsv_mapstat.merge(lookup_present_accessions,
                                               left_on='Template_Name_Key',
                                               right_on='full_name_key', how="right"
        )

        tsv_mapstat_lookup.to_csv(output.tsv_mapstat_lookup, sep='\t', header=True, index=False)

rule remove_plasmids_and_recalculate:
    """
    Remove plasmid hits.
    Recalculate metrics based on ALL sequences from a genome, not only the found sequences.
    """
    input:
        tsv_mapstat_lookup=rules.merge_kma_tsv_mapstat_and_lookup.output.tsv_mapstat_lookup
    output:
        tsv_mapstat_wo_pl_collapsed =WORKING_DIR / 'classification' / 'kma_tsv_mapstat_extended_wo_pl_collapsed.tsv'
    run:
        tsv_mapstat_lookup = pd.read_table(input.tsv_mapstat_lookup)

        tsv_mapstat_lookup_filtered = tsv_mapstat_lookup[~tsv_mapstat_lookup['seq_type'].str.contains('plasmid', case=False, na=False)]

        tsv_mapstat_lookup_filtered_agg = (tsv_mapstat_lookup_filtered
                                           .groupby(['accession', 'species_name'])
                                           .agg(
                                                Genome_length_wo_pl=("length", "sum"),
                                                Read_Count_Aln=("Read_Count_Aln", "sum"),
                                                bpTotal=("bpTotal", "sum"),
                                                refConsensusSum=("refConsensusSum", "sum"),
                                                refCoveredPositions=("refCoveredPositions", "sum"),
                                                conclaveScore=("Score", "sum")
                                                )
                                            .assign(
                                                Template_Identity=lambda d: 100 * d["refConsensusSum"] / d["Genome_length_wo_pl"],
                                                Template_Coverage=lambda d: 100 * d["refCoveredPositions"] / d["Genome_length_wo_pl"],
                                                Template_Depth=lambda d: d["bpTotal"] / d["Genome_length_wo_pl"],
                                                Query_Identity=lambda d: np.where(
                                                    d["refCoveredPositions"] > 0,
                                                    100 * d["refConsensusSum"] / d["refCoveredPositions"],
                                                    0
                                                ),
                                                Query_Depth=lambda d: np.where(
                                                    d["refCoveredPositions"] > 0,
                                                    d["bpTotal"] / d["refCoveredPositions"],
                                                    0
                                                )
                                            )
                                           .astype({
                                                    "Read_Count_Aln": "int",
                                                    "bpTotal": "int",
                                                    "refConsensusSum": "int",
                                                    "refCoveredPositions": "int",
                                                    "conclaveScore": "int"
                                                })
                                            .reset_index()
                                           )
        # Remove if Read_Count_Aln is 0 (i.e., only plasmid hits)
        tsv_mapstat_lookup_filtered_agg = tsv_mapstat_lookup_filtered_agg[tsv_mapstat_lookup_filtered_agg["Read_Count_Aln"] != 0]

        tsv_mapstat_lookup_filtered_agg = (tsv_mapstat_lookup_filtered_agg
                                           .merge(tsv_mapstat_lookup
                                                  .groupby(['accession', 'species_name'])
                                                  .agg(Genome_length=("length", "sum"))
                                                  .reset_index(),
                                                  on=['accession', 'species_name']))

        tsv_mapstat_lookup_filtered_agg.rename(columns={"accession": "genome"}, inplace=True)

        tsv_mapstat_lookup_filtered_agg.to_csv(output.tsv_mapstat_wo_pl_collapsed, sep='\t', header=True, index=False)

rule highest_id_per_species:
    """
    Per species, select the genome with the highest template ID, discard others.
    """
    input:
        tsv_mapstat_wo_pl_collapsed = rules.remove_plasmids_and_recalculate.output.tsv_mapstat_wo_pl_collapsed
    output:
        tsv_mapstat_wo_pl_collapsed_highest_id =WORKING_DIR / 'classification' / 'kma_tsv_mapstat_extended_wo_pl_collapsed_highest_id.tsv'
    run:
        tsv_mapstat = pd.read_table(input.tsv_mapstat_wo_pl_collapsed)

        tsv_mapstat_highest_id = tsv_mapstat.iloc[tsv_mapstat.groupby('species_name', sort=False)["Template_Identity"].idxmax(), :]

        tsv_mapstat_highest_id.to_csv(output.tsv_mapstat_wo_pl_collapsed_highest_id, sep='\t', header=True, index=False)

rule query_information:
    """
    Collect per query file the necessary metadata for the report 
    """
    input:
        N_BASES = rules.get_nb_of_bases.output.N_BASES,
        MAPPING_RESULTS = rules.highest_id_per_species.output.tsv_mapstat_wo_pl_collapsed_highest_id
    output:
        QUERY_INFORMATION =WORKING_DIR / 'query_information_{fdr}.json'
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        sample_yield = lambda wildcards: wildcards.sample_yield,
        fdr = lambda wildcards: wildcards.fdr,
        config= config
    run:
        import json

        with open(input.N_BASES) as f:
            n_bases = int(f.readline().strip('\n'))
        thresholds = pd.read_csv(Path(params.config["thresholds"])).set_index(["yield", "selection_strategy"])
        # Use the right threshold of the interval
        thresholds["right"] = thresholds['threshold'].str.extract(r',\s*([0-9.]+)\]').astype(float)
        thresholds = thresholds["right"].unstack(sort=False).round(2)
        if params.sample_yield != "full":
            closest_yield = params.sample_yield
            closest_threshold = thresholds.at[params.sample_yield, f'FDR {params.fdr}%']
        else:
            sample_yield_int = [int(sample_yield[:-1]) for sample_yield in thresholds.index.get_level_values("yield").unique()]
            closest_yield_int = min(sample_yield_int, key=lambda v: abs(v - n_bases / 1_000_000))
            closest_yield = f'{closest_yield_int}M'
            closest_threshold = thresholds.at[closest_yield, f'FDR {params.fdr}%']

        data = {
            "query": params.sample_name,
            "query_file_path":input.MAPPING_RESULTS,
            "n_bases": n_bases,
            "subsampling": params.sample_yield if params.sample_yield != "full" else None,
            "closest_yield": closest_yield,
            "closest_threshold": closest_threshold,
            "fdr": params.fdr,
        }

        # Write to a JSON file
        with open(output.QUERY_INFORMATION, "w") as f:
            json.dump(data, f, indent=4)



def get_query_combinations(config):
    """
    Flatten the queries dict into aligned lists for expand().
    Returns a dict: {'sample_name': [...], 'sample_yield': [...], 'fdr': [...]}
    """
    sample_name, sample_yield, fdr = [], [], []

    for name, qinfo in config["queries"].items():
        for combo in qinfo["combinations"]:
            sample_name.append(name)
            sample_yield.append(combo["subsampling"] or "full")  # handle None
            fdr.append(combo["fdr"])

    return {"sample_name": sample_name, "sample_yield": sample_yield, "fdr": fdr}

rule create_report:
    """
    Creates the output report of the workflow.
    """
    input:
        QUERY_INFORMATION = lambda wildcards: expand(Path(config['wd']) / '{sample_name}' / '{sample_yield}' / "query_information_{fdr}.json",
                                            zip,
                                            **get_query_combinations(config)

        )
    output:
        HTML = Path(config['output']['html'])
    params:
        dir_out = Path(config['output']['dir']),
        config = config
    run:
        from kaptain.app.utils import reportutils


        thresholds = pd.read_csv(Path(params.config["thresholds"])).set_index(["yield", "selection_strategy"])
        # Use the right threshold of the interval
        thresholds["right"] = thresholds['threshold'].str.extract(r',\s*([0-9.]+)\]').astype(float)
        thresholds = thresholds["right"].unstack(sort=False).round(2)


        # Initialize the report
        report = HtmlReport(Path(output.HTML), Path(params.dir_out))
        if not Path(params.dir_out).exists():
            Path(params.dir_out).mkdir(parents=True)
        report.initialize('KAPTAIN report')
        report.add_pipeline_header('KAPTAIN')

        # Add sections
        report.add_html_object(reportutils.create_analysis_info_section(params.config))
        report.add_html_object(reportutils.create_parameter_section(input.QUERY_INFORMATION))
        report.add_html_object(reportutils.create_thresholds_table(thresholds))
        report.add_html_object(reportutils.create_mapping_section(input.QUERY_INFORMATION))
        report.add_html_object(reportutils.create_citations_section())

        # Save report
        report.save()