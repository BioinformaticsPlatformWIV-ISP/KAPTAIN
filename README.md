# KAPTAIN

KAPTAIN (KMA-bAsed Pipeline for meTAgenomic specIes ideNtification) is an optimized taxonomic classification workflow. It uses [KMA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2336-6) 
as its core tool and adds the following post-processing steps:
- Remove plasmid hits
- Recalculate mapping metrics in light of full genomes
- Select per species the genome with the highest template identity

The template ID is a metric provided by KMA, calculated based on the number of identical bases between the consensus sequence and the template (i.e., the reference) divided by the template’s length: 

**Template ID** = (# identical bases between consensus and template) / (template length)  

in which the consensus sequence reflects the majority vote of the aligned reads.
This way, the template ID is determined by two factors: the breadth of reference coverage and the ID of the covered part.  

A template identity threshold is set according to the user-specified precision, and it generally achieves that precision based on data from ten prior broad samples.
These thresholds are indicative, and all results are included in the final report.

### KAPTAIN is also available on our public [Galaxy instance](https://galaxy.sciensano.be/) (registration required).

## Installation

The KAPTAIN workflow has the following dependencies:
- [Rasusa 2.1.1](https://github.com/mbhall88/rasusa)
- [seqtk 1.4](https://https://github.com/lh3/seqtk)
- [KMA 1.6.2](https://bitbucket.org/genomicepidemiology/kma/src/master/README.md)

The corresponding binaries should be in your PATH to run the workflow. 
Other versions of these tools may work, but have not been tested.

To install, first download the repository and navigate to the directory.
```bash
git clone https://github.com/BioinformaticsPlatformWIV-ISP/KAPTAIN.git
cd KAPTAIN/
```
Optionally create a virtual environment:
```bash
virtualenv kaptain_env --python=python3.10;
. kaptain_env/bin/activate;
```
Once ready, you can install the package with Python >=3.10: 
```bash
pip install . 
```

## Database

In order for the tool to work, a database is needed (`--db`).  
The default database can be downloaded through one of the two options, alongside a lookup file with metadata of the entries in the database (`--db-lookup`):
- Zenodo: https://doi.org/10.5281/zenodo.17435875
- FTP: https://bioit-ftp.sciensano.be/kma/

This database contains (complete) genomes.
See [Custom database](#custom-database) for information on how to provide a custom database.
<!-- TODO provide link and command once FTP is online -->

## USAGE

```bash
usage: kaptain [-h] --ont-in ONT_IN [ONT_IN ...] --db DB --db-lookup DB_LOOKUP [--dir-working DIR_WORKING] --output OUTPUT [--output-html OUTPUT_HTML] [--subsampling {200M,500M,1000M,1500M,2000M,None} [{200M,500M,1000M,1500M,2000M,None} ...]] [--fdr {15,10,5,1} [{15,10,5,1} ...]]
               [--threads THREADS] [--version] [--log]


options:
  -h, --help            show this help message and exit
  --ont-in ONT_IN [ONT_IN ...]
                        ONT input FASTA/Q file (default: None)
  --db DB               Prefix of KMA database. Database exists of four files named *.{comp.b, length.b, name, seq.b} (default: None)
  --db-lookup DB_LOOKUP
                        Lookup file of KMA database (default: None)
  --dir-working DIR_WORKING
                        Working directory (default: /home/alvanuffelen/ProjectsAlvan/AVU_KAPTAIN/kaptain)
  --output OUTPUT       Output directory (default: None)
  --output-html OUTPUT_HTML
                        Output report name (default: report.html)
  --subsampling {200M,500M,1000M,1500M,2000M,None} [{200M,500M,1000M,1500M,2000M,None} ...]
                        Subsample input to number of bases before classification. Leave empty or use None for no downsampling. (default: [None])
  --fdr {15,10,5,1} [{15,10,5,1} ...]
                        FDR setting. (default: [5])
  --threads THREADS     Number of threads (default: 4)
  --version             Print version and exit
  --log                 Write out log information to file (default: False)
```

### Usage examples

1. Basic Classification:
```
kaptain --ont-in query.fq --db my_database --db_lookup my_database.lookup
        --output results/ --subsampling 500M --fdr 5
```

2. Multiple Queries:
```
kaptain --ont-in query.fq query2.fq --db my_database --db_lookup my_database.lookup
        --output results/ --fdr 15
```

3. Multiple FDR and Subsampling Settings:
```
kaptain --ont-in query.fq query.fq query.fq --db my_database --db_lookup my_database.lookup
        --output results/ --subsampling 500M 500M None --fdr 15 10 1


   This example will run the queries as:
    | Query     | Subsample | FDR |
    |-----------|-----------|-----|
    | query.fq  | 500M      | 15  |
    | query.fq  | 500M      | 10  |
    | query2.fq | None      | 1   |

   Both --fdr and --subsampling must either have:
     - a single value (applied to all queries), or
     - as many values as there are query files (applied pairwise).
```

## Custom database
A custom database should contain (complete) genomes.
KMA needs an indexed version (`--db`) of your custom database.  
To build an index from your reference FASTA/Q:
```
kma index -i database.fasta -o kma_db
```
This will make four files: `kma_db.length.b`, `kma_db.comp.b`, `kma_db.name` and `kma_db.seq.b`.  
Additionally, a lookup file is needed (`--db_lookup`) with metadata on each sequence in the database.  
The file should be without headers and tab-separated, that contains the following columns:

| Column        | Description                                          |
|---------------|------------------------------------------------------|
| Header ID\*   | Header ID of sequence in the database                |
| Length        | Length sequence                                      |
| Accession     | Accession of genome/assembly the sequence belongs to |
| Sequence Type | Either 'Chromosome' or 'Plasmid'                     |
| Species Name  | E.g., Escherichia coli                               |

\* Note that the header ID is used, and not the full header line. 
The header ID are all characters before the first space. 
**Each header ID should be unique**.

## CONTACT
[Create an issue](https://github.com/BioinformaticsPlatformWIV-ISP/KAPTAIN/issues) to report bugs, propose new functions or ask for help.

## CITATION
TBD

-----

Copyright - 2025 Alexander Van Uffelen <alexander.vanuffelen@sciensano.be>
