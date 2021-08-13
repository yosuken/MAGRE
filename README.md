
# MAGRE - a tool for MAG quality REfinement by removing contigs that are likely 'contamination' of the MAG.

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)
[![size](https://img.shields.io/github/size/webcaetano/craft/build/phaser-craft.min.js.svg)]()

MAGRE is a tool for quality refinement of prokaryotic MAGs by removing contigs that are likely 'contamination' of the MAG.
Decontamination is performed based on three types of filters (these could be divided into six categories), without depending on single-copy marker genes.
The original form of this tool was developed for quality improvement of a large-scale ocean MAG collection (Nishimura and Yoshizawa, bioRxiv, 2021).
Later, the original scripts were organized as an easy-to-use tool, extended with some additional function (e.g., user can select categories of filter out condition).
For details, see original publication (Nishimura and Yoshizawa, bioRxiv, 2021)

## install (use conda environment)

### [1] make conda environment and install packages
```
$ conda create -n MAGRE -y && conda activate MAGRE
```

### [2] install packages: other versions can be used, but virsorter should be v1.0.6
```
$ conda install -y -c conda-forge r-base=4.1.0 ruby=2.7.2 parallel=20210622
$ conda install -y -c bioconda virsorter=1.0.6 prodigal=2.6.3 hhsuite=3.3.0 cat=5.0.3 blast=2.12.0
$ conda install -y -c biobuilds fasta=36.3.8e
```

### [3] make copy reformat.pl (in hhsuite) to bin directory
- File path (i.e. ~/miniconda) should be modified to fit your environment.
```
$ cp ~/miniconda/envs/MAGRE/scripts/reformat.pl ~/miniconda/envs/MAGRE/bin
```

### [4] download virsorter data (reference: https://github.com/simroux/VirSorter)
- First, move to a directory to put virsorter data
```
$ wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
$ md5sum virsorter-data-v2.tar.gz # md5sum should return dd12af7d13da0a85df0a9106e9346b45
$ tar -xvzf virsorter-data-v2.tar.gz # --> generate 'virsorter-data' directory
```

### [5] download pfam hhsuite db (reference: https://github.com/soedinglab/hh-suite)
- First, move to a directory to put hhsuite db
- download pfamA db (e.g., pfamA_32.0.tar.gz) from http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs
```
$ wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_32.0.tar.gz
$ tar xvf pfamA_32.0.tar.gz # --> generate 7 files: ./pfam{_a3m.ffindex,_a3m.ffdata,_hhm.ffindex,_hhm.ffdata,_cs219.ffindex,_cs219.ffdata,.md5sum}
```

### [6] set up CAT/BAT database (reference: https://github.com/dutilh/CAT)
- First, move to a directory to put CAT database
```
$ CAT prepare --help
$ CAT prepare --fresh -d CAT_database/ -t CAT_taxonomy/
```

### Be careful for diamond version
diamond v0.9.14 will be installed to the environment, due to specification of virsorter=1.0.6.
When the execution, CAT will use this version of diamond.
Please confirm to use the same version of diamond that builds a CAT database. (see https://github.com/dutilh/CAT)
If you have already built a database with a different verion of diamond,
a possible solution could be overwriting the diamond binary within conda environment (see followings).
A straight-ahead solution like 'conda install diamond=x.x.x' will be failed because of conflict of dependencies.
Alternatively, a bolder approach can be taken as workaround. Here, suppose you need diamond v0.9.29.
File path below (i.e. ~/miniconda) should be modified to fit your environment.
```
$ conda create -n diamond_v0.9.29 -y && conda activate diamond_v0.9.29 && conda install -y -c bioconda diamond=0.9.29
$ cp ~/miniconda/envs/MAGRE/bin/diamond ~/miniconda/envs/MAGRE/bin/diamond.original   ## make backup
$ cp ~/miniconda/envs/diamond_v0.9.29/bin/diamond ~/miniconda/envs/MAGRE/bin/diamond  ## replace
```

### [7] modify db_config.txt
- edit 4 file paths of databases

### [8] test run
- Double quatation (") is needed to use * as a wildcard for specifying input files.
```
$ ./MAGRE -n 12 -i "data/testdata/*.fa" -c db_config.txt -o test_out --fcov data/testdata/depth_in_jgi_format.txt --overwrite
```

## requirement
- ruby (ver >= 2.0)
- GNU parallel
- wrapper_phage_contigs_sorter_iPlant.pl (virsorter, v1.0.6)
- prodigal    (prodigal)
- blast       (blast)
- ssearch36   (fasta)
- jackhmmer   (hmmer)
- hhsearch    (hhsuite)
- reformat.pl (hhsuite)
- ccfind      (bundled in this package. https://github.com/yosuken/ccfind)
- pipeline_for_high_sensitive_domain_search  (bundled in this package. https://github.com/yosuken/pipeline_for_high_sensitive_domain_search)

## filters and categories
Currently, three filters and coresponding six categories are available. Users can choose categories used to remove contigs.
For details of each category, see the publication (Nishimura and Yoshizawa, bioRxiv, 2021)

- taxonomy filter (remove taxonomically inconsistent contigs, based on taxonomic assignment of CAT/BAT)
  - taxonomy: remove taxonomically inconsistent contigs based on CAT/BAT result. See the original publication for datails.
- outlier filter (remove outlier contigs in terms of coverage and tetranucleotide composition)
  - coverage: remove outlier of coverage. If PC1 of PCA analysis is <-2.5 or >2.5, contigs are removed.
  - tetranuc: remove outlier of tetranucleotide composition. If PC1 of PCA analysis is <-2.5 or >2.5, contigs are removed.
- mobile element filter (viral contigs detected by virsorter or high sensitive terL gene identification, and circular contigs)
  - virsorter: remove predicted viral contigs based on virsorter result. This category is applied to >=3kb contigs. 
  - terL: remove predicted viral contigs based on the terL gene detection through pipeline_for_high_sensitive_domain_search. This category is applied to <10kb contigs.
  - circular: remove predicted virus or plasmid based on the circular contig detection through ccfind.

## usage 
```
### MAGRE ver 0.1.0 (2021-08-11) ###

[description]
MAGRE is a tool for MAG quality refinement, by automatic removal of contigs that are likely 'contamination' of a given MAG.
The original form of this tool was developed to improve quality of genomes that are reconstructed from a large-scale ocean metagenome collection (Nishimura and Yoshizawa, bioRxiv, 2021).
The original scripts were here re-organized as an easy-to-use tool, with some additional functions (e.g., user can select categories to filter out).
This tool is applicable for any kind of prokaryotic genomes (MAGs/SAGs/isolates).
More functions (e.g., choose thresholds for each category) are planned to develop.

[filters]
1. taxonomy filter (remove taxonomically inconsistent contigs, based on taxonomic assignment of CAT/BAT)
   category for this filter is 'taxonomy'
2. outlier filter (remove outlier contigs in terms of coverage and tetranucleotide composition)
   categories for this filter are 'coverage' and 'tetranuc'
3. mobile element filter (viral contigs detected by virsorter or high sensitive terL gene identification, and circular contigs)
   categories for this filter are 'virsorter', 'terL', and 'circular'

[categories]
1.   taxonomy  - remove taxonomically inconsistent contigs based on CAT/BAT result. See the original publication for datails.
2-a. coverage  - remove outlier of coverage. If PC1 of PCA analysis is <-2.5 or >2.5, contigs are removed.
2-b. tetranuc  - remove outlier of tetranucleotide composition. If PC1 of PCA analysis is <-2.5 or >2.5, contigs are removed.
3-a. virsorter - remove predicted viral contigs. This category is applied to >=3kb contigs.
3-b. terL      - remove predicted viral contigs. This category is applied to <10kb contigs.
3-c. circular  - remove predicted virus or plasmid.

[usage]
$ MAGRE [options] -q <genome fasta file(s)> -c <configuration file> -o <output dir>

[options]
  (general)
    -h, --help
    -v, --version

  (file/directory)
    -i, --input      [file(s)] (required)  -- genomic fasta file(s)
                                              Multiple genome files can be specfied with wildcard (e.g., 'dir/*.fa', quote required) or comma separated values (e.g., A.fa,B.fa).
    -c, --config     [file] (required)     -- configuration file (two columns of space separated values; template: db_config.txt)
    -o, --outdir     [path] (required)     -- output directory (should not exist unless '--overwrite' is specified)
    --fcov           [file]                -- coverage of contigs (the format must be either 'TSV with a header line' or 'output of jgi_summarize_bam_contig_depths')
    --category       [str] (default:all)   -- selction of categories to remove contigs
                                              It must be combination of 'taxonomy', 'coverage', 'tetranuc', 'virsorter', 'terL', 'circular'
                                              e.g: 'taxonomy,coverage' will remove taxonomically inconsistent contigs and coverage outlier contigs.
    --overwrite      [bool]                -- overwrite output directory

  (computation)
    -n, --ncpus      [int] (default: 4)    -- the number of CPUs used for calculation

[output files]
  result/<name>/contig_table.tsv -- result of all analysis. Line 1: header. Line 2: result of input MAG. Line 3-: result of contigs.
  result/<name>/refined.fa       -- genomic fasta of decontaminated MAG
  result/<name>/refined.log      -- simple statistics of decontamination. It is recommended to check the percentage of retained fraction. If the percentage is too small (e.g., <50%), it is better to inspect 'contig_table.tsv' and do manual curation.
```

## format of the output TSV file (contig_table.tsv)
- column 1-5:   genome statistics
- column 6-8:   coding potential based on prodigal
- column 9-13:  CAT/BAT output
- column 14-15: polished lineage (after trimming trivial parts of CAT/BAT predicted lineage; see Nishimura and Yoshizawa for details)
- column 16-17: PC1 value of PCA analysis of coverage and tetranucleotide composition
- column 18:    circular contigs predicted by ccfind
- column 19-23: virsorter result (analysis performed only for >=3kb contigs)
- column 24-27: terL detection by pipeline_for_high_sensitive_domain_search (analysis performed only for <10kb contigs)
- column 28-30: flags of decontamination

## note
- To the given CAT_taxonomy database, subdirectory 'MAGRE' and 'parsed2.txt' within the directory will be created on the first run of MAGRE.
- In the original publication (Nishimura and Yoshizawa, bioRxiv, 2021), pre-screening was performed for the 'terL' category using HMMs built from proteins of aquatic (mainly marine) viral MAGs (EVGs). However, for generalization of the tool (considering the other environments, such as human gut), this pre-screening is not implemented. All genes of the given MAGs were searched for terL using pipeline_for_high_sensitive_domain_search.
- If removal of 'circular' contigs is activated (default setting), large circular contigs are removed even if the contig represents a whole chromosome. I am planing to make a new function to set a length threshold (or some other constraint) on removal of circular contigs. As a workaround, it is recommended to check the percentage of retained fraction. If the percentage is too small (e.g., <50%), it is better to inspect 'contig_table.tsv' and do manual curation.
- If removal of 'coverage' and/or 'tetranuc' categories is activated (default setting), large contigs are removed in some cases. For example, suppose a MAG contains a few long contigs (>10kb; possibly derived from organism A) and multiple short contigs (<1kb; possibly derived from organism B). If the number of the short contigs are larger than that of the long contigs, PCA analysis will regard the long contigs as outlier. I am planning to manipulate such a case, but currently it is recommended to check the percentage of retained fraction and better to inspect 'contig_table.tsv' and do manual curation.
