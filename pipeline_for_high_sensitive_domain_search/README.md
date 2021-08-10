
# pipeline_for_high_sensitive_domain_search - a helper tool for high sensitive HMM-HMM search against Pfam using HHsearch and JackHMMER

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)
[![size](https://img.shields.io/github/size/webcaetano/craft/build/phaser-craft.min.js.svg)]()

## requirements
* jackhmmer (HMMER3)
* hhsearch, reformat.pl (HHsuite3)
* Ruby (ver >=2.1)

## usage 
```
### pipeline_for_high_sensitive_domain_search ver 0.1.3 (2019-07-29) ###

[description]
pipeline_for_high_sensitive_domain_search - a helper tool for high sensitive HMM-HMM search against Pfam using HHsearch and JackHMMER

[reference]
HHsearch  - doi:10.1093/bioinformatics/bti125
JackHMMER - doi:10.1186/1471-2105-11-431

[usage]
$ pipeline_for_high_sensitive_domain_search [options]

[dependencies]
    - jackhmmer             -- included in the HMMER3 package;
                               http://hmmer.org/
    - hhsearch, reformat.pl -- included in the HHsuite package;
                               HHsuite3 - https://github.com/soedinglab/hh-suite
                               HHsuite2 - http://wwwuser.gwdg.de/~compbiol/data/hhsuite/releases/all/
    - ruby (ver >=2.1)

[options]
  (general)
    -i, --in      [path]  (required) query protein fasta file. To use multiple files as input, this value could be (1) comma-separeted file paths or (2) path of input file list (one file in one line).
                          If many proteins are included in a fasta file, adding "belongTo:" tag in comment lines to generate easy-to-understand output.  e.g. ">protein1 len:200 belongTo:contig1"
    -o, --out     [path]  (required) output directory. It should not exist.
    -J, --jackdb  [path]  (required) protein fasta file to create query HMM with query proteins, used as database sequences for JackHMMER
    -H, --hhdb    [path]  (required) Pfam database for HHsearch available here: http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
    -h, --help
    -v, --version

  (jackhmmer related)
    -N, --iter    [int]   (default: 5) number of iteration of JackHMMER serach
        --incE    [num]   (default: 0.001) inclusion threshold for sequences
        --incdomE [num]   (default: 0.001) inclusion threshold for domains

  (hhsearch related)
    -n, --nreport [int]   (default: 3)  number of hits for each query shown in a result table
        --evalue  [num]   (default: 10) E-value threshold for inclusion in a result table
        --pvalue  [num]   (default: 1)  P-value threshold for inclusion in a result table
        --prob    [num]   (default: 30) %-Probability threshold for inclusion in a result table

  (use GNU parallel)
        --ncpus   [int]   -- number of simaltaneous jobs

  (for icr user)          for computation in the ICR supercomputer system
        --queue   [JP1]   -- queue to calculate

[output files]
result/<input file name>/<group>/besthit.tsv   -- table of HHserach besthit. See README.md for column description. Group is what is specified by the "belongTo:" tags or "ungrouped" if not provided (see description of the -i/--in option).
result/<input file name>/<group>/top<n>hit.tsv -- table of HHserach top<n>hit. <n> can be altered by -n/--nreport option (default: 3).
work/<input file name>/<group>/<protein>/      -- The directory includes all the output files of jackhmmer and hhsearch
info/stats.txt                                 -- number of proteins in each file
```

## columns of result table (besthit.tsv / top\<n\>hit.tsv)
```
1  protein            - query protein name
2  query_len          - query protein length (same as length of query HMM)
3  query_num_iter     - the number of iteration by jackhmmer
4  query_num_seq      - the number of sequences included in query HMM
5  query_num_seq_used - the number of sequences included in query HMM that passed through hhsearch prefilter.
6  template_acc       - template Pfam accession
7  template_name      - template Pfam name
8  template_desc      - template Pfam description
9  probability        - reported probability
10 e-value            - reported e-value
11 p-value            - reported p-value
12 score              - reported score
13 cols               - length of the HMM-HMM alignment
14 query_pos          - position of query HMM included in the HMM-HMM alignment
15 template_pos       - position of template HMM included in the HMM-HMM alignment
16 template_len       - length of template HMM
17 rank               - rank of HMM-HMM alignment repored by hhsearch
```
Columns 5 to 17 are derived from hhsearch output (.hhr file).

## information
Old name: HHliner.
Prototype of this package was used in environmental virus studies below.
```
- Nishimura, Y. et al. Environmental viral genomes shed new light on virus–host interactions in the ocean.
  mSphere, 2, e00359-16, doi:10.1128/mSphere.00359-16 (2017)
- Yoshida, T. et al. Locality and diel cycling of viral production revealed by a 24 h time course cross-omics analysis in a coastal region of Japan.
  ISME J., doi:10.1038/s41396-018-0052-x (2018)
```
