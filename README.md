# mclUMI
![](https://img.shields.io/badge/mclUMI-executable-519dd9.svg)
![](https://img.shields.io/badge/last_released-Oct._2021-green.svg)
![](https://img.shields.io/github/stars/cribbslab/mclumi?logo=GitHub&color=blue)
![](https://img.shields.io/pypi/v/mclumix?logo=PyPI)
[![Documentation Status](https://readthedocs.org/projects/mclumi/badge/?version=latest)](https://mclumi.readthedocs.io/en/latest/?badge=latest)
[![Downloads](https://pepy.tech/badge/mclumi)](https://pepy.tech/project/mclumi)

###### tags: `UMI deduplication` `PCR deduplication` `scRNA-seq` `bulk-RNA-seq`

## Overview
```
 __  __  ____ _    _   _ __  __ ___   _____           _ _    _ _   
|  \/  |/ ___| |  | | | |  \/  |_ _| |_   _|__   ___ | | | _(_) |_ 
| |\/| | |   | |  | | | | |\/| || |    | |/ _ \ / _ \| | |/ / | __|
| |  | | |___| |__| |_| | |  | || |    | | (_) | (_) | |   <| | |_ 
|_|  |_|\____|_____\___/|_|  |_|___|   |_|\___/ \___/|_|_|\_\_|\__|
```

This repository deposits the mclUMI toolkit developed by Markov clustering (MCL) network-based algorithms for precisely localizing unique UMIs and thus removing PCR duplicates. mclUMI enables a construction of sub-graphs with UMI nodes to be relatively strongly connected.

## Documentation
The API documentation of mclUMI is available at https://mclumi.herokuapp.com and https://mclumi.readthedocs.io/en/latest.

## System requirement
Linux or Mac

## Installation
We tested the software installation on a Linux system, which has the following configuration:
* Distributor ID: Ubuntu
* Description:    Ubuntu 20.04.3
* Release:        20.04
* Codename:       focal

The anaconda is configured as:
* Conda version: 4.11.0

> You can use `conda update conda` and `conda update anaconda` to keep your anaconda up-to-date.

We recommend using a `Python` of version **`3.9.1`** as the base python to create your conda environment because `NumPy` and `Pandas` in a `Python` of higher version `3.9` may require a few dependencies that are not included in the installation of mclUMI or make conflicts with existing packages.

**Step 1**: create a conda environment, e.g., mclumi
  ```angular2html
  conda create --name mclumi python=3.9.1
      
  conda activate mclumi
  ```
  
  <h1>
      <img src="https://github.com/cribbslab/mclumi/blob/main/imgs/conda-setting.png?raw=true">
      <br>
  </h1>

**Step 2**: sourced from https://pypi.org/project/mclumix.
  ```angular2html
  pip install --upgrade mclumix
  ```
After a two-step installation procedure, you should see the following outputs.
  <h1>
      <img src="https://github.com/cribbslab/mclumi/blob/main/imgs/install.png?raw=true">
      <br>
  </h1>

## Usage
To ease the use of mclUMI for multiple groups of users, we have made it usable in both command-line interface (CLI) and inline mode. 

### 1. CLI
1.1 Parameter illustration

By typing `mclumi -h`, you are able to see the package usage as shown below.

```
usage: mclumi [-h] [--read_structure read_structure] [--lens lens]
              [--input input] [--output output] [--method method]
              [--input_bam input_bam] [--edit_dist edit dist]
              [--inflation_value inflation_value]
              [--expansion_value expansion_value]
              [--iteration_number iteration_number]
              [--mcl_fold_thres mcl_fold_thres] [--is_sv is_sv]
              [--output_bam output_bam] [--verbose verbose]
              [--pos_tag pos_tag] [--gene_assigned_tag gene_assigned_tag]
              [--gene_is_assigned_tag gene_is_assigned_tag]
              tool

Welcome to the mclumi toolkit

positional arguments:
  tool                  trim, dedup_basic, dedup_pos, dedup_gene, dedup_sc

optional arguments:
  -h, --help            show this help message and exit
  --read_structure read_structure, -rs read_structure
                        str - the read structure with elements in conjunction
                        with +, e.g., primer_1+umi_1+seq_1+umi_2+primer_2
  --lens lens, -l lens  str - lengths of all sub-structures separated by +,
                        e.g., 20+10+40+10+20 if the read structure is
                        primer_1+umi_1+seq_1+umi_2+primer_2
  --input input, -i input
                        str - input a fastq file in gz format for trimming
                        UMIs
  --output output, -o output
                        str - output a UMI-trimmed fastq file in gz format.
  --method method, -m method
                        str - a dedup method: unique | cluster | adjacency |
                        directional | mcl | mcl_ed | mcl_val
  --input_bam input_bam, -ibam input_bam
                        str - input a bam file curated by requirements of
                        different dedup modules: dedup_basic, dedup_pos,
                        dedup_gene, dedup_sc
  --edit_dist edit dist, -ed edit dist
                        int - an edit distance used for building graphs at a
                        range of [1, l) where l is the length of a UMI
  --inflation_value inflation_value, -infv inflation_value
                        float - an inflation value for MCL, 2.0 by default
  --expansion_value expansion_value, -expv expansion_value
                        int - an expansion value for MCL at a range of (1,
                        +inf), 2 by default
  --iteration_number iteration_number, -itern iteration_number
                        int - iteration number for MCL at a range of (1,
                        +inf), 100 by default
  --mcl_fold_thres mcl_fold_thres, -fthres mcl_fold_thres
                        float - a fold threshold for MCL at a range of (1, l)
                        where l is the length of a UMI.
  --is_sv is_sv, -issv is_sv
                        bool - to make sure if the deduplicated reads writes
                        to a bam file (True by default or False)
  --output_bam output_bam, -obam output_bam
                        str - output UMI-deduplicated summary statistics to a
                        txt file.
  --verbose verbose, -vb verbose
                        bool - to enable if output logs are on console (True
                        by default or False)
  --pos_tag pos_tag, -pt pos_tag
                        str - to enable deduplication on the position tags (PO
                        recommended when your bam is tagged)
  --gene_assigned_tag gene_assigned_tag, -gt gene_assigned_tag
                        str - to enable deduplication on the gene tag (XT
                        recommended)
  --gene_is_assigned_tag gene_is_assigned_tag, -gist gene_is_assigned_tag
                        str - to check if reads are assigned the gene tag (XS
                        recommended)
```

1.2 Example commands

* extracting and attaching umis to names of reads in fastq format
    ```
    mclumi trim -i ./pcr_1.fastq.gz -o ./pcr_trimmed.fastq.gz -rs primer_1+umi_1+seq_1+umi_2+primer_2 -l 20+10+40+10+20
    ```

* deduplication on only one genome position 
    ```
    mclumi dedup_basic -m mcl -ed 1 -infv 1.6 -expv 2 -ibam ./example_bundle.bam -obam ./dedup.bam
    ```

* deduplication per genome position
    ```
   mclumi dedup_pos -m mcl -pt PO -ed 1 -infv 1.6 -expv 2 -ibam ./example_bundle.bam -obam ./basic/dedup.bam
    ```

* deduplication per gene (applicable to bulk RNA-seq data)
    ```
    mclumi dedup_gene -m directional -gt XT -gist XS -ed 1 -ibam ./hgmm_100_STAR_FC_sorted.bam -obam ./dedup.bam
    ```

* deduplication per cell per gene (applicable to single-cell RNA-seq data)
    ```
    mclumi dedup_sc -m directional -gt XT -gist XS -ed 1 -ibam ./hgmm_100_STAR_FC_sorted.bam -obam ./dedup.bam
    ```

### 2. Inline

see Jupyter notebooks
```
./notebooks/
```

## Output
see `./notebooks/results_spelt_out.ipynb` for result format. More types of output format are about to be added.

## Contact
Homepage: https://www.ndorms.ox.ac.uk/team/adam-cribbs  