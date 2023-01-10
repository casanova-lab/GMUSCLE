# GMUSCLE
Genotyping MUltiplexed-Sequencing of CRISPR-Localized Editing

## Introduction
- CRISPR/Cas9 gene editing is widely used to create genetic modified cells. The genetic editing results in cellular heterogeneity, allelic diversity, and  randomness of de novo sequences, which constitute a major challenge in genotyping CRISPR/Cas9-edited cells. 
- We introduced a streamlined and cost-effective protocol along with the computational tool GMUSCLE, to sequence and genotype the products of CRISPR/Cas9 gene editing with great depth, high accuracy, high efficiency, and handy utilization.
- GMUSCLE enables the quantitative and qualitative identification of the major genotypes from sequencing data of CRISPR/Cas9-edited cells.

## Usage
### Dependency and quick installation
GMUSCLE is written in [python3](https://www.python.org/downloads/), and requires python packages [biopython](https://biopython.org/), [pandas](https://pypi.org/project/pandas/), and [plotnine](https://plotnine.readthedocs.io/en/stable/) installed. 
```
pip install biopython
pip install pandas
pip install plotnine
```

The standalone [blast+](https://www.ncbi.nlm.nih.gov/books/NBK569861/) is also needed. Please download it from [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/).
```
tar zxvpf ncbi-blast-xxx-linux.tar.gz
export PATH=$PATH:$HOME/ncbi-xxx/bin
```

### Input
- fastq/fastq.gz files for all samples (demutiplexed)
- reference genome
- genomic region of interest (CHR, START, END)

### Output
- a summary file of the read counts and the major genotypes (in VCF format) detected in each sample
- a plot of major genotypes (deletions/insertions (indels) of different sizes, at different genomic positions) and their proportions in each sample
- a sequence alignment of all major genotypes against the wild-type sequence, at nucleotide level, for all samples
- a matrix of all major genotypes and their proportions in all samples

### Example
We provided four example fastq files from CRISPR-editing experiments that targeted the 3'-end of exon 5 of *IFNAR1*. The genomic region was chr21-34715869-34716164 (GRCh37).

### Command
```
python GMUSCLE.py \
--dir ./fastq/ \
--sample test_sample_list.txt \
--readtype se \
--genome human_GRCh37_reference_genome.fa \
--chr chr21 \
--start 34715869 \
--end 34716164 \
--count_cutoff 30 \
--ratio_cutoff 0.01 \
--prefix test_
```

### Parameters
Parameter | Type | Description | Default
----------|------|-------------|--------------
*--dir*|str|directory of the fastq files|N.A.
*--sample*|filename|file of sample list|N.A.
*--readtype*|str|single-end (se) or paired-end (pe)|se
*--genome*|filename|file of reference genome|N.A.
*--chr*|str|chromosome (consistent format with the reference genome)|N.A.
*--start*|int|start position of the sequencing region|N.A.
*--end*|int|end position of the sequencing region|N.A.
*--count_cutoff*|int|read count cutoff|30
*--ratio_cutoff*|float|frequency cutoff for being a major genotype|0.01 (1%)
*--prefix*|str|prefix for the output files|N.A.

*Notes:*
- *if sample file and genome file are not in the same folder with GMUSCLE.py, include their path in the command*
- *if paired-end, the fastq files have to be suffixed with _R1/R2*

## Extended Utility
Besides the multiplexed-sequencing ability of this protocol, GMUSCLE software alone is also versatile. GMUSCLE is able to analyze the sequencing data:
- from bulk cell populations in organs/tissues, for which the users just need to provide the path of the large fastq files to GMUSCLE; 
- from cells that were treated with multiple sgRNA to edit multiple target sites, for which the users just need to run GMUSCLE command multiple times with different genomic positions of the target sites; 
- from different experimental gene-editing protocols/systems, as GMUSCLE only needs sequencing data, reference genome, and genomic position from the users for its analysis;
- from gene-edited products in other organisms, for which the users only need to provide the organism’s reference genome to run GMUSCLE.

## References
- *Zhang P, Yang R, Abel A, Casanova J-L.* Genotyping MUltiplexed-Sequencing of CRISPR-Localized Editing (GMUSCLE): 
an experimental and computational approach to analyze CRISPR-edited cells. (2022)

## Contact
> **Developer:** Peng Zhang, Ph.D.

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
