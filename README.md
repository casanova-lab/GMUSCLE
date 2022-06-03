# GMUSCLE
Genotyping MUltiplexed-Sequencing of CRISPR-Localized Editing

## Introduction
- CRISPR/Cas9 gene editing is widely used to create genetic modified cells. The genetic editing results in cellular heterogeneity, allelic diversity, and  randomness of de novo sequences, which constitute a major challenge in genotyping CRISPR/Cas9-edited cells. 
- We introduced a streamlined and cost-effective protocol along with the computational tool GMUSCLE, to sequence and genotype the products of CRISPR/Cas9 gene editing with great depth, high accuracy, high efficiency, and handy utilization.
- GMUSCLE enables the quantitative identification and qualitative visualization of the major genotypes from demultiplexed sequencing data of CRISPR/Cas9-edited cells.

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
- demutiplexed fastq/fastq.gz files for all samples
- reference genome
- genomic region of interest (CHR, START, END)

### Output
- a summary of read counts and major genotypes for each sample
- a plot of major genotypes (deletions/insertions of different sizes) and their frequencies on genomic positions, for each sample
- a matrix of major genotypes and their frequencies for all samples

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
*-dir*|str|directory of the fastq files|N.A.
*-sample*|filename|file of sample list|N.A.
*-readtype*|str|single-end (se) or paired-end (pe)|se
*-genome*|filename|file of reference genome|N.A.
*-chr*|str|chromosome (consistent format with the reference genome)|N.A.
*-start*|int|start position of the sequencing region|N.A.
*-end*|int|end position of the sequencing region|N.A.
*-count_cutoff*|int|read count cutoff|30
*-ratio_cutoff*|float|frequency cutoff for being a major genotype|0.01(1%)
*-prefix*|str|prefix for the output files|N.A.

*Notes:*
- *if sample file and genome file are not in the same folder with GMUSCLE.py, include their path in the command*
- *if paired-end, the fastq files have to be suffixed with R1/R2*

## References
- *Zhang P#, Yang R#, Abel A, Casanova J-L.* Genotyping MUltiplexed-Sequencing of CRISPR-Localized Editing (GMUSCLE): 
a computational approach to analyze CRISPR-edited single cell clones. (2022)

## Contact
> **Developer:** Peng Zhang, Ph.D.

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
