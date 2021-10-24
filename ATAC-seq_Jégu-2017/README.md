ATAC-seq data processing pipeline
================

ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing)
is a molecular biology assay used for analysis of genome-wide chromatin
accessibility. Below is bare-bones ATAC-seq pipeline using a Unix
terminal on MacOS for acquiring and processing data from Jégu <i>et
al.</i> (2017) to prepare for downstream analysis. This pipeline largely
follows that by [Siva
Chudalayandi](https://bioinformaticsworkbook.org/dataAnalysis/ATAC-seq/ATAC_tutorial.html#gsc.tab=0),
but includes some additional detail, changes that I find personally
useful, and minor fixes. In a future update, I will compile all these
steps into a Bash script and also add my own steps on how to do the
downstream analysis.

### Installations

#### Homebrew installation

Install [Homebrew](https://brew.sh/):

``` {install-homebrew}
## ---- Unix terminal ---- ##
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

#### wget installation

Install [wget](https://www.gnu.org/software/wget/),
[FastQC](https://github.com/s-andrews/FastQC),
[SAMtools](http://www.htslib.org/),
[BEDtools](https://bedtools.readthedocs.io/en/latest/), Java Development
Kit (JDK), and [IGV](https://software.broadinstitute.org/software/igv/)
using Homebrew:

``` {wget-install}
## ---- Unix terminal ---- ##
brew install wget
brew install fastqc
brew install bowtie2
brew install samtools
brew install bedtools
# Install Java from openjdk, or:
brew install openjdk
# Alternatively install the 'official' Java JDK from Oracle:
brew install --cask oracle-jdk
# IGV requires Java to run:
brew install --cask igv
```

#### Miniconda installation

Install Miniconda. You can see the [latest installation method for
miniconda here](https://bioconda.github.io/user/install.html). Note that
Python3 will automatically be installed together with Miniconda:

``` {miniconda-install}
## ---- Unix terminal ---- ##
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh
```

Add the following line to your .zprofile and save. Use the command
‘touch \~/.zprofile’ if you don’t have one:

``` {zprofile-path-miniconda}
export PATH="/Users/username/miniconda3/bin/conda:$PATH"
```

Either close and re-start your terminal or enter the command below to
effect zprofile changes:

``` {source-zprofile}
## ---- Unix terminal ---- ##
source ~/.zprofile
```

Install relevant modules using Miniconda:

``` {conda-config}
## ---- Unix terminal ---- ##
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### Bioconda packages installation

Install
[ucsc-bigwigaverageoverbed](https://bioconda.github.io/recipes/ucsc-bigwigaverageoverbed/README.html),
[bioawk](https://anaconda.org/bioconda/bioawk),
[ucsc-bedclip](https://anaconda.org/bioconda/ucsc-bedclip), and
[ucsc-bedgraphtobigwig](https://anaconda.org/bioconda/ucsc-bedgraphtobigwig)
using Miniconda on the bioconda channel:

``` {conda-installs}
## ---- Unix terminal ---- ##
conda install -c bioconda ucsc-bigwigaverageoverbed
conda install -c bioconda bioawk
conda install -c bioconda ucsc-bedclip
conda install -c bioconda ucsc-bedgraphtobigwig
```

#### MACS2 installation

[MACS2](https://pypi.org/project/MACS2/) requires Python3 in order to
run (which was installed together with Miniconda). So we can install it
using pip3:

``` {macs2-install}
pip3 install MACS2
```

### Data acquisition and processing

Create a directory in the home folder called ‘ATAC-seq_Jégu-2017’ and
change into it:

``` {mkdir-home}
## ---- Unix terminal ---- ##
mkdir ATAC-seq_Jégu-2017
cd ATAC-seq_Jégu-2017
```

Download the experimental fasta files using wget:

``` {wget-fastq}
## ---- Unix terminal ---- ##
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR473/002/SRR4733912/SRR4733912_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR473/002/SRR4733912/SRR4733912_2.fastq.gz
```

Perform FastQC analysis on the forward and reverse reads:

``` {fastqc-analysis}
## ---- Unix terminal ---- ##
fastqc -o ~/ATAC-seq_Jégu-2017 SRR4733912_1.fastq.gz
fastqc -o ~/ATAC-seq_Jégu-2017 SRR4733912_2.fastq.gz
```

Download the <i>arapidopsis thaliana</i> genome using wget:

``` {download-genome}
## ---- Unix terminal ---- ##
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gzcat Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz > Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
```

Use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) to
build a genome index called “a_thaliana_TAIR10_build”:

``` {genome-build}
## ---- Unix terminal ---- ##
bowtie2-build --threads 12 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz a_thaliana_TAIR10_build
```

Use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) to
align the fastq files (reads) to the genome index:

``` {bowtie2-alignment}
## ---- Unix terminal ---- ##
bowtie2 --threads 12 -x a_thaliana_TAIR10_build -q -1 SRR4733912_1.fastq.gz -2 SRR4733912_2.fastq.gz -S SRR4733912.sam
```

Use [samtools](http://www.htslib.org/doc/samtools.html) view to convert
sam to bam, and pipe the output to samtools sort to sort the alignments
by leftmost coordinates:

``` {samtools-view}
## ---- Unix terminal ---- ##
samtools view --threads 12 -bS SRR4733912.sam | samtools sort --threads 12 > SRR4733912.sorted.bam
```

Use [samtools](http://www.htslib.org/doc/samtools.html) to view the
header of the sorted bam file:

``` {samtools-head}
## ---- Unix terminal ---- ##
samtools view -H SRR4733912.sorted.bam
```

The above command will produce the following output (which may change
depending on what options used previously):

``` {samtools-head-output}
## ---- Unix terminal ---- ##
@HD VN:1.0  SO:coordinate
@SQ SN:1    LN:30427671
@SQ SN:2    LN:19698289
@SQ SN:3    LN:23459830
@SQ SN:4    LN:18585056
@SQ SN:5    LN:26975502
@SQ SN:Mt   LN:366924
@SQ SN:Pt   LN:154478
@PG ID:bowtie2  PN:bowtie2  VN:2.4.4        CL:"/Users/nick/miniconda3/bin/bowtie2-align-s --wrapper basic-0 --threads 12 -x a_thaliana_TAIR10_build -q -S SRR4733912.sam -1 SRR4733912_1.fastq.gz -2 SRR4733912_2.fastq.gz"
@PG ID:samtools PN:samtools PP:bowtie2      VN:1.13 CL:samtools view --threads 12 -bS SRR4733912.sam
@PG ID:samtools.1   PN:samtools PP:samtools     VN:1.13 CL:samtools sort --threads 12
@PG ID:samtools.2   PN:samtools PP:samtools.1   VN:1.13 CL:samtools view -H SRR4733912.sorted.bam
```

Index the bam file using
[samtools](http://www.htslib.org/doc/samtools.html):

``` {index-bam}
## ---- Unix terminal ---- ##
samtools index SRR4733912.sorted.bam
```

It is important to identify and remove reads aligning to Mt and
Chloroplast in case of ATAC-seq because this can cause confounding
signals. Use [samtools](http://www.htslib.org/doc/samtools.html) do to
this:

``` {remove-mt-chloroplast-dna}
## ---- Unix terminal ---- ##
samtools idxstats SRR4733912.sorted.bam | cut -f1 | grep -v Mt | grep -v Pt | xargs samtools view --threads 7 -b SRR4733912.sorted.bam > SRR4733912.sorted.noorg.bam
```

One can compare the index stats of the BAM file without organelle DNA
alignments and the file with all alignments. You can use this to bench
mark the alignment and the procedure used to make the libraries. Use
[samtools](http://www.htslib.org/doc/samtools.html) do to this:

``` {idx-stats}
## ---- Unix terminal ---- ##
samtools idxstats SRR4733912.sorted.noorg.bam
samtools idxstats SRR4733912.sorted.bam
```

Call peaks using [MACS2](https://pypi.org/project/MACS2/):

``` {macs2-peak-calling}
## ---- Unix terminal ---- ##
macs2 callpeak -t SRR4733912.sorted.noorg.bam -q 0.05 --broad -f BAMPE -n SRR4733912 -B --trackline --outdir . &>SRR4733912.peak.log&
```

The following files are output: \* SRR4733912_control_lambda.bdg \*
SRR4733912_peaks.broadPeak \* SRR4733912_peaks.gappedPeak \*
SRR4733912_peaks.xls \* SRR4733912_treat_pileup.bdg

Make a chromosome sizes file using bioawk:

``` {chr-sizes-file}
bioawk -c fastx '{print $name, length($seq)}' Arabidopsis_thaliana.TAIR10.dna.toplevel.fa > chr.sizes
```

The contents of chromosome sizes file should look like this:

``` {chr-sizes-contents}
1   30427671
2   19698289
3   23459830
4   18585056
5   26975502
Mt  366924
Pt  154478
```

Use [bedtools
slop](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html)
to clip the bed graph files to the correct coordinates:

``` {clip-bed-files}
bedtools slop -i SRR4733912_treat_pileup.bdg -g chr.sizes -b 0 | bedClip stdin chr.sizes SRR4733912_treat_pileup.clipped.bdg
```

Sort the clipped files:

``` {sort-clipped-bed-graph}
sort -k1,1 -k2,2n SRR4733912_treat_pileup.clipped.bdg > SRR4733912_treat_pileup.clipped.sorted.bdg
```

Convert to bigwig files:

``` {bdg-to-bw}
bedGraphToBigWig SRR4733912_treat_pileup.clipped.sorted.bdg chr.sizes SRR4733912_treat_pileup.clipped.sorted.bw
```

Now we can directly view the bigwig file in the IGV Genome browser or
calculate peak scores over a genomic region/interval using
[ucsc-bigwigaverageoverbed](https://bioconda.github.io/recipes/ucsc-bigwigaverageoverbed/README.html).

Open IGV and from the top menu bar load the
‘Arabidopsis_thaliana.TAIR10.dna.toplevel.fa’ genome file: Genomes \>
Load Genome from file

Now, load the ‘SRR4733912_treat_pileup.clipped.sorted.bw’ file in IGV
using the top menu: File \> Load from File…

------------------------------------------------------------------------

### References:

Jégu, Teddy, et al. “The Arabidopsis SWI/SNF protein BAF60 mediates
seedling growth control by modulating DNA accessibility.” Genome biology
18.1 (2017): 114.
