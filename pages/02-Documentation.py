import streamlit as st
st.set_page_config(
    page_title="Adaptive Sampler Check Documentation",
    page_icon=":scissors:",
    layout="wide",
)

st.write("""
# Adaptive Sampler Check Documentation

## Table of Contents
        
1. [Introduction](#introduction)
2. [Requirements](#requirements)
3. [(In)frequently Asked Questions](#in-frequently-asked-questions)
    1. [How do I generate a BED file?](#how-do-i-generate-a-bed-file)
    2. [How do I use custom reference genomes?](#how-do-i-use-custom-reference-genomes)
    3. [I am getting warnings, what do I do?](#i-am-getting-warnings-what-do-i-do)
    4. [What is an AGP file? And do I need it?](#what-is-an-agp-file-and-do-i-need-it)
    5. [What did you do to my BED file?](#what-did-you-do-to-my-bed-file)


## Introduction

### What is Adaptive Sampler Check?

This is a tool to help you design and validate ROIs (Regions of Interest) for your [adaptive sampling](https://archive.is/wy7eG) experiments to be 
using with the [Oxford Nanopore](https://nanoporetech.com/) range of sequencers.

Some of the features include:
- Easy to use web interface served via [Streamlit](https://streamlit.io/)
- Checking formatting errors of the BED file
- Use a preset reference genome or upload your own, to the BED entries
- Check for overlapping regions
- Adjust the buffer size of regions to optimise the enrichment of the regions
- Generate a custom fasta reference file, via the AGP file format
- Rich in metadata, to help you keep track of your experiments


### Requirements

1. A [BED file](https://grch37.ensembl.org/info/website/upload/bed.html) with your regions of interest. This file should have only contain the first 3 columns, 
however additional columns will be disregarded.
2. (Optionally) A custom reference genome where the scaffolds names and sizes are listed in a tab-delimited file: column1=scaffold name, column2=scaffold size.

## (In)frequently Asked Questions

### How do I generate a BED file?

You might generate this file using a tool like [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables), 
shown here is a simple example of selecting a handful of genes from the human genome.

1. Go to the UCSC Table Browser. Select the human genome assembly, then "genes and predictions" and GENCODE. 
""")
st.image("media/1.png", width=600)
st.write('2. Select "paste list". Then enter your list of genes')
st.image("media/2.png", width=600)
st.write("""
3. Select the output format as tsv and plain text then click "get output".

### How do I use custom reference genomes?
         
For this tool you do not need to provide the sequence it self (e.g. a fasta file), but only a file that lists the scaffold names and sizes.
This file should be a tab-delimited file with the scaffold name in the first column and the scaffold size in the second column. There are several 
ways to generate this file, for instance using samtools on the command-line:

    samtools faidx reference.fasta
    cut -f1,2 reference.fasta.fai > reference.sizes


### I am getting warnings, what do I do?

It all depends on the warning, but here are some common ones:
         
**Overlapping regions**: This means that some of the regions in your BED file overlap with each other. By default the tools will try to merge these regions,
although it is unclear to me if this is required by the Nanopore software. It does however improve the readability of the BED file. The overlap might be the result of 
the result of a mistake, or it might be intentional. If it is intentional, you can ignore this warning.
    
**The BED file contains sex chromosomes**:. This is a warning that the BED file contains regions that are in chromosomes X, Y, Z or W. 
This is only to warn you that these regions might seemingly receive less enrichment because of their lesser copy numbers, feel free to ignore.
         
**Region is overlapping with blacklist**: The so-called blacklist is based on the exclude lists from [Delly](https://github.com/dellytools/delly)
where available. These regions are telomeres, centromeres and heterochromatin. Sequencing of these regions might be problematic, but certainly
mapping uniquely to these regions is. It is uncertain to me how much of a problem this is for adaptive sampling, but it is good to be aware of it 
especially if you are exclusively targeting these regions.

### What is an AGP file? And do I need it?

The short answer is that, no you would probably not need it. file is generated from the modified BED file, where ROIs are represented as contigs and gaps as genomic locations.
This is optional and can be used to track enrichment of the ROIs in real time during the sequencing run. To do this you have generate a new fasta reference using e.g. [agptools](https://warrenlab.github.io/agptools/):

         agptools assemble reference.fasta modifications.agp > modified_reference.fasta

""")
st.image("media/minknow_enrichment.jpg", width=800, caption="Screenshot of sequencing coverage per chromosome in MinKNOW using a modified reference genome, showing enrichment of one ROI.")
st.write("""
Also note that if your ROI is a small \% of the genome (<0.01) it might take some time before you see any enrichment, so use with caution.

### What did you do to my BED file?

If you are in doubt, you can inspect resulting BED file by importing it as a track to a genome browser like [IGV](https://igv.org/).
And if you spot anything wrong feel free to either create a [GitHub issue](https://github.com/remiolsen/adaptive_sampler_check/issues) or contact me directly at <remi-andre.olsen@scilifelab.se>.
         
""")

