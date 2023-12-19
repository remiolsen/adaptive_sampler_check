import streamlit as st
import pandas as pd
import numpy as np
from io import StringIO

version = "0.1.0"
genomes = pd.DataFrame({"https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes.txt": "human (T2T CHM13v2.0/hs1)", 
                        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes": "human (GRCh38/hg38)", 
                        "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes": "mouse (GRCm39/mm39)", 
                        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes": "mouse (GRCm38/mm10)"
                        },
                        index=["URL", "Assembly"]).T
accepted_bed_types = ['object', 'int64', 'int64', 'int64', 'int64', 'object']
minimum_ROI_size = 0.001


@st.cache_data
def fetch_genome(url):
    # Fetch genome size from UCSC.
    return pd.read_csv(url, sep='\t', header=None, names=["chrom", "size"])

def write_bed_checks():
    bed_results = st.container(border=True)
    if st.session_state['bed_columns']:
        bed_results.write(f":x: The file should have at least 6 columns, but found only {st.session_state['bed_columns']}")
    else:
        bed_results.write(":white_check_mark: The file has at least 6 columns.")
    if st.session_state['bed_types']:
        bed_results.write(":x: The bed values do not match the expected types:")
        bed_results.write(accepted_bed_types)
    else:
        bed_results.write(":white_check_mark: The file has the correct columns types.")

def write_assembly_checks():
    assembly_results = st.container(border=True)
    if st.session_state['chrom_error']:
        assembly_results.write(f":x: The file contains the following chromosomes not present in the assembly: {st.session_state['chrom_error']}")
    else:
        assembly_results.write(":white_check_mark: All chromosomes in the file are present in the assembly.")
    if st.session_state['size_error']:
        assembly_results.write(f":x: The file contains the following coordinates larger than the chromosome size: {st.session_state['size_error']}")
    else:
        assembly_results.write(":white_check_mark: All coordinates in the file are smaller than the chromosome size.")

# Session flow
if 'state' not in st.session_state:
    st.session_state['state'] = 0
    st.session_state['bed_columns'] = False
    st.session_state['bed_types'] = False
    st.session_state['chrom_error'] = False
    st.session_state['size_error'] = False

'''
# Adaptive Sampler Check
---
This is a sanity check for the adaptive sampling input bed file.

## 1. :page_facing_up: Upload a BED file
'''

uploaded_file = st.file_uploader("Upload a BED file", key="bed_uploader", disabled=st.session_state['state'] > 0, label_visibility='hidden')

# To read file as bytes:
if uploaded_file is not None:
    failed = False
    bytes_data = uploaded_file.getvalue()
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    string_data = stringio.read()
    bed_df = pd.read_csv(uploaded_file, header=None, sep='\t')
    bed_exp = st.expander("See content of the file")
    bed_exp.write(bed_df)
    if bed_df.shape[1] < 6:
        st.session_state['bed_columns'] = bed_df.shape[1] 
    res = bed_df.dtypes
    if list(res[:6]) != accepted_bed_types:
        st.session_state['bed_types'] = res[:6].any()
    st.session_state['state'] = 1

'''
## 2. :dna: Choose an assembly
'''

assembly = st.selectbox("Choose and assembly", genomes, index=1, key="assembly_select", label_visibility='hidden')
fetch_button = st.button('Fetch', disabled=st.session_state['state'] > 1, key="fetch_button")
assembly_df = pd.DataFrame()
if st.session_state['state'] == 1 and fetch_button:
    assembly_df = fetch_genome(genomes.index[genomes["Assembly"] == assembly][0])
    assembly_exp = st.expander("See content of the assembly")
    assembly_exp.write(assembly_df)

    for chrom in bed_df[0].unique():
        if chrom not in assembly_df["chrom"].values:
            st.session_state['chrom_error'] = chrom
    for i, row in bed_df.iterrows():
        if row[1] > assembly_df.loc[assembly_df["chrom"] == row[0], "size"].values[0]:
            st.session_state['size_error'] = row[1]

    st.session_state['state'] = 2

'''
## 3. :sleuth_or_spy: QC results
'''

if st.session_state['state'] > 1:
    ''' :page_facing_up: BED file checks '''
    write_bed_checks()
    ''' :dna: Assembly checks '''
    write_assembly_checks()
    st.session_state['state'] = 3

'''
## 4. :gear: Parameters
'''

'''
## 5. :arrow_down: Download 
'''