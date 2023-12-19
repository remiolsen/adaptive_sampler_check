import streamlit as st
import pandas as pd
import numpy as np
from io import StringIO

version = "0.1.0-dev"
genomes = pd.DataFrame({"https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes.txt": "human (T2T CHM13v2.0/hs1)", 
                        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes": "human (GRCh38/hg38)", 
                        "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes": "mouse (GRCm39/mm39)", 
                        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes": "mouse (GRCm38/mm10)",
                        "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes": "fruit fly (BDGP6/dm6)",
                        },
                        index=["URL", "Assembly"]).T
accepted_bed_types = ['object', 'int64', 'int64', 'int64', 'int64', 'object']
minimum_ROI_size = 0.001
maximum_ROI_size = 0.1


@st.cache_data
def fetch_genome(url):
    # Fetch genome size from UCSC.
    return pd.read_csv(url, sep='\t', header=None, names=["chrom", "size"])

def write_bed_checks():
    bed_results = st.container(border=True)
    if st.session_state['bed_columns']:
        bed_results.write(f":x: The file should have at least 3 columns, but found only `{st.session_state['bed_columns']}`")
    else:
        bed_results.write(":white_check_mark: The file has at least 3 columns.")
    if st.session_state['bed_types']:
        bed_results.write(":x: The bed values do not match the expected types:")
        bed_results.write(accepted_bed_types)
    else:
        bed_results.write(":white_check_mark: The file has the correct columns types.")
    # To be implemented
    bed_results.write(":grey_question: The BED file entries are mostly evenly sized.")
    bed_results.write(":grey_question: The BED does not contain any overlapping regions.")

def write_assembly_checks():
    assembly_results = st.container(border=True)
    if st.session_state['chrom_error']:
        assembly_results.write(f":x: The file contains the following chromosomes not present in the assembly: `{st.session_state['chrom_error']}`")
    else:
        assembly_results.write(":white_check_mark: All chromosomes in the file are present in the assembly.")
    if st.session_state['size_error']:
        assembly_results.write(f":x: The file contains the following coordinates larger than the chromosome size: `{st.session_state['size_error']}`")
    else:
        assembly_results.write(":white_check_mark: All coordinates in the file are smaller than the chromosome size.")
    # To be implemented
    assembly_results.write(":grey_question: The BED does not contain any sex chromosomes.")

def modify_bed(in_bed, assembly, min_size, min_buffer):
    bed = in_bed.copy()
    messages = []
    for i, row in bed.iterrows():
        r_size = row[2] - row[1]
        length_cutoff = assembly.loc[assembly["chrom"] == row[0], "size"].values[0]
        missing = np.round((min_size - r_size) / 2, 0)
        if row[1] - missing > 0 or row[2] + missing < length_cutoff:
            bed.loc[i, 1] = row[1] - missing
            bed.loc[i, 2] = row[2] + missing
        elif row[1] - min_buffer > 0 and row[2] + missing*2 - min_buffer < length_cutoff:
            bed.loc[i, 1] = 0
            bed.loc[i, 2] = row[2] + missing*2 - row[1]
        elif row[1] - missing*2 + min_buffer > 0 and row[2] + min_buffer < length_cutoff:
            bed.loc[i, 1] = row[1] - missing*2 - (length_cutoff - row[2])
            bed.loc[i, 2] = length_cutoff
        else:
            messages.append(f":x: Bed modification failed - new interval illegal `{row[0]} : {row[1]-missing} - {row[2]+missing}`")

    if len(messages) == 0:
        messages.append(":white_check_mark: Bed modification successful")
    return bed, messages            



# Session flow
if 'state' not in st.session_state:
    st.session_state['state'] = 0
    st.session_state['bed_columns'] = False
    st.session_state['bed_types'] = False
    st.session_state['assembly_df'] = pd.DataFrame()
    st.session_state['chrom_error'] = False
    st.session_state['size_error'] = False
    st.session_state['QC_failed'] = False
    st.session_state['BED_mod_failed'] = False


'''
# Adaptive Sampler Check
---
This performs a series of checks on the input BED file and selected assembly for the Oxford Nanopore adaptive sampling method. 
Additionally, it it provides the ability to adjust the region of interest (ROI) size to be used to be used by adding a buffer to each site,
whilst keeping keeping the total and individual ROI sizes within a recommended range.

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
    num_cols = bed_df.shape[1]
    if num_cols < 3:
        st.session_state['bed_columns'] = bed_df.shape[1] 
    res = bed_df.dtypes
    if list(res[:num_cols]) != accepted_bed_types[:num_cols]:
        st.session_state['bed_types'] = res[:6].any()
    st.session_state['state'] = 1

'''
## 2. :dna: Choose an assembly
'''

assembly = st.selectbox("Choose and assembly", genomes, index=1, key="assembly_select", label_visibility='hidden')
fetch_button = st.button('Fetch', disabled=st.session_state['state'] > 1, key="fetch_button")
assembly_df = pd.DataFrame()
assembly_exp = st.expander("See content of the assembly")
if st.session_state['state'] == 1 and fetch_button:
    assembly_df = fetch_genome(genomes.index[genomes["Assembly"] == assembly][0])
    assembly_exp.write(assembly_df)

    for chrom in bed_df[0].unique():
        if chrom not in assembly_df["chrom"].values:
            st.session_state['chrom_error'] = chrom
    for i, row in bed_df.iterrows():
        if not st.session_state['chrom_error'] and row[2] > assembly_df.loc[assembly_df["chrom"] == row[0], "size"].values[0]:
            st.session_state['size_error'] = row[1]
        elif st.session_state['chrom_error']:
            st.session_state['size_error'] = "Not found"
    st.session_state['assembly_df'] = assembly_df
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
    
experiment_name = st.text_input("Experiment name", key="experiment_name", label_visibility='hidden', max_chars=100, placeholder="adaptive_sampling_01")
''' Experiment name '''
operator_name = st.text_input("Operator name", key="operator_name", label_visibility='hidden', max_chars=100, placeholder="N. Hvermannsen")
''' Operator name '''
minimum_buffer_size = st.number_input("Minimum buffer size (bp)", min_value=10000, max_value=10000000, key="min_buffer_size", value=40000, step=10000, format="%d", label_visibility='hidden')
''' Minimum buffer size (bp) '''
roi_size = st.slider("ROI size", min_value=minimum_ROI_size, max_value=maximum_ROI_size, key="ROI_slider",value=0.001, step=0.001, format="%.3f", label_visibility='hidden')
''' ROI size '''
size_override = st.toggle("Override size limitations", key="size_override", value=False, label_visibility='hidden', disabled=True)
''' Override size limitations '''
chrom_prune = st.toggle("Prune missing chromosomes", key="chrom_prune", value=False, label_visibility='hidden', disabled=True)
''' Prune missing chromosomes '''

if st.button('Generate', disabled=st.session_state['state'] < 2, key="generate_button"):
    genome_size = st.session_state['assembly_df']["size"].sum()
    modded_bed, messages = modify_bed(bed_df, st.session_state['assembly_df'], np.round(roi_size * genome_size, 0), minimum_buffer_size)
    msg_container = st.container(border=True)
    for msg in messages:
        msg_container.write(msg)
    mod_exp = st.expander("See modified BED")
    mod_exp.write(modded_bed)


'''
## 5. :arrow_down: Download 
'''

