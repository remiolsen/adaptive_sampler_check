import streamlit as st
import pandas as pd
import numpy as np
from natsort import index_natsorted
import hashlib
import time
import subprocess
from io import StringIO

version = "0.1.0"
program_name = "Adaptive Sampler Check"
program_url = "https://github.com/remiolsen/adaptive_sampler_check/"
program_gh_api = (
    "https://api.github.com/repos/remiolsen/adaptive_sampler_check/git/refs/heads/main"
)
st.set_page_config(
    page_title=f"{program_name} v{version}",
    page_icon=":scissors:",
    initial_sidebar_state="expanded",
)
genomes = pd.DataFrame(
    {
        "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes.txt": "human (T2T CHM13v2.0/hs1)",
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes": "human (GRCh38/hg38)",
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes": "mouse (GRCm39/mm39)",
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes": "mouse (GRCm38/mm10)",
        "https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/rn7.chrom.sizes": "rat (mRatBN/7.2)",
        "https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.chrom.sizes": "rat (RGSC 6.0/rn6)",
        "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes": "fruit fly (BDGP6/dm6)",
        "Custom assembly": "Custom assembly",
    },
    index=["URL", "Assembly"],
).T
accepted_bed_types = ["object", "int64", "int64", "int64", "int64", "object"]
minimum_ROI_size = 0.001
maximum_ROI_size = 0.1


@st.cache_data(persist=True)
def fetch_genome(url):
    # Fetch genome size from UCSC.
    return pd.read_csv(url, sep="\t", header=None, names=["chrom", "size"])


@st.cache_data(persist=True)
def fetch_github_sha(url):
    try:
        return (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
            .decode("ascii")
            .strip()
        )
    except subprocess.CalledProcessError:
        return None
    except subprocess.SubprocessError:
        return None


def find_overlaps(bed_df):
    # Find overlaps in BED file
    # returns a dataframe with the following columns: chrom, start1, end2, overlap, bed1, bed2
    overlaps = []
    for i, row in bed_df.iterrows():
        for j, row2 in bed_df.iterrows():
            if row[0] == row2[0] and i < j:
                ovl = min(row[2], row2[2]) - max(row[1], row2[1])
                if ovl > 0:
                    overlaps.append([row[0], min(row[1],row2[1]), max(row[2],row2[2]), ovl, i, j])

    return pd.DataFrame(overlaps, columns=["chrom", "start1", "end2", "overlap", "bed1", "bed2"])


# Session initialization
def default_state():
    st.session_state["state"] = 0
    st.session_state["bed_df"] = pd.DataFrame()
    st.session_state["bed_columns"] = False
    st.session_state["bed_types"] = False
    st.session_state["bed_sex"] = False
    st.session_state["bed_order"] = False
    st.session_state["bad_input_bed"] = False
    st.session_state["bed_overlaps"] = pd.DataFrame()
    st.session_state["assembly_df"] = pd.DataFrame()
    st.session_state["chrom_error"] = False
    st.session_state["size_error"] = False
    st.session_state["QC_failed"] = False
    st.session_state["BED_mod_failed"] = False
    st.session_state["mod_bed"] = pd.DataFrame()
    st.session_state["metadata"] = pd.DataFrame()
    st.session_state["out_bed"] = ""
    st.session_state["time_hash"] = ""


if "state" not in st.session_state:
    default_state()


"""
# Adaptive Sampler Check
---
This performs a series of checks on the input BED file and selected assembly for the Oxford Nanopore adaptive sampling method. 
Additionally, it it provides the ability to adjust the region of interest (ROI) size to be used to be used by adding a buffer to each site,
whilst keeping keeping the total and individual ROI sizes within a recommended range.
"""

if st.button(
    ":rewind: Start over!",
    key="restart_button",
    help="This will reset the app to its initial state.",
):
    default_state()
    st.rerun()

"""
## 1. :page_facing_up: Upload a BED file
"""

uploaded_file = st.file_uploader(
    "Upload a BED file",
    key="bed_uploader",
    disabled=st.session_state["state"] > 0,
    help="""
                                 A minimum of 3 columns are required: chromosome, start, end. 
                                 Specification of the BED format can be found in the UCSC genome browser FAQ, [here](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
                                 """,
)

# To read file as bytes:
if uploaded_file is not None:
    failed = False
    bytes_data = uploaded_file.getvalue()
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    string_data = []
    for line in stringio.readlines():
        if not line.startswith("track"):
            string_data.append(line)
    st.session_state["bed_df"] = pd.read_csv(
        StringIO("\n".join(string_data)), header=None, sep="\t"
    )
    bed_exp = st.expander("See content of the file")
    bed_exp.write(st.session_state["bed_df"])
    num_cols = st.session_state["bed_df"].shape[1]
    st.session_state["bed_overlaps"] = find_overlaps(st.session_state["bed_df"]).empty

    for i, row in st.session_state["bed_df"].iterrows():
        if row[0].lower() in ["x", "y", "z", "w", "chrx", "chry", "chrz", "chrw"]:
            st.session_state["bed_sex"] = True
        if row[2] < row[1]:
            st.session_state["bed_order"] = True
            st.session_state["bad_input_bed"] = True
    if num_cols < 3:
        st.session_state["bed_columns"] = True
        st.session_state["bad_input_bed"] = True
    res = st.session_state["bed_df"].dtypes
    if list(res[:num_cols]) != accepted_bed_types[:num_cols]:
        st.session_state["bed_types"] = res[:6].any()
        st.session_state["bad_input_bed"] = True
    if st.session_state["state"] < 1:
        st.session_state["state"] = 1

"""
## 2. :dna: Choose an assembly
"""

assembly = st.selectbox(
    "Choose an assembly",
    genomes,
    index=1,
    key="assembly_select",
    disabled=st.session_state["state"] > 1,
    help="This will fetch the size information for the chromosomes in your BED file.",
)
assembly_df = pd.DataFrame()
assembly_exp = st.expander("See content of the assembly")

custom_assembly = st.file_uploader(
    "Or upload a custom assembly .sizes file",
    key="custom_assembly",
    disabled=st.session_state["state"] > 1 or assembly != "Custom assembly",
    help="This will fetch the size information for the chromosomes in your BED file.",
)

if st.button("Fetch", disabled=st.session_state["state"] > 2, key="fetch_button"):
    if custom_assembly is not None:
        assembly_df = pd.read_csv(custom_assembly, sep="\t", header=None, names=["chrom", "size"])
        st.session_state["assembly_df"] = assembly_df
        st.session_state["state"] = 2
    else:
        assembly_df = fetch_genome(genomes.index[genomes["Assembly"] == assembly][0])
        st.session_state["assembly_df"] = assembly_df
        st.session_state["state"] = 2
if st.session_state["state"] < 3 and not st.session_state["assembly_df"].empty:
    assembly_exp.write(st.session_state["assembly_df"])
    for chrom in st.session_state["bed_df"][0].unique():
        if chrom not in st.session_state["assembly_df"]["chrom"].values:
            st.session_state["chrom_error"] = chrom
    for i, row in st.session_state["bed_df"].iterrows():
        if row[0] not in st.session_state["assembly_df"]["chrom"].values:
            pass
        elif (
            row[2]
            > st.session_state["assembly_df"]
            .loc[st.session_state["assembly_df"]["chrom"] == row[0], "size"]
            .values[0]
        ):
            st.session_state["size_error"] = row[1]
            st.session_state["bad_input_bed"] = True

"""
## 3. :sleuth_or_spy: QC results
"""

if st.session_state["state"] > 1:
    """ :page_facing_up: BED file checks """
    bed_results = st.container(border=True)
    if st.session_state["bed_columns"]:
        bed_results.write(
            f":x: The file should have at least 3 columns, but found only `{st.session_state['bed_columns']}`"
        )
    else:
        bed_results.write(":white_check_mark: The file has at least 3 columns.")
    if st.session_state["bed_types"]:
        bed_results.write(":x: The bed values do not match the expected types:")
        bed_results.write(accepted_bed_types)
    else:
        bed_results.write(":white_check_mark: The file has the correct columns types.")

    if not st.session_state["bed_overlaps"]:
        bed_results.write(":waning: The BED file contains overlapping regions")
    if st.session_state["bed_sex"]:
        bed_results.write(":warning: The BED contains sex chromosomes.")
    if st.session_state["bed_order"]:
        bed_results.write(
            ":x: The BED file contains entries where the start coordinate is larger than the end coordinate."
        )

    """ :dna: Assembly checks """
    assembly_results = st.container(border=True)
    if st.session_state["chrom_error"]:
        assembly_results.write(
            f":x: The file contains the following chromosomes not present in the assembly: `{st.session_state['chrom_error']}`"
        )
    else:
        assembly_results.write(
            ":white_check_mark: All chromosomes in the file are present in the assembly."
        )
    if st.session_state["size_error"]:
        assembly_results.write(
            f":x: The file contains the following coordinates larger than the chromosome size: `{st.session_state['size_error']}`"
        )
    else:
        assembly_results.write(
            ":white_check_mark: All coordinates in the file are smaller than the chromosome size."
        )

"""
## 4. :gear: Parameters
"""

col1, col2 = st.columns(2)

with col1:
    experiment_name = st.text_input(
        "Experiment name",
        key="experiment_name",
        max_chars=100,
        placeholder="adaptive_sampling_01",
    )
    minimum_buffer_size = st.number_input(
        "Minimum buffer size (bp)",
        min_value=10000,
        max_value=10000000,
        key="min_buffer_size",
        value=40000,
        step=10000,
        format="%d",
        help="""
                                            In cases where expansion of the ROI is restricted (e.g, beginning and ends of chromosomes) the program requires a minimum buffer to expand into (default: 40000bp).
                                            If this buffer is not available the sequencing run might not be able to cover the entire ROI.
                                            """,
    )
with col2:
    operator_name = st.text_input(
        "Operator name",
        key="operator_name",
        max_chars=100,
        placeholder="N. Hvermannsen",
    )
    roi_size = st.slider(
        "ROI size",
        min_value=0.0,
        max_value=maximum_ROI_size,
        key="ROI_slider",
        value=minimum_ROI_size,
        step=0.001,
        format="%.3f",
        help="This will expand the size of the ROIs of the input bed file to a selected minimum size (default: 0.1pct of the genome size)",
    )

scol1, scol2, scol3 = st.columns(3)
with scol1:
    size_override = st.toggle(
        "Override size limitations",
        key="size_override",
        value=False,
        help="This will override the allowable recommended minimum and maximum ROI sizes. Use with caution.",
    )
with scol2:
    chrom_prune = st.toggle(
        "Prune missing chroms",
        key="chrom_prune",
        value=False,
        help="This will remove any rows in the BED file that contain chromosomes not present in the selected assembly.",
    )
with scol3:
    no_sort = st.toggle(
        "No sorting",
        key="no_sort",
        value=False,
        help="This will try to keep the original order of the BED file instead of sorting it by chromosome and start position.",
    )
merge_overlaps = st.toggle(
    "Merge overlaps",
    key="merge_overlaps",
    value=False,
    help="This will merge overlapping regions in the BED file into a single region.",
)


@st.cache_data
def modify_bed(bed_df, assembly_df, min_size, minimum_buffer_size):
    mod_bed = bed_df.copy()
    messages = []
    bad_bed = False
    for i, row in mod_bed.iterrows():
        r_size = row[2] - row[1]
        length_cutoff = assembly_df.loc[assembly_df["chrom"] == row[0], "size"].values[
            0
        ]
        missing = np.round((min_size - r_size) / 2, 0)
        # The new ROI fits within the chromosome
        if row[1] - missing > 0 and row[2] + missing < length_cutoff:
            mod_bed.loc[i, 1] = row[1] - missing
            mod_bed.loc[i, 2] = row[2] + missing
        # The new ROI limited by the chromosome start. Using minimum buffer, adding buffer downstream
        elif (
            row[1] - minimum_buffer_size > 0
            and row[2] + missing * 2 - minimum_buffer_size < length_cutoff
        ):
            mod_bed.loc[i, 1] = 0
            mod_bed.loc[i, 2] = row[2] + missing * 2 - row[1]
        # The new ROI limited by the chromosome end. Using minimum buffer, adding buffer upstream
        elif (
            row[1] - missing * 2 + minimum_buffer_size > 0
            and row[2] + minimum_buffer_size < length_cutoff
        ):
            mod_bed.loc[i, 1] = row[1] - missing * 2 - (length_cutoff - row[2])
            mod_bed.loc[i, 2] = length_cutoff
        else:
            messages.append(
                f":x: Bed modification failed - new interval illegal `{row[0]} : {row[1]-missing} - {row[2]+missing}`"
            )
            bad_bed = True


    return bad_bed, mod_bed, messages


if st.button("Generate", disabled=st.session_state["state"] < 2 or st.session_state["bad_input_bed"] or (st.session_state["chrom_error"] and not chrom_prune), key="generate_button"):
    bad_bed = False
    genome_size = st.session_state["assembly_df"]["size"].sum()
    roi_size = st.session_state["ROI_slider"]
    input_bed = st.session_state["bed_df"]
    if st.session_state["chrom_prune"]:
        input_bed = input_bed[
            input_bed[0].isin(st.session_state["assembly_df"]["chrom"])
        ]
    bad_bed, mod_bed, messages = modify_bed(
        input_bed,
        st.session_state["assembly_df"],
        np.round(roi_size * genome_size, 0),
        minimum_buffer_size,
    )
    if not st.session_state["no_sort"]:
        mod_bed = mod_bed.sort_values(by=1)
        mod_bed = mod_bed.sort_values(
            by=0, key=lambda x: np.argsort(index_natsorted(mod_bed[0]))
        )
    total_fraction = (mod_bed[2].sum() - mod_bed[1].sum()) / genome_size
    if total_fraction > maximum_ROI_size and not st.session_state["size_override"]:
        messages.append(
            f":x: Total ROI size is larger than the recommended maximum of `{maximum_ROI_size}`"
        )
        bad_bed = True
    # Find overlaps in the modified bed
    mod_overlaps = find_overlaps(mod_bed)
    if not mod_overlaps.empty and not merge_overlaps:
        for i, row in mod_overlaps.iterrows():
            messages.append(
                f":x: Found {row[3]} bp overlap in chromosome {row[0]}, between coordinates {row[1]} and {row[2]}"
            )
        bad_bed = True
    elif not mod_overlaps.empty and merge_overlaps:
        messages.append(
            f":warning: Found {mod_overlaps.shape[0]} overlapping regions, merging them into a single region"
        )
        for i, row in mod_overlaps.iterrows():
            # Remove last entry of overlap and replace the first with merged region
            mod_bed.iloc[row[4]] = [row[0], row[1], row[2]]
            mod_bed = mod_bed.drop(row[5])
    if len(messages) == 0 and not bad_bed:
        messages.append(":white_check_mark: Bed modification successful")

    total_size = mod_bed[2].sum() - mod_bed[1].sum() / 1000000
    median_size = (mod_bed[2] - mod_bed[1]).median() / 1000000
    messages.append(
        f":information_source: total size: `{int(np.round(total_size))} mbp` (`{np.round(total_fraction * 100,2)} %` of genome), median size: `{np.round(median_size,2)} mbp`"
    )
    msg_container = st.container(border=True)
    for msg in messages:
        msg_container.write(msg)
    mod_exp = st.expander("See modified BED")
    mod_exp.write(mod_bed)

    st.session_state["mod_bed"] = mod_bed
    if not st.session_state["operator_name"]:
        st.error(":exclamation: Operator name is required")
        bad_bed = True
    if not st.session_state["experiment_name"]:
        st.error(":exclamation: Experiment name is required")
        bad_bed = True

    if not bad_bed:
        st.session_state["state"] = 4
    else:
        st.session_state["state"] = 3


"""
## 5. :arrow_down: Download 
"""

# Fallback names for output files
mod_bed_name = f"{''.join(i for i in experiment_name if i.isalnum())}.bed"
metadata_name = f"{st.session_state['experiment_name']}_metadata.csv"

gh_hash = fetch_github_sha(program_gh_api)

if st.session_state["state"] > 3:
    local_time = time.localtime()
    tzname_local = local_time.tm_zone
    ts = pd.Timestamp.now(tz=tzname_local).isoformat(timespec="seconds")
    if st.session_state["time_hash"] == "":
        st.session_state["time_hash"] = hashlib.md5(ts.encode("utf-8")).hexdigest()
    mod_bed_name = f"{''.join(i for i in experiment_name if i.isalnum())}_{st.session_state['time_hash'][:10]}.bed"
    st.session_state["out_bed"] = st.session_state["mod_bed"].to_csv(
        index=False, sep="\t", header=False
    )
    st.session_state["metadata"] = pd.DataFrame(
        {
            "key": [
                "operator",
                "experiment",
                "assembly",
                "assembly_url",
                "custom_assembly",
                "custom_assembly_md5",
                "original_bed",
                "original_bed_md5",
                "adaptive_sampling_bed",
                "adaptive_sampling_bed_md5",
                "date_generated",
                "date_generated_hash",
                "version",
                "program_url",
                "program_latest_commit",
            ],
            "value": [
                st.session_state["operator_name"],
                st.session_state["experiment_name"],
                st.session_state["assembly_select"],
                genomes.index[
                    genomes["Assembly"] == st.session_state["assembly_select"]
                ][0],
                custom_assembly.name if custom_assembly else "NA",
                hashlib.md5(custom_assembly.getvalue()).hexdigest() if custom_assembly else "NA",
                uploaded_file.name,
                hashlib.md5(bytes_data).hexdigest(),
                mod_bed_name,
                hashlib.md5(st.session_state["out_bed"].encode("utf-8")).hexdigest(),
                ts,
                st.session_state["time_hash"],
                version,
                program_url,
                gh_hash,
            ],
            "description": [
                "Name of the operator",
                "Name of the experiment",
                "Name of the assembly used. Check that this corresponds to the assembly found on the sequencer",
                "URL of the assembly chromosome sizes",
                "Name of the custom assembly used",
                "MD5 hash of the custom assembly file",
                "Name of the original BED file",
                "MD5 hash of the original BED file",
                "Name of the modified BED file, i.e. the one generated by adaptive sampler check program and to be used on the sequencer",
                "MD5 hash of the modified BED file",
                "Date and time of generation",
                "MD5 hash of the date and time of generation",
                "Version of the program",
                "URL of the program",
                "Latest commit hash",
            ],
        }
    )
    st.table(st.session_state["metadata"])
    metadata_name = f"{st.session_state['experiment_name']}_metadata_{st.session_state['time_hash'][:10]}.csv"

bcol1, bcol2 = st.columns(2)
with bcol1:
    st.download_button(
        "Download BED",
        st.session_state["out_bed"],
        mod_bed_name,
        "text/plain",
        key="dl_bed",
        disabled=st.session_state["state"] < 4,
    )
with bcol2:
    st.download_button(
        "Download Metadata",
        st.session_state["metadata"].to_csv(index=False),
        metadata_name,
        "text/plain",
        key="dl_metadata",
        disabled=st.session_state["state"] < 4,
    )


"""
---
"""
st.markdown(
    f"Created using **{program_name}** *v{version} ({gh_hash})* - [GitHub]({program_url})"
)
