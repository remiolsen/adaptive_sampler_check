# Adaptive Sampler Check

This performs a series of checks on the input BED file and selected assembly for the Oxford Nanopore adaptive sampling method. 
Additionally, it it provides the ability to adjust the region of interest (ROI) size to be used to be used by adding a buffer to each site,
whilst keeping keeping the total and individual ROI sizes within a recommended range.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)

## Installation

This app is deployed in the Streamlit cloud click >[here](https://adaptive-sampler-check.streamlit.app/)< to run it now!
Or take a look at the available deployment options available in the [documentation](https://docs.streamlit.io/deploy).

This software requires `git` `python` v.3.10 and `pip` to install locally:

```
git clone https://github.com/remiolsen/adaptive_sampler_check/
cd adaptive_sampler_check
pip install -r requirements.txt
```

## Usage

To run this app locally use:

```
streamlit run 01-App.py
```

