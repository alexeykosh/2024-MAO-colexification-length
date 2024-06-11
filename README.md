# **Online supplement**: Short words are more likely to refer to several concepts across 192 language families

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11355636.svg)](https://doi.org/10.5281/zenodo.11355636) -->

Authors: 

- Alexey Koshevoy
- Olivier Morin
- Marie Hallo

Preregistration: [https://osf.io/8d4vh](https://osf.io/8d4vh/?view_only=40c3e016cfc64d5491e7306e7590967b)


## Reproduction 

### Downloading the code & requirements:

The code provided in this repository was executed using Python 3.11.7 and R version 4.2.2 (2022-10-31). First, clone the repository:

```bash
git clone https://github.com/alexeykosh/2024-MAO-colexification-length/
```

Then, navigate to the repository:

```bash
cd 2024-MAO-colexification-length/
```

All the required packages are listed in the `requirements.txt` file. To install the required packages, run the following command:

```bash
pip install -r requirements.txt
```

### Data:

[Lexibank](https://lexibank.clld.org/) data used in this study needs to be downloaded from [here](https://zenodo.org/records/7836668). The downloaded zip file needs to be placed in the `data/` directory. 

After downloading the data, run the following command to extract the data:

```bash
python3 preprocessing.py 
```


### Analysis:

- [analysis.R](https://github.com/alexeykosh/2024-MAO-colexification-length/blob/main/analysis.R) -- the main script for the analysis of the data. 