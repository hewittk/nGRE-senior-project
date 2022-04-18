# nGRE Senior Project: A Sequence Searching Tool

A DNA sequence searching tool to parse for DNA sequences for nuclear hormone receptor response elements such as nGREs or other subsequences

[Access the tool's website here](https://share.streamlit.io/hewittk/ngre-senior-project/streamlit_site/streamlit_sharing.py)

## What the sequence searching tool does

The tool uses regular expression based searching to find matches to the nGRE consensus sequence or custom sequences that can contain variability in length, IUPAC codes, and/or mutations within a given gene sequence. The user is also able to input different `penalty` values for each type of mutation that are subtracted from the match's `score` that indicates how closely it resembles the sequence of interest. All scored matches are then returned in a dataframe to the user containing information about each matches sequence, location, and number of mutations.
