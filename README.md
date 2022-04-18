# nGRE Senior Project: A Sequence Searching Tool

A DNA sequence searching tool to parse for DNA sequences for nuclear hormone receptor response elements such as nGREs or other subsequences

[Access the tool's website here](https://share.streamlit.io/hewittk/ngre-senior-project/streamlit_site/streamlit_sharing.py)

## What the sequence searching tool does

The tool uses regular expression based searching to find matches to the nGRE consensus sequence or custom sequences that can contain variability in length, IUPAC codes, and/or mutations within a given gene sequence. The user is also able to input different `penalty` values for each type of mutation that are subtracted from the match's `score` that indicates how closely it resembles the sequence of interest. All scored matches are then returned in a dataframe to the user containing information about each matches sequence, location, and number of mutations.

## Why is the tool useful?

Despite the potential it brings, there is no other web browser based tool for the searching and scoring of imperfect matches to variable length sequences within gene sequences. This tool allows for the searching of nGREs and other subsequences that have variable length elements within gene sequences and returns comprehensive information about all matches to the subsequence found. The tool being browser-based, only relying on users clicking on a URL, simplifies the utilization process of the tool because there is no necesity of downloading software, launching command-line instructions, and/or using programming languages. The browser-based nature of the tool therefore simplifies a user's experience with the tool and increases the accessibility of the tool to biologists that may not have computational background and/or hardware optimized for running programming-based tools.
