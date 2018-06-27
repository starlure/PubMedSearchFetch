This package assume "biopython" is pre-installed. Otherwise it can be downloaded from https://biopython.org/wiki/Download.
biopython tutorial: http://biopython.org/DIST/docs/tutorial/Tutorial.html.

There are two ways to run the "PubMed_SearchFetch.py" code:

1. python PubMed_SearchFetch.py -i inputs.txt
For this way, all the paper PMID are saved line by line in "inputs.txt".

2. python PubMed_SearchFetch.py -s keyword
For this way, the package automatically search for the given keyword and output the info for the first 20 (default) papers.
The number of 20 can be modified in function "search_term()" as a parameter input of "Entrez.esearch()".

The outputs are saved in "outputs.csv" file.
