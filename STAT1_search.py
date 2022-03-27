import regex

stat1_file = open("test_gene_sequences/STAT1.txt")
stat1_sequence = stat1_file.read()
nGRE_start = regex.compile("CACC")

for match_location in nGRE_start.finditer(stat1_sequence):
    print([match_location.start(), match_location.end()])

print(stat1_sequence[27177:27190])
