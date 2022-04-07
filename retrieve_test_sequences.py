import requests

# retrieve TSLP sequence from 30000 before start site to 5000 after end site
TSLP_request = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom=chr18;start=32785383;end=32824799")
print("Status Code: ", TSLP_request.status_code)
TSLP_response = TSLP_request.json()
TSLP_sequence = TSLP_response['dna']
with open("test_gene_sequences/TSLP.txt", "w") as file:
    file.write(TSLP_sequence)
    file.close()

# retrieve STAT1 sequence from 30000 before start site to 5000 after end site
STAT1_request = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm9;chrom=chr1;start=52034034;end=52110508")
print("Status Code: ", STAT1_request.status_code)
STAT1_response = STAT1_request.json()
STAT1_sequence = STAT1_response['dna']
print(STAT1_sequence)
with open("test_gene_sequences/STAT1.txt", "w") as file:
    file.write(STAT1_sequence)
    file.close()
