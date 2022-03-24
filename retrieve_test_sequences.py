import requests

# retrieve TSLP sequence from 30000 before start site to 5000 after end site
TSLP_request = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom=chr18;start=32785383;end=32824799")
print("Status Code: ", TSLP_request.status_code)
TSLP_response = TSLP_request.json()
TSLP_sequence = TSLP_response['dna']
with open("test_gene_sequences/TSLP.txt", "w") as file:
    file.write(TSLP_sequence)
    file.close()
