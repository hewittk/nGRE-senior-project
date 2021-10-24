import requests

response = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom=chr1;start=20542952;end=20618064")
print(response.status_code)
print(response.json())
