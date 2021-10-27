import requests

# get gene name and chromosome, txStart, txEnd for each gene


# request sequence for gene
response = requests.get("http://api.genome.ucsc.edu/getData/sequence?genome=mm10;chrom=chr1;start=20542952;end=20618064")
print(response.status_code)
print(response.json())

# mine json for response

# save to file named after gene


def main():
    print("In main")

if __name__ == "__main__":
    main()
