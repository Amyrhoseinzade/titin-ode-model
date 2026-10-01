import urllib.request

# Download titin sequence from UniPort
url = "https://rest.uniprot.org/uniprotkb/Q8WZ42.fasta"
response = urllib.request.urlopen(url)
fasta = response.read().decode('utf-8')
# seq and header splitting
lines = fasta.split('\n')
seq = ''
for line in lines:
    if not line.startswith('>'):
        seq += line.strip()

groups = {
    'x1_Hyd+': ['I', 'L', 'V', 'W'],
    'x2_Hyd-': ['A', 'M', 'P'],
    'x3_Pos+': ['K', 'R', 'H'],
    'x4_Neg-': ['D', 'E'],
    'x5_Polar': ['S', 'T', 'N', 'Q'],
    'x6_Gly': ['G'],
    'x7_Cys': ['C'],
    'x8_Aro': ['Y', 'F']
}
# count & calculate Abundance
result = {}
for name, aa_list in groups.items():
    count = 0
    for aa in seq:
        if aa in aa_list:
            count += 1
    result[name] = round(count / len(seq), 4)

print(result)