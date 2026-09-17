import re

def natural_sort_key(chrom):
    match = re.search(r'\d+', chrom)
    if match:
        return (0, int(match.group()))
    if 'X' in chrom.upper():
        return (1, 0)
    elif 'Y' in chrom.upper():
        return (1, 1)
    elif 'M' in chrom.upper() or 'MT' in chrom.upper():
        return (1, 2)
    else:
        return (2, 0, chrom)

chroms = ['chr10', 'chr2', 'chr20', 'chr1', 'chr3', 'chrX', 'chrY', 'chrM', 'chr11', 'chr19']
print("Before:", chroms)
chroms.sort(key=natural_sort_key)
print("After:", chroms)
