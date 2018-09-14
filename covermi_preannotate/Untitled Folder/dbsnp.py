import pdb, csv, io, time



def firstwhere(f, lower, upper, func):
        firstrowstarts = lower
        f.seek(lower)
        row = f.readline()
        if func(row):
            return lower
        firstrowends = f.tell()
        lastrowstarts = upper
        while firstrowends != lastrowstarts:
            mid = (lower+upper)//2
            f.seek(mid)
            f.readline()
            testrowstarts = f.tell()
            row = f.readline()
            testrowends = f.tell()
            if func(row):
                upper = mid
                lastrowstarts = testrowstarts
            elif testrowends == upper:
                return upper + 1
            else:
                lower = mid
                firstrowstarts = testrowstarts
                firstrowends = testrowends
        return lastrowstarts

t = time.time()
with open("", "rb") as f:


with open("/home/ed/Desktop/macbook_air_backup_270316/Desktop/snps/All_snps.vcf", "rU") as f:
    row = "#"
    while row.startswith("#"):
        lower = f.tell()
        row = f.readline()
    f.seek(0, 2)
    vcfend = f.tell()

    index = {}
    while lower < vcfend:
        f.seek(lower)
        row = f.readline()
        chrom = row[0:row.find("\t")]
        upper = firstwhere(f, lower, vcfend, lambda row: row[0:row.find("\t")]!=chrom)
        index[chrom] = (lower, upper - 1)
        lower = upper
    
    for chrom, start, stop in [("1", 158647421, 158647449), ("1", 158647661, 158647684)]: 
        offset = firstwhere(f, index[chrom][0], index[chrom][1], lambda row: int(row.split("\t")[1]) >= start)
        f.seek(offset)
        while f.tell() < index[chrom][1]:
            row = f.readline()
            pos = int(row.split("\t")[1])
            if pos > stop:
                break
            print row

