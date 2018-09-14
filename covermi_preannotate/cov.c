
struct readstruct {
    int chrom;
    int start;
    int stop;
    int strand;
    } Read;




        for chrom, start, stop, strand in data:
            depth = 1
            cov = self.data[chrom]
            upperx = -1
            lowerx = -1
            for x in xrange(len(cov) - 1, -1, -1):
                if stop < cov[x][0]:
                    continue
                if stop >= cov[x][1] and start <= cov[x][0]:
                    cov[x][2] += depth
                    if start == cov[x][0]:
                        break
                else:
                    if stop < cov[x][1]:
                        upperx = x
                    if start > cov[x][0]:
                        lowerx = x
                        break
            if upperx != -1:
                cov.insert(upperx, [cov[upperx][0], stop, cov[upperx][2] + depth])
                cov[upperx+1][0] = stop + 1
            if lowerx != -1:
                cov.insert(lowerx+1, [start, cov[lowerx][1], cov[lowerx][2] + depth * (upperx!=lowerx)])
                cov[lowerx][1] = start - 1
                if upperx == lowerx:
                    cov[lowerx][2] -= depth

            if amplicons:
                try:
                    fr_depth[(start, PLUS)].f_depth += 1
                    self.ontarget += 1
                except KeyError:
                    try:
                        fr_depth[(stop, MINUS)].r_depth += 1
                        self.ontarget += 1
                    except KeyError:
                        self.offtarget += 1                        

        self.amplicon_info.sort()
        if bam:
            self.unmapped = bam.unmapped
            self.mapped = bam.mapped
            self.reads = self.mapped + self.unmapped
            self.percent_unmapped = self.unmapped *100.0 / (self.reads or 1)
        if amplicons:
            self.percent_offtarget = self.offtarget * 100.0 / (self.ontarget + self.offtarget or 1)

