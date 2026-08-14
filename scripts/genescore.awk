# Per-cell marker-gene fragment counts from a cellranger-atac fragments.tsv stream.
#
# Restricted to the barcodes the authors called as cells (BCFILE, taken from the
# deposited cell_heteroplasmic_df), so the cell set here is exactly the cell set
# that carries mtDNA data. Windows (BEDFILE) are gene body + 2 kb promoter.
#
# Interval lookup is bucketed at 100 kb so the common case (a fragment nowhere
# near a marker) costs one hash probe rather than a scan over all 59 windows.
BEGIN {
    FS = "\t"
    while ((getline line < BCFILE) > 0) keep[line] = 1
    n = 0
    while ((getline line < BEDFILE) > 0) {
        split(line, a, "\t")
        n++
        chrom[n] = a[1]; s[n] = a[2] + 0; e[n] = a[3] + 0; gname[n] = a[4]
        b0 = int(s[n] / 100000); b1 = int(e[n] / 100000)
        for (b = b0; b <= b1; b++) {
            k = a[1] ":" b
            bin[k] = bin[k] " " n
        }
    }
}
/^#/ { next }
{
    bc = $4
    if (!(bc in keep)) next
    tot[bc]++
    fs = $2 + 0; fe = $3 + 0
    k1 = $1 ":" int(fs / 100000)
    k2 = $1 ":" int(fe / 100000)
    hit = ""
    if (k1 in bin) hit = bin[k1]
    if (k2 != k1 && (k2 in bin)) hit = hit " " bin[k2]
    if (hit == "") next
    m = split(hit, idx, " ")
    delete seen
    for (j = 1; j <= m; j++) {
        i = idx[j] + 0
        if (i == 0 || (i in seen)) continue
        seen[i] = 1
        if (fs < e[i] && fe > s[i]) mark[gname[i] "|" bc]++
    }
}
END {
    for (x in tot) print "T\t" x "\t" tot[x]
    for (x in mark) {
        split(x, p, "|")
        print "M\t" p[2] "\t" p[1] "\t" mark[x]
    }
}
