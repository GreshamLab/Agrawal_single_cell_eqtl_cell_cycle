#!/usr/bin/env python3
"""
Pairwise nucleotide-diversity matrices from the 96-strain all-sites GATK VCF
(no MAF filter), so that pi can be computed for any group of strains (e.g. a
clade) afterwards in R.

For every pair of strains (i, j) the script accumulates, over all callable sites:
  D[i, j]  = number of differing allele pairs between i and j (4 comparisons per site)
  C[i, j]  = number of allele comparisons between i and j (4 per site both called)
and for every strain i (within-strain, i.e. heterozygosity):
  W[i]     = number of heterozygous callable sites
  Wc[i]    = number of callable sites
Pi for a set S of strains (2n alleles) is then
  (sum_{i<j in S} D_ij + sum_{i in S} W_i) / (sum_{i<j in S} C_ij + sum_{i in S} Wc_i)
which is the standard per-site pi with missing data handled per site (as in pixy).
D_sing / W_sing hold the contribution from panel-wide singleton sites (minor
allele count 1 across the 96 strains), so a singleton-excluded pi can be computed
by subtracting them.

Site filters
  * non-variant sites (ALT == '.'): kept; genotype callable if DP >= MIN_DP
  * SNP sites (REF and all ALT alleles single bases, no '*'): kept unless FILTER is
    LowQual (no GATK hard filters; the VCF is raw GenotypeGVCFs output);
    genotype callable if DP >= MIN_DP and GQ >= MIN_GQ
  * InDels / spanning deletions / LowQual SNP sites: dropped entirely
"""
import sys, os, numpy as np
from multiprocessing import Pool

VCF = sys.argv[1]
OUTDIR = sys.argv[2]
NPROC = int(sys.argv[3]) if len(sys.argv) > 3 else 16
MIN_DP, MIN_GQ = 10, 20
BATCH = 100_000
NA = -1

def index_chromosomes(path):
    """Byte offset of the first line of each chromosome, plus header sample names."""
    offsets, order, samples = {}, [], None
    with open(path, "rb") as f:
        pos, last = 0, None
        for line in f:
            if line[:1] == b"#":
                if line.startswith(b"#CHROM"):
                    samples = [s.decode() for s in line.rstrip(b"\n").split(b"\t")[9:]]
            else:
                tab = line.find(b"\t")
                chrom = line[:tab]
                if chrom != last:
                    offsets[chrom] = pos
                    order.append(chrom)
                    last = chrom
            pos += len(line)
    return offsets, order, samples

def accumulate(A, acc):
    """A: sites x n x 2 int8 allele matrix (-1 = not callable)."""
    called = (A[:, :, 0] >= 0)
    N = (called * 2).astype(np.float64)                     # sites x n
    same = np.zeros((N.shape[1], N.shape[1]))
    for k in range(int(A.max()) + 1):
        Ck = ((A[:, :, 0] == k).astype(np.float64) + (A[:, :, 1] == k))
        same += Ck.T @ Ck
    NN = N.T @ N
    acc["C"] += NN
    acc["D"] += NN - same
    het = called & (A[:, :, 0] != A[:, :, 1])
    acc["W"] += het.sum(0)
    acc["Wc"] += called.sum(0)

    # singleton sites: minor allele count == 1 among called alleles
    alt_counts = []
    total = (called.sum(1) * 2)
    maxk = int(A.max()) + 1
    counts = np.stack([((A[:, :, 0] == k).sum(1) + (A[:, :, 1] == k).sum(1)) for k in range(maxk)], 1)
    nonzero = (counts > 0).sum(1)
    sing = (nonzero == 2) & (counts.min(1, where=counts > 0, initial=10**9) == 1)
    if sing.any():
        As = A[sing]
        cs = (As[:, :, 0] >= 0)
        Ns = (cs * 2).astype(np.float64)
        same_s = np.zeros_like(same)
        for k in range(maxk):
            Ck = ((As[:, :, 0] == k).astype(np.float64) + (As[:, :, 1] == k))
            same_s += Ck.T @ Ck
        acc["D_sing"] += Ns.T @ Ns - same_s
        acc["W_sing"] += (cs & (As[:, :, 0] != As[:, :, 1])).sum(0)
    acc["n_sites"] += A.shape[0]
    acc["n_snp_sites"] += int((nonzero >= 2).sum())
    acc["n_singleton_sites"] += int(sing.sum())

def process_chrom(args):
    chrom, start, n = args
    acc = dict(C=np.zeros((n, n)), D=np.zeros((n, n)), D_sing=np.zeros((n, n)),
               W=np.zeros(n), Wc=np.zeros(n), W_sing=np.zeros(n),
               n_sites=0, n_snp_sites=0, n_singleton_sites=0,
               n_dropped_indel=0, n_dropped_filter=0)
    buf = np.full((BATCH, n, 2), NA, dtype=np.int8)
    b = 0
    fmt_cache = {}
    with open(VCF, "rb") as f:
        f.seek(start)
        for line in f:
            fields = line.rstrip(b"\n").split(b"\t")
            if fields[0] != chrom:
                break
            ref, alt = fields[3], fields[4]
            variant = alt != b"."
            if variant:
                alts = alt.split(b",")
                if len(ref) != 1 or any(len(a) != 1 or a == b"*" for a in alts):
                    acc["n_dropped_indel"] += 1
                    continue
                if fields[6] == b"LowQual":
                    acc["n_dropped_filter"] += 1
                    continue
            fmt = fields[8]
            if fmt not in fmt_cache:
                keys = fmt.split(b":")
                fmt_cache[fmt] = (keys.index(b"DP") if b"DP" in keys else None,
                                  keys.index(b"GQ") if b"GQ" in keys else None)
            idp, igq = fmt_cache[fmt]
            row = buf[b]
            row.fill(NA)
            for s, field in enumerate(fields[9:]):
                parts = field.split(b":")
                gt = parts[0]
                if len(gt) != 3 or gt[0:1] == b".":
                    continue
                try:
                    dp = int(parts[idp]) if idp is not None else 0
                except (ValueError, IndexError):
                    continue
                if dp < MIN_DP:
                    continue
                if variant:
                    try:
                        if igq is None or int(parts[igq]) < MIN_GQ:
                            continue
                    except (ValueError, IndexError):
                        continue
                row[s, 0] = gt[0] - 48
                row[s, 1] = gt[2] - 48
            b += 1
            if b == BATCH:
                accumulate(buf[:b], acc)
                b = 0
    if b:
        accumulate(buf[:b], acc)
    return chrom.decode(), acc

if __name__ == "__main__":
    os.makedirs(OUTDIR, exist_ok=True)
    offsets, order, samples = index_chromosomes(VCF)
    n = len(samples)
    print(f"{n} samples; chromosomes: {[c.decode() for c in order]}", flush=True)
    with Pool(NPROC) as pool:
        results = pool.map(process_chrom, [(c, offsets[c], n) for c in order], chunksize=1)
    total = None
    log = []
    for chrom, acc in results:
        log.append([chrom] + [acc[k] for k in ("n_sites", "n_snp_sites", "n_singleton_sites",
                                               "n_dropped_indel", "n_dropped_filter")])
        if total is None:
            total = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in acc.items()}
        else:
            for k, v in acc.items():
                total[k] += v
    np.savetxt(os.path.join(OUTDIR, "samples.txt"), np.array(samples), fmt="%s")
    for k in ("C", "D", "D_sing"):
        np.savetxt(os.path.join(OUTDIR, f"{k}.tsv"), total[k], fmt="%.0f", delimiter="\t")
    np.savetxt(os.path.join(OUTDIR, "within_strain.tsv"),
               np.column_stack([total["W"], total["Wc"], total["W_sing"]]),
               fmt="%.0f", delimiter="\t", header="W\tWc\tW_sing", comments="")
    with open(os.path.join(OUTDIR, "site_counts_per_chromosome.tsv"), "w") as fh:
        fh.write("chrom\tcallable_sites_kept\tsnp_sites\tsingleton_sites\tdropped_indel\tdropped_lowqual\n")
        for r in log:
            fh.write("\t".join(map(str, r)) + "\n")
    print("done", flush=True)
