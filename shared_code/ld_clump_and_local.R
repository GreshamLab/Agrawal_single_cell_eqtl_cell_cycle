# =============================================================================
# Per-gene LD clumping and local/distant eQTL classification
# (sourced by the Figure 3 / Figure 4 Step1 classification Rmds)
#
# Added 2026-10-16. Replaces the earlier PLINK clumping, which
#   (1) clumped all genes together (the PLINK input had no gene column, so variants
#       of different genes were merged and an eQTL could be removed because a nearby
#       variant was a stronger eQTL for ANOTHER gene), and
#   (2) used a 2-kb window, shorter than the LD in the 96 strains.
#
# Updated 2026-10-21: clumping is done with PLINK 1.9 --clump, run separately for every gene
# (plink_clump_per_gene). The R implementation (clump_per_gene) is kept as a fallback; on the
# YPD.1 data it gave identical clumps for 1394 of 1422 genes with >1 variant (the rest differ
# only for pairs with r2 near 0.6: PLINK --clump estimates r2 from haplotype frequencies).
#
# Parameters (justified by the LD decay of the 96-strain SNPs, see
# Figure.3/.../codes/LD_decay_96_strains.Rmd: half-decay 0.7 kb, but pairs with
# r2 > 0.6 stay above background up to ~50 kb):
#   CLUMP_R2 = 0.6    variants of the same gene with r2 > 0.6 are one signal
#   CLUMP_KB = 50     ... if within 50 kb on the same chromosome
#   LOCAL_WINDOW_BP = 25000  local eQTL = any variant of the clump within 25 kb of
#                            the gene (either side), as in published yeast eQTL work
# Local / distant rules per variant type:
#   SNPs, InDels : any clump member within 25 kb of the gene
#   SVs          : the SV span (any clump member) within 25 kb of the gene;
#                  translocations (TRA): EITHER breakpoint within 25 kb
#                  (partner breakpoints from the per-strain calls in the SV VCF)
#   CNVs         : any clump member (placed at its gene start) within 25 kb of the gene,
#                  as for SNPs (since 2026-10-23; was: local only when CNV and gene are the same ORF)
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(BEDMatrix) })

CLUMP_R2        <- 0.6
CLUMP_KB        <- 50
LOCAL_WINDOW_BP <- 25000
SV_VCF <- "/projects/rps/cgsb/gresham/Akriti/wildyeast1/the_vcfs/structural_variant_vcf/SVcalls2/SVcalls/SV.1087Samples.vcf"
# PLINK 1.9 from the bioinformatics/20260810 module. Its `plink` is a wrapper that starts a
# container (~3 s per call), so all PLINK commands are written to one bash script and run in a
# single container session: PLINK_SHELL runs a bash script in which `plink` is on the PATH.
PLINK_SHELL <- "/share/apps/bioinformatics/20260810/run-bioinformatics.bash"

ROMAN16 <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")

# standardised genotypes (strains x variants) for the given variant IDs; missing -> variant mean
geno_z <- function(geno_prefix, ids) {
  bim <- fread(paste0(geno_prefix, ".bim"), header = FALSE)
  fam <- fread(paste0(geno_prefix, ".fam"), header = FALSE)
  G   <- BEDMatrix(paste0(geno_prefix, ".bed"), n = nrow(fam), p = nrow(bim))
  idx <- match(ids, bim$V2)
  if (anyNA(idx)) warning(sum(is.na(idx)), " variant(s) not in ", basename(geno_prefix), ".bim")
  ids <- ids[!is.na(idx)]; idx <- idx[!is.na(idx)]
  M <- as.matrix(G[, idx, drop = FALSE])
  Z <- apply(M, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE); s <- sd(x)
    if (is.na(s) || s == 0) rep(0, length(x)) else (x - mean(x)) / s })
  Z <- matrix(Z, nrow = nrow(M)); colnames(Z) <- ids
  Z
}

# Greedy clumping, separately for every gene (like PLINK --clump, but per gene):
# take the most significant remaining variant as index, merge every remaining variant of
# the same gene on the same chromosome within kb with r2 > r2, repeat.
# assoc needs Pheno, SNP, Chr, ChrPos, PValue. Returns one row per clump (index variant).
clump_per_gene <- function(assoc, Z, r2 = CLUMP_R2, kb = CLUMP_KB) {
  n <- nrow(Z)
  out <- lapply(split(assoc, assoc$Pheno), function(g) {
    g <- g[order(g$PValue), ]; g <- g[!duplicated(g$SNP), ]
    ok <- g$SNP %in% colnames(Z)
    R2 <- matrix(0, nrow(g), nrow(g)); diag(R2) <- 1
    if (sum(ok) > 1) {
      z <- Z[, g$SNP[ok], drop = FALSE]
      R2[ok, ok] <- (crossprod(z) / (n - 1))^2
    }
    remaining <- rep(TRUE, nrow(g)); res <- vector("list", nrow(g)); k <- 0
    for (i in seq_len(nrow(g))) {
      if (!remaining[i]) next
      mem <- which(remaining & g$Chr == g$Chr[i] & abs(g$ChrPos - g$ChrPos[i]) <= kb * 1000 & R2[i, ] > r2)
      mem <- union(i, mem); remaining[mem] <- FALSE; k <- k + 1
      res[[k]] <- data.frame(Pheno = g$Pheno[1], SNP = g$SNP[i], PValue = g$PValue[i],
                             n_members = length(mem), members = paste(g$SNP[mem], collapse = ";"))
    }
    do.call(rbind, res[seq_len(k)])
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  out
}

# PLINK 1.9 --clump, run separately for every gene (same output as clump_per_gene).
# Every significant variant may be an index or a member (--clump-p1 1 --clump-p2 1).
# The significant variants are first extracted into a small fileset; all PLINK commands run
# from one script (run_plink_clump.sh in work_dir, with the per-gene .assoc/.clumped/.log files).
# Genes with one variant need no PLINK run. Variants missing from the .bim cannot be clumped by
# PLINK; they are kept as their own clump.
plink_clump_per_gene <- function(assoc, geno_prefix, work_dir, r2 = CLUMP_R2, kb = CLUMP_KB,
                                 plink_shell = PLINK_SHELL) {
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  a <- as.data.table(assoc)[order(PValue)][!duplicated(paste(Pheno, SNP)), .(Pheno, SNP, PValue)]
  bim <- fread(paste0(geno_prefix, ".bim"), header = FALSE, select = 2)$V2
  miss <- unique(a$SNP[!a$SNP %in% bim])
  if (length(miss)) warning(length(miss), " variant(s) not in ", basename(geno_prefix),
                            ".bim; kept as single-variant clumps")
  a_in   <- a[SNP %in% bim]
  single <- a[!SNP %in% bim, .(Pheno, SNP, n_members = 1L, members = SNP)]
  multi  <- a_in[, .N, by = Pheno][N > 1, Pheno]
  single <- rbind(single, a_in[!Pheno %in% multi, .(Pheno, SNP, n_members = 1L, members = SNP)])

  sig <- file.path(work_dir, "signif_variants")
  writeLines(unique(a_in$SNP), paste0(sig, ".ids"))
  q <- shQuote
  cmds <- paste("plink --bfile", q(geno_prefix), "--extract", q(paste0(sig, ".ids")), "--make-bed",
                "--allow-extra-chr --allow-no-sex --out", q(sig), "> /dev/null || exit 1")
  for (ph in multi) {
    f_out <- file.path(work_dir, ph)
    fwrite(a_in[Pheno == ph, .(SNP, P = PValue)], paste0(f_out, ".assoc"), sep = "\t")
    unlink(paste0(f_out, ".clumped"))
    cmds <- c(cmds, paste("plink --bfile", q(sig), "--clump", q(paste0(f_out, ".assoc")),
                          "--clump-snp-field SNP --clump-field P --clump-p1 1 --clump-p2 1",
                          "--clump-r2", r2, "--clump-kb", kb,
                          "--allow-extra-chr --allow-no-sex --threads 1 --out", q(f_out), "> /dev/null"))
  }
  script <- file.path(work_dir, "run_plink_clump.sh")
  writeLines(c("#!/usr/bin/env bash", cmds), script)
  st <- system2(plink_shell, c("bash", script))
  if (st != 0) stop("PLINK run failed, see the logs in ", work_dir)

  clumped <- rbindlist(lapply(multi, function(ph) {
    f <- file.path(work_dir, paste0(ph, ".clumped"))
    if (!file.exists(f)) stop("PLINK --clump failed for ", ph, ", see ", sub(".clumped$", ".log", f))
    cl  <- fread(f)
    # SP2: clump members other than the index, "ID(1),ID(1),..." or "NONE"
    mem <- ifelse(cl$SP2 == "NONE", "", gsub("\\(1\\)", "", cl$SP2))
    data.table(Pheno = ph, SNP = cl$SNP, n_members = cl$TOTAL + 1L,
               members = ifelse(mem == "", cl$SNP, paste(cl$SNP, gsub(",", ";", mem), sep = ";")))
  }))
  out <- rbind(clumped, single)
  # full-precision p-values from the association file (PLINK rounds them)
  out <- merge(out, a, by = c("Pheno", "SNP"), sort = FALSE)
  as.data.frame(out[, .(Pheno, SNP, PValue, n_members, members)])
}

# Partner breakpoints of translocations: the per-strain calls merged into one TRA record
# (IDLIST_EXT) end with the partner breakpoint "chromosome<N>_<pos>"; it differs between
# strains, so all distinct partners are kept.
read_tra_partners <- function(ids, vcf = SV_VCF) {
  ids <- unique(ids[grepl("_TRA_", ids)])
  if (!length(ids)) return(data.table(SNP = character(), chr = integer(), pos = numeric()))
  rec <- fread(cmd = paste0("grep -v '^#' ", shQuote(vcf), " | awk -F'\\t' '$5==\"<TRA>\"' | cut -f3,8"),
               header = FALSE, sep = "\t", col.names = c("SNP", "INFO"))
  rec <- rec[SNP %in% ids]
  rbindlist(lapply(seq_len(nrow(rec)), function(i) {
    idl <- sub(".*IDLIST_EXT=([^;]*).*", "\\1", rec$INFO[i])
    ent <- strsplit(idl, ",")[[1]]
    m <- regmatches(ent, regexpr("chromosome[0-9]+(\\.[0-9]+)?_[0-9]+$", ent))
    if (!length(m)) return(NULL)
    unique(data.table(SNP = rec$SNP[i], chr = as.integer(sub("chromosome([0-9]+).*", "\\1", m)),
                      pos = as.numeric(sub(".*_([0-9]+)$", "\\1", m))))
  }))
}

# genomic interval(s) of a variant, from its ID and association-file position
variant_intervals <- function(snp, chr, pos, vt, tra = NULL) {
  if (vt == "InDels") {
    ref <- vapply(strsplit(snp, "_"), function(x) nchar(x[3]), 1L)
    return(data.table(SNP = snp, chr = chr, start = pos, end = pos + pmax(ref, 1L) - 1))
  }
  if (vt == "SVs") {
    tk  <- strsplit(snp, "_")
    sub <- vapply(tk, function(x) if (length(x) >= 2) x[2] else "", "")
    len <- vapply(tk, function(x) { i <- which(x == "len"); if (length(i)) abs(as.numeric(x[i + 1])) else NA_real_ }, 0)
    end <- ifelse(sub %in% c("INS", "TRA") | is.na(len), pos, pos + len)
    iv <- data.table(SNP = snp, chr = chr, start = pos, end = end)
    if (!is.null(tra) && nrow(tra)) {          # second breakpoint(s) of translocations
      p <- tra[SNP %in% snp[sub == "TRA"], .(SNP, chr, start = pos, end = pos)]
      iv <- rbind(iv, p)
    }
    return(iv)
  }
  # SNPs, and gene-based CNVs: ChrPos of a CNV marker is the start of its gene (the end is
  # not known for the non-reference YX genes, so every CNV is placed at its gene start)
  data.table(SNP = snp, chr = chr, start = pos, end = pos)
}

# Adds local_distant (+ clump size/members and the closest distance of any clump member
# to the gene). df: annotated index eQTLs with Pheno, Pheno_chr (roman), Pheno_start, Pheno_end.
# assoc: the full association table of that variant type (positions of clump members).
add_local_distant <- function(df, clumps, assoc, vt, tra = NULL, window = LOCAL_WINDOW_BP) {
  cl <- as.data.table(clumps)[, .(Pheno, SNP, clump_n_members = n_members, clump_members = members)]
  df <- merge(as.data.table(df), cl, by = c("Pheno", "SNP"), all.x = TRUE, sort = FALSE)
  # CNVs now use the same 25-kb rule as the other variant types (2026-10-23); was:
  # if (vt == "CNVs") {
  #   same_orf <- mapply(function(mem, ph) {
  #     g <- sub("_Cov[0-9.]+$", "", strsplit(mem, ";")[[1]])
  #     any(gsub("-", ".", g) == ph)
  #   }, df$clump_members, df$Pheno)
  #   df[, `:=`(local_distance_bp = NA_real_, local_distant = ifelse(same_orf, "Local", "Distant"))]
  #   return(as.data.frame(df))
  # }
  pos <- unique(as.data.table(assoc)[, .(Pheno, SNP, Chr, ChrPos)])
  long <- df[, .(member = strsplit(clump_members, ";")[[1]]), by = .(Pheno, SNP)]
  long <- merge(long, pos[, .(Pheno, member = SNP, Chr, ChrPos)], by = c("Pheno", "member"))
  iv <- variant_intervals(long$member, long$Chr, long$ChrPos, vt, tra)
  iv <- merge(unique(long[, .(Pheno, SNP, member)]), unique(iv), by.x = "member", by.y = "SNP",
              allow.cartesian = TRUE)
  gl <- unique(df[, .(Pheno, SNP, Pheno_chr, Pheno_start, Pheno_end)])
  iv <- merge(iv, gl, by = c("Pheno", "SNP"))
  iv[, var_chr := ROMAN16[chr]]
  iv[, d := ifelse(var_chr == Pheno_chr, pmax(0, Pheno_start - end, start - Pheno_end), NA_real_)]
  best <- iv[, .(local_distance_bp = if (all(is.na(d))) NA_real_ else min(d, na.rm = TRUE)), by = .(Pheno, SNP)]
  df <- merge(df, best, by = c("Pheno", "SNP"), all.x = TRUE, sort = FALSE)
  df[, local_distant := ifelse(!is.na(local_distance_bp) & local_distance_bp <= window, "Local", "Distant")]
  as.data.frame(df)
}
