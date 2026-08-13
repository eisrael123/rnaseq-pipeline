# Project: Database for Claude — Genomics Lab File Database

## Background
This project catalogs lab file outputs into a queryable DuckDB database. The lab studies Epstein-Barr virus (EBV) biology using cell line models including Akata (Burkitt lymphoma), DG75, SNU719, and HEK293-derived lines.

Data types present in the lab's file collection:
- RNA-seq splicing analysis (rMATS output: SE, A5SS, A3SS, MXE, RI event tables in JCEC format, plus per-comparison summary.txt files)
- ChIP-seq / CUT&RUN-style signal tracks (.wig, .bw bigwig files, stranded str1/str2) and transcription-factor motif scans (HOCOMOCO-style .bed motif files — AIRE_HUMAN, BACH1_HUMAN, CTCF_HUMAN, etc.)
- Differential motif density (md_scores) comparisons between test and control conditions, with accompanying plots (.png)
- WGBS methylation data (Bismark .cov coverage files)
- Splice junction tracks (SJ.out.tab, converted to .bed, split by strand)
- DESeq2 / differential expression results (.xlsx)

Directory roots (prefer local, fall back to the lab mount):
- Local (copied here for Docker; mounts are slow): `/Volumes/TUNGSACore3/rnaseq_runs`
    - Output folder of analysis (very impportant): '/Volumes/TUNGSACore3/rnaseq_runs/Flemington_RNAseq_analysis_output'
- Lab mount (source of truth if not copied locally): `/Volumes/FlemingtonLabMain1/2b_Flemington_Lab_Experiments/1_RNA_seq`

Subfolders are per experiment, e.g. `RNAseq/EBV_reactivation/Akata/Akata_BCR/`, `RNAseq/EBV_early_gene_transfections_in_DG75/BMRF2/`.

Goal: build and maintain a DuckDB (.duckdb) file that indexes these files (path, experiment, file type, cell line, condition, timepoint) so they're queryable with SQL instead of manually grepping directories.

Database skeleton (quick draft, not final): `/Volumes/TUNGSACore3/database_skeleton`. Use it as the starting point to walk processed final data under the directory roots above and build the DuckDB database.

## Working conventions
<!-- Add specifics here as they come up: naming conventions for tables/columns, how experiments should be grouped, which file types to index vs. ignore, etc. Keep this section updated so future sessions don't relitigate decisions already made. -->


