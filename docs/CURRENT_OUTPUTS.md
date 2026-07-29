# Current pipeline outputs (pre-change inventory)

Snapshot of what the pipeline wrote **before** the warehouse plumbing change, recorded as
required by Task 1 of `PIPELINE_CHANGES.md`. Paths are relative to `<results_dir>` unless
stated otherwise. "Fate" records what happens to the path later in the same run, because
`organize_output()` moves nearly everything and then deletes four top-level directories.

Line numbers refer to `rnaseq.py` at commit `d34fe78`.

## Run driver

| Stage | Code | Writes | Fate |
|---|---|---|---|
| stdout/stderr capture | `setup_logging` (L45) | `processing_log.tsv` | Stays at root. Also receives per-sample failure rows from `log_processing_error` inside `process_sample` (L2088). |
| total runtime | `log_total_time` (L50) | appends to `processing_log.tsv` | — |

## Per-sample stages (`process_sample`, L2085)

| Stage | Code | Writes | Fate |
|---|---|---|---|
| fastp (no trimming) | `run_fastp_no_trimming_checked` (L2183) | `fastp/<sample>/<sample>.html`, `fastp/<sample>/<sample>.json`, `terminal_output_<sample>.tsv` | `fastp/` → `1_Preprocessing/02_quality_control/fastp`; `terminal_output_*.tsv` → `3_Results/09_documents`. **Never parsed into a results table**, though `run_rmats` reads the JSON for read length. |
| FastQC | `run_fastqc` (L78) | `fastqc/<sample>/<fastq_base>_fastqc.{html,zip}` | → `1_Preprocessing/02_quality_control/fastqc`. Module PASS/WARN/FAIL never extracted. |
| STAR | `run_star` (L98) | `star/<sample>/<sample>_Aligned.sortedByCoord.out.bam`, `_Log.final.out`, `_Log.out`, `_Log.progress.out`, `_SJ.out.tab`, `_Chimeric.out.junction`, `_Signal.{Unique,UniqueMultiple}.str{1,2}.out.wig`, `__STARgenome/`, `__STARtmp/` | See below. |
| STAR post-process | `post_process_alignment` (L164) | `.bam.bai` | `__STARgenome/` **deleted here**. |
| RSeQC strandedness | `infer_strandedness` (L195) | `rseqc/<sample>_Aligned.sortedByCoord.out_strandedness.txt` | → `1_Preprocessing/02_quality_control/rseqc`. Free-text file; the strandedness string is re-parsed by string matching in four separate places (`check_and_prepare_rmats`, `run_rmats`, `convert_wig_to_bigWig`, `generate_results_html`). |
| kallisto (species) | `run_kallisto` (L253) | `kallisto/<sample>/<sample>_<species>/{abundance.h5,abundance.tsv,run_info.json,<sample>.kallisto.log,<sample>.kallisto.err}` | → `2_Analysis/04_abundance_estimation/kallisto`. `abundance.tsv` survives on disk but **never reaches a results table**. |
| gene-level summarization | `summarize_to_gene_level` (L312) | `…/<sample>_<species>/abundance_gene_est_counts.tsv`, `abundance_gene_tpm.tsv` | Same move. Two columns only (`gene`, value); no `sample_id`. |
| kallisto (contaminants) | `run_kallisto(is_contaminant_run=True)` | `kallisto/<sample>/<sample>_Foreign_sequences/…` | Same move; read by `mycoplasma_check.py`. |

STAR fates in `organize_output` (L1851): `<sample>*.bw` → `2_Analysis/03_alignment/1_bw/{1_unique,2_unique_multiple}`,
`<sample>*.ba*` → `2_bam`, `<sample>*.tab` → `3_bed/1_tab`, `<sample>*.bed` → `3_bed/2_bed`,
`<sample>*.out` → `4_logs`. Then `star/` is **deleted** (L1893), which also discards
`_Chimeric.out.junction` (not matched by any glob) and `__STARtmp/`.

All `*.wig` (including the `.processed.wig` intermediates) are **deleted** by `remove_wig` (L1492).

## Differential expression

| Stage | Code | Writes | Fate |
|---|---|---|---|
| DESeq2 inputs | `generate_deseq2_metadata` (L448), `generate_counts_matrix` (L473) | `deseq2/deseq2_metadata.tsv`, `deseq2/deseq2_counts_matrix.tsv` | → `05_differential_expression/deseq2/prep_files`. Counts matrix is wide, one column per sample, indexed by biomart `gene` (a **symbol**, not an Ensembl ID). |
| DESeq2 | `deseq2_analysis_ercc.R` | `deseq2/deseq2_results.csv`, `deseq2/deseq2_countsNormalized_matrix.tsv`, `deseq2/PCA_plot.png`, `deseq2/ERCC_counts_{before,after}_normalization.png` | `.csv`/`.tsv` → `prep_files`; `.png` → `3_Results/08_figures`. `deseq2_results.csv` is the only place the raw statistics exist. |
| TPM collation | `collate_tpms` (L1222) | `deseq2/deseq2_collated_tpms.tsv` | → `prep_files`. Wide: one column per sample. |
| merge | `merge_deseq2_results` (L1260) | `deseq2/deseq2_results_genes.tsv` | → `prep_files`. **This is the problem file**: per-sample TPM columns renamed to `<sample> (TPMs)` and concatenated with the statistics columns. |
| XLSX | `reform_deseq2.py` | `deseq2/deseq2_results_genes.xlsx` | → `05_differential_expression/deseq2`. The de facto primary artifact. |
| volcano | `volcano.py` | `deseq2/volcano_plot.png` | → `08_figures`. |

`deseq2/` is **deleted** at L1893 after the moves. `generate_tpms_matrix` (L1192) writes
`deseq2/collated_tpms.tsv` but is **never called**.

## Transcript-level DE

| Stage | Code | Writes | Fate |
|---|---|---|---|
| Sleuth input | `create_sleuth_metadata` (L542) | `sleuth/sleuth_metadata.tsv` | → `sleuth/prep_files`. |
| Sleuth | `sleuth_analysis_ercc.R` | `sleuth/sleuth_DETranscripts_{lrt,wt}.csv`, `sleuth/sleuth_DEGenes_{lrt,wt}.csv` | → `prep_files`. The two gene-level files are written and then never read by anything. |
| TPM collation | `collate_transcript_tpms` (L1322) | `sleuth/sleuth_collated_transcripts_tpms.tsv` | → `prep_files`. Wide, one column per sample. |
| merge | `merge_sleuth_{wald,lrt}_results` (L1360, L1395) | `sleuth/sleuth_results_{wald,lrt}_transcripts.tsv` | → `prep_files`. Same wide-column problem as DESeq2; keeps only `gene`, `b`, `qval` (Wald) or `gene`, `qval` (LRT) — `pval`, `se_b`, `mean_obs`, `test_stat` are dropped. |
| XLSX | `reform_sleuth.py` | `sleuth/sleuth_results_{wald,lrt}_transcripts.xlsx` | → `05_differential_expression/sleuth`. |

`sleuth/` is **deleted** at L1893.

## Enrichment

| Stage | Code | Writes | Fate |
|---|---|---|---|
| GSEA input (≥3 reps) | `create_gsea_input_file` (L617), `create_gsea_cls_file` (L720) | `gsea/test.txt`, `gsea/test.cls` | → `07_enrichment_analysis/1_prep_files`. `test.txt` header is one column per sample. |
| GSEA preranked (<3 reps) | `create_rnk_file` (L791) | `gsea/deseq2_ranked_genes.rnk` | Moves with `gsea/`. |
| GSEA | `run_gsea` (L757) / `run_gsea_preranked` (L833) | `gsea/<gmt_stem>.Gsea.<timestamp>/` or `gsea/preranked_<gmt_stem>.GseaPreranked.<timestamp>/` per MSigDB collection | Whole `gsea/` dir → `2_Analysis/07_enrichment_analysis/gsea`. **No statistics are ever extracted** from `gsea_report_for_*.tsv`; the HTML report is the only consumer. |

## Splicing

| Stage | Code | Writes | Fate |
|---|---|---|---|
| BAM lists | `check_and_prepare_rmats` (L872) | `rmats/testBAMs.txt`, `rmats/controlBAMs.txt` | → `3_Results/09_documents`. `--b1` is **test**, `--b2` is **control**, so rMATS `SAMPLE_1`/`IncLevel1` is the test group. |
| rMATS | `run_rmats` (L950) | `rmats/rmats_out/{SE,A5SS,A3SS,MXE,RI}.MATS.{JC,JCEC}.txt`, `fromGTF.*.txt`, `{JC,JCEC}.raw.input.*.txt`, `summary.txt`, `rmats/rmats_tmp/` | `*.JCEC.txt` + `summary.txt` → `06_alternative_splicing/1_JCEC`; remaining `*.txt` (including all five **JC** files) → `2_otherOutput`; `rmats_tmp/` → `2_otherOutput`. Then `rmats/` is **deleted**. |
| bar charts | `generate_splicing_bar_chart` (L1056) | `rmats/rmats_out/barChart_rmatsASE_fdr{0.05,0.0005}.png` | → `08_figures`. Uses JCEC only, FDR-only significance (no `IncLevelDifference` cut-off), so its counts do not match `summary.txt`. |
| SJ → BED | `convert_SJtab_to_bed` (L1475) + `splice_junction_strand_separator.pl` | `star/<sample>/<sample>_SJ.out.tab.{positive,negative}.bed` | → `3_bed/2_bed`. Perl script emits `start-1` (0-based) and drops strand-0 (undefined) junctions entirely. |
| bigWig | `convert_wig_to_bigWig` (L1130) | `star/<sample>/…Signal.{Unique,UniqueMultiple}.str{1,2}.out.bw` (+ `.processed.wig`) | `.bw` → `1_bw/*`; `.processed.wig` deleted. **Strand is not in the filename**; `str1`/`str2` is preserved as-is and the only strand handling is sign inversion of the wig values. No normalization (`--outWigNorm None`). |

## Reports and summaries

| Stage | Code | Writes | Fate |
|---|---|---|---|
| mycoplasma | `mycoplasma_check.py` | `mycoplasma_report.tsv` | → `09_documents`. |
| top hits | `report_top_genes_and_transcripts` (L1509) | `topGenes.tsv`, `topTranscripts.tsv` | → `09_documents`. Inherits the wide `(TPMs)` columns. |
| top TPMs | `report_top_gene_tpms` (L1540) | `top_gene_tpms.tsv` | → `09_documents`. |
| alignment summary | `summarize_alignments` (L1570) | `alignmentSummary.tsv` | → `09_documents`. Wide: one row per sample, 7 hardcoded STAR metrics. |
| kallisto summary | `summarize_transcript_coverage` (L1624) | `transcriptCoverage.tsv` | → `09_documents`. |
| HTML | `generate_results_html` (L1644) | `results.html` | → `09_documents`. |
| R side effects | DESeq2/Sleuth scripts | `Rplots.pdf` | → `09_documents`. |
| metadata log | `metadata.py` | `metadata_log_<timestamp>.txt` | → `09_documents`. |
| metadata TSV | `metadata.py` | `<investigator>_metadata_<timestamp>.tsv` | **Not moved and not copied** — it sits in `<results_dir>` only because that is where the user was told to put it, and nothing in `rnaseq.py` guarantees it stays. |

## Post-organization renaming

`rename_rmats_files` (L1896), `rename_figures_files` (L1924), `rename_diffExp_output` (L1958)
and `rename_sample_files_and_folders` (L2004) rewrite filenames **after** the moves, prefixing
experiment name / condition-replicate and de-duplicating repeated sample prefixes. This is the
source of the archive's inconsistent naming: the final name depends on the order the renamers
ran and on which prefixes happened to already be present.

## Summary of gaps this change has to close

1. Statistics and per-sample expression are fused into wide, per-experiment column names
   (`deseq2_results_genes.tsv`/`.xlsx`, `sleuth_results_*_transcripts.tsv`, `test.txt`).
2. kallisto `abundance.tsv`, GSEA reports, Sleuth `pval`/`se_b`/`mean_obs`, the five rMATS
   **JC** tables, fastp/FastQC/STAR/RSeQC QC values and the metadata TSV never become
   queryable tables.
3. `str1`/`str2` survives into the final bigWig filenames with no strand resolution and no
   normalization.
4. Nothing records tool versions, genome build, annotation version or reference identity.
