# DEGexplorer 1.0.0

First public release.

## Explore results
- Load DESeq2, limma or edgeR result tables (including nf-core/differentialabundance
  output) or any table with a fold change and a p-value; DEGs at your thresholds.
- Volcano, MA, expression-scatter, PCA + scree and heatmap figures with full
  typography and geometry control; PNG / PDF / TIFF / SVG export.
- Gene explorer: every gene of the count matrix, plotted per sample or per group.
- Enrichment: GO (OrgDb or a bundled official GO release), KEGG and MSigDB with
  clusterProfiler, for human and 12 other organisms; an enrichment gallery of
  ready-made figures; gene intelligence (database links and literature).

## Run analyses (Analysis tab)
- DESeq2, limma-voom, edgeR and EBSeq from raw counts, with designs, covariates,
  comparisons, design checks, filtering, shrinkage, LRT, gene-length offsets and QC.
- Transcript-level DTE from Salmon / kallisto / RSEM: effective-length offsets
  (DESeq2), counts divided by quantification overdispersion (edgeR, limma),
  EBSeq Ng vectors, stageR stage-wise testing.
- Differential transcript usage (DTU): edgeR / limma diffSplice, DEXSeq, DRIMSeq
  and satuRn, each following its authors' workflow; isoform-usage plots.
- Every run exports a bundle with results, QC, a run record and an R script that
  reproduces the numbers exactly.

## Reproducibility
- `renv.lock` records the package versions DEGexplorer 1.0.0 was tested with
  (R 4.6.1, Bioconductor 3.23).
- `tests/smoke_test.R` checks the analysis code; it runs on every push.
