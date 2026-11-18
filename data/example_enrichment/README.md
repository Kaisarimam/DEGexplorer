# Example enrichment tables

Four tables in the layout the **Enrichment gallery → Import my enrichment
tables** option reads. Use them as a template for your own files, or import
them as they are to try the figures.

| File | Database |
|---|---|
| `go.BP.all.csv` | GO: Biological process |
| `go.CC.all.csv` | GO: Cellular component |
| `go.MF.all.csv` | GO: Molecular function |
| `pathway.all.csv` | KEGG pathway |

## Columns

| Column | Meaning |
|---|---|
| `ID` | term accession (`GO:0051607`, `hsa05160`) |
| `Description` | term name |
| `GeneRatio` | DEGs in the term / DEGs annotated in this database |
| `BgRatio` | background genes in the term / background genes annotated |
| `pvalue` | hypergeometric (over-representation) p-value |
| `p.adjust` | Benjamini–Hochberg adjusted p-value |
| `qvalue` | set equal to `p.adjust` in these files |
| `geneID` | the DEGs in the term, as symbols separated by `/` |
| `Count` | number of DEGs in the term |

This is clusterProfiler's layout, so tables saved from `enrichGO()` /
`enrichKEGG()` import without any editing. Each file holds every tested term
that contains at least one DEG; the gallery applies its own significance
cutoff.

## Where the numbers come from

Everything is computed from this repository's own demo data,
`data/example_deseq2_results.tsv.gz` (the bundled example): its 202 DEGs at adj. P ≤
0.05 and |log2FC| ≥ 1, tested against every gene with an adjusted p-value,
with a one-sided hypergeometric test (gene sets of 10–500 genes, BH
adjustment).

**`SOURCE.txt` in this folder states exactly which annotation release the
tables were built from.** Read it before quoting any number from them.

### Rebuilding from the current official release

```sh
Rscript data-raw/make_example_enrichment.R            # run from the repository root
```

The script needs an internet connection and base R only. It:

- downloads `go-basic.obo` and the human annotation file (GAF) from the **same
  dated Gene Ontology Consortium release** (`current.geneontology.org`, or an
  archived one with `--release YYYY-MM-DD`);
- **propagates annotations up the ontology** (is_a + part_of by default;
  `--relations` adds the regulates relations), so a gene annotated to a child
  term counts for every ancestor; NOT-qualified annotations are dropped,
  alternate ids are mapped to their primary term and obsolete terms removed;
- applies the **evidence-code choice** you give it: `--evidence all` (default,
  IEA and IBA included), `no-iea`, or `experimental`;
- takes **KEGG** pathways, gene links and names from the live KEGG REST API
  (academic use; `--no-kegg` skips this);
- **records the release** — the ontology data-version, the annotation file's
  date, the evidence and relation settings, and the KEGG release — in
  `SOURCE.txt`.

Already have the two GO files? Pass them with `--go-obo` and `--gaf`; they must
come from one release.

Because the tables come from the same demo study that loads on the Data tab,
the gene colours in the network figures (log2 fold change, matched by symbol)
line up with them.
