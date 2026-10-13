# Comparative genomics of e14-like phages and prophages

This pipeline compares the cryptic e14 prophage of *Escherichia coli* K-12 with related prophages and phage genomes. It downloads genome records, estimates prophage boundaries using e14-associated *icd* homology, annotates extracted sequences, compares nucleotide similarity and gene organisation, and reconstructs relationships using protein families that include e14.

The workflow also compares the original MG1655 e14 CDS annotation with PHANOTATE calls on the same nucleotide sequence. The original annotation is included as a separate control track in a Clinker comparison, but is excluded from the main nucleotide comparisons and phylogeny.

## Requirements

- Linux with Bash, Conda and standard command-line utilities.
- Internet access for Conda installation, database downloads and NCBI retrieval.
- A tab-separated `manifest.tsv` in the repository root.
- Space for the Pharokka and Phold databases and intermediate outputs. Installation and annotation requirements depend on the database and software versions.

Run all commands from the repository root.

## TL;DR

1. Run `./00_setup.sh && ./01_download.sh  && ./02_find_boundaries.sh  && ./03_annotate.sh && ./04_compare.sh && ./05_finalize.sh`
2. Manually inspect gene homology groups identified by Clinker and resolve any naming conflicts.
3. Run `MODE=apply ./05_finalize.sh && ./06_phylogeny`.

## Scripts

| Script | Purpose |
| --- | --- |
| `00_setup.sh` | Install or update Conda environments, install annotation databases and record versions |
| `01_download.sh` | Download nucleotide FASTA and GenBank records from NCBI |
| `02_find_boundaries.sh` / `02_find_boundaries.py` | Infer candidate boundaries and extract sequences |
| `03_annotate.sh` | Annotate sequences with Pharokka/PHANOTATE and Phold |
| `04_compare.sh` / `04_pairwise_similarity.py` | Compare nucleotide similarity and generate Clinker plots |
| `compare_e14_annotations.py` | Compare original MG1655 CDS models with pipeline CDS calls |
| `05_finalize.sh` / `05_finalize.py` | Review homology-group products and write harmonised annotations |
| `05_augment_manifest_MG1655.py` | Add the original MG1655 annotation track to the stage-05 manifest |
| `05_relabel_clinker_reviewed.py` | Apply reviewed product labels to existing Clinker homology groups |
| `06_phylogeny.sh` / `06_phylogeny.py` | Discover reference phages and construct an e14-inclusive protein phylogeny |

## Input manifest

Use a header row with these columns, in this order. Several shell stages read columns by position.

| Column | Meaning |
| --- | --- |
| `label` | Unique identifier used in filenames; use simple identifiers without spaces or path separators |
| `organism` | Organism or strain description |
| `accession` | NCBI nucleotide accession, preferably including its version |
| `type` | `prophage` for host-embedded candidates or `free_phage` for complete phage records |
| `seed_start` | Approximate prophage start, using 1-based inclusive genomic coordinates |
| `seed_end` | Approximate prophage end, using 1-based inclusive genomic coordinates |
| `flank_bp` | Sequence included on each side of the seed for boundary searches; blank defaults to 10,000 bp |
| `boundary_mode` | `known_icd` for canonical K-12 e14; `infer` for candidates; `none` for complete phage records |

Exactly one row must have `boundary_mode=known_icd`. The default e14 label is `Ecoli_K12_e14`, accession `NC_000913.3`, with coordinates 1,196,209–1,211,412.

For complete phage genomes, leave seed coordinates blank and use `type=free_phage` and `boundary_mode=none`. Do not add `MG1655_curated_e14` to the main manifest: stage 05 creates a separate augmented manifest for this annotation control.

## Running the workflow

### 1. Set up and download

```bash
bash 00_setup.sh
bash 01_download.sh
```

Setup manages the `prophage_compare`, `pharokka`, `phold` and `clinker_py` environments. It reuses populated database directories unless `REFRESH_DBS=1` is supplied. Refreshing removes the existing annotation database directories before reinstalling them.

### 2. Infer and review boundaries

```bash
bash 02_find_boundaries.sh
```

Inspect `boundaries/boundaries.tsv`, candidate-hit tables and the extracted context before proceeding. Candidate sequences include the outer edges of the selected attachment-associated homology blocks. Canonical e14 retains its supplied coordinates; complete phage records retain their full sequence.

The method first searches for paired K-12 e14/*icd*-associated homologs, then falls back to direct boundary homology and finally to the seed interval. These are conservative boundary estimates, not experimentally confirmed attachment sites. A seed-interval fallback requires manual inspection of the genome.

### 3. Annotate and compare

```bash
bash 03_annotate.sh
bash 04_compare.sh
```

Stage 03 sanitises FASTA headers, calls genes with PHANOTATE through Pharokka, and uses Phold for structure-informed annotation. Stable downstream GenBank copies are placed in `annotation/final_gbk/`.

Stage 04 runs nucleotide BLAST in both directions for every pair. It greedily retains non-overlapping query High-scoring Segment Pairs (HSPs) by score and reports directional query coverage and alignment-length-weighted identity, together with symmetric summaries. The coverage calculation is conservative: an overlapping HSP is discarded in full. These summaries are descriptive, rather than phylogenetic distances.

Clinker comparisons include automatic and manifest ordering. The augmented comparison adds the original MG1655 e14 annotation beside the pipeline annotation. The annotation comparison reports coordinate matches, alternative starts/stops and possible splits or merges.

### 4. Review and harmonise product labels

```bash
bash 05_finalize.sh
```

Open `results/final/group_label_review.tsv`. Inspect rows with `review_required=1`, choose an appropriate `final_product`, and set `approved=yes`. Then run:

```bash
MODE=apply bash 05_finalize.sh
```

Apply mode writes harmonised GenBanks and final Clinker plots. Changed product annotations are retained in CDS notes. Product labels describe existing sequence-defined groups; assigning the same product to different groups does not establish homology.

Keep `CLINKER_IDENTITY` identical in stages 04 and 05. Regenerate and recheck the review table whenever the source session or input annotations change. Preserve a copy of any manually edited table before regenerating it.

### 5. Construct the protein phylogeny

```bash
bash 06_phylogeny.sh
```

Stage 06 creates a separate `prophage_phylogeny` Conda environment if absent. It:

1. Searches the INPHARED Mash sketch supplied with Pharokka using manifest entries of type `prophage` by default.
2. Downloads non-focal reference genomes within the Mash distance threshold.
3. Ranks references by their minimum discovery-seed distance, removes near-duplicates against already retained genomes, and caps the retained set. All manifest taxa are retained during this selection.
4. Predicts proteins consistently with PHANOTATE and clusters them with MMseqs2.
5. Selects single-copy families containing e14 that meet occupancy and median-length thresholds.
6. Aligns each family with MAFFT, trims with ClipKIT and concatenates the alignments, representing missing families as gaps.
7. Filters taxa by the number of represented families after trimming and estimates a partitioned maximum-likelihood tree with IQ-TREE, model selection/partition merging, 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates.

The final tree is unrooted. e14 is a family-selection anchor, not an outgroup. Some selected genomes may be excluded from the tree by the post-trimming coverage filter.

## Main outputs

| Directory or file | Contents |
| --- | --- |
| `genomes/` | Downloaded FASTA/GenBank records, download summary and checksums |
| `context/` | Host-genome context for boundary inspection |
| `boundaries/` | Attachment-associated reference, candidate hits and boundary summary |
| `prophages/` | Extracted sequences used downstream |
| `annotation/final_gbk/` | Pipeline-annotated GenBanks |
| `comparisons/blast/` | Pairwise similarity table and coverage/identity matrices |
| `comparisons/annotation_check/` | MG1655 versus PHANOTATE comparison, SVG tracks and curated e14 GenBank |
| `comparisons/clinker/` | Standard and augmented Clinker plots and sessions |
| `results/final/` | Review table, membership and presence/copy-count tables, harmonised annotations and final plots |
| `results/phylogeny/` | Reference-hit tables, retained/dropped taxa and provenance |
| `results/phylogeny/backbone/e14_inclusive/` | Selected families, alignments, supermatrix, partitions, taxon coverage and IQ-TREE outputs |
| `logs/` and `versions/` | Execution logs and software-version records |

The principal tree is:

```text
results/phylogeny/backbone/e14_inclusive/backbone_iqtree.treefile
```

Inspect `backbone_groups.tsv` and `backbone_taxon_coverage.tsv` alongside the tree.

## Parameters

Shell settings are supplied as environment variables, for example:

```bash
THREADS=16 bash 03_annotate.sh
MAX_MASH_DISTANCE=0.10 MAX_TOTAL=100 bash 06_phylogeny.sh
```

| Setting | Default | Purpose |
| --- | --- | --- |
| `THREADS` | `8` | Threads for tools that support the setting |
| `CLINKER_IDENTITY` | `0.3` | Minimum protein identity for Clinker links; use the same value in stages 04 and 05 |
| `PHOLD_GPU` | `0` | CPU mode by default; `1` requests the GPU/Foldseek-GPU mode |
| `DISCOVERY_TYPES` | `prophage` | Comma-separated manifest types used to recruit INPHARED references |
| `MAX_MASH_DISTANCE` | `0.15` | Maximum reference-discovery Mash distance |
| `REDUNDANCY_MASH` | `0.002` | References below this distance to a retained genome are removed |
| `MAX_TOTAL` | `150` | Target cap for focal plus reference genomes; focal genomes are always retained |
| `MASH_K` / `MASH_SKETCH` | `21` / `1000` | Query and candidate sketch settings; k must be compatible with the database sketch |
| `PROTEIN_MIN_ID` / `PROTEIN_COV` | `0.30` / `0.50` | MMseqs2 identity and bidirectional coverage thresholds |
| `E14_CORE_FRACTION` | `0.25` | Minimum family occupancy across retained genomes before tree-taxon filtering |
| `MIN_PROTEIN_AA` | `50` | Minimum median family length in amino acids |
| `MIN_CORE_GROUPS` | `4` | Minimum number of selected families |
| `E14_MIN_TAXON_GROUPS` | `3` | Minimum represented families per tree taxon after trimming |
| `ANCHOR_LABEL` | `Ecoli_K12_e14` | Required taxon in every selected family |

Boundary thresholds are command-line arguments:

```bash
bash 02_find_boundaries.sh --min-icd-identity 80 --min-icd-hit-length 100
```

Defaults are a 2,000-bp K-12 reference window, 75% minimum identity and 70-bp minimum reference-hit length. The generic fallback considers repeats up to 1,000 bp.

## Rerunning and reproducibility

Stages reuse existing outputs. File presence alone does not establish that an output matches the current sequences, parameters or database versions. After changing boundaries, annotation inputs or comparison settings, explicitly rebuild the affected stages:

```bash
FORCE=1 bash 03_annotate.sh
FORCE=1 bash 04_compare.sh
bash 05_finalize.sh
# Reinspect and approve the regenerated review table.
MODE=apply bash 05_finalize.sh
FORCE=1 bash 06_phylogeny.sh
```

If an accession changes under an existing label, refresh its downloads with `FORCE=1 bash 01_download.sh` before rerunning stage 02.

Record the manifest, software environments, database releases and analysis settings used for each analysis. Setup may update packages in existing environments, so version records alone do not pin a reproducible environment. Stage 06 records software versions, INPHARED filenames and its principal parameters in `results/phylogeny/software_versions.txt`.
