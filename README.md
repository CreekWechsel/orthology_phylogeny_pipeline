# orthology_phylogeny_pipeline
Pipeline script to run orthology analysis and phylogenetic tree construction

Name: phylogenetic_pipeline.py

Description:
This script implements an automated analysis pipeline for constructing species phylogenetic trees. Starting from protein sequence data of multiple species, it performs orthology analysis, single-copy gene selection, multiple sequence alignment, alignment filtering, gene tree and supergene tree construction, and finally uses ASTRAL to infer the species tree.

Functional Steps:
1. **Orthology Analysis**: Use OrthoFinder to identify orthologous gene families and single-copy genes.
2. **Single-copy Gene Extraction**: Extract the protein sequences corresponding to single-copy genes and save them as separate FASTA files.
3. **Multiple Sequence Alignment**: Perform multiple sequence alignment for each gene's sequences using PRANK.
4. **Alignment Filtering**: Use trimAl to remove gaps and low-quality regions from the alignments.
5. **Gene Tree Construction**: Construct gene trees for each gene using RAxML.
6. **Supergene Tree Construction**: Merge single-copy gene sequences to construct a supergene tree.
7. **Species Tree Inference**: Use ASTRAL to infer the species tree based on the gene trees.

Input:
- **Protein Sequence Directory** (`--pep_dir`): Contains the protein sequence files for each species (FASTA format).
- **CDS Sequence Directory** (`--cds_dir`): Contains the CDS sequence files for each species (optional).
- **Species List** (`--species`): A list of species identifiers participating in the analysis.
- **Reference Species** (`--ref`): The reference species used as the root for tree construction.
- **ASTRAL Program Path** (`--astral`): Path to the ASTRAL JAR file.

Output:
- **Supergene Tree**: The supergene tree constructed by RAxML, located in the `3.phylogenetic_tree/3.3.pep_merged` directory.
- **Gene Trees**: Gene trees for each single-copy gene, located in the `3.phylogenetic_tree/3.5.out_pep_tree` directory.
- **Species Tree**: The species tree inferred by ASTRAL, `all_pep_astral.nwk`, located in the `3.phylogenetic_tree` directory.

Usage:
```bash
python phylogenetic_pipeline.py \
    --pep_dir ./0.genomic_data/pep \
    --cds_dir ./0.genomic_data/cds \
    --species GGA MGA PHO SMI PCO CPI CMA LNY LSW \
    --ref GGA \
    --cpu 32 \
    --astral /path/to/astral.jar
