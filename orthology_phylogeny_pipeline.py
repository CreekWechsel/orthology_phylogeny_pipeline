#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
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
"""

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(level=logging.INFO,
                        format='[%(levelname)s] %(message)s')

def run_command(cmd, cwd=None, exit_on_failure=True):
    """Run a shell command and handle errors."""
    logging.info(f"Running command: {' '.join(map(str, cmd))}")
    if cwd is not None:
        cwd = str(cwd)
    try:
        subprocess.run(cmd, cwd=cwd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}: {' '.join(map(str, cmd))}")
        if exit_on_failure:
            sys.exit(1)
        else:
            logging.warning("Continuing execution despite the error.")

def check_file_exists(file_path):
    """Check if a file exists."""
    if not Path(file_path).exists():
        logging.error(f"File not found: {file_path}")
        sys.exit(1)

def check_directory_exists(dir_path):
    """Check if a directory exists, if not, create it."""
    dir_path = Path(dir_path)
    if not dir_path.exists():
        logging.info(f"Creating directory: {dir_path}")
        dir_path.mkdir(parents=True, exist_ok=True)

def run_orthofinder(pep_dir, orth_cpu, orth_work):
    """Step 1.0: Run OrthoFinder."""
    logging.info("Step 1.0: Run OrthoFinder")
    orthofinder_ok = orth_work / "0.orthofinder.ok"
    if not orthofinder_ok.exists():
        cmd = ["orthofinder", "-f", str(pep_dir), "-M", "msa", "-a", str(orth_cpu), "-o", str(orth_work)]
        run_command(cmd)
        orthofinder_ok.touch()
    else:
        logging.info("OrthoFinder has already been run. Skipping this step.")

def select_single_copy_genes(orth_work):
    """Step 1.1: Select single-copy genes."""
    logging.info("Step 1.1: Select single-copy genes")
    single_copy_ok = orth_work / "1.single_copy.ok"
    if not single_copy_ok.exists():
        results_dirs = list(orth_work.glob("Results_*"))
        if not results_dirs:
            logging.error("No Results_* directories found in orthologs directory.")
            sys.exit(1)
        results_dir = results_dirs[0]
        orthogroups_txt = results_dir / "Orthogroups" / "Orthogroups_SingleCopyOrthologues.txt"
        orthogroups_tsv = results_dir / "Orthogroups" / "Orthogroups.tsv"
        check_file_exists(orthogroups_txt)
        check_file_exists(orthogroups_tsv)

        output_file = orth_work / "single.copy.genes.txt"
        single_copy_orthogroups = set()
        with open(orthogroups_txt, 'r') as f:
            for line in f:
                og = line.strip()
                if og:
                    single_copy_orthogroups.add(og)

        with open(orthogroups_tsv, 'r') as infile, open(output_file, 'w') as outfile:
            header = infile.readline()
            outfile.write(header)
            for line in infile:
                parts = line.strip().split('\t')
                og_id = parts[0]
                if og_id in single_copy_orthogroups:
                    outfile.write(line)

        single_copy_ok.touch()
    else:
        logging.info("Single-copy genes have already been selected. Skipping this step.")

def split_protein_sequences(pep_dir, orth_work, align_dir, orth_cpu):
    """Step 2: Split peptides with single-copy gene sets and run PRANK."""
    logging.info("Step 2: Split peptides and run PRANK")
    check_directory_exists(align_dir)

    split_pep_dir = align_dir / "2.1.split_protein"
    split_pep_ok = align_dir / "split_protein.ok"
    if not split_pep_ok.exists():
        logging.info("Step 2.1: Split protein sequences")
        check_directory_exists(split_pep_dir)

        protein_files = list(pep_dir.glob("*.fasta"))
        protein_dict = {}
        for pep_file in protein_files:
            species_name = pep_file.stem
            with open(pep_file, 'r') as f:
                seq_id = None
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        seq_id = line[1:].split()[0]
                        protein_dict[seq_id] = ''
                    else:
                        protein_dict[seq_id] += line

        single_copy_file = orth_work / "single.copy.genes.txt"
        check_file_exists(single_copy_file)

        with open(single_copy_file, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split('\t')
                og_id = parts[0]
                genes = parts[1:]
                with open(split_pep_dir / f"{og_id}.fasta", 'w') as outfile:
                    for gene in genes:
                        gene_ids = gene.split(', ')
                        for gid in gene_ids:
                            gid = gid.strip()
                            if gid in protein_dict:
                                outfile.write(f">{gid}\n{protein_dict[gid]}\n")
                            else:
                                logging.warning(f"Gene ID {gid} not found in protein sequences.")

        split_pep_ok.touch()
    else:
        logging.info("Protein sequences have already been split. Skipping this step.")

    aligned_dir = align_dir / "2.2.aligned_AA_prank"
    aligned_ok = align_dir / "aligned_AA.ok"
    if not aligned_ok.exists():
        logging.info("Step 2.2: Run PRANK on split protein sequences")
        check_directory_exists(aligned_dir)

        prank_cmd_list = []
        split_pep_files = list(split_pep_dir.glob("*.fasta"))
        for fasta_file in split_pep_files:
            output_prefix = aligned_dir / f"{fasta_file.stem}_Prank_aligned"
            cmd = f"prank -protein -d=\"{fasta_file}\" -o=\"{output_prefix}\""
            prank_cmd_list.append(cmd)

        parafly_cmd_file = align_dir / "ParaFly.run.prank.command.list.txt"
        with open(parafly_cmd_file, 'w') as f:
            for cmd in prank_cmd_list:
                f.write(cmd + '\n')

        cmd = ["ParaFly", "-c", str(parafly_cmd_file), "-CPU", str(orth_cpu)]
        run_command(cmd)

        aligned_ok.touch()
    else:
        logging.info("PRANK alignment has already been completed. Skipping this step.")

    filtered_dir = align_dir / "2.3.filter_aligned_AA_prank"
    filter_ok = align_dir / "filter_aligned.ok"
    if not filter_ok.exists():
        logging.info("Step 2.3: Filter PRANK results")
        check_directory_exists(filtered_dir)

        aligned_files = list(aligned_dir.glob("*_Prank_aligned.best.fas"))
        logging.info(f"Found {len(aligned_files)} PRANK alignment files to filter.")

        for aln_file in aligned_files:
            seq_dict = {}
            with open(aln_file, 'r') as f:
                seq_id = None
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        seq_id = line[1:]
                        seq_dict[seq_id] = ''
                    else:
                        seq_dict[seq_id] += line
            flag = False
            for seq_id, seq in seq_dict.items():
                if seq.startswith(' ') or seq.startswith('\t') or seq.startswith('\n'):
                    flag = True
                    break
            if not flag and len(seq_dict) > 0:
                with open(filtered_dir / aln_file.name, 'w') as outfile:
                    for seq_id, seq in seq_dict.items():
                        outfile.write(f">{seq_id}\n{seq}\n")
        filter_ok.touch()
    else:
        logging.info("Filtered PRANK results already exist. Skipping this step.")

    gene_align_ok = align_dir / "2.gene_align.ok"
    if not gene_align_ok.exists():
        if split_pep_ok.exists() and aligned_ok.exists() and filter_ok.exists():
            gene_align_ok.touch()
            logging.info("All steps in gene alignment completed. '2.gene_align.ok' created.")
        else:
            logging.warning("Not all steps in gene alignment are completed yet.")
    else:
        logging.info("'2.gene_align.ok' already exists.")

def construct_phylogenetic_tree(construct_tree_dir, align_dir, spe_num, cpu, species_list, ref, astral):
    """Step 3: Construct phylogenetic tree."""
    logging.info("Step 3: Construct phylogenetic tree")
    check_directory_exists(construct_tree_dir)

    trimal_dir = construct_tree_dir / "3.1.trimal_pep"
    trimal_ok = construct_tree_dir / "trimal_pep.ok"
    if not trimal_ok.exists():
        logging.info("Step 3.1: Run trimAl to filter alignments")
        check_directory_exists(trimal_dir)

        filter_aligned_dir = align_dir / "2.3.filter_aligned_AA_prank"
        aligned_files = list(filter_aligned_dir.glob("*.fas"))
        trimal_cmd_list = []
        for aln_file in aligned_files:
            output_file = trimal_dir / f"trimal_{aln_file.name}"
            cmd = f"trimal -in {aln_file} -out {output_file} -nogaps"
            trimal_cmd_list.append(cmd)

        parafly_cmd_file = construct_tree_dir / "para_trim.txt"
        with open(parafly_cmd_file, 'w') as f:
            for cmd in trimal_cmd_list:
                f.write(cmd + '\n')

        cmd = ["ParaFly", "-c", str(parafly_cmd_file), "-CPU", str(cpu)]
        run_command(cmd)

        trimal_ok.touch()
    else:
        logging.info("trimAl filtering has already been done. Skipping this step.")

    trimal_filter_dir = construct_tree_dir / "3.2.trimal_pep_filter"
    trimal_filter_ok = construct_tree_dir / "trimal_pep_filter.ok"
    if not trimal_filter_ok.exists():
        logging.info("Step 3.2: Filter trimAl results")
        check_directory_exists(trimal_filter_dir)

        trimal_files = list(trimal_dir.glob("*.fas"))
        for trimal_file in trimal_files:
            seq_count = 0
            with open(trimal_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        seq_count += 1
            if seq_count == spe_num:
                run_command(["cp", str(trimal_file), str(trimal_filter_dir)])

        trimal_filter_ok.touch()
    else:
        logging.info("Filtered trimAl results already exist. Skipping this step.")

    pep_mer_dir = construct_tree_dir / "3.3.pep_merged"
    pep_merged_ok = construct_tree_dir / "pep_merged.ok"
    if not pep_merged_ok.exists():
        logging.info("Step 3.3: Merge single-copy genes")
        check_directory_exists(pep_mer_dir)

        merged_sequences = {sp: '' for sp in species_list}
        trimal_filter_files = list(trimal_filter_dir.glob("*.fas"))
        for trimal_file in trimal_filter_files:
            seq_dict = {}
            with open(trimal_file, 'r') as f:
                seq_id = None
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        seq_id = line[1:].split('_')[0][:10]
                        seq_dict[seq_id] = ''
                    else:
                        seq_dict[seq_id] += line
            for sp in species_list:
                sp_name = sp[:10]
                if sp_name in seq_dict:
                    merged_sequences[sp] += seq_dict[sp_name]
                else:
                    merged_sequences[sp] += '-' * len(next(iter(seq_dict.values())))

        merged_fasta = pep_mer_dir / "merged.aln.fasta"
        with open(merged_fasta, 'w') as f:
            for sp in species_list:
                sp_name = sp[:10]
                f.write(f">{sp_name}\n{merged_sequences[sp]}\n")

        trimmed_fasta = pep_mer_dir / "trim.merged.aln.fasta"
        cmd = ["trimal", "-in", str(merged_fasta), "-out", str(trimmed_fasta), "-nogaps"]
        run_command(cmd)

        phylip_file = pep_mer_dir / "trim.merged.aln.phy"
        seq_dict = {}
        with open(trimmed_fasta, 'r') as f:
            seq_id = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    seq_id = line[1:]
                    seq_dict[seq_id] = ''
                else:
                    seq_dict[seq_id] += line
        seq_count = len(seq_dict)
        seq_length = len(next(iter(seq_dict.values())))
        with open(phylip_file, 'w') as f:
            f.write(f"{seq_count} {seq_length}\n")
            for sp in species_list:
                sp_name = sp[:10]
                seq = seq_dict.get(sp_name, '-' * seq_length)
                f.write(f"{sp_name.ljust(10)} {seq}\n")

        pep_mer_dir_str = str(pep_mer_dir.resolve())
        if not os.path.exists(pep_mer_dir_str):
            logging.error(f"Directory does not exist: {pep_mer_dir_str}")
            sys.exit(1)
        os.chdir(pep_mer_dir_str)
        cmd = ["raxmlHPC-PTHREADS-AVX", "-T", str(cpu), "-f", "a", "-N", "100",
               "-m", "PROTGAMMAWAG", "-x", "12345", "-p", "12345", "-o", ref,
               "-s", "trim.merged.aln.phy", "-n", "trim.merged.aln.nwk"]
        run_command(cmd)
        os.chdir(str(construct_tree_dir.resolve()))

        pep_merged_ok.touch()
    else:
        logging.info("Merged single-copy genes have already been processed. Skipping this step.")

    supergene_tree_ok = construct_tree_dir / "3.supergene_tree.ok"
    if not supergene_tree_ok.exists():
        if trimal_ok.exists() and trimal_filter_ok.exists() and pep_merged_ok.exists():
            supergene_tree_ok.touch()
            logging.info("Supergene tree construction steps completed. '3.supergene_tree.ok' created.")
        else:
            logging.warning("Not all steps in supergene tree construction are completed yet.")
    else:
        logging.info("'3.supergene_tree.ok' already exists.")

    pep_phy_dir = construct_tree_dir / "3.4.aligned_pep_phy"
    aligned_pep_phy_ok = construct_tree_dir / "aligned_pep_phy.ok"
    if not aligned_pep_phy_ok.exists():
        logging.info("Step 3.4: Convert individual alignments to Phylip format")
        check_directory_exists(pep_phy_dir)

        trimal_filter_files = list(trimal_filter_dir.glob("*.fas"))
        for trimal_file in trimal_filter_files:
            seq_dict = {}
            with open(trimal_file, 'r') as f:
                seq_id = None
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        seq_id = line[1:].split('_')[0][:10]
                        seq_dict[seq_id] = ''
                    else:
                        seq_dict[seq_id] += line
            seq_count = len(seq_dict)
            seq_length = len(next(iter(seq_dict.values())))
            phylip_file = pep_phy_dir / (trimal_file.stem + ".phy")
            with open(phylip_file, 'w') as f:
                f.write(f"{seq_count} {seq_length}\n")
                for sp in species_list:
                    sp_name = sp[:10]
                    seq = seq_dict.get(sp_name, '-' * seq_length)
                    f.write(f"{sp_name.ljust(10)} {seq}\n")
        aligned_pep_phy_ok.touch()
    else:
        logging.info("Aligned peptides have already been converted to Phylip format.")

    out_tree_dir = construct_tree_dir / "3.5.out_pep_tree"
    out_pep_tree_ok = construct_tree_dir / "out_pep_tree.ok"
    if not out_pep_tree_ok.exists():
        logging.info("Step 3.5: Build gene trees using RAxML")
        check_directory_exists(out_tree_dir)

        pep_phy_files = list(pep_phy_dir.glob("*.phy"))
        parafly_cmd_file = construct_tree_dir / "parafly_run.construct.pep.txt"
        with open(parafly_cmd_file, 'w') as f:
            for idx, phy_file in enumerate(sorted(pep_phy_files)):
                gene_name = phy_file.stem
                output_name = f"{gene_name}_{idx}.nwk"
                work_dir = out_tree_dir / f"work_{gene_name}_{idx}"
                check_directory_exists(work_dir)

                cmd = (
                    f'raxmlHPC-PTHREADS-AVX -T 2 -f a -N 100 '
                    f'-m PROTGAMMAWAG -x 12345 -p 12345 -o {ref} '
                    f'-s "{phy_file.resolve()}" -n {output_name} -w "{work_dir.resolve()}"'
                )
                f.write(cmd + '\n')

        os.chdir(str(construct_tree_dir.resolve()))
        cmd = ["ParaFly", "-c", str(parafly_cmd_file.resolve()), "-CPU", str(cpu)]
        run_command(cmd, exit_on_failure=False)

        combined_tree_file = construct_tree_dir / "combined.all.tree.txt"
        with open(combined_tree_file, 'w') as outfile:
            for idx, phy_file in enumerate(sorted(pep_phy_files)):
                gene_name = phy_file.stem
                output_name = f"{gene_name}_{idx}.nwk"
                work_dir = out_tree_dir / f"work_{gene_name}_{idx}"
                tree_file = work_dir / f"RAxML_bestTree.{output_name}"
                if tree_file.exists():
                    with open(tree_file, 'r') as infile:
                        outfile.write(infile.read())
                else:
                    logging.warning(f"Tree file not found: {tree_file}")

        out_pep_tree_ok.touch()
    else:
        logging.info("Gene trees have already been constructed. Skipping this step.")

    astral_ok = construct_tree_dir / "astral.ok"
    if not astral_ok.exists():
        logging.info("Step 3.6: Run ASTRAL to infer species tree")

        astral_dir = construct_tree_dir / "3.6.astral_tree"
        check_directory_exists(astral_dir)

        combined_tree_file = construct_tree_dir / "combined.all.tree.txt"
        astral_output = astral_dir / "all_pep_astral.nwk"

        cmd = [
            "java", "-jar", str(astral),
            "-i", str(combined_tree_file),
            "-o", str(astral_output)
        ]
        run_command(cmd)
        astral_ok.touch()
    else:
        logging.info("ASTRAL analysis has already been completed. Skipping this step.")

def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="""Pipeline script to run orthology analysis and phylogenetic tree construction.

    Author: Chen Zijie
    Email: a604249194@126.com
    QQ: 2264177403
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pep_dir", required=True, help="Directory containing peptide sequences.")
    parser.add_argument("--cds_dir", required=True, help="Directory containing CDS sequences.")
    parser.add_argument("--species", required=True, nargs='+', help="List of species identifiers.")
    parser.add_argument("--ref", required=True, help="Reference species identifier.")
    parser.add_argument("--cpu", required=False, type=int, default=32, help="Number of CPUs to use.")
    parser.add_argument("--astral", required=False, default="/path/to/astral.jar", help="Path to ASTRAL JAR file.")
    args = parser.parse_args()

    main_dir = Path.cwd()
    species_list = args.species
    spe_num = len(species_list)
    ref = args.ref
    cpu = args.cpu
    orth_cpu = 3 * cpu

    args.pep_dir = Path(args.pep_dir).resolve()
    args.cds_dir = Path(args.cds_dir).resolve()
    orth_work = (main_dir / "1.orthologs").resolve()
    align_dir = (main_dir / "2.gene_align").resolve()
    construct_tree_dir = (main_dir / "3.phylogenetic_tree").resolve()
    args.astral = Path(args.astral).resolve()

    run_orthofinder(args.pep_dir, orth_cpu, orth_work)

    select_single_copy_genes(orth_work)

    split_protein_sequences(args.pep_dir, orth_work, align_dir, orth_cpu)

    construct_phylogenetic_tree(construct_tree_dir, align_dir, spe_num, cpu, species_list, ref, args.astral)

if __name__ == "__main__":
    main()
