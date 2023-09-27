#!/usr/bin/env python3

# insert_gene.py - a script to insert sequences into a single strand genome
# (C) 2023

import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint


def read_sequence(fasta_file):
    logging.debug(f"loading sequence from {fasta_file}")
    with open(fasta_file, "r") as input_fh:
        for record in SeqIO.parse(input_fh, "fasta"):
            return record.seq, record.description
        
def insert_sequence(genome, gene, pos):
    len_genome = len(genome)
    len_gene = len(gene)
    logging.debug(f"Inserting gene ({len_gene} bp) in genome ({len_genome} bp) at position {pos}")
    return genome[0:pos] + gene + genome[pos:]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Path to Fasta file with sequence", required=True)
    parser.add_argument("-g", "--genome", help="Path to Fasta file with genome", required=True)
    parser.add_argument("-o", "--output", help="Mutated Fasta output file", required=True)
    parser.add_argument("-i", "--id", help="New sequence id", required=True)
    parser.add_argument("-p", "--pos", type=int, help="Insert location, if missing a random position is used", required=False)
    parser.add_argument("-v", "--verbose", help="Verbose mode", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    
    genome, header = read_sequence(args.genome)
    gene, gene_head = read_sequence(args.fasta)

    if args.pos:
        pos = args.pos
    else:
        pos = randint(0, len(genome))

    new_genome = insert_sequence(genome, gene, pos) 

    with open(args.output, "w") as output_fh:
        new_record = SeqRecord(
                        Seq(new_genome),
                        id=args.id,
                        name=args.id,
                        description=f"{header} | insert: {pos} | {gene_head}",
                    )
        SeqIO.write(new_record, output_fh, "fasta")
    
    logging.debug("Completed")

if __name__ == "__main__":
    main()
