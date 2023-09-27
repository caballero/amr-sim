#!/usr/bin/env python3

# mutator.py - a script to mutate sequences (single strand genomes)
# (C) 2023

import argparse
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint

def mutate(base):
    if base == 'A':
        opt = ['c', 'g', 't', 't']
    elif base == 'T':
        opt = ['c', 'g', 'a', 'a']
    elif base == 'C':
        opt = ['a', 't', 'g', 'g']
    elif base == 'G':
        opt = ['a', 't', 'c', 'c']
    else:
        opt = ['a', 't', 'g', 'c']
    return opt[ randint(0,3) ]

def add_mutations(seq, mut):
    logging.debug(f"Adding {mut} mutations")
    aseq = list(seq)
    n = 0
    mutated = ['a', 'c', 'g', 't']
    while True:
        pos = randint(0, len(aseq) -1)
        old = aseq[pos]
        
        if old in mutated: # skip already mutated positions
            continue

        new = mutate(base=old)
        logging.debug(f"  iter {n}: changing {pos} ({old} -> {new})")
        aseq[pos] = new
        n += 1

        if n == mut:
            break

    return ''.join(aseq)

def main():
    rate = 0.01
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Path to Fasta file", required=True)
    parser.add_argument("-o", "--output", help="Mutated Fasta output file", required=True)
    parser.add_argument("-r", "--rate", type=float, help=f"Mutation rate [0.0 - 1.0], default = {rate}", required=False)
    parser.add_argument("-v", "--verbose", help="Verbose mode", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    
    if args.rate:
        rate = args.rate

    with open(args.fasta) as input_fh, open(args.output, "w") as output_fh:
        for record in SeqIO.parse(input_fh, "fasta"):
            mutations = int(rate * len(record.seq))
            mut_seq = add_mutations(seq=record.seq, mut=mutations)
            mut_record = SeqRecord(
                            Seq(mut_seq),
                            id=record.id,
                            name=record.name,
                            description=record.description + f" MUT_RATE={rate}",
                        )
            SeqIO.write(mut_record, output_fh, "fasta")
    
    logging.debug("Completed")

if __name__ == "__main__":
    main()
