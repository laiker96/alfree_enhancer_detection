from Bio import SeqIO
import sys
import os

fasta_file = sys.argv[1]


def translate_sequences(fasta_file):
    output_file = os.path.splitext(os.path.basename(fasta_file))[0]
    with open('translated_' + output_file + '.fa', 'at') as F:

        with open(fasta_file, 'rt') as f:
            
            seq = SeqIO.parse(f, 'fasta')
            t_seq = [s.translate(id = True, description = True) for s in seq]
            [SeqIO.write(t_s, F, 'fasta-2line') for t_s in t_seq]
            

    return None

translate_sequences(fasta_file)




