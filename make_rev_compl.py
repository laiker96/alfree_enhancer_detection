from Bio import SeqIO
import sys

def reverse_complement(fasta_file):

    fasta_file = SeqIO.parse(fasta_file, 'fasta')
    with open('reverse_complement.fa', 'ta') as F:
        for j in fasta_file:
            rev_comp_seq = j.reverse_complement().seq.__str__()
            seq_id = j.id
            F.write('>' + seq_id + '_RC' + '\n')
            F.write(rev_comp_seq + '\n')

    return None

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    reverse_complement(fasta_file)

