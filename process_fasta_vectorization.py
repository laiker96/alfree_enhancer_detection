#!/usr/bin/env python

import lzma
import pickle
from Bio import SeqIO
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from sklearn.feature_extraction.text import CountVectorizer as Cvec
from itertools import product
from scipy.stats import poisson
from scipy.special import softmax


dmel_bkg = np.load('dmel_bkg.npz.npy').reshape(1, 4096)
dmel_bkg_rev = np.load('dmel_bkg_revcomp.npz.npy').reshape(1, 4096)
dvir_bkg = np.load('dvir_bkg.npz.npy').reshape(1, 4096)
dvir_bkg_rev = np.load('dvir_bkg_revcomp.npz.npy').reshape(1, 4096)

def get_kmers (seq, k = 6) : 
    '''
    Compute kmer spectrum on DNA sequence
    '''
    return  [seq[x:x + k].lower() for x in range (len(seq) - k + 1)]

def tokenizer(kmers, single_seq=False):
    '''
    Create count table for kmer spectrum
    '''
    table = Cvec(vocabulary = [''.join(i) for i in product('acgt', repeat = 6)])
    if single_seq:
        table.fit([' '.join(kmers)])
        rv = table.transform([' '.join(kmers)]).toarray()
    else:
        table.fit(kmers)
        rv = table.transform(kmers).toarray()
    return rv

def compute_kmer_array(fasta_file, k = 6, window_size = 500, stride = 50, return_array = False):
    '''
    Compute kmer tables for fasta file
    '''
        
    file_basename = os.path.splitext(fasta_file)[0]
    out_name = file_basename + '.xz'
        
    region_keys = dict()

    with lzma.open(out_name, "ab") as F:
        
        with open(fasta_file, 'rt') as f:
            
            seqs = SeqIO.parse(f, 'fasta')
            for index, seq in enumerate(seqs):
                
                seq_name, seq_length, seq_string = seq.id, seq.seq.__len__(), seq.seq.__str__()
                region_keys[index] = (seq_name, seq_length)
                rv = []
                for  L in range(0, seq_length, stride):
                    print(seq_name)
                    if (L + window_size < seq_length):
                        r = seq_string[L:L+window_size]
                        r = get_kmers(r, k = k)
                        rv.append(r)
                for i in range (len(rv)): 
                    rv [i]  =  ' '.join(rv[i])
            
                x = tokenizer(rv)
                if return_array:
                    return x
                else:
                    pickle.dump(x, F)
                    id_table = open(file_basename + '.idTables', 'bw')
                    pickle.dump(region_keys, id_table)
                    id_table.close()

def calc_poisson_cdf_common(kmer_table, mu_kmer_array):
    times = kmer_table.shape[0]
    return np.where(kmer_table > 0, 1 - poisson.cdf(kmer_table - 1, np.tile(mu_kmer_array, (times, 1))), 1)
    # return np.where(kmer_table > 0, 1 - poisson.cdf(kmer_table, np.tile(mu_kmer_array, (times, 1))), 1)

def calc_poisson_cdf_all(kmer_table, mu_kmer_array):
    times = kmer_table.shape[0]
    return poisson.cdf(kmer_table, np.tile(mu_kmer_array, (times, 1)))
    # return poisson.cdf(kmer_table, np.tile(mu_kmer_array - 1, (times, 1)))


def test_func(id_tables = 'rregions.idTables', enhancer_dna = 'annot_Ubiquitous_enhancers_S10_hg38.fa'
        , regions_pickle = 'rregions.pickle', k = 6
        , dmel_bkg = None, dvir_bkg = None, dmel_bkg_rev = None):
    
    dmel_bkg = dmel_bkg * (500 - k + 1)
    dmel_bkg_rev = dmel_bkg_rev * (500 - k + 1)
    dvir_bkg = dvir_bkg * (500 - k + 1)

    basename = os.path.splitext(regions_pickle)[0]
    out_file = basename + '.bed'
    SO = open('predicted_enhancers.bed', 'at')

    with open(id_tables, 'br') as F:
        regions_keys = pickle.load(F)

    with open(enhancer_dna, 'rt') as F:
        bio_seqs = list(SeqIO.parse(F, 'fasta'))

    with lzma.open(regions_pickle, 'rb') as F:

        for i,j in regions_keys.items():

            seq_id = j[0]
            seq_length = j[1]
            rv = pickle.load(F)
            tmp = [seq for seq in bio_seqs if seq.id.split("_")[1] == seq_id]
            start_pos = int(seq_id.split(':')[1].split('-')[0])
            chrom = seq_id.split(':')[0]


            for seq in tmp:
                

                rc = seq.reverse_complement(id = True, name = True, description = True).seq.__str__()
                x = seq.seq.__str__()
                x = get_kmers(x, k = 6)
                x = tokenizer(x, single_seq = True)
                rc = get_kmers(rc, k = 6)
                rc = tokenizer(rc, single_seq = True)
                x_minima = np.minimum(rv, x)
                rc_minima = np.minimum(rv, rc)
                poisson_scores_x = calc_poisson_cdf_common(x_minima, dmel_bkg)
                poisson_scores_rc = calc_poisson_cdf_common(rc_minima, dmel_bkg_rev)
                # poisson_scores_region_x = calc_poisson_cdf_common(x_minima, dvir_bkg)
                # poisson_scores_region_rc = calc_poisson_cdf_common(rc_minima, dvir_bkg)

                # poisson_scores_region = calc_poisson_cdf_common(rv, dvir_bkg)
                
                poisson_scores_region_x = calc_poisson_cdf_common(x_minima, dvir_bkg)
                poisson_scores_region_rc = calc_poisson_cdf_common(rc_minima, dvir_bkg)
                ss_x = (poisson_scores_region_x) * (poisson_scores_x) 
                ss_rc = (poisson_scores_region_rc) * (poisson_scores_rc)

                poisson_scores_rc_all = calc_poisson_cdf_all(rc, dmel_bkg_rev)
                poisson_scores_x_all = calc_poisson_cdf_all(x, dmel_bkg)
                poisson_scores_rv_all = calc_poisson_cdf_all(rv, dvir_bkg)
                
                diff_x = np.abs(poisson_scores_x_all - poisson_scores_rv_all)
                diff_x = np.sum(diff_x, axis = 1)/4096
                diff_rc = np.abs(poisson_scores_rc_all - poisson_scores_rv_all)
                diff_rc = np.sum(diff_rc, axis = 1)/4096
                out_x = np.sum(1 - ss_x, axis = 1)/4096
                out_rc = np.sum(1 - ss_rc, axis = 1)/4096
                score_x = out_x - (diff_x * 0.1) 
                score_rc = out_rc - (diff_rc * 0.1)
                np.savetxt('fw_scores', score_x)
                np.savetxt('rv_scores', score_rc)

                # score_x = softmax(out_x - diff_x)
                # score_rc = softmax(out_rc - diff_rc)
                minima_x = np.nanargmax(score_x)
                minima_rc = np.nanargmax(score_rc)
                plt.plot(score_x)
                plt.plot(score_rc, 'r--')
                plt.show()
                if (score_x[minima_x] >= score_rc[minima_rc]):
                    minima = minima_x
                    out_score = score_x[minima_x]
                else:
                    minima = minima_rc
                    out_score = score_rc[minima_rc]
                ortho_position = (chrom + '\t' + f'{start_pos + minima*50}' + '\t' + f'{start_pos + minima*50 + 500}' + '\t' + f'{out_score}')
                
                enhancer_position = '_'.join(seq.id.split('_')[0:2])
                
                out_write = ortho_position + '\t' + enhancer_position + '\t' f'{minima}'
                SO.write(out_write + '\n')


    SO.close()
            

if __name__ == "__main__":
    compute_kmer_array(sys.argv[1])
    basename = os.path.splitext(sys.argv[1])[0]
    test_func(regions_pickle = basename + '.xz', id_tables = basename + '.idTables', dmel_bkg = dmel_bkg, dvir_bkg = dvir_bkg, dmel_bkg_rev=dmel_bkg_rev)



