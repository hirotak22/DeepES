import argparse
import os
import numpy as np
import pandas as pd
import itertools
import pickle

# Argument parser
def arguments():
    parser = argparse.ArgumentParser(description='Evaluate DeepES output')
    parser.add_argument('--input_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--rclass_list', nargs='*', required=True)
    parser.add_argument('--window_size', default=10)
    parser.add_argument('--threshold', default=0.99)
    parser.add_argument('--duplication', default=False)

    args = parser.parse_args()
    return args

# Make probability matrix
def make_probability_matrix(output_dir, fasta_name, rclass_list):
    probability_matrix = []
    for rclass_set in rclass_list:
        synthetic_probability_array = None
        for rclass in rclass_set:
            predicted_probability_array = np.load(f'{output_dir}/inference/{fasta_name}_{rclass}.npy')
            if synthetic_probability_array is None:
                synthetic_probability_array = predicted_probability_array
            else:
                synthetic_probability_array *= predicted_probability_array
        if len(rclass_set) > 1:
            synthetic_probability_array **= 1.0/len(rclass_set)
        probability_matrix.append(synthetic_probability_array)
    return np.array(probability_matrix)

# Pick up genes to maximize window score
def maximize_window_score(probability_matrix_in_window, window_size, duplication):
    best_score = -1.0
    gene_idx_within_window = None
    reaction_num = len(probability_matrix_in_window)
    
    if duplication:
        score = 1.0
        idx_list = []
        for probability_array in probability_matrix_in_window:
            idx = np.nanargmax(probability_array)
            score *= probability_array[idx]
            idx_list.append(idx)
        score **= 1.0/reaction_num
        
        if best_score < score:
            best_score = score
            gene_idx_within_window = idx_list
    
    else:
        for permutation in itertools.permutations(list(range(window_size)), reaction_num):
            score = 1.0
            for idx, probability_array in zip(permutation, probability_matrix_in_window):
                score *= probability_array[idx]
            score **= 1.0/reaction_num
            
            if best_score < score:
                best_score = score
                gene_idx_within_window = permutation
    
    return best_score, gene_idx_within_window

# Carry out window-mapping from one end od the input to the other
def window_mapping(probability_matrix, window_size, duplication):
    input_length = len(probability_matrix[0])
    if input_length < window_size:
        print('Number of input genes < window_size')
        return None
    
    result = []
    for i in range(input_length-window_size+1):
        probability_matrix_in_window = probability_matrix[:,i:i+window_size]
        if np.isnan(probability_matrix_in_window[0]).sum() == window_size:
            result.append([np.nan, i, np.nan])
        else:
            best_score, gene_idx_within_window = maximize_window_score(
                                            probability_matrix_in_window,
                                            window_size,
                                            duplication
                                        )
            result.append([best_score, i, gene_idx_within_window])
    
    return result

# Get candidate genes with window score above threshold
def get_candidate_genes(window_mapping_result, threshold):
    if window_mapping_result is None:
        return None
    
    result_filtered = []
    for content in window_mapping_result:
        window_score, window_idx, gene_idx_within_window = content
        if not np.isnan(window_score) and threshold < window_score:
            result_filtered.append(content)
    return result_filtered

# Summarize candidate gene information
def summarize_candidate_genes(gene_table, rclass_list, window_mapping_result, threshold):
    result_filtered = get_candidate_genes(window_mapping_result, threshold)
    if result_filtered is None or result_filtered == []:
        return None
    else:
        summary = []
        for window_score, window_idx, gene_idx_within_window in result_filtered:
            gene_idx_list = [window_idx+idx for idx in gene_idx_within_window]
            gene_list = gene_table['gene_id'][gene_idx_list].to_list()
            summary.append([window_score, window_idx] + gene_list)
        summary_columns = ['score', 'window_idx'] + [','.join(rclass_set) for rclass_set in rclass_list]
        df_summary = pd.DataFrame(summary, columns=summary_columns)
        return df_summary

def main():
    args = arguments()
    input_dir = args.input_dir
    output_dir = args.output_dir
    rclass_list = args.rclass_list
    window_size = args.window_size
    threshold = args.threshold
    duplication = args.duplication
    
    rclass_list = [rclass_set.split(',') for rclass_set in rclass_list]
    
    for input_fname in os.listdir(input_dir):
        fasta_name = input_fname.split('.')[0]
        probability_matrix = make_probability_matrix(output_dir, fasta_name, rclass_list)
        result = window_mapping(probability_matrix, window_size, duplication)
        pickle.dump(result,open(f'{output_dir}/mapping_result/{fasta_name}.pkl', mode='wb'))
        
        gene_table = pd.read_table(f'{output_dir}/gene_table/{fasta_name}.tsv')
        df_summary = summarize_candidate_genes(gene_table, rclass_list, result, threshold)
        if df_summary is None:
            print(f'No candidate genes were found in {fasta_name}')
        else:
            df_summary.to_csv(f'{output_dir}/candidate_genes/{fasta_name}.tsv', sep='\t', header=True, index=False)

if __name__ == '__main__':
    main()