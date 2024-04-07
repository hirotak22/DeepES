import argparse
import os
import pandas as pd
import torch
import esm
from Bio import SeqIO

# Argument parser
def arguments():
    parser = argparse.ArgumentParser(description='Embed protein sequences using ESM-2')
    parser.add_argument('--input_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--reuse', default=True)
    parser.add_argument('--batch_size', default=1)
    parser.add_argument('--cuda', default=False)
    parser.add_argument('--cpu_num', default=1)

    args = parser.parse_args()
    return args

# Configure settings for GPU and CPU
def set_device(cuda, cpu_num):
    if torch.cuda.is_available() and cuda:
        device = 'cuda:0'
    else:
        torch.set_num_threads(cpu_num)
        device = 'cpu'
    return device

# Load fasta file
def load_fasta(fpath):
    with open(fpath) as file:
        sequence_list = [[seq_record.id, seq_record.seq] for seq_record in SeqIO.parse(file, 'fasta')]
    return sequence_list

# Ignore sequences that are too long or contain "J"
def preprocess_sequences(sequence_list, max_len=4000):
    sequence_list_filtered = []
    gene_table = []
    for gene_id, seq in sequence_list:
        if len(seq) > max_len or 'J' in seq:
            gene_table.append([gene_id, True])
        else:
            sequence_list_filtered.append([gene_id, seq])
            gene_table.append([gene_id, False])
    return sequence_list_filtered, pd.DataFrame(gene_table, columns=['gene_id', 'ignore'])

# Embed input sequences using ESM-2
def embed(sequence_list, model, alphabet, batch_converter, batch_size, device):
    model.to(device)
    embedding_vectors = []
    with torch.no_grad():
        for i in range(0, len(sequence_list), batch_size):
            data = [sequence_list[j] for j in range(i, min(i+batch_size, len(sequence_list)))]
            batch_labels, batch_strs, batch_tokens = batch_converter(data)
            batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
            batch_tokens.to(device)
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
            token_representations = results["representations"][33]
            for j, tokens_len in enumerate(batch_lens):
                embedding_vectors.append(token_representations[j, 1 : tokens_len - 1].cpu().mean(0).clone())
    return torch.stack(embedding_vectors)

def main():
    
    args = arguments()
    input_dir = args.input_dir
    output_dir = args.output_dir
    reuse = args.reuse
    batch_size = args.batch_size
    cuda = args.cuda
    cpu_num = args.cpu_num
    
    device = set_device(cuda, cpu_num)
    
    # Load ESM-2 model 
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    
    for input_fname in os.listdir(input_dir):
        fasta_name = input_fname.split('.')[0]
        if reuse and os.path.isfile(f'{output_dir}/{fasta_name}.pt'):
            continue
        sequence_list = load_fasta(f'{input_dir}/{input_fname}')
        sequence_list_filtered, gene_table = preprocess_sequences(sequence_list)
        embedding_vectors = embed(sequence_list_filtered, model, alphabet, batch_converter, batch_size, device)
        torch.save(embedding_vectors, f'{output_dir}/embedding_vector/{fasta_name}.pt')
        gene_table.to_csv(f'{output_dir}/gene_table/{fasta_name}.tsv', sep='\t', header=True, index=False)

if __name__ == '__main__':
    main()