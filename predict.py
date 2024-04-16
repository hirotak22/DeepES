import argparse
import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

# Argument parser
def arguments():
    parser = argparse.ArgumentParser(description='Predict probability using RClass binary classifiers')
    parser.add_argument('--input_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--model_dir', required=True)
    parser.add_argument('--rclass_list', nargs='*', required=True)
    parser.add_argument('--batch_size', default=16)
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

# DNN architecture
class DNN(nn.Module):

    def __init__(self, input_dim, hidden_dim, output_dim, dropout):
        super().__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.ln = nn.LayerNorm(hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, output_dim)
        self.dropout = nn.Dropout(dropout)
        
    def forward(self, vec):
        hidden = self.dropout(self.ln(self.fc1(vec)))
        hidden = F.relu(hidden, inplace=True)
        hidden = self.fc2(hidden)
        return hidden

# Load pre-trained model parameters
def load_model(model_dir, rclass, input_dim=1280, hidden_dim=64, output_dim=1, dropout=0.1):
    model = DNN(input_dim, hidden_dim, output_dim, dropout)
    model.load_state_dict(torch.load(f'{model_dir}/model_{rclass}.pt'))
    return model

# Inference
def infer(model, embedding_vectors, batch_size, device):
    model.to(device)
    predicted_probability_list = []
    with torch.no_grad():
        for i in range(0, len(embedding_vectors), batch_size):
            vec = embedding_vectors[i:min(i+batch_size, len(embedding_vectors))]
            vec.to(device)
            logits = model(vec).reshape(-1)
            predicted_probability = torch.sigmoid(logits)
            predicted_probability_list.extend(predicted_probability.cpu().tolist())
    return np.array(predicted_probability_list)

def main():
    args = arguments()
    input_dir = args.input_dir
    output_dir = args.output_dir
    model_dir = args.model_dir
    rclass_list = args.rclass_list
    batch_size = args.batch_size
    cuda = args.cuda
    cpu_num = args.cpu_num
    
    os.makedirs(f'{output_dir}/inference', exist_ok=True)
    
    device = set_device(cuda, cpu_num)
    
    unique_rclass_list = list(set(sum([rclass_set.split(',') for rclass_set in rclass_list], [])))
    
    for input_fname in os.listdir(input_dir):
        fasta_name = input_fname.split('.')[0]
        embedding_vectors = torch.load(f'{output_dir}/embedding_vector/{fasta_name}.pt')
        gene_table = pd.read_table(f'{output_dir}/gene_table/{fasta_name}.tsv')
        ignore_index = gene_table[gene_table['ignore']].index.tolist()
        for rclass in unique_rclass_list:
            model = load_model(model_dir, rclass)
            model.eval()
            predicted_probability_list = infer(model, embedding_vectors, batch_size, device)
            if ignore_index != []:
                for idx in ignore_index:
                    np.insert(predicted_probability_list, idx, np.nan)
            np.save(f'{output_dir}/inference/{fasta_name}_{rclass}', predicted_probability_list)

if __name__ == '__main__':
    main()