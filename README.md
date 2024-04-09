# DeepES
## Installation
hoge

## Program Usage
### Embed input protein sequences
```
python embed.py \
    --input_dir data \
    --output_dir output
```
- --input_dir: Please specify the path of your directory where the input data is located.
- --output_dir: Please specify the path to save outputs.
- --reuse: Whether to use the results of previous runs.
- --batch_size: Batch size specified for loading parameters.
- --cuda: Whether to use a GPU.
- --cpu_num: Number of threads.

### Calculate probabilities
```
python predict.py \
    --input_dir data \
    --output_dir output \
    --model_dir model \
    --rclass_list RC01053 RC00004,RC00014 RC01923
```

### Evaluate results
```
python evaluate.py \
    --input_dir data \
    --output_dir output \
    --rclass_list RC01053 RC00004,RC00014 RC01923
```