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
- --reuse: Whether to use the results of previous runs. default:`True`
- --batch_size: Batch size in embedding protein sequences. default:`1`
- --cuda: Whether to use a GPU. default:`False`
- --cpu_num: Number of threads. default:`1`

### Calculate probabilities
```
python predict.py \
    --input_dir data \
    --output_dir output \
    --model_dir model \
    --rclass_list RC01053 RC00004,RC00014 RC01923
```
- --input_dir: Same as above.
- --output_dir: Same as above.
- --model_dir: Please specify the path of your directory where the pre-trained model weights are located.
- --rclass_list: Please specify list of RClass corresponding to sequential enzyme reactions of interest.
- --batch_size: Batch size in calculating probabilities. default:`16`
- --cuda: Same as above. default:`False`
- --cpu_num: Same as above. default:`1`

### Evaluate results
```
python evaluate.py \
    --input_dir data \
    --output_dir output \
    --rclass_list RC01053 RC00004,RC00014 RC01923
```
- --input_dir: Same as above.
- --output_dir: Same as above.
- --rclass_list: Same as above.
- --window_size: Range of contiguous genes to be evaluated at a time. default:`10`
- --threshold: Threshold to obtain candidate genes. default:`0.99`
- --duplication: Whether to allow a single gene to be associated with multiple enzyme reactions. default:`False`