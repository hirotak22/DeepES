# DeepES
[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

DeepES is a deep learning-based framework for enzyme screening.
DeepES can identify orphan enzyme candidate protiens by focusing on biosynthetic gene clusters and KEGG RClass.
DeepES uses protein sequences as inputs and evaluate whether the input genes contain biosynthetic gene clusters of interest.

## Requirements
- python 3.9.19 (with following packages)
  - numpy 1.26.4
  - pandas 2.2.1
  - pytorch 1.13.0
  - biopython 1.83
  - fair-esm 2.0.0

By using `environment.yml`, you can build an anaconda environment exactly the same as this research.
```
conda env create -f environment.yml
conda activate deepes
```

## Pretrained parameters
The model weights are available at [Zotero](https://doi.org/10.5281/zenodo.11123900).
Please download files and unzip.

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

## License
DeepES is released under the [MIT License](LICENSE).