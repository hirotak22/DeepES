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

### Calculate probabilities
```
python predict.py \
    --input_dir data \
    --output_dir output \
    --model_dir /data2/shared/hirota/DREFONG/model/latest \
    --rclass_list RC01053 RC00004,RC00014 RC01923
```

### Evaluate results
```
python evaluate.py \
    --input_dir data \
    --output_dir output \
    --rclass_list RC01053 RC00004,RC00014 RC01923
```

## Options
hoge