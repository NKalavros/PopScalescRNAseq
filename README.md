# Pop Scale scRNA-seq

This repository aims to evaluate popularly used methods for atlas integration of scRNA-seq data and juxtpose them to zero-shot (or slightly finetuned) scRNA-seq based foundation models for the same task.

It is split into three parts:

## Datasets:

A collection of at least 5 studies from multiple disease conditions + healthy samples in 3 tissues:

1. Breast
2. Kidney
3. Endometrium


## Conventional methods:

1. Seurat Sketch 
2. Symphony
3. scVI
4. scPoli
5. CellHint (?)

## scRNA-seq foundation model methods:

1. scGPT
2. scFoundation
3. CellFM
4. Geneformer
5. CellFM (This is quite hard to implement)

## Structure:

```
environments/
```

Conda envs to set up in BigPurple, they cross compile, but you need CUDA to actually run stuff

```
datasets/
```

Scripts that download and minimally process the data, needs only basic deps, such as cURL and scanpy. Checks for rapids to make your life faster and easier.

```
fm_scripts/
```

Scripts that generate embeddings per study with foundation models, in a zero-shot fashion for now. Work on anndata objects.

```
method_scripts/
```

Scripts that generate embeddings per study with conventional methods, with sensible defaults for now.


More notes to come.

