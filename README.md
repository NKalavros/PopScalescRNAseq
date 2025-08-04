# Pop Scale scRNA-seq

This repository aims to evaluate popularly used methods for atlas integration of scRNA-seq data and juxtpose them to zero-shot (or slightly finetuned) scRNA-seq based foundation models for the same task.

It is split into three parts:

## Datasets:

A collection of at least 5 studies from multiple disease conditions + healthy samples in 3 tissues:

1. Breast (Wu, Pal, Kumar, Gray, ) (Should integrate https://www.nature.com/articles/s41597-022-01236-2 and https://www.biorxiv.org/content/10.1101/2025.03.13.643025v2.full)
2. Kidney (Bi, Krishna, Braun, Abedini, Lake)
3. Endometrium (Waiting on data processed by Fenyo lab)


## Conventional methods:

### Inter study:

1. Seurat Sketch (Removed since this remained a Python pipeline for now)
2. Symphony (Through symphonypy)
3. scANVI (Through scvi-tools)
4. scPoli (Through scarches)
5. CellHint (Through the teichlab implementation)

### Intra study:

1. Harmony (through harmonypy)
2. scVI (through scvi-tools)
3. scanorama (through scanorama)
4. ComBat (through ComBat py)

## scRNA-seq foundation model methods:

1. scGPT
2. scFoundation
3. CellFM
4. Geneformer
5. UCE
6. CellFM (This is quite hard to implement)

## Evaluation:

1. UMAP visualizations
2. scIB evaluation

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

