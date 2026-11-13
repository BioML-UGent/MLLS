# Changelog, notes, to-dos

## 2026-2027
- New PC lab replacing the intro to neural networks (now covered in the separate deep learning course).
  Accompanies lecture 11 (t-SNE and UMAP).
- Dataset: human pancreas scRNA-seq of Baron et al. (2016), GEO GSE84133, 14 cell types, at most
  400 cells per type. Built by `make_pancreas.py` (needs scanpy; not for students).
- Chosen over PBMC3k because there the first two PCs already separate the cell types almost as well
  as t-SNE/UMAP (k-NN accuracy in 2D: 0.91 vs 0.93). On the pancreas data: PCA 0.54, t-SNE/UMAP 0.97.
- The SOLVED notebook runs in about 1.5 min on a multi-core laptop; expect a few minutes on Colab
  (t-SNE with perplexity 300 is the slowest step).
