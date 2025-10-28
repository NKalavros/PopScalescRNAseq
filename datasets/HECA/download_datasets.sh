mkdir /gpfs/scratch/nk4167/EndometriumAtlas
cd /gpfs/scratch/nk4167/EndometriumAtlas
wget https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_cells_with_counts.h5ad
wget https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_nuclei.h5ad
mkdir HECA_cells
mkdir HECA_nuclei
# Use python to go back to counts from log-counts
python -c "import scanpy as sc; adata = sc.read_h5ad('endometriumAtlasV2_cells_with_counts.h5ad'); adata.X = adata.layers['counts']; adata.write('HECA_cells/data_cells.h5ad')"
# Nuclei do not have the counts layer, so reconstruct counts using stored totals.
python -c "
import numpy as np
import scanpy as sc
from scipy import sparse

adata = sc.read_h5ad('endometriumAtlasV2_nuclei.h5ad')
counts_key = 'n_counts' if 'n_counts' in adata.obs else 'total_counts'
if counts_key not in adata.obs:
	raise KeyError('AnnData.obs is missing n_counts/total_counts needed to recover raw counts.')

matrix = adata.X.copy()
if sparse.issparse(matrix):
	matrix.data = np.expm1(matrix.data)
else:
	matrix = np.expm1(matrix)

row_sums = np.asarray(matrix.sum(axis=1)).ravel()
total_counts = adata.obs[counts_key].to_numpy()
scale = np.divide(total_counts, row_sums, out=np.zeros(total_counts.shape, dtype=np.float64), where=row_sums != 0)

if sparse.issparse(matrix):
	matrix = matrix.multiply(scale[:, None])
else:
	matrix *= scale[:, None]

adata.X = matrix
print(adata.X.max(), adata.X.min(), adata.X.mean())
# Check that they are all sufficiently close to integers
diff = np.abs(adata.X.data - np.round(adata.X.data))
assert np.all(diff < 1e-2)
# convert to integers
adata.X.data = np.round(adata.X.data).astype(np.int32)
# CSR matrix conversion
adata.X = sparse.csr_matrix(adata.X)
adata.write('HECA_nuclei/data_nuclei.h5ad')"