import json
import sys

import anndata as ad
import numpy as np
import scipy.sparse as sp


def describe_matrix(matrix):
    result = {
        "shape": list(matrix.shape),
        "sparse": bool(sp.issparse(matrix)),
        "dtype": str(matrix.dtype),
    }
    if sp.issparse(matrix):
        values = matrix.data
        result["nnz"] = int(matrix.nnz)
    else:
        values = np.asarray(matrix).ravel()
    if values.size:
        result.update(
            min=float(np.nanmin(values)),
            max=float(np.nanmax(values)),
            nonnegative=bool(np.nanmin(values) >= 0),
            integer_like=bool(np.allclose(values[: min(values.size, 1_000_000)], np.rint(values[: min(values.size, 1_000_000)]))),
        )
    return result


path = sys.argv[1]
adata = ad.read_h5ad(path)
report = {
    "path": path,
    "shape": list(adata.shape),
    "obs_columns": list(adata.obs.columns),
    "var_columns": list(adata.var.columns),
    "layers": list(adata.layers.keys()),
    "obsm": list(adata.obsm.keys()),
    "uns_keys": list(adata.uns.keys()),
    "x": describe_matrix(adata.X),
    "obs_fields": {},
    "var_examples": [
        {
            "feature_id": str(index),
            **{str(column): str(value) for column, value in row.items()},
        }
        for index, row in adata.var.head(30).iterrows()
    ],
}
for column in adata.obs.columns:
    series = adata.obs[column]
    if series.nunique(dropna=False) <= 50:
        counts = series.astype(str).value_counts(dropna=False)
        report["obs_fields"][column] = {str(k): int(v) for k, v in counts.items()}
for layer in adata.layers.keys():
    report.setdefault("layer_details", {})[layer] = describe_matrix(adata.layers[layer])
print(json.dumps(report, ensure_ascii=False, indent=2))
