import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


h5ad_path = sys.argv[1]
count_paths = sys.argv[2:]
adata = ad.read_h5ad(h5ad_path, backed="r")
obs_names = set(map(str, adata.obs_names))
report = {
    "h5ad_obs_examples": list(map(str, adata.obs_names[:10])),
    "replicate_examples": {},
    "replicate_counts": adata.obs["replicate"].astype(str).value_counts().sort_index().to_dict(),
    "count_matrices": [],
}
for field in ("cell_type", "leiden_res0_25"):
    if field in adata.obs:
        table = pd.crosstab(
            adata.obs[field].astype(str),
            adata.obs["replicate"].astype(str),
            dropna=False,
        )
        report[f"{field}_by_replicate"] = {
            str(label): {str(rep): int(count) for rep, count in row.items()}
            for label, row in table.iterrows()
        }
for replicate in sorted(adata.obs["replicate"].astype(str).unique()):
    names = adata.obs_names[adata.obs["replicate"].astype(str) == replicate]
    report["replicate_examples"][replicate] = list(map(str, names[:5]))

for path in count_paths:
    counts = sc.read_10x_h5(path)
    values = counts.X.data if sp.issparse(counts.X) else np.asarray(counts.X).ravel()
    raw_names = list(map(str, counts.obs_names))
    transforms = {
        "exact": set(raw_names),
        "append_-0": {f"{x}-0" for x in raw_names},
        "append_-1": {f"{x}-1" for x in raw_names},
        "append_-2": {f"{x}-2" for x in raw_names},
        "replace_terminal_1_with_0": {x[:-1] + "0" if x.endswith("1") else x for x in raw_names},
        "replace_terminal_1_with_2": {x[:-1] + "2" if x.endswith("1") else x for x in raw_names},
    }
    overlaps = {key: len(value & obs_names) for key, value in transforms.items()}
    best_transform = max(overlaps, key=overlaps.get)
    report["count_matrices"].append({
        "path": path,
        "shape": list(counts.shape),
        "obs_examples": raw_names[:5],
        "var_examples": list(map(str, counts.var_names[:5])),
        "feature_types": counts.var.get("feature_types", []).astype(str).value_counts().to_dict() if "feature_types" in counts.var else {},
        "min_nonzero": float(values.min()) if values.size else None,
        "max": float(values.max()) if values.size else None,
        "integer_like": bool(np.allclose(values, np.rint(values))),
        "barcode_overlap_by_transform": overlaps,
        "best_transform": best_transform,
        "best_overlap": overlaps[best_transform],
        "best_overlap_fraction_of_matrix": overlaps[best_transform] / len(raw_names),
    })
print(json.dumps(report, ensure_ascii=False, indent=2))
