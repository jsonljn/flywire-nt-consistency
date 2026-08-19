import pandas as pd

from paths import DATA, FAFB_MERGED

print("Loading files...")
neurons = pd.read_csv(DATA / "neurons.csv")
cell_types = pd.read_csv(DATA / "cell_types.csv")
classification = pd.read_csv(DATA / "classification.csv")

print(f"neurons: {neurons.shape}")
print(f"cell_types: {cell_types.shape}")
print(f"classification: {classification.shape}")

print("\nneurons columns:", list(neurons.columns))
print("cell_types columns:", list(cell_types.columns))
print("classification columns:", list(classification.columns))

# Merge on root_id
merged = neurons.merge(cell_types, on='root_id', how='left')
merged = merged.merge(classification, on='root_id', how='left')

print(f"\nMerged shape: {merged.shape}")
print(f"Columns: {list(merged.columns)}")

# Check missing data
print(f"\nMissing primary_type: {merged['primary_type'].isna().sum()}")
print(f"Missing nt_type: {merged['nt_type'].isna().sum()}")
print(f"Missing hemilineage: {merged['hemilineage'].isna().sum()}")

merged.to_csv(FAFB_MERGED, index=False)
print(f"\nSaved {FAFB_MERGED}")

# Quick preview
print("\nSample rows:")
print(merged[['root_id', 'primary_type', 'nt_type', 'hemilineage', 'class']].head(10).to_string(index=False))
