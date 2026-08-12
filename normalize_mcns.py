"""
Normalize MCNS neuron attributes to match the FAFB schema so analysis.py
can be reused without modification.
"""
import pandas as pd

print("Loading MCNS neurons.csv...")
df = pd.read_csv('data_mcns/neurons.csv')
print(f"Shape: {df.shape}")
print(f"Original columns: {list(df.columns)}")

# Map MCNS column names -> FAFB-style column names
rename_map = {
    'Root ID': 'root_id',
    'Predicted NT type': 'nt_type',
    'Predicted NT confidence': 'nt_type_score',
    'Primary Cell Type': 'primary_type',
    'Alternative Cell Type(s)': 'additional_type(s)',
    'Flow': 'flow',
    'Super Class': 'super_class',
    'Class': 'class',
    'Sub Class': 'sub_class',
    'Hemilineage': 'hemilineage',
    'Nerve': 'nerve',
    'Soma side': 'side',
}

df = df.rename(columns=rename_map)
print(f"\nRenamed columns: {list(df.columns)}")

# Check key columns
print(f"\nMissing root_id: {df['root_id'].isna().sum()}")
print(f"Missing nt_type: {df['nt_type'].isna().sum()}")
print(f"Missing primary_type: {df['primary_type'].isna().sum()}")

print(f"\nnt_type value counts:")
print(df['nt_type'].value_counts())

print(f"\nSample primary_type values:")
print(df['primary_type'].dropna().head(10).tolist())

df.to_csv('data_mcns/merged_annotations.csv', index=False)
print("\nSaved data_mcns/merged_annotations.csv")
