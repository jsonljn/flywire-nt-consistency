import pandas as pd

from paths import GENERAL_SCAN_FULL, THREE_PATTERNS, ensure_output_dirs

ensure_output_dirs()

df = pd.read_csv(GENERAL_SCAN_FULL)

# Pattern 1: Categorical blind spot (HIST)
p1 = df[(df["is_categorical_blindspot"]) & (~df["agrees_with_mcns"])]

# Pattern 2: ORN SER-vs-ACH confusion
orn = df[df["fafb_cell_type"].astype(str).str.startswith("ORN_")]
p2 = orn[(orn["mcns_confirmed_nt"] == "ACH") & (orn["fafb_dominant_nt"] == "SER")]

# Pattern 3: Dm GABA/ACH-vs-GLUT confusion
dm = df[df["fafb_cell_type"].astype(str).str.match(r"^Dm\d", na=False)]
p3 = dm[(dm["mcns_confirmed_nt"] == "GLUT") & (dm["fafb_dominant_nt"].isin(["GABA", "ACH"]))]

print("PATTERN 1: Categorical blind spot (histamine)")
print(f"  {len(p1)} types: {p1['fafb_cell_type'].tolist()}")
print()
print("PATTERN 2: ORN serotonin confusion")
print(f"  {len(p2)} of {len(orn)} ORN types: {p2['fafb_cell_type'].tolist()}")
print()
print("PATTERN 3: Dm glutamate confusion")
print(f"  {len(p3)} of {len(dm[dm['mcns_confirmed_nt'] == 'GLUT'])} GLUT-confirmed Dm types: {p3['fafb_cell_type'].tolist()}")

p1 = p1.copy()
p1["pattern"] = "categorical_blindspot_HIST"
p2 = p2.copy()
p2["pattern"] = "ORN_SER_confusion"
p3 = p3.copy()
p3["pattern"] = "Dm_GLUT_confusion"
combined = pd.concat([p1, p2, p3])
combined.to_csv(THREE_PATTERNS, index=False)
print(f"\nSaved {THREE_PATTERNS} ({len(combined)} total flagged cell types)")
