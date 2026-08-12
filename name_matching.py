"""
Robust cell type name matching.

Different datasets name the same cell type differently -- e.g. FAFB's "R1-6"
vs MCNS's "R1-R6" (abbreviated vs. spelled-out range notation), or simple
case/whitespace differences. This module normalizes names so equivalent
types can be matched reliably, without resorting to risky fuzzy/edit-distance
matching that could silently produce a wrong correction.

Only two kinds of equivalence are treated as safe to auto-match:
1. Case/whitespace-insensitive exact match
2. Range notation: "R1-6" <-> "R1-R6" (abbreviated vs. full range)

Anything else is left unmatched rather than guessed at.
"""
import re


def canonical_name(name):
    """Return a canonical form for case/whitespace-insensitive matching."""
    return str(name).strip().lower()


def expand_range_notation(name):
    """
    'R1-6' (abbreviated range) -> ('R', '1', '6')
    'R1-R6' (full range)       -> ('R', '1', '6')
    Anything else -> None
    """
    name = str(name).strip()

    # Abbreviated: PREFIX + start_num + '-' + end_num, e.g. R1-6
    m = re.match(r'^([A-Za-z]+)(\d+)-(\d+)$', name)
    if m:
        prefix, start, end = m.groups()
        return (prefix.lower(), start, end)

    # Full: PREFIX + start_num + '-' + PREFIX + end_num, e.g. R1-R6
    m = re.match(r'^([A-Za-z]+)(\d+)-([A-Za-z]+)(\d+)$', name)
    if m:
        prefix1, start, prefix2, end = m.groups()
        if prefix1.lower() == prefix2.lower():
            return (prefix1.lower(), start, end)

    return None


def build_match_index(names):
    """
    Given a list of cell type names, build lookup indexes for:
    - canonical (case/whitespace-insensitive) form -> original name
    - range-notation canonical form -> original name
    """
    canon_idx = {}
    range_idx = {}
    for n in names:
        canon_idx[canonical_name(n)] = n
        r = expand_range_notation(n)
        if r:
            range_idx[r] = n
    return canon_idx, range_idx


def find_match(query_name, canon_idx, range_idx):
    """
    Try to find query_name in the target index. Returns (matched_name, method)
    or (None, None) if no safe match found.
    """
    c = canonical_name(query_name)
    if c in canon_idx:
        return canon_idx[c], 'exact_case_insensitive'

    r = expand_range_notation(query_name)
    if r and r in range_idx:
        return range_idx[r], 'range_notation'

    return None, None


if __name__ == "__main__":
    # Self-test with the known R1-6 / R1-R6 case and a few edge cases
    tests = [
        ('R1-6', ['R1-R6', 'R7', 'R8'], 'R1-R6'),
        ('R1-R6', ['R1-6', 'R7', 'R8'], 'R1-6'),
        ('r7', ['R7', 'R8'], 'R7'),
        ('Lai', ['lai', 'Lawf1'], 'lai'),
        ('Dm12', ['Dm1', 'Dm12', 'Dm16'], 'Dm12'),
        ('XYZ123', ['R1-6', 'R7'], None),
    ]
    print("Self-test:")
    for query, targets, expected in tests:
        canon_idx, range_idx = build_match_index(targets)
        matched, method = find_match(query, canon_idx, range_idx)
        status = "PASS" if (matched == expected or (matched is None and expected is None)) else "FAIL"
        print(f"  [{status}] '{query}' in {targets} -> matched='{matched}' method={method} (expected='{expected}')")
