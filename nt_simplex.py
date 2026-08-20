"""
NT probability simplex utilities.

Cell types are points on the 6-simplex of FAFB classifier outputs
(ACH, GABA, GLUT, DA, SER, OCT). Distance on this simplex is the right
geometry for "does this type look like a known confusion fingerprint?"
— which entropy (a scalar) cannot answer.
"""
from __future__ import annotations

import ast
import json
import re
from typing import Iterable

import numpy as np

# FAFB classifier output categories (histamine is not among them)
NT_ORDER = ("ACH", "GABA", "GLUT", "DA", "SER", "OCT")
FAST_TRANSMITTERS = ("ACH", "GABA", "GLUT")
MONOAMINES = ("DA", "SER", "OCT")

_NUMPY_INT = re.compile(r"np\.int(?:64|32)\((\d+)\)")


def parse_nt_distribution(raw) -> dict[str, int]:
    """Parse a count dict from analysis.py output (JSON or legacy repr)."""
    if raw is None:
        return {}
    if isinstance(raw, dict):
        return {str(k): int(v) for k, v in raw.items()}
    if isinstance(raw, float) and np.isnan(raw):
        return {}
    text = str(raw).strip()
    if not text:
        return {}
    text = _NUMPY_INT.sub(r"\1", text)
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        data = ast.literal_eval(text)
    if not isinstance(data, dict):
        raise ValueError(f"NT distribution is not a dict: {raw!r}")
    return {str(k): int(v) for k, v in data.items()}


def counts_to_vector(counts: dict[str, int], order: Iterable[str] = NT_ORDER) -> np.ndarray:
    """Integer counts aligned to `order` (unknown labels dropped)."""
    return np.array([int(counts.get(nt, 0)) for nt in order], dtype=np.float64)


def normalize(vec: np.ndarray) -> np.ndarray:
    total = float(vec.sum())
    if total <= 0:
        return np.zeros_like(vec, dtype=np.float64)
    return vec / total


def shannon_entropy_from_probs(probs: np.ndarray) -> float:
    p = np.asarray(probs, dtype=np.float64)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    return float(-np.sum(p * np.log2(p)))


def js_divergence(p: np.ndarray, q: np.ndarray, eps: float = 1e-12) -> float:
    """Jensen-Shannon divergence in bits (range 0..1)."""
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    p = np.clip(p, 0.0, None)
    q = np.clip(q, 0.0, None)
    if p.sum() <= 0 or q.sum() <= 0:
        return 1.0
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)

    def _kl(a, b) -> float:
        mask = a > eps
        return float(np.sum(a[mask] * np.log2((a[mask] + eps) / (b[mask] + eps))))

    return 0.5 * _kl(p, m) + 0.5 * _kl(q, m)


def js_distance(p: np.ndarray, q: np.ndarray) -> float:
    """Jensen-Shannon distance = sqrt(divergence)."""
    return float(np.sqrt(max(js_divergence(p, q), 0.0)))


def _kl_batch(A: np.ndarray, B: np.ndarray, eps: float) -> np.ndarray:
    """Row-wise sum_k A_k * log2(A_k / B_k), treating 0*log(0/x) as 0."""
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(A > eps, (A + eps) / (B + eps), 1.0)
        terms = np.where(A > eps, A * np.log2(ratio), 0.0)
    return terms.sum(axis=1)


def batch_js_divergence(P: np.ndarray, q: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """
    Vectorized Jensen-Shannon divergence (bits): each row of P (shape (B, K))
    against a single reference q (shape (K,)). Numerically identical to
    calling `js_divergence(row, q)` for each row -- see
    tests/test_core.py::TestNtSimplex::test_batch_matches_scalar.
    """
    P = np.clip(np.asarray(P, dtype=np.float64), 0.0, None)
    P = P / np.maximum(P.sum(axis=1, keepdims=True), eps)
    q = np.clip(np.asarray(q, dtype=np.float64), 0.0, None)
    q = q / max(float(q.sum()), eps)
    Q = np.broadcast_to(q, P.shape)
    M = 0.5 * (P + Q)
    return 0.5 * _kl_batch(P, M, eps) + 0.5 * _kl_batch(Q, M, eps)


def nearest_seed_js_batch(P: np.ndarray, seed_vectors: list[np.ndarray]) -> np.ndarray:
    """Per-row minimum JS divergence from P (shape (B, K)) to any seed vector."""
    if not seed_vectors:
        return np.full(len(P), np.nan)
    dists = np.stack([batch_js_divergence(P, seed) for seed in seed_vectors], axis=1)
    return dists.min(axis=1)


def mean_prototype(vectors: list[np.ndarray]) -> np.ndarray:
    stacked = np.vstack(vectors)
    return normalize(stacked.mean(axis=0))


def entropy_from_count_matrix(counts: np.ndarray) -> np.ndarray:
    """Shannon entropy in bits for each row of a (G, K) count matrix."""
    counts = np.asarray(counts, dtype=np.float64)
    n = counts.sum(axis=1, keepdims=True)
    n_safe = np.maximum(n, 1.0)
    p = counts / n_safe
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.zeros_like(p)
        mask = p > 0
        logp[mask] = np.log2(p[mask])
    ent = -np.sum(p * logp, axis=1)
    ent[n.ravel() <= 0] = 0.0
    return ent
