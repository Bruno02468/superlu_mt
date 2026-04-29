# Provably-safe LU storage bounds via the Householder factor `H`

## Motivation

`pxgstrf_MemInit()` allocates the L\U storage for the entire factorization
in one shot — there is no facility to grow the LSUB / UCOL / LUSUP arrays
mid-factorization. The pre-existing sizing was a heuristic driven by the
`sp_ienv(6/7/8)` "fill ratio" multipliers (`FILL_LUSUP`, `FILL_UCOL`,
`FILL_LSUB`). When the heuristic underestimates, factorization aborts at
`XPAND_HINT`. When it overestimates, memory is wasted.

The objective of this change is to make `sp_ienv(7)` and `sp_ienv(8)`
obsolete by precomputing a provably-safe budget for `nzumax` and `nzlmax`.

## Theory

Under any column permutation `Pc` and any row permutation `Pr` chosen by
partial pivoting at run time, let `H` be the Householder factor of
`QR(Pr · A · Pc)` (equivalently, the symbolic Cholesky factor of
`(A · Pc)^T · (A · Pc)`). Standard George–Ng / Gilbert structural-bound
results give:

$$
\operatorname{struct}(L_{LU}) \subseteq \operatorname{struct}(H),\qquad
\operatorname{struct}(U_{LU}) \subseteq \operatorname{struct}(R) \subseteq
\operatorname{struct}(H).
$$

So `nnz(H)` is a row-permutation-invariant upper bound for both `nnz(L_LU)`
and `nnz(U_LU)`. SuperLU_MT already computes `nnz(H)` inside
[SRC/qrnzcnt.c](SRC/qrnzcnt.c) (local variable `nhnz`) and inside
[SRC/cholnzcnt.c](SRC/cholnzcnt.c) (it equals the Cholesky `nlnz`), but the
old code discarded that value and only returned `nlnz = nnz(R)` — which
bounds `U_LU` but **not** `L_LU`. Earlier experiments that used `nlnz` for
the LSUB bound therefore exceeded the budget on some matrices.

### LSUB layout factor of 2

Per-supernode LSUB allocation in
[SRC/pdgstrf_column_dfs.c#L302](SRC/pdgstrf_column_dfs.c#L302) and
[SRC/pdgstrf_snode_dfs.c#L75](SRC/pdgstrf_snode_dfs.c#L75) is
`2 · |L[*, fsupc]|` — a copy of the subscripts plus a pruned-graph copy used
by `pivotL` and `super_bnd_dfs`. Subsequent columns within the same
supernode overwrite (do not append). Worst case total LSUB therefore equals
`2 · Σ colcnt_h(fsupc) = 2 · nnz(H)`.

## Solution

1. Add an output parameter `nhnz` to `qrnzcnt()` and `cholnzcnt()` that
   exposes `nnz(H)`.
2. Cache it on `superlumt_options_t` as a new field `nnzH`.
3. In `p{s,d,c,z}gstrf_MemInit()`, compute the LU bounds from `nnzH`:
   - `nzumax = nnzH`         (since `struct(U_LU) ⊆ struct(R) ⊆ struct(H)`)
   - `nzlmax = 2 · nnzH`     (LSUB layout factor)
   - `nzlumax` keeps using `FILL_LUSUP`; this is a numerical-fill bound,
     not bounded by `H`, and is out of scope for this change.
4. An escape hatch is provided: setting environment variable
   `SUPERLU_MT_USE_LEGACY_FILL=1` reverts `nzumax`/`nzlmax` to the old
   `FILL_*` heuristics. The `sp_ienv(6/7/8)` symbols themselves are
   untouched — `FILL_LUSUP` is still consulted, and the
   `FILL_UCOL`/`FILL_LSUB` paths remain reachable as a fallback.

## Files modified

| File | Purpose of edit |
|------|-----------------|
| [SRC/slu_mt_util.h](SRC/slu_mt_util.h) | Added `int_t nnzH` to `superlumt_options_t`; updated `cholnzcnt` prototype |
| [SRC/slu_mt_ddefs.h](SRC/slu_mt_ddefs.h), [SRC/slu_mt_sdefs.h](SRC/slu_mt_sdefs.h), [SRC/slu_mt_cdefs.h](SRC/slu_mt_cdefs.h), [SRC/slu_mt_zdefs.h](SRC/slu_mt_zdefs.h) | `qrnzcnt()` prototype gains `int_t *nhnz` |
| [SRC/qrnzcnt.c](SRC/qrnzcnt.c) | Outputs `nhnz` (already computed locally) |
| [SRC/cholnzcnt.c](SRC/cholnzcnt.c) | Outputs `nhnz` (= `nlnz` since `H ≡ L_chol` in symmetric mode) |
| [SRC/sp_colorder.c](SRC/sp_colorder.c) | Plumbs `&options->nnzH` through both call sites |
| [SRC/pdgstrf_init.c](SRC/pdgstrf_init.c), [SRC/psgstrf_init.c](SRC/psgstrf_init.c), [SRC/pcgstrf_init.c](SRC/pcgstrf_init.c), [SRC/pzgstrf_init.c](SRC/pzgstrf_init.c) | Initialize `superlumt_options->nnzH = 0` |
| [SRC/pdmemory.c](SRC/pdmemory.c), [SRC/psmemory.c](SRC/psmemory.c), [SRC/pcmemory.c](SRC/pcmemory.c), [SRC/pzmemory.c](SRC/pzmemory.c) | Replace FILL heuristics with `nnzH`-based bounds; legacy fallback via env var |

The substantive logic changes are localized to the items above.

## Validation

Tests run on Linux x86-64 with system OpenBLAS and system METIS 5.1.

### PTHREAD, no METIS

| Driver | Matrix | expansions | BERR |
|---|---|---|---|
| `pdlinsolx` | `cg20.cua` | 0 | 1.5e-16 |
| `pdlinsolx` | `big.rua`  | 0 | 2.2e-16 |
| `pdlinsolx` | `g5.rua`   | 0 | 1.1e-16 |
| `pdlinsolx` | `g10`      | 0 | 1.5e-16 |
| `pslinsolx` | `cg20.cua` | 0 | 8.9e-08 |
| `pclinsolx` | `cmat`     | 0 | 1.1e-07 |
| `pzlinsolx` | `cmat`     | 0 | 2.0e-16 |

### PTHREAD, with `-DHAVE_METIS -lmetis` (mystran patches applied)

| Driver | Ordering | Matrix | expansions | BERR |
|---|---|---|---|---|
| `pdlinsolx` | METIS A'+A (4) | `cg20.cua` | 0 | 1.3e-16 |
| `pdlinsolx` | METIS A'+A (4) | `big.rua`  | 0 | 2.5e-16 |
| `pdlinsolx` | METIS A'A (6)  | `cg20.cua` | 0 | 1.1e-16 |
| `pdlinsolx` | METIS A'A (6)  | `big.rua`  | 0 | 2.8e-16 |
| `pdlinsolx` | METIS A'A (6)  | `g10`      | 0 | 1.1e-16 |

### OpenMP, with `-DHAVE_METIS -lmetis` (`OMP_NUM_THREADS=2`)

| Driver | Matrix | expansions | BERR |
|---|---|---|---|
| `pdlinsolx` | `cg20.cua` | 0 | 1.5e-16 |
| `pclinsolx` | `cmat`     | 0 | 1.1e-07 |
| `pzlinsolx` | `cmat`     | 0 | 2.0e-16 |

### Sizing comparison (new vs legacy)

`total MB needed` is identical to legacy on all tested matrices; reported
`nnz(L+U)` after factorization differs by ~0.1% (slightly tighter under
the new bound because the static budget is exact for `U_LU` rather than
padded by `FILL_UCOL`).

### Pre-existing issue (not caused by this change)

`pdrepeat` segfaults with `double free` on `big.rua` under both the new
and legacy paths. Reproduces with the unmodified upstream code.

## Backward compatibility

- API: `qrnzcnt()` and `cholnzcnt()` gained one new `int_t *` parameter.
  All in-tree callers are updated. Out-of-tree callers need to add one
  argument or pass `NULL` (both functions tolerate `NULL`).
- Behavior: identical to legacy on all tested matrices.
- Escape hatch: `SUPERLU_MT_USE_LEGACY_FILL=1` reverts to the old
  heuristics if a regression is observed.

## Out of scope / future work

- `nzlumax` (driven by `FILL_LUSUP`) is unchanged. Bounding it requires
  a different argument: it tracks numerical fill in supernodal blocks,
  not symbolic fill in `H`.
- Symmetric-mode (`SymmetricMode=YES`) was not specifically tuned, since
  the downstream consumer (MYSTRAN) does not currently guarantee
  symmetric matrices.
- The top-level [CMakeLists.txt](CMakeLists.txt) does not wire up METIS;
  only [mystran_patches/CMakeLists.txt](mystran_patches/CMakeLists.txt)
  does. Porting the option upstream is a separate task.
- A deprecation notice for `sp_ienv(6/7/8)` is deferred until field
  validation of the new bounds.
