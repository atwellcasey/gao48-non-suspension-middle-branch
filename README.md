# A Non-Suspension Polynomial Solution in Gao's Five-Dimensional Middle Branch

**Author:** Casey Atwell  
**Manuscript version:** 1.3  
**Manuscript date:** 18 August 2026

This repository accompanies the v1.3 manuscript proving that there exists an exact polynomial
solution of Gao's five-dimensional middle-branch equation whose associated unpadded sweep is not
polynomially equivalent, under source-target polynomial automorphisms, to any triangular
one-coordinate suspension of a four-dimensional Gao sweep.

The theorem is stated independently of Gao's undefined phrase "equivalent to a suspension." Under
the precise source-target polynomial-equivalence and triangular one-coordinate suspension class
defined in the manuscript, the result gives the stated conditional consequence for Gao Open
Problem 4.8.

## Contents

- `paper/` - manuscript PDF and LaTeX source
- `verification/` - exact SymPy verifier and captured output
- `audits/` - source, build, and PDF preflight audits
- `REVISION_NOTES_v1.3.txt` - release-level exposition changes
- `CITATION.cff` - citation metadata
- `SHA256SUMS.txt` - repository integrity hashes

## Exact verification

Run:

```bash
python verification/verify_gao48_middle_branch_v1.3.py
```

Expected first line:

```text
PASS: exact witness verification
```

The verifier checks the explicit witness identities, Gao coefficients, middle-branch equation,
sweep determinant, a rank-three minor, the elimination relation, and displayed singular-gradient
containments. The normalization, UFD, critical-value, and non-cylinder arguments are theoretical
proofs in the manuscript rather than claims of computational certification.

## Scope

Primary theorem: the formally defined triangular one-coordinate Gao-suspension class is excluded
under source-target polynomial equivalence.

Gao Open Problem 4.8: the consequence remains conditional on interpreting Gao's undefined
suspension terminology by the precise class defined in the paper.

## License

A license must be added after author approval.
