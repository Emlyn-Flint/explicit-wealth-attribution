# Explicit Wealth Attribution

Companion code for the August 2026 Peresec Quant Nugget:

[**Multi-Period Attribution: Who Needs Linking Functions?**](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7345225)
Emlyn Flint, Florence Chikurunhe and Luke van Schaik, 2026.

The scripts reproduce the numerical examples in the note and compare:

- Explicit Wealth Attribution (EWA) --> new addition
- Frongello (2002)
- Berg (2014)
- Cariño (1999)
- Reverse GRAP (1997)

EWA preserves the ordinary single-period Brinson-Fachler Allocation, Selection and Interaction effects in money terms and reports benchmark growth on prior active wealth separately as Active-Wealth Carry.

## Files

- `ewa_companion.py` — Python implementation; standard library only.
- `ewa_companion.R` — R implementation; base R only.

Both scripts contain the same four-stock, two-sector, four-period example used in the Quant Nugget.

## Running the examples

### Python

Set:

```python
CASE = "base"
```

or:

```python
CASE = "pure_selection"
```

then run:

```bash
python ewa_companion.py
```

### R

Set:

```r
CASE <- "base"
```

or:

```r
CASE <- "pure_selection"
```

then run:

```r
source("ewa_companion.R")
```

The `base` case reproduces the main numerical example. The `pure_selection` case reproduces the controlled example in which Allocation and Interaction are zero.

## Reference

Flint, E., Chikurunhe, F. and van Schaik, L. (2026),  
*Multi-Period Attribution: Who Needs Linking Functions?*, Peresec Quant Nugget.

[Read the Quant Nugget](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7345225)

## Scope

These scripts are provided to reproduce the calculations in the accompanying research note. They are illustrative research code rather than a general-purpose performance-attribution library.
