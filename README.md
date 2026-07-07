# LocusBlend

LocusBlend is a Streamlit app for visualizing locus-level summary statistics with one, two, or three index variants and blended linkage disequilibrium (LD) coloring.

## Installation

```bash
pip install -r requirements.txt
```

## Run Streamlit

```bash
streamlit run app.py
```

## Load The Example Files

The app loads the bundled synthetic examples by default:

- `examples/example_gwas_dataset.csv`
- `examples/example_qtl_dataset.csv`

To reload them manually, use the sidebar upload controls and select those two CSV files. Keep **LD source** set to **Precomputed LD columns**. The default locus center is `500000`, with index variants `rsDemo001`, `rsDemo002`, and `rsDemo003`.

## Required Input Columns

Summary-statistic CSV files must include:

- `CHR`
- `BP`
- `P`
- `A1`
- `A2`

Recommended columns:

- `rsid`
- `BETA`
- `SE`
- `A1FREQ`
- `N`

For the public demo and path-free use, include precomputed LD columns:

- `r2_1`
- `r2_2`
- `r2_3`

## Example Index Info

`examples/example_index_info.csv` uses this format:

```csv
BP,rs1,rs2,rs3
500000,rsDemo001,rsDemo002,rsDemo003
```

`BP` is the locus center position. `rs1` is the first index SNP, `rs2` is the second index SNP, and `rs3` is the third index SNP.

Two-index mode uses `rs1` and `rs2`. Three-index mode uses `rs1`, `rs2`, and `rs3`.

## LD References

LD reference files are not included. The demo uses precomputed LD columns in the example summary-statistic files. For user data, include `r2_1`, `r2_2`, and `r2_3` directly or generate LD externally and upload an LD matrix or long-format LD table.

## Citation

TBD

## License

TBD
