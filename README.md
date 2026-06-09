# Are Choice Experiments Incentive Compatible? — A Replication

A reproducible **R / Quarto** replication of

> Lusk, J. L., & Schroeder, T. C. (2004). *Are Choice Experiments Incentive
> Compatible? A Test with Quality Differentiated Beef Steaks.* American Journal
> of Agricultural Economics, **86**(2), 467–482.

Replication project for **AGEC 5045 — Applied Discrete Choice Analysis (NTU, 114-2)**.
The original analysis was done in LIMDEP; this reproduces it in R.

## 📖 Read the book

**<https://JunJie-Chang.github.io/ADCA-replication/>**

## What's reproduced

- Data preparation (`harrison.xls` → choice-format panel)
- Tables 1–2 (demographics, choice frequencies)
- Table 3 (multinomial logit + Swait–Louviere scale test)
- IIA test via the universal "mother" logit
- Table 4 (HEV, RPL; MNP discussed)
- Table 5 (willingness to pay + Poe combinatorial test) and Figure 2 (HEV market shares)
- A chapter documenting and interpreting discrepancies vs. the published numbers

## Repository layout

```
_quarto.yml          Quarto book config (renders to docs/)
_common.R            shared setup: packages, helpers, data prep, load results
index.qmd            preface
01–09 *.qmd          chapters
R/                   numbered estimation scripts (01_data_prep … 06_iia)
harrison.xls         raw experimental data
output/              pre-estimated model objects (*.rds) and CSV tables
docs/                rendered website (served by GitHub Pages)
```

## Rebuild

```bash
# 1. (Optional) re-estimate everything, in RStudio, in order:
#    R/01 → R/02 → R/03 → R/06 → R/04 → R/05   (writes output/*.rds)
# 2. Render the site:
quarto render
```

The book loads `output/*.rds` rather than re-running the slow models, so it
renders in seconds.

## License / data note

Code is provided for educational use. The dataset originates from the study's
authors and is included here solely for course replication.
