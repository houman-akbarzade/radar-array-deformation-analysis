# Public modeling report

- [`array_deformation_report.pdf`](array_deformation_report.pdf): English scientific report for the synthetic C-band portfolio model.
- [`array_deformation_report.tex`](array_deformation_report.tex): editable LaTeX source.
- Figures and underlying CSVs: [`../public_results/`](../public_results/).

To rebuild with a LaTeX distribution providing `report`, `amsmath`, `graphicx`, `subcaption`, `booktabs`, `hyperref`, `lmodern`, and `microtype`:

```bash
cd report
pdflatex -interaction=nonstopmode -halt-on-error array_deformation_report.tex
pdflatex -interaction=nonstopmode -halt-on-error array_deformation_report.tex
```

Run from `report/` because the source resolves figures through `../public_results/`. The report uses the public Python configuration and stored numerical summaries. The original engineering project's names, parameters, and measurements are not included; the public model has not been validated against hardware or a complete structural model. MATLAB source is provided separately, but was not executed in a MATLAB runtime for this release.
