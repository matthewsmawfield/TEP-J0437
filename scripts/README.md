# TEP-J0437 Analysis Scripts

Reproducible analysis pipeline for TEP-J0437 (PSR J0437-4715 scintillation holonomy).

## Layout

```
scripts/
├── run_pipeline.py              # Orchestrates all steps (canonical entry point)
├── sync_manuscript_headlines.py # Sync site/components numerics from results/*.json
├── audit_manuscript_consistency.py # Fail if stale headline literals reappear
├── steps/                       # Numbered pipeline steps → results/step_*.json
└── utils/
    ├── headline_stats.py        # Canonical headline values from frozen JSON
    ├── data_loader.py           # Cached J0437 closure loaders (prefers *_j0437.json)
    └── ...
```

## Run the pipeline

```bash
PYTHONPATH=. python scripts/run_pipeline.py
```

## After a pipeline run (manuscript maintenance)

```bash
PYTHONPATH=. python scripts/sync_manuscript_headlines.py
PYTHONPATH=. python scripts/audit_manuscript_consistency.py
npm run build:markdown
```

Canonical headline statistics are read from:

- `results/step_003_closure_final_summary_j0437.json` — weighted epoch ψ, |H|, sample counts
- `results/step_048_cmb_dipole_frame_analysis.json` — unweighted triplet-phase frame tests
- `results/step_007_independent_statistical_validation_results.json` — CV tables and null tests

The generic `step_003_closure_final_summary.json` is a legacy alias written alongside the J0437-specific file; downstream code prefers `*_j0437.json` when present.

## Manuscript source

Edit `site/components/*.html` only (not the generated `16-TEP-J0437-v0.3-Sintra.md`). Rebuild with `npm run build:markdown`.
