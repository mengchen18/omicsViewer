# omicsViewer — Developer Notes

## Git Structure

- **`master`** — active development branch, currently v1.15.x (most recent)
- **`RELEASE_3_22`** — Bioconductor release branch, currently v1.14.x
- **`remotes/origin/devel`** — upstream Bioconductor devel
- **`remotes/upstream/RELEASE_3_*`** — upstream Bioconductor release snapshots

Bug fixes should generally be applied to **`master`** first. If the fix is critical and release-relevant, also apply to `RELEASE_3_22`.

There is a saved stash (`stash@{0}`) on `RELEASE_3_22` from the initial scatter plot bug fix session — drop with `git stash drop` if no longer needed.

## Install Package Locally

```bash
Rscript -e "install.packages('/media/chen/Elements/github/omicsViewer', repos = NULL, type = 'source')"
```

- `devtools` is **not** installed on this machine — use `install.packages(..., type = 'source')` instead
- `pkgload` is available if you need `load_all()` in an R session: `pkgload::load_all('.')`

## Start the Shiny App

```bash
Rscript -e "
  options(shiny.port = 7775, shiny.host = '0.0.0.0')
  omicsViewer::omicsViewer('/media/chen/Elements/github/omicsViewer/inst/extdata/')
" > /tmp/omicsviewer_shiny.log 2>&1 &
```

- Default port: **7775** (chosen to avoid conflict with RStudio Server on 8787)
- Access at: `http://localhost:7775`
- Demo data: `inst/extdata/demo.RDS` (2702 features, 60 samples)
- Log file: `/tmp/omicsviewer_shiny.log`
- Kill the app: `kill $(ss -tlnp | grep 7775 | grep -oP 'pid=\K[0-9]+')`

## Debugging

**Check if the app is running:**
```bash
ss -tlnp | grep 7775
cat /tmp/omicsviewer_shiny.log
```

**Common startup warning** — `Warning: stack imbalance in '::', 2 then 4` — harmless, caused by `omicsViewer::omicsViewer()` double-colon at the top level.

**WebGL:** The remote desktop environment (accessed via gnome-remote-desktop / RDP on ports 3389/3390) has no GPU/WebGL support. The app handles this via browser-side detection in `plotly_scatter_ui()` — the JS sets `input$webgl_supported` and `renderPlotly` conditionally calls `toWebGL()`.

**Chrome browser automation:** Claude Code is running with `--chrome` flag (PID visible via `ps aux | grep "claude --chrome"`). Chrome native messaging host is at `~/.claude/chrome/chrome-native-host` → calls `claude --chrome-native-host`. Extension ID: `fcoeoabgfenejglbffodgkkbkcdhcgfn` (v1.0.63).

## Key Source Files

| File | Role |
|------|------|
| `R/module_scatter.R` | Scatter/beeswarm plot module — includes regression line logic and WebGL fallback |
| `R/00_export_runExpressionSetViewer.R` | `omicsViewer()` entry point |
| `R/00_export_prepEsetViewer.R` | `prepOmicsViewer()` data preparation pipeline |
| `R/L0_module_app.R` | Top-level Shiny app (file loading, snapshots, Excel export) |
| `R/L1_module_data_space.R` | Left panel (scatter, heatmap, tables) |
| `R/L1_module_result_space.R` | Right panel (enrichment, STRING, survival, etc.) |
| `R/auxi_stats.R` | `multi.t.test()`, `correlationAnalysis()` |
| `R/auxi_databaseFormat.R` | SQLite read/write (`saveOmicsViewerDb`, `readESVObj`) |

## Data Format Notes

- Input: `ExpressionSet` or `SummarizedExperiment`
- Column naming convention: `Category|Subcategory|Variable` (e.g. `ttest|KO_vs_WT|log.fdr`)
- Gene sets stored as 3-column data.frame (`featureId`, `gsId`, `weight`) or sparse `dgCMatrix`
- Storage formats: `.RDS` or `.sqlite`/`.db` (SQLite)

## Bugs Fixed (this session)

In `R/module_scatter.R` — regression line (`regressionLine = TRUE`) block:

1. **`is.finite()` filter** — replaced `!is.na()` with `is.finite()` so that `Inf`/`-Inf` are also excluded before `lm()` and `cor.test()`
2. **`cor.test()` on cleaned data** — changed `cor.test(x, y)` to `cor.test(df$x, df$y)` to use the filtered, sorted data rather than original input vectors
3. **WebGL browser fallback** — moved `toWebGL()` out of `plotly_scatter()` into `renderPlotly`; JS in `plotly_scatter_ui()` detects browser WebGL support and sets `input$webgl_supported`; standard plotly is used as fallback
