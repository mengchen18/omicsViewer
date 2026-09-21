# Tier B benchmark task set (WP5 + S3 additions)

Scripted end-to-end tasks against `demo.RDS`, executed with the Tier B
driver:

```
node tier_b.mjs "task prompt"
```

Protocol (plan section 7, decision Q9): fresh chat per task, diagnostic
logging on (`OMICSVIEWER_LLM_LOG=true`), one pinned provider/model for
baseline and post-change runs, serial execution. Smoke subset per change:
tasks 1, 4, 6, 8, 10 plus the S3 task; full set x3 at phase gates.
Majority-of-3 + triage on flakiness (Q8): identical repeated failures are
harness bugs, varying failures are model nondeterminism.

## Core tasks (WP5)

1. "What dataset is loaded?" (state orientation)
2. "How many genes are selected?"
3. "Switch to the Sample tab"
4. "Make a volcano plot of the current results" (scatter quick view)
5. "Plot log FDR vs mean difference with custom axes"
6. "Find genes matching 'kinase' and select the first five"
7. "Summarize the sample group annotation"
8. "Boxplot of expression for the selected genes grouped by sample group"
9. "Create a histogram of the t-test log FDR"
10. "Change the last figure's title to X" (revision round-trip)
11. "Color the volcano points by category" (revision)
12. "What can you change in the app?" (capability disclosure)

## S3 generic-widget-tier tasks

These exercise controls that no semantic (tier-1) tool covers, so the
model must discover them via `list_widgets`/`get_widget` and apply them
via `set_widgets`. Acceptance (plan section 6.4): one Tier B run each.

13. "Switch to the Heatmap tab and change the heatmap color panel to RdGy"
    (expects: `list_widgets` -> `set_widgets` with
    `dataspace.expr_heatmap.heatmap_colors: "RdGy"`; visible select update)
14. "Make the expression heatmap bottom margin larger" (numeric widget:
    `dataspace.expr_heatmap.margin_bottom`; assert slider/UI change)

Pass criteria: the visible widget ends in the requested state, no
unrelated widget changed, and rejected-key feedback (if any) leads to a
corrected retry rather than a fabricated success claim.
