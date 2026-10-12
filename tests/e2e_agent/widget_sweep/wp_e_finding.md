# WP-E finding — dead-end chains & conditional widgets (F5)

**Outcome: implemented as planned for all three legs, with two real
product fixes found on the way; one acceptance number lands differently
than the plan assumed (no-probe 11 = 8 benign + 3 conditional, not 8) —
see the accounting section.**

## Leg 1 — stringdb selected_row end-to-end (criterion 1: MET)

The WP-C finding's "gene count exceeds STRING_MAX_GENES" theory was only
part of the story. The full picture, established with live experiments:

1. The demo's sweep probe genes (X32.TYW5→A2RUC4 etc.) simply have **no
   STRING interactions** — the API returns header-only/empty for them at
   ANY score (curl-verified, including via plain R httr; my earlier
   "curl works" evidence in wp_c_finding.md was a header-only response I
   misread). Known genes (TP53/EGFR/...) return data normally.
2. With a STRING-connected selection (5 demo features:
   X1452.EGFR/X1477.MYC/X1570.TP53/X2033.BRAF/X2371.CCND1 →
   P00533/P01106/P04637/P15056/P24385; 374 enrichment rows via curl),
   the WHOLE chain works through the agent surface
   (`verify_wp_e.mjs`, 18/18):
   - `set_omics_viewer_state` features (5) → `run_control
     resultspace.stringdb.run` (epoch 1) → **50 selected_row choices
     materialize** (first: "GO:0009416"; the enrichment eventReactive is
     pulled lazily by the capability queries themselves — no tab needed)
   - the StringDB tab renders the enrichment table: **"Showing 1 to 10
     of 373 entries"**
   - `set_widgets selected_row=["GO:0009416"]` → applied, **store
     read-back "GO:0009416"**, and the table row is **actually selected
     (tr.selected = 1)** — DOM-level proof.
3. New sweep semantic case `tool:stringdb_chain_run_then_select` encodes
   exactly this chain (select → run → terms → set → DOM row check), with
   a two-attempt retry (the STRING endpoint is external and can
   transiently fail). Standalone gate run: **success in 6.2 s**; in the
   full sweep it needs the baseline entry to be success so a transient
   failure becomes auto-retriable.

## Leg 2 — batch *_selected_row (probe-verified; structurally conditional)

The rows are COMPUTED results: `phenotype_results`/`features_results`
req the show-toggle, a sample selection (two groups), and (R-H4) panel
visibility — a WP-A-style "make providers activation-independent" fix is
not applicable (the values do not exist before the comparison runs).
Treatment = precondition notes (capability records) + end-to-end probe:
samples selected (3) + show toggle on + panel visible → **table computes
("Showing 1 to 8 of 16 entries")** → **16 choices materialize** →
`set_widgets batch_phenotype_selected_row=["General|All|Dose"]` →
applied + **tr.selected = 1**.

## Leg 3 — survival_censor (criterion 2: MET; two real fixes)

The plan assumed "the module already handles late-render pushes". It did
not, quite — two genuine bugs found and fixed:

1. **`module_survival.R` push ordering (R-M7 contract broken)**: the
   epoch-push observer recorded `.sv_desired(cv)` AFTER the
   `reactive_checkpoint()` gate — a push landing before the survival
   view rendered was dropped and the slider rebuilt at the range max
   (observed live: written 20 → slider/store 990→106 paths). The
   recording (+clamp) now happens BEFORE the gate; the gate only guards
   the immediate `updateSliderInput`.
2. Probe-side (no product bug): the first regex matched a `.clientdata`
   output-visibility key instead of the slider input.

Evidence (probe, 18/18): cascade write to `Surv|all|OS` (3 keys applied)
→ slider input EXISTS (`...-censor`, initial 990 = range max) →
`set_widgets survival_censor=20` → applied, **slider 990 → 106 and store
106 in agreement** (20 clamps into the rendered time range 106..990 —
R-M7 clamping working as designed). The pre-existing testServer case
(test_agentWidgets S4: Surv column set → censor push → slider message)
remains green. The PRECONDITION note is in the binding help text.

Also fixed en passant: the result-space navbar materializes only after
the dataset settles, and the tab is labeled **"StringDB"** — probes must
scope to `#app-resultspace-analyst` and match case-insensitively.

## Acceptance arithmetic (criterion 3: DEVIATION, reported)

Baseline no-probe = **11**, decomposing exactly into:
- **8 benign** (F7, single-valid-value): feature_space x/y_subset,
  attr4.tooltip_subset, sample_space x/y_subset, ptm xax ×3
- **3 F5-conditional** (the WP-E legs): stringdb.selected_row,
  batch_phenotype_selected_row, batch_feature_selected_row

The plan's target "8 benign" assumed the three would become probeable in
the sweep after the WP-C run. Mechanically they cannot: the binding-case
loop runs BEFORE any action/semantic case, and these choices are
computed results that do not exist at binding-case time (they need a
successful run / visible panel / Surv variable). Each of the three is
now individually justified: PRECONDITION note in its capability record +
an end-to-end probe-verified chain (and, for stringdb, a sweep semantic
case). Their mechanical sweep statuses stay `no-probe` (a fresh sweep
cannot derive a probe for them), so the baseline keeps 11 — re-curating
them to a synthetic status would make every future run's no-probe count
exceed the baseline and fail the gate spuriously. **The plan's "8" is
therefore not reached mechanically; the honest number is 11 = 8 + 3
documented-conditional.** Flagged for Chen's accept/redirect.

## Guard board (before → after)

widgetStore 69→69, agentWidgets 104→104*, snapshotRoundTrip 63→63
(direct `Rscript`, exit 0, plan lines quoted in the WP report),
agentUiEffects 113 + Tier A 112/112 (post-change).
\* test_agentWidgets was ALREADY FAILING before WP-E (exit 10) on two
exact-key-inventory assertions that WP-D's deliberate search_value
binding had shifted (23→24, 26→27) — reported, and the two fixtures
updated with the WP-D key (the same deliberate-fixture class as the
WP-C allowlists).

## Artifacts

- `tests/e2e_agent/verify_wp_e.mjs` — the 18/18 probe (all three
  chains, store + DOM level)
- `tests/e2e_agent/widget_sweep/wp_e_results.json` — its durable output
- sweep: new semantic case `tool:stringdb_chain_run_then_select`
- product fixes: `R/module_survival.R` (push ordering),
  help-text PRECONDITION notes in `R/module_string.R`,
  `R/module_batch_comparison.R`, `R/module_survival.R`
