# HANDOVER — agent accuracy & widget control plane

Written 2026-09-22 after the S4 heatmap-completion session. Read this
first in a fresh context, then `AGENT_ACCURACY_PLAN.md` (canonical plan +
status table) and `AGENTS.md` (environment + commands). Delete or trim
this file once absorbed.

## Where things stand

Branch `agent-driven-exploration`, DESCRIPTION 2.1.1. S3 (generic widget
tier) shipped previously; **S4 steps 1–2 are done**: heatmap completion
(multi_select kind + full 13-widget parameter panel) and the data-space
table migration (multi-selection switch + shown-columns set on all three
tables).

Commit map (newest first, on top of the S3 history in git log):

| Change | What |
|---|---|
| S4-2 | dataTable store migration (multi_selection + columns, min=1 bound) |
| S4-1 | multi_select kind + heatmap sorting/clustering/annotation registration + tests |

Verified this session (do not re-litigate without evidence):

- Unit: widgetStore 52, agentWidgets 51, aiAssistantTools 19, appState 30,
  agentAssistant 38, agentFigures 26, agentLogging 21, triselectorCascade 6,
  quickViews 17, scatterSelection 6, shinyAuxi 10, tableWidgetState 5
- Tier A browser: **55/55 × 2 consecutive runs** (4c heatmap section:
  five-key sorting/clustering/annotation apply incl. JSON-array
  multi_select, wholesale selection replacement, per-key rejection with
  "Closest matches: General|All|Cell.line" quality suggestions; 4d tables
  section: multi_selection + columns apply with visible header changes,
  min-count rejection, cross-table apply)

## What S4-1 added (carry-forward rules)

- **`multi_select` kind** (auxi_widgetStore.R): value = character vector;
  empty vector / empty JSON array / `""` all clear the selection; entries
  validated against choices_provider/values; **min/max bound the entry
  count** (tables use min=1 so an empty column set is rejected instead of
  blanking the table); literal sentinel strings (`"[]"`, `"null"`) remain
  omitted optionals per the shared AGENT_SENTINEL_STRINGS invariant.
- **Heatmap full registration** (heatmapshinyApp.R): 13 keys per instance
  under `dataspace.{cor,expr,dyn}_heatmap.*`. UI→store sync (clearing a
  multi-select syncs `character(0)`), seeding, and store→UI push follow
  the proven S3 shape; server-side selectize ids (`rowSortBy`,
  `annotRow`) push through `updateSelectizeInput`, everything else
  `updateSelectInput`/`updateSliderInput`.
- **New invariant — choices providers must never `req()`**: a
  `shiny.validation` condition raised inside a choices_provider aborts
  the whole `store_apply` transaction (`tryCatch(error=)` cannot catch
  it). Heatmap providers read module reactives through defensive no-req
  helpers (`.heatmap_hcl_names`, `.heatmap_fd_cols`, …) that mirror the
  UI choices exactly. Apply this pattern to every future registration
  whose choices come from module reactives.
- **Observer retention** still mandatory: every store-glue observer goes
  through `.heatmap_keep(...)` (or the module's equivalent).
- **Snapshot status-path duplication intentionally kept for heatmaps**:
  `observeEvent(status())` still restores heatmap keys directly (legacy
  pre-S3 snapshots need it; seeding-vs-status ordering makes
  store-held-key skipping unsafe). Both paths write identical values,
  idempotently. The full re-route (store authoritative, status path
  retired for migrated keys) stays the LAST S4 step, after every module
  is migrated — see plan §6.4.

## Next steps (S4 continuation, plan §6.4 order)

1. **Result-space modules** (fgsea, ora, survival, string, geneshot,
   feature/sample_general, attr4, PTMotif, …) — the remaining migrations;
   each is a small S2-style move: register bindings (no-req choices
   providers!), UI→store sync, seeding, push with observer retention,
   per-key-resilient restore. The dataTableDownload instances embedded in
   those modules carry no state (buttons only) and need nothing.
3. Re-route snapshot save/restore fully through the store; retire
   status-path duplication for migrated keys.
4. Then WP1 (`sections` param on get_omics_viewer_state), WP3 (figure
   spec round-trip), WP4 (log summarizer); re-run the Tier B task set
   (tests/e2e_agent/tier_b_tasks.md) with heatmap-annotation and
   table-column tasks added.

## Known open threads (none blocking)

1. Re-saving provider settings mid-stream kills the in-flight ellmer
   turn — cosmetic, backlog
2. WP4 log summarizer not built yet; logs are clean and archived per run
3. WP1 + WP3 pending
4. Demo's truncated sample defaults (`PCA|All|PC1(`) still cleanly rejected
5. `store_apply` leaves `pending` entries for system-origin seeds whose
   pushed value equals the current input (never acked) — harmless today,
   revisit if ack bookkeeping ever becomes user-visible. NOTE: multi_select
   seeds (empty annotation selections) hit this regularly since pushing
   `character(0)` onto an already-empty select fires no input change.

## Session workflow reminders

- Reinstall after every source change
  (`Rscript -e "install.packages('<repo>', repos=NULL, type='source')"`)
- Kill orphaned R processes on ports 7775–7786 before browser runs
- NOT_CRAN=true needed for the browser-effect test
- provider.env is gitignored and holds the live key — never commit, never
  echo; artifacts under tests/e2e_agent/artifacts/ are gitignored
- Log timestamps are UTC (user is UTC+2)
