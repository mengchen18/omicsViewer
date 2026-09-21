# HANDOVER — agent accuracy & widget control plane

Written 2026-09-22 at the end of the S3 session. Read this first in a
fresh context, then `AGENT_ACCURACY_PLAN.md` (canonical plan + status) and
`AGENTS.md` (environment + commands). Delete or trim this file once absorbed.

## Where things stand

Branch `agent-driven-exploration`, working tree clean through the S3
commit. DESCRIPTION is 2.1.1. The assistant works end-to-end with a real
provider, and the **generic widget tier (S3)** is live: registry-driven
`list_widgets` / `get_widget` / `set_widgets` tools over the canonical
store, heatmap parameters registered on all three heatmap instances, .ESS
snapshots carrying `widget_store`.

Commit map (S3 session, on top of the S2 history in git log):

| Change | What |
|---|---|
| S3 | generic tier + heatmap registration + snapshot bridge + observer-GC fix |

Verified working (do not re-litigate without evidence):

- Unit suites: widgetStore 37, agentWidgets 31 (new), aiAssistantTools 16,
  triselectorCascade 6, scatterSelection 6, quickViews 17, appState 30,
  agentAssistant 35, agentFigures 22, agentLogging 21
- Tier A browser: **40/40 × 2 consecutive runs** (incl. the new S3 section:
  heatmap palette + margin via the `widgets` hook, per-key rejection
  feedback, unknown-id suggestions, tab restoration afterwards)
- Live Tier B (glm-5.3-flash), task "Switch to the Heatmap tab and change
  the heatmap color panel to RdGy": model chained get_state →
  `list_widgets(section="dataspace")` (26 widgets) → semantic
  `set_omics_viewer_state` for the TAB → `set_widgets` with
  `{"dataspace.expr_heatmap.heatmap_colors": "RdGy"}` → applied cleanly,
  zero rejections — exactly the tier contract. Provider `"null"` sentinels
  normalized at the boundary. Log archived under tests/e2e_agent/artifacts/.

## The two bugs found this session (both latent, both fixed)

1. **Child-store read**: `store_read(child, ids = NULL)` returned empty —
   `names(store$values)` was evaluated on the child env before resolving
   the root. Fixed; unit-tested.
2. **Observer GC** (the important one): Shiny holds observer dependencies
   weakly in the reactive graph. Store-glue observers that (a) call
   `store_epoch(store)()` inline (creating an ephemeral reactive each run)
   and (b) are not referenced anywhere are **garbage collected between
   flushes** — later store→UI pushes silently never fire. Reproduced
   deterministically by calling `gc()` between applies. This caused the
   sporadic Tier A stress/restore flakes attributed to timing before.
   Fix pattern (now in heatmap + meta_scatter glue): create the epoch
   reactive ONCE at module level and keep every store-glue observer in a
   module-level list (`.heatmap_keep(...)` / `.scatter_keep(...)`).
   **Rule for S4: every new store-glue observer must be retained.**

## Architecture rules to carry forward

- `reactive_selector1() %||% input$analysis` cascades: store-backed
  (meta_scatter) vs store-less (feature_general, fgsea, geneshot, tables,
  attr4) regimes both matter; `tests/test_triselectorCascade.R` guards both.
- Generic tier wire contract: ROOT store + full canonical ids; patch is a
  JSON-object string (schema can't express dynamic keys); ellmer
  tibble/array coercions normalized in `.agent_widget_normalize_patch`.
- Tool bodies wrap store access in `isolate(withReactiveDomain(...))` —
  choices providers read module reactives (triset()).
- Diff-only writes: snapshot restore through the store is additive and
  idempotent; panel-status restoration stays authoritative until S4.
- Registry/discovery show ONLY user-editable widgets (agent visibility ==
  user editability). `store_registry_view` prefix match is
  whole-component (`dataspace` ≠ `dataspace2`).

## Known open threads (none blocking S4)

1. Re-saving provider settings mid-stream kills the in-flight ellmer turn —
   cosmetic, backlog
2. WP4 log summarizer not built yet; logs are clean and archived per run
3. WP1 (`sections` param) + WP3 (figure spec round-trip) pending
4. Demo's truncated sample defaults (`PCA|All|PC1(`) still cleanly rejected
5. `store_apply` leaves `pending` entries for system-origin seeds whose
   pushed value equals the current input (never acked) — harmless today
   (pending is overwritten on the next write), revisit if ack bookkeeping
   ever becomes user-visible

## Next stage: S4 (per plan §6.4)

Migrate remaining modules onto the store and make snapshots fully
store-round-tripped:

1. Tables (dataTable modules) then remaining heatmap widgets (sorting,
   clustering, annotations), then result-space modules — each a small
   S2-style migration (register bindings, UI→store sync, store→UI push
   with observer retention, per-key-resilient restore)
2. Re-route snapshot save/restore entirely through `store_snapshot`/
   `store_restore` (the S3 bridge already embeds `widget_store` in .ESS;
   make it authoritative and retire the status-path duplication for
   migrated keys)
3. Acceptance: snapshot round-trip equals store state exactly; Tier A
   green across modules; each migration registers its widgets so the
   generic tier covers them automatically
4. Then WP1 (sections), WP3 (spec round-trip), WP4 (log summarizer);
   re-run the Tier B task set (tests/e2e_agent/tier_b_tasks.md — now
   includes the S3 generic-tier tasks 13–14)

## Session workflow reminders

- Reinstall after every source change
  (`Rscript -e "install.packages('<repo>', repos=NULL, type='source')"`)
- Kill orphaned R processes on ports 7775–7786 before browser runs
- NOT_CRAN=true needed for the browser-effect test
- provider.env is gitignored and holds the live key — never commit, never
  echo; artifacts under tests/e2e_agent/artifacts/ are gitignored
- Log timestamps are UTC (user is UTC+2)
