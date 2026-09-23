# HANDOVER — agent accuracy & widget control plane

Written 2026-09-25 at the close of the WP8–WP12 session. Read this first
in a fresh context, then `AGENT_ACCURACY_PLAN.md` (per-WP status rows are
the source of truth) and `AGENTS.md` (environment, commands, quirks).
Delete or trim this file once absorbed.

## Where we are

Branch `agent-driven-exploration`, DESCRIPTION 2.2.0. **The whole
AGENT_ACCURACY_PLAN is complete**: WP0–WP9, WP11, and WP12 landed;
WP7 patch-mode gated NO (2026-09-24); WP10 resolved by construction
(store epoch). This session:

| What landed | Evidence |
|---|---|
| **WP8 semantic tools** | `set_enrichment_parameters` (ORA/fGSEA ranking/collapse + optional pathway row; opens the panel) and `set_table_view` (feature/sample/expression tables: columns, multi-selection, per-column filters, page; opens the tab). Thin validate + `store_apply` wrappers; per-key resilience with suggestions; L0 applies `apply_agent_enrichment` / `apply_agent_table_view`; test-hook ops `enrichment` / `tableview` |
| **Data-table DT state on the store** | `dataspace.tab_*.{page,column_filters}` — new `mapping` widget kind (named character vector; empty clears; keys validated against colnames), pushed via the DT proxy (`selectPage`/`updateSearch`, filters BEFORE page so paging targets the filtered set), acknowledged by the DT state report. Column ORDERING stays on tab_status (no DT proxy API). **En route fix: the module's dormant `tabproxy` used `ns("table")` — DT double-prefixes, every proxy message targeted a nonexistent table** (dead code until WP8; `dataTableProxy("table")` is correct) |
| **WP9 capability registry** | `R/auxi_agentCapabilities.R` — records GENERATED from store bindings + curated tool metadata (the shared help-text structure); `search_ui_capabilities` / `get_ui_capability` tools; overview carries capability COUNTS only. Widget ids map to covering semantic tools. En route: `.agent_suggest` now also scores whole-string edit distance (full-column typos got no suggestions before) |
| **WP11 conversation-in-snapshot** | Opt-in checkbox in the snapshot modal; assistant module returns an API (`snapshot_payload`/`restore_history`/`has_conversation`); payload = slim turns (base64 previews stripped, credential-like strings redacted, tool VALUES kept as inert context) + text transcript + figure registry (≤20; update_figure re-validates specs at use). Byte budget by real serialization (`object.size` over-counts S7 objects ~170 KB per EMPTY turn — use `serialize()` length). **Deliberate deviation: shinychat `history=TRUE` stores NOT enabled** (file-based production storage would violate the opt-in guardrail); .ESS is the single persistence path |
| **WP12 budgets** | `OMICSVIEWER_LLM_MAX_COST_USD` / `OMICSVIEWER_LLM_MAX_TOKENS` cumulative ceilings checked in `on_request_start`; usage accumulated in `on_request_end` (turn tokens = c(input, output, …), sum all non-NA); `budget_limit` joined the log taxonomy (9 classes). Deputy evaluation stays log-gated |

## Validation snapshot (all green at HEAD)

- Unit: agentCapabilities **45** (new), agentHistory **26** (new),
  aiAssistantTools **43**, agentAssistant **73**, agentWidgets 105,
  agentFigures 67, agentLogSummary 33, agentLogging 21, appState 40,
  widgetStore 52, quickViews 17, scatterSelection 6, triselectorCascade 6,
  tableWidgetState 10, shinyAuxi 10, stats 7, ora 23
- Tier A browser: **109/109 ×2** (new 4h: capability search/get, the
  enrichment tool + typo suggestions, table columns/filters/page through
  the real DT UI); agentUiEffects 110/110
- Live Tier B: none this session (tool surface is unit+Tier-A covered; the
  next full ×3 gate run should fold a set_enrichment_parameters /
  set_table_view task into tier_b_tasks.md)

## Flakes / rules discovered (cumulative; do not re-learn)

- **DT proxies**: `dataTableProxy(session$ns(...))` double-prefixes inside
  modules — pass the module-local id. Proxy messages defer until flush end
  (after output re-renders), so columns+filters pushes in one transaction
  land on the NEW table; `selectPage` beyond the filtered page count logs a
  client-side "Selected page is out of range" (harmless; the sync observer
  mirrors the clamped reality back into the store — override contract).
- **Mapping-kind sentinel semantics**: empty JSON `{}` / `"null"` are
  treated as OMITTED (AGENT_SENTINEL_STRINGS convention); an explicit
  empty map (length-0 named character) clears all filters — the tool
  exposes `clear_filters=true` for the model.
- **jsonlite auto_unbox** serializes length-1 vectors (receipt
  `applied`/`unchanged`) as SCALARS — JS test assertions must wrap with an
  `asArr()` normalizer before `.includes()`.
- **object.size on S7 objects** counts class metadata per instance
  (~170 KB for an empty UserTurn) — byte budgets must measure
  `serialize(x, NULL)` length instead.
- **ellmer 0.5.0**: `AssistantTurn(contents=, tokens=c(in,out,…), cost=,
  duration=, finish_reason=)`; `ContentToolResult(value=, error=)` drops
  `extra` (display payloads) — that is the WP11 slimming. Turns are
  saveRDS-safe.
- Settle on LOG truth for Tier B turns; request cap is load-bearing;
  provider stream hangs happen (infra, not model); selection state has
  multiple mirrors — patch them all; `store_read` is a plain snapshot
  (reactive consumers must read through `store_watch`); seeding from live
  inputs needs `mark_pending = FALSE`; UI→store sync mirrors only
  coherent triples; Tier A stress flake / MockShinySession quirks /
  observer GC / WebGL: AGENTS.md (unchanged).

## Session workflow reminders (unchanged)

Reinstall after every source change; roxygen after touching roxygen blocks
(regenerated man/ for capabilities/assistant/aiAssistant/testHooks/
widgetStore this session); kill orphaned R processes on 7775–7795 before
browser runs; Tier A ×2 convention; Tier B smoke subsets per change, full
×3 at gates; log timestamps are UTC (user is UTC+2).

## Next session candidates

The plan is complete; remaining ideas are log-gated or new work:
- Fold WP8 tool tasks into the Tier B benchmark (enrichment + table-view
  prompts) and re-run the full ×3 gate.
- WP7 patch-mode: reopen ONLY if a future full ×3 run shows revision
  failures persisting (merge base = WP3 registry spec).
- deputy::Agent migration evaluation when it stabilizes (WP12 step 2).
- STATEFUL_UI_GAPS.md still tracks the URL-restore links (not started) and
  SQLite metadata parity (intentionally deferred).
