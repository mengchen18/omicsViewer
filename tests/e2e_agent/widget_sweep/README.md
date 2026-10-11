# Exhaustive widget sweep — every user-facing interactive widget, LLM-tool-pathway, UI-side validation

## Purpose

Establish, with zero gaps, which interactive widgets of the omicsViewer UI an
LLM assistant can operate end-to-end, and record per-widget evidence
(success / rejection / error / no-effect / timeout). No product code is
changed; this directory contains only the harness, per-widget case files,
per-widget result files and inventory artifacts, for later analysis.

## How "the LLM drives the widget"

The app is launched with `OMICSVIEWER_TEST_HOOKS=true`. The hidden hook
module (`R/module_agentTestHooks.R`) exposes the exact server-side pathways
of the assistant tools — same callbacks, same validation, same receipts:

| hook op | assistant tool | server pathway |
|---|---|---|
| `widgets` | `set_widgets` | `agent_widget_apply()` → `store_apply(origin="agent")` |
| `scatter` | `set_scatter_view` | `apply_agent_scatter_view()` |
| `state` | `set_omics_viewer_state` | `apply_agent_state()` |
| `enrichment` | `set_enrichment_parameters` | `apply_agent_enrichment()` |
| `tableview` | `set_table_view` | `apply_agent_table_view()` |
| `overview` | `get_omics_viewer_state` | `agent_state()` |
| `capabilities` | `find_controls` / `describe_control` | `agent_capability_search/get()` |

Every case therefore executes precisely what a model tool call executes.

## Enumeration (no missing, no gap)

1. **DOM-authoritative**: every Shiny input binding (`.shiny-bound-input`)
   enumerated after visiting EVERY tab (both navbars + nested tabsets),
   opening EVERY dropdown panel (gear boxes, Add-column, Map-ID), the
   snapshot modal, and the assistant settings modal — plus DataTables
   global-search boxes, download controls and (listed only) plotly /
   network gesture surfaces.
2. **Store-authoritative**: full store snapshot + capability registry.
3. Cross-map: which DOM inputs each agent write actually moved
   (empirical before/after diff of `Shiny.shinyapp.$inputValues`).
4. Every DOM widget NOT moved by any store binding becomes a `dom-only`
   case: the agent's `set_widgets` is attempted with plausible canonical
   ids and the exact rejection receipts are recorded.

## Case categories

- `binding` — every agent-writable store widget (~1xx cases), probe value
  chosen per kind (flip / first-allowed / +1 / marker string / mapping
  probe), numeric range rejections auto-retried once.
- `semantic-tool` — the five write tools + read/discovery tools, incl.
  deliberately invalid inputs (bogus pathway id) to capture the
  self-correction channel.
- `dom-only` — user-facing widgets with no store binding (buttons, search
  boxes, snapshot modal controls, assistant UI, dataset selector, …).
- `out-of-scope` — plotly/forceNetwork/iheatmap gesture surfaces. NOT
  tested (explicit instruction: no mouse-gesture mimicry); listed so the
  inventory itself is complete.

## Failure rule

Any observable UI effect wait exceeding **30 s** ⇒ status `timeout`
(recorded as failed).

## Outputs

- `cases/<NNN>-<widget>.json` — one test definition per widget
- `results/<NNN>-<widget>.json` — one result per widget: raw tool receipt,
  UI diff (changed inputs before/after), store read-back, duration,
  console errors, status, note
- `artifacts/dom_inventory.json`, `artifacts/store_snapshot.json`,
  `artifacts/capability_registry.json`, `artifacts/mapping_findings.json`,
  `artifacts/overview_state.json`, `artifacts/internal_store_ids.json`
- `summary.json` / `summary.md`

## Run

```bash
cd tests/e2e_agent/widget_sweep
node sweep_harness.mjs        # ~20-40 min; kills its own app on exit
```
