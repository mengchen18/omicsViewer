# HANDOVER — scatter quick-view switch multi-flash fix (IN PROGRESS, BROKEN TREE)

**Bug 1 (volcano corner not selected on load) is FIXED and confirmed by
the user — do not touch it.** The only remaining problem is the
multiple-render/flash on quick-view switches, and the current working
tree (an attempt at that fix) hangs in an infinite reactive loop.

Written 2026-09-24 mid-fix. **The working tree is intentionally left in a
BROKEN state (infinite reactive loop) for the next session to debug or
bisect.** Read this first; then `HANDOVER.md` §"Scatter volcano corner"
for the already-committed background; `AGENT_ACCURACY_PLAN.md` for the
control-plane architecture. Nothing in this file is committed except
this document itself.

## Reproducing the hang (current tree)

```bash
Rscript -e '
  options(shiny.port = 7770, shiny.host = "127.0.0.1")
  eset <- readRDS("inst/extdata/demo.RDS")
  omicsViewer::omicsViewer(dir = "inst/extdata/", ESVObj = eset)
' > /tmp/hang.log 2>&1 &
# open http://127.0.0.1:7770 in a browser (or Playwright with
# --disable-webgl --disable-webgl2)
```

Signature (measured 2026-09-24, port 7770):
- the UI shell loads (tabs, dataset-loaded banner render);
- `Shiny.shinyapp.$inputValues['app-dataspace-eset']` DOES appear ("Feature");
- the shiny busy indicator stays up forever; **no plotly container ever
  renders**;
- the R process pins one core (~95 % CPU) — an infinite reactive loop;
- `/tmp/hang.log` shows NO error, just "Listening on…" + "Using model…".

## Git state

- Last good commit: `aef2796` (fix part 1) + `b6b9250` (docs). That state:
  load corner FIXED (verified), cor→volcano exactly 1 render, volcano→cor
  ≤2 renders, **but volcano→volcano (e.g. RE_vs_ME → RE_vs_LE) flashed
  multiple times** — user-reported regression of the convergence gate
  (pre_volcano dipped TRUE→FALSE→TRUE on every switch; see below).
- Dirty files (the broken second iteration): `R/module_figureAttr4.R`,
  `R/module_meta_scatter.R` (89 insertions / 40 deletions total).
- To get back to the last-good state: `git stash` (keep the stash! it is
  the material to re-apply incrementally).

## Design intent of the broken iteration (all three points are REQUIRED)

1. **Corner state changes ONLY on volcano-ness transitions of the ATOMIC
   store axes.** `pre_vol` must be the pure store-watcher predicate (no
   convergence gate inside it) — otherwise volcano→volcano dips
   FALSE→TRUE and re-fires the corner chain (the user-visible regression
   of commit aef2796).
2. **The corner applies at display convergence, INSIDE the same reactive
   recompute as the axes render** (one paint, rects included). An
   observer-based apply lands one flush LATER (verified empirically:
   volcano→cor painted `Cor|r2` then `Cor|r0`; cor→volcano painted
   `volcano|r0` then `volcano|r2` — 2 renders each, no oscillation but
   the corner was late). Hence the reactive-consumption design:
   attr4 exposes `params$cutoff_reactive` (a pure reactive resolving
   intent+gate+applied corner), and meta_scatter's `rectval` consumes it
   instead of the `params$cutoff` mirror.
3. **Content-deduped side effects** (scorner widget updateSelectInput,
   store_apply scorner, params$cutoff mirror): the historical
   `attr(l,"seed") <- Sys.time()` forced a full plot redraw on every
   input event even when the corner/cutoffs were identical.

## What the broken tree contains (git diff summary)

`R/module_meta_scatter.R`:
- `pre_vol`: pure store-watcher predicate (GOOD — keep this).
- `.scatter_axes_converged`: new reactive (v1/v2 triples == store
  watchers) passed as `corner_apply_gate` to attr4 (GOOD concept).
- `rectval` reads `attr4select$cutoff_reactive` (function) if present,
  else falls back to `attr4select$cutoff` (SUSPECT — see loop analysis).

`R/module_figureAttr4.R`:
- new param `corner_apply_gate = reactive(TRUE)` (roxygen updated).
- `pendingCorner` reactiveVal, set ONLY by `observeEvent(pre_volcano())`
  transitions ("volcano"/"None" intent).
- `corner_effective` reactive: `intent + gate` → intent, else
  `params$cutoff$corner %||% "None"` (**reads params$cutoff**).
- `cutoff_effective` reactive: `list(x=val_xcut(), y=val_ycut(),
  corner=corner_effective())`.
- side-effect `observe`: reads `cutoff_effective()`, computes a content
  key, and on change writes `.a4_cutoff_key`, **`params$cutoff <- l`**,
  `store_apply(scorner, mark_pending=FALSE)`,
  `updateSelectInput("scorner", selected=corner)` (**SUSPECT**).
- the generic user-input observer also writes `params$cutoff` under the
  same content key.
- `params$cutoff_reactive <- cutoff_effective` exported before return.

## Loop analysis (where to look first)

The graph contains a read/write cycle: `corner_effective` READS
`params$cutoff`; the side-effect observe WRITES `params$cutoff`. A
reactiveValues write always invalidates readers (no diffing), so each
apply re-triggers `corner_effective → cutoff_effective → observe`; the
content key is supposed to stop it after one extra pass. It doesn't —
suspects, in order:

1. **The key never stabilizes**: `l$x`/`l$y` come from the DEBOUNCED
   `val_xcut()`/`val_ycut()`. While inputs are initializing these emit
   NULL→value transitions; but more importantly `paste()` of a numeric
   should be stable… verify by logging the key inside the observe.
2. **`updateSelectInput` ping-pong**: every apply sends a scorner update;
   if the browser (selectize) fires an input event even for an unchanged
   value, the generic observer re-runs; it reads `val_xcut()` etc. —
   cheap — but if anything in that chain rewrites `params$cutoff` with a
   fresh list under a key that differs by representation (e.g.
   "0.301029995663981" vs "0.30102999566398120"), the cycle never ends.
   float→paste formatting is a classic source.
3. **`store_apply` → epoch bump → .a4_epoch push observer → scorner
   push**: the attr4 push observer runs on ANY epoch bump and pushes
   scorner when `pending[[scorner]]` is set. `mark_pending=FALSE` means
   no pending — but CHECK the `.a4_epoch` observer's condition actually
   respects that for xcut/ycut too (it pushes those when pending; the
   SEED no longer sets pending, but the RESTORE path does — during init
   no restore runs, so probably fine).
4. **rectval reading a function out of reactiveValues**:
   `attr4select$cutoff_reactive` — reading the member registers a
   dependency on that member; nothing writes it, so it should be inert —
   but confirm Shiny doesn't choke on function-valued reactiveValues
   members during serialization of `params` (the status observer
   snapshots `params` fields — check nothing iterates ALL params members
   into a snapshot; `params$status` is built explicitly, so probably
   fine).

**Fastest debug**: add `message()` counters (with `format(Sys.time(),
"%OS3")` and the computed key) at the top of (a) the side-effect observe,
(b) the generic input observer, (c) the choices-rebuild observer, then
load the page and read the R log — the repeating message names the loop.
Alternatively `profvis::profvis` a few seconds of the spinning process.
Alternatively bisect: comment out (1) the `updateSelectInput` in the
side-effect observe, then (2) the `store_apply`, then (3) the
`params$cutoff <- l` write, reload after each — the loop's driver will
reveal itself.

## Verified-correct facts (do not re-derive)

- `store_read` is a PLAIN snapshot (never invalidates) — reactive
  consumers MUST use `store_watch` (per-key reactives) or `store_epoch`.
- Plotly.react in the app's bundled plotly version DOES clear layout
  shapes when the new layout omits the key (tested in-page). "No rects"
  figures are fine; the visible rects always come from rectval itself.
- `plotly::layout(shapes=NULL/empty)` drops the key entirely (both
  cases) — irrelevant but explains rect serialization.
- The attr4 SEED must keep `mark_pending = FALSE` (a value copied FROM a
  live input can never be acked — its re-assert clobbered the corner at
  load; that fix is committed and verified).
- The triselector UI→store sync must mirror only COHERENT triples
  (`.scatter_triple_coherent`, committed) — mixed cascade echoes were
  treated as user overrides, cleared pending, and oscillated the store
  after every switch (the ROOT of the original user-reported flash).
- Commit aef2796's Tier A additions assert: load corner state
  (scorner=volcano, 2 shapes, badge active), volcano→cor ≤2 renders
  ending `|0`, cor→volcano "exactly 1 render" — NOTE: under the
  intended final design cross-paradigm switches are also 1 render
  (corner resolves in-reactive), so keep/adjust the expectation to
  `exactly 1` for both, and ADD the volcano→volcano case (exactly 1
  render, rects stay 2 — the user's regression report).

## Acceptance for the finished fix (verify in a live browser, then Tier A)

- load: volcano badge active, scorner=volcano, 2 rects, corner features
  selected;
- every quick-view switch (volcano↔cor, volcano↔volcano): exactly ONE
  distinct plot render, no intermediate garbage frames (no rects on a
  cor plot, no rects-flicker), corner state unchanged between two
  volcano views;
- unit board + Tier A ×2 green.

## Environment reminders

- ports 7775–7799: kill orphans before/after browser runs
  (`ss -tlnp | grep -E ':77[0-9][0-9]'`); an orphan serving an OLD build
  poisons verification runs (bit us twice today).
- Chrome headless flags: `--disable-webgl --disable-webgl2` (mandatory).
- Reinstall after every source change:
  `Rscript -e 'install.packages("<repo>", repos = NULL, type = "source")'`.
