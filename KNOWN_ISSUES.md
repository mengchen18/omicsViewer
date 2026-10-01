# Known issues — snapshot save/reload (Stage 1, todo.md §2; 2026-09-30)

## 1. Plotly-internal events are not restored — by design

Lasso/box **shapes** drawn in the browser and other plotly-internal event
state are browser-owned display state and are deliberately NOT part of a
snapshot. Only the **semantic selection** (feature/sample ids, origin,
anchor) rides the snapshot through the selection bus. After a restore the
selection is emphasized via the anchor path; the drawn shape itself is not
recreated. This matches the display contract in
`R/module_meta_scatter.R` ("Plotly owns the visible box/lasso ... emphasis
comes from the settled selVal").

**Headless verification limit:** plotly event inputs cannot be simulated
under plotly >= 4.12 because event registration happens client-side
(`event_register` on the rendered figure); `event_data()` returns NULL
server-side with a "not registered" warning. This is the same root cause
behind the **pre-existing** `tests/test_scatterSelection.R` failures on
this machine (14 not-ok on unmodified master). Figure-origin selections
are therefore covered indirectly (golden v2 fixture carries origin
"figure" for the sample space); live verification happens in the Tier A
browser suite.

## 2. Hand-built v2 fixture selection records — RESOLVED (todo 4.1, 2026-09-30)

The two formerly commented-out assertions in `tests/test_snapshotRoundTrip.R`
case 9 ("v2 selection records restore", "v2 record origin survives the round
trip") are re-enabled and green. Root cause was the restore applying the
selection in the SAME flush as the axes cascade, so a corner-echo clear ate
the records of a hand-built v2 fixture (no feature-fig panel status). The
phased restore controller (todo 4.1 / deferred 2.8) now applies the
selection through the selection bus only after the restored axes have
converged, which sequences the corner re-engagement after the record write.
Suite now 63 assertions.

## 3. Transient "--select--" sentinel may rest in the widget store

When a restore clears an attr4 colour mapping that was never picked in
this session, the canonical unset sentinel ("--select--") can remain as
the stored value of the `attr4.<group>_variable` key (the unset marker
cannot commit because the cascade's analysis/subset were never settled).
The widget and figure behave correctly (placeholder shown, no mapping)
and the panel status reports the mapping as unset; only a registry view
would show the sentinel string. Harmless; revisit if the agent-facing
registry should never expose it.

## 4. Snapshot function is experimental (surfaced in the UI)

The save/load modal and the restore-confirm dialog both display the
disclaimer: *"The snapshot function is experimental. Saved states may not
restore exactly across package versions or revised datasets.
Plotly-internal events (lasso/box shapes) are not restored — only the
semantic selection of features/samples is. Do not rely on snapshots as
the only record of an analysis."*

## 5. Pre-existing (not Stage 1): test_scatterSelection fails on master

14 of 20 assertions fail identically on unmodified master
(v2.1.57) — plotly 4.12 event registration (see §1). Not caused by the
Stage 1 changes; needs its own fix (e.g. registering events server-side
in the test harness or pinning plotly).
