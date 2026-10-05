# Known issues — snapshot save/reload

## 1. Plotly-internal events are not restored — by design

Lasso/box **shapes** drawn in the browser and other plotly-internal event
state are browser-owned display state and are deliberately NOT part of a
snapshot. Only the **semantic selection** (feature/sample ids, origin,
anchor) rides the snapshot. After a restore the selection is emphasized via
the anchor path; the drawn shape itself is not recreated.

## 2. Snapshot function is experimental (surfaced in the UI)

The save/load modal and the restore-confirm dialog both display the
disclaimer: *"The snapshot function is experimental. Saved states may not
restore exactly across package versions or revised datasets.
Plotly-internal events (lasso/box shapes) are not restored — only the
semantic selection of features/samples is. Do not rely on snapshots as
the only record of an analysis."*
