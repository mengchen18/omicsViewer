#!/usr/bin/env bash
# =============================================================================
# run_regression.sh — WP-0: the widget-sweep regression gate (FIX_PLAN.md).
# =============================================================================
# Boots the omicsViewer app (demo.RDS), runs the exhaustive agent widget sweep
# (sweep_harness.mjs — every user-facing widget driven through the exact LLM
# tool surface, validated browser-side), then diffs the status matrix against
# the committed baseline_summary.json and fails on any regression:
#   * binding timeout / no-ui-effect / no-probe / rejected / unchanged / error
#     count increase
#   * dom-only unaddressable count increase
#   * any case transitioning success -> non-success
# Improvements and vanished/new cases are reported but never fail the gate.
#
# The sweep writes to a scratch directory (WIDGET_SWEEP_OUT), so the committed
# evidence artifacts (cases/, results/, summary.json) are never clobbered.
# Regressed cases are automatically re-run once in isolation (congestion
# guard) before the gate declares failure.
#
# Usage:
#   NOT_CRAN=true tests/e2e_agent/widget_sweep/run_regression.sh [options]
# Options:
#   --install       reinstall the package from source first (AGENTS.md: the
#                   sweep exercises the INSTALLED library, not the work tree)
#   --retries N     isolation re-runs of regressed cases (default 1)
#   --only IDS      comma-separated caseId substrings; run only those cases
#                   (gated against the baseline exactly like a full run)
#   --baseline F    alternative baseline file (e.g. a candidate refreshed
#                   baseline under review); default: baseline_summary.json
#   --keep          keep the scratch directory even on success
#   --selftest      run only the offline diff-engine selftest (no browser)
# Environment:
#   NOT_CRAN=true   required (same gating convention as tests/test_agentUiEffects.R)
# Needs: node >= 18, playwright module (npm install in tests/e2e_agent),
#        system Chrome, R with the omicsViewer package installed.
# =============================================================================
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
PORT=7794

say() { printf '\n== %s ==\n' "$*"; }
die() { printf '\nGATE FAIL: %s\n' "$*" >&2; exit 1; }

# --- offline selftest --------------------------------------------------------
if [[ "${1:-}" == "--selftest" ]]; then
  node "$HERE/check_regression.mjs" --selftest || die "diff-engine selftest"
  echo "GATE selftest PASS"
  exit 0
fi

# --- NOT_CRAN gate (never inside R CMD check) --------------------------------
if [[ "${NOT_CRAN:-}" != "true" ]]; then
  echo "widget-sweep regression gate: skipped (set NOT_CRAN=true to run)"
  exit 0
fi

# --- option parse ------------------------------------------------------------
DO_INSTALL=0
RETRIES=1
ONLY=""
KEEP=0
BASELINE="$HERE/baseline_summary.json"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --install) DO_INSTALL=1 ;;
    --retries) RETRIES="${2:?--retries needs a value}"; shift ;;
    --only)    ONLY="${2:?--only needs a value}"; shift ;;
    --baseline) BASELINE="${2:?--baseline needs a value}"; shift ;;
    --keep)    KEEP=1 ;;
    *) die "unknown option: $1" ;;
  esac
  shift
done
[[ -f "$BASELINE" ]] || die "baseline file not found: $BASELINE"

# --- environment checks (skip cleanly, like tests/test_agentUiEffects.R) ----
skip() { echo "widget-sweep regression gate: skipped ($*)"; exit 0; }
command -v node >/dev/null 2>&1 || skip "node unavailable"
[[ -f "$REPO/tests/e2e_agent/node_modules/playwright/package.json" ]] || skip "playwright module unavailable (npm install in tests/e2e_agent)"
CHROME=""
for c in /usr/bin/google-chrome /usr/bin/google-chrome-stable /usr/bin/chromium /usr/bin/chromium-browser; do
  [[ -x "$c" ]] && CHROME="$c" && break
done
[[ -n "$CHROME" ]] || skip "no Chrome/Chromium found"
command -v Rscript >/dev/null 2>&1 || skip "Rscript unavailable"
Rscript -e 'invisible(packageVersion("omicsViewer"))' >/dev/null 2>&1 || skip "omicsViewer package not installed"

# The sweep drives the INSTALLED library: warn on version drift.
INST_VER=$(Rscript -e 'cat(as.character(packageVersion("omicsViewer")))' 2>/dev/null | tail -1)
SRC_VER=$(sed -n 's/^Version: *//p' "$REPO/DESCRIPTION")
if [[ "$INST_VER" != "$SRC_VER" ]]; then
  echo "WARNING: installed omicsViewer $INST_VER != source $SRC_VER — the gate tests the INSTALLED library (use --install)" >&2
fi

# --- port hygiene (orphaned R processes poison runs; see AGENTS.md) ----------
sweep_port() {
  local pids
  pids=$(ss -tlnp 2>/dev/null | grep ":$PORT " | grep -oP 'pid=\K[0-9]+' | sort -u || true)
  if [[ -n "$pids" ]]; then
    echo "killing stale process(es) on port $PORT: $(echo "$pids" | tr '\n' ' ')"
    kill -9 $pids 2>/dev/null || true
    sleep 1
  fi
}
sweep_port
trap 'sweep_port' EXIT

# --- optional package (re)install --------------------------------------------
if [[ "$DO_INSTALL" == 1 ]]; then
  say "installing omicsViewer from source"
  Rscript -e "install.packages('$REPO', repos = NULL, type = 'source')" || die "package install failed"
fi

# --- run the sweep into a scratch dir ----------------------------------------
WORK=$(mktemp -d /tmp/omicsviewer-sweep-gate.XXXXXX)
echo "scratch dir: $WORK"
run_sweep() { # $1 = out dir, $2 = ONLY filter ("" = full)
  say "running sweep_harness.mjs (out: $1, filter: ${2:-<full>})"
  local start=$SECONDS
  local env_args=(WIDGET_SWEEP_OUT="$1")
  if [[ -n "${2:-}" ]]; then env_args+=(WIDGET_SWEEP_ONLY="$2"); fi
  if ! env "${env_args[@]}" node "$HERE/sweep_harness.mjs" 2>&1 | tee "$1/sweep_run.log"; then
    tail -30 "$1/sweep_run.log" || true
    die "sweep_harness.mjs failed (log: $1/sweep_run.log)"
  fi
  echo "sweep finished in $((SECONDS - start))s"
}

run_sweep "$WORK" "$ONLY"

DIRS=("$WORK/results")
# RETRIES = number of isolation re-rounds; the loop performs RETRIES+1
# checks (first full check, then one check per re-run).
for attempt in $(seq 0 "$RETRIES"); do
  say "diff against baseline (check $((attempt + 1))/$((RETRIES + 1)))"
  REG="$WORK/regressions-$((attempt + 1)).json"
  CHECK_ARGS=(--baseline "$BASELINE" --write-regressions "$REG")
  for d in "${DIRS[@]}"; do CHECK_ARGS+=(--results "$d"); done
  set +e
  node "$HERE/check_regression.mjs" "${CHECK_ARGS[@]}"
  RC=$?
  set -e
  if [[ "$RC" == 0 ]]; then
    say "GATE PASS — no regression against $BASELINE"
    if [[ "$KEEP" != 1 ]]; then rm -rf "$WORK"; else echo "scratch kept: $WORK"; fi
    exit 0
  fi
  if [[ "$RC" != 1 ]]; then
    die "check_regression.mjs usage/IO error (rc=$RC)"
  fi
  IDS=$(node -e "const j=require('$REG'); process.stdout.write(j.caseIds.join(','))")
  if [[ -z "$IDS" ]]; then
    die "count-level regression with no retriable case ids — inspect $WORK"
  fi
  if [[ "$attempt" == "$RETRIES" ]]; then
    echo "regressed case(s), still failing after $RETRIES isolation re-run(s): $IDS"
    break
  fi
  say "retrying regressed cases in isolation (round $((attempt + 1))/$RETRIES): $IDS"
  RETRY="$WORK/retry-$((attempt + 1))"
  mkdir -p "$RETRY"
  run_sweep "$RETRY" "$IDS"
  DIRS+=("$RETRY/results")
done

echo
tail -40 "$WORK/sweep_run.log" || true
die "regression detected — scratch dir kept for inspection: $WORK (results: ${DIRS[*]}; diff details above)"
