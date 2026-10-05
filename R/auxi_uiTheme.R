#' omicsViewer UI theme: iconography + namespaced CSS
#'
#' Returns a single \code{<style>} tag with the omicsViewer visual theme:
#' accent-colored navigation bars, tab icons, button/selector polish and
#' slim scrollbars.
#'
#' @section CSS scoping contract:
#' Every rule is selected through names this package owns - the root class
#' \code{.omicsviewer-app} (attached by \code{\link{app_ui}}) or custom
#' \code{.omicsviewer-*} classes - never through bare Bootstrap class
#' names. This keeps the theme stable when the viewer is embedded into a
#' host Shiny application: host stylesheets may redefine generic classes
#' such as \code{.nav-tabs} or \code{.btn-default}, but the more specific
#' namespaced selectors keep applying inside the viewer subtree only, and
#' nothing leaks out into the host UI. The snapshot modal lives outside
#' the subtree (rendered at body level), so it carries its own
#' \code{.omicsviewer-modal} hook class.
#'
#' @return A \code{shiny.tag} \code{<style>} element.
#' @keywords internal
omicsviewer_ui_theme <- function() {
  tags$style(HTML("
/* ================================================================
   omicsViewer theme - all selectors are namespaced:
   .omicsviewer-app ... (in-app subtree) or .omicsviewer-modal ...
   (app-owned modal hook). No bare Bootstrap class rules.
   ================================================================ */
.omicsviewer-app {
  --ov-accent: #0e7490;
  --ov-accent-strong: #155e75;
  --ov-ink: #16394f;
  --ov-ink-soft: #52606d;
  --ov-line: #e4ebf1;
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
    'Helvetica Neue', Arial, sans-serif;
  color: #243b53;
}

/* ---- main navigation bars (Data / Analysis panels) ---- */
.omicsviewer-app .navbar-default {
  background-color: #fff;
  border: none;
  border-bottom: 2px solid var(--ov-line);
  box-shadow: 0 1px 2px rgba(16, 42, 67, 0.07);
  margin-bottom: 10px;
  min-height: 38px;
}
/* two-row header: the brand block (title + one-line guidance) takes its
   own row, the tab bar flows directly underneath with a tight gap.
   Bootstrap floats .navbar-header left at >=768px, which would let the
   tabs sit BESIDE the title again on wide screens - unfloat the header
   so the split holds at every width. */
.omicsviewer-app .navbar-default .navbar-header {
  float: none;
}
.omicsviewer-app .navbar-default .navbar-brand {
  display: block;
  float: none;
  height: auto;
  color: var(--ov-accent-strong);
  font-size: 15px;
  font-weight: 600;
  line-height: 1.3;
  /* left padding 17px lines the brand block up with the first tab
     label: container edge (15px) + 17px = 32px, the same x-position
     the first .navbar-nav icon/text starts at (nav margin -6px +
     link padding 8px) - keeps title and tabset names vertically
     aligned as one column. */
  padding: 7px 17px 3px 17px;
}
.omicsviewer-app .navbar-default .navbar-brand:hover,
.omicsviewer-app .navbar-default .navbar-brand:focus {
  color: #0c4a5e;
}
.omicsviewer-app .navbar-default .navbar-brand i {
  margin-right: 5px;
}
.omicsviewer-app .navbar-default .navbar-brand .omicsviewer-nav-sub {
  font-size: 12px;
  font-weight: 400;
  color: #627d98;
}
.omicsviewer-app .navbar-default .navbar-nav {
  margin: 0 -15px 0 -6px;
}
.omicsviewer-app .navbar-default .navbar-nav > li > a {
  color: var(--ov-ink-soft);
  font-size: 13px;
  padding: 9px 8px;
  transition: color 0.15s ease, background-color 0.15s ease;
}
.omicsviewer-app .navbar-default .navbar-nav > li > a:hover,
.omicsviewer-app .navbar-default .navbar-nav > li > a:focus {
  color: var(--ov-accent);
  background-color: rgba(14, 116, 144, 0.07);
}
.omicsviewer-app .navbar-default .navbar-nav > .active > a,
.omicsviewer-app .navbar-default .navbar-nav > .active > a:hover,
.omicsviewer-app .navbar-default .navbar-nav > .active > a:focus {
  color: var(--ov-accent);
  font-weight: 600;
  background-color: transparent;
  box-shadow: inset 0 -3px 0 0 var(--ov-accent);
}
.omicsviewer-app .navbar-default .navbar-nav > li > a i {
  font-size: 12px;
  margin-right: 4px;
  opacity: 0.85;
}
.omicsviewer-app .navbar-default .navbar-nav > .active > a i {
  opacity: 1;
}
/* On narrower screens drop only the guidance sentence so the header
   stays one short line (the tab bar wraps as before). */
@media (max-width: 1499px) {
  .omicsviewer-app .navbar-default .navbar-brand .omicsviewer-nav-sub {
    display: none;
  }
}

/* ---- inner tabsets (heatmap controls, dropdown panels) ---- */
.omicsviewer-app .nav-tabs {
  border-bottom: 1px solid var(--ov-line);
}
.omicsviewer-app .nav-tabs > li > a {
  color: var(--ov-ink-soft);
  padding: 7px 12px;
  border-radius: 7px 7px 0 0;
}
.omicsviewer-app .nav-tabs > li > a:hover {
  color: var(--ov-accent);
  background-color: #f7fafc;
  border-color: var(--ov-line) var(--ov-line) #b9c6d1;
}
.omicsviewer-app .nav-tabs > li.active > a,
.omicsviewer-app .nav-tabs > li.active > a:hover,
.omicsviewer-app .nav-tabs > li.active > a:focus {
  color: var(--ov-accent);
  font-weight: 600;
  border: 1px solid var(--ov-line);
  border-bottom-color: #fff;
}

/* ---- buttons ----
   Spacelab ships a CHARCOAL .btn-default (dark bg, white text); a host
   theme may ship something else again. The viewer therefore owns the
   COMPLETE color pair (background + gradient flattened + text + border)
   for every button variant it renders, so contrast can never depend on
   whichever .btn-* rules win outside. */
.omicsviewer-app .btn-default {
  background-color: #fff;
  background-image: none;
  color: var(--ov-ink-soft);
  border-color: #d3dde6;
}
.omicsviewer-app .btn-default:hover,
.omicsviewer-app .btn-default:focus,
.omicsviewer-app .btn-default:active,
.omicsviewer-app .btn-default.active {
  background-color: #eef5f8;
  background-image: none;
  color: var(--ov-accent);
  border-color: var(--ov-accent);
}
.omicsviewer-app .btn-primary {
  background-color: var(--ov-accent);
  background-image: none;
  color: #fff;
  border-color: #0b5f78;
}
.omicsviewer-app .btn-primary:hover,
.omicsviewer-app .btn-primary:focus,
.omicsviewer-app .btn-primary:active,
.omicsviewer-app .btn-primary.active {
  background-color: var(--ov-accent-strong);
  background-image: none;
  color: #fff;
  border-color: #0c4a5e;
}
.omicsviewer-app .btn-info {
  background-color: var(--ov-accent);
  background-image: none;
  color: #fff;
  border-color: #0b5f78;
}
.omicsviewer-app .btn-info:hover,
.omicsviewer-app .btn-info:focus,
.omicsviewer-app .btn-info:active,
.omicsviewer-app .btn-info.active {
  background-color: var(--ov-accent-strong);
  background-image: none;
  color: #fff;
  border-color: #0c4a5e;
}
.omicsviewer-app .btn-danger {
  background-color: #d9230f;
  background-image: none;
  color: #fff;
  border-color: #ba1f10;
}
.omicsviewer-app .btn-danger:hover,
.omicsviewer-app .btn-danger:focus,
.omicsviewer-app .btn-danger:active,
.omicsviewer-app .btn-danger.active {
  background-color: #ba1f10;
  background-image: none;
  color: #fff;
  border-color: #98190d;
}
.omicsviewer-app .btn-warning {
  background-image: none;
}

/* ---- compact widget sizing ----------------------------------------
   The viewer's tool widgets share ONE small size (matching the
   quick-view badges and the AI drawer buttons): the table toolbar
   (Show all / multi-selection switch / Save table / Add column), every
   shinyWidgets dropdown toggle (the figure-attribute gear, the heatmap
   Controls buttons), and all buttons inside the snapshot modals. The
   containers are selected through app-owned names only: the
   .omicsviewer-toolbar class set in dataTable_ui(), the
   .omicsviewer-modal hook from L0_module_snapshot.R, and the
   shinyWidgets .sw-dropdown wrapper scoped under our root. */
.omicsviewer-app .omicsviewer-toolbar .btn,
.omicsviewer-app .sw-dropdown > .btn,
.omicsviewer-modal .btn {
  padding: 1px 8px;
  font-size: 12px;
  line-height: 1.6;
  border-radius: 7px;
}
.omicsviewer-app .sw-dropdown > .btn i {
  font-size: 11px;
}

/* top-right floating tools (dataset export, snapshots): full color
   pair + the compact footprint */
.omicsviewer-app .omicsviewer-tool-btn {
  background-color: #fff;
  background-image: none;
  color: var(--ov-ink-soft);
  border-radius: 9px;
  font-size: 12px;
  padding: 3px 9px;
  box-shadow: 0 1px 2px rgba(16, 42, 67, 0.10);
}
.omicsviewer-app .omicsviewer-tool-btn:hover {
  background-color: #eef5f8;
  box-shadow: 0 2px 6px rgba(14, 116, 144, 0.25);
}
.omicsviewer-app .omicsviewer-tool-btn i {
  margin-right: 3px;
}

/* multi-selection switch inside the table toolbar: soften the frame to
   match the compact buttons */
.omicsviewer-app .omicsviewer-toolbar .bootstrap-switch {
  border-color: #d3dde6;
  border-radius: 7px;
}
.omicsviewer-app .omicsviewer-toolbar .bootstrap-switch
.bootstrap-switch-handle-on,
.omicsviewer-app .omicsviewer-toolbar .bootstrap-switch
.bootstrap-switch-handle-off,
.omicsviewer-app .omicsviewer-toolbar .bootstrap-switch
.bootstrap-switch-label {
  font-size: 11px;
  color: var(--ov-ink-soft);
}

/* ---- app title block (dataset summary line) ---- */
/* padding-left matches the column(6) inset used by the Data/Analysis
   panels below, so the app title lines up with the left edge of the
   panel boxes instead of touching the page border. */
.omicsviewer-app .omicsviewer-titleblock {
  padding-left: 15px;
}
.omicsviewer-app .omicsviewer-app-title {
  display: inline;
  font-size: 26px;
  font-weight: 700;
  letter-spacing: 0.3px;
  color: var(--ov-ink);
  margin: 0;
}
.omicsviewer-app .omicsviewer-app-title i {
  color: var(--ov-accent);
  font-size: 22px;
  margin-right: 7px;
  vertical-align: 2px;
}
.omicsviewer-app .omicsviewer-app-sub {
  display: inline;
  font-size: 16px;
  font-weight: 400;
  color: #627d98;
  margin-left: 4px;
}

/* ---- dataset selector (selectize) ---- */
.omicsviewer-app .selectize-input {
  border-radius: 9px;
  border-color: #d3dde6;
  box-shadow: none;
}
.omicsviewer-app .selectize-input.focus {
  border-color: var(--ov-accent);
  box-shadow: 0 0 0 3px rgba(14, 116, 144, 0.15);
}
.omicsviewer-app .selectize-dropdown {
  border-color: #d3dde6;
  border-radius: 9px;
  box-shadow: 0 10px 28px rgba(16, 42, 67, 0.16);
  overflow: hidden;
}
.omicsviewer-app .selectize-dropdown .active {
  background-color: #e4f3f7;
  color: var(--ov-accent);
}

/* ---- keyboard focus ring ---- */
.omicsviewer-app a:focus-visible,
.omicsviewer-app .btn:focus-visible {
  outline: 3px solid rgba(14, 116, 144, 0.45);
  outline-offset: 1px;
}

/* ---- slim scrollbars inside the viewer subtree ---- */
.omicsviewer-app ::-webkit-scrollbar {
  width: 9px;
  height: 9px;
}
.omicsviewer-app ::-webkit-scrollbar-thumb {
  background-color: #c3ced9;
  border-radius: 8px;
}
.omicsviewer-app ::-webkit-scrollbar-thumb:hover {
  background-color: #9fb1c1;
}
.omicsviewer-app ::-webkit-scrollbar-track {
  background: transparent;
}

/* ---- snapshot modals (rendered at body level, outside the root
   subtree; selected via the app-owned .omicsviewer-modal hook class
   that modalDialog() calls in L0_module_snapshot.R append) ---- */
.omicsviewer-modal .modal-content {
  border: none;
  border-radius: 12px;
  box-shadow: 0 18px 50px rgba(16, 42, 67, 0.28);
}
.omicsviewer-modal .modal-header {
  border-bottom: 1px solid var(--ov-line, #e4ebf1);
  padding: 14px 18px;
}
.omicsviewer-modal .modal-title {
  color: var(--ov-ink, #16394f);
  font-weight: 600;
}
.omicsviewer-modal .modal-title i {
  color: var(--ov-accent, #0e7490);
  margin-right: 6px;
}
.omicsviewer-modal .modal-body {
  padding: 16px 18px;
}
/* modal buttons live outside the .omicsviewer-app subtree: give them the
   same owned color pairs so the spacelab charcoal (or a host theme)
   cannot darken them */
.omicsviewer-modal .btn-default {
  background-color: #fff;
  background-image: none;
  color: #52606d;
  border-color: #d3dde6;
}
.omicsviewer-modal .btn-default:hover,
.omicsviewer-modal .btn-default:focus,
.omicsviewer-modal .btn-default:active,
.omicsviewer-modal .btn-default.active {
  background-color: #eef5f8;
  background-image: none;
  color: #0e7490;
  border-color: #0e7490;
}
.omicsviewer-modal .btn-primary {
  background-color: #0e7490;
  background-image: none;
  color: #fff;
  border-color: #0b5f78;
}
.omicsviewer-modal .btn-primary:hover,
.omicsviewer-modal .btn-primary:focus,
.omicsviewer-modal .btn-primary:active,
.omicsviewer-modal .btn-primary.active {
  background-color: #155e75;
  background-image: none;
  color: #fff;
  border-color: #0c4a5e;
}
.omicsviewer-modal .btn-danger {
  background-color: #d9230f;
  background-image: none;
  color: #fff;
  border-color: #ba1f10;
}
.omicsviewer-modal .btn-danger:hover,
.omicsviewer-modal .btn-danger:focus,
.omicsviewer-modal .btn-danger:active,
.omicsviewer-modal .btn-danger.active {
  background-color: #ba1f10;
  background-image: none;
  color: #fff;
  border-color: #98190d;
}
"))
}
