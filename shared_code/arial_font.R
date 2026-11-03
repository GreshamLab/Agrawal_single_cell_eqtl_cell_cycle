# Publication font for Figures 3-4: Arial, matching Figures 2 and 5 (PUB_FONT <- "Arial").
#
# Plots in the Figure 3/4 scripts call theme_classic()/theme_bw()/theme_minimal()/theme_void()
# without base_family, so they rendered in the cluster's default sans (Noto Sans). Sourcing
# this file once makes Arial the default everywhere without editing individual plots:
#   - ggplot2 complete themes get base_family = "Arial" (geom_text/geom_label/annotate
#     inherit it through from_theme(family) in ggplot2 >= 4.0)
#   - ggrepel and ggsignif text, which ignore the theme, get family = "Arial"
# The cluster has no Arial; systemfonts resolves "Arial" to Arimo (metric-compatible), the
# same font Figures 2 and 5 render in.
local({
  env <- globalenv()
  for (fn in c("theme_bw", "theme_classic", "theme_minimal", "theme_void", "theme_light",
               "theme_linedraw", "theme_gray", "theme_grey", "theme_dark")) {
    f <- get(fn, envir = asNamespace("ggplot2"))
    formals(f)$base_family <- "Arial"
    assign(fn, f, envir = env)
  }
  ggplot2::theme_set(ggplot2::theme_grey(base_family = "Arial"))
  # pdf()/postscript() only know their built-in font families ("invalid font type" otherwise);
  # map Arial to Helvetica, which has the same metrics
  if (!"Arial" %in% names(grDevices::pdfFonts()))
    grDevices::pdfFonts(Arial = grDevices::pdfFonts()$Helvetica)
  if (!"Arial" %in% names(grDevices::postscriptFonts()))
    grDevices::postscriptFonts(Arial = grDevices::postscriptFonts()$Helvetica)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggplot2::update_geom_defaults(ggrepel::GeomTextRepel, list(family = "Arial"))
    ggplot2::update_geom_defaults(ggrepel::GeomLabelRepel, list(family = "Arial"))
  }
  if (requireNamespace("ggsignif", quietly = TRUE)) {
    gs <- ggsignif::geom_signif
    formals(gs)$family <- "Arial"
    assign("geom_signif", gs, envir = env)
  }
})
