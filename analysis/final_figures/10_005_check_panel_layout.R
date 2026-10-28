# Measure the actual drawn plot areas and run the portable alignment gate.
# This helper does not change figure data, scales or geometry.
require_patchwork_panel_alignment <- function(plot, manifest_path, report_path,
    width_in, height_in, audit_script, python, strict=TRUE) {
  if (!file.exists(audit_script) || !nzchar(python))
    stop("Provide the repository alignment checker and a Python executable.")
  for (pkg in c("patchwork", "jsonlite"))
    if (!requireNamespace(pkg, quietly=TRUE)) stop("Missing R package: ", pkg)
  grob <- patchwork::patchworkGrob(plot)
  cells <- grob$layout[grepl("^panel(?:-[0-9]+)?$", grob$layout$name, perl=TRUE), ]
  cells <- cells[order(cells$t, cells$l), ]
  if (nrow(cells)<2) stop("At least two final panel viewports are required.")
  probe <- tempfile(fileext=".pdf")
  grDevices::pdf(probe, width=width_in, height=height_in, useDingbats=FALSE)
  opened <- TRUE
  on.exit({ if(opened) grDevices::dev.off(); unlink(probe) }, add=TRUE)
  grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  names <- grid::grid.ls(viewports=TRUE, grobs=FALSE, print=FALSE)$name
  panels <- lapply(seq_len(nrow(cells)), function(i) {
    cell <- cells[i, ]
    target <- paste0(cell$name,".",cell$t,"-",cell$l,"-",cell$b,"-",cell$r)
    if (!target %in% names) {
      matches <- names[startsWith(names,paste0(cell$name,"."))]
      if (length(matches)!=1) stop("Cannot uniquely identify panel viewport: ",target)
      target <- matches[1]
    }
    grid::seekViewport(target)
    lo <- grid::deviceLoc(grid::unit(0,"npc"),grid::unit(0,"npc"),valueOnly=TRUE)
    hi <- grid::deviceLoc(grid::unit(1,"npc"),grid::unit(1,"npc"),valueOnly=TRUE)
    grid::upViewport(0)
    list(id=letters[i],bbox_pt=c(lo$x,lo$y,hi$x,hi$y)*72,
         grid_id="patchwork-grid-1",row_start=as.integer(cell$t-1),
         row_stop=as.integer(cell$b),col_start=as.integer(cell$l-1),
         col_stop=as.integer(cell$r))
  })
  grDevices::dev.off(); opened <- FALSE
  manifest <- list(schema_version=1L,backend="r-patchwork",
    figure=list(width_pt=width_in*72,height_pt=height_in*72),
    panels=panels,exemptions=list())
  jsonlite::write_json(manifest,manifest_path,auto_unbox=TRUE,pretty=TRUE,digits=NA)
  cmd <- c(shQuote(audit_script),shQuote(manifest_path),"--json-out",shQuote(report_path),
           "--tolerance-pt","1.5")
  output <- suppressWarnings(system2(python,args=cmd,stdout=TRUE,stderr=TRUE))
  status <- attr(output,"status");if(is.null(status)) status <- 0L
  if(length(output)) message(paste(output,collapse="\n"))
  if(status!=0L && strict) stop("Panel alignment gate failed; inspect ",report_path)
  invisible(manifest)
}
