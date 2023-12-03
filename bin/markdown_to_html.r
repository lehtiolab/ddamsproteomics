#!/usr/bin/env Rscript

# Command line argument processing
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: markdown_to_html.r <input.md> <output.html>", call.=FALSE)
}
markdown_fn <- args[1]
output_fn <- args[2]

# Load / install packages
if (!require("markdown")) {
  install.packages("markdown", dependencies=TRUE, repos='http://cloud.r-project.org/')
  library("markdown")
}

mark_html(
  file = markdown_fn,
  output = output_fn,
  options = c('toc', 'base64_images', 'highlight_code')
)
