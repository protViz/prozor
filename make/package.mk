# Build the documentation site with altdoc instead of pkgdown.
SITE_CMD := Rscript -e "altdoc::render_docs(freeze = FALSE)"
