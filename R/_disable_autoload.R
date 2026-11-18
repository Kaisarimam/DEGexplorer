# The presence of this file stops shiny::runApp() from sourcing every file in R/
# in alphabetical order before app.R. app.R sources the modules itself, in the
# order their top-level definitions need (organisms.R before annotation_sources.R,
# for example), so Shiny's autoloading would both run them twice and run them in
# the wrong order.
