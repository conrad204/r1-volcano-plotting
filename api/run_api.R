# run_api.R — entrypoint for the volcano-plot plumber API
library(plumber)

# Plumb the API definition (path relative to WORKDIR, i.e. /app)
pr <- plumb("api/plumber.R")

# Serve the static frontend from web/
pr_static(pr, "/", "web")

port <- as.integer(Sys.getenv("PORT", "8000"))
host <- Sys.getenv("HOST", "0.0.0.0")  # CRITICAL for Docker
message(sprintf("Starting plumber on %s:%d  (static files from web/)", host, port))
pr$run(host = host, port = port)