options(rsconnect.http.trace = TRUE)
options(rsconnect.http.verbose = TRUE)

source("ui.R")
source("server.R")
options(rsconnect.max.bundle.size = 1 * 1024^3)  # 1 GB

shinyApp(ui, server)
