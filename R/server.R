#' StartHttpServer
#'
#' Create HTTP server
#' @param server.host Server host
#' @param server.port Server port
#' @return Returns HTTP SERVER instance
#'
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom httpuv startDaemonizedServer
#'
#' @export
#'
StartHttpServer <- function(server.host = "127.0.0.1", server.port = 9009) {
  options(encoding="UTF-8")
  options(verbose=T)

  is.defined <- function(sym) {
    sym <- deparse(substitute(sym))
    env <- parent.frame()
    exists(sym, env)
  }

  app <- list(
    onWSOpen = function(ws) {
      cat("Opened conenction with WS client => Start to receive message\n")
      ws$onMessage(function(binary, message) {
        if(message == "stop_server") {
          StopHttpServer()
          return(NA)
        }
        tryCatch({
          resultData <- fromJSON('{}')
          message <- fromJSON(message)
          response <- eval(parse(text=message$command))
          resultData$code <- 0
          resultData$result <- response
          resultData$commandID <- message$commandID
          ws$send(jsonlite::toJSON(resultData, auto_unbox=TRUE))
        }, error = function(error_message) {
          resultData$code <- 1
          resultData$result$error <- error_message[[1]]
          if(is.object(message) == FALSE) {
            if(is.defined(message$commandID) == TRUE) {
              resultData$commandID <- message$commandID
            }
          }
          ws$send(jsonlite::toJSON(resultData, auto_unbox=TRUE))
        })
      })
    }
  )

  server <- httpuv::runServer(server.host, server.port, app)
  return(server)
}

#' StopHttpServer
#'
#' Stop HTTP server
#' @param server.host Server host
#' @param server.port Server port
#' @return Returns HTTP SERVER instance
#'
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom httpuv startDaemonizedServer
#'
#' @export
#'
StopHttpServer <- function() {
  tryCatch({
    httpuv::stopAllServers();
    quit(save = "no")
  }, error = function(error_message) {
    print(error_message)
  })
}
