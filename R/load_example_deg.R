#' #' Load Example Dataset for Differential Expression Analysis
#' #'
#' #' @param id A dataset identifier
#' #' @param id i the compound number in GSE85871
#' #' @return a `list`.
#' #' @export
#' #'
#' #' @examples
#' #' deg <- load_example_deg(id = "GSE85871", i)
#' #' deg.all <- lapply(c(1:30, 32:64, 66:105), function(x) load_example_deg(id = "GSE85871", x))
#' load_example_deg <- function(id = "GSE85871", i) {
#'   if (id == "GSE85871") {
#'     data <- load_example_dataset()
#'     group_list_all <- stringr::str_split(as.character(data$pdata$title), "_", simplify = T)[, 2]
#'     if (i %in% 1:30) {
#'       sel <- c(2 * i - 1, 2 * i, 61, 62)
#'     } else if (i %in% 32:64) {
#'       sel <- c(2 * i - 1, 2 * i, 129, 130)
#'     } else if (i %in% 66:105) {
#'       sel <- c(2 * i - 1, 2 * i, 211, 212)
#'     }
#'     exp <- data$expr[, sel]
#'     group_list <- group_list_all[sel]
#'     message("Now calculate is ", group_list[1], "  VS  ", group_list[3])
#'     result <- deg.cal(exp, group = group_list, level = group_list[c(3, 1)], method = "limma")
#'     return(result)
#'     message("Done.")
#'   }
#' }
