# Utilities

# Escape numerics in a string -----

#' Escape Entrez ID numbers in a string.
#'
#' @description
#' Finds Entrez ID in a string and escapes them with '``.
#'
#' @return a sting.
#'
#' @param str a string.

  escape_numbers <- function(str) {

    stopifnot(is.character(str))

    numbs <- sort(unique(unlist(stri_extract_all(str, regex = "\\d+"))))

    for(i in numbs) {

      search_rex <- paste0("(?<!\\d{1})", i, "(?!\\d+)")

      str <-
        stri_replace_all(str,
                         regex = search_rex,
                         replacement = paste0(" `", i, "` "))

    }

    return(str)

  }

# END ------
