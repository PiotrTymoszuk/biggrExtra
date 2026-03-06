# Extraction of features from Recon data frames with reaction annotation

# Mapping of reactions to gene identifiers -------

#' Extract gene identifiers from reaction annotations.
#'
#' @description
#' The functions extract Entrez ID gene identifiers mapped to all or
#' user-specified reactions present in a data frame with reaction annotation
#' and gene - reaction association rules.
#' Entrez IDs of genes associated with reactions are listed, and character
#' strings with the gene - reaction association rules are translated to R
#' expressions.
#' `as_reactDB()` offers a simplified tools for generation of
#' \code{\link{reactDB}} data frames,
#' while `extract_genes()` allows for selection of reactions of interest and
#' parsing error diagnostics.
#'
#' @details
#' Association rules with unrecognized Entrez IDs are removed from the output
#' and a parsing warning is raised.
#' The gene association rules have to operate with "bare" Entrez IDs without
#' the version information (i.e. numbers after a dot).
#'
#' @references
#' King ZA, Lu J, Dräger A, Miller P, Federowicz S, Lerman JA, Ebrahim A,
#' Palsson BO, Lewis NE. BiGG Models: A platform for integrating,
#' standardizing and sharing genome-scale models.
#' Nucleic Acids Res (2016) 44:D515–D522. doi:10.1093/NAR/GKV1049
#'
#' @return a data frame of class \code{\link{reactDB}} containing reaction
#' IDs (`id`), reaction names (`name`), subsystem information (`subsystem`),
#' list of gene identifiers (`entrez_id`) and a list of expressions to be
#' evaluated by at calculating reaction activity estimates.
#'
#' @param x a data frame with the following obligatory columns: `id` with reaction
#' identifiers, `name` with reaction names, `subsystem` with subsystem assignment,
#' and `gene_association` with character strings with gene association rules.
#' @param react_id BiGG reaction ID, with or without the leading 'R_' string.
#' Defaults to `NULL`, which means that all reactions are mapped to genes.
#' @param inspect_errors logical. If `TRUE`, the function returns parsing errors.
#' @param ... additional arguments passed to `extract_genes()`.
#'
#' @export

  extract_genes <- function(x,
                            react_id = NULL,
                            inspect_errors = FALSE) {

    ## entry control ----------

    if(!is.data.frame(x)) stop("`x` has to be a data frame.", call. = FALSE)

    fix_cols <- c("id", "name", "subsystem", "gene_association")

    missing_cols <- setdiff(fix_cols, names(x))

    if(length(missing_cols) > 0) {

      stop(paste("The following obligatory columns are missing from `x`:",
                 paste(missing_cols, collapse = ", ")),
           call. = FALSE)

    }

    x <- x[, fix_cols]

    class_check <- map_lgl(x, is.character)

    if(any(!class_check)) {

      stop(paste("The following obligatory column in `x` must be of character type:",
                 paste(fix_cols[!class_check]), collapse = ", "),
           call. = FALSE)

    }

    if(any(stri_detect(na.omit(x[["gene_association"]]),
                       regex = "\\d+\\.\\d+"))) {

      stop(paste("Gene association rules contain identifiers with",
                 "version information, i.e. number after a dot.",
                 "Please remove it to proceed."),
           call. = FALSE)

    }

    if(!is.null(react_id)) {

      stopifnot(is.character(react_id))

      react_id <-
        ifelse(!stri_detect(react_id, regex = '^R_'),
               paste0('R_', react_id), react_id)

      x <- filter(x, .data[["id"]] %in% react_id)

      if(nrow(x) == 0) {

        stop("No reactions to precess after filtering with `react_id`.",
             call. = FALSE)

      }

    }

    stopifnot(is.logical(inspect_errors))
    inspect_errors <- inspect_errors[1]

    ## filtering the reaction annotation df ---------

    ### getting rid of empty rules and space-only rules in the annotation
    ### data frame

    proc_df <-
      filter(x[, c("id", "gene_association")],
             !is.na(.data[["id"]]),
             .data[["id"]] != "",
             !is.na(.data[["gene_association"]]),
             .data[["gene_association"]] != "",
             !stri_detect(.data[["gene_association"]],
                          regex = "^\\s+$"))

    ### removal of leading and trailing spaces in the gene association rules

    proc_df[["gene_association"]] <-
      stri_replace_all(proc_df[["gene_association"]],
                       regex = "^\\s+",
                       replacement = "")

    proc_df[["gene_association"]] <-
      stri_replace_all(proc_df[["gene_association"]],
                       regex = "\\s+$",
                       replacement = "")

    ## extraction of the Entrez IDs --------

    proc_df[["entrez_id"]] <-
      map(proc_df[["gene_association"]],
          stri_extract_all, regex = "\\d+")

    proc_df[["entrez_id"]] <- map(proc_df[["entrez_id"]],
                                  unlist)

    proc_df[["entrez_id"]] <- map(proc_df[["entrez_id"]],
                                  unique)

    ## translation of gene assignment rules to R expressions ------

    exp_lst <-
      map(set_names(proc_df[["gene_association"]],
                    proc_df[["id"]]),
          escape_numbers)

    exp_lst <-
      map(exp_lst,
          stri_replace_all,
          regex = "and|AND",
          replacement = "%AND%")

    exp_lst <-
      map(exp_lst,
          stri_replace_all,
          regex = "or|OR",
          replacement = "%OR%")

    exp_lst <- map(exp_lst, safely(str2lang))

    parse_errors <- compact(map(exp_lst, ~.x$error))

    if(length(parse_errors) > 0) {

      warning(paste('There were', length(parse_errors), 'parsing errors.'),
              call. = FALSE)

    }

    if(inspect_errors) return(parse_errors)

    exp_lst <- compact(map(exp_lst, ~.x$result))

    ## output -------

    id <- NULL
    exprs <- NULL

    map_df <- tibble(id = names(exp_lst),
                     exprs = exp_lst)

    out_df <- left_join(proc_df, map_df, by = "id")

    out_tbl <- left_join(x[, c("id", "name", "subsystem")],
                         out_df,
                         by = "id")

    return(reactDB(out_tbl))

  }

#' @rdname extract_genes
#' @export

  as_reactDB <- function(x, ...) extract_genes(x, ...)

# END ------
