# Extraction of features from SBML Document class databases.

# Mapping of reactions to gene identifiers -------

#' Extract gene identifiers from a database.
#'
#' @description Extracts Entrez ID gene identifiers mapped to all or
#' user-specified reactions present in the given database SBML database.
#' @details A wrapper around \code{\link[BiGGR]{extractGeneAssociations}}.
#' Entrez ID identifier version is discarded silently. Reactions with
#' non-Entrez ID features mapped are removed from the output and
#' a parsing warning is raised.
#' @return a data frame containing reaction ID ('react_id'),
#' list of gene identifiers ('entrez_id') and a list of expressions to be
#' evaluated by at calculating reaction regulation estimates.
#' @param database an object of class SBMLDocument.
#' @param react_id BiGG reaction ID, with or without the leading 'R_' string.
#' Defaults to NULL, which means that all reactions are mapped to genes.
#' @export

  extract_genes <- function(database, react_id = NULL) {

    ## the entry control is done in part by the BiGGR function

    gene_lst <- extractGeneAssociations(database)

    if(!is.null(react_id)) {

      react_id <-
        ifelse(!stringi::stri_detect(react_id, regex = '^R_'),
               paste0('R_', react_id), react_id)

    }

    ## getting rid of empty strings and

    gene_lst <-
      purrr::map(gene_lst,
                 function(x) if(stringi::stri_detect_regex(x, pattern = '^\\s{1}$')) NULL else x)

    gene_lst <- purrr::compact(gene_lst)

    ## removing the leading spaces

    gene_lst <- purrr::map(gene_lst,
                           stringi::stri_replace_all_regex,
                           pattern = '^\\s{1}',
                           replacement = '')

    ## extraction of the Entrez IDs

    entrez_lst <- purrr::map(gene_lst,
                             stringi::stri_extract_all_regex,
                             pattern = '\\d+\\.\\d+')

    entrez_lst <- purrr::map(entrez_lst, unlist)

    entrez_lst <- purrr::map(entrez_lst,
                             stringi::stri_replace_all_regex,
                             pattern = '\\.\\d+',
                             replacement = '')

    ## generating an expression list, warning at possible parsing errors
    ## escaping the numbers

    exp_lst <- purrr::map(gene_lst,
                          stringi::stri_replace_all_fixed,
                          pattern = 'and',
                          replacement = '%AND%')

    exp_lst <- purrr::map(exp_lst,
                          stringi::stri_replace_all_fixed,
                          pattern = 'or',
                          replacement = '%OR%')

    exp_lst <- purrr::map(exp_lst,
                          escape_numbers)

    exp_lst <- purrr::map(exp_lst,
                          purrr::safely(str2lang))

    parse_errors <- purrr::map(exp_lst, ~.x$error)

    parse_errors <- purrr::compact(parse_errors)

    if(length(parse_errors) > 0) {

      warning(paste('There were', length(parse_errors), 'parsing errors.'),
              call. = FALSE)

    }

    exp_lst <- purrr::map(exp_lst, ~.x$result)

    exp_lst <- purrr::compact(exp_lst)

    ## output

    map_tbl <- tibble::tibble(react_id = names(exp_lst),
                              entrez_id = entrez_lst[names(exp_lst)],
                              exprs = exp_lst)

    if(is.null(react_id)) {

      return(map_tbl)

    }

    dplyr::filter(map_tbl,
                  .data[['react_id']] %in% .env[['react_id']])

  }

# Mapping reactions to subsystem ------

#' Extract subsystem assignment from a database.
#'
#' @description Extracts mapping of all or user-specified reactions to
#' subsystems defines in the given database SBML database.
#' @return if 'as_list' is set to FALSE a data frame with the reaction ID
#' ('react_id') and subsystem variable ('subsystem').
#' Otherwise a subsystem-named list of reaction ID vectors.
#' @param database an object of class SBMLDocument.
#' @param react_id BiGG reaction ID, with or without the leading 'R_' string.
#' Defaults to NULL, which means that all reactions are mapped.
#' @param as_list logical, should a list be returned? Defaults to FALSE.
#' @export

  extract_subsystems <- function(database,
                                 react_id = NULL,
                                 as_list = FALSE) {

    ## entry control -----

    if(!'SBML' %in% class(database)) {

      stop("'database' has to be a valid SBML Document class object.",
           call. = FALSE)

    }

    if(!is.null(react_id)) {

      react_id <-
        ifelse(!stringi::stri_detect(react_id, regex = '^R_'),
               paste0('R_', react_id), react_id)

    }

    ## mapping -------

    model_notes <- purrr::map(database@model@reactions, ~.x@notes)

    if(!is.null(react_id)) {

      model_notes <- model_notes[react_id]

      if(length(model_notes) == 0) {

        warning('No subsystems for the specified reaction were found.',
                call. = FALSE)

        return(NULL)

      }

    }

    model_notes <- purrr::map(model_notes,
                              stringi::stri_extract,
                              regex = '<p>SUBSYSTEM:\\s{1}.*</p>')

    model_notes <- purrr::map(model_notes,
                              stringi::stri_extract,
                              regex = '\\s{1}.*</p>$')

    model_notes <- purrr::map(model_notes,
                              stringi::stri_replace_all,
                              regex = '(^\\s{1})|(</p>$)',
                              replacement = '')

    sub_tbl <- purrr::map2_dfr(model_notes, names(model_notes),
                               ~tibble::tibble(react_id = .y,
                                               subsystem = .x))

    if(!as_list) {

      return(sub_tbl)

    }

    plyr::dlply(sub_tbl, 'subsystem', function(x) x$react_id)

  }
