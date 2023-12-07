# Extraction of features from SBML Document class databases.

# Mapping of reactions to gene identifiers -------

#' Extract gene identifiers from a database.
#'
#' @description
#' Extracts Entrez ID gene identifiers mapped to all or
#' user-specified reactions present in the given database SBML database.
#'
#' @details
#' A wrapper around \code{\link[BiGGR]{extractGeneAssociations}}.
#' Entrez ID identifier version is discarded silently. Reactions with
#' non-Entrez ID features mapped are removed from the output and
#' a parsing warning is raised.
#'
#' @return a data frame containing reaction ID ('react_id'),
#' list of gene identifiers ('entrez_id') and a list of expressions to be
#' evaluated by at calculating reaction regulation estimates.
#'
#' @param database an object of class SBMLDocument.
#'
#' @param react_id BiGG reaction ID, with or without the leading 'R_' string.
#' Defaults to NULL, which means that all reactions are mapped to genes.
#'
#' @export

  extract_genes <- function(database, react_id = NULL) {

    ## the entry control is done in part by the BiGGR function -------

    gene_lst <- extractGeneAssociations(database)

    if(!is.null(react_id)) {

      react_id <-
        ifelse(!stri_detect(react_id, regex = '^R_'),
               paste0('R_', react_id), react_id)

    }

    ## getting rid of empty strings ----

    non_empty_rules <- map_lgl(gene_lst, ~!stri_detect(.x, regex = '^\\s{1}$'))

    gene_lst <- gene_lst[non_empty_rules]

    ## removing the leading spaces -------

    gene_lst <- map(gene_lst,
                    stri_replace_all,
                    regex = '^\\s{1}',
                    replacement = '')

    ## extraction of the Entrez IDs --------

    entrez_lst <- map(gene_lst,
                      stri_extract_all,
                      regex = '\\d+\\.\\d+')

    entrez_lst <- map(entrez_lst, unlist)

    entrez_lst <- map(entrez_lst,
                      stri_replace,
                      regex = '\\.\\d+',
                      replacement = '')

    ## generating an expression list ------
    ## warning at possible parsing errors,
    ## escaping the numbers gene entry versions

    exp_lst <- map(gene_lst,
                   stri_replace_all,
                   fixed = 'and',
                   replacement = '%AND%')

    exp_lst <- map(exp_lst,
                   stri_replace_all,
                   fixed = 'or',
                   replacement = '%OR%')

    exp_lst <- map(exp_lst, escape_numbers)

    exp_lst <- map(exp_lst, safely(str2lang))

    parse_errors <- compact(map(exp_lst, ~.x$error))

    if(length(parse_errors) > 0) {

      warning(paste('There were', length(parse_errors), 'parsing errors.'),
              call. = FALSE)

    }

    exp_lst <- compact(map(exp_lst, ~.x$result))

    ## output -------

    map_tbl <- tibble(react_id = names(exp_lst),
                      entrez_id = entrez_lst[names(exp_lst)],
                      exprs = exp_lst)

    if(is.null(react_id)) {

      return(map_tbl)

    }

    filter(map_tbl, .data[['react_id']] %in% .env[['react_id']])

  }

# Mapping reactions to subsystem ------

#' Extract subsystem assignment from a database.
#'
#' @description
#' Extracts mapping of all or user-specified reactions to
#' subsystems defines in the given database SBML database or geneSBML model.
#'
#' @return if `as_list` is set to FALSE a data frame with the reaction ID
#' ('react_id') and subsystem variable ('subsystem').
#' Otherwise a subsystem-named list of reaction ID vectors.
#'
#' @param database an object of class SBMLDocument or class geneSBML.
#' @param react_id BiGG reaction ID, with or without the leading 'R_' string.
#' Defaults to NULL, which means that all reactions are mapped.
#' @param as_list logical, should a list be returned? Defaults to FALSE.
#'
#' @export

  extract_subsystems <- function(database,
                                 react_id = NULL,
                                 as_list = FALSE) {

    ## entry control -----

    if(!is_SBML(database)) {

      if(!is_geneSBML(database)) {

        stop(paste("'database' has to be a valid SBML Document or",
                   "a geneSBML class object."),
             call. = FALSE)

      }

    }

    if(!is.null(react_id)) {

      react_id <-
        ifelse(!stri_detect(react_id, regex = '^R_'),
               paste0('R_', react_id), react_id)

    }

    ## mapping -------

    if(is_geneSBML(database)) {

      reaction_lst <- database$model@reactions

    } else {

      reaction_lst <- database@model@reactions

    }

    model_notes <- map(reaction_lst, ~.x@notes)

    if(!is.null(react_id)) {

      model_notes <- model_notes[react_id]

      if(length(model_notes) == 0) {

        warning('No subsystems for the specified reaction were found.',
                call. = FALSE)

        return(NULL)

      }

    }

    model_notes <- map(model_notes,
                       stri_extract,
                       regex = '<p>SUBSYSTEM:\\s{1}.*</p>')

    model_notes <- map(model_notes,
                       stri_extract,
                       regex = '\\s{1}.*</p>$')

    model_notes <- map(model_notes,
                       stri_replace_all,
                       regex = '(^\\s{1})|(</p>$)',
                       replacement = '')

    subsystem <- NULL

    sub_tbl <- map2_dfr(model_notes, names(model_notes),
                        ~tibble(react_id = .y,
                                subsystem = .x))

    if(!as_list) {

      return(sub_tbl)

    }

    split(sub_tbl$react_id, sub_tbl[[subsystem]])

  }

# END ------
