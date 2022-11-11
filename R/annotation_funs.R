# Functions used for gene, metabolite and reaction annotation

# Annotate metabolite or reaction ------

#' Retrieve metabolite or reaction features from a local BiGG database.
#'
#' @description Retrieves the requested metabolite or reaction name from
#' a local database given the BiGG ID as a key.
#' @details You may use the 'metabolites' or 'reactions' data sets provided
#' by the package or download the latest annotation data sets
#' from http://bigg.ucsd.edu/data_access.
#' @param bigg_id BiGG ID identifier of reactions or metabolites.
#' A character vector, with or without trailing 'R_' or 'M_' strings.
#' @param value the name of the annotation table output variable.
#' @param annotation_db annotation data frame or a SBML database.
#' @param id_type type of the identifier: reaction (default) or metabolite.
#' Ignored if annotation_db is a data frame.
#' @return a vector with the requested features.
#' @export

  annotate_bigg <- function(bigg_id,
                            value = 'name',
                            annotation_db = reactions,
                            id_type = c('reaction', 'metabolite')) {

    ## entry control

    if(!is.character(bigg_id)) {

      stop('bigg_id has to be a character vector.', call. = FALSE)

    }

    if(!is.data.frame(annotation_db)) {

      if(!any(c('SBML', 'geneSBML') %in% class(annotation_db))) {

        stop('annotation_db hast to be a data frame, a SBML database or a geneSBML object.',
             call. = FALSE)

      }

    }

    if(is.data.frame(annotation_db)) {

      if(!value %in% names(annotation_db)) {

        stop('value absent from annotation_db', call. = FALSE)

      }

      new_id <- stringi::stri_replace(bigg_id,
                                      regex = '(R_)|(M_)',
                                      replacement = '')

      annot_vec <- rlang::set_names(annotation_db[[value]],
                                    annotation_db$bigg_id)

      return(rlang::set_names(annot_vec[new_id],
                              bigg_id))

    } else {

      if(is_geneSBML(annotation_db)) {

        annotation_db <- annotation_db$model

      } else {

        annotation_db <- annotation_db@model

      }

      id_type <- match.arg(id_type[1], c('reaction', 'metabolite'))

      class_type <- switch(id_type,
                           reaction = 'Reaction',
                           metabolite = 'Species')

      slot_name <- switch(id_type,
                          reaction = 'reactions',
                          metabolite = 'species')

      all_slots <- names(getSlots(class_type))

      if(!value %in% all_slots) {

        stop('value absent from annotation_db', call. = FALSE)

      }

      annot_vec <- names(slot(annotation_db, slot_name))

      annot_vec <- annot_vec[annot_vec %in% bigg_id]

      val_vec <- slot(annotation_db, slot_name)[annot_vec]

      val_vec <- purrr::map(val_vec, ~slot(.x, value))

      return(val_vec)

    }


  }

# Check if reactions are present in the the geneSBML model ------

#' Check if reactions or metabolites are present in the geneSBML model -------
#'
#' @description Checks if the given reaction can be found in the geneSBML model.
#' @return a list of logical vectors logical vectors, each for the reactions
#' and metabolite IDs, if both the reactions and metabolites are queried.
#' Otherwise a named vector with logical values.
#' @param react_id BiGG reaction ID, with or without the leading 'R_' string.
#' @param metab_id BiGG metabolite ID, with or without the leading 'M_' string.
#' @param geneSBML a geneSBML object.
#' @export

  check_geneSBML <- function(react_id = NULL, metab_id = NULL, geneSBML) {

    ## entry control

    if(!is_geneSBML(geneSBML)) {

      stop('A geneSBML object is required.', call. = FALSE)

    }

    if(is.null(react_id) & is.null(metab_id)) {

      stop('At least one is required: react_id or metab_id vector.',
           call. = FALSE)

    }

    ## checking

    if(!is.null(react_id)) {

      react_id <-
        purrr::map(react_id, function(x) if(!stringi::stri_detect(x, regex = '^R_')) paste0('R_', x) else x)

      react_tst <-
        rlang::set_names(unlist(react_id) %in% geneSBML$reg$react_id,
                         unlist(react_id))

    } else {

      react_tst <- NULL

    }

    if(!is.null(metab_id)) {

      metab_id <-
        purrr::map(metab_id, function(x) if(!stringi::stri_detect(x, regex = '^M_')) paste0('M_', x) else x)

      metab_tst <-
        rlang::set_names(unlist(metab_id) %in% names(geneSBML$model@species),
                         unlist(metab_id))

    } else {

      metab_tst <- NULL

    }

    res <- purrr::compact(list(react_id = react_tst,
                               metab_id = metab_tst))

    if(length(res) == 1) res <- res[[1]]

    res

  }

# Map reactions to genes ------

#' Map reactions to gene EntrezID.
#'
#' @description Maps reaction BiGG IDs to the Entrez IDs of the associated
#' genes.
#' @param react_id BiGG reaction ID, with or without the leading 'R_' string.
#' @param geneSBML a geneSBML object.
#' @return a data frame with reaction ids and a list of the associated
#' gene Entrez IDs.
#' @export

  react_to_gene <- function(react_id, geneSBML) {

    ## entry control

    if(!is_geneSBML(geneSBML)) {

      stop('A geneSBML object is required.', call. = FALSE)

    }

    react_id <-
      ifelse(stringi::stri_detect(react_id, regex = '^R_'),
             react_id,
             paste0('R_', react_id))

    ## mapping

    genes <- dplyr::filter(components(geneSBML, type = 'gene_map'),
                           .data[['react_id']] %in% .env[['react_id']])

    genes <- dplyr::mutate(genes,
                           entrez_id = purrr::map(entrez_id, unique))

    genes[c('react_id', 'entrez_id')]

  }

# Map genes to reactions --------

#' Map gene EntrezIDs to reactions.
#'
#' @description Maps gene Entrez IDs to the associated reactions.
#' @param entrez_id a vector of of Entrez IDs.
#' @param geneSBML a geneSBML object.
#' @param lead_str logical, should the reaction ID contain the leading 'R_'
#' string? Defaults to TRUE.
#' @return a data frame with Entrez IDs and a list of the associated
#' reaction BiGG IDs.
#' @export

  gene_to_react <- function(entrez_id, geneSBML, lead_str = TRUE) {

    ## entry control

    if(!is_geneSBML(geneSBML)) {

      stop('A geneSBML object is required.', call. = FALSE)

    }

    stopifnot(is.logical(lead_str))

    ## mapping

    react_map <- components(geneSBML, 'gene_map')

    sel_lst <- purrr::map(entrez_id,
                          function(gene_id) purrr::map_lgl(react_map[['entrez_id']],
                                                           ~any(.x == gene_id)))

    sel_lst <- rlang::set_names(sel_lst,
                                entrez_id)

    reacts <- purrr::map(sel_lst,
                         ~react_map[['react_id']][.x])

    react_tbl <- tibble::tibble(entrez_id = entrez_id,
                                react_id = reacts)

    if(!lead_str) {

      react_tbl <- dplyr::mutate(react_tbl,
                                 react_id = purrr::map(react_id,
                                                       stringi::stri_replace,
                                                       regex = '^R_',
                                                       replacement = ''))

    }

    react_tbl

  }

# Map reactions to metabolites -------

#' Map reactions to metabolites.
#'
#' @description Maps reaction BiGG IDs to BiGG IDs of
#' the associated metabolites.
#' @param react_id BiGG reaction ID, with or without the leading 'R_' string.
#' @param annotation_db an annotation data frame with
#' the 'bigg_id' and 'reaction_string' variables.
#' @param lead_str logical, should the metabolite ID contain the leading 'M_'
#' string? Defaults to TRUE.
#' @param exc_regex a regular expression to exclude
#' some of metabolites (such as water) from the output. Ignored if NULL.
#' @param detailed logical, should the output contain information on substrates
#' and products? Dafaluts to FALSE.
#' @return a data frame with reaction IDs and a list
#' of the associated metabolite IDs. If 'detailed' is set to TRUE,
#' the colums 'lhs' and 'rhs' contain the left and right hand side metabolites
#' of the reaction.
#' @export

  react_to_metab <- function(react_id,
                             annotation_db = reactions,
                             lead_str = TRUE,
                             exc_regex = NULL,
                             detailed = FALSE) {

    ## entry control ------

    if(!is.data.frame(annotation_db)) {

      stop('annotation_db should be a data frame with bigg_id and reaction_string variables.',
           call. = FALSE)

    }

    stopifnot(is.logical(lead_str))
    stopifnot(is.logical(detailed))

    ## metabolite database handling -------

    metabs <- annotate_bigg(react_id,
                            value = 'reaction_string',
                            annotation_db = annotation_db)

    metabs <- rlang::set_names(metabs,
                               react_id)

    ## simple output -------

    if(!detailed) {

      metabs <- purrr::map(metabs,
                           stringi::stri_split,
                           regex = '(\\s{1}<->\\s{1})|(\\s{1}\\+\\s{1})')

      metabs <- purrr::map(metabs,
                           unlist)

      metabs <- purrr::map(metabs,
                           unique)

      metabs <- purrr::map(metabs,
                           stringi::stri_replace_all,
                           regex = '^\\d+\\.\\d+\\s{1}',
                           replacement = '')

      metabs <- purrr::map(metabs,
                           stringi::stri_replace_all,
                           fixed = '__',
                           replacement = '_')

      metabs <- purrr::map(metabs,
                           unique)

      metabs <- purrr::map(metabs,
                           ~.x[!is.na(.x)])

      if(!is.null(exc_regex)) {

        metabs <- purrr::map(metabs,
                             ~.x[!stringi::stri_detect(.x, regex = exc_regex)])

      }

      if(lead_str) {

        metabs <- purrr::map(metabs,
                             ~ifelse(stringi::stri_detect(.x, regex = '^M_'),
                                     .x,
                                     paste0('M_', .x)))

      }

      return(tibble::tibble(react_id = names(metabs),
                            metab_id = metabs))

    }

    ## substrates and metabolites -----

    metabs <- purrr::map(metabs,
                         stringi::stri_split,
                         fixed = ' <-> ')

    metab_tbl <-
      tibble::tibble(react_id = names(metabs),
                     lhs = purrr::map_chr(metabs, ~.x[[1]][1]),
                     rhs = purrr::map_chr(metabs, ~.x[[1]][2]))

    metab_tbl <-
      purrr::map_dfc(metab_tbl,
                     stringi::stri_replace_all,
                     regex = '\\d+\\.\\d+\\s{1}',
                     replacement = '')

    metab_tbl <-
      purrr::map_dfc(metab_tbl,
                     stringi::stri_replace_all,
                     fixed = '__',
                     replacement = '_')

    metab_tbl <-
      dplyr::mutate(metab_tbl,
                    lhs = stringi::stri_split(lhs,
                                              fixed = ' + '),
                    rhs = stringi::stri_split(rhs,
                                              fixed = ' + '))

    return(metab_tbl)

  }

# END ------
