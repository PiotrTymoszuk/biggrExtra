# Functions used for gene, metabolite and reaction annotation

# Annotate metabolite or reaction ------

#' Retrieve metabolite or reaction features from a local BiGG database.
#'
#' @description Retrieves the requested metabolite or reaction feature from
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
#' @param react_id BiGG reaction ID, with or without the leading 'R_' string
#' @param metab_id BiGG metabolite ID, with or without the leading 'M_' string
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
