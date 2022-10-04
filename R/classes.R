# S3 class definitions

# geneSBML class ------

#' Create a geneSBML object.
#'
#' @description Creates a geneSBML object given a SBML model object and a
#' data frame with reaction regulation estimates and, optionally, regulation
#' errors, and a data frame with gene - reaction mapping.
#' @return an instance of geneSBML class.
#' @param model a SBML model.
#' @param reg a data frame with reaction IDs ('react_id'), regulation estimates
#' ('fold_reg') and, optionally, regulation errors ('error').
#' @param gene_map a data frame with reaction IDs ('react_id'),
#' list of Entrez ID identifiers ('entrez_id')
#' and evaluation expressions ('exprs').
#' @param mc optional, a matrix with reaction regulation estimates in each step
#' of Monte Carlo simulation. Columns are reactions, rows are subsequent runs.
#' @export

  geneSBML <- function(model, reg, gene_map, mc = NULL) {

    ## entry control

    if(!'Model' %in% class(model)) {

      stop('model has to be a valid SBML model.', call. = FALSE)

    }

    if(!is.data.frame(reg)) {

      stop("reg has to be a data frame with 'react_id' and 'fold_reg' columns.",
           call. = FALSE)

    }

    if(any(!c('react_id', 'fold_reg') %in% names(reg))) {

      stop("reg has to be a data frame with 'react_id' and 'fold_reg' columns.",
           call. = FALSE)

    }

    if(!is.data.frame(gene_map)) {

      stop("gene_map has to be a data frame with 'react_id', 'entrez_id' and 'exprs' columns.",
           call. = FALSE)

    }

    if(any(!c('react_id', 'entrez_id', 'exprs') %in% names(gene_map))) {

      stop("gene_map has to be a data frame with 'react_id', 'entrez_id' and 'exprs' columns.",
           call. = FALSE)

    }

    if(!is.null(mc)) {

      if(!is.numeric(mc)) {

        stop('mc has to be a numeric matrix.', call. = FALSE)

      }

      if(!is.matrix(mc)) {

        stop('mc has to be a numeric matrix.', call. = FALSE)

      }

    }

    ## object

    structure(list(model = model,
                   reg = reg,
                   gene_map = gene_map,
                   mc = mc),
              class = 'geneSBML')

  }

# END ------

