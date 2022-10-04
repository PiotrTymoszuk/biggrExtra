# Hyperdraw utils

# Generate a customized hypergraph from an SBML model ---------

#' Generate a hypergraph of the SMBL model.
#'
#' @description Produces a hypergraph (an instance of 'RagraphBPH' class).
#' The code was accommodated from \code{\link[BiGGR]{sbml2hyperdraw}} with
#' few modifications allowing for labeling of confidence intervals,
#' errors and significance.
#' @return a hypergraph object of the 'RagraphBPH' class. Plotted
#' conveniently with the plot() method.
#' @param model a SBML model.
#' @param rates a named numeric vector with reaction regulation estimates
#' (linear, not log2). The names should be reaction IDs compatible with the
#' SBML model.
#' @param rate_sep a character separating the reaction identifier
#' and its regulation rate in the graph.
#' @param suffixes additional information to be displayed in the graph along
#' with the reactions and their rates. Omitted if NULL.
#' @param suffix_sep a character separating the reaction identifier/rate
#' from the suffix in the in the graph.
#' @param colors a vector of colors for the reactions. Should be named with
#' reaction IDs.
#' @param node_color default node color.
#' @param signif_digits significant digits used for rounding of the
#' reaction rate.
#' @param fontsize size of the font, in cex-units.
#' @param relevant.species metabolites to be plotted. If not specified by
#' the user, all metabolites participating in the reactions specified by
#' 'relevant.reactions' are displayed.
#' @param relevant.reactions BiGG ID of the reactions to be plotted.
#' @param layoutType is a character string representing the layout engine to
#' be used for visualization. Current supported layouts are "dot", "twopi",
#' "neato","fdp","sfdp" and "circo". Defaults to "dot".
#' @param lwd.max a numeric given the maximum edge width. Defaults to 3.
#' @param lwd.min a numeric given the minimum edge width. Defaults to 0.5
#' @param plt.margins A numerical vector of the form
#' c(bottom, left, top, right) giving additional white space around the graph
#' (in case long node or edge labels fall outside the plotting region).
#' Defaults to c(150,150,150,150).

  sbml_to_hd <- function (model,
                          rates,
                          rate_sep = ':',
                          suffixes = NULL,
                          suffix_sep = ' ',
                          colors = NULL,
                          node_color = 'gray30',
                          signif_digits = 2,
                          fontsize = 0.75,
                          relevant.species = names(model@species),
                          relevant.reactions = names(model@reactions),
                          layoutType = 'dot',
                          lwd.max = 3,
                          lwd.min = 0.5,
                          plt.margins = c(150, 150, 150, 150)) {

    ## entry control -----

    if(!'Model' %in% class(model)) {

      stop('model has to be a valid SBML model.', call. = FALSE)

    }

    if(is.null(names(rates))) {

      stop('rates has to be a named numeric vector.', call. = FALSE)

    }

    if(all(!names(rates) %in% names(model@reactions))) {

      stop('no rates were found for the model reactions. Wrong rates vector?',
           call. = FALSE)

    }

    if(!is.null(suffixes)) {

      if(is.null(names(suffixes))) {

        stop('suffixes has to be a named character vector.',
             call. = FALSE)

      }

    }

    if(!is.null(colors)) {

      if(is.null(names(colors))) {

        stop('colors has to be a named character vector.',
             call. = FALSE)

      }

    }

    stopifnot(is.numeric(signif_digits))

    signif_digits <- as.integer(signif_digits)

    stopifnot(is.numeric(fontsize))
    stopifnot(is.numeric(lwd.max))
    stopifnot(is.numeric(lwd.min))
    stopifnot(is.numeric(plt.margins))

    miss_reactions <- relevant.reactions[!relevant.reactions %in% names(rates)]

    if(length(miss_reactions) > 0) {

      warning('Some of the relevant.reactions have no rates provided.
              They will be skipped from the graph.',
              call. = FALSE)

    }

    relevant.reactions <-
      relevant.reactions[relevant.reactions %in% names(rates)]

    ## defining the hyperedges, their labels and nodes -------

    react2edge <- function(r) {

      if (r@id %in% relevant.reactions) {

        reactants <- intersect(sapply(r@reactants, species),
                               relevant.species)

        products <- intersect(sapply(r@products, species),
                              relevant.species)

        if (length(reactants) > 0 & length(products) > 0) {

          my.label <- ifelse(r@id %in% names(rates),
                             paste(r@id,
                                   signif(rates[r@id], signif_digits),
                                   sep = rate_sep),
                             r@id)

          my.label <- ifelse(r@id %in% names(suffixes),
                             paste(my.label,
                                   suffixes[r@id],
                                   sep = suffix_sep),
                             my.label)

          if (!is.null(rates) && rates[r@id] < 0) {

            hypergraph::DirectedHyperedge(products,
                                          reactants,
                                          label = my.label)

          } else {

            hypergraph::DirectedHyperedge(reactants,
                                          products,
                                          label = my.label)

          }

        }
      }

    }

    hyperedges <- unlist(sapply(model@reactions,
                                react2edge,
                                simplify = TRUE))

    if (hasArg(rates)) {

      old_hyper_names <- names(hyperedges)
      old_rate_names <- names(rates)

      names(hyperedges) <-
        paste0(names(hyperedges),
               rate_sep,
               signif(rates, signif_digits)[names(hyperedges)])

      names(rates) <-
        paste0(names(rates),
               rate_sep,
               signif(rates, signif_digits))

      if(hasArg(suffixes)) {

        names(hyperedges) <- paste(names(hyperedges),
                                   suffixes[old_hyper_names],
                                   sep = suffix_sep)

        names(rates) <- paste(names(rates),
                              suffixes[old_rate_names],
                              sep = suffix_sep)

      }

    } else {

      rates <- sapply(model@reactions, function(r) 1)

    }

    node.names <- unique(unlist(c(lapply(hyperedges, function(x) x@head),
                                  lapply(hyperedges, function(x) x@tail))))

    ## base graph -----

    hg <- hypergraph::Hypergraph(node.names, hyperedges)

    testbph <- hyperdraw::graphBPH(hg)

    my.graph <- hyperdraw::graphLayout(testbph, layoutType = layoutType)

    nodeDataDefaults(my.graph, 'shape') <- 'box'
    nodeDataDefaults(my.graph, 'margin') <- 'unit(3, \"mm\")'
    nodeDataDefaults(my.graph, 'color') <- node_color
    nodeDataDefaults(my.graph, 'cex') <- fontsize

    edgeDataDefaults(my.graph, 'lwd') <- 1
    edgeDataDefaults(my.graph, 'color') <- node_color

    graphDataDefaults(my.graph, 'arrowLoc') <- 'end'

    ## customizing the graph -----

    ### scaling the line width - it should correspond to the regulation rate

    lwds <- abs(log2(rates[my.graph@edgeNodes]))

    lwds <- (lwds - min(lwds))/(max(lwds) - min(lwds))

    lwds <- ifelse(is.nan(lwds), 0, lwds)

    lwds <- lwd.min + lwds * lwd.max

    ## modifying the line width and color
    ## the edge color should correspond to the colors provided
    ## in a vector

    if(!is.null(colors)) {

      colors <- rlang::set_names(colors[old_rate_names],
                                 names(rates))

    }

    for(rxn.id in names(rates)) {

      if(!is.null(colors)) {

        lapply(my.graph@edgeNodeIO$outgoing[[rxn.id]],
               function(x) edgeData(my.graph, rxn.id, x, "color") <-
                 colors[[rxn.id]])

        lapply(my.graph@edgeNodeIO$incoming[[rxn.id]],
               function(x) edgeData(my.graph, x, rxn.id, "color") <-
                 colors[[rxn.id]])

      }

      lapply(my.graph@edgeNodeIO$outgoing[[rxn.id]],
             function(x) edgeData(my.graph, rxn.id, x, "lwd") <-
               as.character(lwds[rxn.id]))

      lapply(my.graph@edgeNodeIO$incoming[[rxn.id]],
             function(x) edgeData(my.graph, x, rxn.id, "lwd") <-
               as.character(lwds[rxn.id]))

    }

    my.graph@graph@boundBox@botLeft@y <-
      my.graph@graph@boundBox@botLeft@y - plt.margins[1]

    my.graph@graph@boundBox@botLeft@x <-
      my.graph@graph@boundBox@botLeft@x - plt.margins[2]

    my.graph@graph@boundBox@upRight@y <-
      my.graph@graph@boundBox@upRight@y + plt.margins[3]

    my.graph@graph@boundBox@upRight@x <-
      my.graph@graph@boundBox@upRight@x + plt.margins[4]


    return(my.graph)

  }

# END -----
