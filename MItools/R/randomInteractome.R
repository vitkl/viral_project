##' Retrieve molecular interactions for the random set of proteins (of a particular taxon)
##' @name randomInteractome
##' @author Vitalii Kleshchevnikov
##' @param MITABdata data.table containing pre-loaded molecular interaction data as returned by \code{\link{queryPSICQUICrlib}}, usefull for taking multiple samples, the default in NULL
##' @param degree_data data.table containing pre-calculated (using \code{\link{edgelist2degree}}) degree for each node in MITABdata, usefull for taking multiple samples, the default in NULL
##' @param n_inter integer (1L), the number of proteins for which to retrieve the random set of interactions
##' @param degree integer (\code{n_inter} L), specify the degree for \code{n_inter} number of proteins to produce the network with the specific degree distribution, if set to NULL (default) the degree distribution will correspond to that of \code{taxid} interactome
##' @details Random network can be specified to have specific degree distribution. If the (\code{degree} parameter is set \code{taxid} proteins will be split by degree and from each degree group a sample of the size specified by how many times specific degree number is repeated in \code{degree} will be taken.
##' @details If the degree distribution is not specified a sample of \code{n_inter} is taken from all proteins which have interaction data available in the \code{database} for \code{taxid}. In this case, the degree distribution of the resulting set of proteins will be similar to the degree distribution in the interactome of \code{taxid} in \code{database}.
##' @details \code{randomInteractome} retrieves molecular interactions using \code{\link{fullInteractome}}
##' @import data.table
##' @export randomInteractome
##' @examples
##' # retrive the interactome using PSICQIUC servise (or by reading local copy) from IMEx databases for a list of 200 random human (9606) proteins, not specifying their degree distribution
##' set.seed(1)
##' random = randomInteractome(n_inter = 200, degree = NULL, taxid = "9606", database = "imex", protein_only = TRUE)
##' # retrive the interactome from MITABdata for a list of 200 random human (9606) proteins, not specifying their degree distribution
##' full = fullInteractome(taxid = "9606", database = "imex", format = "tab25", clean = TRUE, protein_only = TRUE)
##' set.seed(1)
##' random = randomInteractome(MITABdata = full, n_inter = 200, degree = NULL)
randomInteractome = function(MITABdata = NULL, degree_data = NULL, n_inter, degree = NULL, ...){
  # if MITABdata is NULL retrive the cleaned full Interactome
  if(is.null(MITABdata)) full_interactome_clean = fullInteractome(..., format = "tab25", clean = TRUE)
  # if MITABdata is supplied use that data
  if(!is.null(MITABdata)) full_interactome_clean = MITABdata

  # check if the data has necessary columns:
  if(mean(c("IDs_interactor_A","IDs_interactor_B", "pair_id") %in% colnames(full_interactome_clean)) != 1) stop("MITABdata is supplied in the wrong format")

  # get interactors
  interactors = full_interactome_clean[, unique(c(IDs_interactor_A, IDs_interactor_B))]

  # if the degree distribution is not specified - sample n_inter of interactors
  if(is.null(degree)){
    random_interactors = sample(interactors, n_inter)
  }

  # if the degree distribution is specified - calculate degree of each interactor (the number of interacting partners per interactor)
  if(!is.null(degree) & (length(degree) == n_inter)){
    # calculate degree data if not provided
    if(is.null(degree_data)) degree_data = edgelist2degree(full_interactome_clean[,.(pair_id)])
    # stop if degree_data doesn't match MITABdata
    if(!is.null(degree_data) & mean(degree_data$ID %in% interactors) != 1) stop("degree_data doesn't contain degree information for all nodes in MITABdata")
    # find the degree distribution
    degree_frequency = ecdf(degree)
    unique_degrees = unique(degree_data[,.(N)])
    unique_degrees[, ecdf_freq := degree_frequency(N[1:nrow(unique_degrees)]) - degree_frequency(c(0,N[1:(nrow(unique_degrees)-1)]))]
    degree_data = degree_data[unique_degrees, on = "N"]
    random_interactors = sample(x = degree_data[ecdf_freq != 0, ID],
                                size = n_inter,
                                prob = degree_data[ecdf_freq != 0, ecdf_freq])
  }
  # if the degree distribution is specified but doesn't match the number of interactions n_inter
  if(!is.null(degree) & !(length(degree) == n_inter)) stop(paste0("number of interactions n_inter doesn't match the length of the degree vector, N: ",n_inter," while the length of degree vector: ", length(degree)))

  # retrive interactions for a set of random proteins
  interactors2interactions(full_interactome_clean, random_interactors)

}
