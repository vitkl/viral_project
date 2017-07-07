##' Retrieve molecular interactions for the random set of proteins
##' @name randomInteractome
##' @param taxid character (1L), taxonomy id of the species which interaction participants should belong to
##' @param n_inter integer (1L), the number of proteins for which to retrieve the random set of interactions
##' @param degree integer (\code{n_inter} L), specify the degree for \code{n_inter} number of proteins to produce the network with the specific degree distribution, if set to NULL (default) the degree distribution will correspond to that of \code{taxid} interactome
##' @param database argument for \code{\link{queryPSICQUIC}}, default is "imex"
##' @param seed integer (1L), set seed (\code{\link{set.seed}}) for random sampling (default is 1)
##' @param protein_only logical (1L), if TRUE the interaction participants are restricted to proteins (exclude other types of molecules such as RNA or small molecules)? By default is TRUE.
##' @details \code{taxid} is used to query specified database using PSICQUIC client, only interactions in which both participants belong the \code{taxid} are retured (\code{"taxidA:9606 AND taxidB:9606"}, not \code{"species:9606"})
##' @details Random network can be specified to have specific degree distribution. If the (\code{degree} parameter is set \code{taxid} proteins will be split by degree and from each degree group a sample of the size specified by how many times specific degree number is repeated in \code{degree} will be taken.
##' @details If the degree distribution is not specified a sample of \code{n_inter} is taken from all proteins which have interaction data available in the \code{database} for \code{taxid}. In this case, the degree distribution of the resulting set of proteins will be similar to the degree distribution in the interactome of \code{taxid} in \code{database}.
##' @import data.table
##' @export randomInteractome
randomInteractome = function(taxid = "9606", n_inter, degree = NULL, database = "imex", seed = 1, protein_only = TRUE){
  # if the interaction data for species taxid and from database is not saved in the library - queryPSICQUIC for interaction data for taxid interactions in the database and in MITAB2.5 format, save results to the library
  queryPSICQUIC(query = "id:P74565 AND detmethod:\"MI:0018\"",
                                format = "tab27",
                                database = "imex",
                                file = "P74565_2H_interactions_imex_tab27.tsv")
}
