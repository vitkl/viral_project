##' \code{cleanMITAB} extracts interactor Uniprot (or else, specified by interactor_IDs_databases) IDs, interactor taxonomy IDs, Publication Identifiers, Confidence values and generates unique pair ID
##' @name cleanMITAB
##' @author Vitalii Kleshchevnikov
##' @param mitab data.table containing molecular interaction data in MITAB 2.5 or 2.7 formats. Details: \code{\link{queryPSICQUIC}}
##' @return data.table for MITAB 2.5: containing the interactor Uniprot IDs, interactor taxonomy IDs, Publication Identifiers, Confidence values and unique pair ID(alphanumerically sorted)
##' @details Output column description:
##' @details pair_id - unique identifier of the undirected interaction: ordered alphabetically and concatenated interacting molecule IDs
##' @details IDs_interactor_A, IDs_interactor_B - interacting molecule ID
##' @details interactor_IDs_databases_A, interactor_IDs_databases_B - database that provides interacting molecule ID such as UniProt, ChEMBL or IntAct (IntAct: when the molecule cannot be mapped to ID in any other resource)
##' @details Taxid_interactor_A, Taxid_interactor_B - taxonomic species that interacting molecule belongs to
##' @details Publication_Identifiers - pubmed ID of the papers that has reported the interaction
##' @details Confidence_values - MIscore. Details: \link{https://psicquic.github.io/MITAB25Format.html}, \link{http://www.ebi.ac.uk/intact/pages/faq/faq.xhtml}, \link{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4316181/}
##' @import data.table
cleanMITAB = function(mitab){
  miformat = NA
  if(ncol(mitab) == 15) miformat = "2.5"
  if(ncol(mitab) == 42) miformat = "2.7"
  mitab = copy(mitab)
  if(is.na(miformat)) stop("the table is not in MITAB 2.5 or 2.7 format, check if any columns were added or deleted from the original query output")

  if(miformat == "2.5"){
    # cleaning Taxid "taxid:9606(human)|taxid:9606(Homo sapiens)" to 9606
    {
      mitab[, Taxid_interactor_A := gsub("taxid:|\\(.*$","",V10)]
      mitab[, Taxid_interactor_B := gsub("taxid:|\\(.*$","",V11)]
      # saving identifier types and cleaning interactor ids
      mitab[, interactor_IDs_databases_A := gsub(":.*$","",V1)]
      mitab[, interactor_IDs_databases_B := gsub(":.*$","",V2)]
      mitab[, IDs_interactor_A := gsub("^.*:","",V1)]
      mitab[, IDs_interactor_B := gsub("^.*:","",V2)]
      # isoform "-1" is a canonical sequence, IntAct uses isoform "-1" when it's clear that the isoform is "-1" and a canonical identifier if it's not clear which isoform was used in the experiment. Removing isoform sign "-1":
      mitab[, IDs_interactor_A := gsub("-1$", "", IDs_interactor_A)]
      mitab[, IDs_interactor_B := gsub("-1$", "", IDs_interactor_B)]
      # cleaning other information
      mitab[, Publication_Identifiers := gsub("^.*pubmed:|\\|.*$","",V9)]
      mitab[, Confidence_values := gsub(".*intact-miscore:","",V15)]
      mitab[, Confidence_values := gsub("-","NA",Confidence_values)]
      mitab[, Confidence_values := as.numeric(Confidence_values)]
      #mitab[, Interaction_identifiers := unlist(gsubfn::strapplyc(Interaction_identifiers,"EBI-[[:digit:]]+",simplify = T)), by =Interaction_identifiers]
      # generating unique identifier for interacting pairs
      mitab[, pair_id := apply(data.table(IDs_interactor_A,IDs_interactor_B,stringsAsFactors = F), 1,
                               function(a) { z = sort(a)
                               paste0(z[1],"|",z[2]) })]
      mitab = mitab[,.(pair_id,
               IDs_interactor_A, IDs_interactor_B,
               interactor_IDs_databases_A, interactor_IDs_databases_B,
               Taxid_interactor_A, Taxid_interactor_B,
               Publication_Identifiers, Confidence_values)]
    }
  }
  mitab = unique(mitab)
  return(mitab)
}
