
#' Main functions
#'
#' The Main function of MDBrainT
#' @param sig_matrix  sig_matrix file path to gene expression from isolated cells, or a matrix of expression profile of cells. Beta values are required.
#'
#' @param mixture_file mixture_file file path to heterogenous mixed expression file, or a matrix of heterogenous mixed expression
#'
#' @import utils
#' @importFrom nnls 
#' @export
#' @examples
#' \dontrun{
#'   ## example 1
#'   sig_matrix <- system.file("extdata", "brain_sig_matrix.txt", package = "MDBrainT")
#'   mixture_file <- system.file("extdata", "example_brain.txt", package = "MDBrainT")
#'   celltype_anno <- system.file("extdata", "celltype_anno.txt", package = "MDBrainT")
#'   results <- MDBrainT(sig_matrix, mixture_file, cell_annotation = celltype_anno)
#'   ## example 2
#'   data(brain_sig_matrix)
#'   data(example_brain)
#'   data(celltype_anno)
#'   results <- MDBrainT(sig_matrix = brain_sig_matrix, mixture_file = example_brain, cell_annotation = celltype_anno)
#' }

require("nnls")


MDBrainT <- function(sig_matrix, mixture_file, cell_annotation){
  
  #read in data
  if (is.character(sig_matrix)) {
    X <- read.delim(sig_matrix, header=T, sep="\t", row.names=1, check.names = F)
    X <- data.matrix(X)
  } else {
    X <- sig_matrix
  }
  
  if (is.character(mixture_file)) {
    Y <- read.delim(mixture_file, header=T, sep="\t", row.names=1, check.names = F)
    Y <- data.matrix(Y)
  } else {
    Y <- mixture_file
  }
  
  if (is.character(mixture_file)) {
    sig.df <- read.delim(cell_annotation, header=T, sep="\t", row.names=1, check.names = F)
    sig.df <- as.data.frame(sig.df)
  } else {
    sig.df <- mixture_file
  }
  
  
  #order
  #X <- X[order(rownames(X)),]
  #Y <- Y[order(rownames(Y)),]

  
  #intersect markers
  keep <- intersect(rownames(X), rownames(Y))
  X <- X[keep, ]
  bVals_decon <- Y[keep, ]
  rownames(sig.df) <- sig.df$id
  keep <- intersect(rownames(X), rownames(bVals_decon))
  sig.df <- sig.df[keep,]

  
  decon.df <- c()
  for (i in 1:ncol(bVals_decon)) {
    #print(i)
    tmp <- bVals_decon[,i]
    #print(colnames(bVals_decon)[i])
    tmp.df <- as.matrix(bVals_decon[,i])
    
    ref.mat <- as.matrix(X[,])
    mod1 <- nnls(ref.mat, tmp.df)
    decon <- mod1$x/sum(mod1$x)
    #colnames(tissue.hyper.mat)
    decon.df <- cbind(decon.df, decon)
    #print("merge Done")
  }
  rownames(decon.df) <- colnames(X[,])
  colnames(decon.df)  <- colnames(bVals_decon)
  decon.df <- round(decon.df, 4)
  decon.df

}



