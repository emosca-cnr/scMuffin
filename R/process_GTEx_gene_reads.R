#' Process GTEx expression data
#'
#' @param geneReads GTEx gene reads "gct" file, e.g. "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
#' @param GTEx_annot GTEx Sample annotation file, e.g. "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
#' @param tissue tissue name, see param "SMTS" in annotation 
#' @param id_type whether the output should be specified as gene_id or symbol
#' @param nrows limit the read to the top nrows. Useful for testing.
#' @param norm.type Method for normalization. See Seurat::NormalizeData.
#' @param ... further arguments passed to functin Seurat::NormalizeData.
#' @description GTEx gene reads file is processed to obtain the average normalized gene expression in a tissue
#' @return average normalized gene expression values
#' @export
#' @importFrom Seurat NormalizeData
#' @importFrom org.Hs.eg.db org.Hs.egENSEMBL
#' @importFrom utils read.delim read.table


process_GTEx_gene_reads <- function(geneReads=NULL, GTEx_annot=NULL, tissue=NULL, id_type=c("gene_id", "symbol"), nrows=-1, norm.type="LogNormalize", ...) {
    
    id_type <- match.arg(id_type, c("gene_id", "symbol"))
    
    # load data and annotations
    cat("Reading input files: this could take a while...\n")
    GTEx <- read.table(geneReads, skip=2, header = TRUE, sep = "\t", nrows = nrows) 
    annot <- read.delim(GTEx_annot)
 
    cat("Processing data...\n")
    # filter data by tissue
    annot <- annot[ which(annot$SMTS==tissue), ]
    
    if(nrow(annot)<1){
        cat(tissue, " not found. Available tissues:\n")
        print(sort(unique(annot$SMTS)))
    }
    
    # merge GTEx data with annotation via 'SAMPID'
    GTEx <- GTEx[, colnames(GTEx) %in% c("Name", gsub('-', '.', annot$SAMPID))]
    if(nrow(GTEx)<1){
        cat("None of the samples in the annotation file were found in the gene expression dataset.\n")
    }
    
    GTEx$Name <- gsub("\\.[^\\.]+$", "", GTEx$Name) #remove ENSG version to enable the mapping with Entrez
    
    #  calculate mean to retrieve one-column vector + symbol and entrezid
    GTEx_mean <- rowMeans(GTEx[, -1])
    GTEx <- cbind(GTEx, GTEx_mean)
    
    # ANNOTATION
    entrez_ensembl <- as.data.frame(org.Hs.egENSEMBL)
 
    GTEx <- merge(entrez_ensembl, GTEx, by.x="ensembl_id", by.y="Name", sort=FALSE)
    
    #this will keep the highest mean after removing duplicated

    if(id_type=="symbol"){

        entrez_symbol <- as.data.frame(org.Hs.egSYMBOL)
        
        GTEx <- merge(entrez_symbol, GTEx, by="gene_id", sort=FALSE)
        
        GTEx$gene_id <- GTEx$ensembl_id <- NULL
        
        GTEx <- unique(GTEx)
        GTEx <- GTEx[order(-GTEx$GTEx_mean), ]
        GTEx <- GTEx[!duplicated(GTEx$symbol), ]
        
        rownames(GTEx) <- GTEx$symbol
        GTEx$symbol <- GTEx$GTEx_mean <- NULL
        
    }else{
        
        GTEx$ensembl_id <- NULL
        GTEx <- unique(GTEx)
        
        GTEx <- GTEx[order(-GTEx$GTEx_mean), ]
        GTEx <- GTEx[!duplicated(GTEx$gene_id),]
        
        rownames(GTEx) <- GTEx$gene_id
        GTEx$gene_id <- GTEx$GTEx_mean <- NULL
    }
    

    #normalization
    GTEx <- NormalizeData(GTEx, normalization.method = norm.type, scale.factor = 10000, ...)
    
    #calculate the average
    GTEx <- rowMeans(as.data.frame(GTEx))
    GTEx <- as.data.frame(GTEx)
    
    return(GTEx)
    
}
