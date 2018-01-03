## summarize imputed metrics
## optional argument to provide "snp keeplist":
## 2 col df of variant ID and chrom for all variants initially selected into imputation

#' Create counts of imputation basis and target variants

#' @param start_chr Start chromosome (1:23)
#' @param end_chr End chromosome (1:23)
#' @param parX logical whether to include PAR (not implemented)
#' @param out_dir Output directory
#' @param mets_all Dataframe of imputation metrics created by \code{\link{readMetrics}}
#' @param keep_list File path of an initial variant keep list, with columns variant name and chromosome (integer)

#' @return Writes a csv file of variant counts
#' @export

impSummary <- function(start_chr, end_chr, parX=FALSE, out_dir,
                       mets_all, keep_list=NA) {
    message("\nPrinting variant summary by chrom\n")
    
    # get counts for mets_all: imputation basis, study only, imputed
    # then user can combine in anyway desired
    imp.cnts <- unclass(table(mets_all$chr[mets_all$Genotyped %in% "Imputed"]))
    geno.cnts <- unclass(table(mets_all$chr[mets_all$Genotyped %in% "Genotyped"]))
    studyonly.cnts <- unclass(table(mets_all$chr[mets_all$Genotyped %in% "Typed_Only"]))
    
    cnts.dat <- data.frame(chromosome = names(imp.cnts), imp.basis.vars = geno.cnts, 
        imputed.vars = imp.cnts, study.only.vars = studyonly.cnts)
    
    # read in keeplist of snps, optionally
    if(!is.na(keep_list)){
      if(!file.exists(keep_list)) {
        message(keep_list," doesn't exist")
      } else
      message("\nReading in ", keep_list, " to compare varcounts with info file\n")
            kept <- read.table(file = keep_list, comment.char = "",
                               stringsAsFactors = FALSE)
            names(kept)[1:2] <- c("snp", "chr")

      if (sum(is.element(1:22, kept$chr)) %in% 0) {
       message("\tdoesn't look like your keeplist has a 'chr' column - not using")
            } else {

                # for now combine all X (par and non-par) into one
                kept$chr.cnt <- kept$chr
                kept$chr.cnt[kept$chr > 22] <- "X"

                ## count by chrom
                bychr <- as.data.frame(table(kept$chr.cnt))

                ## change chrom from factor to character
                bychr$chromosome <- as.character(bychr$Var1)
                names(bychr)[2] <- "nsnp"

                # add counts from keeplist onto 'cnts.dat'
                cnts.dat$keeplist.vars <- bychr$nsnp[match(cnts.dat$chromosome, bychr$chromosome)]

          } 
    }
    
    # write out csv snp summary
    filo <- paste0(out_dir, "/SNPSummary_chr", start_chr, "-", end_chr, ".csv")
    write.csv(cnts.dat, file = filo, eol = "\n", quote = FALSE, row.names = FALSE)
    
    # report on fraction of imputed variants passing various Rsq thresholds
    
    mets.imp <- mets_all[mets_all$Genotyped %in% "Imputed", ]
    
    message("\nSummaries of imputed variants passing various 'Rsq' score thresholds:")
    
    ## loop through different Rsq thresholds (keeping same as info thresholds)
    for (t in c(0.3, 0.5, 0.8)) {
        message("\nFraction of all imputed variants with Rsq >", t, " is ",
                round(sum(mets.imp$Rsq > t)/nrow(mets.imp), 4))
    }
    
}  # close 'impSummary' function


#' Summarize masked variant metrics 

#' @param maf_thresh Minor allele frequency threshold for grouping masked variants 
#' @param mets_all Dataframe of imputation metrics created by \code{\link{readMetrics}}
#' @param out_dir Output directory
#' @param start_chr Start chromosome (1:23)
#' @param end_chr End chromosome (1:23)

#' @return Writes a csv file of masked metric summaries
#' @export

maskedSummary <- function(maf_thresh=0.05, mets_all, out_dir, start_chr, end_chr) {
    
    # subset metrics to vars in the 'Loo' (leave-one-out) masking expirements
    mets.mask <- mets_all[mets_all$Genotyped %in% "Genotyped", ]
    
    # loop through MAF categories to summarize masked metrics
    type <- "SNP"  # can expand to iterate over type=SNP, type=SV/INDEL
    for (m in maf_thresh) message("\nDichotomizing masked ", type, " metrics by study MAF of ", 
        maf_thresh, "\n")
    
    # save mean and median values as csv
    maf.lab1 <- paste("<", maf_thresh)
    maf.lab2 <- paste(">=", maf_thresh)
    
    df <- data.frame(matrix(ncol = 6, nrow = 4))
    met.nms <- c("EmpRsq", "EmpR", "Dose0")
    names(df) <- c("maf", "nsnp", "metric", met.nms)
    df$maf <- c(rep(maf.lab1, 2), rep(maf.lab2, 2))
    df$metric <- rep(c("mean", "median"), 2)
    
    # make logical flag that's is TRUE for vars below maf thresh
    lm.sel <- mets.mask$MAF < maf_thresh
    
    # fill in avg values for metrics - loop through metrics
    for (met in met.nms) {
        met.vals <- mets.mask[, met]
        df[, met] <- c(round(mean(met.vals[lm.sel]), 3), round(median(met.vals[lm.sel]), 
            3), round(mean(met.vals[!lm.sel]), 3), round(median(met.vals[!lm.sel]), 
            3))
    }
    
    df$nsnp <- c(rep(sum(lm.sel), 2), rep(sum(!lm.sel), 2))
    
    # save as csv
    filo <- paste(out_dir, "/AllMaskedSNPMetrics_chr", start_chr, "-", end_chr, sep = "")
    write.csv(df, file = paste0(filo, "_summary.csv"), quote = FALSE, row.names = FALSE)
    
}  
