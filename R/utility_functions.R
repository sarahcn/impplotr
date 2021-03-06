#' Read in imputation metrics files
#' 
#' Reads in metrics files generated by University of Michigan Imputation Server. 
#' Note: optimized for use with reference HRC 1.1 (SNP only)

#' @param imp_dir Imputation directory containing metrics files
#' @param start_chr Start chromosome (1:23)
#' @param end_chr End chromosome (1:23)
#' @param parX logical whether to include PAR (not implemented)

#' @return Dataframe of imputation metrics, from \code{start_chr} to \code{end_chr}

#' @examples
#' \dontrun{
#' mets <- readMetrics("mymetsdir", "proj", start_chr=1, end_chr=23)
#' }

#' @export

readMetrics <- function(imp_dir, start_chr, end_chr, parX=FALSE) {
    
    chr.ranges <- start_chr:end_chr
    chr.ranges.fn <- chr.ranges

    if (parX) {message("\n Note - plotting of PAR region metrics not yet implemented\n")} 
    
    # note that if we have end_chr=23, the server output has 'X'
    # plotting of PAR regions not yet implemented 
    
    if (end_chr %in% 23) {
        # end_chr <- "X"
         message("\nFor chromX, using female metrics\n")
        chr.ranges <- sub(23, "X", chr.ranges)
       #  chrX.no.auto_female.info.gz
        chr.ranges.fn <- sub("X", "X.no.auto_female", chr.ranges)
    }
    
    message("\nReading in imputation metrics on chroms ", start_chr, " through ", 
        end_chr, "...")


    # get metrics file names
    mets.fnames <- paste0(imp_dir, "/chr", chr.ranges.fn, ".info.gz")
    
    # read in the first metrics file to check for expected columns
    # note 'REF(0)' 'ALT(1)' get read in as 'REF.0.'  'ALT.1.'
    minimac.mets.cols <- c("SNP", "REF.0.", "ALT.1.", "ALT_Frq", "MAF", "AvgCall",
                           "Rsq", "Genotyped", "LooRsq", "EmpR", "EmpRsq", "Dose0", "Dose1")
    
    mets.test <- read.table(gzfile(mets.fnames[1]), nrow = 1, header = TRUE, as.is = TRUE)
    
    missing_columns <- setdiff(minimac.mets.cols, names(mets.test))
    
    if (length(missing_columns) > 0) 
        stop("It does not look like you've given me minimac metrics files - they're missing ", 
            paste(missing_columns, collapse = "; "))
    
    ## use 'count.fields' to determine dimension of each per-chr metrics files and
    ## thus dimension of total combined chrom metrics file (rather than rbinding)
    dim.dat <- rep(NA, length(chr.ranges))
    names(dim.dat) <- chr.ranges
    for (i in 1:length(chr.ranges)) {
        mets.fn <- mets.fnames[i]
        Cnt <- system2(command = "zcat", args = paste(mets.fn, "| wc -l"), stdout = TRUE)
        nvars <- as.integer(Cnt) - 1  # subtract header row
        dim.dat[names(dim.dat) == chr.ranges[i]] <- nvars
    }  # close chrom loop to determine variant dimensions
    
    ## loop through selected chroms, combining metrics files add on column for chrom
    mets.all <- data.frame(matrix(nrow = sum(dim.dat),
                                  ncol = length(minimac.mets.cols) +  1))
    
    # start counter for keeping place in combined, annotated metrics files
    cnter <- 1
    for (i in 1:length(chr.ranges)) {
        # report which chrom is being read in, in case of errors reading the combined
        # metrics file
        chr <- chr.ranges[i]
        message("\tReading in metrics for chrom ", chr, "...")
        nvars <- dim.dat[names(dim.dat) == chr]
        mfil <- mets.fnames[i]
        met <- read.table(gzfile(mfil), header = TRUE, as.is = TRUE, comment.char = "", 
            nrow = nvars, na.strings = "-")
        
        # add chrom as final column
        met$chr <- chr
        
        ## add to other chroms' metrics files
        mets.all[cnter:(cnter + nvars - 1), ] <- met
        
        # update counter
        cnter <- cnter + nvars
        
    }  # end loop through chrom columns
    
    names(mets.all) <- c(minimac.mets.cols, "chr")
    
    # additional processing steps we are not doing:
    # flagging SNPs vs indels (HRC r1.1 is only SNPs)
    # make MAF column: minimac already has
    
    # report on number of genotyped vs imputed snps
    nvar.tot <- nrow(mets.all)
    nvar.imp <- sum(mets.all$Genotyped %in% "Imputed")
    message("\nTotal of ", prettyNum(nvar.tot, big.mark = ","),
            " variants; ", prettyNum(nvar.imp, big.mark = ","),
            " of which are imputed\n")
    
    return(mets.all)
    
}


#' Wrapper function to run all plotting and summary functions

#' @param imp_dir Imputation directory
#' @param out_dir Output directory
#' @param start_chr Start chromosome (1:23)
#' @param end_chr End chromosome (1:23)
#' @param sets Sets of imputation summaries
#' @param plots Sets of plots
#' @param fmt Plot output format (only pdf supported)
#' @param summary Logical indicator of whether to run \code{\link{impSummary}}
#' @param rsq_filt Rsq threshold used in \code{\link{maskedPlots}}
#' @param verbose Logical indicator for verbose reporting
#' @param maf_thresh Minor allele frequency threshold for grouping masked variants in \code{\link{maskedSummary}}
#' @param parX logical whether to include PAR (not implemented)
#' @param keep_list File path of an initial variant keep list, with columns variant name and chromosome (integer)

#' @export

makeAllPlots <- function(imp_dir, out_dir, start_chr = 1, end_chr = 23,
                      sets = c("imputed", "masked"),
                      plots = c("boxplot", "info_by_maf", "masked_check_strand"),
                      fmt = "pdf", summary = TRUE, rsq_filt = 0.8,
                      verbose = TRUE, maf_thresh = 0.05,  parX = FALSE,
                         keep_list = NA) {
  
    # I. Read in metrics and annotate with alleles
    mets_all <- readMetrics(imp_dir, start_chr, end_chr, parX)
    
    ## II. Summarize imputed metrics
    if (is.element("imputed", sets) & summary) {
        impSummary(start_chr, end_chr, parX, out_dir,
                   mets_all = mets_all, keep_list = keep_list)
    }
    
    ## III. Graph imputed metrics: info by MAF
    if (is.element("imputed", sets) & is.element("info_by_maf", plots)) {
        impPlotsInfo(out_dir, start_chr, end_chr, parX, fmt,
                     mets_all = mets_all)
    }
    
    ## IV. Graph imputed metrics: by-chrom boxplots
    if (is.element("imputed", sets) & is.element("boxplot", plots)) {
        impPlotsByChr(out_dir, start_chr, end_chr, parX, mets_all = mets_all)
    }
    
    ## V. Graph masked metrics
    if (is.element("masked", sets)) {
        maskedPlots(out_dir, start_chr, end_chr, parX, rsq_filt, 
            mets_all = mets_all)
    }
    
    ## VI. Summarize masked metrics
    if (is.element("masked", sets) & verbose) {
        maskedSummary(maf_thresh, mets_all = mets_all, out_dir, start_chr, end_chr)
    }
    
    ## VII. Masked strand check
    if (is.element("masked", sets) & is.element("masked_check_strand", plots)) {
        maskedStrandCheck(out_dir, start_chr, end_chr, mets_all = mets_all)
    }
}  



