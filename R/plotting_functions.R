# graph imputed metrics: info by MAF

impPlotsInfo <- function(out_dir, start_chr, end_chr,
                         parX, fmt, maf, mets_all) {
    
    # excerpt imputed vars (all snps)
    mets_all_imp <- mets_all[mets_all$Genotyped %in% "Imputed", ]
    
    nsnps <- nrow(mets_all_imp)
    
    message("\nRsq by MAF plots will represent a total of ",
            prettyNum(nsnps, big.mark = ","), " imputed variants (all SNPs)\n")
    
    # I would never expect a negative Rsq, but build in a safety net
    neg_sel <- with(mets_all_imp, Rsq < 0)
    if (sum(neg_sel) > 0) 
        warning("\nNote there are ", neg_sel, " imputed variants with negative Rsq! These are being included in the graphs (namely because it's unexpected so the code was not built to exclude them)\n")
    
    # PAR not built into pipeline 
    # if (parX)  message("\t including PAR1 and PAR2")
    
    # create output file name
    filo <- paste0(out_dir, "/AllImpMetrics_chr", start_chr, "-", end_chr)
    
    mafs_save <- NULL
    
    mets_curr <- mets_all_imp
    
    set <- "all"
    message("\nplotting info by MAF bin for ", set, " imputed variants\n")
    
    ####### MAF data frame loop
    mafs_list <- seq(0, 0.5, length.out = 51)
    mafs <- as.data.frame(mafs_list[-51])
    names(mafs) <- "maf_min"
    mafs$nsnps <- mafs$mean_rsq <- rep(NA, times = length(mafs$maf_min))
    
    for (i in 1:length(mafs$maf.min)) {
        # bins are defined as [lower, upper)
        maf_sel <- with(mets_curr, MAF > mafs_list[i] & MAF <= mafs_list[i + 1])
        
        # if it's the first bin, grab everything on the left: (0, 0.01)
        if (i %in% 1) {
            maf_sel <- with(mets_curr, MAF >= mafs_list[i] & MAF <= mafs_list[i + 
                1])
        }
        
        ## Count of variants per MAF bin
        mafs$nsnps[i] <- sum(maf_sel)
        ## Mean Rsq
        mafs$mean_rsq[i] <- mean(mets_curr$Rsq[maf_sel])
    }
    ####### close MAF data frame loop
    
    ## plot values at the mid point of maf bin
    mafs$xval <- mafs$maf_min + 0.005
    
    p1 <- ggplot(mafs, aes(x = xval, y = mean_rsq)) + 
                geom_point(shape = 24, fill = "black", size = 3) + 
                xlab("MAF of imputed variants") + 
                ylab("Rsq (mean per MAF bin)") + 
                theme_bw() + ylim(c(0, 1))
    
    p2 <- ggplot(mafs, aes(x = xval, y = nsnps)) + geom_point(shape = 25, fill = "gray47", 
        colour = "gray47", size = 3) + xlab("MAF of imputed variants") + ylab("variant count") + 
        theme_bw()
    
    
    if (fmt == "png") {
        png.fn <- paste(filo, set, "rsq_by_studyMAF.png", sep = ".")
        png(png.fn, width = 900, units = "px")
    }
    if (fmt == "pdf") {
        pdf.fn <- paste(filo, set, "rsq_by_studyMAF.pdf", sep = ".")
        pdf(pdf.fn, width = 12)  # width, height, in inches (default=7)
    }
    
    grid.arrange(p1, p2, ncol = 2)
    dev.off()
    
    # save data frame
    mafs$xval <- NULL
    save(mafs, file = paste0(filo, "rsq_by_MAF.RData"))
    
}

## graph imputed metrics: by-chrom boxplots
impPlotsByChr <- function(out_dir, start_chr, end_chr, parX, mets_all) {
    
    # excerpt imputed vars (all snps)
    mets_all_imp <- mets_all[mets_all$Genotyped %in% "Imputed", ]
    
    nsnps <- nrow(mets_all_imp)
    
    message("\nRsq by-chrom box plots will represent a total of ", prettyNum(nsnps, 
        big.mark = ","), " imputed variants (all SNPs)\n")
    
    # I would never expect a negative Rsq, but build in a safety net
    neg.sel <- with(mets_all_imp, Rsq < 0)
    if (sum(neg.sel) > 0) 
        warning("\nNote there are ", neg.sel, " imputed variants with negative Rsq! These are being included in the graphs (namely because it's unexpected so the code was not built to exclude them)\n")
    
    mets_curr <- mets_all_imp
    
    # create a chrom col that will order 1-X along X axis will likely want to adjust
    # for PARs vs nonPARs
    mets_curr$chr.int.use <- as.integer(mets_curr$chr)
    mets_curr$chr.int.use[mets_curr$chr %in% "X"] <- 23
    
    mets_curr$chr.fac.use <- factor(mets_curr$chr.int.use, ordered = TRUE, levels = 1:23, 
        labels = c(1:22, "X"))
    
    # using ggplot excluding outliers from plot
    
    set <- "all"
    message("making boxplots for ", set, " variants")
    
    # create output file name
    filo <- paste0(out_dir, "/AllImpMetrics_chr", start_chr, "-", end_chr)
    set <- "all"
    png.fn <- paste(filo, set, "boxplot_%d.png", sep = ".")
    png(png.fn)
    
    # Rsq
    print(ggplot(mets_curr, aes(x = chr.fac.use, y = Rsq)) + geom_boxplot(fill = "seashell2", 
        outlier.shape = NA) + xlab("chromosome") + ylab("Rsq") + theme_bw())
    
    # MAF. log10 scale the y axis - will exclude monomorphs
    nmono <- sum(mets_curr$MAF %in% 0)
    subt <- paste("*excluding", prettyNum(nmono, big.mark = ","), "monomorphic variants")
    
    print(ggplot(mets_curr, aes(x = chr.fac.use, y = MAF)) + geom_boxplot(fill = "seashell2", 
        outlier.shape = NA) + scale_y_log10() + xlab("chromosome") + ylab("imputed MAF (log10 scaling)*") + 
        theme_bw() + annotate("text", label = subt, x = unique(mets_curr$chr.fac.use)[1], 
        y = 0, colour = "red", fontface = "italic", vjust = 0, hjust = 0))
    
    dev.off()
    
}  

## graph masked metrics
maskedPlots <- function(out_dir, start_chr, end_chr, parX, maf, rsq_filt, mets_all) {
    
    # subset metrics to vars in the 'Loo' (leave-one-out) masking expirements
    mets.mask <- mets_all[mets_all$Genotyped %in% "Genotyped", ]
    
    nsnps <- nrow(mets.mask)
    message("\n\nTotal number of masked SNPs will be ", prettyNum(nsnps, big.mark = ","))
    
    # may need to edit this stmt once we understand chrX par and non-par output
    # if (parX & 23 %in% start_chr:end_chr) 
    #    message("\t including PAR1 and PAR2")
    
    # create MAF bins and define mean masked SNP metrics per bin
    mafs_list <- seq(0, 0.5, length.out = 51)
    mafs <- as.data.frame(mafs_list[-51])
    names(mafs) <- "maf.min"
    nas <- rep(NA, times = length(mafs$maf.min))
    mafs$nsnps <- mafs$meanEmpRsq <- mafs$meanEmpR <- mafs$meanDose0 <- nas
    mafs.blank <- mafs  # perserve empty copy
    
    for (i in 1:length(mafs$maf.min)) {
        
        # bins are defined as [lower, upper) (I think that's the correct notation!)
        maf.sel <- with(mets.mask, MAF > mafs_list[i] & MAF <= mafs_list[i + 1])
        
        # if it's the first bin, grab everything on the left: (0, 0.01)
        if (i %in% 1) {
            maf.sel <- with(mets.mask, MAF >= mafs_list[i] & MAF <= mafs_list[i + 
                1])
        }
        
        ## Count of variants per MAF bin
        mafs$nsnps[i] <- sum(maf.sel)
        
        ## Mean of masked snp metrics per bin
        mafs$meanEmpRsq[i] <- mean(mets.mask$EmpRsq[maf.sel])
        mafs$meanEmpR[i] <- mean(mets.mask$EmpR[maf.sel])
        mafs$meanDose0[i] <- mean(mets.mask$Dose0[maf.sel])
        
    }  # close MAF data frame loop
    
    ## remove bins with no snps
    mafs <- mafs[mafs$nsnps > 0, ]
    mafs.all <- mafs
    
    ## repeat thresholding on LooRsq (like Rsq, but from internal masking)
    mets.filt <- mets.mask[mets.mask$LooRsq >= rsq_filt, ]
    psnps <- round((nrow(mets.filt)/nsnps), 5) * 100
    
    message("\n\nPercent of masked SNPs passing LooRsq threshold: ", prettyNum(psnps, 
        big.mark = ","), "%")
    
    mafs <- mafs.blank
    for (i in 1:length(mafs$maf.min)) {
        
        ## bins are defined as [lower, upper) (I think that's the correct notation!)
        maf.sel <- with(mets.filt, MAF > mafs_list[i] & MAF <= mafs_list[i + 1])
        
        ## if it's the first bin, grab everything on the left: (0, 0.01)
        if (i %in% 1) {
            maf.sel <- with(mets.filt, MAF >= mafs_list[i] & MAF <= mafs_list[i + 
                1])
        }
        
        ## Count of variants per MAF bin
        mafs$nsnps[i] <- sum(maf.sel)
        
        ## Mean of masked snp metrics per bin
        mafs$meanEmpRsq[i] <- mean(mets.filt$EmpRsq[maf.sel])
        mafs$meanEmpR[i] <- mean(mets.filt$EmpR[maf.sel])
        mafs$meanDose0[i] <- mean(mets.filt$Dose0[maf.sel])
        
    }  # close MAF data frame loop
    
    ## keep filtered maf bins in overall maf bins
    mafs.filt <- mafs
    mafs.filt <- mafs.filt[mafs.filt$maf.min %in% mafs.all$maf.min, ]
    mafs.filt$frac.snps <- round(mafs.filt$nsnps/mafs.all$nsnps, 3)
    
    # create base name for all masked metrics plots
    filo <- paste(out_dir, "/AllMaskedSNPMetrics_chr", start_chr, "-", end_chr, sep = "")
    
    ## combine two data frames of average metrics by MAF bin
    mafs.all$frac.snps <- 1
    mafs.comb <- rbind(mafs.all, mafs.filt)
    save(mafs.comb, file = paste(filo, "mets_by_MAF.RData", sep = ""))
    
    ## plot values at the mid point of maf bin
    mafs.filt$xval <- mafs.filt$maf.min + 0.005
    mafs.all$xval <- mafs.all$maf.min + 0.005
    
    ## use grid arrange for 4x4 plot: (1) SNP count; (2) %pct passing Rsq filt; (3)
    ## EmpR2; (4) Dose1
    
    p1 <- ggplot(mafs.all, aes(x = xval, y = nsnps)) + geom_point(shape = 16, colour = "black", 
        size = 2) + theme_bw() + xlab("study MAF") + ylab("SNP count") + ggtitle("(A)")
    
    p2 <- ggplot(mafs.filt, aes(x = xval, y = frac.snps)) + geom_point(shape = 17, 
        colour = "gray47", size = 2) + theme_bw() + xlab("study MAF") + ylab(paste("fraction of SNPs with LooRsq >=", 
        rsq_filt)) + ggtitle("(B)")
    
    # EmpRsq - LooRsq filtered as second series of points
    p3 <- ggplot(mafs.all, aes(x = xval, y = meanEmpRsq)) + geom_point(shape = 16, 
        colour = "black", size = 3) + theme_bw() + xlab("study MAF") + ylab("empirical r2 (EmpRsq)") + 
        geom_point(data = mafs.filt, colour = "gray47", shape = 17, size = 2, aes(x = xval, 
            y = meanEmpRsq)) + ggtitle("(C)")
    
    # 'Dose1' metric - LooRsq filtered as second series of points
    p4 <- ggplot(mafs.all, aes(x = xval, y = meanDose0)) + geom_point(shape = 16, 
        colour = "black", size = 3) + theme_bw() + xlab("study MAF") + ylab("Dose0") + 
        geom_point(data = mafs.filt, colour = "gray47", shape = 17, size = 2, aes(x = xval, 
            y = meanDose0)) + ggtitle("(D)")
    
    
    pdf.fn <- paste(filo, "comboplot.pdf", sep = ".")
    pdf(pdf.fn, width = 10, height = 10)  # need to specify width and height?
    
    grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    
    dev.off()
    
} 

## plot masked SNP strand check
maskedStrandCheck <- function(out_dir, start_chr, end_chr, mets_all) {

    message("\nMaking masked SNP strand check plot\n")
    
    # limit to masked snps
    mets.mask <- mets_all[mets_all$Genotyped %in% "Genotyped", ]
    
    # flag strand ambiguous SNPs
    mets.mask$alleles <- with(mets.mask, pasteSorted(REF.0., ALT.1.))
    mets.mask$strand.amb <- is.element(mets.mask$alleles, c("A/T", "C/G"))
    
    ## ggplot: High Density Scatterplot with Color Transparency
    cols.transp <- rgb(0, 100, 0, 50, maxColorValue = 255)
    leg.name <- "Strand Ambiguous (A/T or C/G) SNPs"
    
    # make scatterplot unlike IMPUTE2, we don't have masked concordance so using
    # EmpRsq instead you get line of EmpRsq near 0 mostly b/c of rare variants could
    # try reducing to MAF > .1% - or just explain in narrative [mets.mask$MAF >
    # 0.001,],
    
    plot_scatter <- ggplot(mets.mask, aes(x = LooRsq, y = EmpRsq, color = strand.amb, 
        shape = strand.amb)) + geom_point() + scale_shape_manual(values = c(16, 8), 
        name = leg.name) + scale_color_manual(values = c(cols.transp, "black"), name = leg.name) + 
        # xlim(0,1) + ylim(0,1) +
    theme(legend.position = "bottom") + geom_abline(intercept = 0, slope = 1, linetype = "dashed")
    
    # see http://www.r-bloggers.com/ggplot2-cheatsheet-for-visualizing-distributions/
    # placeholder plot - prints nothing at all
    empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        panel.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
    
    # marginal density of x - plot on top
    plot_top <- ggplot(mets.mask, aes(LooRsq)) + geom_density() + theme(legend.position = "none") + 
        xlab("") + ylab("")
    
    # marginal density of y - plot on the right
    plot_right <- ggplot(mets.mask, aes(EmpRsq)) + geom_density() + coord_flip() + 
        theme(legend.position = "none") + xlab("") + ylab("")
    
    # arrange the plots together, with appropriate height and width for each row and
    # column note ggsave doesn't work with grid.arrange
    filo <- paste(out_dir, "/AllMaskedSNPMetrics_chr", start_chr, "-", end_chr, sep = "")
    png.fn <- paste(filo, "check_strand.png", sep = ".")
    png(png.fn)
    grid.arrange(plot_top, empty, plot_scatter, plot_right, ncol = 2, nrow = 2, widths = c(4, 
        1), heights = c(1, 4))
    dev.off()
    
    # report number of snps in 'problem quadrant'
    probsel <- with(mets.mask, EmpRsq < 0.5 & LooRsq > 0.5)
    nprob <- sum(probsel)
    nprob.amb <- sum(probsel & with(mets.mask, strand.amb & MAF > 0.001))
    message("Of ", prettyNum(nrow(mets.mask), big.mark = ","), " total masked SNPs, ", 
        prettyNum(nprob, big.mark = ","), " are in problem quadrant of strand check (LooRsq > 0.5 and EmpRsq < 0.5) - ", 
        prettyNum(nprob.amb, big.mark = ","), " of which are also strand ambiguous and MAF > 0.1%\n\n")
    
    # subset to strand ambiguous, excluding very rare, which seem to have low EmpRsq
    # by nature of being rare....?
    
}  # end masked strand check function
