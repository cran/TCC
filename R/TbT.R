#####################################
###  TbT  normalization methods   ###
#####################################
do_TbT <- function(data, data.cl, sample_num = 10000){
    RAW <- data
    ###  Step 1: first normalization  ###
    d <- DGEList(counts=data, group=data.cl)
    d <- calcNormFactors(d)
    norm_f_TMM <- d$samples$norm.factors
    names(norm_f_TMM) <- colnames(data)

    ###  Step 2: DEG identification  ###
    groups <- list(NDE=rep(1, length(data.cl)), DE=data.cl)
    RPM <- sweep(data, 2, 1000000/colSums(data), "*")
    data <- round(RPM)
    once_normalized <- new("countData", data=as.matrix(data), replicates=data.cl, libsizes=colSums(data)*norm_f_TMM, groups=groups)
    once_normalized.NB <- getPriors.NB(once_normalized, samplesize=sample_num, estimation="QL", cl=NULL)
    out <- getLikelihoods.NB(once_normalized.NB, pET="BIC", cl=NULL)
    PDEG <- out@estProps[2] #proportion of differentially expressed genes
    rank_bayseq <- rank(-out@posteriors[,2])
    NDEG <- (nrow(data) * PDEG)# number of differentially expressed genes

    ###  Step 3: second normalization  ###
    obj_DEGy <- (rank_bayseq < NDEG)
    obj_DEGn <- (rank_bayseq >= NDEG)
    data <- RAW[obj_DEGn,]
    d <- DGEList(counts=data, group=data.cl)
    d <- calcNormFactors(d)
    norm_f_TbTorg <- d$samples$norm.factors*colSums(data)/colSums(RAW)
    norm_f_TbT <- norm_f_TbTorg/mean(c(mean(norm_f_TbTorg[data.cl==1]),mean(norm_f_TbTorg[data.cl==2])))

    ###  calculation of PA value (degree of biased expression)  ###
    RPM_TMM <- sweep(RPM, 2, 1/norm_f_TMM, "*")
    data <- RPM_TMM
    logratio <- log2(apply(data[,data.cl==2], 1, mean)) - log2(apply(data[,data.cl==1], 1, mean))
    PA <- sum(logratio[rank_bayseq < NDEG] < 0)/NDEG

    retval <- list(norm_f_TbT, PDEG, PA, obj_DEGn, obj_DEGy, norm_f_TMM, norm_f_TbTorg)
    names(retval) <- c("norm_f_TbT", "PDEG", "PA", "nonDEG_posi", "DEG_posi", "norm_f_TMM", "norm_f_TbTorg")
    return(retval)
}
exactTestafterTbT <- function(names, counts, group, sample_num = 10000, prop.used = 0.5, grid.length=500){
    tbtout <- do_TbT(counts, group, sample_num)
    d <- DGEList(counts=counts, group=group)
    d$samples$norm.factors <- tbtout$norm_f_TbT 
    d <- estimateCommonDisp(d)
    d <- estimateTagwiseDisp(d, prop.used=prop.used, grid.length=grid.length) 
    out <- exactTest(d)
    if(is.vector(out$table$PValue)){#for current edgeR
        FDR <- p.adjust(out$table$PValue, method="BH")
    }else if(is.vector(out$table$p.value)){#for older edgeR
        FDR <- p.adjust(out$table$p.value, method="BH")
    }else{#something strange
        warning("PValue was not available")
    }
    rank_edgeR <- rank(FDR)
    retval <- cbind(names, counts, out$table, FDR, rank_edgeR) 
    return(retval)
}

