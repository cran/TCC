# generate negative binomial distributed datasets with different frequencies
NBsample <- function(DEG_foldchange = 4, repA = 3, repB = 3, Ngene = 3000, PDEG = 0.15, PA = 0.2){
    arab <- NULL;rm(arab) # to avoid note by R CMD check
    data(arab) #arab dataset from NBPseq
    data.cl <- c(rep(1, 3), rep(2, 3))
    RPM <- sweep(arab, 2, 1000000/colSums(arab), "*")
    RPM_A <- RPM[,data.cl == 1]
    RPM_B <- RPM[,data.cl == 2]
    RPM_A <- RPM_A[apply(RPM_A, 1, var) > 0,]
    RPM_B <- RPM_B[apply(RPM_B, 1, var) > 0,]
    MEAN <- c(apply(RPM_A, 1, mean), apply(RPM_B, 1, mean))
    VARIANCE <- c(apply(RPM_A, 1, var), apply(RPM_B, 1, var))
    DISPERSION <- (VARIANCE - MEAN)/(MEAN*MEAN)
    mean_disp_tmp <- cbind(MEAN, DISPERSION)
    mean_disp_tmp <- mean_disp_tmp[mean_disp_tmp[,2] > 0,]
    resampling_vector <- sample(1:nrow(mean_disp_tmp), Ngene, replace=TRUE)
    mean_disp <- mean_disp_tmp[resampling_vector,]
    mu <- mean_disp[,1]
    DEG_degree_A <- rep(1, Ngene)
    DEG_degree_A[1:(Ngene*PDEG*PA)] <- DEG_foldchange
    mu_A <- mu*DEG_degree_A
    DEG_degree_B <- rep(1, Ngene)
    DEG_degree_B[(Ngene*PDEG*PA+1):(Ngene*PDEG)] <- DEG_foldchange
    mu_B <- mu*DEG_degree_B
    DEG_posi_org <- (DEG_degree_A*DEG_degree_B) > 1
    nonDEG_posi_org <- (DEG_degree_A*DEG_degree_B) == 1
    outA <- NULL
    colnamev <-NULL
    for(i in 1:repA){
        outA <- cbind(outA, rnbinom(n=length(mu_A), mu=mu_A, size=1/mean_disp[,2]))
        colnamev <-cbind(colnamev, paste("A", as.character(i), sep=""))
    }
    outB <- NULL
    for(i in 1:repB){
        outB <- cbind(outB, rnbinom(n=length(mu_B), mu=mu_B, size=1/mean_disp[,2]))
        colnamev <-cbind(colnamev, paste("B", as.character(i), sep=""))
    }
    out <- cbind(outA, outB)
    colnames(out) <- colnamev
    obj <- rowSums(out) > 0
    RAW <- out[obj,]
    DEG_posi <- DEG_posi_org[obj]
    nonDEG_posi <- nonDEG_posi_org[obj]
    retval <- list(RAW, DEG_posi, nonDEG_posi)
    names(retval) <- c("data", "DEG_posi", "nonDEG_posi")
    return(retval)
}

