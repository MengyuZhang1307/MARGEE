#' Conduct robust single-variant gene-environment interaction tests and joint test for common variants
#'
#' @param null.obj The 'glmmkin' class object. Please use sparse kinship matrix (e.g. dsTMatrix class) when fit null model with 'glmmkin' function in GMMAT package
#' @param interaction a vector contains indices or variable name of the environmental factors
#' @param related.id a numeric or a character vector contains family id
#' @param geno.file the name of a GDS file (including the suffix .gds)
#' @param outfile the directory name of outfile res.txt
#' @param interaction.covariates a vector contains indices or variable name of the interaction covariates
#' @param meta.output boolean value to modiy the output file. If TRUE, the GxE effect estimate and variance and covariance associated with the effect estimate are included in the output file. (default = FALSE)
#' @param center If TRUE, genotypes will be centered before tests. Otherwise, original values will be used in the tests (default = TRUE)
#' @param MAF.range a numeric vector of length 2 defining the minimum and maximum minor allele frequencies of variants that should be included in the analysis (default = c(0.01, 0.5), i.e. common variants).
#' @param miss.cutoff the maximum missing rate allowed for a variant to be included (default = 1, including all variants)
#' @param missing.method method of handling missing genotypes.Either "impute2mean" or "omit" (default = "impute2mean")
#' @param nperbatch an positive integer for how many SNPs should be tested in a batch (default = 100). 
#' @param ncores a positive integer indicating the number of cores used in parallel computing (default = 1)
#' @param is.dosage whether imputed dosage should be used from a GDS infile (default = FALSE).
#' @param verbose  whether failed matrix inversions should be written to outfile.err for debugging (default = FALSE).
#' @export

rom = function(null.obj, outfile, 
    interaction, geno.file, interaction.covariates=NULL, related.id,
    meta.output=F, center=T, MAF.range = c(0.01, 0.5), 
    miss.cutoff = 1, missing.method = "impute2mean", nperbatch = 100, ncores = 1, is.dosage = FALSE, verbose = FALSE){

    if(Sys.info()["sysname"] == "Windows" && ncores > 1) {
        warning("The package doMC is not available on Windows... Switching to single thread...")
        ncores <- 1
    }

    if(!grepl("\\.gds$|\\.bgen$", geno.file[1])) stop("Error: only .gds and .bgen format is supported in geno.file!")
    if(grepl("\\.bgen$", geno.file[1])) stop("Error: currently geno.file must be .gds, .bgen format not yet supported...")
    if(!inherits(null.obj, c("glmmkin", "glmmkin.multi"))) stop("Error: null.obj must be a class glmmkin or glmmkin.multi object!")
    if(inherits(null.obj,"glmmkin.multi")) stop("Error: currently null.obj must be a class glmmkin object, glmmkin.multi not yet supported...")
    #if(meta.output) stop("Error: currently meta output not yet supported...")
    #if(!is.na(null.obj$P)) stop("Error: Please use sparse kinship matrix when fitting the null model")
  
    #if (!dir.exists(outfile)){
    #    dir.create(outfile)
    #} else {
    #    print("Dir already exists!")
    #}
    
    # data manipulation and cleaning
    n.pheno <- null.obj$n.pheno # types of phenotype

    missing.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
    if(inherits(missing.method,"try-error")) stop("Error: \"missing.method\" must be \"impute2mean\" or \"omit\".")
    if(!inherits(interaction, c("integer", "numeric", "character"))) stop("Error: \"interaction\" should be an integer, numeric, or character vector.")
    if(length(related.id) != length(null.obj$id_include)) stop("Error: Cluster ID does not match the length of individual ID")
    

    residuals <- null.obj$scaled.residuals
    ei <- length(interaction)
    qi <- length(interaction.covariates)
    ei1 <- ei+1
    
    if(inherits(interaction,"character")) {
        if(!is.null(interaction.covariates)) {
            if(any(interaction.covariates %in% interaction)) {stop("there are interaction.covariates also specified as interaction.")}
            interaction <- c(interaction, interaction.covariates)
        }
        if (!all(interaction %in% colnames(null.obj$X))) {stop("there are interactions not in column name of covariate matrix.")}
        E <- as.matrix(null.obj$X[,interaction])
    } else { # indices
        if(!is.null(interaction.covariates)) {
          if(any(interaction.covariates %in% interaction)) {stop("there are interaction.covariates also specified as interaction.")}
          interaction <- c(interaction, interaction.covariates)
        }
        E <- as.matrix(null.obj$X[,interaction+1])
    }

    ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
    J <- NULL

    if (!inherits(geno.file, "SeqVarGDSClass")) {
      gds <- SeqArray::seqOpen(geno.file) 
    } else {
      gds <- geno.file
    }
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    variant.idx.all <- SeqArray::seqGetData(gds, "variant.id")
    if(any(is.na(match(null.obj$id_include, sample.id)))) warning("Check your data... Some individuals in null.obj$id_include are missing in sample.id of geno.file!")
    sample.id <- sample.id[sample.id %in% null.obj$id_include]
    if(length(sample.id) == 0) stop("Error: null.obj$id_include does not match sample.id in geno.file!")
    #J <- NULL
    if(any(duplicated(null.obj$id_include))) {
      match.id <- null.obj$id_include %in% sample.id
      null.obj$id_include <- null.obj$id_include[match.id]
      J <- t(sparseMatrix(i=1:length(null.obj$id_include), j=match(null.obj$id_include, unique(null.obj$id_include)[match(sample.id, unique(null.obj$id_include))]), x=1))
    } else match.id <- match(sample.id, null.obj$id_include)

    
    E <- as.matrix(E[match.id, , drop = FALSE])  ###
    
    residuals <- residuals[match.id]
    if(!is.null(null.obj$P)) {
      null.obj$P <- null.obj$P[match.id, match.id]
    } else { # for sparse matrices ???
      null.obj$Sigma_iX <- null.obj$Sigma_iX[match.id, , drop = FALSE]
      null.obj$Sigma_i <- null.obj$Sigma_i[match.id, match.id]
    }
    null.obj$X <- null.obj$X[match.id, , drop = FALSE]
    related.id <- related.id[match.id]
    nfam = length(unique(related.id))
    print(paste0("There are ", nfam, " families"))

    Ebin<-apply(as.matrix(E),2,function(x) length(unique(x))<=20)

    if (any(Ebin)){
      Ecat<-as.matrix(E[,Ebin])
      strata <- apply(Ecat, 1, paste, collapse = "_")
      
      uni.strata<-unique(rev(strata))
      uni.strata<-sort(uni.strata)
      cat_inter<-paste(interaction, collapse = '_')
      tmp<-apply(as.matrix(uni.strata),1,function(x) paste(x, collapse = '_'))
      tmp1<-paste0(cat_inter,"_",tmp)
      tmp2<-c("N","AF")
      bin_header<-c(apply(as.matrix(tmp1),1, function(x) paste0(tmp2,"_",x)))
    }else {
      bin_header=NULL
    }
    strata.cat<-apply(E,2,function(x) length(unique(x))<=20)
    strata <- if (any(strata.cat))  as.numeric(as.factor(apply(as.matrix(E[,strata.cat]), 1, paste, collapse = "_"))) else NULL

   
    E <- scale(E, scale = FALSE)
    E <- cbind(1, E)
    ncolE <- ncol(E)
    if(!is.null(strata)) {
      strata.list <- lapply(sort(unique(strata)), function(x) which(strata==x))
    } else {
      strata.list <- NULL
    }
    if (!inherits(geno.file, "SeqVarGDSClass")) {
      SeqArray::seqClose(gds)
    }
    
    sparsej <- sparseMatrix(i=1:length(related.id), j=match(related.id,unique(related.id)), x=1)

    p.all <- length(variant.idx.all)

    if(ncores>1){
        doMC::registerDoMC(cores = ncores)
        p.percore <- (p.all-1) %/% ncores + 1
        n.p.percore_1 <- p.percore * ncores - p.all
    
        if (meta.output) {
            interaction2 <- c("G", paste0("G-", interaction))
            cov.header = matrix(paste(rep(paste0("Cov_Beta_", interaction2), each = ncolE), interaction2, sep = "_"),ncolE, ncolE)
            meta.header = c(paste0("Beta_", interaction2), paste0("SE_Beta_", interaction2),cov.header[lower.tri(cov.header)],paste0("Robust_SE_BETA_", interaction2[2:(1+ei)]))
        
            if (is.null(strata.list))
            {
                totalCol =  11 + 3*(ei+qi) + ((ei+qi) * ((ei+qi) - 1) / 2 + ei)
            }
            else {
                totalCol =  11 +2*length(unique(strata))+ 3*(ei+qi) + ((ei+qi) * ((ei+qi) - 1) / 2 + ei)
            }
        } else {
            interaction2 <- paste0("G-", interaction[1:ei])
            if (ei != 1) {
                cov.header = matrix(paste(rep(paste0("Cov_Beta_", interaction2), each = ei), interaction, sep = "_G-"), ei, ei)
                meta.header = c(paste0("Beta_", interaction2), paste0("SE_Beta_", interaction2),cov.header[lower.tri(cov.header)], paste0("Robust_SE_BETA_", interaction2))
            } else {
                meta.header = c(paste0("Beta_", interaction2), paste0("SE_Beta_", interaction2),paste0("Robust_SE_BETA_", interaction2))
            }
        
            if (is.null(strata.list)) {totalCol = 9 + ei + ei + ei * (ei - 1) / 2 + ei }
            else {totalCol = 9 +2*length(unique(strata))+ ei + ei + ei * (ei - 1) / 2 + ei }
        }

        if (is.null(strata.list)){
            write.table(t(data.frame(n = c("SNPID","CHR","POS","Non_Effect_Allele","Effect_Allele", "N_Sample", "AF", "Beta_Marginal", "SE_Beta_Marginal", meta.header, "P_Value_Marginal", "P_Value_Interaction", "P_Value_Joint", "Robust_P_Value_Main", "Robust_P_Value_Interaction", "Robust_P_Value_Joint"))), outfile, quote = F, col.names = F, row.names = F, sep="\t")
        } else {
            write.table(t(data.frame(n = c("SNPID","CHR","POS","Non_Effect_Allele","Effect_Allele", "N_Sample", "AF", bin_header, "Beta_Marginal", "SE_Beta_Marginal", meta.header, "P_Value_Marginal",  "P_Value_Interaction", "P_Value_Joint", "Robust_P_Value_Main", "Robust_P_Value_Interaction", "Robust_P_Value_Joint"))), outfile, quote = F, col.names = F, row.names = F, sep="\t")
        }
# b = 1
        foreach(b=1:ncores, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
                debug_file <- paste0(outfile, "_tmp.", b, ".err")
                file.create(debug_file)
            
                variant.idx <- if(b <= n.p.percore_1) variant.idx.all[((b-1)*(p.percore-1)+1):(b*(p.percore-1))] else variant.idx.all[(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1-1)*p.percore+1):(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1)*p.percore)]
                p <- length(variant.idx)
                if (class(geno.file)[1] != "SeqVarGDSClass") {
                gds <- SeqArray::seqOpen(geno.file)
                }
                SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
                #rm(sample.id)
                nbatch.flush <- (p-1) %/% 100000 + 1
                ii <- 0
                for(i in 1:nbatch.flush) {
                    gc()
                    #### each 100000 a batch
                    tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*100000+1):p] else variant.idx[((i-1)*100000+1):(i*100000)]
                    SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
                    MISSRATE <- if(is.dosage) SeqArray::seqApply(gds, "annotation/format/DS", function(xx) mean(is.na(xx)), margin = "by.variant", as.is = "double") else SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
                    AF <- if(is.dosage) SeqArray::seqApply(gds, "annotation/format/DS", mean, margin = "by.variant", as.is = "double", na.rm = TRUE)/2 else 1 - SeqVarTools::alleleFrequency(gds)
                    include <- (MISSRATE <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
                    if(sum(include) == 0) {
                      next
                    }
                    ii <- ii + 1
                    tmp.variant.idx <- tmp.variant.idx[include]
                    tmp.p <- length(tmp.variant.idx)
                    SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
                    SNP <- SeqArray::seqGetData(gds, "annotation/id")
                    SNP[SNP == ""] <- NA
                    out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
                    rm(SNP)
                    alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
                    out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
                    out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
                    out$MISSRATE <- MISSRATE[include]
                    out$AF <- AF[include]
                    include <- include[include]
                    rm(alleles.list)
                    tmp_idx <- 1
                    tmp.out <- lapply(1:((tmp.p-1) %/% nperbatch + 1), function(j){
                        # j = 1
                        if(j == (tmp.p-1) %/% nperbatch + 1) {
                            tmp2.variant.idx = tmp.variant.idx[((j-1)*nperbatch+1):tmp.p]
                        }  else {
                            tmp2.variant.idx = tmp.variant.idx[((j-1)*nperbatch+1):(j*nperbatch)]
                        }
                        SeqArray::seqSetFilter(gds, variant.id = tmp2.variant.idx, verbose = FALSE)
                        geno <- if(is.dosage) SeqVarTools::imputedDosage(gds, use.names = FALSE) else SeqVarTools::altDosage(gds, use.names = FALSE)
                        ng <- ncol(geno)
                        freq <- colMeans(geno, na.rm = TRUE)/2
                        if(any(duplicated(null.obj$id_include))) geno <- crossprod(J, geno)
                        N <- nrow(geno) - colSums(is.na(geno))
                        if(!is.null(strata.list)) { # E is not continuous
                            freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2)
                            if(is.null(ncol(freq.tmp))) freq.tmp <- matrix(freq.tmp, nrow = 1, dimnames = list(NULL, names(freq.tmp)))
                            n.tmp <- sapply(strata.list, function(x) colSums(!is.na(geno[x, , drop = FALSE])))
                            if(is.null(ncol(n.tmp))) n.tmp <- matrix(n.tmp, nrow = 1, dimnames = list(NULL, names(n.tmp)))
                            freq.tmp.rev<-freq.tmp[,order(ncol(freq.tmp):1),drop=FALSE]
                            n.tmp.rev<-n.tmp[,order(ncol(n.tmp):1),drop=FALSE]
                            # combine the freq.tmp and n.tmp by alterating columns
                            rows.freq_N<-nrow(freq.tmp)
                            cols.freq_N<-ncol(freq.tmp)+ncol(n.tmp)
                            freq_N<-matrix(NA,nrow =rows.freq_N,ncol=cols.freq_N)
                            freq_N[,seq(1,cols.freq_N,2)]<-n.tmp
                            freq_N[,seq(2,cols.freq_N,2)]<-freq.tmp
                        }
                    
                        miss.idx <- which(is.na(geno))
                        if(length(miss.idx)>0) {
                            geno[miss.idx] <- if(missing.method == "impute2mean") 2*freq[ceiling(miss.idx/nrow(geno))] else NA
                        }
                        if(center) geno <- scale(geno, scale = FALSE)
                        miss.idx <- which(is.na(geno))
                        if(length(miss.idx)>0) {
                            geno[miss.idx] <- 0
                        }
                    
                        K <- do.call(cbind, sapply(1:ncolE, function(xx) geno*E[,xx], simplify = FALSE), envir = environment()) 
                        U <- as.vector(crossprod(geno, residuals))
                        if(!is.null(null.obj$P)) {
                            PG <- crossprod(null.obj$P, geno)
                        } else {
                            GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
                            PG <- crossprod(null.obj$Sigma_i, geno) - tcrossprod(null.obj$Sigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
                        }

                        GPG <- as.matrix(crossprod(geno, PG)) * (matrix(1, 1, 1) %x% diag(ng)) 
                        GPG_i <- try(solve(GPG), silent = TRUE)
                        if(class(GPG_i)[1] == "try-error") GPG_i <- MASS::ginv(GPG)
                        V_i <- diag(GPG_i)

                        BETA.MAIN <- V_i * U 
                        SE.MAIN   <- sqrt(V_i)
                        STAT.MAIN <- BETA.MAIN * U
                        PVAL.MAIN <- ifelse(V_i>0, pchisq(STAT.MAIN, df=1, lower.tail=FALSE), NA)




                        if(!is.null(null.obj$P)) {
                            KPK <- crossprod(K,crossprod(null.obj$P,K))
                        } else {
                            KSigma_iX <- crossprod(K, null.obj$Sigma_iX)
                            KPK <- crossprod(K, crossprod(null.obj$Sigma_i, K)) - tcrossprod(KSigma_iX, tcrossprod(KSigma_iX, null.obj$cov))
                        }
                        KPK <- as.matrix(KPK) * (matrix(1, ncolE, ncolE) %x% diag(ng))

                        IV.V_i <- try(solve(KPK), silent = TRUE)
                        if(class(IV.V_i)[1] == "try-error") IV.V_i <- try(MASS::ginv(KPK), silent = TRUE)
                        if (class(IV.V_i)[1] == "try-error") {
                            fix_out <- fix.dgesdd(gds, out, debug_file, null.obj, J, residuals, tmp2.variant.idx, meta.output, center, missing.method, strata.list, ncolE, E, ei, meta.header, totalCol, tmp_idx, include)
                            tmp_idx <<- fix_out[[1]]
                            include <<- fix_out[[2]]
                            return(fix_out[[3]])
                        }

                        IV.U <- (rep(1, ncolE) %x% diag(ng)) * as.vector(crossprod(K,residuals))
                        BETA.INT <- crossprod(IV.V_i, IV.U) 

                        ng1   <- ng+1
                        ngei1 <- ng*ei1

                        IV.E_i <- try(solve(IV.V_i[ng1:ngei1, ng1:ngei1]), silent = TRUE)
                        if(class(IV.E_i)[1] == "try-error") IV.E_i <-try(MASS::ginv(IV.V_i[ng1:ngei1, ng1:ngei1]), silent = TRUE)
                        if(class(IV.E_i)[1] == "try-error") {
                            fix_out <- fix.dgesdd(gds, out, debug_file, null.obj, J, residuals,tmp2.variant.idx, meta.output, center, missing.method, strata.list,ncolE, E, ei, meta.header, totalCol, tmp_idx, include)
                            tmp_idx <<- fix_out[[1]]
                            include <<- fix_out[[2]]
                            return(fix_out[[3]])
                        }
                        STAT.INT   <- diag(crossprod(BETA.INT[ng1:ngei1,], crossprod(IV.E_i, BETA.INT[ng1:ngei1,])))

                        IV.GE_i <- try(solve(IV.V_i[1:ngei1, 1:ngei1]), silent = TRUE)
                        if(class(IV.GE_i)[1] == "try-error") IV.GE_i <-try(MASS::ginv(IV.V_i[1:ngei1, 1:ngei1]), silent = TRUE)
                        if(class(IV.GE_i)[1] == "try-error") {
                            fix_out <- fix.dgesdd(gds, out, debug_file, null.obj, J, residuals,tmp2.variant.idx, meta.output, center, missing.method, strata.list,ncolE, E, ei, meta.header, totalCol, tmp_idx, include)
                            tmp_idx <<- fix_out[[1]]
                            include <<- fix_out[[2]]
                            return(fix_out[[3]])
                        }
                        STAT.JOINT <- diag(crossprod(BETA.INT[1:ngei1,], crossprod(IV.GE_i,BETA.INT[1:ngei1,])))
                        PVAL.INT   <- pchisq(STAT.INT, df=ei, lower.tail=FALSE)
                        PVAL.JOINT <- ifelse(is.na(PVAL.MAIN), NA, pchisq(STAT.JOINT, df=1+ei, lower.tail=FALSE))
                    
                    

                    # Adjusted for covariates X
                    
                    
                        GK_X = K - tcrossprod(null.obj$X,tcrossprod(crossprod(K,null.obj$Sigma_iX),null.obj$cov))
                        if(!is.null(null.obj$P)) {
                            PK <- crossprod(null.obj$P,K)
                        } else {
                            PK <- crossprod(null.obj$Sigma_i, K) - tcrossprod(null.obj$Sigma_iX,tcrossprod(KSigma_iX,  null.obj$cov))
                        } 
                    
            
                        true_Residuals = residuals - crossprod(t(PK),(rep(1, ncolE) %x% diag(ng)) * as.vector(tcrossprod(crossprod(residuals,K),IV.V_i))) # residual adjusted for X, G and K. Dim: N x nperbatch (np)
                        true_Score = as.matrix(matrix(GK_X, ncol = (ei1+qi)*ng) * true_Residuals) # N * (pxei+p)
                    
                        M = tcrossprod(crossprod(true_Score, sparsej))
                    
                        batchindex = rep(((0:((ei1+qi)*ng-1)) * (ei1+qi)*ng),each=(ei1+qi))+rep(rep((0:(ng-1)),each=(ei1+qi)),(ei1+qi))+rep(seq(1,(ei1+qi)*ng,by=ng),(ei1+qi)*ng)
                    
                        M[setdiff(1:((ei1+qi)*ng)^2, batchindex)] = 0
                    
                    
                        joint_cov = as.matrix(crossprod(IV.V_i, crossprod(M, IV.V_i)))
                        
                        SE.SW.MAIN = if (ei == 1 & ng == 1) sqrt(joint_cov[1:ng,1:ng]) else sqrt(diag(joint_cov[1:ng,1:ng]))
                    
                        SE.SW.INT = if (ei == 1 & ng == 1) sqrt(joint_cov[ng1:ngei1,ng1:ngei1]) else sqrt(diag(joint_cov[ng1:ngei1,ng1:ngei1]))
                        
                        joint_cov.G_i <- try(solve(joint_cov[1:ng,1:ng]), silent = TRUE) #4x4
                        if(class(joint_cov.G_i)[1] == "try-error") joint_cov.G_i <- MASS::ginv(joint_cov[1:ng,1:ng])
                        
                        SW.STAT.MAIN = diag(crossprod(BETA.INT[1:ng,],crossprod(joint_cov.G_i, BETA.INT[1:ng,])))
                        SW.PVAL.MAIN = pchisq(SW.STAT.MAIN, df=ei, lower.tail=FALSE)
                        SW.PVAL.MAIN = ifelse(SW.PVAL.MAIN == 0, NA, SW.PVAL.MAIN)
                    
                        
                        joint_cov.E_i <- try(solve(joint_cov[ng1:ngei1,ng1:ngei1]), silent = TRUE) #4x4
                        if(class(joint_cov.E_i)[1] == "try-error") joint_cov.E_i <- MASS::ginv(joint_cov[ng1:ngei1,ng1:ngei1])
                        
                        SW.STAT.INT = diag(crossprod(BETA.INT[ng1:ngei1,],crossprod(joint_cov.E_i, BETA.INT[ng1:ngei1,])))
                        SW.PVAL.INT = pchisq(SW.STAT.INT, df=ei, lower.tail=FALSE)
                        SW.PVAL.INT = ifelse(SW.PVAL.INT == 0, NA, SW.PVAL.INT)
                    
                        
                        joint_cov.GE_i <- try(solve(joint_cov[1:ngei1,1:ngei1]), silent = TRUE) #4x4
                        if(class(joint_cov.GE_i)[1] == "try-error") joint_cov.GE_i <- MASS::ginv(joint_cov[1:ngei1,1:ngei1])
                       
                        SW.STAT.JOINT = diag(crossprod(BETA.INT[1:ngei1,], crossprod(joint_cov.GE_i, BETA.INT[1:ngei1,])))
                        SW.PVAL.JOINT = ifelse(is.na(PVAL.MAIN), NA, pchisq(SW.STAT.JOINT, df=1+ei,lower.tail=FALSE))
                        SW.PVAL.JOINT = ifelse(SW.PVAL.JOINT == 0, NA, SW.PVAL.JOINT)
                    
                    
                        split_mat <- matrix(1:(ncolE*ncolE), ncolE, ncolE)
                        if (ng > 1) {
                            IV.V_i <- split(IV.V_i, split_mat %x% diag(ng))[-1]
                            joint_cov <- split(joint_cov, split_mat %x% diag(ng))[-1]
                        }else {
                            IV.V_i <- split(IV.V_i, split_mat %x% diag(ng))
                            joint_cov <- split(joint_cov, split_mat %x% diag(ng))
                        }
                        tmp_idx <<- tmp_idx + ng


                        if (!is.null(strata.list)){
                            if (meta.output) {
                                return(as.matrix(rbind(N, t(freq_N), BETA.MAIN, SE.MAIN, 
                                                diag(as.matrix(BETA.INT[1:ng,])), # Beta G;
                                                t(do.call(cbind, lapply(2:ncolE, function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE and then Beta Covariates
                                                t(sqrt(do.call(cbind, lapply(seq(1,ncolE*ncolE, ncolE+1), function(x) {IV.V_i[[x]]})))),
                                                t(do.call(cbind, lapply(split_mat[lower.tri(split_mat)], function(x) {IV.V_i[[x]]}))),
                                                t(matrix(SE.SW.INT, ncol = ei)), # GxE Robust SE
                                                    PVAL.MAIN, # G Main pval
                                                    PVAL.INT, # GxE pval
                                                    PVAL.JOINT,
                                                    SW.PVAL.MAIN,
                                                    SW.PVAL.INT, # GxE Robust pval
                                                    SW.PVAL.JOINT)))
                            } else {
                                    split_mat <- as.matrix(split_mat[2:(ei+1),2:(ei+1)])
                                    if (length(split_mat) == 1) {
                                        return(as.matrix(rbind(N, t(freq_N), BETA.MAIN, SE.MAIN,
                                                        t(do.call(cbind, lapply(2:(ei+1), function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE only
                                                        t(sqrt(do.call(cbind,lapply(diag(split_mat), function(x) {IV.V_i[[x]]})))),   # SE Beta GxE only
                                                        t(matrix(SE.SW.INT, ncol = ei)),
                                                        PVAL.MAIN,
                                                        PVAL.INT,
                                                        PVAL.JOINT,
                                                        SW.PVAL.MAIN,
                                                        SW.PVAL.INT, # Robust pval
                                                        SW.PVAL.JOINT)))
                                    } else {
                                        return(as.matrix(rbind(N, t(freq_N), BETA.MAIN, SE.MAIN,
                                                        t(do.call(cbind, lapply(2:(ei+1), function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE only
                                                        t(sqrt(do.call(cbind,lapply(diag(split_mat), function(x) {IV.V_i[[x]]})))),   # SE Beta GxE only
                                                        t(do.call(cbind, lapply(split_mat[lower.tri(split_mat)], function(x) {IV.V_i[[x]]}))),
                                                        t(matrix(SE.SW.INT, ncol = ei)),
                                                        PVAL.MAIN,
                                                        PVAL.INT,
                                                        PVAL.JOINT,
                                                        SW.PVAL.MAIN,
                                                        SW.PVAL.INT, # Robust pval
                                                        SW.PVAL.JOINT)))
                                    }
                                }
                        } else {
                            if (meta.output) {
                                    return(as.matrix(rbind(N,  BETA.MAIN, SE.MAIN, 
                                                        diag(as.matrix(BETA.INT[1:ng,])), # Beta G;
                                                        t(do.call(cbind, lapply(2:ncolE, function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE and then Beta Covariates
                                                        t(sqrt(do.call(cbind, lapply(seq(1,ncolE*ncolE, ncolE+1), function(x) {IV.V_i[[x]]})))),
                                                        t(do.call(cbind, lapply(split_mat[lower.tri(split_mat)], function(x) {IV.V_i[[x]]}))),
                                                        t(matrix(SE.SW.INT, ncol = ei)),
                                                        PVAL.MAIN,
                                                        PVAL.INT,
                                                        PVAL.JOINT,
                                                        SW.PVAL.MAIN,
                                                        SW.PVAL.INT, # Robust pval
                                                        SW.PVAL.JOINT)))
                                } else {
                                    split_mat = as.matrix(split_mat[2:(ei+1),2:(ei+1)])
                                    if (length(split_mat) == 1) {
                                        return(as.matrix(rbind(N,BETA.MAIN, SE.MAIN,
                                                    t(do.call(cbind, lapply(2:(ei+1), function(x){diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxEonly 101:200
                                                    t(sqrt(do.call(cbind,lapply(diag(split_mat),function(x) {IV.V_i[[x]]})))),   # SE Beta GxEonly
                                                    t(matrix(SE.SW.INT, ncol = ei)), # GxE Robust SE
                                                    PVAL.MAIN, # G Main pval
                                                    PVAL.INT, # GxE pval
                                                    PVAL.JOINT,
                                                    SW.PVAL.MAIN,
                                                    SW.PVAL.INT, # GxE Robust pval
                                                    SW.PVAL.JOINT)))
                                        } else {
                                            return(as.matrix(rbind(N, BETA.MAIN, SE.MAIN,
                                                        t(do.call(cbind, lapply(2:(ei+1), function(x){diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxEonly
                                                        t(sqrt(do.call(cbind,lapply(diag(split_mat),function(x) {IV.V_i[[x]]})))),   # SE BetaGxE only
                                                        t(do.call(cbind,lapply(split_mat[lower.tri(split_mat)],function(x) {IV.V_i[[x]]}))),
                                                        t(matrix(SE.SW.INT, ncol = ei)),
                                                        PVAL.MAIN,
                                                        PVAL.INT,
                                                        PVAL.JOINT,
                                                        SW.PVAL.MAIN,
                                                        SW.PVAL.INT, # Robust pval
                                                        SW.PVAL.JOINT)))
                                        }
                                }
                        }
                    })
                
                
                    if (!is.null(strata.list)){ 
                        if (any(include)) {
                            out <- out[include,]
                            tmp.out <- matrix(unlist(tmp.out), ncol = totalCol, byrow = TRUE, dimnames = list(NULL,   c("N", bin_header, "BETA.MARGINAL", "SE.MARGINAL", meta.header, "PVAL.MARGINAL", "PVAL.INT",   "PVAL.JOINT", "SW.PVAL.MAIN", "SW.PVAL.INT", "SW.PVAL.JOINT")))
                            out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N", drop = F], out[,"AF",drop=F], tmp.out[,c(bin_header, "BETA.MARGINAL", "SE.MARGINAL", meta.header,"PVAL.MARGINAL", "PVAL.INT", "PVAL.JOINT", "SW.PVAL.MAIN", "SW.PVAL.INT", "SW.PVAL.JOINT"), drop=F])
                            write.table(out, paste0(outfile, "_tmp.", b, ".txt"), quote=FALSE, row.names=FALSE,  col.names=FALSE, sep="\t", append=TRUE, na=".")
                        }
                    } else {
                        if (any(include)) {
                            out <- out[include,]
                            tmp.out <- matrix(unlist(tmp.out), ncol = totalCol, byrow = TRUE, dimnames = list(NULL, c("N",   "BETA.MARGINAL", "SE.MARGINAL", meta.header, "PVAL.MARGINAL", "PVAL.INT",   "PVAL.JOINT", "SW.PVAL.MAIN", "SW.PVAL.INT", "SW.PVAL.JOINT")))
                            out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N", drop = F], out[,"AF",drop=F],   tmp.out[,c( "BETA.MARGINAL", "SE.MARGINAL", meta.header, "PVAL.MARGINAL", "PVAL.INT",   "PVAL.JOINT", "SW.PVAL.MAIN", "SW.PVAL.INT", "SW.PVAL.JOINT"), drop = F])
                            write.table(out, paste0(outfile, "_tmp.", b, ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE, na=".")
                        }
                    }
                    #print(object.size(out))
                    rm(tmp.out)
                    rm(out)
                }
                SeqArray::seqClose(gds)
        }
        if (class(geno.file)[1] == "SeqVarGDSClass") {
                SeqArray::seqClose(geno.file)
        }
        for(b in 1:ncores) {
            system(paste0("cat ", outfile, "_tmp.", b, ".txt", " >> ", outfile))
            system(paste0("cat ", outfile, "_tmp.", b , ".err", " >> ", outfile, ".err"))
            unlink(paste0(outfile, "_tmp.", b, ".txt"))
            unlink(paste0(outfile, "_tmp.", b , ".err"))
        }
    } else { # use single core
            variant.idx <- variant.idx.all
        rm(variant.idx.all)
        p <- length(variant.idx)
        if (!inherits(geno.file, "SeqVarGDSClass")) {
            gds <- SeqArray::seqOpen(geno.file)
        }
        SeqArray::seqSetFilter(gds, sample.id = sample.id, verbose = FALSE)
        rm(sample.id)
        nbatch.flush <- (p-1) %/% 100000 + 1
        ii <- 0
        if (meta.output) {
         interaction2 <- c("G", paste0("G-", interaction))
         cov.header = matrix(paste(rep(paste0("Cov_Beta_", interaction2), each = ncolE), interaction2, sep = "_"),ncolE, ncolE)
         meta.header = c(paste0("Beta_", interaction2), paste0("SE_Beta_", interaction2),cov.header[lower.tri(cov.header)], ,paste0("Robust_SE_BETA_", interaction2[2:(1+ei)]))
         
         if (is.null(strata.list))
         {
             totalCol =  11 + 3*(ei+qi) + ((ei+qi) * ((ei+qi) - 1) / 2 + ei)
         }
         else {
             totalCol =  11 +2*length(unique(strata))+ 3*(ei+qi) + ((ei+qi) * ((ei+qi) - 1) / 2 + ei)
         }
        } else {
         interaction2 <- paste0("G-", interaction[1:ei])
         if (ei != 1) {
             cov.header = matrix(paste(rep(paste0("Cov_Beta_", interaction2), each = ei), interaction, sep = "_G-"), ei, ei)
             meta.header = c(paste0("Beta_", interaction2), paste0("SE_Beta_", interaction2),cov.header[lower.tri(cov.header)], paste0("Robust_SE_BETA_", interaction2))
         } else {
             meta.header = c(paste0("Beta_", interaction2), paste0("SE_Beta_", interaction2),paste0("Robust_SE_BETA_", interaction2))
         }
         
         if (is.null(strata.list)) {totalCol = 9 + ei + ei + ei * (ei - 1) / 2 + ei}
         else {totalCol = 9 +2*length(unique(strata))+ ei + ei + ei * (ei - 1) / 2 + ei}
        }  
 
     if (is.null(strata.list)){
         write.table(t(data.frame(n = c("SNPID","CHR","POS","Non_Effect_Allele","Effect_Allele", "N_Sample", "AF", "Beta_Marginal", "SE_Beta_Marginal", meta.header, "P_Value_Marginal", "P_Value_Interaction", "P_Value_Joint", "Robust_P_Value_Main", "Robust_P_Value_Interaction", "Robust_P_Value_Joint"))), outfile, quote = F, col.names = F, row.names = F, sep="\t")
     } else {
         write.table(t(data.frame(n = c("SNPID","CHR","POS","Non_Effect_Allele","Effect_Allele", "N_Sample", "AF", bin_header, "Beta_Marginal", "SE_Beta_Marginal", meta.header, "P_Value_Marginal",  "P_Value_Interaction", "P_Value_Joint", "Robust_P_Value_Main", "Robust_P_Value_Interaction", "Robust_P_Value_Joint"))), outfile, quote = F, col.names = F, row.names = F, sep="\t")
     }

     debug_file <- paste0(outfile, "_tmp.err")
     file.create(debug_file)
     for(i in 1:nbatch.flush) {
        gc()
        tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*100000+1):p] else variant.idx[((i-1)*100000+1):(i*100000)]
        SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
        MISSRATE <- if(is.dosage) SeqArray::seqApply(gds, "annotation/format/DS", function(xx) mean(is.na(xx)), margin = "by.variant", as.is = "double") else SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
        AF <- if(is.dosage) SeqArray::seqApply(gds, "annotation/format/DS", mean, margin = "by.variant", as.is = "double", na.rm = TRUE)/2 else 1 - SeqVarTools::alleleFrequency(gds)

        include <- (MISSRATE <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))

        if(sum(include) == 0) {
          next
        }
        ii <- ii + 1
        tmp.variant.idx <- tmp.variant.idx[include]
        tmp.p <- length(tmp.variant.idx)
        SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
        SNP <- SeqArray::seqGetData(gds, "annotation/id")
        SNP[SNP == ""] <- NA
        out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
        rm(SNP)
        alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
        out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
        out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
        out$MISSRATE <- MISSRATE[include]
        out$AF <- AF[include]

        include <- include[include]
        rm(alleles.list)
        tmp_idx <- 1
        tmp.out <- lapply(1:((tmp.p-1) %/% nperbatch + 1), function(j) {
          if(j == (tmp.p-1) %/% nperbatch + 1) {
            tmp2.variant.idx = tmp.variant.idx[((j-1)*nperbatch+1):tmp.p]
          }  else {
            tmp2.variant.idx = tmp.variant.idx[((j-1)*nperbatch+1):(j*nperbatch)]
          }   
        SeqArray::seqSetFilter(gds, variant.id = tmp2.variant.idx, verbose = FALSE)
        geno <- if(is.dosage) SeqVarTools::imputedDosage(gds, use.names = FALSE) else SeqVarTools::altDosage(gds, use.names = FALSE)
        ng <- ncol(geno)
        freq <- colMeans(geno, na.rm = TRUE)/2

        if(any(duplicated(null.obj$id_include))) geno <- crossprod(J, geno)

        N <- nrow(geno) - colSums(is.na(geno))
        if(!is.null(strata.list)) { # E is not continuous
            freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2)
            if(is.null(ncol(freq.tmp))) freq.tmp <- matrix(freq.tmp, nrow = 1, dimnames = list(NULL, names(freq.tmp)))
                n.tmp <- sapply(strata.list, function(x) colSums(!is.na(geno[x, , drop = FALSE])))
            if(is.null(ncol(n.tmp))) n.tmp <- matrix(n.tmp, nrow = 1, dimnames = list(NULL, names(n.tmp)))
                freq.tmp.rev<-freq.tmp[,order(ncol(freq.tmp):1),drop=FALSE]
                n.tmp.rev<-n.tmp[,order(ncol(n.tmp):1),drop=FALSE]
            #combine the freq.tmp and n.tmp by alterating columns
                rows.freq_N<-nrow(freq.tmp)
                cols.freq_N<-ncol(freq.tmp)+ncol(n.tmp)
                freq_N<-matrix(NA,nrow =rows.freq_N,ncol=cols.freq_N)
                freq_N[,seq(1,cols.freq_N,2)]<-n.tmp
                freq_N[,seq(2,cols.freq_N,2)]<-freq.tmp
            }else {freq_N<-NA}
          
    
          miss.idx <- which(is.na(geno))
          if(length(miss.idx)>0) {
            geno[miss.idx] <- if(missing.method == "impute2mean") 2*freq[ceiling(miss.idx/nrow(geno))] else NA
          }

          if(center) geno <- scale(geno, scale = FALSE)
          miss.idx <- which(is.na(geno))
          if(length(miss.idx)>0) { # omit
            geno[miss.idx] <- 0
          }

          K <- do.call(cbind, sapply(1:ncolE, function(xx) geno*E[,xx], simplify = FALSE), envir = environment())
          U <- as.vector(crossprod(geno, residuals))

          if(!is.null(null.obj$P)) {
            PG <- crossprod(null.obj$P, geno)
          } else {
            GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
            PG <- crossprod(null.obj$Sigma_i, geno) - tcrossprod(null.obj$Sigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
          }
          
          GPG <- as.matrix(crossprod(geno, PG)) * (matrix(1, 1, 1) %x% diag(ng))
          
          GPG_i <- try(solve(GPG), silent = TRUE)
          if(inherits(GPG_i, "try-error")) GPG_i <- MASS::ginv(GPG)
          V_i <- diag(GPG_i)

          BETA.MAIN <- V_i * U
          
          SE.MAIN   <- sqrt(V_i)
          STAT.MAIN <- BETA.MAIN * U
          PVAL.MAIN <- ifelse(V_i>0, pchisq(STAT.MAIN, df=1, lower.tail=FALSE), NA)

          if(!is.null(null.obj$P)) {
            KPK <- crossprod(K,crossprod(null.obj$P,K))
          } else {
            KSigma_iX <- crossprod(K, null.obj$Sigma_iX)
            KPK <- crossprod(K, crossprod(null.obj$Sigma_i, K)) - tcrossprod(KSigma_iX, tcrossprod(KSigma_iX, null.obj$cov))
          }
   
          KPK <- as.matrix(KPK) * (matrix(1, ncolE, ncolE) %x% diag(ng))

          
          IV.V_i <- try(solve(KPK), silent = TRUE)
     
          if(inherits(IV.V_i, "try-error")) IV.V_i <- try(MASS::ginv(KPK), silent = TRUE)
          if (inherits(IV.V_i, "try-error")) {
            fix_out <- fix.dgesdd(gds, out, debug_file, null.obj, J, residuals, tmp2.variant.idx, meta.output, center, missing.method, strata.list, ncolE, E, ei, meta.header, totalCol, tmp_idx, include)
            tmp_idx <<- fix_out[[1]]
            include <<- fix_out[[2]]
            return(fix_out[[3]])
          }

          IV.U <- (rep(1, ncolE) %x% diag(ng)) * as.vector(crossprod(K,residuals))
          
          BETA.INT <- crossprod(IV.V_i, IV.U)

          
          ng1   <- ng+1
          ngei1 <- ng*ei1
          
          IV.E_i <- try(solve(IV.V_i[ng1:ngei1, ng1:ngei1]), silent = TRUE)
          if(inherits(IV.E_i, "try-error")) IV.E_i <- try(MASS::ginv(IV.V_i[ng1:ngei1, ng1:ngei1]), silent = TRUE)
          if(inherits(IV.E_i, "try-error")) {
            fix_out <- fix.dgesdd(gds, out, debug_file, null.obj, J, residuals, tmp2.variant.idx, meta.output, center, missing.method, strata.list, ncolE, E, ei, meta.header, totalCol, tmp_idx, include)
            tmp_idx <<- fix_out[[1]]
            include <<- fix_out[[2]]
            return(fix_out[[3]])
          }
          STAT.INT   <- diag(crossprod(BETA.INT[ng1:ngei1,], crossprod(IV.E_i, BETA.INT[ng1:ngei1,])))
         
          IV.GE_i <- try(solve(IV.V_i[1:ngei1, 1:ngei1]), silent = TRUE)
          if(inherits(IV.GE_i, "try-error")) IV.GE_i <- try(MASS::ginv(IV.V_i[1:ngei1, 1:ngei1]), silent = TRUE)
          if(inherits(IV.GE_i, "try-error")) {
            fix_out <- fix.dgesdd(gds, out, debug_file, null.obj, J, residuals, tmp2.variant.idx, meta.output, center, missing.method, strata.list, ncolE, E, ei, bin_header, meta.header, totalCol, tmp_idx, include)
            tmp_idx <<- fix_out[[1]]
            include <<- fix_out[[2]]
            return(fix_out[[3]])
          }
          STAT.JOINT <- diag(crossprod(BETA.INT[1:ngei1,], crossprod(IV.GE_i, BETA.INT[1:ngei1,])))
          
          PVAL.INT   <- pchisq(STAT.INT, df=ei, lower.tail=FALSE)
          PVAL.JOINT <- ifelse(is.na(PVAL.MAIN), NA, pchisq(STAT.JOINT, df=1+ei, lower.tail=FALSE))





          GK_X = K - tcrossprod(null.obj$X,tcrossprod(crossprod(K,null.obj$Sigma_iX),null.obj$cov)) # joint Nx(pxei+p)
          if(!is.null(null.obj$P)) {
              PK <- crossprod(null.obj$P,K)
          } else {
              PK <- crossprod(null.obj$Sigma_i, K) - tcrossprod(null.obj$Sigma_iX,tcrossprod(KSigma_iX,  null.obj$cov))
          } # N x (ei+1)*np
          
  
          true_Residuals = residuals - crossprod(t(PK),(rep(1, ncolE) %x% diag(ng)) * as.vector(tcrossprod(crossprod(residuals,K),IV.V_i))) # residual adjusted for X, G and K. Dim: N x nperbatch (np)
          true_Score = as.matrix(matrix(GK_X, ncol = (ei1+qi)*ng) * true_Residuals) # N * (pxei+p)
         
          M = tcrossprod(crossprod(true_Score, sparsej))
          
          batchindex = rep(((0:((ei1+qi)*ng-1)) * (ei1+qi)*ng),each=(ei1+qi))+rep(rep((0:(ng-1)),each=(ei1+qi)),(ei1+qi))+rep(seq(1,(ei1+qi)*ng,by=ng),(ei1+qi)*ng)
          
          M[setdiff(1:((ei1+qi)*ng)^2, batchindex)] = 0
          
          
          joint_cov = as.matrix(crossprod(IV.V_i, crossprod(M, IV.V_i)))
          
          SE.SW.MAIN = if (ei == 1 & ng == 1) sqrt(joint_cov[1:ng,1:ng]) else sqrt(diag(joint_cov[1:ng,1:ng]))             
          SE.SW.INT = if (ei == 1 & ng == 1) sqrt(joint_cov[ng1:ngei1,ng1:ngei1]) else sqrt(diag(joint_cov[ng1:ngei1,ng1:ngei1]))


          #SE.SW.JOINT = sqrt(diag(joint_cov))

          joint_cov.G_i <- try(solve(joint_cov[1:ng,1:ng]), silent = TRUE) #4x4
          if(class(joint_cov.G_i)[1] == "try-error") joint_cov.G_i <- MASS::ginv(joint_cov[1:ng,1:ng])
                        
          SW.STAT.MAIN = diag(crossprod(BETA.INT[1:ng,],crossprod(joint_cov.G_i, BETA.INT[1:ng,])))
          SW.PVAL.MAIN = pchisq(SW.STAT.MAIN, df=ei, lower.tail=FALSE)
          SW.PVAL.MAIN = ifelse(SW.PVAL.MAIN == 0, NA, SW.PVAL.MAIN)
                    

          joint_cov.E_i <- try(solve(joint_cov[ng1:ngei1,ng1:ngei1]), silent = TRUE) #4x4
          if(class(joint_cov.E_i)[1] == "try-error") joint_cov.E_i <- MASS::ginv(joint_cov[ng1:ngei1,ng1:ngei1])
          
          SW.STAT.INT = diag(crossprod(BETA.INT[ng1:ngei1,],crossprod(joint_cov.E_i, BETA.INT[ng1:ngei1,])))
          SW.PVAL.INT = pchisq(SW.STAT.INT, df=ei, lower.tail=FALSE)
          SW.PVAL.INT = ifelse(SW.PVAL.INT == 0, NA, SW.PVAL.INT)
          #interaction covariates 30 min
          
          joint_cov.GE_i <- try(solve(joint_cov[1:ngei1,1:ngei1]), silent = TRUE) #4x4
          if(class(joint_cov.GE_i)[1] == "try-error") joint_cov.GE_i <- MASS::ginv(joint_cov[1:ngei1,1:ngei1])
          
          SW.STAT.JOINT = diag(crossprod(BETA.INT[1:ngei1,], crossprod(joint_cov.GE_i, BETA.INT[1:ngei1,])))                        
          SW.PVAL.JOINT = ifelse(is.na(PVAL.MAIN), NA, pchisq(SW.STAT.JOINT, df=1+ei,lower.tail=FALSE))
          SW.PVAL.JOINT = ifelse(SW.PVAL.JOINT == 0, NA, SW.PVAL.JOINT)
          
          
          split_mat <- matrix(1:(ncolE*ncolE), ncolE, ncolE)
          if (ng > 1) {
          IV.V_i <- split(IV.V_i, split_mat %x% diag(ng))[-1]
          joint_cov <- split(joint_cov, split_mat %x% diag(ng))[-1]
          }else {
          IV.V_i <- split(IV.V_i, split_mat %x% diag(ng))
          joint_cov <- split(joint_cov, split_mat %x% diag(ng))
          }
          tmp_idx <<- tmp_idx + ng
            
          if (!is.null(strata.list)){
                            if (meta.output) {
                                return(as.matrix(rbind(N, t(freq_N), BETA.MAIN, SE.MAIN, 
                                                diag(as.matrix(BETA.INT[1:ng,])), # Beta G;
                                                t(do.call(cbind, lapply(2:ncolE, function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE and then Beta Covariates
                                                t(sqrt(do.call(cbind, lapply(seq(1,ncolE*ncolE, ncolE+1), function(x) {IV.V_i[[x]]})))),
                                                t(do.call(cbind, lapply(split_mat[lower.tri(split_mat)], function(x) {IV.V_i[[x]]}))),
                                                t(matrix(SE.SW.INT, ncol = ei)), # GxE Robust SE
                                                    PVAL.MAIN, # G Main pval
                                                    PVAL.INT, # GxE pval
                                                    PVAL.JOINT,
                                                    SW.PVAL.MAIN,
                                                    SW.PVAL.INT, # GxE Robust pval
                                                    SW.PVAL.JOINT)))
                            } else {
                                    split_mat <- as.matrix(split_mat[2:(ei+1),2:(ei+1)])
                                    if (length(split_mat) == 1) {
                                        return(as.matrix(rbind(N, t(freq_N), BETA.MAIN, SE.MAIN,
                                                        t(do.call(cbind, lapply(2:(ei+1), function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE only
                                                        t(sqrt(do.call(cbind,lapply(diag(split_mat), function(x) {IV.V_i[[x]]})))),   # SE Beta GxE only
                                                        t(matrix(SE.SW.INT, ncol = ei)),
                                                        PVAL.MAIN,
                                                        PVAL.INT,
                                                        PVAL.JOINT,
                                                        SW.PVAL.MAIN,
                                                        SW.PVAL.INT, # Robust pval
                                                        SW.PVAL.JOINT)))
                                    } else {
                                        return(as.matrix(rbind(N, t(freq_N), BETA.MAIN, SE.MAIN,
                                                        t(do.call(cbind, lapply(2:(ei+1), function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE only
                                                        t(sqrt(do.call(cbind,lapply(diag(split_mat), function(x) {IV.V_i[[x]]})))),   # SE Beta GxE only
                                                        t(do.call(cbind, lapply(split_mat[lower.tri(split_mat)], function(x) {IV.V_i[[x]]}))),
                                                        t(matrix(SE.SW.INT, ncol = ei)),
                                                        PVAL.MAIN,
                                                        PVAL.INT,
                                                        PVAL.JOINT,
                                                        SW.PVAL.MAIN,
                                                        SW.PVAL.INT, # Robust pval
                                                        SW.PVAL.JOINT)))
                                    }
                                }
                        } else {
                            if (meta.output) {
                                    return(as.matrix(rbind(N,  BETA.MAIN, SE.MAIN, 
                                                        diag(as.matrix(BETA.INT[1:ng,])), # Beta G;
                                                        t(do.call(cbind, lapply(2:ncolE, function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE and then Beta Covariates
                                                        t(sqrt(do.call(cbind, lapply(seq(1,ncolE*ncolE, ncolE+1), function(x) {IV.V_i[[x]]})))),
                                                        t(do.call(cbind, lapply(split_mat[lower.tri(split_mat)], function(x) {IV.V_i[[x]]}))),
                                                        t(matrix(SE.SW.INT, ncol = ei)),
                                                        PVAL.MAIN,
                                                        PVAL.INT,
                                                        PVAL.JOINT,
                                                        SW.PVAL.MAIN,
                                                        SW.PVAL.INT, # Robust pval
                                                        SW.PVAL.JOINT)))
                                } else {
                                    split_mat = as.matrix(split_mat[2:(ei+1),2:(ei+1)])
                                    if (length(split_mat) == 1) {
                                        return(as.matrix(rbind(N,BETA.MAIN, SE.MAIN,
                                                    t(do.call(cbind, lapply(2:(ei+1), function(x){diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxEonly 101:200
                                                    t(sqrt(do.call(cbind,lapply(diag(split_mat),function(x) {IV.V_i[[x]]})))),   # SE Beta GxEonly
                                                    t(matrix(SE.SW.INT, ncol = ei)), # GxE Robust SE
                                                    PVAL.MAIN, # G Main pval
                                                    PVAL.INT, # GxE pval
                                                    PVAL.JOINT,
                                                    SW.PVAL.MAIN,
                                                    SW.PVAL.INT, # GxE Robust pval
                                                    SW.PVAL.JOINT)))
                                        } else {
                                            return(as.matrix(rbind(N, BETA.MAIN, SE.MAIN,
                                                        t(do.call(cbind, lapply(2:(ei+1), function(x){diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxEonly
                                                        t(sqrt(do.call(cbind,lapply(diag(split_mat),function(x) {IV.V_i[[x]]})))),   # SE BetaGxE only
                                                        t(do.call(cbind,lapply(split_mat[lower.tri(split_mat)],function(x) {IV.V_i[[x]]}))),
                                                        t(matrix(SE.SW.INT, ncol = ei)),
                                                        PVAL.MAIN,
                                                        PVAL.INT,
                                                        PVAL.JOINT,
                                                        SW.PVAL.MAIN,
                                                        SW.PVAL.INT, # Robust pval
                                                        SW.PVAL.JOINT)))
                                        }
                                }
                        }
          
        })

        

        if (!is.null(strata.list)){ 
                    if (any(include)) {
                            out <- out[include,]
                            tmp.out <- matrix(unlist(tmp.out), ncol = totalCol, byrow = TRUE, dimnames = list(NULL,   c("N", bin_header, "BETA.MARGINAL", "SE.MARGINAL", meta.header, "PVAL.MARGINAL", "PVAL.INT",   "PVAL.JOINT", "SW.PVAL.MAIN", "SW.PVAL.INT", "SW.PVAL.JOINT")))
                            out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N", drop = F], out[,"AF",drop=F], tmp.out[,c(bin_header, "BETA.MARGINAL", "SE.MARGINAL", meta.header,"PVAL.MARGINAL", "PVAL.INT", "PVAL.JOINT", "SW.PVAL.MAIN", "SW.PVAL.INT", "SW.PVAL.JOINT"), drop=F])
                            write.table(out, paste0(outfile, "_tmp.", b, ".txt"), quote=FALSE, row.names=FALSE,  col.names=FALSE, sep="\t", append=TRUE, na=".")
                        }
                } else {
                    if (any(include)) {
                        out <- out[include,]
                        tmp.out <- matrix(unlist(tmp.out), ncol = totalCol, byrow = TRUE, dimnames = list(NULL, c("N",   "BETA.MARGINAL", "SE.MARGINAL", meta.header, "PVAL.MARGINAL", "PVAL.INT",   "PVAL.JOINT", "SW.PVAL.MAIN", "SW.PVAL.INT", "SW.PVAL.JOINT")))
                        out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N", drop = F], out[,"AF",drop=F],   tmp.out[,c( "BETA.MARGINAL", "SE.MARGINAL", meta.header, "PVAL.MARGINAL", "PVAL.INT",   "PVAL.JOINT", "SW.PVAL.MAIN", "SW.PVAL.INT", "SW.PVAL.JOINT"), drop = F])
                        write.table(out, outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE, na=".")
                    }
                }
        
        rm(tmp.out)
        rm(out)
      }
      SeqArray::seqClose(gds)

    }
    if(!verbose) unlink(paste0(outfile, ".err"))
    return(invisible(NULL))
    
}

    
    
fix.dgesdd <- function(gds, out, debug_file, null.obj, J, residuals, tmp2.variant.idx, meta.output, center, missing.method, strata.list, ncolE, E, ei, bin_header, meta.header, totalCol, tmp_idx, include) {
  ei1 <- ei+1
  tmp_idx0 <- tmp_idx
  tmp.out <- lapply(tmp2.variant.idx, function(j) {
   
    SeqArray::seqSetFilter(gds, variant.id = j, verbose = FALSE)
    geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
    ng <- ncol(geno)
    freq <- colMeans(geno, na.rm = TRUE)/2
    if(any(duplicated(null.obj$id_include))) geno <- crossprod(J, geno)
    N <- nrow(geno) - colSums(is.na(geno))
    if(!is.null(strata.list)) { # E is not continuous
      freq.tmp <- sapply(strata.list, function(x) colMeans(geno[x, , drop = FALSE], na.rm = TRUE)/2)
      if(is.null(ncol(freq.tmp))) freq.tmp <- matrix(freq.tmp, nrow = 1, dimnames = list(NULL, names(freq.tmp)))
      n.tmp <- sapply(strata.list, function(x) colSums(!is.na(geno[x, , drop = FALSE])))
      if(is.null(ncol(n.tmp))) n.tmp <- matrix(n.tmp, nrow = 1, dimnames = list(NULL, names(n.tmp)))
      freq.tmp.rev<-freq.tmp[,order(ncol(freq.tmp):1),drop=FALSE]
      n.tmp.rev<-n.tmp[,order(ncol(n.tmp):1),drop=FALSE]
      rows.freq_N<-nrow(freq.tmp)
      cols.freq_N<-ncol(freq.tmp)+ncol(n.tmp)
      freq_N<-matrix(NA,nrow =rows.freq_N,ncol=cols.freq_N)
      freq_N[,seq(1,cols.freq_N,2)]<-n.tmp
      freq_N[,seq(2,cols.freq_N,2)]<-freq.tmp
    }
    
    miss.idx <- which(is.na(geno))
    if(length(miss.idx)>0) {
      geno[miss.idx] <- if(missing.method == "impute2mean") 2*freq[ceiling(miss.idx/nrow(geno))] else NA
    }
    if(center) geno <- scale(geno, scale = FALSE)
    miss.idx <- which(is.na(geno))
    if(length(miss.idx)>0) { # omit
      geno[miss.idx] <- 0
    }
            
    K <- do.call(cbind, sapply(1:ncolE, function(xx) geno*E[,xx], simplify = FALSE), envir = environment())

    U <- as.vector(crossprod(geno, residuals))
    if(!is.null(null.obj$P)) {
      PG <- crossprod(null.obj$P, geno)
    } else {
      GSigma_iX <- crossprod(geno, null.obj$Sigma_iX)
      PG <- crossprod(null.obj$Sigma_i, geno) - tcrossprod(null.obj$Sigma_iX, tcrossprod(GSigma_iX, null.obj$cov))
    }

    GPG <- as.matrix(crossprod(geno, PG)) * (matrix(1, 1, 1) %x% diag(ng))
    GPG_i <- try(solve(GPG), silent = TRUE)
    if(inherits(GPG_i, "try-error")) GPG_i <- MASS::ginv(GPG)
    V_i <- diag(GPG_i)
    
    BETA.MAIN <- V_i * U
    SE.MAIN   <- sqrt(V_i)
    STAT.MAIN <- BETA.MAIN * U
    PVAL.MAIN <- ifelse(V_i>0, pchisq(STAT.MAIN, df=1, lower.tail=FALSE), NA)

    if(!is.null(null.obj$P)) {
      KPK <- crossprod(K,crossprod(null.obj$P,K))
    } else {
      KSigma_iX <- crossprod(K, null.obj$Sigma_iX)
      KPK <- crossprod(K, crossprod(null.obj$Sigma_i, K)) - tcrossprod(KSigma_iX, tcrossprod(KSigma_iX, null.obj$cov))
    }
    
    KPK <- as.matrix(KPK) * (matrix(1, ncolE, ncolE) %x% diag(ng))
    IV.V_i <- try(solve(KPK), silent = TRUE)
    if(inherits(IV.V_i, "try-error")) IV.V_i <- try(MASS::ginv(KPK), silent = TRUE)
    if (inherits(IV.V_i, "try-error")) {
      write.table("Variant Info: ", debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table(out[tmp_idx, ], debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table("KPK: ", debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table(KPK, debug_file, row.names = F, col.names = F, quote = F, append = T)
      include[tmp_idx] <<- F
      tmp_idx <<- tmp_idx + 1
      return(NULL)
    }
    
    IV.U <- (rep(1, ncolE) %x% diag(ng)) * as.vector(crossprod(K,residuals))
    BETA.INT <- crossprod(IV.V_i, IV.U)
    
    ng1   <- ng+1
    ngei1 <- ng*ei1
    
    IV.E_i <- try(solve(IV.V_i[ng1:ngei1, ng1:ngei1]), silent = TRUE)
    if(inherits(IV.E_i, "try-error")) IV.E_i <- try(MASS::ginv(IV.V_i[ng1:ngei1, ng1:ngei1]), silent = TRUE)
    if(inherits(IV.E_i, "try-error")) {
      write.table("Variant Info: ", debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table(out[tmp_idx, ], debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table("KPK: ", debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table(KPK, debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table("IV.E_i: ", debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table(IV.V_i[ng1:ngei1, ng1:ngei1], debug_file, row.names = F, col.names = F, quote = F, append = T)
      include[tmp_idx] <<- F
      tmp_idx <<- tmp_idx + 1
      return(NULL)
    }
    STAT.INT   <- diag(crossprod(BETA.INT[ng1:ngei1,], crossprod(IV.E_i, BETA.INT[ng1:ngei1,])))
    
    IV.GE_i <- try(solve(IV.V_i[1:ngei1, 1:ngei1]), silent = TRUE)
    if(inherits(IV.GE_i, "try-error")) IV.GE_i <- try(MASS::ginv(IV.V_i[1:ngei1, 1:ngei1]), silent = TRUE)
    if(inherits(IV.GE_i, "try-error")) {
      write.table("Variant Info: ", debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table(out[tmp_idx, ], debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table("KPK: ", debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table(KPK, debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table("IV.GE_i: ", debug_file, row.names = F, col.names = F, quote = F, append = T)
      write.table(IV.V_i[1:ngei1, 1:ngei1], debug_file, row.names = F, col.names = F, quote = F, append = T)
      include[tmp_idx] <<- F
      tmp_idx <<- tmp_idx + 1
      return(NULL)
    }
    STAT.JOINT <- diag(crossprod(BETA.INT[1:ngei1,], crossprod(IV.GE_i, BETA.INT[1:ngei1,])))
    
    PVAL.INT   <- pchisq(STAT.INT, df=ei, lower.tail=FALSE)
    PVAL.JOINT <- ifelse(is.na(PVAL.MAIN), NA, pchisq(STAT.JOINT, df=1+ei, lower.tail=FALSE))
    
    
    split_mat <- matrix(1:(ncolE*ncolE), ncolE, ncolE)
    if (ng > 1) {
      IV.V_i <- split(IV.V_i, split_mat %x% diag(ng))[-1]
    }else {
      IV.V_i <- split(IV.V_i, split_mat %x% diag(ng))
    }
    tmp_idx <<- tmp_idx + 1
    if (meta.output) {
      return(rbind(N, t(freq_N), BETA.MAIN, SE.MAIN, 
                   diag(as.matrix(BETA.INT[1:ng,])), # Beta G;
                   t(do.call(cbind, lapply(2:ncolE, function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE and then Beta Covariates
                   t(do.call(cbind, lapply(split_mat[lower.tri(split_mat)], function(x) {IV.V_i[[x]]}))),
                   PVAL.MAIN, STAT.INT, PVAL.INT, PVAL.JOINT))
    } else {
      split_mat <- as.matrix(split_mat[2:(ei+1),2:(ei+1)])
      if (length(split_mat) == 1) {
        return(rbind(N, t(freq_N), BETA.MAIN, SE.MAIN,
                     t(do.call(cbind, lapply(2:(ei+1), function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE only
                     t(sqrt(do.call(cbind,lapply(diag(split_mat), function(x) {IV.V_i[[x]]})))),   # SE Beta GxE only
                     PVAL.MAIN, STAT.INT, PVAL.INT, PVAL.JOINT))
      } else {
        return(rbind(N, t(freq_N), BETA.MAIN, SE.MAIN,
                     t(do.call(cbind, lapply(2:(ei+1), function(x) {diag(as.matrix(BETA.INT[(((x-1)*ng)+1):(ng*x),]))}))), # Beta GxE only
                     t(sqrt(do.call(cbind,lapply(diag(split_mat), function(x) {IV.V_i[[x]]})))),   # SE Beta GxE only
                     t(do.call(cbind, lapply(split_mat[lower.tri(split_mat)], function(x) {IV.V_i[[x]]}))),
                     PVAL.MAIN, STAT.INT, PVAL.INT, PVAL.JOINT)) 
      }
    }
    
  })
  if (any(include[tmp_idx0:(tmp_idx-1)])) {
    tmp.out <- matrix(unlist(tmp.out), nrow = totalCol, byrow = F, dimnames = list(c("N", bin_header, "BETA.MARGINAL", "SE.MARGINAL", meta.header, "PVAL.MARGINAL", "STAT.INT", "PVAL.INT", "PVAL.JOINT"), NULL))
    return(list(tmp_idx, include, tmp.out))
  } else {
    return(list(tmp_idx, include, NULL))
  }
  
}











