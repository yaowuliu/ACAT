
#'
#' Aggregated Cauchy Assocaition Test
#'
#' A p-value combination method using the Cauchy distribution.
#'
#'
#'
#' @param Weights a numeric vector of non-negative weights for the combined p-values. When it is NULL, the equal weights are used.
#' @param Pvals a numeric vector of p-values to be combined by ACAT.
#' @return p-value of ACAT.
#' @author Yaowu Liu
#' @examples p.values<-c(2e-02,4e-04,0.2,0.1,0.8)
#' @examples ACAT(Pvals=p.values)
#' @export
ACAT<-function(Pvals,Weights=NULL){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
        stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
        stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(sum(Pvals==0)>=1)
    is.one<-(sum(Pvals==1)>=1)
    if (is.zero && is.one){
        stop("Cannot have both 0 and 1 p-values!")
    }
    if (is.zero){
        return(0)
    }
    if (is.one){
        warning("There are p-values that are exactly 1!")
        return(1)
    }

    #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(Weights)){
        Weights<-rep(1/length(Pvals),length(Pvals))
    }else if (length(Weights)!=length(Pvals)){
        stop("The length of weights should be the same as that of the p-values")
    }else if (sum(Weights<0)>0){
        stop("All the weights must be positive!")
    }else{
        Weights<-Weights/sum(Weights)
    }


    #### check if there are very small non-zero p values
    is.small<-(Pvals<1e-16)
    if (sum(is.small)==0){
        cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
    }else{
        cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
        cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
    }
    #### check if the test statistic is very large.
    if (cct.stat>1e+15){
        pval<-(1/cct.stat)/pi
    }else{
        pval<-1-pcauchy(cct.stat)
    }
    return(pval)
}

#'
#' A set-based test that uses ACAT to combine the variant-level p-values.
#'
#'
#' @param G a numeric matrix or dgCMatrix with each row as a different individual and each column as a separate gene/snp. Each genotype should be coded as 0, 1, 2.
#' @param obj an output object of the NULL_Model function.
#' @param weights.beta a numeric vector of parameters for the beta weights for the weighted kernels. If you want to use your own weights, please use the “weights” parameter. It will be ignored if “weights” parameter is not null.
#' @param weights a numeric vector of weights for the SNPs. When it is NULL, the beta weight with the “weights.beta” parameter is used.
#' @param mac.thresh a threshold of the minor allele count (MAC). The Burden test will be used to aggregate the SNPs with MAC less than this thrshold.
#' @return p-value of ACAT-V.
#' @author Yaowu Liu
#' @examples  data(Geno)
#' @examples  G<-Geno[,1:100] # Geno is a dgCMatrix of genotypes
#' @examples  Y<-rnorm(nrow(G)); X<-matrix(rnorm(nrow(G)*4),ncol=4)
#' @examples  obj<-NULL_Model(Y,X)
#' @examples  ACAT_V(G,obj)
#' @export
ACAT_V<-function(G,obj,weights.beta=c(1,25),weights=NULL,mac.thresh=10){
    ### check weights
    if (!is.null(weights) && length(weights)!=ncol(G)){
        stop("The length of weights must equal to the number of variants!")
    }

    mac<-Matrix::colSums(G)
    ### remove SNPs with mac=0
    if (sum(mac==0)>0){
        G<-G[,mac>0,drop=FALSE]
        weights<-weights[mac>0]
        mac<-mac[mac>0]
        if (length(mac)==0){
            stop("The genotype matrix do not have non-zero element!")
        }
    }
    ### p and n
    p<-length(mac)
    n<-nrow(G)
    ###

    if (sum(mac>mac.thresh)==0){  ## only Burden
        pval<-Burden(G,obj, weights.beta = weights.beta, weights = weights)
    }else if (sum(mac<=mac.thresh)==0){ ## only cauchy method
        if (is.null(weights)){
            MAF<-mac/(2*n)
            W<-dbeta(MAF,weights.beta[1],weights.beta[2])/dbeta(MAF,0.5,0.5)
        }else{
            W<-weights
        }

        Mpvals<-Get.marginal.pval(G,obj)
        pval<-ACAT(Mpvals,W)
    }else{  ## Burden + Cauchy method
        is.very.rare<-mac<=mac.thresh
        weights.sparse<-weights[is.very.rare]
        weights.dense<-weights[!is.very.rare]
        pval.dense<-Burden(G[,is.very.rare,drop=FALSE],obj, weights.beta = weights.beta, weights = weights.sparse)

        Mpvals<-Get.marginal.pval(G[,!is.very.rare,drop=FALSE],obj)

        Mpvals<-c(Mpvals,pval.dense)
        if (is.null(weights)){
            MAF<-mac/(2*n)
            mafs<-c(MAF[!is.very.rare],mean(MAF[is.very.rare])) ## maf for p-values
            W<-(dbeta(mafs,weights.beta[1],weights.beta[2])/dbeta(mafs,0.5,0.5))^2
        }else{
            W<-c(weights.dense,mean(weights.sparse))
        }


        is.keep<-rep(T,length(Mpvals))
        is.keep[which(Mpvals==1)]<-F  ## remove p-values of 1.
        pval<-ACAT(Mpvals[is.keep],W[is.keep])
    }
    return(pval)
}

#'
#'
#' Get parameters and residuals from the NULL model
#'
#' Compute model parameters and residuals for ACAT-V
#'
#'
#' @param Y a numeric vector of outcome phenotypes.
#' @param X a numeric matrix of covariates. X must be full-rank.
#' @return This function returns an object that has model parameters and residuals of the NULL model of no association between genetic variables and outcome phenotypes. After obtaining it, please use ACAT-V function to conduct the association test.
#' @author Yaowu Liu
#' @examples  Y<-rnorm(10000)
#' @examples  X<-matrix(rnorm(10000*4),ncol=4)
#' @examples  obj<-NULL_Model(Y,X)
#' @export
NULL_Model<-function(Y,X=NULL){
    n<-length(Y)
    #### check the type of Y
    if ((sum(Y==0)+sum(Y==1))==n){
        out_type<-"D"
    }else{
        out_type<-"C"
    }
    #### Add intercept
    X.tilde<-cbind(rep(1,length(Y)),X)
    #colnames(X.tilde)[1]<-"Intercept"
    if (out_type=="C"){
        #### estimate of sigma square
        X.med<-X.tilde%*%solve(chol(t(X.tilde)%*%X.tilde))   ## X.med%*%t(X.med) is the projection matrix of X.tilde
        Y.res<-as.vector(Y-(Y%*%X.med)%*%t(X.med))
        sigma2<-sum(Y.res^2)/(n-ncol(X.med))
        #### output
        res<-list()
        res[["out_type"]]<-out_type
        res[["X.med"]]<-X.med
        res[["Y.res"]]<-Y.res
        res[["sigma2"]]<-sigma2
    }else if (out_type=="D"){
        #### fit null model
        g<-glm(Y~0+X.tilde,family = "binomial")
        prob.est<-g[["fitted.values"]]
        #### unstandarized residuals
        Y.res<-(Y-prob.est)
        ### Sigma when rho=0
        sigma2.Y<-prob.est*(1-prob.est)  ### variance of each Y_i
        ### output
        res<-list()
        res[["out_type"]]<-out_type
        res[["X.tilde"]]<-X.tilde
        res[["Y.res"]]<-Y.res
        res[["sigma2.Y"]]<-sigma2.Y
    }
    return(res)
}




Get.marginal.pval<-function(G,obj){
    ### check obj
    if (names(obj)[1]!="out_type"){
        stop("obj is not calculated from MOAT_NULL_MODEL!")
    }else{
        out_type<-obj[["out_type"]]
        if (out_type=="C"){
            if (!all.equal(names(obj)[2:length(obj)],c("X.med","Y.res","sigma2"))){
                stop("obj is not calculated from MOAT_NULL_MODEL!")
            }else{
                X.med<-obj[["X.med"]]
                Y.res<-obj[["Y.res"]]
                n<-length(Y.res)
                SST<-obj[["sigma2"]]*(n-ncol(X.med))
            }
        }else if (out_type=="D"){
            if (!all.equal(names(obj)[2:length(obj)],c("X.tilde","Y.res","sigma2.Y"))){
                stop("obj is not calculated from MOAT_NULL_MODEL!")
            }else{
                X.tilde<-obj[["X.tilde"]]
                Y.res<-obj[["Y.res"]]
                sigma2.Y<-obj[["sigma2.Y"]]
                n<-length(Y.res)
            }
        }
    }

    if (class(G)!="matrix" && class(G)!="dgCMatrix"){
        stop("The class of G must be matrix or dgCMatrix!")
    }

    if (out_type=="C"){
        G_tX.med<-as.matrix(Matrix::crossprod(X.med,G))
        ### Sigma^2 of G
        Sigma2.G<-Matrix::colSums(G^2)-Matrix::colSums(G_tX.med^2)
        SSR<-as.vector((Y.res%*%G)^2/Sigma2.G)
        SSR[Sigma2.G<=0]<-0
        df.2<-n-1-ncol(X.med)
        t.stat<-suppressWarnings(sqrt(SSR/((SST-SSR)/df.2)))
        marginal.pvals<-2*pt(t.stat,(n-1-ncol(X.med)),lower.tail = FALSE)
    }else if (out_type=="D"){
        Z.stat0<-as.vector(Y.res%*%G)
        ### Sigma when rho=0
        tG_X.tilde_sigma2<-as.matrix(Matrix::crossprod(G,X.tilde*sigma2.Y))
        Sigma2.G<-Matrix::colSums(G^2*sigma2.Y)-diag(tG_X.tilde_sigma2%*%solve(t(X.tilde)%*%(X.tilde*sigma2.Y))%*%t(tG_X.tilde_sigma2))
        marginal.pvals<-2*pnorm(abs(Z.stat0)/sqrt(Sigma2.G),lower.tail = FALSE)
    }

    return(marginal.pvals)
}


Burden<-function(G,obj,kernel="linear.weighted",weights.beta=c(1,25),weights=NULL){
    ### check obj
    if (names(obj)[1]!="out_type"){
        stop("obj is not calculated from NULL_MODEL!")
    }else{
        out_type<-obj[["out_type"]]
        if (out_type=="C"){
            if (!all.equal(names(obj)[2:length(obj)],c("X.med","Y.res","sigma2"))){
                stop("obj is not calculated from NULL_MODEL!")
            }else{
                X.med<-obj[["X.med"]]
                Y.res<-obj[["Y.res"]]/sqrt(obj[["sigma2"]])  ## rescaled residules
                n<-length(Y.res)
            }
        }else if (out_type=="D"){
            if (!all.equal(names(obj)[2:length(obj)],c("X.tilde","Y.res","sigma2.Y"))){
                stop("obj is not calculated from NULL_MODEL!")
            }else{
                X.tilde<-obj[["X.tilde"]]
                Y.res<-obj[["Y.res"]]
                sigma2.Y<-obj[["sigma2.Y"]]
                n<-length(Y.res)
            }
        }
    }
    ### MAF
    MAF<-Matrix::colSums(G)/(2*dim(G)[1])
    p<-length(MAF)
    #### weights
    if (kernel=="linear.weighted"){
        if (is.null(weights)){
            W<-dbeta(MAF,weights.beta[1],weights.beta[2])
        }else{
            if (length(weights)==p){
                W<-weights
            }else{
                stop("The length of weights must equal to the number of variants!")
            }
        }

    }else if (kernel=="linear"){
        W<-rep(1,p)
    }else{
        stop("The kernel name is not valid!")
    }

    ###### if G is sparse or not
    if (class(G)=="matrix" || class(G)=="dgCMatrix"){
        if (out_type=="C"){
            Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
            Gw<-G%*%W
            sigma.z<-sqrt(sum(Gw^2)-sum((t(X.med)%*%(Gw))^2))
        }else if (out_type=="D"){
            Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
            Gw<-as.vector(G%*%W)
            sigma.z<-sum(Gw^2*sigma2.Y)-((Gw*sigma2.Y)%*%X.tilde)%*%solve(t(X.tilde)%*%(X.tilde*sigma2.Y))%*%t((Gw*sigma2.Y)%*%X.tilde)
            sigma.z<-as.vector(sqrt(sigma.z))
        }
    }else{
        stop("The class of G must be matrix or dgCMatrix!")
    }

    V<-Z.stat.sum/sigma.z
    Q<-V^2   ## Q test statistic
    pval<-1-pchisq(Q,df=1)
    return(pval)
}
