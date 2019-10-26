#' @importFrom grpreg grpreg
estimate.plbpsm <- function (G,criterion,method,ind_c,VS,control,group0,...) {
  ## Do plbpsm estimation and penalty parameter selection...
  if (!method%in%c("SCAD","ALASSO")) stop("unknown variable selection method")
  if (!criterion%in%c("GCV","CV")) stop("unknown smoothing penalty parameter selection criterion")
  
  family <- G$family;
  
  # SMILE #
  if (G$m_b>0){
    Gb=G$basis_info[[length(G$basis_info)]]
  } else {Gb=NULL}
  if (all(group0==0)){MI=FALSE} else {MI=TRUE}
  if (MI){
    XX=G$XX
    for (i in 1: length(G$basis_info_MI)){
      basis_info_MI=G$basis_info_MI[[i]]
      XX <- cbind(XX,basis_info_MI$Bx0)
    }
    group0=c(rep(0,length(G$basis_info_MI)),group0)
    if (!is.null(Gb)){
      BQ2=Gb$B%*%Gb$Q2
      BQ2=as.matrix(BQ2)
      dim_BQ2=dim(BQ2)[2]
      group0=c(rep(0,dim(BQ2)[2]),group0)
    } else {
      BQ2=NULL
      dim_BQ2=0
    }
    p.x=length(G$basis_info_MI)
    colnames(XX)[1:p.x]=colnames(G$Xorig)
    X=G$X
    y=G$y
    fit_smil=LM_ALASSO(y,as.matrix(cbind(BQ2,XX)),group=group0,
                       family=family,...)
    theta.hat=fit_smil$theta.hat
    # change!
    if(length(which(G$group0==0))!=0){
      theta.hat=theta.hat[-which(G$group0==0)]
    }
    ind.l0=which(theta.hat[1+dim_BQ2+1:p.x]!=0,arr.ind=TRUE)
    ind.nl=unique((group0[theta.hat[-1]!=0])[-(1:(dim_BQ2+p.x))])
    G$ind.nl=ind.nl
    # ind.l0=which(theta.hat[tail(1:length(theta.hat),p.x)]!=0,arr.ind=TRUE)
    # ind.nl=unique((group0[theta.hat[-1]!=0])[-c(1:dim(BQ2)[2],tail(1:length(theta.hat[-1]),p.x))])
    ind.l=setdiff(ind.l0,ind.nl)
    G$ind.nl=ind.nl
    G$ind.l=ind.l
    
    cat("Linear selected:",ind.l,"\n")
    cat("Nonlinear selected:",ind.nl,"\n")
    
    ### if there are known linear part, the index should be changed! later needed to be changed.
    linearNULL=as.matrix(G$Xorig[,1:p.x])
    linear_cov=as.matrix(linearNULL[,ind.l])
    ub=c()
    if (length(ind.nl)>0){
      for (j in ind.nl){
        ubb=G$basis_info[[j]]$B
        ub=cbind(ub,ubb)
      }
    } else {ub=NULL}
    G$UB=ub
    if(length(ub)>0){
      ub=ub
    } else {ub=NULL}
    G$X=cbind(linear_cov,ub)
    #G$X=if(G$intercept) cbind(1,G$X)
    G$linear_cov=linear_cov
    linear.names=names(ind.l0)[which(ind.l0 %in% ind.l)]
    nonlinear.names=names(ind.l0)[which(ind.l0 %in% ind.nl)]
    G$linear.names=linear.names
    G$nonlinear.names=nonlinear.names
    if (length(G$nonlinear.names)==0) {G$backfitting<-FALSE}
  } else {
    if (length(G$smooth.univariate)==0){G$backfitting<-FALSE}
    ind.nl=NULL
    ind.l=NULL
  }
  
  object <-
    if (G$am){
    grplsfit(G, criterion=criterion, method=method, family = family, ind_c = ind_c, VS=VS,MI=MI,...)} else {
    ggrplsfit(G, criterion=criterion, method=method, family = family, ind_c = ind_c, VS=VS,control=control,MI=MI)
   }
  object$MI=MI
  object$ind.nl=ind.nl
  object$ind.l=ind.l
  object$scale <- G$scale
  object$criterion <- criterion
  object$method <- method
  object$basis_info<-G$basis_info
  object$basis_info_MI<-G$basis_info_MI
  object$term.names=G$term.names
  object$linear.names=G$linear.names
  object$nonlinear.names=G$nonlinear.names
  
  if (G$intercept &&  G$m_b>0){
    term.names_ni=G$term.names[-1]
  } else {
    term.names_ni=G$term.names
  }
  if(length(term.names_ni)>0){
    if (!MI){
      names(object$coefficients) <- term.names_ni
    } else {
      a=term.names_ni
      if (length(linear.names)>0){
        for (k in 1:length(linear.names)){
          a=a[-grep(linear.names[k],a)]
        }
        #print(a)
        term.names_ni2=if(!is.null(Gb)) {
          c(linear.names,a)
        } else if (G$intercept){
          c(a[1],linear.names,a[2:length(a)])
        } else {c(linear.names,a)}
        names(object$coefficients) <- term.names_ni2
      } else {
        names(object$coefficients) <- term.names_ni
      }
    }
  }
  ## term.names include non as well???
  ## name Ve and se_beta
  if (length(object$Ve)>1){
    object$Ve=as.matrix(object$Ve)
    rownames(object$Ve)[1:length(object$ind_c)] <- colnames(object$Ve)[1:length(object$ind_c)] <- term.names_ni[object$ind_c]
    names(object$se_beta) <- term.names_ni[object$ind_c]
  }
  object
}## end estimate.plbpsm
