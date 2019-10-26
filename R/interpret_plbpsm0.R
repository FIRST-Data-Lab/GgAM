#' @importFrom stats reformulate as.formula terms.formula
interpret.plbpsm0 <- function(gf,textra=NULL,extra.special =NULL)
  # interprets a plbpsm formula of the generic form:
  #   y~x0+x1+x2 + g(w)+ f(u,v) ....
  # and returns:
  # 1. a model formula for the parametric part: pf (and pfok indicating whether it has terms)
  # 2. a list of descriptors for the non-para
  # this is function does the work, and is called by in interpret.plbpsm
{ p.env <- environment(gf) # environment of formula
tf <- terms.formula(gf,specials=c("u","b")) # specials attribute indicates which terms are non-para

terms <- attr(tf,"term.labels") # labels of the model terms
nt <- length(terms) # how many terms?

if (attr(tf,"response") > 0) {  # start the replacement formulae
  response <- as.character(attr(tf,"variables")[2])
} else {
  response <- NULL
}
sp <- c(attr(tf,"specials")$u,attr(tf,"specials")$b)     # array of indices of bivarate smooth terms
#sp2<-attr(tf,"specials")$u # array of univariate smooth terms
off <- attr(tf,"offset") # location of offset in formula

## have to translate sp so that they relate to terms,
## rather than elements of the formula...
vtab <- attr(tf,"factors") # cross tabulation of vars to terms
if (length(sp)>0) for (i in 1:length(sp)) {
  ind <- (1:nt)[as.logical(vtab[sp[i],])]
  sp[i] <- ind # the term that smooth relates to
}
# if (length(sp2)>0) for (i in 1:length(sp2)) {
#   ind2 <- (1:nt)[as.logical(vtab[sp2[i],])]
#   sp2[i] <- ind2 # the term that smooth relates to
# }
## re-referencing is complete

k <-  ks <- kp <- 1 # counters for terms in the 2 formulae
len.sp <- length(sp)
# len.sp2 <- length(sp2)

ns <- len.sp #+ len.sp2 # number of smooths
#ns2 <- len.sp2 # number of smooths

pav <- av <- rep("",0)
smooth.spec <- list()
## Once build our package we may change this
mgcvat <- "package:ggam" %in% search() ## is ggam in search path?
if (nt) for (i in 1:nt) { # work through all terms
  if ( ( k <= ns&&((ks<=len.sp&&sp[ks]==i)) )) { # it's a smooth
    ## have to evaluate in the environment of the formula or you can't find variables
    ## supplied as smooth arguments, e.g. k <- 5;plbpsm(y~s(x,k=k)), fails,
    ## but if you don't specify namespace of mgcv then stuff like
    ## loadNamespace('ggam'); k <- 10; ggam::interpret.plbpsm(y~s(x,k=k)) fails (can't find s)
    ## eval(parse(text=terms[i]),envir=p.env,enclos=loadNamespace('ggam')) fails??
    ## following may supply namespace of ggam explicitly if not on search path...
    if (mgcvat) st <- eval(parse(text=terms[i]),envir=p.env) else {
      st <- try(eval(parse(text=terms[i]),envir=p.env),silent=TRUE)
      if (inherits(st,"try-error")) st <-
          eval(parse(text=terms[i]),enclos=p.env,envir=loadNamespace('GgAM'))
    }
    if (!is.null(textra)) { ## modify the labels on smooths with textra
      pos <- regexpr("(",st$lab,fixed=TRUE)[1]
      st$label <- paste(substr(st$label,start=1,stop=pos-1),textra,
                        substr(st$label,start=pos,stop=nchar(st$label)),sep="")
    }
    smooth.spec[[k]] <- st
    if (ks<=len.sp&&sp[ks]==i) ks <- ks + 1
    k <- k + 1
    # if (ks2<=len.sp2&&sp2[ks2]==i) {
    #     ks2 <- ks2 + 1
    #     k <- k + 1      # counts non-para terms
    #   }
  } else {          # parametric
    av[kp] <- terms[i] ## element kp on rhs of parametric
    kp <- kp+1    # counts parametric terms
  }
}
if (!is.null(off)) { ## deal with offset
  av[kp] <- as.character(attr(tf,"variables")[1+off])
  kp <- kp+1
}

pf <- paste(response,"~",paste(av,collapse=" + "))
if (attr(tf,"intercept")==0) {
  pf <- paste(pf,"-1",sep="")
  if (kp>1) pfok <- 1 else pfok <- 0
} else {
  pfok <- 1;if (kp==1) {
    pf <- paste(pf,"1");
  }
}

fake.formula <- pf
## combine para and non-para
if (length(smooth.spec)>0)
  for (i in 1:length(smooth.spec)) {
    nt <- length(smooth.spec[[i]]$term)
    ff1 <- paste(smooth.spec[[i]]$term[1:nt],collapse="+")
    fake.formula <- paste(fake.formula,"+",ff1)
    # if (smooth.spec[[i]]$by!="NA") {
    #   fake.formula <- paste(fake.formula,"+",smooth.spec[[i]]$by)
    #   av <- c(av,smooth.spec[[i]]$term,smooth.spec[[i]]$by)
    # } else
    av <- c(av,smooth.spec[[i]]$term)
  }
fake.formula <- as.formula(fake.formula,p.env)
if (length(av)) {
  pred.formula <- as.formula(paste("~",paste(av,collapse="+")))
  pav <- all.vars(pred.formula) ## trick to strip out 'offset(x)' etc...
  pred.formula <- reformulate(pav)
} else  pred.formula <- ~1
ret <- list(pf=as.formula(pf,p.env),pfok=pfok,smooth.spec=smooth.spec,
            fake.formula=fake.formula,response=response,fake.names=av,
            pred.names=pav,pred.formula=pred.formula)
class(ret) <- "split.plbpsm.formula"
ret
} ## interpret.bpst0
