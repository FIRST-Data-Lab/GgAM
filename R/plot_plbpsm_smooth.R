#' @importFrom plotly plot_ly add_contour layout %>%
#' @importFrom ggplot2 ggplot geom_raster aes guide_colourbar geom_contour scale_fill_gradientn ggtitle theme
#' geom_polygon element_text element_blank element_line element_rect stat_density2d
#' @importFrom colorRamps matlab.like
#' @importFrom directlabels direct.label
#' @importFrom graphics abline
#' @importFrom utils globalVariables

### The plots information for the bivariate smooth term.
plot.bivariate.smooth=function(x,P=NULL,data=NULL,main,xlim,ylim,xlab,ylab,n1,n2,...){
  ..level..=NULL
  if (is.null(P)) {
    ## later we may adjust this to one dimension and so on.
    if (x$dim!=2) {
      return(NULL)## shouldn't or can't plot
    } else {
      ## get plotting information...
      xterm <- x$term[1]
      if (is.null(xlab)) xlabel <- xterm else xlabel <- xlab
      yterm <- x$term[2]
      if (is.null(ylab)) ylabel <- yterm else ylabel <- ylab
      # raw <- data.frame(x=as.numeric(data[xterm][[1]]),
      #                   y=as.numeric(data[yterm][[1]]))

      ### Set up the fine grid for prediction surface
      n2 <- max(10,n2)
      if (is.null(xlim)) xm <- seq(min(x$V[,1]),max(x$V[,1]),length=n1) else
        xm <- seq(xlim[1],xlim[2],length=n1)
      if (is.null(ylim)) ym <- seq(min(x$V[,2]),max(x$V[,2]),length=n2) else
        ym <- seq(ylim[1],ylim[2],length=n2)
      xx <- rep(xm,n2)
      yy <- rep(ym,each=n1)
      if (is.null(main)) {
        main <- x$label
      }
      if (is.null(ylim)) ylim <- range(ym)
      if (is.null(xlim)) xlim <- range(xm)
      return(list(xx=xx,yy=yy,xm,ym=ym,xlab=xlabel,ylab=ylabel,
             main=main,xlim=xlim,ylim=ylim))
    }
  } else { ## produce plot
    if (x$dim==2) {
      Xpred=cbind(P$xx,P$yy)
      B0=basis(as.matrix(x$V),as.matrix(x$Tr),x$d,x$r,as.matrix(Xpred))
      Bpred=B0$B
      ind=B0$Ind.inside
      u=unique(P$xx)
      v=unique(P$yy)
      n1=length(u)
      n2=length(v)
      mpred=matrix(NA,n1*n2,1)
      mpred[ind,]=as.matrix(Bpred)%*%x$coefficients
      mpred.mtx=matrix(mpred,n1,n2)
      # image(x = u,y = v,
      #                z=(mpred.mtx))
      # contour(x = u,y = v,
      #       z=(mpred.mtx),add = TRUE)
      # ax <- list(
      #   zeroline = FALSE,
      #   showline = FALSE,
      #   showticklabels = FALSE,
      #   showgrid = FALSE,
      #   autotick = FALSE      )
      grid <- data.frame(xx=P$xx,yy=P$yy,mpred=mpred)

      p.beta=ggplot(grid, aes(xx, yy, z = mpred)) +
        geom_raster(aes(fill = mpred),hjust = -0.1, show.legend=TRUE) +# show.lengend=TRUE/lengend=true#geom_line(data = as.data.frame(horseshoe), aes(x=V1, y=V2))  +
        geom_contour(aes(colour = ..level..), color='black', na.rm=T)+ #scale_fill_continuous(type='viridis',guide="colorbar",na.value="white")+
        scale_fill_gradientn(colours = matlab.like(200),na.value = 'white',guide = guide_colourbar(ticks.colour = "black",label.position = "right", label.hjust = 1.5))+
        ggtitle(paste('Predicted Surface for',P$main,sep=' '))+
        theme(plot.title = element_text(hjust = 0.55))+
      theme( plot.background = element_rect(fill = "transparent"),    # Background of the entire plot
               axis.ticks.y = element_blank(),        ## <- this line
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),        ## <- this line
              axis.text.x = element_blank(),
              axis.title.x = element_blank(), ## <- and this line
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background =  element_rect(fill = "white"),
              panel.border = element_blank(),
              legend.title = element_blank(),text=element_text(size=12),legend.key.size = unit(1.5, "cm"),
              legend.key.width = unit(0.6,"cm"),
              legend.box.spacing = unit(0.4, "cm"),
              legend.spacing.x = unit(0.2, 'cm'),
              legend.text=element_text(size=12))#+
        ## change boundaries here horseshoe/geobry/SydneyRealEstateBdry(x=longitude,y=latitude)/bd(x,y)
        # geom_polygon(data=bd,aes(x=x,y=y),inherit.aes=F,
        #              colour='black', fill=NA, lwd=1)
      p1 = direct.label(p.beta, list("bottom.pieces", colour='black'))

      # p.beta=plot_ly(x=u,y=v,z=t(mpred.mtx),
      #                # marker=list(
      #                #   colorscale='Viridis',
      #                #   reversescale =T
      #                # ),
      #                colorscale = "Greys", type = "contour") %>% layout(
      #     xaxis = ax, yaxis = ax,
      #     title = paste('Predicted Surface for',P$main,sep=' ')
        # )
      ### Triangulation plot
      # triplot=create.MESH.2D(nodes=x$V,triangles = x$Tr,...)
      # #plot(triplot)
      # ## Creates the basis
      # FEMbasis = create.FEM.basis(triplot)
    } else {
      warning("no automatic plotting for smooths of variables other than two.")
    }
   ## end of plot production
  return(list(p.beta=p1,triplot=NULL,FEMbasis=NULL,x=P$xm,y=P$ym,p.xlab=P$xlabel,p.ylab=P$ylabel,
              xlab=xlab,ylab=ylab,p.main=P$main,p.ylim=P$ylim,p.xlim=P$xlim,main=main,ylim=ylim,xlim=xlim))
  }
}
