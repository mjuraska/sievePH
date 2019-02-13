plot.summary.sievePH <- function(object, mark=NULL, tx=NULL, xlim=NULL, ylim=NULL, xtickAt=NULL, xtickLab=NULL, ytickAt=NULL, ytickLab=NULL, xlab=NULL, ylab=NULL, txLab=c("Placebo", "Treatment"), title=NULL){
  contrast <- names(object)[length(names(object))]

  cexAxis <- 1.3
  cexLab <- 1.4
  cexTitle <- 1.6
  cexText <- 1.2
  cexLegend <- 1.2
  parMar <- c(5, 6, 2, 1)

  par(mar=parMar, oma=rep(0, 4), cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexTitle)

  # a 2-dimensional plot only when the mark is univariate
  if (NCOL(object[[contrast]])==4){
    if (is.null(xlim)){
      xlim <- range(object[[contrast]][, 1])
    }

    if (is.null(ylim)){
      ylim <- range(object[[contrast]][, 2], na.rm=TRUE)
    }

    # need extra room for box plots on the top
    if (!any(c(is.null(mark), is.null(tx)))){
      ylim <- c(ylim[1], ylim[2] + 0.15 * (ylim[2] - ylim[1]))
    }

    if (is.null(xlab)){ xlab <- colnames(object[[contrast]])[1] }
    if (is.null(ylab)){ ylab <- switch(colnames(object[[contrast]])[2], TE="Treatment Efficacy", HR="Hazard Ratio", LogHR="Log Hazard Ratio") }

    plot(object[[contrast]][, 1], object[[contrast]][, 2], xlim=xlim, ylim=ylim, type="n", xlab="", ylab="", xaxt=ifelse(is.null(xtickAt), "s", "n"), yaxt="n", bty="l")

    if (!is.null(xtickAt)){
      if (is.null(xtickLab)){ xtickLab <- xtickAt }
      axis(side=1, at=xtickAt, labels=xtickLab, cex.axis=cexAxis)
    }

    if (!is.null(ytickAt)){
      if (is.null(ytickLab)){ ytickLab <- ytickAt }
      axis(side=2, at=ytickAt, labels=ytickLab, las=1, cex.axis=cexAxis)
    } else {
      # to avoid overlapping tickmarks
      if (!any(c(is.null(mark), is.null(tx)))){
        axis(side=2, at=axTicks(2)[axTicks(2) <= 1.1 * max(object[[contrast]][, 4], na.rm=TRUE)])
      }
    }

    mtext(xlab, side=1, line=3, cex=cexLab)
    mtext(ylab, side=2, line=3.5, las=3, cex=cexLab)

    if (!is.null(title)){ mtext(title, side=3, font=2, line=1, cex=cexTitle, outer=TRUE, at=0, adj=0) }

    abline(h=ifelse(colnames(object[[contrast]])[2]=="HR", 1, 0), col="gray70", lwd=2)

    ### point and interval estimates of VE(v)
    # colCI <- "darkgoldenrod2"
    # colRGB <- c(col2rgb(colCI))
    # colRGB <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha=255*0.55, maxColorValue=255)
    # polygon(c(out$v, rev(out$v)), c(ifelse(out$LBve>-1, out$LBve, -1), rev(ifelse(out$UBve<=1, out$UBve, 1))), col=colRGB, border=NA)

    lines(object[[contrast]][, 1], ifelse(object[[contrast]][, 2] >= ylim[1], object[[contrast]][, 2], NA), lwd=4)
    lines(object[[contrast]][, 1], ifelse(object[[contrast]][, 3] >= ylim[1], object[[contrast]][, 3], NA), lwd=3.5, lty="dashed")
    lines(object[[contrast]][, 1], ifelse(object[[contrast]][, 4] >= ylim[1], object[[contrast]][, 4], NA), lwd=3.5, lty="dashed")

    # text(min(out$v), -0.6, paste0("Marginal Sieve Test P ",ifelse(pMarginalSieve<0.001,"< 0.001",paste0("= ",format(pMarginalSieve,digits=2))),marginalSignifMark), pos=4, cex=cexText)

    #legend("bottomleft", fill=colCI, border=colCI, legend="95% Pointwise CI", cex=cexLegend, bty="n")
    #legend("bottomleft", lwd=c(3.5,2), lty=c("dashed","longdash"), legend=c("95% Pointwise CI","Overall Hazard-Ratio PE"), col=c("black","darkorange"), cex=cexLegend, bty="n")
    legend("bottomleft", lwd=3.5, lty="dashed", legend="95% Pointwise CI", col="black", cex=cexLegend, bty="n")

    # add scatter/box plots of the observed mark values by treatment
    if (!any(c(is.null(mark), is.null(tx)))){
      par(fig=c(0,1,0.85,1), new=TRUE)

      data <- na.omit(cbind(mark, tx))
      plotMarkHoriz(data[, 1], data[, 2], parMar=c(0, parMar[-1]), yLim=xlim, txLab=txLab)

      par(fig=c(0,1,0,1))
    }

  # a 3-dimensional plot (a surface) when the mark is bivariate
  } else if (NCOL(object[[contrast]])==5){
    # STEPHANIE: please fill in using the 'persp' function or something better
  } else {
    stop("Plotting of results is available for univariate and bivariate marks only.")
  }
}

plotMarkHoriz <- function(mark, tx, parMar, yLim, txLab=c("Placebo", "Treatment")){
  cexAxis <- 1.3
  cexLab <- 1.4

  par(mar=parMar, oma=c(0,0,0,0), cex.axis=cexAxis, cex.lab=cexLab)
  boxplot(mark ~ as.factor(tx), at=c(0.5,1.5), xlim=c(0,2), ylim=c(yLim[1], yLim[2]), frame.plot=FALSE, xaxt="n", yaxt="n",
          xlab="", ylab="", boxwex=0.8, outline=FALSE, border="black", lwd=2.5, horizontal=TRUE)
  axis(side=2, at=c(0.5,1.5), labels=txLab, cex.axis=cexAxis, las=1)
  points(mark, jitter(tx + 0.5, factor=0.9), col=ifelse(tx==1, "red3", "blue"), pch=ifelse(tx==1, 24, 21), lwd=2, cex=ifelse(tx==1, 1.2, 1.1))
}
