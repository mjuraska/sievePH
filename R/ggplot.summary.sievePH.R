#' Plotting Mark-Specific Proportional Hazards Model Fits Using ggplot
#'
#' \code{plot} method for class \code{summary.sievePH}. For univariate marks, it plots point and interval estimates of the mark-specific treatment effect parameter specified by \code{contrast} in \code{\link{summary.sievePH}}, and,
#' optionally, scatter/box plots of the observed mark values by treatment. 
#' @param x an object returned by \code{\link{summary.sievePH}}
#' @param mark either a numeric vector specifying a univariate continuous mark or a data frame specifying a multivariate continuous mark.
#' For subjects with a right-censored time-to-event, the value(s) in \code{mark} should be set to \code{NA}.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0 if placebo)
#' @param xlim a numeric vector of length 2 specifying the x-axis range (\code{NULL} by default)
#' @param ylim a numeric vector of length 2 specifying the y-axis range (\code{NULL} by default)
#' @param xtickAt a numeric vector specifing the position of x-axis tickmarks (\code{NULL} by default)
#' @param xtickLab a numeric vector specifying labels for tickmarks listed in \code{xtickAt}. If \code{NULL} (default), the labels are determined by \code{xtickAt}.
#' @param ytickAt a numeric vector specifing the position of y-axis tickmarks (\code{NULL} by default)
#' @param ytickLab a numeric vector specifying labels for tickmarks listed in \code{ytickAt}. If \code{NULL} (default), the labels are determined by \code{ytickAt}.
#' @param xlab a character string specifying the x-axis label (\code{NULL} by default)
#' @param ylab a character string specifying the y-axis label (\code{NULL} by default)
#' @param axis.title.size a numeric number specifying the font size of axis title in the bottom panel (\code{15} by default)
#' @param axis.text.size a numeric number specifying the font size of axis text in the bottom panel (\code{14} by default)
#' @param legend.text.size a numeric number specifying the font size of legend text in the bottom panel (\code{11} by default)
#' @param txLab a character vector of length 2 specifying the placebo and treatment labels (in this order). The default labels are \code{placebo} and \code{treatment}.
#' @param txLab.size a numeric number specifying the size of placebo and treatment labels (\code{5} by default). 
#' @param boxplot.width a numeric number specifying the width of the boxplot (\code{0.8}) by default
#' @jitter.height a numeric number specifying the amount of vertical jitter (\code{NULL} by default: this means the jitter values will occupy 80% of the implied bins.  
#' @jitter.seed a numeric number setting the seed of R's random number generator for jitter (\code{0} by default)
#' @param title a character string specifying the plot title (\code{NULL} by default)
#' @param subtitle a character string specifying a plot subtitle (\code{NULL} by default)
#' @param title.size a numeric number specifying font of the plot title (\code{16} by default)

#' @param ... other arguments to be passed to plotting functions
#' @details
#'
#' @return a ggplot object. 
#'
#' @examples
#' n <- 500
#' tx <- rep(0:1, each=n/2)
#' tm <- c(rexp(n/2, 0.2), rexp(n/2, 0.2 * exp(-0.4)))
#' cens <- runif(n, 0, 15)
#' eventTime <- pmin(tm, cens, 3)
#' eventInd <- as.numeric(tm <= pmin(cens, 3))
#' mark <- ifelse(eventInd==1, c(rbeta(n/2, 2, 5), rbeta(n/2, 2, 2)), NA)
#' markRng <- range(mark, na.rm=TRUE)
#'
#' # fit a model with a univariate mark
#' fit <- sievePH(eventTime, eventInd, mark, tx)
#' sfit <- summary(fit, markGrid=seq(markRng[1], markRng[2], length.out=10))
#' print(ggplot.summary.sievePH(sfit, mark, tx))
#'
#' @seealso \code{\link{sievePH}}, \code{\link{sievePHipw}}, \code{\link{sievePHaipw}} and \code{\link{summary.sievePH}}
#'
#' @export
ggplot.summary.sievePH <- function(x, mark=NULL, tx=NULL, xlim=NULL, ylim=NULL, xtickAt=NULL, xtickLab=NULL, ytickAt=NULL, 
                                   ytickLab=NULL, xlab=NULL, ylab=NULL, axis.title.size = 15, axis.text.size = 14, legend.text.size = 10,
                                   txLab=c("Placebo", "Treatment"), txLab.size = 5, boxplot.width = 0.8, jitter.height = NULL, jitter.seed = 0, 
                                   title=NULL, title.size = 16, subtitle = NULL, subtitle.size = 10, ...){
  library(ggpubr)
  library(ggplot2)
  
  contrast <- names(x)[length(names(x))]
  if (is.null(xlab)){ xlab <- colnames(x[[contrast]])[1] }
  if (is.null(ylab)){ ylab <- switch(colnames(x[[contrast]])[2], TE="Treatment Efficacy", HR="Hazard Ratio", LogHR="Log Hazard Ratio") }
  fit <- x[[contrast]]
  colnames(fit) <- c("mark","EST", "LB", "UB")
  if (is.null(xlim)){
    xlim <- range(fit[, 1])
  }
  
  if (is.null(ylim)){
    ylim <- range(fit[, -1], na.rm=TRUE)
  }

  p1 <- ggplot(fit)+
    geom_line(aes(x = mark, y = EST), linetype = "solid", size = 1.6)+
    geom_line(aes(x = mark, y = LB, linetype = '95% Pointwise CI', size = '95% Pointwise CI'))+
    geom_line(aes(x = mark, y = UB, linetype = '95% Pointwise CI', size = '95% Pointwise CI'))+
    theme_bw()+
    scale_linetype_manual(name="", labels = c('95% Pointwise CI'), values = c('95% Pointwise CI'= "dashed"))+
    scale_size_manual(name="", labels = c('95% Pointwise CI'), values = c('95% Pointwise CI'= 1.2))+
    xlab (xlab)+
    ylab (ylab)+
    theme(legend.key.size = unit(0.65, "cm"),
          legend.margin=margin(grid::unit(0,"cm")),
          legend.text=element_text(size=legend.text.size),
          legend.position = c(0,0.02),
          legend.justification = c(0, 0),
          legend.key = element_blank(),
          legend.key.width = unit(1.4,"cm"),
          legend.background = element_blank(),
          axis.title.x = element_text(size = axis.title.size, margin = margin(t = 10 )),
          axis.title.y = element_text(size = axis.title.size, margin = margin(r = 10 )),
          axis.text.x = element_text(size = axis.text.size, colour = "black"),
          axis.text.y = element_text(size = axis.text.size, colour = "black"),
          plot.margin=unit(c(-0.5,1,1,1), "cm"))
  
  #Add x-axis limits, ticks and labels
  if(is.null(xtickAt) & is.null(xtickLab)){
    p1 <- p1 +scale_x_continuous(limits = xlim, breaks = scales::pretty_breaks(n = 5))
    
  }else if (!is.null(xtickAt) & is.null(xtickLab)){
    p1 <- p1 +scale_x_continuous(limits = xlim, breaks = xtickAt, labels = xtickAt)
    
  }else if (is.null(xtickAt) & !is.null(xtickLab)){
    p1 <- p1 +scale_x_continuous(limits = xlim, breaks = xtickLab, labels = xtickLab)
    
  }else {
    p1 <- p1 +scale_x_continuous(limits = xlim, breaks = xtickAt, labels = xtickLab)
  }
  
  #Add y-axis limits, ticks and labels
  if(is.null(ytickAt) & is.null(ytickLab)){
    p1 <- p1 +scale_y_continuous(limits = ylim, breaks = scales::pretty_breaks(n = 5))
    
  }else if (!is.null(ytickAt) & is.null(ytickLab)){
    p1 <- p1 +scale_y_continuous(limits = ylim, breaks = ytickAt, labels = ytickAt)
    
  }else if (is.null(ytickAt) & !is.null(ytickLab)){
    p1 <- p1 +scale_y_continuous(limits = ylim, breaks = ytickLab, labels = ytickLab)
    
  }else {
    p1 <- p1 +scale_y_continuous(limits = ylim, breaks = ytickAt, labels = ytickLab)
  }
  
  data <- data.frame("tx" = as.factor(tx), "mark" = mark)
  set.seed(jitter.seed)
  p2 <- ggplot(data)+
    geom_boxplot(aes(x = mark, y = tx), width = boxplot.width, fill = "gray80",color = "black", lwd = 0.8, outlier.shape = NA, data = data)
  if(is.null(jitter.height)){
    p2 <- p2 + geom_jitter(aes(x = mark, y = tx, shape = as.factor(tx), color = as.factor(tx)), alpha = 1,size = 1.7,  
                           width = 0, fill = "white",  stroke = 1, data = data)
  }else{
    p2 <- p2 + geom_jitter(aes(x = mark, y = tx, shape = as.factor(tx), color = as.factor(tx)), alpha = 1,size = 1.7, height = jitter.height, 
                width = 0, fill = "white",  stroke = 1, data = data)
  }
  p2 <- p2 +  scale_shape_manual(name = "", labels = c("0" = txLab[1], "1" = txLab[2]), values = c(21,24))+
    scale_color_manual(name="", values = c("red","blue"), breaks = c("1","0"))
  p2 <- p2 + geom_text(x = min(xlim)-0.07*(xlim[2]-xlim[1]), y = "0", label = txLab[1], size = txLab.size, hjust = 1) 
  p2 <- p2 + geom_text(x = min(xlim)-0.07*(xlim[2]-xlim[1]), y = "1", label = txLab[2], size = txLab.size, hjust = 1)
  
  if(is.null(xtickAt)){
    p2 <- p2 +scale_x_continuous(limits = xlim, breaks = scales::pretty_breaks(n = 5),labels = NULL)
  }else{
    p2 <- p2 +scale_x_continuous(limits = xlim, breaks = xtickAt, labels = NULL)
  }
  if(!is.null(subtitle)){
    p2 <- p2 + labs(subtitle = subtitle)
  }
  
  if(!is.null(title)){
    p2 <- p2 + labs(title = title)
    
  }
  # if(!is.null(title) & !is.null(subtitle)){
  #   p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(title, size = title.size, face = "bold", vjust = 1.2, hjust = 0.2))
  #   
  # }
  p2 <- p2 +
    xlab("")+
    ylab("")+
    scale_y_discrete(labels = c("1" = txLab[2], "0" = txLab[1]),breaks = c("1", "0"))+
    theme_bw()+
    guides(shape = "none")+
    guides(color="none")+
    theme(plot.title = element_text(size = title.size, hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size = subtitle.size, hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin=unit(c(1,1,-0.12,2.7), "lines")
    )+coord_cartesian(clip = "off")
 
  
  
  
  p <- ggpubr::ggarrange(p2,p1,heights = c(0.33,0.67),ncol=1,nrow=2,align = "v")
  #p <- p2+p1+patchwork::plot_layout(ncol=1,heights = c(1,2))
  
  
  
  
  
  
  return(p)
}

