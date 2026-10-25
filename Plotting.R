# Shared figure style. Colour is accompanied by line type or a text label.
library(ggplot2)

zone_colours <- c(red="#C35A54", yellow="#EBD89C", green="#4A9984")
direction_colours <- c("Standard better"="#B34440", "Internal better"="#237B65")

paper_theme <- function(base_size=10) {
  theme_minimal(base_size=base_size, base_family="sans") +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(colour="#E8E8E8", linewidth=.3),
          plot.title=element_text(face="bold",size=base_size+1,hjust=0),
          plot.subtitle=element_text(colour="#555555",size=base_size-1),
          axis.title=element_text(size=base_size),
          legend.position="top",legend.title=element_blank(),
          legend.key.width=grid::unit(1.6,"lines"),
          strip.text=element_text(face="bold"),
          strip.background=element_rect(fill="#F4F4F4",colour=NA),
          plot.margin=margin(7,10,6,7))
}

save_figure <- function(plot, path, width=7.2, height=3.6) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  ggsave(path,plot,width=width,height=height,device=grDevices::pdf,
         family="Helvetica",useDingbats=FALSE,bg="white",limitsize=FALSE)
  ggsave(sub("\\.pdf$",".png",path),plot,width=width,height=height,
         dpi=180,bg="white",limitsize=FALSE)
}

plot_comparison <- function(test, x=NULL, title=NULL, subtitle=NULL,
                            change=NULL, show_zones=TRUE) {
  p <- test$path
  if (is.null(x)) x <- seq_len(nrow(p))
  stopifnot(length(x)==nrow(p))
  df <- rbind(data.frame(x=x,value=pmax(p$logMminus/log(10),-2),direction="Standard better",episode=p$episode),
              data.frame(x=x,value=pmax(p$logMplus/log(10),-2),direction="Internal better",episode=p$episode))
  threshold <- log10(test$config$threshold)
  # Values below 10^-2 are drawn at the labelled floor; saved paths stay exact.
  top <- max(threshold*1.12,min(max(df$value),threshold*1.6),.5)
  plot <- ggplot(df,aes(x,value,colour=direction,linetype=direction,
                        group=interaction(direction,episode))) +
    geom_hline(yintercept=0,colour="#B8B8B8",linewidth=.25) +
    geom_hline(yintercept=threshold,colour="#5E5E5E",linetype="dashed",linewidth=.4) +
    geom_line(linewidth=.55) +
    scale_colour_manual(values=direction_colours) +
    scale_linetype_manual(values=c("Standard better"="solid","Internal better"="22")) +
    scale_y_continuous(breaks=c(-2,-1,0,1,2,3,4),
                       labels=c("<= -2","-1","0","1","2","3","4")) +
    coord_cartesian(ylim=c(-2,top),clip="on") +
    labs(title=title,subtitle=subtitle,x=if(inherits(x,"Date")) NULL else "Observation",
         y=expression(log[10]~"e-value")) + paper_theme()
  if (inherits(x,"Date")) plot <- plot + scale_x_date(
    breaks=seq(as.Date("2008-01-01"),as.Date("2024-01-01"),by="4 years"),date_labels="%Y")
  if (!is.null(change)) plot <- plot +
    geom_vline(xintercept=change,colour="#8C8C8C",linetype="dotted",linewidth=.4)
  if (nrow(test$reports)) {
    reports <- test$reports
    reports$x <- x[match(reports$end,p$time)]
    reports$xstart <- x[match(reports$start,p$time)]
    reports$direction <- ifelse(reports$color=="red","Standard better","Internal better")
    # Completed intervals are shaded retrospectively in their report colour.
    shading <- geom_rect(data=reports,aes(xmin=xstart,xmax=x,ymin=-Inf,ymax=Inf,fill=color),
                          inherit.aes=FALSE,alpha=.07,colour=NA,show.legend=FALSE)
    plot$layers <- c(list(shading),plot$layers)
    plot <- plot + scale_fill_manual(values=zone_colours) +
      geom_point(data=reports,aes(x=x,y=threshold,colour=direction),
                              inherit.aes=FALSE,size=1.7,shape=16)
  }
  if (show_zones && any(p$color!="yellow")) {
    zones <- data.frame(x=x,value=-1.95,color=p$color)
    # A narrow strip identifies the current evidence without filling the plot.
    plot <- plot + geom_point(data=zones,aes(x,value,fill=color),shape=22,size=1.2,
                              colour="transparent",inherit.aes=FALSE,show.legend=FALSE) +
      NULL
    if(!nrow(test$reports)) plot <- plot + scale_fill_manual(values=zone_colours)
  }
  plot
}

plot_heatmap <- function(cells, title=NULL) {
  model_order <- unique(cells$internal)
  cells$internal <- factor(cells$internal,levels=model_order)
  cells$standard <- factor(cells$standard,levels=rev(model_order))
  ggplot(cells,aes(internal,standard,fill=color)) +
    geom_tile(colour="white",linewidth=.4) +
    geom_text(aes(label=ifelse(color=="red","R",ifelse(color=="green","G",""))),
              size=2.1,colour="#333333") +
    scale_fill_manual(values=c(zone_colours,diagonal="#ECECEC"),
                      limits=c("red","yellow","green","diagonal"),
                      breaks=c("red","yellow","green"),
                      labels=c("Standard better","Inconclusive","Internal better"),drop=FALSE) +
    coord_equal() + facet_wrap(~level,nrow=1) +
    labs(title=title,x="Internal model",y="Standard model") +
    paper_theme(9) + theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1),
                          legend.position="bottom",axis.ticks=element_blank())
}
