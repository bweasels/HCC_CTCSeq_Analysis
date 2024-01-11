homogenizePlots <- function(plot){
  plot <- plot + theme_bw() + theme(panel.grid = element_blank(), 
                                    panel.border = element_blank(), 
                                    axis.line=element_line(color='black'), 
                                    axis.text = element_text(size=12),
                                    axis.title= element_text(size=14))
  return(plot)
}