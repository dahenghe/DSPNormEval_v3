#' Evaluate Performance of Normlaization Methods by Comparing Biological Variation
#' @description The function generates PCA plots of normalized DSP data of the ROIs that are treated as biological replicates by a list of normalization methods. A better normalization method is expected to yield larger biological variations between two known biological subgroups.
#' @return ggplot2 object of PCA plot.
#' @export
#' @import graphics stats ggplot2 ggfortify
#' @param normed.data.list a list object containing the normalized DSP expression matrix under all normalization methods, each normalized matrix can be either a data.frame or matrix.
#' @param norm.methods.list a numerical vector of the method IDs to be compared, for instance, norm.methods.list = 1:17.
#' @param grouping.info a character vector describing grouping info of each ROI, must be arranged in the same order of the columns of DSP expression matrix. Each unique character in the vector is to be treated as a unique biological subgroup. For optimal visualization, it is strongly recommended to consider only two unique biological groups of interest in each round of comparison.
#' @param main title to be added to the output figure.
#' @param title.size font size of title.
#' @examples
#' data(COVID.19.DSP)
#'
#' data.normed.list <- COVID.19.DSP$Expr.normed.DSP
#'
#' group.info <- COVID.19.DSP$Meta.DSP
#'
#' controls <- c("Ms IgG1", "Ms IgG2a", "Rb IgG",
#'               "S6", "GAPDH", "Histone H3")
#'
#' pick.id.1 <- group.info$organ.short=="Alveolar"
#'
#' pick.id.2 <- group.info$patient=="D9"
#'
#' group.info.pick <- group.info[pick.id.1 & pick.id.2, ]
#'
#' norm.data.list.s <- list()
#' for (method_id in c(1:5, 8:16)) {
#'   data.pick <- log(data.normed.list[[method_id]])
#'   norm.data.list.s[[method_id]] <-
#'      data.pick[!row.names(data.pick)%in%controls, pick.id.1 & pick.id.2]
#' }
#'
#' DSPNorm.eval.bio.var(normed.data.list = norm.data.list.s,
#'                      norm.methods.list = c(1:5, 8:16),
#'                      grouping.info = group.info.pick[ ,"PanCK.Status"],
#'                      main = "Sample D9",
#'                      title.size = 6.5)
DSPNorm.eval.bio.var <- function(normed.data.list, norm.methods.list, grouping.info, main="", title.size=6){

  # library(ggplot2)
  # library("ggfortify")
  # library(gridExtra)


  plots.list = list()
  plot_id = 1


  # collect PCA plots for all normalization methods on the same page:
  if(length(unique(grouping.info)) > 1){
    for (method_id in norm.methods.list) {

      data.scaled.bg_corrected.normed.pt.i = normed.data.list[[method_id]]


      for.pca = as.data.frame(t(data.scaled.bg_corrected.normed.pt.i))

      for.pca$Group = grouping.info

      library(ggfortify)  # After loading {ggfortify}, you can use ggplot2::autoplot function for stats::prcomp and stats::princomp objects.

      myplot.keep = ggplot2::autoplot(stats::prcomp(for.pca[ ,-dim(for.pca)[2]]), data=for.pca, label=F, shape = T, colour = 'Group', frame=T, frame.type = 'norm')
      myplot.keep = myplot.keep + ggplot2::theme(axis.text=ggplot2::element_text(size=5),axis.title=ggplot2::element_text(size=5), legend.position = "none", plot.title = ggplot2::element_text(size=title.size, face = "bold"), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
      if(main==""){
        myplot.keep = myplot.keep + ggplot2::ggtitle(paste(method_id, sep = ""))
      } else{
        myplot.keep = myplot.keep + ggplot2::ggtitle(paste(main, ", ", method_id, sep = ""))
      }


      plots.list[[plot_id]] = myplot.keep


      plot_id = plot_id + 1

    }

  } else{
    cat("   Found only one group, PCA plot is therefore ignored!\n")
  }


return(plots.list)
}











