#' Check Quality of Data of NanoString GeoMx Protein Assays
#' @description QC.DSP provides several measures to check quality of NanoString GeoMx Protein DSP data, revealing potential data quality issues before proceeding to downstream analysis, also providing important clues in deciding optimal normalization methods.
#' @return graphs for QC diagnosis.
#' @export
#' @import graphics stats
#' @param data.raw a raw DSP data matrix that has only been ERCC-normalized. The raw data matrix may be either a matrix or data.frame object with each row for a protein and each column for a ROI.
#' @param Controls a character vector listing both the negative control and housekeeping proteins that appear in the row names of data.raw, preferably in the order of negative controls first followed by housekeepers.
#' @param colData a data.frame object providing annotation info associated with the DSP data, such as area/nuclei count of each ROI, grouping info of each ROI, and subject to which each ROI belongs. The rows of colData and the columns of data.raw must be arranged in exactly the same order.
#' @param Comparison.Pair a character vector of length two, which can be found in the grouping info column of colData, as specified by the argument colData.grouping.index.
#' @param colData.grouping.index a character specifying the name of the column of colData containing the grouping info of the study.
#' @param colData.ROI.area.index a character specifying the name of the column of colData containing the area info of ROIs.
#' @param colData.ROI.nuclei.count.index a character specifying the name of the column of colData containing the nuclei count info of ROIs. If missing, then the QC steps concerning nuclei count are skipped.
#' @param colData.Subject.index a character specifying the name of the column of colData containing the subject info of ROIs. if missing, then single subject is assumed.
#' @examples
#' data(COVID.19.DSP)
#'
#' data.raw <- COVID.19.DSP$Expr.raw.DSP
#'
#' group.info <- COVID.19.DSP$Meta.DSP
#'
#' control.SignalToNoise.use <- c("Ms IgG1", "Ms IgG2a", "Rb IgG")
#' control.HK.use <- c("Histone H3", "S6", "GAPDH")
#' control.both <- c(control.SignalToNoise.use, control.HK.use)
#'
#' group.info.Alveolar <- group.info[group.info$organ.short=="Alveolar", ]
#'
#' panck.uniq <- unique(group.info.Alveolar$PanCK.Status)
#'
#' data.raw.Alveolar <- data.raw[ ,group.info$organ.short=="Alveolar"]
#'
#' QC.DSP(data.raw = data.raw.Alveolar,
#'        Controls = control.both,
#'        colData = group.info.Alveolar,
#'        Comparison.Pair = panck.uniq,
#'        colData.grouping.index = "PanCK.Status",
#'        colData.ROI.area.index = "area",
#'        colData.Subject.index = "patient")
QC.DSP <- function(data.raw, Controls, colData, Comparison.Pair, colData.grouping.index, colData.ROI.area.index, colData.ROI.nuclei.count.index=NULL, colData.Subject.index=NULL){

  # only use those obsn of the interested groups:
  compare.pick.IDs <- colData[ ,colData.grouping.index]%in%Comparison.Pair
  colData <- colData[compare.pick.IDs, ]
  data.raw <- data.raw[ ,compare.pick.IDs]


  colData$color.use <- "red"

  colData$color.use[colData[ ,colData.grouping.index]==Comparison.Pair[2]] <- "blue"


  if(!is.null(colData.Subject.index)){
    subject.uniq <- unique(colData[ ,colData.Subject.index])
  } else{

    cat("The colData.Subject.index is not provided, thus assuming single subject\n")

    subject.uniq <- NULL

  }



  # check area vs nuclei count relationship:
  if(!is.null(colData.ROI.nuclei.count.index)){

    if(is.null(subject.uniq)){
      par(mfcol=c(1,1))

      plot(x=colData[ ,colData.ROI.area.index], y=colData[ ,colData.ROI.nuclei.count.index], pch=16, cex=0.75, xlab = "Area", ylab = "Nuclei Count", main=paste("Assuming Single Subject", sep = ""), col=colData$color.use)
      nuclei.area <- lm(colData[ ,colData.ROI.nuclei.count.index] ~ colData[ ,colData.ROI.area.index])
      abline(nuclei.area, lwd=2, col="black")

      #plot.new()

    } else{

      par(mfrow=c(2,2))

      subjects.all <- colData[ ,colData.Subject.index]

      for (subject_id in 1:length(subject.uniq)) {

        colData.i <- colData[subjects.all==subject.uniq[subject_id], ]

        plot(x=colData.i[ ,colData.ROI.area.index], y=colData.i[ ,colData.ROI.nuclei.count.index], pch=16, cex=0.75, xlab = "Area", ylab = "Nuclei Count", main=paste("Subject: ", subject.uniq[subject_id], sep = ""), col=colData.i$color.use)
        nuclei.area <- lm(colData.i[ ,colData.ROI.nuclei.count.index] ~ colData.i[ ,colData.ROI.area.index])
        abline(nuclei.area, lwd=2, col="black")

      }

      #plot.new()

    }


  } else{
    cat("Nuclei count index not provided, the area vs nuclei relationship check is therefore ingnored\n")
  }


  # check confounding between group of interest and area:
  subtype <- rep(NA, length(unique(colData[ ,colData.grouping.index])))
  for(stg_id in 1:length(unique(colData[ ,colData.grouping.index])))
  {
    subtype.i <- unique(colData[ ,colData.grouping.index])[order(unique(colData[ ,colData.grouping.index]),decreasing = F)][stg_id]
    size.i <- length(which(colData[ ,colData.grouping.index]==subtype.i))
    subtype[stg_id] <- paste(subtype.i, "\n(n = ", size.i, ")", sep = "")
  }

  par(mfcol=c(1,1))

  myboxplot.DSP(y=colData[, colData.ROI.area.index],
                x=colData[ ,colData.grouping.index], title="", type = subtype,
                ylab = paste("Area", sep = ""))

  if(!is.null(colData.ROI.nuclei.count.index)){

    myboxplot.DSP(y=colData[, colData.ROI.nuclei.count.index],
                  x=colData[ ,colData.grouping.index], title="", type = subtype,
                  ylab = paste("Nuclei Count", sep = ""))

  }

  #plot.new()


  # check correlation between group of interest and HK/Negative controls:
  if(!is.null(colData.ROI.nuclei.count.index)){
    par(mfrow=c(3,3))
  } else{
    par(mfrow=c(3,2))
  }

  for(ctrl_id in 1:length(Controls)){

    expr.ctrl.i.raw <- as.numeric(data.raw[Controls[ctrl_id], ])
    expr.ctrl.i.raw[expr.ctrl.i.raw==0] <- min(expr.ctrl.i.raw[expr.ctrl.i.raw > 0], na.rm = T)

    myboxplot.DSP(y=log(expr.ctrl.i.raw),
                  x=colData[ ,colData.grouping.index], title="", type = subtype,
                  ylab = paste("Log(Raw ", Controls[ctrl_id], ")", sep = ""))

    expr.ctrl.i.unit.area <- expr.ctrl.i.raw/as.numeric(colData[, colData.ROI.area.index])

    myboxplot.DSP(y=log(expr.ctrl.i.unit.area),
                  x=colData[ ,colData.grouping.index], title="", type = subtype,
                  ylab = paste("Log(", Controls[ctrl_id], "/Area)", sep = ""))

    if(!is.null(colData.ROI.nuclei.count.index)){
      expr.ctrl.i.unit.nuclei <- expr.ctrl.i.raw/as.numeric(colData[, colData.ROI.nuclei.count.index])

      myboxplot.DSP(y=log(expr.ctrl.i.unit.nuclei),
                    x=colData[ ,colData.grouping.index], title="", type = subtype,
                    ylab = paste("Log(", Controls[ctrl_id], "/Nuclei)", sep = ""))
    }

  }
  #plot.new()


  # check correlation among HK and Neg. controls, as well as the controls vs area:
  # in this measure, we use log-transformed values
  data.raw.ctrls.all <- as.data.frame(log(t(data.raw[Controls, ])))

  if(is.null(subject.uniq)){

    # correlation among controls:
    par(mfrow=c(1,1))
    plot(data.raw.ctrls.all, pch=16, cex=0.75, main=paste("Assuming Single Subject", sep = ""), col=colData$color.use)
    #plot.new()

    # correlation between controls and area:
    par(mfcol=c(3,2))
    for (ctrl_id in 1:length(Controls)) {

      expr.ctrl.i <- as.numeric(data.raw[Controls[ctrl_id], ])

      cor.area <- cor.test(x=as.numeric(colData[ ,colData.ROI.area.index]), y=expr.ctrl.i)

      plot(x=colData[ ,colData.ROI.area.index], y=expr.ctrl.i,
           #xlim = c(0, max(group.info.pt.i$surface_area.num)),
           #ylim = c(0,max(ctrl.val)),
           xlab = "Area",
           pch=16,
           ylab = paste("Raw ", Controls[ctrl_id], sep = ""),
           main = paste("Assuming Single Subject", "\n", "Corr = ", signif(cor.area$estimate,3),
                        ", P-Value = ", signif(cor.area$p.value, 3), sep = ""), col=colData$color.use)
      model.area <- lm(expr.ctrl.i~as.numeric(colData[ ,colData.ROI.area.index]))
      abline(model.area, lwd=2, col="black")
      # abline(h=0, lty="dashed", lwd=1.5)
      # abline(v=0, lty="dashed", lwd=1.5)



    }
    #plot.new()

    # correlation between controls and nuclei counts:
    if(!is.null(colData.ROI.nuclei.count.index)){
      par(mfcol=c(3,2))
      for (ctrl_id in 1:length(Controls)) {

        expr.ctrl.i <- as.numeric(data.raw[Controls[ctrl_id], ])

        cor.nuclei <- cor.test(x=as.numeric(colData[ ,colData.ROI.nuclei.count.index]), y=expr.ctrl.i)

        plot(x=colData[ ,colData.ROI.nuclei.count.index], y=expr.ctrl.i,
             #xlim = c(0, max(group.info.pt.i$surface_area.num)),
             #ylim = c(0,max(ctrl.val)),
             xlab = "Nuclei Count",
             pch=16,
             ylab = paste("Raw ", Controls[ctrl_id], sep = ""),
             main = paste("Assuming Single Subject", "\n", "Corr = ", signif(cor.nuclei$estimate,3),
                          ", P-Value = ", signif(cor.nuclei$p.value, 3), sep = ""), col=colData$color.use)
        model.nuclei <- lm(expr.ctrl.i~as.numeric(colData[ ,colData.ROI.nuclei.count.index]))
        abline(model.nuclei, lwd=2, col="black")
        # abline(h=0, lty="dashed", lwd=1.5)
        # abline(v=0, lty="dashed", lwd=1.5)



      }
      #plot.new()

    }



  } else{

    subjects.all <- colData[ ,colData.Subject.index]

    for (subject_id in 1:length(subject.uniq)) {

      # correlation among controls:
      data.raw.ctrls.j <- data.raw.ctrls.all[subjects.all==subject.uniq[subject_id], ]

      colData.j <- colData[subjects.all==subject.uniq[subject_id], ]

      par(mfrow=c(1,1))
      plot(data.raw.ctrls.j, pch=16, cex=0.75, main=paste("Subject: ", subject.uniq[subject_id], sep = ""), col=colData.j$color.use)
      #plot.new()

      par(mfcol=c(3,2))
      for (ctrl_id in 1:length(Controls)) {

        expr.ctrl.i.j <- as.numeric(data.raw[Controls[ctrl_id], subjects.all==subject.uniq[subject_id]])

        cor.area <- cor.test(x=as.numeric(colData.j[ ,colData.ROI.area.index]), y=expr.ctrl.i.j)

        plot(x=colData.j[ ,colData.ROI.area.index], y=expr.ctrl.i.j,
             #xlim = c(0, max(group.info.pt.i$surface_area.num)),
             #ylim = c(0,max(ctrl.val)),
             xlab = "Area",
             pch=16,
             ylab = paste("Raw ", Controls[ctrl_id], sep = ""),
             main = paste("Subject: ", subject.uniq[subject_id], "\n", "Corr = ", signif(cor.area$estimate,3),
                          ", P-Value = ", signif(cor.area$p.value, 3), sep = ""), col=colData.j$color.use)
        model.area <- lm(expr.ctrl.i.j~as.numeric(colData.j[ ,colData.ROI.area.index]))
        abline(model.area, lwd=2, col="black")
        # abline(h=0, lty="dashed", lwd=1.5)
        # abline(v=0, lty="dashed", lwd=1.5)

      }
      #plot.new()

      if(!is.null(colData.ROI.nuclei.count.index)){
        par(mfcol=c(3,2))
        for (ctrl_id in 1:length(Controls)) {

          expr.ctrl.i.j <- as.numeric(data.raw[Controls[ctrl_id], subjects.all==subject.uniq[subject_id]])

          cor.nuclei <- cor.test(x=as.numeric(colData.j[ ,colData.ROI.nuclei.count.index]), y=expr.ctrl.i.j)

          plot(x=colData.j[ ,colData.ROI.nuclei.count.index], y=expr.ctrl.i.j,
               #xlim = c(0, max(group.info.pt.i$surface_area.num)),
               #ylim = c(0,max(ctrl.val)),
               xlab = "Nuclei Count",
               pch=16,
               ylab = paste("Raw ", Controls[ctrl_id], sep = ""),
               main = paste("Subject: ", subject.uniq[subject_id], "\n", "Corr = ", signif(cor.nuclei$estimate,3),
                            ", P-Value = ", signif(cor.nuclei$p.value, 3), sep = ""), col=colData.j$color.use)
          model.nuclei <- lm(expr.ctrl.i.j~as.numeric(colData.j[ ,colData.ROI.nuclei.count.index]))
          abline(model.nuclei, lwd=2, col="black")
          # abline(h=0, lty="dashed", lwd=1.5)
          # abline(v=0, lty="dashed", lwd=1.5)



        }
        #plot.new()

      }

    }

  }






}
