#' Normalize Data of NanoString GeoMx Protein Assays
#' @description NanostringNorm provides flexible options of normalization methods for raw data (only ERCC-normalized) of NanoString GeoMx protein assays for downstream performance evaluation.
#' @description  Note that certain normalization methods may be combined in a serial manner to constitute an addtional normalization method that is not on the list of available options in the 'normalization.method' parameter.
#' @return a matrix or data.frame object, same class and data structure as the input of raw DSP data matrix.
#' @export
#' @importFrom EnvStats geoMean
#' @param data.raw a raw DSP data matrix that has only been ERCC-normalized. The raw data matrix may be either a matrix or data.frame object with each row for a protein and each column for a ROI.
#' @param control.AreaToPixelSquare a numerical vector of ROI areas, must be arranged in the same order of the columns of data.raw. If missing, the normalization methods concerning area cannot be performed.
#' @param control.Nuclei a numerical vector of ROI nuclei counts, must be arranged in the same order of the columns of data.raw. If missing, the normalization methods concerning nuclei count cannot be performed.
#' @param control.HK a character vector listing all the housekeeping proteins that appear in the row names of data.raw.
#' @param control.SignalToNoise a character vector listing all the negative control proteins that appear in the row names of data.raw.
#' @param normalization.method a character specifying which normalization method is applied to normalize raw DSP data. Available normalization methods include: "None", "Area_geoMean", "Nuclei_geoMean", "STNR_geoMean", "STNS_mean2SD", "HK_geoMean".
#' @param control.rm a logical value indicating whether the housekeep/negative controls should be removed from the output of normalized data matrix.
#' @examples
#' data(COVID.19.DSP)
#'
#' data.raw <- COVID.19.DSP$Expr.raw.DSP
#'
#' group.info <- COVID.19.DSP$Meta.DSP
#'
#' # normalize the raw data by area of each ROI:
#' data.norm <- DSPNorm(data.raw = data.raw,
#'                        control.AreaToPixelSquare = group.info$area,
#'                        control.HK = c("S6", "GAPDH", "Histone H3"),
#'                        control.SignalToNoise = c("Ms IgG1", "Ms IgG2a", "Rb IgG"),
#'                        normalization.method = "Area_geoMean",
#'                        control.rm = FALSE)
#'
#' # normalize the raw data by three housekeeping controls:
#' data.norm <- DSPNorm(data.raw = data.raw,
#'                        control.HK = c("S6", "GAPDH", "Histone H3"),
#'                        control.SignalToNoise = c("Ms IgG1", "Ms IgG2a", "Rb IgG"),
#'                        normalization.method = "HK_geoMean",
#'                        control.rm = FALSE)
DSPNorm <- function(data.raw, control.AreaToPixelSquare, control.Nuclei, control.HK, control.SignalToNoise, normalization.method, control.rm=TRUE){



  if(normalization.method=="None")
  {
    data.norm.f = data.raw

    if(control.rm){
      data.norm.f.ctrl.rm = data.norm.f[!row.names(data.norm.f)%in%c(control.SignalToNoise, control.HK), ]
    } else{
      data.norm.f.ctrl.rm = data.norm.f[ , ]
    }

  }


  # normalize data by Area with geomean of area included:
  if(normalization.method=="Area_geoMean")
  {
    data.norm.f = data.raw
    # control.AreaToPixelSquare.align = control.AreaToPixelSquare[ ,colnames(data.raw)]
    # control.AreaToPixelSquare must be a vector following the same order of the raw data columns exactly!
    # AreaToPixelSquare.avg = mean(as.numeric(control.AreaToPixelSquare))
    Area.geomean = geoMean(as.numeric(control.AreaToPixelSquare), na.rm = T)
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id]*(Area.geomean/control.AreaToPixelSquare[col_id])
    }

    if(control.rm){
      data.norm.f.ctrl.rm = data.norm.f[!row.names(data.norm.f)%in%c(control.SignalToNoise, control.HK), ]
    } else{
      data.norm.f.ctrl.rm = data.norm.f[ , ]
    }

  }



  # normalize data by nuclei, i.e., per cell:
  if(normalization.method=="Nuclei_geoMean")
  {
    data.norm.f = data.raw
    # control.AreaToPixelSquare.align = control.AreaToPixelSquare[ ,colnames(data.raw)]
    # control.AreaToPixelSquare must be a vector following the same order of the raw data columns exactly!
    # AreaToPixelSquare.avg = mean(as.numeric(control.AreaToPixelSquare))
    control.Nuclei.rm0 = as.numeric(control.Nuclei)
    if(sum(control.Nuclei.rm0 == 0) > 0) {
      cat("Warning: in scaling by nuclei counts, found ROIs with 0 count, replaced 0 with 1!\n")
      control.Nuclei.rm0[control.Nuclei.rm0 == 0] = 1
    }

    Nuclei.geomean = geoMean(control.Nuclei.rm0, na.rm = T)
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id]*(Nuclei.geomean/control.Nuclei.rm0[col_id])
    }

    if(control.rm){
      data.norm.f.ctrl.rm = data.norm.f[!row.names(data.norm.f)%in%c(control.SignalToNoise, control.HK), ]
    } else{
      data.norm.f.ctrl.rm = data.norm.f[ , ]
    }

  }




  # normalize data by HK:
  if(normalization.method=="HK_geoMean")
  {
    if(!all(c(control.SignalToNoise, control.HK)%in%row.names(data.raw)))
    {
      cat("Warning: in doing HK normalization, some HK-Ctrl feature(s) are NOT contained in your raw data!\n")
      return(NULL)
    }

    data.norm.f = data.raw
    #control.HK.align = control.HK[ ,colnames(data.raw)]
    control.HK.align = data.raw[control.HK, colnames(data.raw)]
    mean.geometric = rep(NA, dim(control.HK.align)[2])
    for(col_id in 1:dim(control.HK.align)[2])
    {
      mean.geometric[col_id] = geoMean(control.HK.align[ ,col_id], na.rm = T)
    }
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id]*(mean(mean.geometric, na.rm = T)/mean.geometric[col_id])
    }

    if(control.rm){
      data.norm.f.ctrl.rm = data.norm.f[!row.names(data.norm.f)%in%c(control.SignalToNoise, control.HK), ]
    } else{
      data.norm.f.ctrl.rm = data.norm.f[ , ]
    }

  }





  # normalize data by Signal/Noise by ratio:
  if(normalization.method=="STNR_geoMean")
  {
    if(!all(c(control.SignalToNoise, control.HK)%in%row.names(data.raw)))
    {
      cat("Warning: in doing SignalToNoise normalization, some SignalToNoise-Ctrl feature(s) are NOT contained in your raw data!\n")
      return(NULL)
    }

    data.norm.f = data.raw
    # control.SignalToNoise.align = control.SignalToNoise[ ,colnames(data.raw)]
    control.SignalToNoise.align = data.raw[control.SignalToNoise, colnames(data.raw)]
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id]/geoMean(control.SignalToNoise.align[ ,col_id], na.rm = T)
    }

    if(control.rm){
      data.norm.f.ctrl.rm = data.norm.f[!row.names(data.norm.f)%in%c(control.SignalToNoise, control.HK), ]
    } else{
      data.norm.f.ctrl.rm = data.norm.f[ , ]
    }
  }



  # normalize data by Signal/Noise by subtraction of mean +2*SD:
  if(normalization.method=="STNS_mean2SD")
  {
    if(!all(c(control.SignalToNoise, control.HK)%in%row.names(data.raw)))
    {
      cat("Warning: in doing SignalToNoise normalization, some SignalToNoise-Ctrl feature(s) are NOT contained in your raw data!\n")
      return(NULL)
    }

    data.norm.f = data.raw
    # control.SignalToNoise.align = control.SignalToNoise[ ,colnames(data.raw)]
    control.SignalToNoise.align = data.raw[row.names(data.raw)%in%control.SignalToNoise, colnames(data.raw)]
    for(col_id in 1:dim(data.raw)[2])
    {
      data.norm.f[ ,col_id] = data.raw[ ,col_id] - (mean(control.SignalToNoise.align[ ,col_id], na.rm = T) + 2*sd(control.SignalToNoise.align[ ,col_id], na.rm = T))
    }

    data.norm.f[data.norm.f <= 0] = 1.0

    if(control.rm){
      data.norm.f.ctrl.rm = data.norm.f[!row.names(data.norm.f)%in%c(control.SignalToNoise, control.HK), ]
    } else{
      data.norm.f.ctrl.rm = data.norm.f[ , ]
    }
  }




  return(data.norm.f.ctrl.rm)
}

