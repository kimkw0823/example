#'
#'
#' MSA gage R&R test
#'@param data
#'@return msa.gageRR
#'@export

# var                 측정값
# part                부품
# appr                측정자
# lsl                 하한
# usl                 상한
# sigma               표준편차 (default =  6)
# data                활용데이터셋
# main                대제목 (default = "Gage R&R Study")
# sub                 소제목 (default = "My MSA project")
# alphaLim            유의수준 (default =0.05)
# digits              최소단위 (default = 4)

msa.gageRR <- function(var, part, appr, lsl = NA, usl = NA, sigma = 6, data,
                       main = "Gage R&R Study", sub = "My MSA project",
                       alphaLim = 0.05, digits = 4){

  library(MSA)
  library(grid)
  library(gridExtra)
  library(gplots)
  library(ggplot2)
  library(lattice)
  library(SixSigma)

  if (is.data.frame(data)) {
    if (deparse(substitute(var)) %in% names(data)) {
      var <- deparse(substitute(var))
    }
    if (!(var %in% names(data))) {
      stop(var, "is not a valid column name for", deparse(substitute(data)))
    }
    if (deparse(substitute(part)) %in% names(data)) {
      part <- deparse(substitute(part))
    }
    if (deparse(substitute(appr)) %in% names(data)) {
      appr <- deparse(substitute(appr))
    }
    if (part %in% names(data)) {
      data[[part]] <- factor(data[[part]])
    }
    else {
      stop(part, "is not a valid column name for", data)
    }
    if (appr %in% names(data)) {
      data[[appr]] <- factor(data[[appr]])
    }
    else {
      stop(appr, "is not a valid column name for", data)
    }
  }

  msa.gageRR2(var = var, part = part, appr = appr, lsl = lsl, usl = usl, sigma = sigma, data = data,
              main = main, sub = sub,
              alphaLim = alphaLim, digits = digits)

  ss.rr(var =var, part = part, appr = appr, lsl = lsl, usl = usl, sigma = sigma, data = data, main = main, sub = sub)


}
