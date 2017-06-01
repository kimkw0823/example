#'
#'
#' MSA distructive Testing
#'@param data
#'@return msa.distTest
#'@export
#'
#' @export

# data 데이터
# tol = LSL - USL
# LSL = 하한
# USL = 상한
# resol = 최소단위
# target = 타겟
# alpha = (0.05)
# mainTitle = 대제목
# subTitle = 소제목

msa.destGageTable <- function (gdo, method="nested", sigma=6, alpha=0.05,tolerance ,dig = 5,
                               main = "Gage R&R table",sub = "My MSA project")
{

  library(ggplot2)
  library(gridExtra)
  library(qualityTools)
  if (method %in% c("crossed", "nested"))
    method = method
  else method = gdo@method
  yName = names(gdo)[5]
  aName = names(gdo)[3]
  bName = names(gdo)[4]
  if (method == "crossed")
    abName = paste(aName, ":", bName, sep = "")
  if (method == "nested")
    abName = paste(bName, "(", aName, ")", sep = "")
  bTobName = paste(bName, "to", bName, sep = " ")
  a = gdo@X[, aName]
  b = gdo@X[, bName]
  y = gdo@X[, yName]
  nestedFormula = as.formula(paste(yName, "~", aName, "/",
                                   bName))
  crossedFormula = as.formula(paste(yName, "~", aName, "*",
                                    bName))
  reducedFormula = as.formula(paste(yName, "~", aName, "+",
                                    bName))
  if (!is.null(tolerance))
    tolerance(gdo) = tolerance
  if (is.na(y) || !is.numeric(y))
    stop("Measurements need to be numeric")
  if (method == "nested") {
    numA <- nlevels(a[, drop = T])
    numB <- nlevels(b[, drop = T])
    numMPP <- length(y)/((numB) * numA)
    gdo@numO = numA
    gdo@numP = numB
    gdo@numM = numMPP
    fit = aov(nestedFormula, data = gdo)
    meanSq <- anova(fit)[, 3]
    gdo@ANOVA = fit
    gdo@method = "nested"
    MSa = meanSq[1]
    MSab = meanSq[2]
    MSe = meanSq[3]
    Cerror = MSe
    Cb = (MSab - MSe)/numMPP
    Ca = (MSa - MSab)/(numB * numMPP)
    if (Ca <= 0)
      Ca = 0
    if (Cb <= 0)
      Cb = 0
    Cab = 0
    totalRR = Ca + Cab + Cerror
    repeatability = Cerror
    reproducibility = Ca
    bTob = Cb
    totalVar = Cb + Ca + Cab + Cerror
    estimates = list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
    varcomp = list(totalRR = totalRR, repeatability = repeatability,
                   reproducibility = reproducibility, bTob = bTob, totalVar = totalVar)
    gdo@Estimates = estimates
    gdo@Varcomp = varcomp
  }
  if (method == "crossed") {
    numA <- nlevels(a[, drop = T])
    numB <- nlevels(b[, drop = T])
    numMPP <- length(a)/(numA * numB)
    gdo@numO = numA
    gdo@numP = numB
    gdo@numM = numMPP
    fit = aov(crossedFormula, data = gdo)
    model <- anova(fit)
    gdo@ANOVA = fit
    gdo@method = "crossed"
    MSb = MSa = MSab = MSe = 0
    if (bName %in% row.names(model))
      MSb = model[bName, "Mean Sq"]
    else warning(paste("missing factor", bName, "in model"))
    if (aName %in% row.names(model))
      MSa = model[aName, "Mean Sq"]
    else warning(paste("missing factor", aName, "in model"))
    if (abName %in% row.names(model))
      MSab = model[abName, "Mean Sq"]
    else warning(paste("missing interaction", abName, "in model"))
    if ("Residuals" %in% row.names(model))
      MSe = model["Residuals", "Mean Sq"]
    else warning("missing Residuals in model")
    Cb = Ca = Cab = Cerror = 0
    Cb = (MSb - MSab)/(numA * numMPP)
    Ca = (MSa - MSab)/(numB * numMPP)
    Cab = (MSab - MSe)/(numMPP)
    Cerror = (MSe)
    gdo@RedANOVA = gdo@ANOVA
    if ((Cab < 0) || (model[abName, "Pr(>F)"] >= alpha)) {
      redFit <- aov(reducedFormula, data = gdo)
      model <- anova(redFit)
      MSb = MSa = MSab = MSe = 0
      if (bName %in% row.names(model))
        MSb = model[bName, "Mean Sq"]
      else warning(paste("missing factor", bName, "in model"))
      if (aName %in% row.names(model))
        MSa = model[aName, "Mean Sq"]
      else warning(paste("missing factor", aName, "in model"))
      if ("Residuals" %in% row.names(model))
        MSe = model["Residuals", "Mean Sq"]
      else warning("missing Residuals in model")
      Cb = Ca = Cab = Cerror = 0
      Cb = (MSb - MSe)/(numA * numMPP)
      Ca = (MSa - MSe)/(numB * numMPP)
      Cab = 0
      Cerror = (MSe)
      gdo@RedANOVA = redFit
    }
    gdo@method = "crossed"
    Ca = max(0, Ca)
    Cb = max(0, Cb)
    Cab = max(0, Cab)
    totalRR = Ca + Cab + Cerror
    repeatability = Cerror
    reproducibility = Ca + Cab
    bTob = max(0, Cb)
    totalVar = Cb + Ca + Cab + Cerror
    estimates = list(Cb = Cb, Ca = Ca, Cab = Cab, Cerror = Cerror)
    varcomp = list(totalRR = totalRR, repeatability = repeatability,
                   reproducibility = reproducibility, a = Ca, a_b = Cab,
                   bTob = bTob, totalVar = totalVar)
    gdo@Estimates = estimates
    gdo@Varcomp = varcomp
  }
  cat("\n")
  cat(paste("AnOVa Table - ", gdo@method, "Design\n"))
  print(summary(gdo@ANOVA))
  cat("\n")
  cat("----------\n")
  if (!identical(gdo@RedANOVA, gdo@ANOVA) && gdo@method ==
      "crossed") {
    cat(paste("AnOVa Table Without Interaction - ", gdo@method,
              "Design\n"))
    print(summary(gdo@RedANOVA))
    cat("\n")
    cat("----------\n")
  }
  Source = names(gdo@Varcomp)
  Source[Source == "repeatability"] = " repeatability"
  Source[Source == "reproducibility"] = " reproducibility"
  Source[Source == "a_b"] = paste("  ", abName)
  Source[Source == "a"] = paste("  ", aName)
  Source[Source == "bTob"] = bTobName
  VarComp = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))]),
                  3)
  Contribution = round(as.numeric(gdo@Varcomp[c(1:length(gdo@Varcomp))])/as.numeric(gdo@Varcomp[length(gdo@Varcomp)]),
                       3)
  VarComp = t(data.frame(gdo@Varcomp))
  VarCompContrib = VarComp/gdo@Varcomp$totalVar
  Stdev = sqrt(VarComp)
  StudyVar = Stdev * gdo@Sigma
  StudyVarContrib = StudyVar/StudyVar["totalVar", ]
  SNR = 1
  ptRatio = NULL
  temp = NULL
  if ((length(gdo@GageTolerance) > 0) && (gdo@GageTolerance >
                                          0)) {
    ptRatio = StudyVar/gdo@GageTolerance
    temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar,
                      StudyVarContrib, ptRatio)
    names(temp)[6] = c("P/T Ratio")
    row.names(temp) = c(Source)
  }
  else {
    temp = data.frame(VarComp, VarCompContrib, Stdev, StudyVar,
                      StudyVarContrib)
    row.names(temp) = c(Source)
  }
  cat("\n")
  cat("Gage R&R\n")
  tempout = temp

  msa.prepCanvas(main = main, sub = sub)
  table <- grid.table(round(as.data.frame(tempout),digits = dig))
  print(table, newpage = FALSE)

  invisible(gdo)
}
