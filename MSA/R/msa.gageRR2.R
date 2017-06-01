#'
#'
#' MSA Type1 study analysis
#'@param data
#'@return msa.gageRR2
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
msa.gageRR2 <- function(var, part, appr, lsl = NA, usl = NA, sigma = 6, data,
                        main = "Gage R&R Study(ANOVA)", sub = "My MSA project",
                        alphaLim = 0.05, digits = 5,errorTerm = "interaction"){
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
  }else {
    stop("A data.frame object is needed as data argument")
  }
  a <- nlevels(data[[part]])
  b <- nlevels(data[[appr]])
  n <- nrow(data)/(a * b)
  options(show.signif.stars = FALSE)
  if (b == 1) {
    modelf <- as.formula(paste(var, "~", part))
    model <- aov(modelf, data = data)
    modelm <- summary(model)
    rownames(modelm[[1]])[2] <- "Repeatability"
    modelm[[1]] <- rbind(modelm[[1]], c(colSums(modelm[[1]][,1:2]), rep(NA, 3)))
    rownames(modelm[[1]])[3] <- "Total"
    cat("One-way ANOVA (single appraiser):\n\n")
    print(modelm)
    modelrm <- NULL
  }else {
    modelf <- as.formula(paste(var, "~", part, "*", appr))
    modelfm <- as.formula(paste(var, "~", part, "*", appr,
                                "+ Error(", part, "/", appr, ")"))
    model <- aov(modelf, data = data)
    modelm <- summary(model)
    if (errorTerm == "interaction") {
      modelm[[1]][1:2, 4] <- modelm[[1]][1:2, 3]/modelm[[1]][3,3]
      modelm[[1]][1:2, 5] <- pf(modelm[[1]][1:2, 4], modelm[[1]][1:2,1], modelm[[1]][3, 1], lower.tail = FALSE)
    }
    rownames(modelm[[1]])[4] <- "Repeatability"
    modelm[[1]] <- rbind(modelm[[1]], c(colSums(modelm[[1]][,1:2]), rep(NA, 3)))
    rownames(modelm[[1]])[5] <- "Total"
    cat("Complete model (with interaction):\n\n")
    print(modelm)
    cat("\nalpha for removing interaction:", alphaLim, "\n")
    pint <- modelm[[1]][3, 5]
    if (pint > alphaLim) {
      modelfr <- as.formula(paste(var, "~", part, "+", appr))
      modelr <- aov(modelfr, data = data)
      modelrm <- summary(modelr)
      rownames(modelrm[[1]])[3] <- "Repeatability"
      modelrm[[1]] <- rbind(modelrm[[1]], c(colSums(modelrm[[1]][,1:2]), rep(NA, 3)))
      rownames(modelrm[[1]])[4] <- "Total"
      cat("\n\nReduced model (without interaction):\n\n")
      print(modelrm)
    }
    else modelrm <- NULL
  }
  #b =1 일때 modelm 츨력, else (modelm 출력, pint>alphaLim일때 modelrm 출력)
  varComp <- matrix(ncol = 6, nrow = 7)
  rownames(varComp) <- c("Total Gage R&R", "  Repeatability","  Reproducibility",
                         paste0("    ", appr), paste0(part,":", appr), "Part-To-Part", "Total Variation")
  colnames(varComp) <- c("VarComp", "%Contrib", "StdDev", "StudyVar", "%StudyVar", "%Tolerance")
  if (b == 1) {
    varComp[2, 1] <- modelm[[1]][2, 3]
    varComp[4, 1] <- NA
    varComp[5, 1] <- NA
    varComp[3, 2] <- NA
    varComp[6, 1] <- max(c((modelm[[1]][1, 3] - modelm[[1]][2,3])/(b * n), 0))
    varComp[1, 1] <- varComp[2, 1]
    varComp[7, 1] <- varComp[1, 1] + varComp[6, 1]
  }else {
    if (pint > alphaLim) {
      varComp[2, 1] <- modelrm[[1]][3, 3]
      varComp[4, 1] <- max(c((modelrm[[1]][2, 3] - modelrm[[1]][3, 3])/(a * n), 0))
      varComp[5, 1] <- NA
      varComp[3, 1] <- varComp[4, 1]
      varComp[6, 1] <- max(c((modelrm[[1]][1, 3] - modelrm[[1]][3, 3])/(b * n), 0))
      varComp[1, 1] <- varComp[2, 1] + varComp[3, 1]
      varComp[7, 1] <- varComp[1, 1] + varComp[6, 1]
    }else {
      varComp[2, 1] <- modelm[[1]][4, 3]
      varComp[4, 1] <- max(c((modelm[[1]][2, 3] - modelm[[1]][3,3])/(a * n), 0))
      varComp[5, 1] <- max(c((modelm[[1]][3, 3] - modelm[[1]][4,3])/n, 0))
      varComp[3, 1] <- varComp[4, 1] + varComp[5, 1]
      varComp[6, 1] <- max(c((modelm[[1]][1, 3] - modelm[[1]][3, 3])/(b * n), 0))
      varComp[1, 1] <- varComp[2, 1] + varComp[3, 1]
      varComp[7, 1] <- varComp[1, 1] + varComp[6, 1]
    }
  }
  varComp[, 2] <- round(100 * (varComp[, 1]/varComp[7, 1]), 2)
  varComp[, 3] <- sqrt(varComp[, 1])
  varComp[, 4] <- varComp[, 3] * sigma
  varComp[, 5] <- round(100 * (varComp[, 3]/varComp[7, 3]), 2)
  varComp[, 6] <- round(100 * (varComp[, 4]/(usl - lsl)), 2)
  ncat <- max(c(1, floor((varComp[6, 4]/varComp[1, 4]) * 1.41)))
  if (b == 1) {
    varComp <- varComp[-c(3:5), ]
  }else {
    if (pint > alphaLim) {
      varComp <- varComp[-c(5), ]
    }
  }
  cat(paste("\nGage R&R\n\n"))
  print(varComp[, 1:2])
  cat("\n")
  if (!is.na(usl) && !is.na(lsl)) {
    print(varComp[, c(1, 3:6)])
  }else {
    print(varComp[, 3:5])
  }
  # if(!is.na(usl) && !is.na(lsl)) varComp[,1:6] else varComp[,1:2,3:5]
  cat(paste("\nNumber of Distinct Categories =", ncat, "\n"))


#
  ##그래프 그리기
  msa.prepCanvas(main, sub)


  vp.plots <- grid::viewport(name = "plots", layout = grid::grid.layout(2,2))
  grid::pushViewport(vp.plots)

  #vp.matrix(1:1)
  if(b==1 || (b!=1 && (pint <= alphaLim)) ){
    #vp.matrix(1,1:2)
    vp.matrix <- viewport(name = "matrix", layout.pos.row = 1, layout.pos.col = 1:2)
    pushViewport(vp.matrix)

    grid.text(expression(paste("Complete model (with interaction):")),
              x = 0.20, y = 0.90, just = c("center", "top"),gp = gpar(cex = 1.2))

    result.frame <- round(modelm[[1]],digits = digits)
    table <- grid.table(result.frame)
    print(table, newpage = FALSE)

    grid.text(expression(paste("alpha for removing interaction:")),
              x = 0.80, y = 0.15, just = c("center", "top"),gp = gpar(cex = 1))
    grid.text(alphaLim,x = 0.93, y = 0.15,just = c("left", "top"),gp = gpar(cex = 1, col = "red"))

    grid.lines(x = unit(c(0, 1), "npc"),
               y = unit(0.05, "npc"),gp=gpar(col = "azure3"))

    popViewport()#vp.matrix(1,1:2) END

  }else{
    #vp.matrix(1,1)
    vp.matrix <- viewport(name = "matrix", layout.pos.row = 1, layout.pos.col = 1)
    pushViewport(vp.matrix)

    grid.text(expression(paste("Complete model (with interaction):")),
              x = 0.40, y = 0.90, just = c("center", "top"),gp = gpar(cex = 1.1))

    result.frame <- round(modelm[[1]],digits = digits)
    table <- grid.table(result.frame)
    print(table, newpage = FALSE)

    grid.lines(x = unit(c(0, 1), "npc"),
               y = unit(0.05, "npc"),gp=gpar(col = "azure3"))

    popViewport()#vp.matrix(1,1) END

    #vp.matrix(1,2)
    vp.matrix2 <- viewport(name = "matrix2", layout.pos.row = 1, layout.pos.col = 2)
    pushViewport(vp.matrix2)

    grid.text(expression(paste("Reduced model (without interaction):")),
              x = 0.30, y = 0.90, just = c("center", "top"),gp = gpar(cex = 1.1))

    result.frame <-round(modelrm[[1]],digits = digits)
    table <- grid.table(result.frame)
    print(table, newpage = FALSE)

    grid.text(expression(paste("alpha for removing interaction:")),
              x = 0.65, y = 0.2, just = c("center", "top"),gp = gpar(cex = 1))
    grid.text(alphaLim,x = 0.90, y = 0.2,just = c("left", "top"),gp = gpar(cex = 1, col = "red"))

    grid.lines(x = unit(c(0, 1), "npc"),
               y = unit(0.05, "npc"),gp=gpar(col = "azure3"))

    popViewport()#vp.matrix(1,2) END
  }

  #vp.matrix3(2,1:2)
  vp.matrix3 <- viewport(name = "matrix3", layout.pos.row = 2, layout.pos.col = 1:2)
  pushViewport(vp.matrix3)

  grid.text(expression(paste("Gage R&R:")),
            x = 0.11, y = 1, just = c("center", "top"),gp = gpar(cex = 1.2))

  if(!is.na(usl) && !is.na(lsl)){

    result.frame2 <- round(varComp[,1:6],digits = digits)
    table2 <- grid.table(result.frame2)
    print(table2, newpage = FALSE)

  }else{

    result.frame2 <- round(varComp[,c(1:2,3:5)],digits = digits)
    table2 <- grid.table(result.frame2)
    print(table2, newpage = FALSE)
  }

  popViewport()#vp.matrix2(2:1) END


  grid.text(expression(paste("\nNumber of Distinct Categories: ")),
              x = 0.83, y = 0.025, just = c("center", "top"),gp = gpar(cex = 1))
  grid.text(ncat,x = 0.97, y = 0.025,just = c("left", "top"),gp = gpar(cex = 1, col = "red"))


  # grid.text(expression(paste("alpha for removing interaction:")),
  #           x = 0.80, y = 0.15, just = c("center", "top"),gp = gpar(cex = 1))
  # grid.text(alphaLim,x = 0.93, y = 0.15,just = c("left", "top"),gp = gpar(cex = 1, col = "red"))

  invisible(list(anovaTable = modelm,
                 anovaRed = modelrm,
                 varComp = varComp[,1:2],
                 studyVar = varComp[, 3:6], ncat = ncat))
}
