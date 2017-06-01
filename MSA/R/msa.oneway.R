#'
#'
#' MSA gageR&R oneway study analysis
#' @param data
#' @return msa.oneway
#' @export
#'
#' @export

# part                측정기
# result              측정값
# lsl                 하한
# usl                 상한
# data                사용데이터셋
# alpha               유의수준 (default = 0.05)
# sigma               표준편차 (default = 6)
# digits              최소단위 (default = 5)
# main                대제목(default = "Gage Oneway Test")
# sub                 소제목 (default = "My MSA project")

msa.oneway <- function(part, result, lsl, usl,data,
                      alpha = 0.05, sigma = 6, digits = 5,
                      main = "Gage Oneway Test", sub = "My MSA project"){

  library(MSA)
  library(grid)
  library(gridExtra)
  library(gplots)
  library(ggplot2)
  library(lattice)
  library(SixSigma)

  options(scipen=999) # do not show scientific notation
  options(show.signif.stars = FALSE)

  if (is.data.frame(data)) {
    if (deparse(substitute(result)) %in% names(data)) {
      result <- deparse(substitute(result))
    }
    if (!(result %in% names(data))) {
      stop(result, "is not a valid column name for", deparse(substitute(data)))
    }
    if (deparse(substitute(part)) %in% names(data)) {
      part <- deparse(substitute(part))
    }
    if (part %in% names(data)) {
      data[[part]] <- factor(data[[part]])
    }else {
      stop(part, "is not a valid column name for", data)
    }
  }else {
    stop("A data.frame object is needed as data argument")
  }


  a <- nlevels(data[[part]])
  n <- nrow(data)/a

  model_formula <- as.formula(paste(result, "~", part))
  model <- aov(model_formula, data = data)
  model_Anova <- summary(model)
  rownames(model_Anova[[1]])[2] <- "Repeatability"
  model_Anova[[1]] <- rbind(model_Anova[[1]], c(colSums(model_Anova[[1]][,1:2]), rep(NA, 3)))
  rownames(model_Anova[[1]])[3] <- "Total"

  #b =1 일때 model_Anova 츨력, else (model_Anova 출력, pint>alphaLim일때 modelrm 출력)
  varComp <- matrix(ncol = 6, nrow = 4)
  rownames(varComp) <- c("Total Gage R&R", "  Repeatability",
                         "Part-To-Part", "Total Variation")
  colnames(varComp) <- c("VarComp", "%Contrib", "StdDev", "StudyVar", "%StudyVar", "%Tolerance")

  varComp[1, 1] <- model_Anova[[1]][2,3]
  varComp[2, 1] <- model_Anova[[1]][2,3]
  varComp[3, 1] <- max(c((model_Anova[[1]][1, 3] - model_Anova[[1]][2,3])/n, 0))
  varComp[4, 1] <- varComp[1, 1] + varComp[3, 1]

  varComp[, 2] <- round(100 * (varComp[, 1]/varComp[4, 1]), digits = digits)
  varComp[, 3] <- sqrt(varComp[, 1])
  varComp[, 4] <- varComp[, 3] * sigma
  varComp[, 5] <- round(100 * (varComp[, 3]/varComp[4, 3]), digits = digits)

  if(!missing(lsl) & !missing(usl)){
    varComp[, 6] <- round(100 * (varComp[, 4]/(usl - lsl)), digits = digits)
  }else if(!missing(lsl)){
    varComp[, 6] <- round(100 * (varComp[, 4]/(2*(mean(data[[result]])-lsl))),digits = digits)
  }else if(!missing(usl)){
    varComp[, 6] <- round(100 * (varComp[, 4]/(2*(usl - mean(data[[result]])))),digits = digits)
  }

  ncat <- max(c(1, floor((varComp[3, 3]/varComp[1, 3]) * 1.41)))

  ##그래프 그리기
  msa.prepCanvas(main, sub)

  vp.plots <- grid::viewport(name = "plots", layout = grid::grid.layout(3,2))
  grid::pushViewport(vp.plots)

  #(1:1)
  vp.var <- grid::viewport(name = "variation", layout.pos.row = 1,layout.pos.col = 1)
  grid::pushViewport(vp.var)

  rowstoplot <- c(1,2,4)
  colstoplot <- c(2,5,6)
  klabels <- c("%Contribution", "%Study Var", "%Tolerance")

  databar <- varComp[rowstoplot, colstoplot]

  rownames(databar) <- c("G.R&R", "Repeat", "Part2Part")

  plot <- barchart(databar, freq = FALSE, grid = TRUE,
                   par.settings = list(axis.text = list(cex = 0.6),
                                       par.ylab.text = list(cex = 0.8),
                                       par.main.text = list(cex = 0.85)),
                   ylab = list("Percent",fontsize = 8),
                   panel = function(...){
                     panel.barchart(...)
                     panel.abline(h = 0)
                     panel.abline(h = c(10, 30), lty = 2, col = "gray") },
                   auto.key = list(text = klabels, cex = 0.8, columns = length(colstoplot),
                                   space = "bottom", rectangles = TRUE, points = FALSE,
                                   adj = 1, rep = FALSE),
                   stack = FALSE, horizontal = FALSE,
                   main = list("Components of Variation", fontsize = 14))
  print(plot, newpage = FALSE)
  grid::popViewport()#(1:1) END

  #(2:2)
  vp.varByPart <- grid::viewport(name = "varByPart", layout.pos.row = 2,layout.pos.col = 1)
  pushViewport(vp.varByPart)

  varByPart_plot <- stripplot(as.formula(paste(result, "~", part)),
                              data = data, grid = TRUE,
                              par.settings = list(axis.text = list(cex = 0.6),
                                                  par.xlab.text = list(cex = 0.8),
                                                  par.ylab.text = list(cex = 0.8),
                                                  par.main.text = list(cex = 0.9)),
                              main = paste(result, "by", part), type = c("p", "a"))
  print(varByPart_plot, newpage = FALSE)
  popViewport()#(2:2) END

  data.xrange <- aggregate(as.formula(paste(result, "~", part)),data = data, function(x) {max(x) - min(x)})
  data.xbar <- aggregate(as.formula(paste(result, "~", part)), data = data, mean)
  ar <- mean(data.xrange[[result]])

  #(2:2)
  vp.Xchart <- grid::viewport(name = "Xchart", layout.pos.row = 2,layout.pos.col = 2)
  pushViewport(vp.Xchart)

  xbar <- mean(data[[result]], na.rm = TRUE)
  ucl <- xbar + (3/(ss.cc.getd2(n) * sqrt(n))) * ar
  lcl <- xbar - (3/(ss.cc.getd2(n) * sqrt(n))) * ar
  glimits <- c(min(range(data.xbar[[result]])[1], lcl),
               max(range(data.xbar[[result]])[2],ucl)) + c(-1, 1) * 0.1 * diff(range(data.xbar[[result]]))
  plot <- xyplot(as.formula(paste(result, "~", part)),
                 data = data.xbar, pch = 16, par.settings = list(axis.text = list(cex = 0.6),
                                                                 par.xlab.text = list(cex = 0.8),
                                                                 par.ylab.text = list(cex = 0.8),
                                                                 par.main.text = list(cex = 0.9)),
                 par.strip.text = list(cex = 0.6),
                 main = expression(bold(bar(x) * " Chart by " * part)),
                 grid = TRUE, type = "b", ylim = glimits,
                 panel = function(...) {
                   panel.xyplot(...)
                   panel.abline(h = xbar, lty = 2)
                   panel.abline(h = ucl, col = "red3")
                   panel.abline(h = lcl, col = "red3")
                 })
  print(plot, newpage = FALSE)
  popViewport()#(2:2) END

  #(1:2)
  vp.Rchart <- grid::viewport(name = "Rchart", layout.pos.row = 1,layout.pos.col = 2)
  pushViewport(vp.Rchart)

  this.d3 <- ss.cc.getd3(n)
  this.d2 <- ss.cc.getd2(n)
  rlimits <- c(max(ar * (1 - 3 * (this.d3/(this.d2))), 0),
               ar * (1 + 3 * (this.d3/(this.d2))))
  glimits <- c(min(range(data.xrange[[result]])[1], rlimits[1]),
               max(range(data.xrange[[result]])[2], rlimits[2]))+ c(-1, 1) * 0.1 * diff(range(data.xrange[[result]]))
  plot <- xyplot(as.formula(paste(result, "~", part)),
                 data = data.xrange, pch = 16,
                 par.settings = list(axis.text = list(cex = 0.6),
                                     par.xlab.text = list(cex = 0.8),
                                     par.ylab.text = list(cex = 0.8),
                                     par.main.text = list(cex = 0.9)),
                 par.strip.text = list(cex = 0.6),
                 main = paste("R Chart by", part), grid = TRUE,
                 type = "b", ylim = glimits,
                 panel = function(...) {
                   panel.xyplot(...)
                   panel.abline(h = ar, lty = 2)
                   panel.abline(h = rlimits[1], col = "red3")
                   panel.abline(h = rlimits[2], col = "red3")})
  print(plot, newpage = FALSE)
  popViewport()#(1:2)

  #(3,1:2)
  vp.gagerr <- grid::viewport(name = "gagerr", layout.pos.row = 3,layout.pos.col = 1:2)
  pushViewport(vp.gagerr)

  result.frame <- round(varComp,digits = digits)
  table <- grid.table(result.frame)

  print(table, newpage = FALSE)

  popViewport()#(3,1:2)

  if(ncat >=5){ncatcol = "blue"}else{ncatcol = "red"}

  grid.text(expression(paste("\nNumber of Distinct Categories: ")),
            x = 0.75, y = 0.04, just = c("center", "top"),gp = gpar(cex = 1))
  grid.text(ncat,x = 0.90, y = 0.04,just = c("left", "top"),gp = gpar(cex = 1, col = ncatcol))

  invisible(list(anovaTable = model_Anova,
                 varComp = varComp[,1:2],
                 studyVar = varComp[, 3:6], ncat = ncat))



}
