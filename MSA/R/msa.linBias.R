#'
#'
#' MSA linearity and bias test
#' @param data
#' @return msa.linBias
#' @export

# part                부품
# reference           기준값
# output              측정값
# data                사용데이터셋
# conf.level          신뢰수준 (default = 0.95)
# ylim                y축범위
# main                대제목(default = "Gage Linearity and Bias Study")
# sub                 소제목 (default = "My MSA project")

msa.linBias <- function(part, reference,output, data , conf.level = 0.95, ylim,
                        main = "Gage Linearity and Bias Study", sub = "My MSA project"){

  library(MSA)
  library(grid)
  library(gridExtra)
  library(dplyr)
  library(tidyr)
  library(qualityTools)

  if (is.data.frame(data)) {
    if (deparse(substitute(part)) %in% names(data)) {
      part <- deparse(substitute(part))
    }
    if (!(part %in% names(data))) {
      stop(part, "is not a valid column name for", deparse(substitute(data)))
    }
    if (deparse(substitute(reference)) %in% names(data)) {
      reference <- deparse(substitute(reference))
    }
    if (!(reference %in% names(data))) {
      stop(reference, "is not a valid column name for", deparse(substitute(data)))
    }
    if (deparse(substitute(output)) %in% names(data)) {
      output <- deparse(substitute(output))
    }
    if (!(output %in% names(data))) {
      stop(output, "is not a valid column name for", deparse(substitute(data)))
    }

    if (part %in% names(data)) {
      data[[part]] <- factor(data[[part]])
    }else {
      stop(part, "is not a valid column name for", deparse(substitute(data)))
    }
    if (reference %in% names(data)) {
      data[[reference]] <- factor(data[[reference]])
    }else {
      stop(reference, "is not a valid column name for", deparse(substitute(data)))
    }
  }else {
    stop("A data.frame object is needed as data argument")
  }

  partl <- nlevels(data[[part]])
  stanl <- nlevels(data[[reference]])
  n <- nrow(data) / partl

  uniq <- data.frame(unique(data[,c(part,reference)]))

  if(partl != nrow(uniq)){
    stop("One part need only one reference value. Check the data.")
  }

  reference_list <- list()

  for(i in 1:partl){

    name <- as.character(i)
    num <- as.numeric(as.character(filter(data, part ==i)[[reference]][1]))
    reference_list[[name]] <- num
  }

  result_frame <- matrix(nrow = partl , ncol =n)
  rownames(result_frame) <- levels(data[[part]])

  for(i in 1:partl){
    for(j in 1:n){
      result_frame[i,j] <- filter(data, part ==i)[[output]][j]
    }
  }

  reference_vector <- rapply(reference_list,c)
  linDesign <- gageLinDesign(ref = reference_vector, n = n)
  result_frame <- as.data.frame(result_frame)
  response(linDesign) <- result_frame

  resultLin <- gageLin(linDesign, stats = FALSE, plot = FALSE)

  #visualization
  msa.prepCanvas(main, sub)

  #vp.basic
  vp.basic <- viewport(name = "basic", layout = grid.layout(2,2, widths = c(0.6,0.4), heights = c(0.6,0.4)))
  pushViewport(vp.basic)

  #vp.plot (1,1)
  vp.plot <- viewport(name = "plot", layout.pos.row = 1, layout.pos.col = 1)
  pushViewport(vp.plot)

  g = nrow(resultLin@X[2])
  m = ncol(resultLin@Y)
  bias = resultLin@Y
  mbias = numeric(g)
  for (i in 1:g) bias[i, ] = resultLin@Y[i, ] - resultLin@X[i, 2]
  for (i in 1:g) mbias[i] = mean(as.numeric(bias[i, ]))
  if (missing(ylim))
    ylim = c(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE))

  temp_frame <- data.frame(ref = resultLin@X$Ref, mbias =mbias)

  y <- ggplot(temp_frame, aes(x = ref, y = mbias)) +
       geom_point(color= "red", shape = 18, size = 3) +
       ylim(ylim[1],ylim[2]) +
       xlab("Reference Values") + ylab("Bias")

  temp_frame2 <- data.frame(ref = resultLin@X$Ref, bias = as.vector(t(bias)))

  y <- y  +geom_point(data = temp_frame2, aes(x = ref, y = bias)) +
            geom_hline(yintercept = 0,linetype = 2)

  BIAS = numeric()
  ref = numeric()
  for (i in 1:g) {
    BIAS = c(BIAS, as.numeric(bias[i, ]))
    ref = c(ref, rep(resultLin@X$Ref[i], length = m))
  }
  pre <- predict.lm(resultLin@model, interval = "confidence", level = conf.level)
  uniq2 <- data.frame(unique(pre[,c(1:3)]))
  temp_frame3 <- data.frame(ref, pre)

  finalgraph <- y + geom_line(data = temp_frame3, aes(x = ref, y = fit)) +
                    geom_line(data = temp_frame3, aes(x = ref, y = lwr),col = "blue",linetype=2) +
                    geom_line(data = temp_frame3, aes(x = ref, y = upr),col = "blue",linetype=2)


  print(finalgraph, newpage = FALSE)

  popViewport()#vp.plot

  #vp.plot (2, 1)
  vp.gageLin <- viewport(name = "gageLin", layout.pos.row = 2, layout.pos.col = 1)
  pushViewport(vp.gageLin)

  sigma <- round(summary(resultLin@model)[["sigma"]],digits = 5)
  R <- round(summary(resultLin@model)[["r.squared"]],digits = 5)
  R2 <- round(summary(resultLin@model)[["adj.r.squared"]],digits = 5)

  grid.text(expression(bold(paste("Gage Linearity"))),
            x = 0.55, y = 0.95, just = c("center", "top"),gp = gpar(cex = 1))

  table <- round(summary(resultLin@model)[["coefficients"]], digits = 5)
  rownames(table) <- c("Constant", "Slope")
  colnames(table) <- c("Coef", "SE Coef", "t value", "PR(>|t|)")
  table <- grid.table(table)

  grid.text(expression(paste("Residual standard error: ")),
            x = 0.20, y = 0.2, just = c("center", "top"),gp = gpar(cex = 0.9))
  grid.text(sigma,x = 0.40, y = 0.2,just = c("left", "top"),gp = gpar(cex = 0.9))
  grid.text(expression(paste("Multiple R-squared: ")),
            x = 0.20, y = 0.1, just = c("center", "top"),gp = gpar(cex = 0.9))
  grid.text(R,x = 0.40, y = 0.1,just = c("left", "top"),gp = gpar(cex = 0.9))
  grid.text(expression(paste("Adjusted R-squared: ")),
            x = 0.65, y = 0.1, just = c("center", "top"),gp = gpar(cex = 0.9))
  grid.text(R2,x = 0.80, y = 0.1,just = c("left", "top"),gp = gpar(cex = 0.9))

  print(table, newpage = FALSE)

  popViewport()#vp.plot (2, 1) END

  #vp.names (1:2,2)
  vp.names <- viewport(name = "names", layout.pos.row = 1:2, layout.pos.col = 2,
                       layout = grid.layout(3,1, heights = c(0.5,0.1,0.4)))
  pushViewport(vp.names)

  #vp.legends (1,2,2)의 (1,1)
  vp.legends <- viewport(name = "legends", layout.pos.row = 1, layout.pos.col = 1)
  pushViewport(vp.legends)

  ci <- paste(conf.level*100,"% confidence interval", sep = "")
  grid.text(expression(bold(paste("Legends"))),
            x = 0.50, y = 0.9, just = c("center", "top"),gp = gpar(cex = 1.0))
  grid.lines(x = c(0.1, 0.4), y = c(0.7, 0.7), gp = gpar(lty = 1,lwd = 3))# Regression
  grid.text("Regreesion", x = 0.55, y = 0.71, just = c("left",  "center"), gp = gpar(cex = 0.8))
  grid.lines(x = c(0.1, 0.4), y = c(0.6, 0.6), gp = gpar(lty = 2,lwd = 3, col = "blue"))# confidence interval
  grid.text(ci, x = 0.55, y = 0.61, just = c("left",  "center"), gp = gpar(cex = 0.8))
  grid.text("\u2666", x = 0.25, y = 0.51, just = c("left",  "center"), gp = gpar(cex = 1.3, col = "red"))# Average Bias
  grid.text("Average Bias", x = 0.55, y = 0.51, just = c("left",  "center"), gp = gpar(cex = 0.8))
  grid.text("\u25cf", x = 0.25, y = 0.41, just = c("left",  "center"), gp = gpar(cex = 1.1, col = "black"))# Data
  grid.text("Data", x = 0.55, y = 0.41, just = c("left",  "center"), gp = gpar(cex = 0.8))

  popViewport()#vp.legends(1,2) END

  vp.biasTitle <- viewport(name = "biasTitle", layout.pos.row = 2, layout.pos.col = 1)
  pushViewport(vp.biasTitle)

  grid.text(expression(bold(paste("Gage Bias"))),
            x = 0.50, y = 0.4, just = c("center", "top"),gp = gpar(cex = 1.0))

  popViewport()

  #vp.bias (1,2,2)의 (3,1)
  vp.bias <- viewport(name = "bias", layout.pos.row = 3, layout.pos.col = 1)
  pushViewport(vp.bias)

  biastable <- data.frame(reference = c(round(resultLin@X$Ref,digits = 5),"Average"),
                          bias = round(c(mbias,mean(mbias)),digits = 5))
  biastable <- grid.table(biastable)

  print(biastable, newpage = FALSE)

  popViewport()#vp.bias(1,2,2)의 (1,3) END

  popViewport() #vp.names (1:2,2) END

  popViewport()#vp.basic END

}
