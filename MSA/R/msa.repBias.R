#'
#'
#' MSA Reapeatability and Bias test
#'@param data
#'@return msa.repBias
#'@export
#'
#' @export

# data          데이터
# tol           공차 (LSL - USL)
# LSL           하한
# USL           상한
# target        목표값
# resol         최소단위, 분해능, 해(resolution)
# alpha         유의수준 (default = 0.05)
# mainTitle     대제목 (default = "Measurment System Analysis Graph")
# subTitle      소제목 (default = "My MSA project")

msa.repBias <- function(data, tol, LSL, USL, resol, target, alpha = 0.05
                      ,mainTitle = "Gage Reapeatability and Bias Test"
                      ,subTitle = "My MSA project"){

  options(scipen=999) # do not show in scientific notation

  # get lowest decimal place in data
  decimalplaces <- function(x) {
    if ((x %% 1) != 0) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }
  if(missing(resol)){

    low_deci <- decimalplaces(data[1])
    for(i in 2: length(data)){
      temp_deci <- decimalplaces(data[i])
      if(low_deci < temp_deci){
        low_deci <- temp_deci
      }
    } # get lowest decimal place in data END
    resol <- 10 ^(-1*(low_deci))
  }else{
    low_deci <- decimalplaces(resol)
  }

  if(missing(tol) & missing(LSL) & missing(USL) ){ # lsl, usl, tol 모두 없을때

    stop(print("Input missing variables."))
  }

  if(missing(tol) & (missing(LSL)|missing(USL)) ){ # lsl, usl 한쪽만 있고 tol도 없을때

    tol <- NA
    cg <- NA
    cgk <- NA

    #type1 study
    meanx <- round(mean(data),digits = low_deci)
    sdx <- round(sd(data),digits = low_deci)

    bias <- round(meanx - target, digits = low_deci)

    passQ <- TRUE # 합격 1 불합격 0
    if(!missing(USL)){
      standard <- USL - 4*sdx
      for(i in 1:length(data)){
        if(data[i] > standard) passQ <- FALSE
      }
    }
    if(!missing(LSL)){
      standard <- LSL + 4*sdx
      for(i in 1:length(data)){
        if(data[i] < standard) passQ <- FALSE
      }
    }

    if(passQ){
      print("pass the type1 test")
    }else{
      print("fail the type1 test")
    }
    testT <- t.test(data, alternative = "two.sided",mu=target,conf.level = 1-alpha)[["statistic"]]
    tt <- qt(p = 1-alpha/2, df = length(data)-1)
    pv <- t.test(data, alternative = "two.sided",mu=target,conf.level = 1-alpha)[["p.value"]]

    if(abs(testT) < tt){
      print(paste("passed the bias test : testT = |",testT,"| < ",tt," = t(",length(data)-1,") (",1-alpha/2,"%)",sep = ""))
    }else{
      print(paste("failed the bias test : testT = |",testT,"| >= ",tt," = t(",length(data)-1,") (",1-alpha/2,"%)",sep = ""))
    }

  }else{ # lsl, usl 둘다 있거나 tol이 있을때

    #tolerance
    if(missing(tol)){
      tempdeci <- max(decimalplaces(USL),decimalplaces(LSL))
      tol <- round(USL - LSL,digits=tempdeci)
    }
    #percent RE
    percentRE <- (target / tol) * 100

    meanx <- round(mean(data),digits = low_deci)
    sdx <- round(sd(data),digits = low_deci)

    ##Type1 check(GM, Bosch)

    #bias
    bias <- round(meanx - target,digits = low_deci)

    #cg, cgk (반복성, 편의)
    cg <- round((0.2 * tol)/(6 * sdx),digits = low_deci)
    cgk <- round((0.1 * tol - abs(bias))/(3 * sdx),digits =low_deci)

    #judgment
    if(cg >= 1.33 & cgk >= 1.33){
      print(paste("passed the test : cg = ",cg," , cgk = ",cgk,sep = ""))
    }else{
      print(paste("failed the test : cg = ",cg," , cgk = " ,cgk,sep = ""))
    }

    ##AIAG process
    #repeatability test
    percentEV <- round((6 * sdx) / tol * 100,digits = 3)
    if(percentEV <= 30){
      print(paste("passed the repeatability test : ",percentEV,sep = ""))
    }else{
      print(paste("failed the repeatability test : ",percentEV,sep = ""))
    }
    #bias test
    testT <- t.test(data, alternative = "two.sided",mu=target,conf.level = 1-alpha)[["statistic"]]
    tt <- qt(p = 1-alpha/2, df = length(data)-1)
    pv <- t.test(data, alternative = "two.sided",mu=target,conf.level = 1-alpha)[["p.value"]]

    if(abs(testT) < tt){
      print(paste("passed the bias test : testT = |",testT,"| < ",tt," = t(",length(data)-1,") (",1-alpha/2,"%)",sep = ""))
    }else{
      print(paste("failed the bias test : testT = |",testT,"| >= ",tt," = t(",length(data)-1,") (",1-alpha/2,"%)",sep = ""))
    }
  }
  #결과창
  msa.prepCanvas(main = mainTitle, sub = subTitle)

  vp.visualize <- viewport(name = "visualize", layout = grid.layout(2, 3, heights =  c(0.75, 0.25)))
  pushViewport(vp.visualize)

  # 그래프 그리기(1,1:3)
  vp.hist <- viewport(name = "hist", layout.pos.row = 1, layout.pos.col = 1:3)
  pushViewport(vp.hist)

  xlim <- c(0, length(data))
  histdata <- data.frame(x = seq(1:length(data)), y = data)
  hist <- ggplot(histdata, aes(x =x, y =y))+
    geom_point(shape = 20) + geom_line(color ="steelblue", size = 1)+
    geom_point(data = histdata, aes(x = x, y = y), col = "blue")+
    geom_hline(aes(yintercept= target, linetype = "Standard"), colour= "red")+
    labs(x="observations",y="")
  if(!is.na(tol)){
    hist <- hist +
      geom_hline(aes(yintercept= target+0.1*tol, linetype = ""), colour= 'red',lty = 2)+
      geom_hline(aes(yintercept= target-0.1*tol, linetype = ""), colour= 'red',lty = 2)
  }
  print(hist, newpage = FALSE)

  popViewport() # 그래프 그리기(1,1:3) END

  #result(2,1)
  vp.result1 <- viewport(name = "result1", layout.pos.row =2, layout.pos.col = 1)
  pushViewport(vp.result1)

  grid.text(expression(bold("Basic Statistics")), x = 0.55, y = 0.9, just = c("center", "top"),gp = gpar(cex = 0.9))
  grid.text(expression(bold("Target : ")), x = 0.50, y = 0.70, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(target,x = 0.60, y = 0.70,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Mean : ")), x = 0.50, y = 0.55, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(meanx,x = 0.60, y = 0.55,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("S.D. : ")), x = 0.50, y = 0.40, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(sdx,x = 0.60, y = 0.40,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("tolerance : ")), x = 0.50, y = 0.25, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(tol,x = 0.60, y = 0.25,just = c("left", "top"),gp = gpar(cex = 0.8))

  popViewport() #result(2,1) END

  #result(2,2)
  vp.result2 <- viewport(name = "result2", layout.pos.row =2, layout.pos.col = 2)
  pushViewport(vp.result2)

  grid.text(expression(bold("Bias")), x = 0.55, y = 0.9, just = c("center", "top"),gp = gpar(cex = 0.9))
  grid.text(expression(bold("Bias : ")), x = 0.50, y = 0.70, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(bias,x = 0.60, y = 0.70,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("test T : ")), x = 0.50, y = 0.55, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(round(testT,digits = low_deci),x = 0.60, y = 0.55,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("p value : ")), x = 0.50, y = 0.40, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(round(pv,digits = low_deci),x = 0.60, y = 0.40,just = c("left", "top"),gp = gpar(cex = 0.8))

  popViewport() #result(2,2) END

  #result(2,3)
  vp.result3 <- viewport(name = "result3", layout.pos.row =2, layout.pos.col = 3)
  pushViewport(vp.result3)

  grid.text(expression(bold("Process Capability")), x = 0.55, y = 0.9, just = c("center", "top"),gp = gpar(cex = 0.9))
  grid.text(expression(bold("Cg : ")), x = 0.50, y = 0.70, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(cg,x = 0.60, y = 0.70,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Cgk : ")), x = 0.50, y = 0.55, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(cgk,x = 0.60, y = 0.55,just = c("left", "top"),gp = gpar(cex = 0.8))

  popViewport() #result(2,3) END




}


