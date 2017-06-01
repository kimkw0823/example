#'
#'@param data, LSL = NA,USL = NA,Target = NA, group = 1,alpha = 0.05,f.na.rm = TRUE, f.main = "Six Sigma Capability Analysis Study", f.sub = ""
#'@return pca.normDist
#'@export

#data data$variable
#(LSL = NA) 하한
#(USL = NA) 상한
#(Target = NA) 타겟
#(whatGroup = 1) 부분군 갯수를 무엇으로 구할건지(1: 상수, 2: 데이터 셋 열)
#(group1 = 1) 몇개의 군이 있는지                (whatGroup이 1일 때 사용)
#(group2 = 1) 데이터 셋 열                      (whatGroup이 2일 때 사용)
#(f.na.rm = TRUE) NA 값 있을 시 제외하고 통계 계산
#(withinSd1 = 1) 부분군 갯수 2이상일 때) 군내 표준편차 추정 방식
#                     (1: R_bar, 2: S_bar, 3: 합동표준편차 방식)
#(withinSd2 = 1) 부분군 갯수 1일 때)     군내 표준편차 추정 방식
#                     (1: 이동범위의 평균 방식, 2: 이동범위의 중위수 방식, 3: MSSD의 제곱근 방식)
#(maintTitle) 제목
#(subTitle) 소제목

pca.normdist <- function(data, LSL = NA,USL = NA,Target = NA, whatGroup =1, group1 = 1, group2 = c(),
                         f.na.rm = TRUE, withinSd1 =1, withinSd2 = 1,
                         mainTitle = "Process Capability Analysis Graph", subTitle = "My PCA project") {

  library(qcc)
  library(grid)
  library(ggplot2)
  library(reshape2)
  if(is.na(Target)){
    stop("Target is needed.")
  }
  if(is.na(LSL) & is.na(USL)){
    stop("No specification limits provided.")
  }
  if(!(withinSd1 == 1 ) & !(withinSd1 ==2) & !(withinSd1 ==3)){
    stop("Select 1 or 2 or 3")
  }
  if(!(withinSd2 == 1 ) & !(withinSd2 ==2) & !(withinSd2 ==3)){
    stop("Select 1 or 2 or 3")
  }

  #부분군 갯수 정의
  if(whatGroup ==1){
    group <- group1
  }else{
    group <- nlevels(as.factor(group2))
  }
  #변수 정의 (meanST = 공정평균, N = 전체 데이터 갯수(N), n = 군내 표본 갯수(n), data : vector 형태로 전환)
  meanST <- mean(data, na.rm = f.na.rm)
  N <- length(data[!is.na(data)])
  n <- N / group
  data <- as.vector(data)

  split <- NA
  rowSplit <- NA
  colSplit <- NA
  sr <- list(NA)
  overallC4 <- NA
  sdOver <- NA
  sdW <- NA

  #부분군 크기가 2이상인 경우
  if(n >=2){
    #데이터 군별로 쪼개기
    split <- splitby(data = data, group = group)
    rowSplit <- nrow(split)
    colSplit <- ncol(split)
    #측정 데이터 평균, 범위, 표준편차
    split$mean <- apply(split, 1,mean)
    split$range <- apply(split[1:colSplit],1,rangeMinMax)
    split$sd <- apply(split[1:colSplit] ,1,sd)
    #각 값들의 평균 구하기
    sr <- list(mean = mean(split[,colSplit+1]),
               range = mean(split[,colSplit+2]),sd = mean(split[,colSplit+3]))

    #sdOver = OVERALL S.D., uc>c4 사용
    if(!is.na(uc[N,1])){
      overallC4 <- uc[N,1]
    }else{
      overallC4 <- 4*(N -1)/(4*N -3)
    }
    sdOver <- sdOverAll(split = split, rowSplit = rowSplit, colSplit = colSplit,
                        N = N, mean = sr$mean, na.rm = f.na.rm, c4 = overallC4)

    #군내 표준편차 추정(qcc sd.R, sd.S 방식 사용)
    if(withinSd1 ==1){                        # R 방식
      if(group ==1){
        sdW <- sr$range / uc[n,3]
      }else{
      sdW <- sd.R(data = split[1:colSplit])
      }
    }else if(withinSd1 ==2){                  # S 방식
      if(group ==1){
        sdW <- sr$sd / uc[n,1]
      }else{
      sdW <- sd.S(data = split[1:colSplit])}
    }else{                                    # 합동표준편차방식
      sp <- getSp(split = split, rowSplit= rowSplit, colSplit = colSplit, N = N, group = group)
      tempc4 <- group*(n-1)+1
      sdW <- sp/uc[tempc4,1] # uc>c4 사용
    }

  } # 부분군 크기가 2이상인 경우(end)

  #부분군 크기가 1인 경우 (N = k = group)
  if(n ==1){

    #데이터 군별로 쪼개기
    split <- splitby(data = data, group = group)
    names(split) <- "X1"
    rowSplit <- nrow(split)
    colSplit <- ncol(split)

    #측정 데이터 이동범위(moving range), 연속차(successive differences)
    split$MR <- c(NA,abs(split[2:rowSplit, colSplit] - split[1:(rowSplit-1), colSplit]))
    split$sucdiff <-c(NA,split[2:rowSplit, colSplit] - split[1:(rowSplit-1), colSplit])

    #각 값들의 평균 구하기
    sr <- list(mean = mean(split[,1]),MR = sum(split[,colSplit+1],na.rm = TRUE)/(group-1))

    #sdOver = OVERALL S.D., uc>c4 사용
    overallC4 <- uc[N,1]
    sdOver <- sdOverAll2(split = split, rowSplit = rowSplit, colSplit = colSplit, N = group,
                         mean = sr$mean, na.rm = f.na.rm, c4 = overallC4)


    #군내 표준편차 추정
    if(withinSd2 ==1){                          # 이동범위의 평균 방식
      sdW <- sr$MR/uc[2,3] # uc>d2
    }else if(withinSd2 ==2){                    # 이동범위의 중위수 방식
      sdW <- median(split$MR,na.rm = TRUE)/uc[2,5] # uc>d4
    }else{                                      # MSSD의 제곱근 방식
      disq <- sum((split$sucdiff)^2)
      sdwup <- sqrt(disq/(2*(rowSplit-1)))
      sdW <- sdwup / uc[rowSplit,1] # uc>c4
    }


  }# 부분군 크기가 1인 경우(end)

  #pp, cp, ppm
  cplist <- ppclist(data = data,LSL = LSL, USL = USL, meanST = meanST, sdTotal = sdW) #cp, cpu, cpl, cpk
  pplist <- tpplist(data = data,LSL = LSL, USL = USL, meanST = meanST, sdOver = sdOver) #pp, ppu, ppl, ppk
  obPPMlist <- obPerform(data = data, LSL = LSL, USL = USL, N = N) # 관측 성능
  overPPMlist <- expPerform_overall(data = data, LSL = LSL, USL = USL, mean = meanST, sdOver = sdOver) # 전체 기대 성능
  inPPMlist <- expPerform_ingroup(data = data, LSL = LSL, USL = USL, mean = meanST,sdIngroup = sdW) # 군내 기대 성능

  #결과창
  pca.prepCanvas(main = mainTitle, sub = subTitle)

  vp.visualize <- viewport(name = "visualize", layout = grid.layout(2, 2,widths = c(0.7,0.3),heights =  c(0.7, 0.3)))
  pushViewport(vp.visualize)

  # 그래프 그리기(1:1)
  vp.hist <- viewport(name = "hist", layout.pos.row = 1, layout.pos.col = 1)
  pushViewport(vp.hist)

  grid.text("Histogram & Density", y = 1, just = c("center", "top"))

  binwST <- round(diff(range(data))/sqrt(N),digits = 3)
  ggdata <- melt(data)
  range(data)[2]- range(data)[1]
  qqp <- ggplot(ggdata, aes(x = data))+
    xlab("")+
    ylab("")# histogram 배경화면

  hist <- qqp + geom_histogram(aes(y = ..density..), binwidth = binwST, fill = "royalblue3",
                               colour = "black",stat = "bin")

  xST_density <- density(data, bw = binwST)
  if (!is.na(LSL)) {
    hist <- hist + annotate(geom = "text", x = LSL, y = max(xST_density$y),
                            label = "LSL", hjust = -0.1, size = 5)
  }
  hist <- hist + annotate(geom = "text", x = Target, y = max(xST_density$y),
                          label = "Target", hjust = -0.1, size = 5)
  if (!is.na(USL)) {
    hist <- hist + annotate(geom = "text", x = USL, y = max(xST_density$y),
                            label = "USL", hjust = 1.1, size = 5)
  }
  #hist <- hist  + theme(axis.text.y = element_blank())
  if (!is.na(LSL)) {
    hist <- hist + geom_vline(xintercept = LSL, linetype = 2,
                              size = 1)
  }
  if (!is.na(USL)) {
    hist <- hist + geom_vline(xintercept = USL, linetype = 2,
                              size = 1)
  }
  meandata <- melt(split$mean)
  hist <- hist + geom_vline(xintercept = Target,linetype = 3,size = 1) +
     stat_function(fun = dnorm,args = with(ggdata,c(mean(value), sd(value))),linetype = 1,size = 1)


  print(hist, newpage = FALSE)
  popViewport() # 그래프 그리기(1:1) END

  #오른쪽 표(1,2) 그리기
  vp.basicValue <- viewport(name = "basic value", layout.pos.row = 1, layout.pos.col = 2,layout = grid.layout(2, 2))
  pushViewport(vp.basicValue)

  grid.rect(gp = gpar(col = "gray0", lwd = 2)) # 사각형 틀 그리기

  #오른쪽 표(1,2)의 (1,1:2)부분
  vp.basicTitle <- viewport(name = "basic title", layout.pos.row = 1, layout.pos.col = 1:2)
  pushViewport(vp.basicTitle)

  grid.text(expression(bold("Density Lines Legend \n   &Summary statistics")), y = 0.70, just = c("center", "top"))
  grid.lines(x = c(0.05, 0.3), y = c(0.4, 0.4), gp = gpar(lty = 1,lwd = 3))# Density ST
  grid.text("Overall", x = 0.35, y = 0.40, just = c("left",  "center"), gp = gpar(cex = 0.8))

  popViewport() #오른쪽 표(1,2)의 (1,1:2)부분 END

  #오른쪽 표(1,2)의 (2,1)부분
  vp.input <- viewport(name = "input", layout.pos.row = 2, layout.pos.col = 1)
  pushViewport(vp.input)

  grid.text(expression(bold("LSL: ")), x = 0.50, y = 0.9, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(LSL,x = 0.50, y = 0.9,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Target: ")), x = 0.50, y = 0.7, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(Target, x = 0.50, y = 0.7, just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("USL: ")), x = 0.50, y = 0.5, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(USL, x = 0.50, y = 0.5, just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold("Group: ")), x = 0.50, y = 0.3, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(group, x = 0.50, y = 0.3, just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#오른쪽 표(1,2)의 (2,1)부분 END

  #오른쪽 표(1,2)의 (2,2)부분
  vp.process <- viewport(name = "proc", layout.pos.row = 2, layout.pos.col = 2)
  pushViewport(vp.process)

  grid.text(expression(bold("Sample N: ")), x = 0.50, y = 0.9, # n(데이터)
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(N, x = 0.5, y = 0.90,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold("Mean: ")), x = 0.5, y = 0.7, #평균
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(sprintf("%.3f",round(meanST,digits = 3) ), x = 0.5, y = 0.7,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold("StDev(Within): ")),x = 0.5, y = 0.5, # 군내 표준편차(within s.d.)
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(sprintf("%.3f", sdW), x = 0.5, y = 0.5,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold("StDev(Overall): ")), x = 0.5, y = 0.30, # 전체 표준편차(overall s.d.)
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(sprintf("%.3f", sdOver), x = 0.5, y = 0.30,
            just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#오른쪽 표(1,2)의 (2,2)부분 END

  popViewport()#오른쪽 표(1,2) 그리기 END

  # 아래 결과 표(2,1:2) 그리기
  vp.result <- viewport(name = "result", layout.pos.row = 2, layout.pos.col = 1:2, layout = grid.layout(1,2,widths = c(0.6,0.4)))
  pushViewport(vp.result)

  grid.rect(gp = gpar(col = "gray0", lwd = 2))  # 사각형 틀 그리기

  # 아래 왼쪽 (2:1)
  vp.perform <- viewport(name = "perform", layout.pos.col = 1,layout = grid.layout(1,4))
  pushViewport(vp.perform)

  grid.rect(gp = gpar(col = "gray0", lwd = 2))  # 사각형 틀 그리기

  grid.text(expression(bold("Performance")), y = 0.95, just = c("center", "top"))

  #아래 왼쪽 (2:1) 4구역 1번째
  vp.perform1 <- viewport(name = "perform1", layout.pos.col = 1)
  pushViewport(vp.perform1)

  grid.text(expression("PPM < st. lower limit"), x = 0.08, y = 0.60,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression("PPM < st. upper limit"), x = 0.08, y = 0.40,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression("PPM < st. total sum"),x = 0.08, y = 0.20,
            just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#아래 왼쪽 (2:1) 4구역 1번째 END

  #아래 왼쪽 (2:1) 4구역 2번째
  vp.perform2 <- viewport(name = "perform1", layout.pos.col = 2)
  pushViewport(vp.perform2)

  grid.text(expression("Observation"), x = 0.35, y = 0.75,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(obPPMlist$ppm_down, x = 0.5, y = 0.60,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(obPPMlist$ppm_up, x = 0.5, y = 0.40,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(obPPMlist$ppm_total, x = 0.5, y = 0.20,
            just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#아래 왼쪽 (2:1) 4구역 2번째 END

  #아래 왼쪽 (2:1) 4구역 3번째
  vp.perform3 <- viewport(name = "perform1", layout.pos.col = 3)
  pushViewport(vp.perform3)

  grid.text(expression(" Expected\n perform(total)"), x = 0.40, y = 0.75,
            just = c("left", "top"), gp = gpar(cex = 0.7))
  grid.text(round(overPPMlist$ppm_down,digits=3), x = 0.5, y = 0.60,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(round(overPPMlist$ppm_up,digits=3), x = 0.5, y = 0.40,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(round(overPPMlist$ppm_total,digits=3), x = 0.5, y = 0.20,
            just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#아래 왼쪽 (2:1) 4구역 3번째 END

  #아래 왼쪽 (2:1) 4구역 4번째
  vp.perform4 <- viewport(name = "perform1", layout.pos.col = 4)
  pushViewport(vp.perform4)

  grid.text(expression(" Expected\n perform(within)"), x = 0.30, y = 0.75,
            just = c("left", "top"), gp = gpar(cex = 0.7))
  grid.text(round(inPPMlist$ppm_down,digits=3), x = 0.5, y = 0.60,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(round(inPPMlist$ppm_up,digits=3), x = 0.5, y = 0.40,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(round(inPPMlist$ppm_total,digits=3), x = 0.5, y = 0.20,
            just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#아래 왼쪽 (2:1) 4구역 4번째 END

  popViewport()# 아래 왼쪽 (2:1) END


  #Indices(2:2)
  vp.indices <- viewport(name = "indi", layout.pos.col = 2, layout = grid.layout(1, 2))
  pushViewport(vp.indices)

  grid.rect(gp = gpar(col = "gray0", lwd = 2)) # 사각형 틀 그리기

  grid.text(expression(bold("Indices")), y = 0.95, just = c("center", "top"))

  #Indices left
  vp.indiLeft <- viewport(name = "indiLeft", layout.pos.col = 1)
  pushViewport(vp.indiLeft)

  grid.text(expression(bold(C[p] * ": ")), x = 0.60, y = 0.75, just = c("right", "top"), gp = gpar(cex = 0.8)) #cp
  grid.text(sprintf("%.4f", cplist$cp), x = 0.60, y = 0.75, just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold(C[pl] * ": ")), x = 0.60, y = 0.60, just = c("right", "top"),gp = gpar(cex = 0.8)) #cpl
  grid.text(sprintf("%.4f", cplist$cpl), x = 0.60, y = 0.60, just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold(C[pu] * ": ")), x = 0.60, y = 0.45, just = c("right", "top"), gp = gpar(cex = 0.8)) #cpu
  grid.text(sprintf("%.4f", cplist$cpu), x = 0.60, y = 0.45, just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold(C[pk] *": ")), x = 0.60, y = 0.30, just = c("right", "top"), gp = gpar(cex = 0.8)) #cpk
  grid.text(sprintf("%.4f", cplist$cpk), x = 0.60, y = 0.30, just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#Indices left END

  #Indices right
  vp.indiRight <- viewport(name = "indiRight", layout.pos.col = 2)
  pushViewport(vp.indiRight)

  grid.text(expression(bold(P[p] * ": ")), x = 0.40, y = 0.75, just = c("right", "top"), gp = gpar(cex = 0.8)) #cp
  grid.text(sprintf("%.4f", pplist$pp), x = 0.40, y = 0.75, just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold(P[pl] * ": ")), x = 0.40, y = 0.60, just = c("right", "top"),gp = gpar(cex = 0.8)) #cpl
  grid.text(sprintf("%.4f", pplist$ppl), x = 0.40, y = 0.60, just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold(P[pu] * ": ")), x = 0.40, y = 0.45, just = c("right", "top"), gp = gpar(cex = 0.8)) #cpu
  grid.text(sprintf("%.4f", pplist$ppu), x = 0.40, y = 0.45, just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold(P[pk] *": ")), x = 0.40, y = 0.30, just = c("right", "top"), gp = gpar(cex = 0.8)) #cpk
  grid.text(sprintf("%.4f", pplist$ppk), x = 0.40, y = 0.30, just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#Indices right END

  popViewport()#Indices(2:2) END

  popViewport()#아래 결과 표(2,1:2) 그리기 END

}


