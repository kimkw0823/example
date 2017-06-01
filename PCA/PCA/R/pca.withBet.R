#'
#'@param data xLT= NA, LSL = NA,USL = NA,Target = NA, group = 1,alpha = 0.05, f.na.rm = TRUE, f.main = "Six Sigma Capability Analysis Study", f.sub = ""
#'@return pca.withBet
#'@export

#data data$variable
#(LSL = NA) 하한
#(USL = NA) 상한
#(Target = NA) 타겟
#(whatGroup = 1) 부분군 갯수를 무엇으로 구할건지(1: 상수, 2: 데이터 셋 열)
#(group1 = 1) 몇개의 군이 있는지                (whatGroup이 1일 때 사용)
#(group2 = 1) 데이터 셋 열                      (whatGroup이 2일 때 사용)
#(f.na.rm = TRUE) NA 값 있을 시 제외하고 통계 계산
#(withinSD = 1) 군내 표준편차 추정 방식 (1: R_bar, 2: S_bar, 3: 합동표준편차 방식)
#(betSd = 1) 군간 표준편차 추정 방식 (1: 이동범위의 평균 방식, 2: 이동범위의 중위수 방식)
#(maintTitle) 제목
#(subTitle) 소제목

pca.withBet <- function(data, LSL = NA,USL = NA,Target = NA, whatGroup =1, group1 = 1, group2 = c(),
                        f.na.rm = TRUE, withinSd =1, betSd = 1,
                        mainTitle = "PCA within/Between dist. graph",subTitle = "My PCA graph") {

  library(qcc)
  if(is.na(Target)){
    stop("Target is needed.")
  }
  if(is.na(LSL) & is.na(USL)){
    stop("No specification limits provided.")
  }
  if(!(withinSd == 1 ) & !(withinSd ==2) & !(withinSd ==3)){
    stop("Select 1(R_bar) or 2(S_bar) or 3(Pooled Standard deviation)")
  }
  if(!(betSd == 1 ) & !(betSd ==2)){
    stop("Select 1(Simple moving average) or 2(Simple median average)")
  }

  #부분군 갯수 정의
  if(whatGroup ==1){
    group <- group1
  }else{
    group <- nlevels(as.factor(group2))
  }
  #변수 정의 (meanST = 공정평균, N = 전체 데이터 갯수=N, n = 군내 표본 갯수(n))
  meanST <- mean(data, na.rm = f.na.rm)
  N <- length(data[!is.na(data)])
  n <- N / group
  data <- as.vector(data)

  #데이터 군별로 쪼개기
  split <- splitby(data = data, group = group)
  rowSplit <- nrow(split)
  colSplit <- ncol(split)
  #측정 데이터 평균, 이동범위, 범위, 표준편차
  split$mean <- apply(split, 1,mean)
  split$MR <- c(NA,abs(split[2:rowSplit, colSplit+1] - split[1:(rowSplit-1), colSplit+1]))
  split$range <- apply(split[1:colSplit],1,rangeMinMax)
  split$sd <- apply(split[1:colSplit] ,1,sd)
  #각 값들의 평균 구하기
  sr <- list(mean = mean(split[,colSplit+1]), MR = sum(split[,colSplit+2],na.rm = TRUE)/(group-1),
             range = mean(split[,colSplit+3]),sd = mean(split[,colSplit+4]))

  #군내 표준편차 추정(sdW)
  sdW <- NA
  if(withinSd ==1){#R 방식
    sdW <- sd.R(data = split[1:colSplit])
  }else if(withinSd ==2){#S 방식
    sdW <- sd.S(data = split[1:colSplit])
  }else{#합동표준편차방식
    sp <- getSp(split = split, rowSplit= rowSplit, colSplit = colSplit, N = N, group = group)
    tempc4 <- group*(n-1)+1
    sdW <- sp/uc[tempc4,1] # uc>c4 사용
  }

  #군간 표준편차 추정(sdX)
  sdX <- NA
  if(betSd ==1){#이동범위의 평균 방식
    sdX <- sr$MR/uc[2,3] # uc>d2
  }else{#이동범위의 중위수 방식
    sdX <- median(split$MR,na.rm = TRUE)/uc[2,5] # uc>d4
  }

  sdB <- sqrt((sdX^2)-(sdW^2/n)) # sdB = 군간 표준편차

  sdTotal <- sqrt(sdW^2 + sdB^2) #sdTotal = TOTAL S.D.

  sdOver <- sdOverAll(split,rowSplit = rowSplit, colSplit = colSplit,
                      N = N, mean =sr$mean, na.rm = f.na.rm, c4 = uc[N,1]) #sdOver = OVERALL S.D., uc>c4 사용
  cplist <- ppclist(data = data,LSL = LSL, USL = USL, meanST = meanST, sdTotal = sdTotal) #cp, cpu, cpl, cpk
  pplist <- tpplist(data = data,LSL = LSL, USL = USL, meanST = meanST, sdOver = sdOver) #pp, ppu, ppl, ppk
  obPPMlist <- obPerform(data = data, LSL = LSL, USL = USL, N = N) # 관측 성능
  overPPMlist <- expPerform_overall(data = data, LSL = LSL, USL = USL, mean = meanST, sdOver = sdOver) # 전체 기대 성능
  inPPMlist <- expPerform_ingroup(data = data, LSL = LSL, USL = USL, mean = meanST,sdIngroup = sdTotal) # 군내 기대 성능

  #결과창
  pca.prepCanvas(main = mainTitle, sub = subTitle)

  # 그래프 그리기
  vp.basic <- viewport(name = "basic", layout = grid.layout(2, 2,widths = c(0.7,0.3),heights =  c(0.7, 0.3)))
  pushViewport(vp.basic)

  vp.hist <- viewport(name = "hist", layout.pos.row = 1, layout.pos.col = 1)
  pushViewport(vp.hist)

  grid.text("Histogram & Density", y = 1, just = c("center", "top"))

  binwST <- round(diff(range(data))/sqrt(N),digits = 3)
  ggdata <- melt(data)
  qqp <- ggplot(ggdata, aes(x = data))+
    xlab("")+
    ylab("density")
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
    stat_function(fun = dnorm,args = with(ggdata,c(mean(value), sd(value))),linetype = 1,size = 1) +
    stat_function(fun = dnorm,args = with(meandata,c(mean(value), sd(value))),linetype = 2,size = 1)

  print(hist, newpage = FALSE)
  popViewport()
  # 그래프 그리기 END

  #오른쪽 표(1,2) 그리기
  vp.basicValue <- viewport(name = "basic value", layout.pos.row = 1, layout.pos.col = 2,layout = grid.layout(3, 2,heights = c(0.3,0.2,0.5)))
  pushViewport(vp.basicValue)

  grid.rect(gp = gpar(col = "gray0", lwd = 2)) # 사각형 틀 그리기

  #오른쪽 표(1,2)의 (1,1:2)부분
  vp.basicTitle <- viewport(name = "basic title", layout.pos.row = 1, layout.pos.col = 1:2)
  pushViewport(vp.basicTitle)

  grid.text(expression(bold("Density Lines Legend \n &Summary statistics")), y = 0.70, just = c("center", "top"))
  grid.lines(x = c(0.05, 0.3), y = c(0.30, 0.30), gp = gpar(lty = 1,lwd = 3))# Density ST
  grid.text("Overall", x = 0.35, y = 0.30, just = c("left",  "center"), gp = gpar(cex = 0.8))
  grid.lines(x = c(0.05, 0.3), y = c(0.10, 0.10), gp = gpar(lty = 2,  lwd = 3))#Theoretical ST
  grid.text("Within / Between", x = 0.35, y = 0.10, just = c("left", "center"), gp = gpar(cex = 0.8))

  popViewport() #오른쪽 표(1,2)의 (1,1:2)부분 END

  #오른쪽 표(1,2)의 (2,1:2)부분
  vp.input <- viewport(name = "input", layout.pos.row = 2, layout.pos.col = 1:2)
  pushViewport(vp.input)

  grid.text(expression(bold("LSL: ")), x = 0.30, y = 0.80, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(LSL,x = 0.30, y = 0.80,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Target: ")), x = 0.30, y = 0.50, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(Target, x = 0.30, y = 0.50, just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("USL: ")), x = 0.70, y = 0.80, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(USL, x = 0.70, y = 0.80, just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold("Group: ")), x = 0.70, y = 0.50, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(group, x = 0.70, y = 0.50, just = c("left", "top"), gp = gpar(cex = 0.8))


  popViewport() #오른쪽 표(1,2)의 (2,1:2)부분 END

  #오른쪽 표(1,2)의 (3,1:2)부분
  vp.sumstat <- viewport(name = "sd", layout.pos.row = 3, layout.pos.col = 1:2,layout = grid.layout(2, 2,heights = c(0.4,0.6)))
  pushViewport(vp.sumstat)

  #오른쪽 표(1,2)의 (3,1:2)부분 중 (1,1:2)
  vp.sampMean <- viewport(name = "sampMean",layout.pos.row = 1, layout.pos.col = 1:2)
  pushViewport(vp.sampMean)

  grid.lines(x = c(0, 1), y = c(1, 1), gp = gpar(lty = 3,  lwd = 3,col = "gray48"))

  grid.text(expression(bold("Sample N: ")), x = 0.50, y = 0.80, # n(데이터)
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(N, x = 0.5, y = 0.85,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold("Mean: ")), x = 0.5, y = 0.55, #평균
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(sprintf("%.4f",round(meanST,digits = 3) ), x = 0.50, y = 0.55,
            just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#오른쪽 표(1,2)의 (3,1:2)부분 중 (1,1:2) END

  #오른쪽 표(1,2)의 (3,1:2)부분 중 (2,1)
  vp.allsd <- viewport(name = "allsd", layout.pos.row = 2, layout.pos.col = 1)
  pushViewport(vp.allsd)

  grid.text(expression(bold("StDev\n(Within): ")),x = 0.6, y = 0.8, # 군내 표준편차(within s.d.)
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(sprintf("%.3f", sdW), x = 0.6, y = 0.8,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold("StDev\n(Between): ")), x = 0.65, y = 0.30, # 군간 표준편차(Between s.d.)
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(sprintf("%.3f", sdB), x = 0.65, y = 0.30,
            just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#오른쪽 표(1,2)의 (3,2)부분 즁 (2,1) END

  #오른쪽 표(1,2)의 (3,1:2)부분 중 (2,2)
  vp.allsd2 <- viewport(name = "allsd2", layout.pos.row = 2, layout.pos.col = 2)
  pushViewport(vp.allsd2)

  grid.text(expression(bold("StDev\n(Total): ")),x = 0.6, y = 0.8, # 총 표준편차(Total s.d.)
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(sprintf("%.3f", sdTotal), x = 0.6, y = 0.8,
            just = c("left", "top"), gp = gpar(cex = 0.8))
  grid.text(expression(bold("StDev\n(Overall): ")), x = 0.60, y = 0.30, # 전체 표준편차(Overall s.d.)
            just = c("right", "top"), gp = gpar(cex = 0.8))
  grid.text(sprintf("%.3f", sdOver), x = 0.60, y = 0.30,
            just = c("left", "top"), gp = gpar(cex = 0.8))

  popViewport()#오른쪽 표(1,2)의 (3,2)부분 즁 (2,2) END

  popViewport()#오른쪽 표(1,2)의 (3,1:2)부분 END

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
