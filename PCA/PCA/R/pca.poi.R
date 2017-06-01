#'
#'@param
#'@return pca.poi
#'@export

#data data$variable 불량
#(whatGroup = 1) 표본 크기를 무엇으로 구할건지(1: 상수, 2: 데이터 셋 열)
#(N1 = 1) 표본 갯수                (whatGroup이 1일 때 사용)
#(N2 = 1) 사용될 크기              (whatGroup이 2일 때 사용)
#(Target = NA) 타겟
#(alpha = 0.05)
#(f.na.rm = TRUE) NA 값 있을 시 제외하고 통계 계산
#(maintTitle) 제목
#(subTitle) 소제목

pca.poi <- function(data, whatGroup =1, N1 = 1, N2 = c(),Target = NA,
                    alpha = 0.05, f.na.rm = TRUE,
                    mainTitle = "Process Capability Analysis Graph", subTitle = "My PCA project"){

  k <- NA # 부분군 수
  N <- NA # 총 검사 단위
  D <- NA # 총 부적합수(결점수)

  #k : 부분군 수
  k <- length(data)
  #N : 총 검사갯수
  if(whatGroup ==1){
    N <- N1 * length(data)
  }else if(whatGroup==2){
    N <- sum(N2,na.rm =  f.na.rm)
  }else{
    stop("Use appropriate value.")
  }
  #D : 총 부적합품 수
  D <- sum(data,na.rm =  f.na.rm)

  #부분군당 평균 부적합수(결점수) = 평균 결점
  meanD <- D/k

  #자유도 pi1, pi2
  pi1 <- NA ; pi2 <- NA
  pi1 <- 2*D
  pi2 <- 2*(D+1)

  #부분군당 평균 부적합수(결점수)의 100*(1-alpha)% 신뢰구간 (신뢰하한 : pl, 신뢰상한 : pu)
  perk_confInter <- list(CI_pl = round(0.5*(1/k)*qchisq(alpha/2,pi1),digits = 3),
                         CI_pu = round(0.5*(1/k)*qchisq(alpha/2,pi2,lower.tail = FALSE),digits = 3))



  #단위당 평균 부적합수(결점수) = DPU
  DPU <- D/N

  #단위당 평균 부적합수(결점수)의 100(1-alpha)% 신뢰구간 (신뢰하한 : pl, 신뢰상한 : pu)
  perN_confInter <- list(CI_pl = round(0.5*(1/N)*qchisq(alpha/2,pi1),digits =3),
                         CI_pu = round(0.5*(1/N)*qchisq(alpha/2,pi2,lower.tail = FALSE),digits =3))


  #결과창
  pca.prepCanvas(main = mainTitle, sub = subTitle)

  vp.visualize <- viewport(name = "visualize", layout = grid.layout(1, 2,widths = c(0.7,0.3)))
  pushViewport(vp.visualize)


  # 그래프 그리기(1:1)
  vp.hist <- viewport(name = "hist", layout.pos.row = 1, layout.pos.col = 1)
  pushViewport(vp.hist)

  data1 <-data/100
  binwST <- round(diff(range(data1))/sqrt(D/100),digits = 3)
  ggdata <- melt(data1)

  qqp <- ggplot(ggdata, aes(x = ggdata)) +
    scale_x_continuous(breaks=seq(from = range(data1)[1]-pca.gcd(data1)*2,
                                  to = range(data1)[2]+pca.gcd(data1)*2,
                                  by = pca.gcd(data1)*2))+
    xlab("DPU")+
    ylab("frequency")
  hist <- qqp + geom_histogram(aes(y = ..count..), binwidth = pca.gcd(data1), fill = "royalblue3",
                               colour = "black",stat = "bin")

  x_density <- density(data1, bw = binwST)
  my.ggp.yrange <- ggplot_build(hist)$layout$panel_ranges[[1]]$y.range
  if(is.na(Target)){
    nTarget <- 0
  }else{
    nTarget <- Target/100
  }
  if(!is.na(Target)){
    hist <- hist + annotate(geom = "text", x = nTarget, y = my.ggp.yrange[2]-abs(my.ggp.yrange[1]),
                            label = "Target", hjust = -0.1, size = 5)
    hist <- hist + geom_vline(xintercept = nTarget,linetype = 3,size = 1)
  }
  print(hist, newpage = FALSE)

  popViewport() #그래프 그리기 (1:1)

  #summary statistics (1,2)
  vp.ss <- viewport(name = "ss", layout.pos.row = 1, layout.pos.col = 2)
  pushViewport(vp.ss)

  grid.rect(gp = gpar(col = "gray0", lwd = 2)) # 사각형 틀 그리기
  grid.text(expression(bold("Summary statistics")), y = 0.95, just = c("center", "top"),gp=gpar(cex = 1.0)) #제목
  grid.text("(95% confidence interval)",y = 0.88, just = c("center",  "top"), gp = gpar(cex = 0.7))


  grid.text(expression(bold("Average defects: ")), x = 0.55, y = 0.75, just = c("right", "top"),gp = gpar(cex = 0.8)) # 내용
  grid.text(meanD,x = 0.73, y = 0.75,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Lower CI: ")), x = 0.55, y = 0.68, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(perN_confInter$CI_pl,x = 0.70, y = 0.68,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Upper CI: ")), x = 0.55, y = 0.61, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(perN_confInter$CI_pu,x = 0.73, y = 0.61,just = c("left", "top"),gp = gpar(cex = 0.8))

  grid.text(expression(bold("mean DPU: ")), x = 0.55, y = 0.54, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(DPU,x = 0.73, y = 0.54,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Lower CI: ")), x = 0.55, y = 0.48, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(perN_confInter$CI_pl,x = 0.73, y = 0.48,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Upper CI: ")), x = 0.55, y = 0.42, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(perN_confInter$CI_pu,x = 0.73, y = 0.42,just = c("left", "top"),gp = gpar(cex = 0.8))

  grid.text(expression(bold("Min DPU: ")), x = 0.55, y = 0.36, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(range(data1)[1],x = 0.73, y = 0.36,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Max DPU: ")), x = 0.55, y = 0.30, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(range(data1)[2],x = 0.73, y = 0.30,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Target DPU: ")), x = 0.55, y = 0.24, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(nTarget,x = 0.73, y = 0.24,just = c("left", "top"),gp = gpar(cex = 0.8))

  grid.text(expression(bold("* CI : Confidence Interval ")), x = 0.55, y = 0.05, just = c("left", "top"),gp = gpar(cex = 0.6))

  popViewport()#summary statistics (1,2) END
}
