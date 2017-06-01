#'
#'@param
#'@return pca.bino
#'@export

#data data$variable 불량
#(Target = NA) 타겟
#(whatGroup = 1) 표본 크기를 무엇으로 구할건지(1: 상수, 2: 데이터 셋 열)
#(N1 = 1) 표본 갯수                (whatGroup이 1일 때 사용)
#(N2 = 1) 사용될 크기              (whatGroup이 2일 때 사용)
#(alpha = 0.05)
#(f.na.rm = TRUE) NA 값 있을 시 제외하고 통계 계산
#(maintTitle) 제목
#(subTitle) 소제목

pca.bino <- function(data, Target = NA, whatGroup =1, N1 = 1, N2 = c(),
                     alpha = 0.05,f.na.rm = TRUE,
                     mainTitle = "Binomial distribution", subTitle = "My PCA project") {
  library(grid)
  library(ggplot2)
  library(reshape2)
  N <- NA
  D <- NA
  percentData <- c()

  #총 검사갯수
  if(whatGroup ==1){
    N <- N1 * length(data)
    percentData <- (data / N1) *100
  }else if(whatGroup==2){
    N <- sum(N2,na.rm = f.na.rm)
    percentData <- (data/N2) *100
  }else{
    stop("Use appropriate value.")
  }

  #총 부적합품 수
  D <- sum(data,na.rm = f.na.rm)
  #공정 부적합품률(불량률) p의 추정량 p_hat
  p_hat <- D/N

  #p의 100*(1-alpha)% 신뢰구간(신뢰하한 : pl, 신뢰상한 : pu)
  pl <- NA
  pu <- NA
  pi1 <- NA ; pi2 <- NA ; pi3<-NA ; pi4 <- NA
  pi1 <- 2*(N-D+1)
  pi2 <- 2*D
  pi3 <- 2*(D+1)
  pi4 <- 2*(N-D)
  pl <- pi2 / (pi1*qf(alpha/2,df1 = pi1,df2 = pi2,lower.tail = FALSE)+pi2)
  pu <- (pi3*qf(alpha/2,df1=pi3,df2 = pi4,lower.tail = FALSE))/(pi3*qf(alpha/2,df1=pi3,df2=pi4,lower.tail = FALSE)+pi4)


  # %불량품
  percentD <- p_hat *100
  # %불량품의 100*(1-alpha)% 신뢰구간
  confInter <- list(CI_pl = round(pl*100,digits=3), CI_pu = round(pu *100,digits=3))

  # 불량품 PPM
  ppm <- p_hat * 1000000
  # 불량품 PPM의 100*(1-alpha)% 신뢰구간
  ppm_confInter <- list(CI_pl = round(pl*1000000,digits = 3), CI_pu = round(pu * 1000000,digits=3))

  # 공정 Z
  Z <- round(qnorm(1-p_hat),digits=3)
  # 공정 Z의 100*(1-alpha)% 신뢰구간
  Z_confInter <- list(CI_pl = round(qnorm(1-pl),digits=3),CI_pu = round(qnorm(1-pu),digits=3))

  #결과창
  pca.prepCanvas(main = mainTitle, sub = subTitle)

  vp.visualize <- viewport(name = "visualize", layout = grid.layout(1, 2,widths = c(0.7,0.3)))
  pushViewport(vp.visualize)


  # 그래프 그리기(1:1)
  vp.hist <- viewport(name = "hist", layout.pos.row = 1, layout.pos.col = 1)
  pushViewport(vp.hist)

  binwST <- round(diff(range(percentData))/sqrt(N),digits = 3)
  ggdata <- melt(percentData)

  qqp <- ggplot(ggdata, aes(x = percentData)) +
    scale_x_continuous(breaks=seq(from = range(percentData)[1]-pca.gcd(percentData)*2,
                                  to = range(percentData)[2]+pca.gcd(percentData)*2,
                                  by = pca.gcd(percentData)*2))+
    xlab("%inferior goods")+
    ylab("frequency")
  hist <- qqp + geom_histogram(aes(y = ..count..), binwidth = pca.gcd(percentData), fill = "royalblue3",
                               colour = "black",stat = "bin")
  x_density <- density(percentData, bw = binwST)
  my.ggp.yrange <- ggplot_build(hist)$layout$panel_ranges[[1]]$y.range
  if(is.na(Target)){
    nTarget <- 0.00
  }else{
    nTarget <- Target
  }
  if(!is.na(Target)){
    hist <- hist + annotate(geom = "text", x = Target, y = my.ggp.yrange[2]-abs(my.ggp.yrange[1]),
                            label = "Target", hjust = -0.1, size = 5)
    hist <- hist + geom_vline(xintercept = Target,linetype = 3,size = 1)
  }
  print(hist, newpage = FALSE)

  popViewport() #그래프 그리기 (1:1)

  #summary statistics (1,2)
  vp.ss <- viewport(name = "ss", layout.pos.row = 1, layout.pos.col = 2)
  pushViewport(vp.ss)

  grid.rect(gp = gpar(col = "gray0", lwd = 2)) # 사각형 틀 그리기
  grid.text(expression(bold("Summary statistics")), y = 0.95, just = c("center", "top"),gp=gpar(cex = 1.0)) #제목
  grid.text("(95% confidence interval)",y = 0.88, just = c("center",  "top"), gp = gpar(cex = 0.7))


  grid.text(expression(bold("%inferior goods: ")), x = 0.55, y = 0.75, just = c("right", "top"),gp = gpar(cex = 0.8)) # 내용
  grid.text(percentD,x = 0.73, y = 0.75,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Lower CI: ")), x = 0.55, y = 0.68, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(confInter$CI_pl,x = 0.70, y = 0.68,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Upper CI: ")), x = 0.55, y = 0.61, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(confInter$CI_pu,x = 0.73, y = 0.61,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Target: ")), x = 0.55, y = 0.54, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(nTarget,x = 0.73, y = 0.54,just = c("left", "top"),gp = gpar(cex = 0.8))

  grid.text(expression(bold("inferior goods PPM: ")), x = 0.55, y = 0.48, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(ppm,x = 0.73, y = 0.48,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Lower CI: ")), x = 0.55, y = 0.42, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(ppm_confInter$CI_pl,x = 0.73, y = 0.42,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Upper CI: ")), x = 0.55, y = 0.36, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(ppm_confInter$CI_pu,x = 0.73, y = 0.36,just = c("left", "top"),gp = gpar(cex = 0.8))

  grid.text(expression(bold("Z.bench: ")), x = 0.55, y = 0.30, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(Z,x = 0.73, y = 0.30,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Lower CI: ")), x = 0.55, y = 0.24, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(Z_confInter$CI_pl,x = 0.73, y = 0.24,just = c("left", "top"),gp = gpar(cex = 0.8))
  grid.text(expression(bold("Upper CI: ")), x = 0.55, y = 0.18, just = c("right", "top"),gp = gpar(cex = 0.8))
  grid.text(Z_confInter$CI_pu,x = 0.73, y = 0.18,just = c("left", "top"),gp = gpar(cex = 0.8))

  grid.text(expression(bold("* CI : Confidence Interval ")), x = 0.55, y = 0.05, just = c("left", "top"),gp = gpar(cex = 0.6))

  popViewport()#summary statistics (1,2) END
}
