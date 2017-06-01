#'
#'@param
#'@return pca.nonNorm
#'@export

#data data$variable 불량
#(LSL = NA) 하한
#(USL = NA) 상한
#(Target = NA) 타겟
#(maintTitle) 제목

pca.noNorm <- function(data,LSL = NA,USL = NA, target = NA,
                       mainTitle = "My PCA project"){
  library(qualityTools)

  dist <- c("normal","exponential","weibull","gamma","logistic")
  admax <- 0
  addist <- ""
  cat("Searching best distribution for this data... It may take a few seconds...")
  for(i in 1:5){
    ad <- pca.sim(data,distribution = dist[i])
    adcomp <- ad$p_value
    admax <- max(admax,adcomp)
    if(adcomp ==admax){
      addist <- dist[i]
    }
  }
  if(is.na(LSL)){
    LSL <- NULL
  }
  if(is.na(USL)){
    USL <- NULL
  }
  if(is.na(target)){
    target <- NULL
  }
  pcr(x=data, distribution = addist,lsl = LSL,usl = USL,target = target, main = mainTitle,
      lineCol="gray0",specCol = "red3",lineWidth = 2)
}
