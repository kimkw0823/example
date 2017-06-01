#'
#'This is used for getting total standard deviation.
#'@param c(split, rowSplit, colSplit, N, mean, na.rm = FALSE, c4 = 1)
#'@return total standard devtiation for n!=1
#'@export
sdOverAll <- function(split, rowSplit, colSplit, N, mean, na.rm = FALSE, c4 = 1){
  subpow <- 0
  for(i in 1:rowSplit){
    for(j in 1:colSplit){
      subpow <- subpow + (split[i,j]-mean)^2
    }
  }
  overallBefore <- sqrt(subpow/(N-1))
  result<- overallBefore/c4
  return(result)
}

#`
#'This is used for getting total standard deviation.(when n ==1)
#'@param c(split, rowSplit, colSplit, N, mean, na.rm = FALSE, c4 = 1)
#'@return total standard devtiation for n==1
#'@export
sdOverAll2 <- function(split, rowSplit, colSplit, N, mean, na.rm = FALSE, c4 = 1){
  subpow <- 0
  for(i in 1:rowSplit){
    for(j in 1:colSplit){
      subpow <- subpow + (split[i,j]-mean)^2
    }
  }
  overallBefore <- sqrt(subpow/(N-1))
  result<- overallBefore/c4
  return(result)
}
#'
#'This is used for spliting data by group
#'@param c(data,group)
#'@return splited data
#'@export
splitby <- function (data, group){

  reshape <- function(df, nrow, ncol, byrow = TRUE)
    data.frame(matrix(as.matrix(df), nrow, ncol, byrow = byrow))
  ndata <- length(data)
  ncolumn <- floor(ndata/group)
  result <- reshape(df = data, nrow = group, ncol = ncolumn)
  return(result)
}

#'
#'This is used for subtract max to min
#'@param data
#'@return max - min
#'@export
rangeMinMax <- function (data){

  result <- range(data)[2]-range(data)[1]
  return(result)
}

#'
#'This is used to get Sp
#'@param data
#'@return sp
#'@export
getSp <- function (split, rowSplit, colSplit, N, group){
  n <- N / group
  subpow <- 0
  for(i in 1:rowSplit){
    for(j in 1:colSplit){
      subpow <- subpow + (split[i,j] - split[i,colSplit+1])^2
    }
  }
  sp <- subpow/(group*(n-1))
  return(sp)
}

#'
#'This is used to get prepare visualization
#'@param data
#'@return prep Canvas
#'@export
pca.prepCanvas<-function(main="Process Capability Analysis Graph", sub="My PCA project",
                         pca.col=c("#666666", "#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")){

      grid.newpage()
      grid.rect(gp=gpar(col=pca.col[2], lwd=2, fill=pca.col[5]))
      #vp.canvas
      vp.canvas<-viewport(name="canvas", width=unit(1,"npc")-unit(6,"mm"),
                          height=unit(1,"npc")-unit(6,"mm"),
                          layout=grid.layout(3,1,heights=unit(c(3,1,2), c("lines", "null", "lines"))))
      pushViewport(vp.canvas)
      grid.rect(gp=gpar(col="#FFFFFF", lwd=0, fill="#FFFFFF"))

        #Title
        vp.title<-viewport(layout.pos.col=1, layout.pos.row=1, name="title")
        pushViewport(vp.title)
        grid.text (main, gp=gpar(fontsize=20))
        popViewport()

        #Subtitle
        vp.subtitle<-viewport(layout.pos.col=1, layout.pos.row=3, name="subtitle")
        pushViewport(vp.subtitle)
        grid.text (sub, gp=gpar(col=pca.col[1]))
        popViewport()

        #Container
        vp.container<-viewport(layout.pos.col=1, layout.pos.row=2, name="container")
        pushViewport(vp.container)

      }


#`
#'getting GCD
#'@param data
#'@return GCD
#'@export
gcd <- function(x,y){
  r <- x%%y;
  return(ifelse(r, gcd(y, r), y))
}
#`
#'getting GCD
#'@param data
#'@return GCD
#'@export
pca.gcd <- function(data){

  l <- length(data)
  for(i in 1:l-1){

    a <- gcd(data[i],data[i+1])

  }
  return(a)
}
