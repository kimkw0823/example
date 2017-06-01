#'
#'
#' MSA destructive Test
#'@param data
#'@return msa.destTest
#'@export

# operators           측정자
# parts               부품
# results             측정값
# data                사용데이터셋
# method              gageR&R test 방법(msa.destTest 고정값 : "nested")
# lsl                 하한
# usl                 상한
# alpha               유의수준 (default = 0.05)
# sigma               표준편차 (default = 6)
# digits              최소단위 (default = 5)
# main                대제목(default = "Gage R&R destructive Test")
# sub                 소제목 (default = "My MSA project")

msa.destTest <- function(operators, parts, results, data, method = "nested",
                         lsl, usl, alpha = 0.05, sigma = 6, digits = 5,
                         main = "Gage R&R destructive Test", sub = "My MSA project"){

  library(MSA)
  library(qualityTools)

  if (is.data.frame(data)) {
    if (deparse(substitute(results)) %in% names(data)) {
      results <- deparse(substitute(results))
    }
    if (!(results %in% names(data))) {
      stop(results, "is not a valid column name for", deparse(substitute(data)))
    }
    if (deparse(substitute(parts)) %in% names(data)) {
      parts <- deparse(substitute(parts))
    }
    if (deparse(substitute(operators)) %in% names(data)) {
      operators <- deparse(substitute(operators))
    }
    if (parts %in% names(data)) {
      data[[parts]] <- factor(data[[parts]])
    }else {
      stop(parts, "is not a valid column name for", data)
    }
    if (operators %in% names(data)) {
      data[[operators]] <- factor(data[[operators]])
    }else {
      stop(operators, "is not a valid column name for", data)
    }
  }else {
    stop("A data.frame object is needed as data argument")
  }

  op_level <- nlevels(data[[operators]])
  pa_level <- nlevels(data[[parts]]) / op_level
  me_level <- length(data[[results]]) / nlevels(data[[parts]])

  gdo <- gageRRDesign(Operators = op_level, Parts = pa_level, Measurements = me_level,
                     method = method,sigma = sigma, randomize = FALSE)

  response(gdo) <- data[[results]]


  if(!missing(lsl) && !missing(usl)){
    gdo <- gageRR(gdo, method=method, sigma=sigma, alpha=alpha,tolerance = c(lsl,usl),dig = digits)
    # msa.destGageTable(gdo, method=method, sigma=sigma, alpha = alpha,tolerance = c(lsl, usl) ,
    #               dig = digits,main = main,sub = sub)

  }else{
    gdo <- gageRR(gdo,method=method,sigma=sigma,alpha=alpha,dig = digits)
    # msa.destGageTable(gdo, method=method, sigma=sigma, alpha = alpha,dig = digits,main = main,sub = sub)

  }

  plot(gdo)
}
