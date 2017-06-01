
#'
#'This is used for observation performance
#'
#'@return list of observation performance
#'@export

obPerform <- function(data,LSL, USL, N){

  ppm_down <- round((length(data[data<LSL])/N)*1000000,digits=3)
  ppm_up <- round((length(data[data>USL])/N)*1000000,digits=3)
  ppm_total <- round(ppm_down + ppm_up,digits=3)
  obPerform_list <- list(ppm_down=ppm_down, ppm_up = ppm_up, ppm_total = ppm_total)
  return(obPerform_list)
}

#'
#'This is used for overall expected performance
#'
#'@return list of overall expected performance
#'@export

expPerform_overall <- function(data, LSL, USL,mean,sdOver){

  ppm_down <- round(pnorm((LSL-mean)/sdOver) * 1000000,digits =3)
  ppm_up <- round(pnorm((USL-mean)/sdOver,lower.tail = FALSE) *1000000,digits=3)
  ppm_total <- round(ppm_down + ppm_up,digits=3)
  overPerform_list <- list(ppm_down = ppm_down, ppm_up = ppm_up, ppm_total = ppm_total)
  return(overPerform_list)

}


#'
#'This is used for ingroup expected performance
#'
#'@return ingroup expected performance
#'@export

expPerform_ingroup <- function(data, LSL, USL, mean,sdIngroup){

  ppm_down <- round(pnorm((LSL-mean)/sdIngroup) * 1000000,digits=3)
  ppm_up <- round(pnorm((USL-mean)/sdIngroup,lower.tail =FALSE) * 1000000,digits=3)
  ppm_total <- round(ppm_down + ppm_up,digits=3)
  inPerform_list <- list(ppm_down = ppm_down, ppm_up = ppm_up, ppm_total = ppm_total)
  return(inPerform_list)

}
