#' Resampling of growth increment data (dL/dt) from mark-recapture (tagging) studies
#'
#' @param 'increm' a growth increment object of the class `data.frame`.
#' @param nresamp numeric; the number of permutations to run (Default: `nresamp = 200`)

#' @description The function resamples the `increm` data by rows 
#'   (i.e., by recapture  date) several times ('nresamp' times, default: nresamp = 200).
#'   Then, a VBGF curve is fitted to each resapled data set. The  resulting 
#'   posteriors (e.g., VGBGF function parameters K and Linf) are stored in a data.frame
#'   called '$bootRaw'. 
#'   The '$bootRaw' table also includes the growth performance index Phi', 
#'   seaonal parameters "u" and "w" sensu Francis,1988, whih are equal to "C" 
#'   sensu Pauly and Gaschütz, 1979 and and "ts" sensu Mildenberger et al., 2017.
#'   The '$bootRaw' table also includes seed values and system time. 
#'
#'  @param L1 Vector of length at release of tagged fish (see fishmethods::grotag)
#'  @param L2 Vector of length at recovery of tagged fish (see fishmethods::grotag)
#'  @param T1 Vector of julian time at release of tagged fish (see fishmethods::grotag)
#'  @param T2 Vector of julian time at recovery of tagged fish (see fishmethods::grotag)
#'  @param alpha Numeric value giving an arbitrary length alpha (see fishmethods::grotag)
#'  @param beta Numeric value giving an arbitrary length beta (beta> alpha) (see fishmethods::grotag)
#'  @param design List specifying the design of the model to estimate. Use 1 to designate 
#'  ether a parameter(s) should be estimated.  Type of parameters are:  nu=growth 
#'  variability (1 parameter),  m=bias parameter of measurement error (1 parameter), 
#'  p=outlier probability (1 parameter), and sea=seasonal variation (2 parameters: u and w). (see fishmethods::grotag)
#'  Model 1 of Francis is the default settings of 0 for nu, m, p and sea. (see fishmethods::grotag)
#'  @param stvalue Starting values of sigma (s) and depending on the design argument, nu, m, p, u, and w 
#'  used as input in the nonlinear estimation (function optim) routine.
#'  @param upper Upper limit of the model parameters' (nu, m, p, u, and w) region to be 
#'  investigated.(see fishmethods::grotag)
#'  @param lower Lower limit of the model parameters' (nu, m, p, u, and w) region to be 
#'  investigated. (see fishmethods::grotag)
#'  @param gestimate Logical specifying whether starting values of ga and gb 
#'  (growth increments of alpha and beta) should be estimated automatically. 
#'  Default = TRUE.
#'  @param st.ga If gestimate=FALSE, user-specified starting value for ga. (see fishmethods::grotag)
#'  @param st.gb If gestimate=FALSE, user-specified starting value for gb. (see fishmethods::grotag)
#'  @param st.galow If gestimate=FALSE, user-specified lower limit for st.ga used in optimization. (see fishmethods::grotag)
#'  @param st.gaup If gestimate=FALSE, user-specified upper limit for st.ga used in optimization. (see fishmethods::grotag)
#'  @param st.gblow If gestimate=FALSE, user-specified lower limit for st.gb used in optimization. (see fishmethods::grotag)
#'  @param st.gbup If gestimate=FALSE, user-specified upper limit for st.gb used in optimization. (see fishmethods::grotag)
#'  @param control Additional controls passed to the optimization function optim. (see fishmethods::grotag)
#' 
#' 
#'@return a list of class `lfqBoot` containing 2 levels: 
#'  `$bootRaw` - a data.frame of fitted VBGF parameters (columns) by resampling (rows), 
#'  `$seed` - a vector of seed values set prior to each resampling call to `Increm_Resample`.
#'
#'
#' @examples
#'
# library(fishmethods)
# 
# data(bonito)
# 
# res <- grotag_boot(bonito, nresamp = 10,
#                         L1=L1, L2=L2, T1=T1, T2=T2,alpha=35,beta=55,
#                         design=list(nu=1,m=1,p=1,sea=1),
#                         stvalue=list(sigma=0.9,nu=0.4,m=-1,p=0,u=0.4,w=0.4),
#                         upper=list(sigma=5,nu=1,m=2,p=0.5,u=1,w=1),
#                         lower=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0),control=list(maxit=1e4))
# 
# res$bootRaw
# res$seed
# 
# library(fishboot)
# LinfK_scatterhist(res)
# 


grotag_boot <- function(input.data, nresamp = 200,
                             L1 = NULL, L2 = NULL, T1 = NULL, T2 = NULL, alpha = NULL, beta = NULL,
                             design = list(nu = 0, m = 0, p = 0, sea = 0),
                             stvalue = list(sigma = 0.9, nu = 0.4, m = -1, p = 0.01, u = 0.4, w = 0.4),
                             upper = list(sigma = 5, nu = 1, m = 2, p = 1, u = 1, w = 1),
                             lower = list(sigma = 0, nu = 0, m = -2, p = 0, u = 0, w = 0), gestimate = TRUE,
                             st.ga = NULL, st.gb = NULL, st.galow = NULL, st.gaup = NULL, st.gblow = NULL,
                             st.gbup = NULL, control = list(maxit = 10000))
  
{ 

 
res_tab <- data.frame(rep.int(NA,(nresamp)),
                      rep.int(NA,(nresamp)),
                      rep.int(NA,(nresamp)),
                      rep.int(NA,(nresamp)),
                      rep.int(NA,(nresamp)),
                      rep.int(NA,(nresamp)),
                      rep.int(NA,(nresamp))) # creates an empty results table for
                      

names(res_tab) <- c( "Linf", "K", "PhiL", "u", "w", "seed","time")


for(x in 1:(nresamp))
  
{
  
  tryCatch({
  
  seed <- round((runif(1, min = 0, max = 1000)), 0)      
  seed_final <- (x + seed +.Random.seed[x]) 
  set.seed <- seed_final
  
  sam_size <-  nrow(input.data)
  input.rownames <- row.names(input.data)
  samplerows <- sample(input.rownames, sam_size, replace = TRUE)
  input.data_p <- input.data[samplerows,]
  
  
  
  fitboot <- with(input.data_p,
                  fishmethods::grotag( 
                    L1=L1, L2=L2, T1=T1, T2=T2,alpha=35,beta=55,
                    design=list(nu=1,m=1,p=1,sea=1),
                    stvalue=list(sigma=0.9,nu=0.4,m=-1,p=0,u=0.4,w=0.4),
                    upper=list(sigma=5,nu=1,m=2,p=0.5,u=1,w=1),
                    lower=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0),control=list(maxit=1e4)))
  
  fitboot
  
  Linf <-  as.numeric (- fitboot$VBparms$Estimate[2]) # Linf
  K  <-  as.numeric  (fitboot$VBparms$Estimate[3])    # K
  PhiL <- log10(K) + 2 * log10(Linf)                  # Phi prime
  
  u_estimate <- fitboot$table$Estimate[4]  # seasonal amplitude (sensu Pauly and Gaschütz 1979, Francis, 1988)
  w_estimate  <-  fitboot$table$Estimate[5]  # summer point  (sensu Francis, 1988)
  
  
  res_tab[x,1] <- Linf
  res_tab[x,2] <- K
  res_tab[x,3] <- PhiL
  res_tab[x,4] <- u_estimate
  res_tab[x,5] <- w_estimate
  res_tab[x,6] <- seed_final
  res_tab[x,7] <- as.character(paste(Sys.time()))
  
  
  
  # remove NAs
  res_tab <- res_tab[complete.cases(res_tab), ]
  
  }, error=function(e){})
  
  }

ret <- list()
ret$bootRaw <- res_tab
ret$seed <- seed
class(ret) <- "lfqBoot"

return(ret)


}

