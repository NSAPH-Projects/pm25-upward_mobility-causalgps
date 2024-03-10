percent_change <- function(x) {
  (exp(x) - 1) * 100
}

confint <- function(est, se, level = 0.95) {
  
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  ci = est + se %o% qnorm(a)
  return(ci)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



RR <- function(est) {
  class(est) <- c("RR", "estimate")
  est
}

#' @rdname effect_measures
#' @export
OR <- function(est, rare) {
  class(est) <- c("OR", "estimate")
  attr(est, "rare") <- rare
  est
}

#' @rdname effect_measures
#' @export
HR <- function(est, rare) {
  class(est) <- c("HR", "estimate")
  attr(est, "rare") <- rare
  est
}

#' @rdname effect_measures
#' @export
RD <- function(est) {
  class(est) <- c("RD", "estimate")
  est
}

#' @rdname effect_measures
#' @export
OLS <- function(est, sd) {
  class(est) <- c("OLS", "estimate")
  attr(est, "sd") <- sd
  est
}

#' @rdname effect_measures
#' @export
MD <- function(est) {
  class(est) <- c("MD", "estimate")
  est
}

#' @export
print.estimate <- function(x, ...) {
  attr(x, "sd") <- NULL
  attr(x, "rare") <- NULL
  attr(x, "history") <- NULL
  class(x) <- "numeric"
  print.default(x, ...)
}

#' @export
summary.estimate <- function(object, ...) {
  if (is.null(attr(object, "history"))) return(cat(class(object)[1], "=", object))
  history <- attr(object, "history")
  cat(class(object)[1], "=", object,
      "\nThis is an approximate conversion of the original", 
      history[1,1], "estimate =", history[1,2])
}

#' Convert an effect measure
#'
#' @description These helper functions are mostly used internally to convert
#'   [effect measures][effect_measures] for the calculation of E-values. The
#'   approximate conversion of odds and hazard ratios to risk ratios depends on
#'   whether the rare outcome assumption is made.
#' @name convert_measures
#' @param est The effect estimate; constructed with one of [RR()], [OR()], [HR()],
#'   [MD()], [OLS()].
#' @param rare When converting a [OR()] or [HR()] estimate, a logical indicating
#'   whether the outcome is sufficiently rare to approximate a risk ratio.
#' @param delta When converting an [OLS()] estimate, the contrast of interest 
#'   in the exposure. Defaults to 1 (a 1-unit contrast in the exposure).
#' @param ... Arguments passed to other methods.
#' @return An object of class "estimate" and the desired effect measure. Also
#'   includes as an attribute its conversion history.
#' @details Uses the conversions listed in Table 2 of VanderWeele TJ, Ding P.
#'   *Sensitivity Analysis in Observational Research: Introducing the E-Value.*
#'   Annals of Internal Medicine. 2017;167(4):268–75.
#'
#'   See references.
#'
#'   Regarding the continuous outcome, the function uses the effect-size
#'   conversions in Chinn (2000) and VanderWeele (2017) to approximately convert
#'   the mean difference between these exposure "groups" to the odds ratio that
#'   would arise from dichotomizing the continuous outcome.
#'
#' @references Chinn, S (2000). A simple method for converting an odds ratio to
#' effect size for use in meta-analysis. \emph{Statistics in Medicine}, 19(22),
#' 3127-3131.
#'
#' VanderWeele, TJ (2017). On a square-root transformation of the odds ratio for
#' a common outcome. \emph{Epidemiology}, 28(6), e58.
#'
#' VanderWeele TJ (2020). *Optimal approximate conversions of odds ratios and
#' hazard ratios to risk ratios.* Biometrics.
#' @examples
#' # Both odds ratios are 3, but will be treated differently
#' # depending on whether rare outcome assumption is reasonable
#' OR(3, rare = FALSE)
#' OR(3, rare = TRUE)
#' toRR(OR(3, rare = FALSE))
#' toRR(OR(3, rare = TRUE))
#' attributes(toRR(toMD(OLS(3, sd = 1.2), delta = 1)))
#' @rdname convert_measures
#' @export
toRR <- function(est, rare, delta = 1, ...) {
  UseMethod("toRR", est)
}
#' @rdname convert_measures
#' @export
toMD <- function(est, delta = 1, ...) {
  UseMethod("toMD", est)
}



#' @export
toMD.OLS <- function(est, delta = 1, ... ) {
  sd_attr <- attr(est, "sd")
  
  if (is.null(sd_attr)) 
    stop("Must specify the outcome standard deviation. Use argument sd = in the OLS() function")
  
  MD <- est * delta / sd_attr
  class(MD) <- c("MD", "estimate")
  attr(MD, "history") <- rbind(attr(est, "history"), c("OLS", est))
  MD
}

#' @export
toRR.MD <- function(est, ... ) {
  RR <- exp( 0.91 * est )
  class(RR) <- c("RR", "estimate")
  attr(RR, "history") <- rbind(attr(est, "history"), c("MD", est))
  RR
}

#' @export
toRR.OLS <- function(est, rare = NULL, delta = 1, ... ) {
  toRR(toMD(est, delta = delta))
}

#' @export
toRR.HR <- function(est, rare, ... ) {
  rare_attr <- attr(est, "rare")
  
  if (is.null(rare_attr)) 
    stop("Must specify whether the rare outcome assumption can be made. Use argument rare = in the HR() function.")
  
  if (rare_attr) RR <- est else {
    RR <- ( 1 - 0.5^sqrt( est ) ) / ( 1 - 0.5^sqrt( 1 / est ) )
  }
  class(RR) <- c("RR", "estimate")
  attr(RR, "history") <- rbind(attr(est, "history"), c("HR", est))
  RR
}

#' @export
toRR.OR <- function(est, rare, ... ) {
  rare_attr <- attr(est, "rare")
  
  if (is.null(rare_attr)) 
    stop("Must specify whether the rare outcome assumption can be made. Use argument rare = in the OR() function.")
  
  if (rare_attr) RR <- est else RR <- sqrt(est)
  class(RR) <- c("RR", "estimate")
  attr(RR, "history") <- rbind(attr(est, "history"), c("OR", est))
  RR
}


#' @export
toRR.default <- function(est, ... ) {
  stop("RR conversion is currently available only for estimates of class \"OR\", \"HR\", \"MD\", and \"OLS\"")
}


#' @export
toMD.default <- function(est, ...) {
  stop("MD conversion is currently available only for estimates of class \"OLS\"")
}

#########################

evalues.OLS = function( est, se = NA, sd, delta = 1, true = 0, ... ) {
  
  if ( !is.na( se ) ) {
    if ( se < 0 ) stop( "Standard error cannot be negative" )
  }
  
  if ( delta < 0 ) {
    delta = -delta
    print( "Recoding delta to be positive" )
  }
  
  if ( !inherits(est, "OLS") ) est = OLS( est, sd = sd )
  if ( !inherits(se, "OLS") ) se = OLS( se, sd = attr(est, "sd") )
  if ( !inherits(true, "MD") ) true = MD( true )
  
  # rescale to reflect a contrast of size delta
  est = toMD( OLS( est, sd = sd ), delta = delta )
  se = toMD( OLS(se, se), delta = delta )
  
  return( evalues.MD( est = est, se = se, true = true ) )
}




evalues.MD = function( est, se = NA, true = 0, ... ) {
  
  if ( !is.na( se ) ) {
    if ( se < 0 ) stop( "Standard error cannot be negative" )
  }
  
  if ( !inherits(est, "MD") ) est = MD(est)
  if ( !inherits(true, "MD") ) true = MD(true)
  
  lo = NA
  hi = NA
  if ( !is.na(se) ) {
    lo = exp( 0.91 * est - 1.78 * se )
    hi = exp( 0.91 * est -- 1.78 * se )
    #lo =  exp( log( est ) - 1.96 * log( MDtoRR( se ) )) # ( est converted )
    #hi =  exp( log( est ) + 1.96 * log( MDtoRR( se ) ))
  }
  
  if ( !is.na(lo) ) lo = RR(lo)
  if ( !is.na(hi) ) hi = RR(hi)
  est = toRR(est)
  true = toRR(true)
  
  return( evalues.RR( est = est, lo = lo, hi = hi, true = true ) )
}




evalues.HR = function( est, lo = NA, hi = NA, rare = NA, true = 1, ... ) {
  
  # sanity checks
  if ( est < 0 ) stop( "HR cannot be negative" )
  
  if ( is.na(rare) ) rare = NULL # for compatibility w/ HR constructor
  
  if ( !inherits(est, "HR") ) est = HR( est, rare = rare )
  if ( !is.na(lo) && !inherits(lo, "HR") ) lo = HR( lo, rare = attr(est, "rare") )
  if ( !is.na(hi) && !inherits(hi, "HR") ) hi = HR( hi, rare = attr(est, "rare") )
  if ( !inherits(true, "HR") ) true = HR( true, rare = attr(est, "rare") )
  
  est = toRR(est)
  if ( !is.na(lo) ) lo = toRR(lo)
  if ( !is.na(hi) ) hi = toRR(hi)
  true = toRR(true)
  
  return( evalues.RR( est = est, lo = lo, hi = hi, true = true ) )
}




evalues.OR = function( est, lo = NA, hi = NA, rare = NA, true = 1, ... ) {
  
  # sanity checks
  if ( est < 0 ) stop( "OR cannot be negative" )
  
  if ( is.na(rare) ) rare = NULL # for compatibility w/ OR constructor
  
  if ( !inherits(est, "OR") ) est = OR( est, rare = rare )
  if ( !is.na(lo) && !inherits(lo, "OR") ) lo = OR( lo, rare = attr(est, "rare") )
  if ( !is.na(hi) && !inherits(hi, "OR") ) hi = OR( hi, rare = attr(est, "rare") )
  if ( !inherits(true, "OR") ) true = OR( true, rare = attr(est, "rare"))
  
  est = toRR(est)
  if ( !is.na(lo) ) lo = toRR(lo)
  if ( !is.na(hi) ) hi = toRR(hi)
  true = toRR(true)
  
  return( evalues.RR( est = est, lo = lo, hi = hi, true = true ) )
}




evalues.RR = function( est, lo = NA, hi = NA, true = 1, ... ) {
  
  # organize user's values
  values = c( est, lo, hi )
  
  # sanity checks
  if ( est < 0 ) stop( "RR cannot be negative" )
  if ( true < 0 ) stop( "True value is impossible" )
  
  # warn user if using non-null true value
  if ( true != 1 ) print(c("You are calculating a \"non-null\" E-value,",
                                 "i.e., an E-value for the minimum amount of unmeasured",
                                 "confounding needed to move the estimate and confidence",
                                 "interval to your specified true value rather than to",
                                 "the null value."))
  
  # check if CI crosses null
  null.CI = NA
  if ( est > true & !is.na( lo ) ) {
    null.CI = ( lo < true )
  }
  
  if ( est < true & !is.na( hi ) ) {
    null.CI = ( hi > true )
  }
  
  
  # sanity checks for CI
  if ( !is.na( lo ) & !is.na( hi ) ) {
    # check if lo < hi
    if ( lo > hi ) stop( "Lower confidence limit should be less than upper confidence limit" )
  }
  
  if ( !is.na( lo ) & est < lo ) stop( "Point estimate should be inside confidence interval" )
  if ( !is.na( hi ) & est > hi ) stop( "Point estimate should be inside confidence interval" )
  
  # compute E-values
  E = sapply( values, FUN = function(x) threshold( x, true = as.numeric(true) ) )
  
  
  # clean up CI reporting
  # if CI crosses null, set its E-value to 1
  if ( !is.na(null.CI) & null.CI == TRUE ){
    E[ 2:3 ] = 1
    print("Confidence interval crosses the true value, so its E-value is 1.") 
  }
  
  # if user provides either CI limit...
  if ( !is.na(lo) | !is.na(hi) ) {
    # ...then only report E-value for CI limit closer to null
    if ( est > true ) E[3] = NA
    if ( est < true ) E[2] = NA
    if ( est == true ) {
      E[2] = 1
      E[3] = NA
    }
  }
  
  result = rbind(values, E)
  
  rownames(result) = c("RR", "E-values")
  colnames(result) = c("point", "lower", "upper")
  class(result) = c("evalue", "matrix")
  
  result
}



twoXtwoRR = function( n11, n10, n01, n00, alpha = 0.05 ){
  
  p1     = n11/( n11 + n10 )
  p0     = n01/( n01 + n00 )
  RR     = p1/p0
  logRR  = log( RR )
  
  selogRR  = sqrt( 1/n11 - 1/( n11+n10 ) + 1/n01 - 1/( n01+n00 ) )
  q.alpha  = qnorm( 1 - alpha/2 )
  
  upperRR  = exp( logRR + q.alpha*selogRR )
  lowerRR  = exp( logRR - q.alpha*selogRR )
  
  res         = c( RR, upperRR, lowerRR )
  names(res)  = c( "point", "upper", "lower" )
  
  return(res) 
}

threshold = function( x, true = 1 ) {
  
  if ( is.na(x) ) return(NA)
  
  if( x < 0 ){
    warning("The risk ratio must be non-negative.")
  }  
  
  if( x <= 1 ){
    x = 1 / x
    true = 1 / true
  }
  
  # standard case: causal effect is toward null
  if ( true <= x ) return( ( x + sqrt( x * ( x - true ) ) ) / true )
  
  # causal effect is away from null
  else if ( true > x ) {
    # ratio that is > 1
    rat = true / x 
    return( rat + sqrt( rat * ( rat - 1 ) ) )
  }
  
}

evalues.RD = function( n11, n10, n01, n00,  
                       true = 0, alpha = 0.05, grid = 0.0001, ... ) {
  
  # sanity check
  if ( any( c(n11, n10, n01, n00) < 0 ) ) stop("Negative cell counts are impossible.")
  
  # sample sizes
  N = n10 + n11 + n01 + n00
  N1 = n10 + n11  # total X=1
  N0 = n00 + n01  # total X=0
  
  # compute f = P(X = 1)
  f = N1 / N
  
  # P(D = 1 | X = 1)
  p1  = n11 / N1
  
  # P(D = 1 | X = 0)
  p0  = n01 / N0
  
  if( p1 < p0 ) stop("RD < 0; please relabel the exposure such that the risk difference > 0.")
  
  
  # standard errors
  se.p1 = sqrt( p1 * ( 1-p1 ) / N1 )
  se.p0 = sqrt( p0 * ( 1-p0 ) / N0 )
  
  # back to Peng's code
  s2.f   = f*( 1-f )/N
  s2.p1  = se.p1^2
  s2.p0  = se.p0^2
  diff   = p0*( 1-f ) - p1*f
  
  # bias factor and E-value for point estimate
  est.BF = ( sqrt( ( true + diff )^2 + 4 * p1 * p0 * f * ( 1-f )  ) - ( true + diff ) ) / ( 2 * p0 * f )
  est.Evalue    = threshold(est.BF)   
  if( p1 - p0 <= true ) stop("For risk difference, true value must be less than or equal to point estimate.")
  
  # compute lower CI limit
  Zalpha        = qnorm( 1-alpha/2 )  # critical value
  lowerCI       = p1 - p0 - Zalpha*sqrt( s2.p1 + s2.p0 )
  
  # check if CI contains null
  if ( lowerCI <= true ) {
    
    # warning( "Lower CI limit of RD is smaller than or equal to true value." )
    return( list( est.Evalue = est.Evalue, lower.Evalue = 1 ) )
    
  } else {
    # find E-value for lower CI limit
    # we know it's less than or equal to E-value for point estimate
    BF.search = seq( 1, est.BF, grid )
    
    # population-standardized risk difference
    RD.search = p1 - p0 * BF.search
    f.search  = f + ( 1-f )/BF.search
    
    # using equation for RD^true on pg 376, compute the lower CI limit for these parameters
    # RD.search * f.search is exactly the RHS of the inequality for RD^true ( population )
    Low.search = RD.search * f.search -
      Zalpha * sqrt( ( s2.p1 + s2.p0 * BF.search^2 ) * f.search^2 +
                       RD.search^2 * ( 1 - 1 / BF.search )^2 * s2.f )
    
    # get the first value for BF_u such that the CI limit hits the true value
    Low.ind    = ( Low.search <= true )
    Low.no     = which( Low.ind==1 )[1]
    lower.Evalue = threshold( BF.search[Low.no] )
    
    
    return(list(est.Evalue   = est.Evalue,
                lower.Evalue = lower.Evalue))
  }
  
}

bias_plot = function( RR, xmax ) {
  
  x = seq( 0, xmax, 0.01 )
  
  # MM: reverse RR if it's preventive
  if ( RR < 1 ) RR = 1/RR
  
  plot( x, x, lty = 2, col = "white", type = "l", xaxs = "i", yaxs = "i", xaxt="n", yaxt = "n",
        xlab = expression( RR[EU] ), ylab = expression( RR[UD] ),
        xlim = c( 0,xmax ),
        main = "" )
  
  x = seq( RR, xmax, 0.01 )
  
  y    = RR*( RR-1 )/( x-RR )+RR
  
  lines( x, y, type = "l" )
  
  
  high = RR + sqrt( RR*( RR-1 ) )
  
  
  points( high, high, pch = 19 )
  
  label5 = seq( 5, 40, by = 5 )
  axis( 1, label5, label5, cex.axis = 1 )
  axis( 2, label5, label5, cex.axis = 1 )
  
  g = round( RR + sqrt( RR * ( RR - 1 ) ), 2 )
  label = paste( "( ", g, ", ", g, " )", sep="" )
  
  text( high + 3, high + 1, label )
  
  legend( "bottomleft", expression(
    RR[EU]*RR[UD]/( RR[EU]+RR[UD]-1 )==RR
  ), 
  lty = 1:2,
  bty = "n" )
  
}


#' @export
evalue.RR = function( est, lo = NA, hi = NA, se = NA, delta = NA, true = 1, ... ){
  evalues.RR(est = est, lo = lo, hi = hi, true = true, ...)
}

#' @export
evalue.OR = function(est, lo = NA, hi = NA, se = NA, delta = NA, true = 1, ...){
  evalues.OR(est = est, lo = lo, hi = hi, true = true, ...)
}

#' @export
evalue.HR = function(est, lo = NA, hi = NA, se = NA, delta = NA,true = 1, ...){
  evalues.HR(est = est, lo = lo, hi = hi, true = true, ...)
}

#' @export
evalue.OLS = function(est, lo = NA, hi = NA, se = NA, delta = 1, true = 0, ...){
  evalues.OLS(est, se = se, delta = delta, true = true, ...)
}

#' @export
evalue.MD = function(est, lo = NA, hi = NA, se = NA, delta = NA, true = 0, ...){
  evalues.MD(est, se = se, true = true, ...)
}


#' @export
evalue.default <- function(est, ...) {
  
  if (is.null(measure) && !inherits(est, "estimate")) stop("Effect measure must be specified")
  
  measure <- class(est)[1]
  
  evalues_func = switch(measure,
                        "HR" = evalues.HR,
                        "OR" = evalues.OR,
                        "RR" = evalues.RR,
                        "RD" = evalues.RD,
                        "OLS" = evalues.OLS,
                        "MD" = evalues.MD)
  
  evalues_func(est, ...)
}

evalue = function( est, lo = NA, hi = NA, se = NA, delta = 1, true = c(0, 1), ... ) {
  UseMethod( "evalue")
}

#' @export
summary.evalue = function( object, ... ) {
  if ( !inherits(object, "evalue")) stop('Argument must be of class "evalue"')
  object[2,1]
}

#' @export
print.evalue = function( x, ... ) {
  class(x) <- "matrix" # to suppress attr printing
  print.default(x)
}

# identify US census region from state code
statecode_to_region = function(statecode){
  #' New England census division
  .new_england <- c("CT", "MA", "ME", "NH", "RI", "VT")
  #' Mid-Atlantic census division
  .mid_atlantic <- c("NJ", "NY", "PA")
  
  #' East North Central census division
  .east_north_central <- c("IL", "IN", "MI", "OH", "WI")
  #' West North Central census division
  .west_north_central <- c("IA", "KS", "MN", "MO", "NE", "ND", "SD")
  
  #' South Atlantic census division
  .south_atlantic <- c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA", "WV")
  #' East South Central census division
  .east_south_central <- c("AL", "KY", "MS", "TN")
  #' West South Central census division
  .west_south_central <- c("AR", "LA", "OK", "TX")
  
  #' Mountain census division
  .mountain <- c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY")
  #' Pacific census division
  .pacific <- c("AK", "CA", "HI", "OR", "WA")
  
  #' Northeast census region
  .northeast_region <- c(.new_england, .mid_atlantic)
  #' North-Central census region
  .north_central_region <- c(.east_north_central, .west_north_central)
  #' Midwest census region
  .midwest_region <- .north_central_region
  #' South census region
  .south_region <- c(.south_atlantic, .east_south_central, .west_south_central)
  #' West census region
  .west_region <- c(.mountain, .pacific)
  
  USregions = c("NE", "MW", "S", "W")
  n = length(statecode)
  out = c()
  for(i in 1:n){
    out[i] = USregions[which(c(statecode[i] %in% .northeast_region,
                               statecode[i] %in% .midwest_region,
                               statecode[i] %in% .south_region,
                               statecode[i] %in% .west_region))]
  }
  
  return(out)
}
