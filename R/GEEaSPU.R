#' GEEaSPU test.
#'
#' It gives p-values of the GEESPU tests and GEEaSPU test.
#'
#' @param traits trait matrix. The row for individuals and the column for traits.
#'
#' @param geno A matrix of genetic information.
#'
#' @param Z covariates.
#'
#' @param pow power used in GEEaSPU test. A vector of the powers.
#'
#' @param n.sim number of simulations
#'
#' @return p-values for the GEE-SPU and GEE-aSPU test.
#'
#' @author Junghi Kim, Wei Pan and Il-Youp Kwak
#'
#' @references
#'
#' Il-Youp Kwak, Wei Pan (2015)
#' Adaptive Gene- and Pathway-Trait Association Testing with GWAS Summary Statistics,
#' in revision.
#'
#' Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014)
#' A powerful and adaptive association test for rare variants,
#' Genetics, 197(4), 1081-95
#'
#' Junghi Kim, Jeffrey R Wozniak, Bryon A Mueller, Xiaotong Shen and Wei Pan (2014) Comparison of statistical tests for group differences in brain functional networks, NeuroImage, 1;101:681-694
#'
#' @examples
#'
#' traits <- matrix(rnorm(100*5, 0,1), ncol=5)
#' Z <- rnorm(100, 2, 0.5)
#' geno <- rbinom(100, 2, 0.5)
#' out <- GEEaSPU(traits, geno, Z = NULL,family = "gaussian", 
#'		  gamma = c(1:8,Inf), n.perm = 1000)
#'

GEEaSPU <- function(traits, geno, 
	  Z = NULL, family = c("binomial", "gaussian"), 
	  gamma = c(1:8,Inf), n.perm=10000) {

   Score <- GEE(traits = traits, geno = geno, Z = Z,family = family)
   U <- Score$U
   V <- Score$Cov
   GEEspu.score(U = U, V = V, gamma = gamma, B = n.perm)

}

