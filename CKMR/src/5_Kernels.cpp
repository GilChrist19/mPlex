#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************
 * Helpers
 *****************************************************************************/
/**************************************
 * Helper Constants
 *************************************/

/* limits for truncated distributions */
const static double inf_pos = std::numeric_limits<double>::infinity();
//const static double inf_neg = -std::numeric_limits<double>::infinity();

/**************************************
 * Helper Functions
 *************************************/

/* approximate equality */
template<typename T>
static bool approxEqual(T f1, T f2) {
  return (std::fabs(f1 - f2) <= std::numeric_limits<T>::epsilon() * fmax(std::fabs(f1), std::fabs(f2)));
}

/* truncated exponential distribution */
inline double dtruncExp(double x, double r, double a, double b){
  if(a >= b){
    Rcpp::stop("argument a is greater than or equal to b\n");
  }
  double scale = 1.0/r;
  double Ga = R::pexp(a,scale,true,false);
  double Gb = R::pexp(b,scale,true,false);
  if(approxEqual(Ga,Gb)){
    Rcpp::stop("Truncation interval is not inside the domain of the density function\n");
  }
  double density = R::dexp(x,scale,false) / (R::pexp(b,scale,true,false) - R::pexp(a,scale,true,false));
  return density;
}

/******************************************************************************
 * Kernels
 *****************************************************************************/
/**************************************
 * lognormal
 *************************************/

//' Calculate Lognormal Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calcVinEll}},
//' calculate a stochastic matrix where one step movement probabilities follow a lognormal density.
//'
//' The distribution and density functions for the lognormal kernel are given below:
//' \deqn{
//' F(x)=\frac{1}{2} + \frac{1}{2} \mathrm{erf}[\frac{\mathrm{ln}x-\mu}{\sqrt{2}\sigma}]
//' }
//' \deqn{
//' f(x)=\frac{1}{x\sigma\sqrt{2\pi}}\mathrm{exp}\left( -\frac{(\mathrm{ln}x-\mu)^{2}}{2\sigma^{2}} \right)
//' }
//' where \eqn{\mu} is the mean on the log scale, and \eqn{\sigma} is the standard deviation on the log scale.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calcVinEll}}
//' @param meanlog log mean of \code{\link[stats]{Lognormal}} distribution
//' @param sdlog log standard deviation of \code{\link[stats]{Lognormal}} distribution
//'
//' @examples
//' # setup distance matrix
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Ellipsoid  distance formula
//' distMat = calcVinEll(latLongs = latLong)
//'
//' # calculate lognormal distribution over distances
//' #  mean and standard deviation are just for example
//' kernMat = calcLognormalKernel(distMat = distMat, meanlog = 100, sdlog = 10)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcLognormalKernel(const Rcpp::NumericMatrix& distMat,
                                        const double& meanlog, const double& sdlog){

  size_t n = distMat.nrow();
  Rcpp::NumericMatrix kernMat(n,n);

  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      kernMat(i,j) = R::dlnorm(distMat(i,j),meanlog,sdlog,false);
    }
    kernMat(i,_) = kernMat(i,_) / Rcpp::sum(kernMat(i,_)); /* normalize density */
  }

  return kernMat;
};

/**************************************
 * gamma
 *************************************/

//' Calculate Gamma Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calcVinEll}}, calculate a
//' stochastic matrix where one step movement probabilities follow a gamma density.
//'
//' The distribution and density functions for the gamma kernel are given below:
//' \deqn{
//' F(x)=\frac{1}{\Gamma(\alpha)}\gamma(\alpha,\beta x)
//' }
//' \deqn{
//' f(x)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x}
//' }
//' where \eqn{\Gamma(\alpha)} is the Gamma function, \eqn{\gamma(\alpha,\beta x)} is hte lower incomplete
//' gamma function, and \eqn{\alpha,\beta} are the shape and rate parameters, respectively.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calcVinEll}}
//' @param shape shape parameter of \code{\link[stats]{GammaDist}} distribution
//' @param rate rate parameter of \code{\link[stats]{GammaDist}} distribution
//'
//' @examples
//' # setup distance matrix
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Ellipsoid  distance formula
//' distMat = calcVinEll(latLongs = latLong)
//'
//' # calculate gamma distribution over distances
//' #  shape and rate are just for example
//' kernMat = calcGammaKernel(distMat = distMat, shape = 1, rate = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcGammaKernel(const Rcpp::NumericMatrix& distMat, const double& shape, const double& rate){

  size_t n = distMat.nrow();
  Rcpp::NumericMatrix kernMat(n,n);

  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      kernMat(i,j) = R::dgamma(distMat(i,j),shape,rate,false);
    }
    kernMat(i,_) = kernMat(i,_) / Rcpp::sum(kernMat(i,_)); /* normalize density */
  }

  return kernMat;
};

/**************************************
 * exponential
 *************************************/

//' Calculate Exponential Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calcVinEll}}, calculate a
//' stochastic matrix where one step movement probabilities follow an exponential density.
//'
//' The distribution and density functions for the exponential kernel are given below:
//' \deqn{
//' F(x)=1-e^{-\lambda x}
//' }
//' \deqn{
//' f(x)=\lambda e^{-\lambda x}
//' }
//' where \eqn{\lambda} is the rate parameter of the exponential distribution.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calcVinEll}}
//' @param rate rate parameter of \code{\link[stats]{Exponential}} distribution
//'
//' @examples
//' # setup distance matrix
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Ellipsoid  distance formula
//' distMat = calcVinEll(latLongs = latLong)
//'
//' # calculate exponential distribution over distances
//' #  rate is just for example
//' kernMat = calcExpKernel(distMat = distMat, rate = 10)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcExpKernel(const Rcpp::NumericMatrix& distMat, const double& rate){

  size_t n = distMat.nrow();
  Rcpp::NumericMatrix kernMat(n,n);
  double scale = 1.0/rate;

  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      kernMat(i,j) = R::dexp(distMat(i,j),scale,false);
    }
    kernMat(i,_) = kernMat(i,_) / Rcpp::sum(kernMat(i,_)); /* normalize density */
  }

  return kernMat;
};

/**************************************
 * zero-inflated exponential (point mass at zero + zero-truncated exponential distribution)
 *************************************/

//' Calculate Zero-inflated Exponential Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calcVinEll}}, calculate a
//' stochastic matrix where one step movement probabilities follow an zero-inflated
//' exponential density with a point mass at zero. The point mass at zero represents
//' the first stage of a two-stage process, where mosquitoes decide to stay at
//' their current node or leave anywhere. This parameter can be calculated from
//' lifetime probabilities to stay at the current node with the helper function
//' \code{\link[MGDrivE]{calcZeroInflation}}.
//'
//' If a mosquito leaves its current node, with probability \eqn{1-p_{0}}, it
//' then chooses a destination node according to a standard exponential density
//' with rate parameter \eqn{rate}.
//'
//' The distribution and density functions for the zero inflated exponential kernel are given below:
//' \deqn{
//' F(x)=p_{0}\theta(x) + (1-p_{0})(1-e^{-\lambda x})
//' }
//' \deqn{
//' f(x)=p_{0}\delta(x)+(1-p_{0})\lambda e^{-\lambda x}
//' }
//' where \eqn{\lambda} is the rate parameter of the exponential distribution,
//' \eqn{\theta(x)} is the Heaviside step function and \eqn{\delta(x)} is the
//' Dirac delta function.
//'
//' @param distMat Distance matrix from \code{\link[MGDrivE]{calcVinEll}}
//' @param rate Rate parameter of \code{\link[stats]{Exponential}} distribution
//' @param p0 Point mass at zero
//' @param eps Cutoff for extremely small probabilities, default is 1e-20
//'
//' @examples
//' # setup distance matrix
//' # two-column matrix with latitude/longitude, in degrees
//' latLong = cbind(runif(n = 5, min = 0, max = 90),
//'                 runif(n = 5, min = 0, max = 180))
//'
//' # Vincenty Ellipsoid  distance formula
//' distMat = calcVinEll(latLongs = latLong)
//'
//' # calculate hurdle exponential distribution over distances
//' #  rate and point mass are just for example
//' kernMat = calcHurdleExpKernel(distMat = distMat, rate = 1/1e6, p0 = 0.1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcHurdleExpKernel(const Rcpp::NumericMatrix& distMat,
                                        const double& rate, const double& p0,
                                        const double& eps = 1.0e-20){
  // This is the third hurdle exponential function
  // It protects against numerically 0 probabilities
  // It uses a simplified exponential draw, instead of the truncated version
  //  These are almost numerically identical, unless the distances are beneath
  //  the lower truncation bound

  const double scale(1.0/rate); /* convert rate to scale for R::dexp() */
  const double oneMP0(1.0-p0); /* one minus p0, since it doesn't change */

  double holdVal(0.0), sumVal(0.0); /* values for loop */

  const size_t n = distMat.nrow(); /* dimensions */
  Rcpp::NumericMatrix kernMat(n,n); /* initialize return matrix */

  for(size_t i=0; i<n; i++){ /* loop over rows */
    for(size_t j=0; j<n; j++){ /* loop over columns */
      if(i != j){ /* skip diagonal */

        // draw exponential density
        holdVal = R::dexp(distMat(i,j),scale,false);

        // safety, make sure the probability is large enough
        if(holdVal > eps){
          kernMat(i,j) = holdVal;
          sumVal += holdVal;
        }

      } /* end diagonal check */
    } /* end column loop */

    // normalize
    if(sumVal > eps){
      kernMat(i,_) = (kernMat(i,_) / sumVal) * oneMP0; /* normalize density */
      kernMat(i,i) = p0; /* point mass at zero */
      sumVal = 0.0; /* reset sum */
    } else {
      kernMat(i,i) = 1.0; /* all mass at zero */
    }

  } // end row loop

  return kernMat;
}


// This is the original hurdle exponential function
// It is fine for all normal values
// It returns NaNs if all of the probabilities are 0
//  the probs are 0, then renormalization is a divide by 0, then NaNs
// Rcpp::NumericMatrix calcHurdleExpKernelOLD(const Rcpp::NumericMatrix& distMat,
//                                         const double& rate, const double& p0){
//
//   const double a = 1.0e-10; /* lower truncation bound */
//
//   size_t n = distMat.nrow();
//   Rcpp::NumericMatrix kernMat(n,n);
//
//   for(size_t i=0; i<n; i++){
//     for(size_t j=0; j<n; j++){
//       if(i==j){
//         kernMat(i,j) = 0;
//       } else {
//         kernMat(i,j) = dtruncExp(distMat(i,j),rate,a,inf_pos); /* truncated density */
//       }
//     }
//     kernMat(i,_) = (kernMat(i,_) / Rcpp::sum(kernMat(i,_))) *(1-p0); /* normalize density */
//     kernMat(i,i) = p0; /* point mass at zero */
//   }
//
//   return kernMat;
// }


// This is the second hurdle exponential function
// It protects against numerically 0 probabilities, so instead of returning NaNs,
//  it will converge to an identity matrix
// Rcpp::NumericMatrix calcHurdleExpKernel(const Rcpp::NumericMatrix& distMat,
//                                         const double& rate, const double& p0){
//
//   const double a = 1.0e-10; /* lower truncation bound */
//   const double b = 1.0e-20; /* numerical lower bound, "randomly" chosen, no reason*/
//
//   double holdVal(0.0), sumVal(0.0); /* values for loop */
//
//   size_t n = distMat.nrow();
//   Rcpp::NumericMatrix kernMat(n,n);
//
//   for(size_t i=0; i<n; i++){
//     for(size_t j=0; j<n; j++){
//       // skip diagonal
//       if(i != j){
//
//         holdVal = dtruncExp(distMat(i,j),rate,a,inf_pos); /* truncated density */
//
//         // safety, make sure the probability is large enough
//         if(holdVal > b){
//           kernMat(i,j) = holdVal;
//           sumVal += holdVal;
//         }
//
//       } // end diagonal check
//     } // end column loop
//
//     // normalize
//     if(sumVal > b){
//       kernMat(i,_) = (kernMat(i,_) / sumVal) * (1.0-p0); /* normalize density */
//       kernMat(i,i) = p0; /* point mass at zero */
//       sumVal = 0.0; /* reset sum */
//     } else {
//       kernMat(i,i) = 1.0; /* all mass at zero */
//     }
//
//   } // end row loop
//
//   return kernMat;
// }
