#include <math.h>
#include <limits>
#include <Rcpp.h>

using namespace Rcpp;


/* limits for truncated distributions */
const static double inf_pos = std::numeric_limits<double>::infinity();
const static double inf_neg = -std::numeric_limits<double>::infinity();
const static double R_earth = 6371.0; // Earth mean radius [km]

/* Convert degrees to radians */
inline double deg2rad(const double& deg){return deg*M_PI/180.0;};


/* Calculates the geodesic distance between two points specified by
 * radian latitude/longitude using the Haversine formula
 * Ouputs distance between sites 1 and 2 as meters
*/

inline double gcd_hf(const double& long1, const double& lat1, const double& long2, const double& lat2){
  double deltaLong = (long2 - long1);
  double deltaLat = (lat2 - lat1);
  double a = std::pow(std::sin(deltaLat/2.0),2) + std::cos(lat1) * std::cos(lat2) * std::pow(std::sin(deltaLong/2.0),2);
  double sqrtA = std::sqrt(a);
  double c = 2.0 * std::asin(fmin(1.0,sqrtA));
  double d = (R_earth * c)*1000.0;
  return d;
};

/* Fxn to calculate matrix of distances between each two sites
 * INPUT: a matrix in which longs are in first column and lats in second column
 * OUTPUT: a distance matrix between all pairwise sites
 * Output distances are in meters
*/

//' Calculate Haversine Distance
//'
//' Calculate the great-circle distance (Haversine distance) between sets of longitude - latitude points, see \url{https://en.wikipedia.org/wiki/Haversine_formula}
//'
//' @param latlongs numeric matrix where first column is vector of latitudes and second column is vector of longitudes
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calc_haversine(const Rcpp::NumericMatrix& latlongs){
  size_t n = latlongs.nrow();
  Rcpp::NumericMatrix zz = Rcpp::NumericMatrix(n,n);
  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      zz(i,j) = gcd_hf(deg2rad(latlongs(i,1)), deg2rad(latlongs(i,0)),deg2rad(latlongs(j,1)), deg2rad(latlongs(j,0)));
    }
  }
  return zz;
};

/* calculate lognormal kernel */

//' Calculate Lognormal Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calc_haversine}}, calculate a stochastic matrix where one step movement probabilities follow a lognormal density.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calc_haversine}}
//' @param meanlog log mean of \code{\link[stats]{Lognormal}} distribution
//' @param sdlog log standard deviation of \code{\link[stats]{Lognormal}} distribution
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calc_LognormalKernel(const Rcpp::NumericMatrix& distMat, const double& meanlog, const double& sdlog){

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

/* calculate gamma kernel */

//' Calculate Gamma Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calc_haversine}}, calculate a stochastic matrix where one step movement probabilities follow a gamma density.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calc_haversine}}
//' @param shape shape parameter of \code{\link[stats]{GammaDist}} distribution
//' @param rate rate parameter of \code{\link[stats]{GammaDist}} distribution
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calc_GammaKernel(const Rcpp::NumericMatrix& distMat, const double& shape, const double& rate){

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

/* calculate exponential kernel */

//' Calculate Exponential Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calc_haversine}}, calculate a stochastic matrix where one step movement probabilities follow an exponential density.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calc_haversine}}
//' @param r rate parameter of \code{\link[stats]{Exponential}} distribution
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calc_ExpKernel(const Rcpp::NumericMatrix& distMat, const double& r){

  size_t n = distMat.nrow();
  Rcpp::NumericMatrix kernMat(n,n);
  double scale = 1.0/r;

  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      kernMat(i,j) = R::dexp(distMat(i,j),scale,false);
    }
    kernMat(i,_) = kernMat(i,_) / Rcpp::sum(kernMat(i,_)); /* normalize density */
  }

  return kernMat;
};

/* calculate hurdle exponential kernel (point mass at zero + zero-truncated exponential distribution) */

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

//' Calculate Hurdle Exponential Stochastic Matrix
//'
//' Given a distance matrix from \code{\link[MGDrivE]{calc_haversine}}, calculate a stochastic matrix where one step movement probabilities follow an zero-truncated exponential density with a point mass at zero.
//'
//' @param distMat distance matrix from \code{\link[MGDrivE]{calc_haversine}}
//' @param r rate parameter of \code{\link[stats]{Exponential}} distribution
//' @param pi point mass at zero
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calc_HurdleExpKernel(const Rcpp::NumericMatrix& distMat, double r, double pi){
  const double a = 1.0e-10; /* lower truncation bound */

  size_t n = distMat.nrow();
  Rcpp::NumericMatrix kernMat(n,n);

  for(size_t i=0; i<n; i++){
    for(size_t j=0; j<n; j++){
      if(i==j){
        kernMat(i,j) = 0;
      } else {
        kernMat(i,j) = dtruncExp(distMat(i,j),r,a,inf_pos); /* truncated density */
      }
    }
    kernMat(i,_) = (kernMat(i,_) / Rcpp::sum(kernMat(i,_))) *(1-pi); /* normalize density */
    kernMat(i,i) = pi; /* point mass at zero */
  }

  return kernMat;
}
