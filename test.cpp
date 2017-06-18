#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include <iostream>
#include <vector>
#include "Vasicek.h"


TEST_CASE("Test computeExpectation: should equal long run", "[Vasicek]"){
    std::vector<double> y0(1);
    std::vector<double> alpha(1);
    std::vector<double> beta(1);
    y0[0]=.5;
    alpha[0]=.3;
    beta[0]=y0[0];
    REQUIRE(vasicek::computeExpectation(y0, alpha, beta, 1, 1.0)[0]==Approx(beta[0])); 
}
TEST_CASE("Test computeIntegralExpectation", "[Vasicek]"){
    std::vector<double> y0(1);
    std::vector<double> alpha(1);
    std::vector<double> beta(1);
    double tau=1.0;
    y0[0]=.5;
    alpha[0]=.3;
    beta[0]=y0[0];
    REQUIRE(vasicek::computeIntegralExpectation(y0, alpha, beta, 1, tau)[0]==Approx(tau*beta[0]+(1-exp(-alpha[0]*tau))*(y0[0]-beta[0])/alpha[0])); 
}
TEST_CASE("Test computeIntegralVarianceVasicek", "[Vasicek]"){
    std::vector<double> alpha(1);
    std::vector<double> sigma(1);
    std::vector<std::vector<double>> rho;
    std::vector<double> myTmp(1);

    myTmp[0]=1;
    rho.emplace_back(myTmp);
    double tau=1.0;
    sigma[0]=.2;
    alpha[0]=.3;
    rho[0][0]=1;
    double elem=(1.0-exp(-alpha[0]*tau))/alpha[0];
    double anotherelem=(1.0-exp(-2*alpha[0]*tau))/(2*alpha[0]);
    double coefElem=(sigma[0]*sigma[0])/(alpha[0]*alpha[0]);
    REQUIRE(vasicek::computeIntegralVarianceVasicek(alpha, sigma, rho, 1, tau)[0][0]==Approx(coefElem*(tau-2*elem+anotherelem))); 
}
