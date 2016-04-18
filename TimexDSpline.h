#pragma once
#include <Eigen/Dense>
#include <vector>

// a + b x + c x^2 + d x^3
// a0 + a1 x + a2 x^2 + a3 x^3
struct SplineSet{
    double t1;
    double t2;
    double a0;
    double a1;
    double a2;
    double a3;
};
typedef std::vector<double> vec;

typedef std::vector<SplineSet> SplineVec;

class TimexDSpline {
public:
    TimexDSpline(const std::vector<Eigen::VectorXd>& timePoints);
    std::vector<SplineVec> getCoeffs();
    double pdistance2(const vec& p1, const vec& p2); 
    double splineEval(double x, const SplineVec& coeffVec);

private:
    std::vector<SplineVec> coeffSet;
    void spline(const vec &x, const vec &y, SplineVec& coeff);
    int pieces;
    double tolerance;
};

