
#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>

#include "TimexDSpline.h"

//typedef vector<double> vec;
using namespace std;

// t, x, y , z
TimexDSpline::TimexDSpline(const vector<Eigen::VectorXd>& timePoints):tolerance(1.), 
    pieces(0)
{
    pieces = timePoints.size() - 1;

    int i,j;
    int tsize = timePoints.size();
    int axisSize = timePoints[0].rows()-1;

    vec timeStamp(timePoints.size());
    for (j = 0; j < timePoints.size(); j++ ) {
        timeStamp[j] = timePoints[j][0];
    }
    vector<bool> timeStampUsed(tsize);
    timeStampUsed[0] = true;
    timeStampUsed[tsize-1] = true;
    for ( i = 1; i < tsize-1; i++ )
    {
        timeStampUsed[i] = false;
    }
    // initialized done

    double maxDist = -1;
    int maxIndex = -1;
    while ( maxDist == -1 || maxDist > tolerance ) 
    {
        maxDist = 0;
        maxIndex = -1;

        std::vector<SplineVec> coeffWorkingSet(axisSize);
        vec realTimestamp;
        vector<vec> lines(axisSize);
        for ( i = 0; i < tsize; i++ )
        {
            if ( timeStampUsed[i] ) 
            {
                realTimestamp.push_back(timeStamp[i]);
                for ( j = 0; j < axisSize; j++ ) {
                    lines[j].push_back(timePoints[i][j+1]);
                }
            }
        }
        if ( realTimestamp.size() == timeStamp.size() )
            break;
        for (int k = 0; k < realTimestamp.size(); ++k)
            std::cout << realTimestamp[k] << ' ';
        for ( i = 0; i < axisSize; i++ ) {
            spline(realTimestamp, lines[i], coeffWorkingSet[i]);
        }
        for ( i = 1; i < pieces; i++ ) 
        {
            double t = timePoints[i][0];
            vec originalPoint(axisSize);
            vec newPoint(axisSize);
            std::copy(&timePoints[i](1), &timePoints[i](1)+timePoints[i].size()-1, originalPoint.begin());
            for ( j = 0; j < axisSize; j++ ) 
            {
                newPoint[j] = splineEval(t, coeffWorkingSet[j]);
            }
            double distance2 = pdistance2(originalPoint, newPoint);

            if ( distance2 > maxDist ) {
                maxDist = distance2;
                maxIndex = i;
            }
        }
        assert(timeStampUsed[maxIndex] == false);
        timeStampUsed[maxIndex] = true;
        coeffSet = coeffWorkingSet;
        //cout << "max dist is " << maxDist << " at index " << maxIndex << endl;
    }

}

//vector<SplineSet> spline(vec &x, vec &y)
void TimexDSpline::spline(const vec &x, const vec &y, SplineVec& coeff)
{
    int n = x.size()-1;
    assert(x.size() == y.size());
    vec a;
    a.insert(a.begin(), y.begin(), y.end());
    vec b(n);
    vec d(n);
    vec h;

    for(int i = 0; i < n; ++i)
        h.push_back(x[i+1]-x[i]);

    vec alpha;
    for(int i = 0; i < n; ++i)
        alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );

    vec c(n+1);
    vec l(n+1);
    vec mu(n+1);
    vec z(n+1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for(int i = 1; i < n; ++i)
    {
        l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for(int j = n-1; j >= 0; --j)
    {
        c[j] = z [j] - mu[j] * c[j+1];
        b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
        d[j] = (c[j+1]-c[j])/3/h[j];
    }

    vector<SplineSet> output_set(n);
    for(int i = 0; i < n; ++i)
    {
        output_set[i].a0 = a[i];
        output_set[i].a1 = b[i];
        output_set[i].a2 = c[i];
        output_set[i].a3 = d[i];
        output_set[i].t1 = x[i];
        output_set[i].t2 = x[i+1];
    }
//return output_set;
    coeff = output_set;
}

double TimexDSpline::pdistance2(const vec& p1, const vec& p2) 
{
    double distanceSq, distance;
    assert(p1.size() == p2.size());
    distanceSq = 0;
    for ( int i = 0; i < p1.size(); i++ ) 
    {
        distanceSq += pow((p2[i] - p1[i]),2);
    }
    return distanceSq;
}

double TimexDSpline::splineEval(double x, const SplineVec& coeffVec)
{
    SplineSet coeff;
    int i;
    int found = 0;
    
    assert(coeffVec.size() > 0 );
    for ( i = 0; i< coeffVec.size(); i++)
    {
        if ( coeffVec[i].t1 <= x && x <= coeffVec[i].t2 )
        {
            coeff = coeffVec[i];
            found = 1;
            break;
        }
    }
    assert(found);

    double x0 = x - coeff.t1;
    double value = coeff.a0 + coeff.a1* x0 + coeff.a2 * pow(x0,2) + coeff.a3 * pow(x0,3);
    return value;
}

std::vector<SplineVec> TimexDSpline::getCoeffs()
{
    return coeffSet;
}

