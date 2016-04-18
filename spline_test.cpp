#include "TimexDSpline.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;
int main()
{
    std::ifstream  data("tmp.log");
    std::string line;

    /*
    vector<Eigen::VectorXd> dataPoints(3000);
    int i = 0;
    int j;
    while(std::getline(data,line))
    {
        std::stringstream  lineStream(line);
        std::string        cell;
        j = 0;
        dataPoints[i].resize(7);
        while(std::getline(lineStream,cell,','))
        {
            // You have a cell!!!!
            dataPoints[i][j] = atof(cell.c_str());
            j++;
        }
        i++;
    }*/
    vector<Eigen::VectorXd> dataPoints(4);

    dataPoints[0].resize(4);
    dataPoints[1].resize(4);
    dataPoints[2].resize(4);
    dataPoints[3].resize(4);
    dataPoints[0] << 0,0,0,0;
    dataPoints[1] << 1 ,1 ,1 ,1;
    dataPoints[2] << 1.5, 1,2,3;
    dataPoints[3] << 2, 3,2,3;

    TimexDSpline stuff(dataPoints);
    std::vector<SplineVec> coeffs;
    coeffs = stuff.getCoeffs();
    for ( int i = 0; i < coeffs.size(); i++) 
    {
        cout << "axis " << i << endl;
        for ( int j = 0; j < coeffs[i].size(); j++ )
        {
            cout << coeffs[i][j].t1 << " " << coeffs[i][j].t2 <<" " << coeffs[i][j].a0 <<" " << coeffs[i][j].a1 <<" " << coeffs[i][j].a2 <<" " << coeffs[i][j].a3 ;
            cout << endl;
        }
    }
    //std::cout << coeffs[0];
}
