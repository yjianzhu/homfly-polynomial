#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include "homfly.h"
#include "my_function.h"
#include "links.h"

using namespace std;
using namespace Eigen;
#define epsilon 0.0000001

int main(int argc, char *argv[])
{
    fstream read;
    read.open(argv[1],ios::in);
    vector<array<double,3>> points;
    vector<int> ends;

    int count=0;
    while(!read.eof())
    {   
        if(!read_data_links(points,read,ends)) 
            break;

        // 
        links L(points,ends);
        L.print_links();
        L.cal_number_of_beads();   
        L.KMT();
        L.print_links();
        L.cal_number_of_beads();

        L.get_interSection_Matrix();
        cout<<L.num_beads<<endl;
        L.print_IM();

        L.get_gauss_notation();

        count++;
        points.clear();
        ends.clear();
    }
}