#ifndef MY_FUNCTION
#define MY_FUNCTION

#include<vector>
#include<string>
#include<iostream>
#include<algorithm>
#include<fstream>
#include<map>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>
#include<Eigen/Core>
#include <cstdio>
#include <cstdlib>
#include <iomanip>

#define my_epsilon 0.0000001
using namespace std;

//MC

extern double Rsq, Hhalf; 
//旋转两个端点
void rotate_side_of_chain(vector<double *> &x,vector<double *> &xCtry,int L,int i,double theta,double phi);

void rotate_array_vector(vector<double *> &x,int n);
void rotate_array(vector<double *> &x,int n);
void print_trajectory(vector<double *> &x);

double get_angle(double a[3],double b[3],double c[3]);

double cal_dist_sq(double a[3], double b[3]);

void rotation_around_axis(double axis[3], double theta, double *xin,double *xout);
void rotation_around_axis(double axis[3], double theta, int ibead, int jbead, double xin[][3],double xout[][3]);

void crankshaft(vector<double *> &x,int ibead,int jbead,double theta, vector<double *> &xout);
double check_bead_crash(vector<double *> &x,vector<double *> &xCtry,int ibead,int jbead);

//MC  END



void write_map(map<int,int> &hist_count,fstream &write);

void write_conf(string knottype,vector<double *> &x);

int write_core_tail(vector<int> &knotSize,vector<double *> &point,int n_tail);
int write_core(vector<int> &knotSize,vector<double *> &point);//把相同纽结大小的构象写入同一个文件

int get_segment(vector<int> &segment,int offset,int point);
bool pair_compare(pair<int,double>a,pair<int,double>b);

int judge_triangle(double *a,double *b,double *c,double *plain,double *line_1,double *line_2);

double cal_interSection(double *line1_1,double *line1_2,double *line2_1,double *line2_2,int *up_down);

vector<double> cross(const vector<double> a,const vector<double> b);

void cross_product(double a[3],double b[3], double res[3]);

void cal_normals(double I[3], double J[3], double K[3], double planeijk[4]);

int read_data_links(vector<array<double,3>> &x,fstream &read,vector<int> &ends);

void read_data(vector<double *> & x,char *s);

int read_data_cpp(vector<double *> &x,fstream &read);

void write_data_cpp(vector<double *> &x,string s);

int read_data_cpp_oxdna(vector<double *> &x,fstream &read,int NB);

void clean_pointer(vector<double *> &x);

void print_knot(vector<double *> &points3D);

void error_out(vector<double *> &x,string s="error.txt");

#endif