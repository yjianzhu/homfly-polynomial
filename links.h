#include<vector>
#include<string>
#include<iostream>
#include<algorithm>
#include<cmath>
#include<fstream>
#include<map>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include "homfly.h"
#include "my_function.h"

using namespace std;
using namespace Eigen;
#define epsilon 0.0000001

int judge_triangle(array<double,3> &a,array<double,3> &b,array<double,3> &c,double *plain,array<double,3> &line_1,array<double,3> &line_2);


double cal_interSection(array<double,3> &line1_1,array<double,3> &line1_2,array<double,3> &line2_1,array<double,3> &line2_2,int *up_down);


void cross_product_hull(double a[3],double b[3], double res[3]);

void cal_normals(array<double,3> I, array<double,3> J, array<double,3> K, double planeijk[4]) ;

// class definition
class links
{
private:
    /* data */
public:
    int num_beads=0;
    vector<vector<array<double,3>>> strings;

    vector<vector<double>> interSection_Matrix;

    vector<vector<int>> interSection_Matrix_up_down;

    links(vector<array<double,3>> a,vector<int> b/* args */);

    links(vector<array<double,3>> a)
    {
        strings.push_back(a);
    }
    void print_links(){
        cout<<"there are "<<strings.size()<<"chains!"<<endl;
        for(int i=0;i<strings.size();i++)
        {
            for(int j=0;j<strings[i].size();j++)
                printf("%10.5lf\t%10.5lf\t%10.5lf\t\n",strings[i][j][0],strings[i][j][1],strings[i][j][2]);
        }
    }
    ~links();

    void cal_number_of_beads()
    {
        num_beads=0;
        for(int i=0;i<strings.size();i++)
        {
            num_beads+=strings[i].size();
        }
        cout<<"there are "<<num_beads<<" beads!"<<endl;
    };


    void KMT()
    {
        for(int j=0;j<strings.size();j++)
        {
            // j表示第几条链，i，k表示第几个点。
            while (true)
            {
                int number=strings[j].size();
                int flag;
                for(int i=1;i<=strings[j].size();i++)
                {
                    //cout<<i<<'\t'<<strings.size()<<endl;
                    double plain[4]={0};
                    flag=0;
                    cal_normals(strings[j][i-1],strings[j][i%strings[j].size()],strings[j][(i+1)%strings[j].size()],plain);

                    if(fabs(plain[0])<epsilon and fabs(plain[1])<epsilon and fabs(plain[2])<epsilon)
                    {//如果三个点在一条线上，可以省去
                        strings[j].erase(strings[j].begin()+i%strings[j].size());
                        i--;
                        //cout<<i<<endl;
                        //cout<<strings.size()<<endl;
                        continue;
                    }

                    for(int k=i+1;k<=strings[j].size()+i-2;k++)
                    {
                        if(judge_triangle(strings[j][i-1],strings[j][i%strings[j].size()],strings[j][(i+1)%strings[j].size()],plain,strings[j][k%strings[j].size()],strings[j][(k+1)%strings[j].size()]))
                        {
                            flag=1;
                            break; 
                        }
                    }

                    //对其他链做判断
                    for(int i_chain=0;i_chain<strings.size();i_chain++)
                    {
                        if(i_chain==j)
                            continue;
                        for(int k=0;k<strings[i_chain].size();k++)
                        {
                            if(judge_triangle(strings[j][i-1],strings[j][i%strings[j].size()],strings[j][(i+1)%strings[j].size()],plain,strings[i_chain][k%strings[i_chain].size()],strings[i_chain][(k+1)%strings[i_chain].size()]))
                            {
                                flag=1;
                                break; 
                            }
                        }
                    }

                    if(flag==0)
                    {
                        strings.erase(strings.begin()+i%strings.size());
                        i--;
                        //cout<<i<<endl;
                        //cout<<strings.size()<<endl;
                    }

                }
                if(number==strings[j].size())
                    break;
            }
        }   
    }


    void get_interSection_Matrix()
    {
        int num_chains=strings.size();
        int n=0;
        int up_down;
        for(int i=0;i<num_chains;i++)
            n+=strings[i].size();

        vector<vector<double>> i_Matrix(n,vector<double> (n,0));
        vector<vector<int>> i_Matrix_up_down(n,vector<int> (n,0));

        int i=0,j=0;    // 矩阵下标

        for(int i_chain=0;i_chain<num_chains;i_chain++)
        {
            for(int i_bead=0;i_bead<strings[i_chain].size();i_bead++,i++)
            {
                j=0;
                for(int j_chain=0;j_chain<num_chains;j_chain++)
                {
                    for(int j_bead=0;j_bead<strings[j_chain].size();j_bead++,j++)
                    {
                        if(i==j)    
                            continue;
                        up_down=0;
                        i_Matrix[i][j]=cal_interSection(strings[i_chain][i_bead],strings[i_chain][(i_bead+1)%strings[i_chain].size()],strings[j_chain][j_bead],strings[j_chain][(j_bead+1)%strings[j_chain].size()],&up_down);
                        i_Matrix_up_down[i][j]=up_down; 
                    }
                }
            }
        }


        // 删除
        // for(int i=0;i<n;i++)
        // {
        //     for(int j=0;j<n;j++)
        //     {   
        //         if(i==j)    
        //             continue;
        //         up_down=0;
        //         i_Matrix[i][j]=cal_interSection(strings[i],strings[(i+1)%n],strings[j],strings[(j+1)%n],&up_down);
        //         i_Matrix_up_down[i][j]=up_down; 
        //     }
        //}
        interSection_Matrix=i_Matrix;
        interSection_Matrix_up_down=i_Matrix_up_down;
    }

    void print_IM()
    {
        printf("\n Intersection Matrix:\n");
        int n=num_beads;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                printf("%5.10f\t",interSection_Matrix[i][j]);
            }
            printf("\n");
        }

        printf("\n Intersection Matrix up or down:\n");
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                printf("%5d\t",interSection_Matrix_up_down[i][j]);
            }
            printf("\n");
        }
    }

    void get_gauss_notation()
    {
        int cross_num=0;
        string notation;
        notation+=to_string(strings.size())+"\n\n";
        map<pair<int,int>,pair<int,int>> cross;
        int n=strings.size();
        int count=0;

        for(int i=0;i<num_beads;i++)
        {
            for(int j=i+1;j<num_beads;j++)
            {
                if(fabs(interSection_Matrix[i][j])>epsilon)
                    {
                        int isOver0=interSection_Matrix[i][j]>0.?1:-1;//大于0的记为1，小于0记为-1
                        cross.insert(map<pair<int,int>,pair<int,int>>::value_type (pair<int,int>(i,j), pair<int,int>(count,isOver0)));
                        count++;
                    }
            }
        }

        // 按照链的顺序输出
        int cur_index=0;    //表示键的index
        for(int i_chain=0;i_chain<strings.size();i_chain++)
        {
            int cross_num_single=0;//和temp_count有点重复

            int temp_count=0;// 只用当前链上的交叉点序号
            string notation_single{};
            for(int index=0;index<strings[i_chain].size();index++,cur_index++)
            {
                vector<pair<int,double>> oneLine;
                for(int j=0;j<num_beads;j++)
                {
                    if(fabs(interSection_Matrix[cur_index][j])>epsilon)
                        oneLine.push_back(pair<int,double>(j,fabs(interSection_Matrix[cur_index][j])));
                }
                sort(oneLine.begin(),oneLine.end(),pair_compare);

                
                for(auto it=oneLine.begin();it!=oneLine.end();it++)
                {   
                    pair<int,int> temp;
                    if(cur_index<it->first) {temp.first=cur_index;temp.second=it->first;}
                    else {temp.first=it->first;temp.second=cur_index;}
                    //notation_single+=to_string(cross[temp].first)+' ';
                    notation_single+=to_string(temp_count)+' ';
                    notation_single+=to_string(interSection_Matrix_up_down[cur_index][it->first])+' ';
                    temp_count++;
                }
                cross_num_single+=oneLine.size();
            }
            notation+=to_string(cross_num_single)+'\n';
            notation+=notation_single+'\n';
        }


        // notation+=to_string(cross.size()*2)+' ';
        // cross_num+=cross.size();
        // for(int i=0;i<n;i++)
        // {   
        //     vector<pair<int,double>> oneLine;
        //     for(int j=0;j<n;j++)
        //     {
        //         if(fabs(interSection_Matrix[i][j])>epsilon)
        //             oneLine.push_back(pair<int,double>(j,fabs(interSection_Matrix[i][j])));
        //     }
        //     sort(oneLine.begin(),oneLine.end(),pair_compare);
        //     for(auto it=oneLine.begin();it!=oneLine.end();it++)
        //     {   
        //         pair<int,int> temp;
        //         if(i<it->first) {temp.first=i;temp.second=it->first;}
        //         else {temp.first=it->first;temp.second=i;}
        //         notation+=to_string(cross[temp].first)+' ';
        //         notation+=to_string(interSection_Matrix_up_down[i][it->first])+' ';
        //     }
        // }

        for(auto i:cross)
        {
            notation+=to_string(i.second.first)+' '+to_string(i.second.second)+' ';
        }
        cout<<notation<<endl;

        char * writable = new char[notation.size() + 1];
        copy(notation.begin(), notation.end(), writable);
        writable[notation.size()] = '\0';
        cout<<homfly_str(writable)<<endl;
        delete[] writable;
    }
};



// 待使用
void move3D(vector<double *> a,vector<int> b);

