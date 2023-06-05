#include "links.h"

int judge_triangle(array<double,3> &a,array<double,3> &b,array<double,3> &c,double *plain,array<double,3> &line_1,array<double,3> &line_2)
{
    // a b c 三角形三个顶点。 line_1 and line_2 是向量
    // 输出1 表示有重叠
    double dis1=line_1[0]*plain[0]+line_1[1]*plain[1]+line_1[2]*plain[2]+plain[3];
    double dis2=line_2[0]*plain[0]+line_2[1]*plain[1]+line_2[2]*plain[2]+plain[3];
    if( dis1*dis2>epsilon)
        return 0;
    //上面三行代码判断是否为平面同侧。
    if(fabs(dis1)<epsilon and fabs(dis2)<epsilon  )
    {
        //两个点都在面上，暂时先不给用
        return 1;
    }
    // if( fabs(dis1)<epsilon )
    // {
    //     Matrix2d ab;
    //     ab<<b[0]-a[0],c[0]-a[0],b[1]-a[1],c[1]-a[1];
    //     Vector2d b0,x0;
    //     b0<<line_1[0]-a[0],line_1[1]-a[1];
    //     x0=ab.inverse()*b0;
    //     if(x0(0,0)<0 or x0(1,0)<0 or x0(0,0)+x0(1,0)>1)
    //         return 0;
    //     else return 1;
    // }
    // if( fabs(dis2)<epsilon )
    // {
    //     Matrix2d ab;
    //     ab<<b[0]-a[0],c[0]-a[0],b[1]-a[1],c[1]-a[1];
    //     Vector2d b0,x0;
    //     b0<<line_2[0]-a[0],line_2[1]-a[1];
    //     x0=ab.inverse()*b0;
    //     if(x0(0,0)<-epsilon or x0(1,0)<-epsilon or x0(0,0)+x0(1,0)>1+epsilon)
    //         return 0;
    //     else return 1;
    // }
    
    

    Vector3d T,d,E1,E2,M,K;
    T<<line_1[0]-a[0],line_1[1]-a[1],line_1[2]-a[2];
    d<<line_2[0]-line_1[0],line_2[1]-line_1[1],line_2[2]-line_1[2];
   
    E1<<b[0]-a[0],b[1]-a[1],b[2]-a[2];
    E2<<c[0]-a[0],c[1]-a[1],c[2]-a[2];
    M=d.cross(E2);
    double det=M.dot(E1);
    K=T.cross(E1);
    double t=K.dot(E2)/det;
    double u=M.dot(T)/det;
    double v=K.dot(d)/det;
    
    if(u<=0 or v<=0 or u+v>=1 or t>=1 or t<=0)
        return 0;
    //cout<<det<<"\tuvt\t"<<u<<'\t'<<v<<'\t'<<t<<endl;
    return 1;    
}


double cal_interSection(array<double,3> &line1_1,array<double,3> &line1_2,array<double,3> &line2_1,array<double,3> &line2_2,int *up_down)
{
    Vector2d k,b;
    b<<line2_1[0]-line1_1[0],line2_1[1]-line1_1[1];
    Matrix2d x;
    x<<line1_2[0]-line1_1[0],line2_1[0]-line2_2[0],line1_2[1]-line1_1[1],line2_1[1]-line2_2[1];
    if(x.determinant()==0)  return 0;//平行
    k=x.inverse()*b;
    if(k(0,0)<0+epsilon or k(0,0)>1-epsilon or k(1,0)<0+epsilon or k(1,0)>1-epsilon)    return 0;//用一个小量解决可能出现的共点情况

    if((k(0,0)*(line1_2[2]-line1_1[2])+line1_1[2])>(k(1,0)*(line2_2[2]-line2_1[2])+line2_1[2]))
        {//第一个片段在上方
            *up_down=1;
            if(x(0,0)*(-x(1,1))+x(1,0)*x(0,1)>0)
                return k(0,0);
            else
                return -1*k(0,0);
        }
    else
        {//第二个片段在下方
            *up_down=-1;
            if(x(0,0)*(-x(1,1))+x(1,0)*x(0,1)>0)
                return -1*(k(0,0));
            else
                return k(0,0);
        }
}

void cross_product_hull(double a[3],double b[3], double res[3]) {
   res[0] = a[1]*b[2] - a[2]*b[1];
   res[1] = a[2]*b[0] - a[0]*b[2];
   res[2] = a[0]*b[1] - a[1]*b[0];
}
void cal_normals(array<double,3> I, array<double,3> J, array<double,3> K, double planeijk[4]) {
   int d;
   double  vij[3], vik[3], vtmp[3], vabs;

   for(d=0;d<3;d++) vij[d] = J[d] - I[d];
   for(d=0;d<3;d++) vik[d] = K[d] - I[d];
   cross_product_hull(vij, vik, vtmp);
   vabs = sqrt( vtmp[0]*vtmp[0] + vtmp[1]*vtmp[1] + vtmp[2]*vtmp[2] );
   if(fabs(vabs)<epsilon) return;
   for(d=0;d<3;d++) planeijk[d] = vtmp[d]/vabs;
   planeijk[3]  = 0.;
   for(d=0;d<3;d++)  planeijk[3] += I[d]*planeijk[d];
   planeijk[3]  = -planeijk[3];

}


void move3D(vector<array<double,3>> a,vector<int> b)
{//不能用
    int N_components=b.size()+1;
    b.push_back(a.size());
    vector<array<double,3>> single_knot;
    single_knot.assign(a.begin(),a.begin()+b[0]);
    // knot temp(single_knot,b);
    // temp.print_knot();


    //接口留给links，目前只用单一闭合链。
    for(int component=1;component<N_components;component++)
    {
        vector<array<double,3>> single_knot;
        printf("%d\n",b[component]);
        single_knot.assign(a.begin()+b[component-1],a.begin()+b[component]);
        links temp(single_knot);
        temp.print_links();
    }


}

links::~links()
{
    //printf("class has been killed!\n");
}

links::links(vector<array<double,3>> a,vector<int> b/* args */)
{
    if(b.size()==0)
    {
        strings.push_back(a);
    }
    else
    {
        vector<array<double,3>> single_knot;
        single_knot.assign(a.begin(),a.begin()+b[0]);
        strings.push_back(single_knot);

        for(int i=1;i<b.size();i++)
        {
            vector<array<double,3>> single_knot;
            single_knot.assign(a.begin()+b[i-1],a.begin()+b[i]);
            strings.push_back(single_knot);
        }
        single_knot.assign(a.begin()+b[b.size()-1],a.end());
        strings.push_back(single_knot);
    }

}