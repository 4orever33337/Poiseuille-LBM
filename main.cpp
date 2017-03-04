//泊肃叶流非平衡外推+修正反弹
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
const int Q=9;            //D2Q9模型
const int NX=70;         //X方向网格数
const int NY=70;          //Y方向网格数
const double U=0.1;

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};   //表示D2Q9离散速度模型速度空间的速度配置
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};     //表示权系数
int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

double rho[NX+1][NY+1];   //密度ρ
double u0[NX+1][NY+1][2]; //n时层的速度
double u[NX+1][NY+1][2];  //n+1时层的速度,为宏观速度
double f[NX+1][NY+1][Q];  //演化前的分布函数，256个格子，257个点
double f1[NX+1][NY+1][Q];  //碰撞后的分布函数，中间临时变量
double F[NX+1][NY+1][Q];  //演化后的分布函数
double solid_flag[NX+1][NY+1];//固体流体区域划分标志，固体为1，流体为0

int i,j,k,ip,jp,n;		  //ip,jp代表前一个时层
double c;                 //格子速度
double Re;                //雷诺数
double dx,dy;             //x,y方向的网格步长
double Lx,Ly;             //x,y方向的长度
double dt;                //时间步长
double rho0;              //流场初始密度ρ0
double P0;                //流场初始压力P0
double tau_f;             //无量纲松弛时间τ
double niu;               //运动粘度系数ν
double error;             //两个相邻时层速度的最大相对误差


//函数声明
void init();
double feq(int k,double rho,double u[2]);
void evolution();
void output(int m);
void Error();

int main()
{
	using namespace std;
	cout<<4<<endl;
	init();//初始化
	output(1);
	for(n=0;;n++)
	{
		evolution();//进行演化
		if(n%100==0)
		{
			Error();
			cout<<"The"<<n<<"th computation result:"<<endl<<"The u,v of point (NX/2,NY/2)is:"<<setprecision(6)<<u[NX/2][NY/2][0]<<","<<u[NX/2][NY/2][1]<<endl;
			cout<<"The max relative error of uv is:"<<setiosflags(ios::scientific)<<error<<endl;
			if(n>=1000)
			{
				if(n%100==0) output(n);
				if(error<1.0e-6) break;
			}
		}
	}
	return 1;
}

void init()
{
	dx=1.0;
	dy=1.0;
	Lx=dx*double(NX);
	Ly=dy*double(NY);
	dt=dx;                 
	c=dx/dt;               
	rho0=1.0;
	Re=100;
	niu=U*Lx/Re;           
	tau_f=3.0*niu+0.5;     //由v＝cs^2*(τ-0.5)
//	tau_f=1.0;
//	P0=1.0e-5;             //外力F
	std::cout<<"tau_f="<<tau_f<<endl;
	for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			solid_flag[i][j]=0;
		}
	}
	for(i=0;i<=NX;i++)
	{
		solid_flag[i][0]=1;
		rho[i][0]=0;
		solid_flag[i][NY]=1;
		rho[i][NY]=0;
	}
	for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			if(solid_flag[i][j]!=1)
			{
				if(i==0)
				{
					u[i][j][0]=0.8*j*(NY-j)/(NY*NY);   //网格点初始化x速度为0
				}
				else
				{
					u[i][j][0]=0;
				}
				u[i][j][1]=0;   //网格点初始化y速度为0
				rho[i][j]=rho0; //网格点密度初始化
				for(k=0;k<Q;k++)//计算每个网格点的初始分布函数
				{
					f[i][j][k]=feq(k,rho[i][j],u[i][j]);
				}
			}
		}
	}

}

double feq(int k,double rho,double u[2]) //计算平衡态分布函数
{
	double eu,uv,feq;
	eu= (e[k][0]*u[0]+e[k][1]*u[1]);    //求得e·u
	uv= (u[0]*u[0]+u[1]*u[1]);          //求得u的平方            关于速度u的矩阵处理，矩阵问题
	feq=w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);//Cs的平方=1/3
	return feq;
}

void  evolution()                       //演化方程
{
	for(i=0;i<=NX;i++)                    //碰撞过程
	{
		for(j=0;j<=NY;j++)
		{
			if(solid_flag[i][j]==0)
			{
				for(k=0;k<Q;k++)
				{
					f1[i][j][k]=f[i][j][k]+(feq(k,rho[i][j],u[i][j])-f[i][j][k])/tau_f;//+3*e[k][0]*P0*w[k];	
				}
			}
		}
	}
	///////////////////
	for(i=0;i<=NX;i++)                   //格子点的迁移过程
	{
		for(j=0;j<=NY;j++)
		{
			if(solid_flag[i][j]==0)
			{
				for(k=0;k<Q;k++)
				{
					ip=i-e[k][0];           
					jp=j-e[k][1];
					if(solid_flag[ip][jp]==1)
					{
						F[i][j][k]=f1[i][j][opp[k]];
					}
					if(solid_flag[ip][jp]==0)
					{
						F[i][j][k]=f1[ip][jp][k];
					}
				}
			}
		}
	}
////////////////////////////////
	for(i=1;i<NX;i++)                    ///区域宏观量计算
	{
		for(j=0;j<=NY;j++)
		{
			if(solid_flag[i][j]==0)
			{
				u0[i][j][0]=u[i][j][0];
				u0[i][j][1]=u[i][j][1];
				rho[i][j]=0;
				u[i][j][0]=0;
				u[i][j][1]=0;
				for(k=0;k<Q;k++)
				{
					f[i][j][k]=F[i][j][k];
					rho[i][j]+=f[i][j][k];  //宏观密度的控制方程
					u[i][j][0]+=e[k][0]*f[i][j][k];
					u[i][j][1]+=e[k][1]*f[i][j][k];
				}
				u[i][j][0]=u[i][j][0]/rho[i][j];//宏观速度的控制方程
				u[i][j][1]=u[i][j][1]/rho[i][j];
			}
		}
	}
	for(j=0;j<=NY;j++)                     
	{
		for(k=0;k<Q;k++)
		{   
			if(solid_flag[i][j]==0)
			{
				rho[0][j]=rho[1][j];
				f[0][j][k]=feq(k,rho[0][j],u[0][j])+f[1][j][k]-feq(k,rho[1][j],u[1][j]);//节点分布函数     
				u[NX][j][0]=u[NX-1][j][0];
				u[NX][j][1]=0;
				rho[NX][j]=rho[NX-1][j];
				f[NX][j][k]=feq(k,rho[NX][j],u[NX][j])+f[NX-1][j][k]-feq(k,rho[NX-1][j],u[NX-1][j]);//节点分布函数
			}
		}
	}
	
}

void output(int m)                      //输出函数
{
	ostringstream name;
	name<<"cavity_"<<m<<".dat";
	ofstream out(name.str().c_str());
	out<<"Title=\"LBM Lid Driven Flow\"\n"<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n"<<"ZONE T=\"BOX\",I="<<NX+1<<",J="<<NY+1<<",F=POINT"<<endl;
	for(j=0;j<=NY;j++)
	{
		for(i=0;i<=NX;i++)
		{
			out<<double(i)<<"	"
				<<double(j)<<"	"<<u[i][j][0]<<"	"<<u[i][j][1]<<endl;
		}
	}
}

void Error()
{
	double temp1,temp2;                  //为了进行比较准则的临时变量
	temp1=0;
	temp2=0;
	for(i=1;i<NX;i++)
	{
		for(j=1;j<NY;j++)
		{
			temp1+=(u[i][j][0]-u0[i][j][0])*(u[i][j][0]-u0[i][j][0])+(u[i][j][1]-u0[i][j][1])*(u[i][j][1]-u0[i][j][1]);
			temp2+=(u[i][j][0]*u[i][j][0])+(u[i][j][1]*u[i][j][1]);     
		}
	}
	temp1=sqrt(temp1);
	temp2=sqrt(temp2);
	error=temp1/(temp2+1e-30);
}

