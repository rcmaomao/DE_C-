
#include "stdafx.h"
#include <iostream>
#include <random>
#include <ctime>

using namespace std;
//种群最大个体数
constexpr int MAXPOP = 100;
//每个个体最大基因数
constexpr int MAXDIM = 15;

//用于存储初始种群
double c[MAXPOP][MAXDIM];
//用于存储老种群
double oldarray[MAXPOP][MAXDIM];
//用于存储新种群
double newarray[MAXPOP][MAXDIM];

//初始化函数1
void initial1(int& D,int& NP, long& nfeval, double inibound_h[], double inibound_l[],int n);
//声明评估函数1
double  evaluate(int D, double tmp[], long *nfeval,int n);


//完成数组的复制功能，将b中的内容复制到a中
void CopyVector(double a[], double b[])
{
	//对数组中的每一个元素
	for (int k = 0; k < MAXDIM; k++)
	{
		//b中的内容放到a中
		a[k] = b[k];
	}

}

//完成二维数组的复制功能
void CopyArray(double dest[MAXPOP][MAXDIM], double src[MAXPOP][MAXDIM])
{
	//对每一行
	for (int j = 0; j < MAXPOP; j++)
	{
		//对每一列上的元素
		for (int k = 0; k < MAXDIM; k++)
		{
			//复制
			dest[j][k] = src[j][k];
		}
	}

}

//主函数
int main()
{	
	int   i, j, L, n;      // 循环控制变量                 
	int   r1, r2, r3;	   // 声明用于父个体选择的下标       
	int   D;               // 个体的实际维数      
	int   NP;              // 种群数量      
	int   imin;            // 最小值的索引
	int   gen=0, genmax=300;
	long  nfeval=0;          // 执行评价函数的次数     
	double trial_energy;    // 临时变量，用于存储评价值                    
	double inibound_h[MAXDIM];      // 参数上界              
	double inibound_l[MAXDIM];      // 参数下界             
	double tmp[MAXDIM], best[MAXDIM], bestit[MAXDIM];   
	double energy[MAXPOP];  // 每个个体的评价值存储数组                 
	double F=0.8, CR=0.6;           // 缩放因子与杂交概率            
	double emin;            // 最小值  
	int strategy;
	cout << "请选择序号" << endl;
	cin >> strategy;
	initial1(D, NP, nfeval,inibound_h, inibound_l, strategy);
	double r;//声明double变量
	//设置随机数引擎和种子
	default_random_engine e(time(0));
	//default_random_engine e(3216586);
	//设置0到1之间均匀分布的实数随机数
	uniform_real_distribution<double> randReal(0, 1);
	//设置0到NP-1间均匀分布的整数随机数
	uniform_int_distribution<unsigned> randIntNP(0, NP-1);
	//设置0到D-1之间均匀分布的整数随机数
	uniform_int_distribution<unsigned> randIntD(0, D-1);
	
	//初始化种群，对种群内的每个个体
	for (i = 0; i < NP; i++) 
	{
		//初始化各种的每个维度
		for (j = 0; j < D; j++) 
		{
			//获得0到1的随机数
			r = randReal(e);
			//复制，初始化
			c[i][j] = inibound_l[j] + r * (inibound_h[j] - inibound_l[j]);
		}
		//对这个个体进行评估，并将结果保存
		energy[i] = evaluate(D, c[i], &nfeval, strategy);

	}
	//将最小值临时设置为评估值存储数组的第一个值
	emin = energy[0];
	//将最小值索引设置为0
	imin = 0;
	//遍历种群中的所有个体
	for (i = 1; i < NP; i++)
	{
		//如果这个个体的评估值小于当前最小值
		if (energy[i] < emin)
		{
			//更新当前最小值
			emin = energy[i];
			//更新最小值索引为当前索引
			imin = i;
		}
	}
	//将当前最好的个体保存
	CopyVector(best, c[imin]);
	//将当前最好的个体保存
	CopyVector(bestit, c[imin]);
	//将初始化种群赋为老种群
	CopyArray(oldarray, c);
	//迭代次数初始化为0
	gen = 0; 
	//当迭代次数小于设置的最大迭代次数时
	while ((gen < genmax)) 
	{
		//首先将当前迭代次数加一
		gen++;
		//将当前最优值索引设置为0
		imin = 0;
		//遍历种群中的个体
		for (i = 0; i < NP; i++) 
		{ 
			do {

				//获得0到NP间的一个整数
				r1 = randIntNP(e);
			} while (r1 == i);
			//若r1不等于i则终止

			do {
				//产生一个0到1之间的随机数
				r2 = randIntNP(e);
			} while ((r2 == i) || (r2 == r1));
			//若r2不等于r1且不等于i则终止

			do {
				//产生一个0到1之间的随机数
				r3 = randIntNP(e);
			} while ((r3 == i) || (r3 == r1) || (r3 == r2));
			//若r3不等于i且r3不等于r2和r1则终止

			//将当前个体数据保存到临时数组中
			CopyVector(tmp, oldarray[i]);
			//将n设置为0到D-1之间的一个随机整数
			n = randIntD(e);
			//对于这个个体的所有基因
			for (L = 0; L < D; L++) 
			{
				//当产生的0到1的一个随机数小于杂交概率CR或者当前执行循环次数为基因总数时
				if ((randReal(e) < CR) || L == (D - 1)) 
				{
					//更新这个个体的第n+1个基因
					tmp[n] = oldarray[r1][n] + F * (oldarray[r2][n] - oldarray[r3][n]);
					//如果这个基因值超出了边界
					if (tmp[n] > inibound_h[n] || tmp[n] < inibound_l[n])
					{
						//则随机产生新的边界内的基因值
						tmp[n]= inibound_l[n] + r * (inibound_h[n] - inibound_l[n]);
					}
				}
				//修改索引值，从而修改不同位置的基因
				n = (n + 1) % D;
			}
			// 对当前变异结果评估
			trial_energy = evaluate(D, tmp, &nfeval, strategy);
			//如果这个个体优于父代中对应位置的个体
			if (trial_energy <= energy[i]) 
			{
				//修改评估值存储数组中的对应位置个体的评价值
				energy[i] = trial_energy;
				//将这个个体添加到新种群中
				CopyVector(newarray[i], tmp);
				// 如果这个个体是目前为止最优的
				if (trial_energy < emin)
				{
					//更新最小值
					emin = trial_energy;
					//更新最小值索引
					imin = i;
					//将这个个体存储起来
					CopyVector(best, tmp);
				}
			}
			//如果这个个体不优于老种群中对应位置的个体
			else 
			{
				//将老种群中这个位置的个体加入到新种群中
				CopyVector(newarray[i], oldarray[i]);
			}
		}
		//将经过一次遍历之后最优的个体存储
		CopyVector(bestit, best);
		//用新种群替换旧种群
		CopyArray(oldarray, newarray);

	}

	//输出最佳结果
	cout << "目前为止的最佳值是：" << emin << endl;
	//对最优个体的每个基因
	for (j = 0; j < D; j++) 
	{
		//输出这个基因值
		cout << "best[" <<j << "]:"<<best[j] << endl;
	}
	//输出迭代次数
	cout << "迭代次数是：" << gen << endl;
	//输出评估次数
	cout << "评估次数是：" << nfeval << endl;
	//输出所使用参数值
	cout << " NP: " << NP << "  F: " << F << "  CR: " << CR << endl;

	return(0);
}

//实现评估函数
double evaluate(int D, double tmp[], long *nfeval,int n)
{
	//评估函数执行次数增加
	(*nfeval)++;
	//f(x1,x2)=5cos(x1x2)+4x1+x2
	//return 5 * cos(tmp[0] * tmp[1]) + 4 * tmp[0] + tmp[1];
	if (n == 1)
	{
		double result = 0;
		for (int i = 0; i < 10; i++)
		{
			result += tmp[i] * tmp[i] - 10 * cos(2 * 3.1415926*tmp[i]);
		}
		return result + 100;
	}
	if (n == 2)
	{
		return -20 * exp(-0.2*sqrt(0.5*(tmp[0] * tmp[0] + tmp[1] * tmp[1]))) - exp(0.5*(cos(2 * 3.1415926*tmp[0])
			+ cos(2 * 3.1415926*tmp[1]))) + 20 + exp(1);
	}
	if (n == 3)
	{
		double result = 0;
		for (int i = 0; i < 10; i++)
		{
			result += tmp[i] * tmp[i];
		}
		return result;
	}
	if (n==4)
	{
		double result = 0;
		for (int i = 0; i < 9; i++)
		{
			result += 100 * (tmp[i + 1] - tmp[i] * tmp[1])*(tmp[i + 1] - tmp[i] * tmp[1])
				+ (1 - tmp[i])*(1 - tmp[i]);
		}
		return result;
	}
	if (n == 5)
	{
		double x = tmp[0];
		double y = tmp[1];
		double result = 0;
		result += (1.5 - x - x * y)*(1.5 - x - x * y) + (2.25 - x + x * y*y)*(2.25 - x + x * y*y);
		result += (2.625 - x + x * y*y*y)*(2.625 - x + x * y*y*y);
		return result;
	}
	if (n == 6)
	{
		
		double x = tmp[0];
		double y = tmp[1];
		double div = cos(abs(x*x - y * y))*cos(abs(x*x - y * y))-0.5;
		double dd = (1 + 0.001*(x*x + y * y))*(1 + 0.001*(x*x + y * y));
		return 0.5 + div / dd;
	}
	if (n==7)
	{
		double result = 0;
		for (int i = 0; i < D; i++)
		{
			result += pow(tmp[i], 4) - 16 * tmp[i] * tmp[i] + 5 * tmp[i];
		}
		return result / 2;
	}
	if (n == 8)
	{
		return 100 * sqrt(abs(tmp[1] - 0.01*tmp[0] * tmp[0])) + 0.01*abs(tmp[0] + 10);
	}
	if (n == 9)
	{
		double x = tmp[0];
		double y = tmp[1];
		double result = (x*x + y - 11)*(x*x + y - 11);
		result += (x + y * y - 7)*(x + y * y - 7);
		return result;
	}
	if (n == 10)
	{
		double x = tmp[0];
		double y = tmp[1];
		double exponent = abs(100 - (sqrt(x*x + y * y) / 3.1415926))+1;
		double base = abs(sin(x)*sin(y)*exp(exponent));
		return -0.0001*pow(base, 0.1);
	}

}

void initial1(int& D, int& NP, long& nfeval, double inibound_h[], double inibound_l[], int n)
{
	nfeval = 0;
	if (n == 1)
	{
		D = 10;
		NP = 30;
		for (int i = 0; i < 10; i++)
		{
			inibound_h[i] = 5.12;
			inibound_l[i] = -5.12;
		}
	}
	else if (n==2)
	{
		D = 2;
		NP = 30;
		for (int i = 0; i < D; i++)
		{
			inibound_h[i] = 5;
			inibound_l[i] = -5;
		}
	}
	else if (n==3)
	{
		D = 10;
		NP = 40;
		for (int i = 0; i < D; i++)
		{
			inibound_h[i] = 100;
			inibound_l[i] = -100;
		}
	}
	else if (n==4)
	{
		D = 10;
		NP = 30;
		for (int i = 0; i < D; i++)
		{
			inibound_h[i] = 30;
			inibound_l[i] = -30;
		}
	}
	else if (n==5)
	{
		D = 2;
		NP = 20;
		for (int i = 0; i < D; i++)
		{
			inibound_h[i] = 4.5;
			inibound_l[i] = -4.5;
		}
	}
	else if (n==6)
	{
		D = 2;
		NP = 50;
		for (int i = 0; i < D; i++)
		{
			inibound_h[i] = 100;
			inibound_l[i] = -100;
		}
	}
	else if (n==7)
	{
		D = 10;
		NP = 30;
		for (int i = 0; i < D; i++)
		{
			inibound_h[i] = 5;
			inibound_l[i] = -5;
		}
	}
	else if (n==8)
	{
		D = 2;
		NP = 25;
		inibound_h[0] = -5;
		inibound_l[0] = -15;
		inibound_h[1] = 3;
		inibound_l[1] = -3;
	}
	else if (n==9)
	{
		D = 2;
		NP = 30;
		for (int i = 0; i < D; i++)
		{
			inibound_h[i] = 5;
			inibound_l[i] = -5;
		}
	}
	else if (n==10)
	{
		D = 2;
		NP = 30;
		for (int i = 0; i < D; i++)
		{
			inibound_h[i] = 10;
			inibound_l[i] = -10;
		}
	}

}
