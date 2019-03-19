
#include "stdafx.h"
#include <iostream>
#include <random>
#include <ctime>

using namespace std;
//��Ⱥ��������
constexpr int MAXPOP = 100;
//ÿ��������������
constexpr int MAXDIM = 15;

//���ڴ洢��ʼ��Ⱥ
double c[MAXPOP][MAXDIM];
//���ڴ洢����Ⱥ
double oldarray[MAXPOP][MAXDIM];
//���ڴ洢����Ⱥ
double newarray[MAXPOP][MAXDIM];

//��ʼ������1
void initial1(int& D,int& NP, long& nfeval, double inibound_h[], double inibound_l[],int n);
//������������1
double  evaluate(int D, double tmp[], long *nfeval,int n);


//�������ĸ��ƹ��ܣ���b�е����ݸ��Ƶ�a��
void CopyVector(double a[], double b[])
{
	//�������е�ÿһ��Ԫ��
	for (int k = 0; k < MAXDIM; k++)
	{
		//b�е����ݷŵ�a��
		a[k] = b[k];
	}

}

//��ɶ�ά����ĸ��ƹ���
void CopyArray(double dest[MAXPOP][MAXDIM], double src[MAXPOP][MAXDIM])
{
	//��ÿһ��
	for (int j = 0; j < MAXPOP; j++)
	{
		//��ÿһ���ϵ�Ԫ��
		for (int k = 0; k < MAXDIM; k++)
		{
			//����
			dest[j][k] = src[j][k];
		}
	}

}

//������
int main()
{	
	int   i, j, L, n;      // ѭ�����Ʊ���                 
	int   r1, r2, r3;	   // �������ڸ�����ѡ����±�       
	int   D;               // �����ʵ��ά��      
	int   NP;              // ��Ⱥ����      
	int   imin;            // ��Сֵ������
	int   gen=0, genmax=300;
	long  nfeval=0;          // ִ�����ۺ����Ĵ���     
	double trial_energy;    // ��ʱ���������ڴ洢����ֵ                    
	double inibound_h[MAXDIM];      // �����Ͻ�              
	double inibound_l[MAXDIM];      // �����½�             
	double tmp[MAXDIM], best[MAXDIM], bestit[MAXDIM];   
	double energy[MAXPOP];  // ÿ�����������ֵ�洢����                 
	double F=0.8, CR=0.6;           // �����������ӽ�����            
	double emin;            // ��Сֵ  
	int strategy;
	cout << "��ѡ�����" << endl;
	cin >> strategy;
	initial1(D, NP, nfeval,inibound_h, inibound_l, strategy);
	double r;//����double����
	//������������������
	default_random_engine e(time(0));
	//default_random_engine e(3216586);
	//����0��1֮����ȷֲ���ʵ�������
	uniform_real_distribution<double> randReal(0, 1);
	//����0��NP-1����ȷֲ������������
	uniform_int_distribution<unsigned> randIntNP(0, NP-1);
	//����0��D-1֮����ȷֲ������������
	uniform_int_distribution<unsigned> randIntD(0, D-1);
	
	//��ʼ����Ⱥ������Ⱥ�ڵ�ÿ������
	for (i = 0; i < NP; i++) 
	{
		//��ʼ�����ֵ�ÿ��ά��
		for (j = 0; j < D; j++) 
		{
			//���0��1�������
			r = randReal(e);
			//���ƣ���ʼ��
			c[i][j] = inibound_l[j] + r * (inibound_h[j] - inibound_l[j]);
		}
		//�����������������������������
		energy[i] = evaluate(D, c[i], &nfeval, strategy);

	}
	//����Сֵ��ʱ����Ϊ����ֵ�洢����ĵ�һ��ֵ
	emin = energy[0];
	//����Сֵ��������Ϊ0
	imin = 0;
	//������Ⱥ�е����и���
	for (i = 1; i < NP; i++)
	{
		//���������������ֵС�ڵ�ǰ��Сֵ
		if (energy[i] < emin)
		{
			//���µ�ǰ��Сֵ
			emin = energy[i];
			//������Сֵ����Ϊ��ǰ����
			imin = i;
		}
	}
	//����ǰ��õĸ��屣��
	CopyVector(best, c[imin]);
	//����ǰ��õĸ��屣��
	CopyVector(bestit, c[imin]);
	//����ʼ����Ⱥ��Ϊ����Ⱥ
	CopyArray(oldarray, c);
	//����������ʼ��Ϊ0
	gen = 0; 
	//����������С�����õ�����������ʱ
	while ((gen < genmax)) 
	{
		//���Ƚ���ǰ����������һ
		gen++;
		//����ǰ����ֵ��������Ϊ0
		imin = 0;
		//������Ⱥ�еĸ���
		for (i = 0; i < NP; i++) 
		{ 
			do {

				//���0��NP���һ������
				r1 = randIntNP(e);
			} while (r1 == i);
			//��r1������i����ֹ

			do {
				//����һ��0��1֮��������
				r2 = randIntNP(e);
			} while ((r2 == i) || (r2 == r1));
			//��r2������r1�Ҳ�����i����ֹ

			do {
				//����һ��0��1֮��������
				r3 = randIntNP(e);
			} while ((r3 == i) || (r3 == r1) || (r3 == r2));
			//��r3������i��r3������r2��r1����ֹ

			//����ǰ�������ݱ��浽��ʱ������
			CopyVector(tmp, oldarray[i]);
			//��n����Ϊ0��D-1֮���һ���������
			n = randIntD(e);
			//���������������л���
			for (L = 0; L < D; L++) 
			{
				//��������0��1��һ�������С���ӽ�����CR���ߵ�ǰִ��ѭ������Ϊ��������ʱ
				if ((randReal(e) < CR) || L == (D - 1)) 
				{
					//�����������ĵ�n+1������
					tmp[n] = oldarray[r1][n] + F * (oldarray[r2][n] - oldarray[r3][n]);
					//����������ֵ�����˱߽�
					if (tmp[n] > inibound_h[n] || tmp[n] < inibound_l[n])
					{
						//����������µı߽��ڵĻ���ֵ
						tmp[n]= inibound_l[n] + r * (inibound_h[n] - inibound_l[n]);
					}
				}
				//�޸�����ֵ���Ӷ��޸Ĳ�ͬλ�õĻ���
				n = (n + 1) % D;
			}
			// �Ե�ǰ����������
			trial_energy = evaluate(D, tmp, &nfeval, strategy);
			//�������������ڸ����ж�Ӧλ�õĸ���
			if (trial_energy <= energy[i]) 
			{
				//�޸�����ֵ�洢�����еĶ�Ӧλ�ø��������ֵ
				energy[i] = trial_energy;
				//�����������ӵ�����Ⱥ��
				CopyVector(newarray[i], tmp);
				// ������������ĿǰΪֹ���ŵ�
				if (trial_energy < emin)
				{
					//������Сֵ
					emin = trial_energy;
					//������Сֵ����
					imin = i;
					//���������洢����
					CopyVector(best, tmp);
				}
			}
			//���������岻��������Ⱥ�ж�Ӧλ�õĸ���
			else 
			{
				//������Ⱥ�����λ�õĸ�����뵽����Ⱥ��
				CopyVector(newarray[i], oldarray[i]);
			}
		}
		//������һ�α���֮�����ŵĸ���洢
		CopyVector(bestit, best);
		//������Ⱥ�滻����Ⱥ
		CopyArray(oldarray, newarray);

	}

	//�����ѽ��
	cout << "ĿǰΪֹ�����ֵ�ǣ�" << emin << endl;
	//�����Ÿ����ÿ������
	for (j = 0; j < D; j++) 
	{
		//����������ֵ
		cout << "best[" <<j << "]:"<<best[j] << endl;
	}
	//�����������
	cout << "���������ǣ�" << gen << endl;
	//�����������
	cout << "���������ǣ�" << nfeval << endl;
	//�����ʹ�ò���ֵ
	cout << " NP: " << NP << "  F: " << F << "  CR: " << CR << endl;

	return(0);
}

//ʵ����������
double evaluate(int D, double tmp[], long *nfeval,int n)
{
	//��������ִ�д�������
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
