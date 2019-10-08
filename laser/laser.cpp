// laser.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include <math.h>

#define PI 3.14159265358979323846


const double B_table[6][10] = 
{
	{0.99528393, 0.01311258, -0.00148875, 0.00137257, -0.0121624, 0.00000038, 0.00000026, -0.00042213, 0, 0.00376331},
	{0.99300868, -0.01959749, -0.00227073, 0.00212873, -0.0183182, 0.00000035, 0.00000067, 0.0006653, 0, 0.00570955},
	{1.65651297, -1.92604145, -0.00646772, -0.00008351, 1.88594691, -0.00023356, 0.00014466, 0.00247856, 0.00000045, -0.61624228},
	{1.00426951, -0.01300905, 0.00361072, 0.00342946, 0.01322436, -0.00000698, 0.00000654, -0.00107642, 0, -0.00448256},
	{1.31523023, -0.92958114, -0.00727529, 0.00867352, 0.90633602, -0.00011192, 0.00008907, -0.0030339, 0.00000006, -0.29389689},
	{0.46726199, 1.80965096, -0.21839542, 0.34811215, -1.6162656, -0.00197022, 0.00074627, -0.13346963, 0.00001856, 0.03007356}
};
const double P_0 = 101325;
//0��ľ����¶�
const double T_0 = 273.15;
//t=0��ʱ�ı���ˮ��ѹ
const double E_0 = 611;
//�������ӵ�ƽ���뾶
const double r = 1.8E-4;
const double lambda_0 = 0.55;


double get_a_m(double lambda, double T_1, double I)
{
	double B_lambda[10] = {};
	if (0.83 <= lambda)
	{
		memcpy(B_lambda, B_table[0], sizeof(B_lambda));
	} 
	else if (1.06 <= lambda)
	{
		memcpy(B_lambda, B_table[1], sizeof(B_lambda));
	}
	else if (1.33 <= lambda)
	{
		memcpy(B_lambda, B_table[2], sizeof(B_lambda));
	}
	else if (1.54 <= lambda)
	{
		memcpy(B_lambda, B_table[3], sizeof(B_lambda));
	}
	else if (3.82 <= lambda)
	{
		memcpy(B_lambda, B_table[4], sizeof(B_lambda));
	}
	else
	{
		memcpy(B_lambda, B_table[5], sizeof(B_lambda));
	}

	return -log(B_lambda[0] + B_lambda[1]*T_1 + B_lambda[2]*I + B_lambda[3]*T_1*I + \
		B_lambda[4]*pow(T_1,2) + B_lambda[5]*pow(I,2) + B_lambda[6]*T_1*pow(I,2) + \
		B_lambda[7]*pow(T_1,2)*I + B_lambda[8]*pow(I,3) + B_lambda[1]*pow(T_1,3));
}

double get_a_s(double lambda, double N_g)
{
	//���ӵ�ɢ�����
	double A = 4 * PI * pow(r,2);
	return 0.827 * (N_g * pow(A,3) / pow(lambda,4)) * (1E-3);
}

double get_k_ms(double lambda, double V, double q)
{
	return (3.912 / V) * pow((lambda_0 / lambda), q);
}

double get_V(double V_0, double H, double b)
{
	return V_0 * exp(b * H);
}

double get_q(double V)
{
	if (V > 50)
	{
		return 1.6;
	}
	else if (V > 6 && V <= 50)
	{
		return 1.3;
	}
	else if (V > 1 && V <= 6)
	{
		return 0.16 * V + 0.34;
	}
	else if (V > 0.5 && V <= 1)
	{
		return V - 0.5;
	}
	return 0;
}

double get_t(double t_0, double H)
{
	return t_0 - 5*H;
}

double get_T(double t)
{
	return t + T_0;
}

double get_T_1(double H, double t_0)
{
	double t = t_0 - 5*H;
	double T = t + 273.15;
	double T_1 = T / T_0;
	return T_1;
}

double get_P(double H)
{
	return P_0 * (1 - 0.11975*H + 0.0053788*pow(H,2));
}

// double get_rho(double P, double T_1)
// {
// 	const double rho_0 = 1293;
// 	return P_0 * (P / P_0) * (1 / T_1);
// }

double get_N_g(double P, double T)
{
	return (6.02E23 / 22.4E3) * ((P*T_0) / (P_0*T));
}

double get_e(double f, double E)
{
	return f * E;
}

double get_E(double t, double T)
{
	return E_0 * exp(19.9*t / T);
}

double get_I(double e, double T)
{
	return 2.17*e / T;
}

double get_d(double lambda, double D_0, double R)
{
	//sin(theta) = (1.22*lambda) / D_0 * 1E-6
	//d = theta * R
	return (1.22*lambda) / D_0 * R * 1E-1;
}

double get_P_r(double P_t, double beta, double R, double d)
{
	return (4 * P_t * exp(-beta * R)) / (PI * pow(d,2));
}

int _tmain(int argc, _TCHAR* argv[])
{
	// input_0
	double b = 0.1;//��ɳ���������µ���1.25��������������ԼΪ0.1
	// input_1
	//0���θ߶ȴ����ܼ���
	double V_0 = 10.0;
	//0���θ߶ȵ�����
	double t_0 = 10.0;
	//�������ʪ��
	double f = 80.0;
	// input_2
	//���ⷢ�书��
	double P_t = 100000.0;
	//����������ھ�
	double D_0 = 0.5;
	//���θ߶�
	double H = 10.0;
	//���Ⲩ��
	double lambda = 1.06;
	// input_3
	//���⴫�����
	double R = 5.0;

	//�����ܼ���
	double V = get_V(V_0, H, b);
	//ϵ��
	double q = get_q(V);

	//���϶�
	double t = get_t(t_0, H);
	//������
	double T = get_T(t);
	//�¶���t��ı���ˮ��ѹ
	double E = get_E(t, T);
	//ˮ��ѹ
	double e = get_e(f, E);

	//����ѹǿ
	double P = get_P(H);
	//�����¶�
	double T_1 = get_T_1(H, t_0);
	//�������ӵ�λ�������
	double N_g = get_N_g(P, T);
	//��������ʪ��
	double I = get_I(e, T);

	//������������ϵ��
	double a_m = get_a_m(lambda, T_1, I);
	//��������ɢ��ϵ��
	double a_s = get_a_s(lambda, N_g);
	//���ܽ����պ�ɢ��ϵ��
	double k_ms = get_k_ms(lambda, V ,q);
	//����˥��ϵ��
	double beta = a_m + a_s + k_ms;

	//������ֱ��
	double d = get_d(lambda, D_0, R);

	//�������䵽Ŀ���ϵĹ����ܶ�
	double P_r = get_P_r(P_t, beta, R, d);

	return 0;
}

