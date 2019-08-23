// DLT.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <cstring>

using namespace std;

void TransposeMatrix(double *m1, double *m2, int row, int col); //矩阵转置
void AddMatrix(double *m1, double *m2, int row, int col); //矩阵相加
void SubMatrix(double *m1, double *m2, int row, int col); //矩阵相减
void InverseMatrix(double *a, int n); //矩阵求逆
void MultiplyMatrix(double *m1, double *m2, double *result, int m, int n, int l); //矩阵相乘，其中result指向输出矩阵
struct Coordinate
{
	int num;
	double x;
	double y;
	double X;
	double Y;
	double Z;

};

struct ImageCoordinate
{
	int num;
	double x;
	double y;

};

struct ObjectCoordinate
{
	int num;
	double x;
	double y;
	double z;

};

struct point
{
	int num;
	//两张像片上的像点坐标
	double x1;
	double y1;
	double x2;
	double y2;
	//控制点坐标
	double X;
	double Y;
	double Z;
	//计算出的控制点坐标
	double X_;
	double Y_;
	double Z_;

};

int ReadFile(char *s1, char *s2, Coordinate *pk);
double  *DLT_(char *s1, char *s2);
int main()
{
	vector<double> vec; //像点坐标保存向量
	vector<double> vec1; //物方空间点坐标保存向量
	ifstream infile;
	infile.open("912.txt");
	if (!infile.is_open()) cout << "can't open this file!" << endl;

	double temp;
	while (infile >> temp) {
		vec.push_back(temp);
	}

	infile.close(); //将像点坐标保存在vec中
	ifstream infile1;
	infile1.open("wu91.txt");
	if (!infile1.is_open()) cout << "can't open this file" << endl;

	double temp1;
	while (infile1 >> temp1) {
		vec1.push_back(temp1); //将物方空间点坐标保存在vec1中
	}
	infile1.close();


	int num = vec1.size() / 4; //参与计算的点的个数
	Coordinate *c = new Coordinate[num];
	for (int i = 0; i < num; i++) {
		c[i].num = i;
		c[i].x = vec[3 * i + 1] * 0.0051966;
		c[i].y = vec[3 * i + 2] * 0.0051966;
		c[i].X = vec1[4 * i + 1];
		c[i].Y = vec1[4 * i + 2];
		c[i].Z = vec1[4 * i + 3];
	}
	//初始化li系数
	double l1 = 0, l2 = 0, l3 = 0, l4 = 0, l5 = 0, l6 = 0, l7 = 0, l8 = 0, l9 = 0, l10 = 0, l11 = 0;

	//选取的5.5个点位的像方坐标和物方坐标
	double x1 = vec[1] * 0.0051966, y1 = vec[2] * 0.0051966, x2 = vec[4] * 0.0051966, y2 = vec[5] * 0.0051966,
		x3 = vec[7] * 0.0051966, y3 = vec[8] * 0.0051966, x4 = vec[10] * 0.0051966, y4 = vec[11] * 0.0051966,
		x5 = vec[13] * 0.0051966, y5 = vec[14] * 0.0051966, x6 = vec[16] * 0.0051966;
	/*cout << x1 << endl << y1 << endl << x2 << endl << y2 << endl << x3 << endl << y3 << endl << x4 << endl << y4 <<
		endl << x5 << endl << y5 << endl << x6 << endl;
	system("pause");*/
	
	double Z1 = -1 * vec1[1], X1 = vec1[2], Y1 = vec1[3],
		Z2 = -1 * vec1[5], X2 = vec1[6], Y2 = vec1[7],
		Z3 = -1 * vec1[9], X3 = vec1[10], Y3 = vec1[11],
		Z4 = -1 * vec1[13], X4 = vec1[14], Y4 = vec1[15],
		Z5 = -1 * vec1[17], X5 = vec1[18], Y5 = vec1[19],
		Z6 = -1 * vec1[21], X6 = vec1[22], Y6 = vec1[23];
	/*cout << X1 << endl << Y1 << endl << Z1 << endl
		<< X2 << endl << Y2 << endl << Z2 << endl
		<< X3 << endl << Y3 << endl << Z3 << endl
		<< X4 << endl << Y4 << endl << Z4 << endl
		<< X5 << endl << Y5 << endl << Z5 << endl
		<< X6 << endl << Y6 << endl << Z6 << endl;
	system("pause");*/

	//初始化设定
	double ATL[11] = { 0.0 };
	double A[11 * 12] = { 0.0 };
	double AT[11 * 12] = { 0.0 };
	double ATA[11 * 11] = { 0.0 };

	double L[11] = { 0.0 };
	double XX[15] = { 0.0 };
	double Qii[15] = { 0.0 };

	double x0, y0;
	//对A赋值
	A[0] = X1; A[1] = Y1; A[2] = Z1; A[3] = 1; A[4] = 0; A[5] = 0; 
	A[6] = 0; A[7] = 0; A[8] = x1 * X1; A[9] = x1 * Y1; A[10] = x1 * Z1;

	A[11] = 0; A[12] = 0; A[13] = 0; A[14] = 0; A[15] = X1; A[16] = Y1; 
	A[17] = Z1; A[18] = 1; A[19] = y1 * X1; A[20] = y1 * Y1; A[21] = y1 * Z1;

	A[22] = X2; A[23] = Y2; A[24] = Z2; A[25] = 1; A[26] = 0; A[27] = 0;
	A[28] = 0; A[29] = 0; A[30] = x2 * X2; A[31] = x2 * Y2; A[32] = x2 * Z2;

	A[33] = 0; A[34] = 0; A[35] = 0; A[36] = 0; A[37] = X2; A[38] = Y2; 
	A[39] = Z2; A[40] = 1; A[41] = y2 * X2; A[42] = y2 * Y2; A[43] = y2 * Z2;

	A[44] = X3; A[45] = Y3; A[46] = Z3; A[47] = 1; A[48] = 0; A[49] = 0;
	A[50] = 0; A[51] = 0; A[52] = x3 * X3; A[53] = x3 * Y3; A[54] = x3 * Z3;

	A[55] = 0; A[56] = 0; A[57] = 0; A[58] = 0; A[59] = X3; A[60] = Y3;
	A[61] = Z3; A[62] = 1; A[63] = y3 * X3; A[64] = y3 * Y3; A[65] = y3 * Z3;

	A[66] = X4; A[67] = Y4; A[68] = Z4; A[69] = 1; A[70] = 0; A[71] = 0;
	A[72] = 0; A[73] = 0; A[74] = x4 * X4; A[75] = x4 * Y4; A[76] = x4 * Z4;

	A[77] = 0; A[78] = 0; A[79] = 0; A[80] = 0; A[81] = X4; A[82] = Y4; 
	A[83] = Z4; A[84] = 1; A[85] = y4 * X4; A[86] = y4 * Y4; A[87] = y4 * Z4;

	A[88] = X5; A[89] = Y5; A[90] = Z5; A[91] = 1; A[92] = 0; A[93] = 0;
	A[94] = 0; A[95] = 0; A[96] = x5 * X5; A[97] = x5 * Y5; A[98] = x5 * Z5;

	A[99] = 0; A[100] = 0; A[101] = 0; A[102] = 0; A[103] = X5; A[104] = Y5;
	A[105] = Z5; A[106] = 1; A[107] = y5 * X5; A[108] = y5 * Y5; A[109] = y5 * Z5;

	A[110] = X6; A[111] = Y6; A[112] = Z6; A[113] = 1; A[114] = 0; A[115] = 0;
	A[116] = 0; A[117] = 0; A[118] = x6 * X6; A[119] = x6 * Y6; A[120] = x6 * Z6;

	//L赋值
	L[0] = x1; L[1] = y1; L[2] = x2; L[3] = y2;
	L[4] = x3; L[5] = y3; L[6] = x4; L[7] = y4;
	L[8] = x5; L[9] = y5; L[10] = x6;
	//li系数近似值解算
	//double A1[11 * 12] = { 0.0 };
	double L1[11] = { 0.0 };
	for (int i = 0; i < 11; i++) {
		L1[i] = -1 * L[i];
	}
	TransposeMatrix(A, AT, 12, 11);
	MultiplyMatrix(AT, A, ATA, 11, 12, 11);
	MultiplyMatrix(AT, L, ATL, 11, 12, 1);
	InverseMatrix(ATA, 11);
	MultiplyMatrix(ATA, ATL, XX, 11, 11, 1);

	//内方位元素x0 y0解算
	x0 = -(XX[0] * XX[8] + XX[1] * XX[9] + XX[2] * XX[10]) / (XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10]);
	y0 = -(XX[4] * XX[8] + XX[5] * XX[9] + XX[6] * XX[10]) / (XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10]);

	//计算fx,fy
	double fx, fy;  //x向主距，y向主距
	double tem1 = 0, tem2 = 0, tem3 = 0;
	//double *c = new double[num];

	tem1 = XX[0] * XX[0] + XX[1] * XX[1] + XX[2] * XX[2];
	tem2 = XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10];
	tem3 = XX[4] * XX[4] + XX[5] * XX[5] + XX[6] * XX[6];

	fx = sqrt(tem1 / tem2 - x0*x0);

	double r3 = 0, ds = 0;
	r3 = 1 / (XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10]);
	ds = sqrt((r3*r3*tem1 - x0*x0) / (r3*r3*tem3 - y0*y0)) - 1;
	fy = fx / (1 + ds);

	double fx1, fy1;
	double *M = new double[15 * 2 * num];
	double *MT = new double[15 * 2 * num];
	double *W = new double[2 * num];
	double MTM[225] = { 0.0 };
	double MTW[15] = { 0.0 };
	double aa = 0, bb = 0, cc = 0;

	//解算li系数精确值
	do {
		fx1 = fx;
		fy1 = fy;
		for (int i = 0; i < num; i++) {
			double temp = XX[8] * c[i].X + XX[9] * c[i].Y + XX[10] * c[i].Z + 1;
			double r = sqrt((c[i].x - x0)*(c[i].x - x0) + (c[i].y - y0)*(c[i].y - y0));
			double r2 = pow(r, 2);//r的平方
			double r4 = pow(r, 4);//r的四次方
			M[30 * i + 0] = -c[i].X / temp;
			M[30 * i + 1] = -c[i].Y / temp;
			M[30 * i + 2] = -c[i].Z / temp;
			M[30 * i + 3] = -1 / temp;
			M[30 * i + 4] = 0;
			M[30 * i + 5] = 0;
			M[30 * i + 6] = 0;
			M[30 * i + 7] = 0;
			M[30 * i + 8] = -(c[i].x*c[i].X) / temp;
			M[30 * i + 9] = -(c[i].x*c[i].Y) / temp;
			M[30 * i + 10] = -(c[i].x*c[i].Z) / temp;
			M[30 * i + 11] = -(c[i].x - x0)*r2;
			M[30 * i + 12] = -(c[i].x - x0)*r4;
			M[30 * i + 13] = -r2 - 2 * (c[i].x - x0)*(c[i].x - x0);
			M[30 * i + 14] = -2 * (c[i].x - x0)*(c[i].y - y0);

			M[30 * i + 15] = 0;
			M[30 * i + 16] = 0;
			M[30 * i + 17] = 0;
			M[30 * i + 18] = 0;
			M[30 * i + 19] = -c[i].X / temp;
			M[30 * i + 20] = -c[i].Y / temp;
			M[30 * i + 21] = -c[i].Z / temp;
			M[30 * i + 22] = -1 / temp;
			M[30 * i + 23] = -(c[i].y*c[i].X) / temp;
			M[30 * i + 24] = -(c[i].y*c[i].Y) / temp;
			M[30 * i + 25] = -(c[i].y*c[i].Z) / temp;
			M[30 * i + 26] = -(c[i].y - y0)*r2;
			M[30 * i + 27] = -(c[i].y - y0)*r4;
			M[30 * i + 28] = -2 * (c[i].x - x0)*(c[i].y - y0);
			M[30 * i + 29] = -r2 - 2 * (c[i].y - y0)*(c[i].y - y0);

			W[2 * i + 0] = c[i].x / temp;
			W[2 * i + 1] = c[i].y / temp;
		}
		TransposeMatrix(M, MT, 2 * num, 15);
		MultiplyMatrix(MT, M, MTM, 15, 2 * num, 15);
		MultiplyMatrix(MT, W, MTW, 15, 2 * num, 1);
		InverseMatrix(MTM, 15);
		MultiplyMatrix(MTM, MTW, XX, 15, 15, 1);

		r3 = 1 / (XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10]);
		x0 = -r3*(XX[0] * XX[8] + XX[1] * XX[9] + XX[2] * XX[10]);
		y0 = -r3*(XX[4] * XX[8] + XX[5] * XX[9] + XX[6] * XX[10]);
		aa = r3*(XX[0] * XX[0] + XX[1] * XX[1] + XX[2] * XX[2]) - x0*x0;
		bb = r3*(XX[4] * XX[4] + XX[5] * XX[5] + XX[6] * XX[6]) - y0*y0;
		cc = r3*(XX[0] * XX[4] + XX[1] * XX[5] + XX[2] * XX[6]) - x0*y0;

		fx = sqrt((aa*bb - cc*cc) / bb);
		fy = sqrt((aa*bb - cc*cc) / aa);
	} while (fabs(fx - fx1)>0.01 || fabs(fy - fy1)>0.01);

	//计算单位权中误差
	double *Vii = new double[2 * num];
	MultiplyMatrix(M, XX, Vii, 2 * num, 15, 1);
	SubMatrix(Vii, W, 2 * num, 1);

	double v = 0;
	for (int i = 0; i < 2 * num; i++) {
		v += Vii[i] * Vii[i];
	}
	double m0 = sqrt(v / (2 * num - 15));

	for (int i = 0; i < 15; i++) {
		Qii[i] = MTM[15 * i + i];
	}
	double m[15] = { 0.0 };
	for (int i = 0; i < 15; i++) {
		m[i] = m0 * sqrt(Qii[i]);
	}

	x0 = -r3*(XX[0] * XX[8] + XX[1] * XX[9] + XX[2] * XX[10]);
	y0 = -r3*(XX[4] * XX[8] + XX[5] * XX[9] + XX[6] * XX[10]);
	//11个线性变换系数
	double l[11] = { 0 };
	for (int i = 0; i < 11; i++) {
		l[i] = XX[i];
	}
		//四个畸变系数
	double k1 = XX[11];
	double k2 = XX[12];
	double p1 = XX[13];
	double p2 = XX[14];

		//外方位线元素
	double wai[9] = { 0.0 };
	wai[0] = XX[0]; wai[1] = XX[1]; wai[2] = XX[2];
	wai[4] = XX[3]; wai[1] = XX[4]; wai[5] = XX[6];
	wai[8] = XX[6]; wai[1] = XX[7]; wai[8] = XX[10];

	double waiL[3] = { 0.0 };
	waiL[0] = -1 * XX[3]; waiL[1] = -1 * XX[7]; waiL[2] = -1;

	double w[3] = { 0 };
	InverseMatrix(wai, 3);
	MultiplyMatrix(wai, waiL, w, 3, 3, 1);

	//化为地面摄影测量坐标
	double Xs = w[0], Ys = w[1], Zs = w[2];
	double temz = Zs;
	Zs = Ys;
	Ys = Xs;
	Xs = -1 * temz;

	ds = sqrt(aa / bb) - 1;
	double dp = asin(sqrt(cc*cc / aa / bb));

	double LL[3 * 4] = { XX[0],XX[1],XX[2],XX[3],XX[4],XX[5],XX[6],XX[7],XX[8],XX[9],XX[10],1 };
	double Lnei[3 * 3] = { fx,-fx*tan(dp),-x0,0,fx / cos(dp) / (1 + ds),-y0,0,0,1 };
	double RT[3 * 4] = { 0 };
	InverseMatrix(Lnei, 3);
	MultiplyMatrix(Lnei, LL, RT, 3, 3, 4);

	double phi = atan(-RT[8] / RT[10]);
	double omega = asin(-RT[9]);
	double kappa = atan(RT[1] / RT[5]);
	double pixelsize = 0.0051966;//像素大小，用于单位换算

	cout << "单位权中误差  " << m0 << endl;
	cout << "内方位元素  " << x0 << "  " << y0 << "  " << fx << "  " << fy << endl;
	cout << "外方位元素  " << Xs << "  " << Ys << "  " << Zs << "  " << phi << "  " << omega << "  " << kappa << endl;
	cout << "畸变改正系数和中误差  " << k1 << "  " << m[11] << endl << k2 << "  " << m[12] << endl << p1 << "  " << m[13] << endl << p2 << "  " << m[14] << endl;
	cout << "直接线性变换和中误差  " << endl;
	for (int i = 0; i < 11; i++) {
		cout << i + 1 << "  " << l[i] << "  " << i + 1 << double(m[i] * 0.01) << endl;
	}
	cout << "像点观测值残差  " << endl;
	cout << "点号" << "  " << "x方向残差" << "  " << "y方向残差" << endl;
	for (int i = 0; i < num; i++) {
		cout << vec[3 * i] << "  " << Vii[2 * i] << "  " << Vii[2 * i + 1] << endl;
	}

	ofstream fileout("result1.txt");
	if (!fileout)
		cout << "wrong" << endl;
	fileout << "单位权中误差  " << m0 << endl;
	fileout << "内方位元素  " << x0 << "  " << y0 << "  " << fx << "  " << fy << endl;
	fileout << "外方位元素  " << Xs << "  " << Ys << "  " << Zs << "  " << phi << "  " << omega << "  " << kappa << endl;
	fileout << "畸变改正系数和中误差  " << k1 << "  " << m[11] << endl << k2 << "  " << m[12] << endl << p1 << "  " << m[13] << endl << p2 << "  " << m[14] << endl;
	fileout << "直接线性变换和中误差  " << endl;
	for (int i = 0; i < 11; i++) {
		fileout << i + 1 << "  " << l[i] << "  " << i + 1 << double(m[i] * 0.01) << endl;
	}
	fileout << "像点观测值残差  " << endl;
	fileout << "点号" << "  " << "x方向残差" << "  " << "y方向残差" << endl;
	for (int i = 0; i < num; i++) {
		fileout << vec[3 * i] << "  " << Vii[2 * i] << "  " << Vii[2 * i + 1] << endl;
	}

	//另一张像片的li和内方位元素
	double ll1[21] = { -0.0142034 , -0.0557398 , -0.00024109 , 256.572 , 0.00055014 , 0.000108932 , -0.0575101 , -9.28197 , 0.00202693 ,
		-0.000533327 , 2.35285e-05 , 0.000131276 , -3.09395e-07 , 0.000130023 , 2.33842e-05 , -0.21228 , 0.0674031, 0.000131276 , -3.09395e-07,
		0.000130023 , 2.33842e-05 };
	double ll2[21] = { 0.0 };
	for (int i = 0; i < 15; i++) {
		ll2[i] = XX[i];
	}
	ll2[15] = x0;
	ll2[16] = y0;
	ll2[17] = k1;
	ll2[18] = k2;
	ll2[19] = p1;
	ll2[20] = p2;
	int num1 = 0, num2 = 0;
	char *s4 = "4090.txt"; char *s5 = "4091.txt"; char *s3 = "wu.txt"; char *s2 = "wu92.txt";
	Coordinate *c1 = new Coordinate[500];
	Coordinate *c2 = new Coordinate[500];
	
	vector<double> vect1, vect2, vect3, vect4;
	//num1 = //ReadFile(s4, s3, c1);
		//num2 = //ReadFile(s5, s3, c2);

	ifstream infiles;
	infiles.open(s4);
	if (!infiles.is_open()) cout << "can't open this file!" << endl;
	double temps;
	while (infiles >> temps) {
		vect1.push_back(temps);
	}
	infiles.close(); //将像点坐标保存在vec中
	ifstream infiles1;
	infiles1.open(s5);
	if (!infiles1.is_open()) cout << "can't open this file!" << endl;
	double tempsq;
	while (infiles1 >> tempsq) {
		vect2.push_back(tempsq);
	}
	infiles1.close(); //将像点坐标保存在vec中
	ifstream infiles2;
	infiles2.open(s3);
	if (!infiles2.is_open()) cout << "can't open this file!" << endl;
	double tempsw;
	while (infiles2 >> tempsw) {
		vect3.push_back(tempsw);
	}
	infiles2.close(); //将像点坐标保存在vec中

	for (int i = 0; i < vect1.size() / 3; i++) {
		c1[i].num = int(vect1[3 * i]);
		c1[i].x = vect1[3 * i + 1];
		c1[i].y = vect1[3 * i + 2];
		c1[i].Z = -1 * vect3[4 * i + 1];
		c1[i].X = vect3[4 * i + 2];
		c1[i].Y = vect3[4 * i + 3];
	}

	for (int i = 0; i < vect2.size() / 3; i++) {
		c2[i].num = int(vect2[3 * i]);
		c2[i].x = vect2[3 * i + 1];
		c2[i].y = vect2[3 * i + 2];
		c2[i].Z = -1 * vect3[4 * i + 1];
		c2[i].X = vect3[4 * i + 2];
		c2[i].Y = vect3[4 * i + 3];
	}

	num1 = vect1.size() / 3;
	num2 = vect2.size() / 3;
	point *p = new point[500];
	int nnn = 0;
	for (int i = 0; i < num1; i++) {
		for (int j = 0; j < num2; j++) {
			if (c1[i].num == c2[j].num) {
				p[nnn].num = c1[i].num;
				p[nnn].x1 = c1[i].x * 0.0051966;
				p[nnn].y1 = c1[i].y * 0.0051966;
				p[nnn].x2 = c2[j].x * 0.0051966;
				p[nnn].y2 = c2[j].y * 0.0051966;
				p[nnn].X = -1 * c1[i].Z;
				p[nnn].Y = c1[i].X;
				p[nnn].Z = c1[i].Y;
				p[nnn].X_ = 0;//赋初值
				p[nnn].Y_ = 0;//赋初值
				p[nnn].Z_ = 0;//赋初值

				nnn++;
				break;
			}
		}
	}

	double r1s = 0, r12s = 0, r14s = 0, r2s = 0, r22s = 0, r24s = 0, temp1s = 0, temp2s = 0;
	double XYZ1s[3] = { 0.0 }, XYZ2s[3] = { 0.0 };
	double As[12] = { 0.0 }, ATs[12] = { 0.0 }, ATAs[9], Ls[4] = { 0.0 }, ATLs[3] = { 0.0 };
	double Ns[12] = { 0.0 }, NTs[12] = { 0.0 }, NTNs[9], Qs[4] = { 0.0 }, NTQs[3] = { 0.0 };

	for (int i = 0; i < nnn; i++) {
		//待定点像点坐标改正
		r1s = sqrt((p[i].x1 - ll1[15]) * (p[i].x1 - ll1[15]) + (p[i].y1 - ll1[16]) * (p[i].y1 - ll1[16]));
		r12s = r1s * r1s;
		r14s = r1s * r1s * r1s * r1s;
		p[i].x1 = p[i].x1 + (p[i].x1 - ll1[15]) * (ll1[17] * r12s + ll1[18] * r14s) + ll1[19] * (r12s + 2 * (p[i].x1 - ll1[15]) * (p[i].x1 - ll1[15])) + 2 * ll1[20] * (p[i].x1 - ll1[15]) * (p[i].y1 - ll1[16]);
		p[i].y1 = p[i].y1 + (p[i].y1 - ll1[16]) * (ll1[17] * r12s + ll1[18] * r14s) + ll1[20] * (r12s + 2 * (p[i].y1 - ll1[16]) * (p[i].y1 - ll1[16])) + 2 * ll1[19] * (p[i].x1 - ll1[15]) * (p[i].y1 - ll1[16]);

		r2s = sqrt((p[i].x2 - ll2[15]) * (p[i].x2 - ll2[15]) + (p[i].y2 - ll2[16]) * (p[i].y2 - ll2[16]));
		r22s = r2s * r2s;
		r24s = r2s * r2s * r2s * r2s;
		p[i].x2 = p[i].x2 + (p[i].x2 - ll2[15]) * (ll2[17] * r22s + ll2[18] * r24s) + ll2[19] * (r22s + 2 * (p[i].x2 - ll2[15]) * (p[i].x2 - ll2[15])) + 2 * ll2[20] * (p[i].x2 - ll2[15]) * (p[i].y2 - ll2[16]);
		p[i].y2 = p[i].y2 + (p[i].y2 - ll2[16]) * (ll2[17] * r22s + ll2[18] * r24s) + ll2[20] * (r22s + 2 * (p[i].y2 - ll2[16]) * (p[i].y2 - ll2[16])) + 2 * ll2[19] * (p[i].x2 - ll2[15]) * (p[i].y2 - ll2[16]);
	    
		//待定点物方空间坐标近似值计算
		As[0] = p[i].x1 * ll1[8] + ll1[0];
		As[1] = p[i].x1 * ll1[9] + ll1[1];
		As[2] = p[i].x1 * ll1[10] + ll1[2];
		As[3] = p[i].y1 * ll1[8] + ll1[4];
		As[4] = p[i].y1 * ll1[9] + ll1[5];
		As[5] = p[i].y1 * ll1[10] + ll1[6];
		As[6] = p[i].x2 * ll2[8] + ll2[0];
		As[7] = p[i].x2 * ll2[9] + ll2[1];
		As[8] = p[i].x2 * ll2[10] + ll2[2];
		As[9] = p[i].y2 * ll2[8] + ll2[4];
		As[10] = p[i].y2 * ll2[9] + ll2[5];
		As[11] = p[i].y2 * ll2[10] + ll2[6];

		Ls[0] = -p[i].x1 - ll1[3];
		Ls[1] = -p[i].y1 - ll1[7];
		Ls[2] = -p[i].x2 - ll2[3];
		Ls[3] = -p[i].y2 - ll2[7];

		TransposeMatrix(As, ATs, 4, 3);
		MultiplyMatrix(ATs, As, ATAs, 3, 4, 3);
		MultiplyMatrix(ATs, Ls, ATLs, 3, 4, 1);
		InverseMatrix(ATAs, 3);
		MultiplyMatrix(ATAs, ATLs, XYZ1s, 3, 3, 1);

		
			for (int j = 0; j < 3; j++)
				XYZ2s[j] = XYZ1s[j];

			temp1s = ll1[8] * XYZ1s[0] + ll1[9] * XYZ1s[1] + ll1[10] * XYZ1s[2] + 1;
			temp2s = ll2[8] * XYZ1s[0] + ll2[9] * XYZ1s[1] + ll2[10] * XYZ1s[2] + 1;

			Ns[0] = -(ll1[0] + ll1[8] * p[i].x1) / temp1s;
			Ns[1] = -(ll1[1] + ll1[9] * p[i].x1) / temp1s;
			Ns[2] = -(ll1[2] + ll1[10] * p[i].x1) / temp1s;
			Ns[3] = -(ll1[4] + ll1[8] * p[i].y1) / temp1s;
			Ns[4] = -(ll1[5] + ll1[9] * p[i].y1) / temp1s;
			Ns[5] = -(ll1[6] + ll1[10] * p[i].y1) / temp1s;
			Ns[6] = -(ll2[0] + ll2[8] * p[i].x2) / temp2s;
			Ns[7] = -(ll2[1] + ll2[9] * p[i].x2) / temp2s;
			Ns[8] = -(ll2[2] + ll2[10] * p[i].x2) / temp2s;
			Ns[9] = -(ll2[4] + ll2[8] * p[i].y2) / temp2s;
			Ns[10] = -(ll2[5] + ll2[9] * p[i].y2) / temp2s;
			Ns[11] = -(ll2[6] + ll2[10] * p[i].y2) / temp2s;

			Qs[0] = (ll1[3] + p[i].x1) / temp1s;
			Qs[1] = (ll1[7] + p[i].y1) / temp1s;
			Qs[2] = (ll2[3] + p[i].x2) / temp2s;
			Qs[3] = (ll2[7] + p[i].y2) / temp2s;

			TransposeMatrix(Ns, NTs, 4, 3);
			MultiplyMatrix(NTs, Ns, NTNs, 3, 4, 3);
			MultiplyMatrix(NTs, Qs, NTQs, 3, 4, 1);
			InverseMatrix(NTNs, 3);
			MultiplyMatrix(NTNs, NTQs, XYZ1s, 3, 3, 1);
			//for(int i=0;i<3;i++)
			//	printf("%lf\n",XYZ1[i]-XYZ2[i]); 



		p[i].X_ = XYZ1s[0];
		p[i].Y_ = XYZ1s[1];
		p[i].Z_ = XYZ1s[2];



	}

	double *vs = new double[3 * 13];
	for (int i = 0; i < 13; i++){
		vs[i * 3 + 0] = p[i].X - p[i].X_;
		vs[i * 3 + 1] = p[i].Y - p[i].Y_;
		vs[i * 3 + 2] = p[i].Z - p[i].Z_;
	}
	double mx0s = 0, my0s = 0, mz0s = 0;
	double vxs = 0;
	for (int i = 0; i < 13; i++){
		vxs += vs[i * 3 + 0] * vs[i * 3 + 0];
	}
	double vys = 0;
	for (int i = 0; i < 13; i++){
		vys += vs[3 * i + 1] * vs[3 * i + 1];
	}
	double vzs = 0;
	for (int i = 0; i < 13; i++){
		vzs += vs[i * 3 + 2] * vs[i * 3 + 2];
	}
	mx0s = sqrt(vzs / (4 * 13 - 3 * 13));
	my0s = sqrt(vxs / (4 * 13 - 3 * 13));
	mz0s = sqrt(vys / (4 * 13 - 3 * 13));


	double vvs = 0;
	for (int i = 0; i <  13; i++){
		vvs += vs[i] * vs[i];
	}
	double m0s = 0;
	m0s = (sqrt(vvs / 13));

	cout << "外精度  " << m0s << endl;
	cout << "x方向精度  " << mx0s << endl;
	cout << "y方向精度  " << my0s << endl;
	cout << "z方向精度  " << mz0s << endl;
	//cout << mx0s * 0.0051966 << endl;
	for (int i = 0; i < 13; i++) {
		cout << "解求出来的检查点物方坐标    " << p[i].num << "  " << p[i].X_ << "  " << p[i].Y_ << "  " << p[i].Z_ << endl;
		cout << "实际情况下的检查点物方坐标  " <<p[i].num << "  " << p[i].X << "  " << p[i].Y << "  " << p[i].Z << endl;
	}
	for (int i = 0; i < 13; i++) {
		double chazhi = sqrt((p[i].X - p[i].X_) * (p[i].X - p[i].X_) + (p[i].Y - p[i].Y_) * (p[i].Y - p[i].Y_) + (p[i].Z - p[i].Z_) * (p[i].Z - p[i].Z_));
		cout << "物方差值  " << chazhi << endl;
	}

	for (int i = 13; i < nnn; i++) {
		double chazhi1 = sqrt((p[13].X * 0.0051966 - p[i].X_ * 0.0051966) * (p[13].X * 0.0051966 - p[i].X_ * 0.0051966) + (p[13].Y * 0.0051966 - p[i].Y_ * 0.0051966) * (p[13].Y * 0.0051966 - p[i].Y_* 0.0051966) + (p[13].Z* 0.0051966 - p[i].Z_* 0.0051966) * (p[13].Z * 0.0051966 - p[i].Z_* 0.0051966));
		cout << "待定点坐标  " << p[i].num << "  " << p[i].X_ << "  " << p[i].Y_ << "  " << p[i].Z_ << "  " << endl;
		//cout << "待定点距离  "  << "  " << chazhi1 << endl;
	}

	for (int i = 14; i < nnn; i++) {
		double chazhi1 = sqrt((p[13].X * 0.0051966 - p[i].X_ * 0.0051966) * (p[13].X * 0.0051966 - p[i].X_ * 0.0051966) + (p[13].Y * 0.0051966 - p[i].Y_ * 0.0051966) * (p[13].Y * 0.0051966 - p[i].Y_* 0.0051966) + (p[13].Z* 0.0051966 - p[i].Z_* 0.0051966) * (p[13].Z * 0.0051966 - p[i].Z_* 0.0051966));
		cout << "待定点距离  " << p[i].num << "  " << chazhi1 << endl;
	}

	ofstream fileout1("result2.txt");
	if (!fileout1)
		cout << "wrong" << endl;
	fileout1 << "外精度  " << m0s << endl;
	fileout1 << "x方向精度  " << mx0s << endl;
	fileout1 << "y方向精度  " << my0s << endl;
	fileout1 << "z方向精度  " << mz0s << endl;
	//cout << mx0s * 0.0051966 << endl;
	for (int i = 0; i < 13; i++) {
		fileout1 << "解求出来的检查点物方坐标    " << p[i].num << "  " << p[i].X_ << "  " << p[i].Y_ << "  " << p[i].Z_ << endl;
		fileout1 << "实际情况下的检查点物方坐标  " << p[i].num << "  " << p[i].X << "  " << p[i].Y << "  " << p[i].Z << endl;
	}
	for (int i = 0; i < 13; i++) {
		double chazhi = sqrt((p[i].X - p[i].X_) * (p[i].X - p[i].X_) + (p[i].Y - p[i].Y_) * (p[i].Y - p[i].Y_) + (p[i].Z - p[i].Z_) * (p[i].Z - p[i].Z_));
		fileout1 << "物方差值  " << chazhi << endl;
	}

	for (int i = 13; i < nnn; i++) {
		double chazhi1 = sqrt((p[13].X * 0.0051966 - p[i].X_ * 0.0051966) * (p[13].X * 0.0051966 - p[i].X_ * 0.0051966) + (p[13].Y * 0.0051966 - p[i].Y_ * 0.0051966) * (p[13].Y * 0.0051966 - p[i].Y_* 0.0051966) + (p[13].Z* 0.0051966 - p[i].Z_* 0.0051966) * (p[13].Z * 0.0051966 - p[i].Z_* 0.0051966));
		fileout1 << "待定点坐标  " << p[i].num << "  " << p[i].X_ << "  " << p[i].Y_ << "  " << p[i].Z_ << "  " << endl;
		//fileout1 << "待定点距离  " << p[i + 1].num << "  " << chazhi1 << endl;
	}

	for (int i = 14; i < nnn; i++) {
		double chazhi1 = sqrt((p[13].X * 0.0051966 - p[i].X_ * 0.0051966) * (p[13].X * 0.0051966 - p[i].X_ * 0.0051966) + (p[13].Y * 0.0051966 - p[i].Y_ * 0.0051966) * (p[13].Y * 0.0051966 - p[i].Y_* 0.0051966) + (p[13].Z* 0.0051966 - p[i].Z_* 0.0051966) * (p[13].Z * 0.0051966 - p[i].Z_* 0.0051966));
		fileout1 << "待定点距离  " << p[i].num << "  " << chazhi1 << endl;
	}
	//找到两张像片的同名点

	//待定点像点坐标改正

	//待定点物方空间坐标近似值计算


	//现在我们有了内外方位元素 来算一下待定点的物方空间坐标
	/*vector<double> vec; //像点坐标保存向量
	vector<double> vec1; //物方空间点坐标保存向量
	ifstream infile;
	infile.open("912.txt");
	if (!infile.is_open()) cout << "can't open this file!" << endl;

	double temp;
	while (infile >> temp) {
		vec.push_back(temp);
	}

	infile.close(); //将像点坐标保存在vec中*/
	

	//cout << "外精度单位权中误差  " << m0s << endl;

	system("pause");
    return 0;
}

void TransposeMatrix(double *m1, double *m2, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			m2[j * row + i] = m1[i * col + j];
		}
	}
}

void AddMatrix(double *m1, double *m2, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			m2[i * col + j] = m2[i * col + j] + m1[i * col + j];
		}
	}
}

void SubMatrix(double *m1, double *m2, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			m2[i * col + j] = m2[i * col + j] - m1[i * col + j];
		}
	}
}

void InverseMatrix(double *a, int n) {
	int i, j, k;
	for (k = 0; k<n; k++) {
		for (i = 0; i<n; i++) {
			if (i != k)
				*(a + i*n + k) = -*(a + i*n + k) / (*(a + k*n + k));
		}
		*(a + k*n + k) = 1 / (*(a + k*n + k));
		for (i = 0; i<n; i++) {
			if (i != k) {
				for (j = 0; j<n; j++) {
					if (j != k)
						*(a + i*n + j) += *(a + k*n + j)* *(a + i*n + k);
				}
			}
		}
		for (j = 0; j<n; j++) {
			if (j != k)
				*(a + k*n + j) *= *(a + k*n + k);
		}
	}
}

void MultiplyMatrix(double *m1, double *m2, double *result, int m, int n, int l) {
	for (int i = 0; i<m; i++) {
		for (int j = 0; j<l; j++) {
			result[i*l + j] = 0.0;							//输出矩阵初始化
			for (int k = 0; k<n; k++)
				result[i*l + j] += m1[i*n + k] * m2[j + k*l];		//输出矩阵赋值
		}
	}
}

double  *DLT_(char *s1, char *s2)
{
	ifstream fp1;//像点文件
	ifstream fp2;//控制点文件

	double *temp1 = new double[1000];
	double *temp2 = new double[5000];

	//读入像点文件
	int n1 = 0;
	fp1.open(s1, ios::in | ios::out);
	for (n1 = 0; fp1.eof() == 0; n1++)
	{
		fp1 >> temp1[n1];

	}
	fp1.close();

	//读入控制点文件
	int n2 = 0;
	fp2.open(s2, ios::in | ios::out);
	for (n2 = 0; fp2.eof() == 0; n2++)
	{
		fp2 >> temp2[n2];

	}
	fp2.close();

	//像点个数
	int num1 = (n1 - 1) / 3;

	ImageCoordinate *imagec = new ImageCoordinate[num1];
	for (int i = 0; i<num1; i++)
	{
		imagec[i].num = temp1[3 * i + 0];
		imagec[i].x = temp1[3 * i + 1] - 4272 / 2.0;
		imagec[i].y = 2848 / 2.0 - temp1[3 * i + 2];
		//printf("%lf\n",imagec[i].y);
	}

	//控制点个数
	int num2 = (n2) / 4;

	ObjectCoordinate *objectc = new ObjectCoordinate[num2];
	for (int i = 0; i<num2; i++)
	{
		objectc[i].num = temp2[4 * i + 0];
		objectc[i].x = temp2[4 * i + 1];
		objectc[i].y = temp2[4 * i + 2];
		objectc[i].z = temp2[4 * i + 3];
	}

	int nn = 0;
	Coordinate *c = new Coordinate[num1];
	for (int i = 0; i<num2; i++)
		for (int j = 0; j<num1; j++)
		{
			if (objectc[i].num == imagec[j].num)
			{

				c[nn].num = imagec[j].num;
				c[nn].x = imagec[j].x*0.0051966;
				c[nn].y = imagec[j].y*0.0051966;
				c[nn].X = objectc[i].y;
				c[nn].Y = objectc[i].z;
				c[nn].Z = -1 * objectc[i].x;

				//printf("%lf %lf %lf %lf %lf\n",c[nn].x,c[nn].y,c[nn].X,c[nn].Y,c[nn].Z);
				nn++;
			}

		}

	int num = nn;//像点个数

				 //1.l系数近似值的解算
				 //选址前六个像点及其对应的控制点进行计算
	double *A = new double[11 * 2 * 6]; //大小12*11
	double *AT = new double[11 * 2 * 6]; //大小11*12
	double ATA[11 * 11] = { 0.0 };       //大小11*11
	double *L = new double[2 * 6];
	double ATL[11] = { 0.0 };
	double XX[15] = { 0.0 };
	double Qii[15] = { 0.0 };
	//double *v=new double[2*num];
	////初始化
	//for(int i=0;i<2*num;i++)
	//	v[i]=0;

	//l系数的初始值计算
	for (int i = 0; i<6; i++)
	{
		A[0 + 22 * i] = c[i].X;
		A[1 + 22 * i] = c[i].Y;
		A[2 + 22 * i] = c[i].Z;
		A[3 + 22 * i] = 1;
		A[4 + 22 * i] = 0;
		A[5 + 22 * i] = 0;
		A[6 + 22 * i] = 0;
		A[7 + 22 * i] = 0;
		A[8 + 22 * i] = c[i].x*c[i].X;
		A[9 + 22 * i] = c[i].x*c[i].Y;
		A[10 + 22 * i] = c[i].x*c[i].Z;
		A[11 + 22 * i] = 0;
		A[12 + 22 * i] = 0;
		A[13 + 22 * i] = 0;
		A[14 + 22 * i] = 0;
		A[15 + 22 * i] = c[i].X;
		A[16 + 22 * i] = c[i].Y;
		A[17 + 22 * i] = c[i].Z;
		A[18 + 22 * i] = 1;
		A[19 + 22 * i] = c[i].y*c[i].X;
		A[20 + 22 * i] = c[i].y*c[i].Y;
		A[21 + 22 * i] = c[i].y*c[i].Z;

		L[0 + 2 * i] = -c[i].x;
		L[1 + 2 * i] = -c[i].y;
	}

	TransposeMatrix(A, AT, 12, 11);
	MultiplyMatrix(AT, A, ATA, 11, 12, 11);
	MultiplyMatrix(AT, L, ATL, 11, 12, 1);
	InverseMatrix(ATA, 11);
	MultiplyMatrix(ATA, ATL, XX, 11, 11, 1);

	//l系数的精确值解算

	//计算x0，y0
	double x0, y0;  //内方位元素	
	x0 = -(XX[0] * XX[8] + XX[1] * XX[9] + XX[2] * XX[10]) / (XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10]);
	y0 = -(XX[4] * XX[8] + XX[5] * XX[9] + XX[6] * XX[10]) / (XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10]);
	//计算fx,fy
	double fx, fy;  //x向主距，y向主距
	double tem1 = 0, tem2 = 0, tem3 = 0;

	tem1 = XX[0] * XX[0] + XX[1] * XX[1] + XX[2] * XX[2];
	tem2 = XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10];
	tem3 = XX[4] * XX[4] + XX[5] * XX[5] + XX[6] * XX[6];

	fx = sqrt(tem1 / tem2 - x0*x0);

	double r3 = 0, ds = 0;
	r3 = 1 / (XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10]);
	ds = sqrt((r3*r3*tem1 - x0*x0) / (r3*r3*tem3 - y0*y0)) - 1;
	fy = fx / (1 + ds);

	double fx_, fy_;
	double *M = new double[15 * 2 * num];
	double *MT = new double[15 * 2 * num];
	double *W = new double[2 * num];
	double MTM[225] = { 0.0 };
	double MTW[15] = { 0.0 };
	double aa = 0, bb = 0, cc = 0;
	do
	{
		fx_ = fx;
		fy_ = fy;
		for (int i = 0; i<num; i++)
		{

			////读出数组中一组数据
			//double x,y,X,Y,Z;
			//x=c[i].x;
			//y=c[i].y;
			//X=c[i].X;
			//Y=c[i].Y;
			//Z=c[i].Z;
			double temp = XX[8] * c[i].X + XX[9] * c[i].Y + XX[10] * c[i].Z + 1;
			//计算向径
			double r = sqrt((c[i].x - x0)*(c[i].x - x0) + (c[i].y - y0)*(c[i].y - y0));
			double r2 = pow(r, 2);//r的平方
			double r4 = pow(r, 4);//r的四次方
			M[30 * i + 0] = -c[i].X / temp;
			M[30 * i + 1] = -c[i].Y / temp;
			M[30 * i + 2] = -c[i].Z / temp;
			M[30 * i + 3] = -1 / temp;
			M[30 * i + 4] = 0;
			M[30 * i + 5] = 0;
			M[30 * i + 6] = 0;
			M[30 * i + 7] = 0;
			M[30 * i + 8] = -(c[i].x*c[i].X) / temp;
			M[30 * i + 9] = -(c[i].x*c[i].Y) / temp;
			M[30 * i + 10] = -(c[i].x*c[i].Z) / temp;
			M[30 * i + 11] = -(c[i].x - x0)*r2;
			M[30 * i + 12] = -(c[i].x - x0)*r4;
			M[30 * i + 13] = -r2 - 2 * (c[i].x - x0)*(c[i].x - x0);
			M[30 * i + 14] = -2 * (c[i].x - x0)*(c[i].y - y0);

			M[30 * i + 15] = 0;
			M[30 * i + 16] = 0;
			M[30 * i + 17] = 0;
			M[30 * i + 18] = 0;
			M[30 * i + 19] = -c[i].X / temp;
			M[30 * i + 20] = -c[i].Y / temp;
			M[30 * i + 21] = -c[i].Z / temp;
			M[30 * i + 22] = -1 / temp;
			M[30 * i + 23] = -(c[i].y*c[i].X) / temp;
			M[30 * i + 24] = -(c[i].y*c[i].Y) / temp;
			M[30 * i + 25] = -(c[i].y*c[i].Z) / temp;
			M[30 * i + 26] = -(c[i].y - y0)*r2;
			M[30 * i + 27] = -(c[i].y - y0)*r4;
			M[30 * i + 28] = -2 * (c[i].x - x0)*(c[i].y - y0);
			M[30 * i + 29] = -r2 - 2 * (c[i].y - y0)*(c[i].y - y0);

			W[2 * i + 0] = c[i].x / temp;
			W[2 * i + 1] = c[i].y / temp;
		}
		TransposeMatrix(M, MT, 2 * num, 15);
		MultiplyMatrix(MT, M, MTM, 15, 2 * num, 15);
		MultiplyMatrix(MT, W, MTW, 15, 2 * num, 1);
		InverseMatrix(MTM, 15);
		MultiplyMatrix(MTM, MTW, XX, 15, 15, 1);

		//计算
		r3 = 1 / (XX[8] * XX[8] + XX[9] * XX[9] + XX[10] * XX[10]);
		x0 = -r3*(XX[0] * XX[8] + XX[1] * XX[9] + XX[2] * XX[10]);
		y0 = -r3*(XX[4] * XX[8] + XX[5] * XX[9] + XX[6] * XX[10]);
		aa = r3*(XX[0] * XX[0] + XX[1] * XX[1] + XX[2] * XX[2]) - x0*x0;
		bb = r3*(XX[4] * XX[4] + XX[5] * XX[5] + XX[6] * XX[6]) - y0*y0;
		cc = r3*(XX[0] * XX[4] + XX[1] * XX[5] + XX[2] * XX[6]) - x0*y0;

		fx = sqrt((aa*bb - cc*cc) / bb);
		fy = sqrt((aa*bb - cc*cc) / aa);

	} while (fabs(fx - fx_)>0.01 || fabs(fy - fy_)>0.01);

	//计算单位权中误差
	double *V = new double[2 * num];
	MultiplyMatrix(M, XX, V, 2 * num, 15, 1);
	SubMatrix(V, W, 2 * num, 1);
	double v = 0;
	for (int i = 0; i<2 * num; i++)
	{
		v += V[i] * V[i];
	}
	double m0 = sqrt(v / (2 * num - 15));



	for (int i = 0; i<15; i++)
	{
		Qii[i] = MTM[15 * i + i];
	}
	double m[15] = { 0.0 };
	for (int i = 0; i<15; i++)
	{
		m[i] = m0*sqrt(Qii[i]);
	}


	//内方位元素
	x0 = -r3*(XX[0] * XX[8] + XX[1] * XX[9] + XX[2] * XX[10]);
	y0 = -r3*(XX[4] * XX[8] + XX[5] * XX[9] + XX[6] * XX[10]);
	//11个线性变换系数
	double l[11] = { 0 };
	for (int i = 0; i<11; i++)
		l[i] = XX[i];
	//四个畸变系数
	double k1 = XX[11];
	double k2 = XX[12];
	double p1 = XX[13];
	double p2 = XX[14];

	//外方位线元素,书本150页
	double wai[9] = { 0.0 };
	wai[0] = XX[0]; wai[1] = XX[1]; wai[2] = XX[2];
	wai[4] = XX[3]; wai[1] = XX[4]; wai[5] = XX[6];
	wai[8] = XX[6]; wai[1] = XX[7]; wai[8] = XX[10];

	double waiL[3] = { 0.0 };
	waiL[0] = -1 * XX[3]; waiL[1] = -1 * XX[7]; waiL[2] = -1;

	double w[3] = { 0 };
	InverseMatrix(wai, 3);
	MultiplyMatrix(wai, waiL, w, 3, 3, 1);

	//化为地面摄影测量坐标
	double Xs = w[0], Ys = w[1], Zs = w[2];
	double temz = Zs;
	Zs = Ys;
	Ys = Xs;
	Xs = -1 * temz;

	ds = sqrt(aa / bb) - 1;
	double dp = asin(sqrt(cc*cc / aa / bb));

	double LL[3 * 4] = { XX[0],XX[1],XX[2],XX[3],XX[4],XX[5],XX[6],XX[7],XX[8],XX[9],XX[10],1 };
	double Lnei[3 * 3] = { fx,-fx*tan(dp),-x0,0,fx / cos(dp) / (1 + ds),-y0,0,0,1 };
	double RT[3 * 4] = { 0 };
	InverseMatrix(Lnei, 3);
	MultiplyMatrix(Lnei, LL, RT, 3, 3, 4);

	double phi = atan(-RT[8] / RT[10]);
	double omega = asin(-RT[9]);
	double kappa = atan(RT[1] / RT[5]);
	double pixelsize = 0.0051966;//像素大小，用于单位换算
	
	double *need = new double[21];
	for (int i = 0; i<15; i++)
		need[i] = XX[i];
	need[15] = x0;
	need[16] = y0;
	need[17] = k1;
	need[18] = k2;
	need[19] = p1;
	need[20] = p2;

	return need;
}

int ReadFile(char *s1, char *s2, Coordinate *pk)
{
	//读入像点文件
	int n1 = 0;
	ifstream fp1;//像点文件
	ifstream fp2;//控制点文件

	double *temp1 = new double[1000];
	double *temp2 = new double[5000];
	fp1.open(s1, ios::in | ios::out);
	for (n1 = 0; fp1.eof() == 0; n1++)
	{
		fp1 >> temp1[n1];

	}
	fp1.close();

	//读入控制点文件
	int n2 = 0;
	fp2.open(s2, ios::in | ios::out);
	for (n2 = 0; fp2.eof() == 0; n2++)
	{
		fp2 >> temp2[n2];

	}
	fp2.close();

	//像点个数
	int num1 = (n1 - 1) / 3;

	ImageCoordinate *imagec = new ImageCoordinate[num1];
	for (int i = 0; i<num1; i++)
	{
		imagec[i].num = temp1[3 * i + 0];
		imagec[i].x = temp1[3 * i + 1] - 4272 / 2.0;
		imagec[i].y = 2848 / 2.0 - temp1[3 * i + 2];
	}

	//控制点个数
	int num2 = (n2) / 4;

	ObjectCoordinate *objectc = new ObjectCoordinate[num2];
	for (int i = 0; i<num2; i++)
	{
		objectc[i].num = temp2[4 * i + 0];
		objectc[i].x = temp2[4 * i + 1];
		objectc[i].y = temp2[4 * i + 2];
		objectc[i].z = temp2[4 * i + 3];
	}

	int nn = 0;
	for (int i = 0; i<num2; i++)
		for (int j = 0; j<num1; j++)
		{
			if (objectc[i].num == imagec[j].num)
			{

				pk[nn].num = imagec[j].num;
				pk[nn].x = imagec[j].x*0.0051966;
				pk[nn].y = imagec[j].y*0.0051966;
				pk[nn].X = objectc[i].y;
				pk[nn].Y = objectc[i].z;
				pk[nn].Z = -1 * objectc[i].x;

				nn++;
			}

		}
	return nn;
}