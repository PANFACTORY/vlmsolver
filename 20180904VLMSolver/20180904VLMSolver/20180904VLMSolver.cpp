//*************************************************************************************************
//Name		: VLM Solver
//Purpose	: Culculate AeroDynamicForces with VLM method
//
//Author	: Yuta Tanabe
//
//Created	: 06 / 13 / 2018
//Modified	: 09 / 04 / 2018 
//Copyright : (c)YutaTanabe 2018
//---------------------------修正・検討が必要な個所------------------------------------------------
//①連立方程式のソルバはガウスの消去法よりLU分解などの逆行列を求めるタイプにした方が良い．
//　迎角を変えたときの解析が速くなる
//②キャンバーを求めるときに端は片側差分になるが入れられていない
//　（よほどx方向分割数を増やさない限りは大丈夫か?）
//③座標軸の張り方は主流方向にx軸正，右翼側へy軸正，鉛直上向きにz軸正
//④パネルの格納方法を再検討する必要あり（現状ではx軸方向への積分が面倒）
//⑤VLMにおける粘性考慮の方法は？（XFLR5のようにCLから参照するのでよいのか）
//*************************************************************************************************

#include "stdafx.h"

#include "WingElementClass.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

//*************************************************************************************************
//文字列分割用の関数
//		_str		：分割対称文字列
//		_sep		：区切り文字
//*************************************************************************************************
vector<string> split(const string &_str, char _sep)
{
	vector<string> v;
	stringstream ss(_str);
	string buffer;
	while (getline(ss, buffer, _sep)) {
		v.push_back(buffer);
	}
	return v;
}

//*************************************************************************************************
//ガウスの消去法
//		_n			：要素数
//		_a			：係数行列
//		_b			：係数ベクトル
//*************************************************************************************************
void MSolve(int _n, double **_a, double *_b) {
	// 前進消去
	double p;
	for (int k = 0; k < _n - 1; k++) {
		for (int i = k + 1; i < _n; i++) {
			p = -_a[i][k] / _a[k][k];
			for (int j = k + 1; j < _n; j++) {
				_a[i][j] += p * _a[k][j];
			}
			_b[i] += p * _b[k];
		}
	}
	// 後退代入
	for (int i = _n - 1; i > -1; i--) {
		for (int j = i + 1; j < _n; j++) {
			_b[i] -= _a[i][j] * _b[j];
		}
		_b[i] = _b[i] / _a[i][i];
	}
}

//*************************************************************************************************
//VLMソルバ
//		_element	：パネル格納配列
//		_alpha		：迎え角
//*************************************************************************************************
void vlm_solver(vector<WingElementClass> &_element, double _alpha) {
	//ビオサバールの影響係数を計算
	double** qijz = new double*[_element.size()];
	for (int i = 0; i < _element.size(); i++) {
		qijz[i] = new double[_element.size()];

		for (int j = 0; j < _element.size(); j++) {
			double r0[3], r1[3], r2[3];
			for (int n = 0; n < 3; n++) {
				r0[n] = _element[j].pr[n] - _element[j].pl[n];
				r1[n] = _element[i].pcp[n] - _element[j].pl[n];
				r2[n] = _element[i].pcp[n] - _element[j].pr[n];
			}
			double r1nrel = sqrt(pow(r1[0], 2.0) + pow(r1[1], 2.0) + pow(r1[2], 2.0));
			double r2nrel = sqrt(pow(r2[0], 2.0) + pow(r2[1], 2.0) + pow(r2[2], 2.0));
			//束縛渦
			double qpsirel = 0.0, phidrel[3];
			for (int n = 0; n < 3; n++) {
				qpsirel += r0[n] * (r1[n] / r1nrel - r2[n] / r2nrel);
				phidrel[n] = r1[(n + 1) % 3] * r2[(n + 2) % 3] - r2[(n + 1) % 3] * r1[(n + 2) % 3];
			}
			double phinrel = sqrt(pow(phidrel[0], 2.0) + pow(phidrel[1], 2.0) + pow(phidrel[2], 2.0));
			double phirel[3], vabrel[3];
			for (int n = 0; n < 3; n++) {
				phirel[n] = phidrel[n] / pow(phinrel, 2.0);
				vabrel[n] = qpsirel * phirel[n];
			}
			//後曳渦左
			double vaexrel = 0.0;
			double vaeyrel = r1[2] / (pow(r1[2], 2.0) + pow(r1[1], 2.0)) * (1 + r1[0] / r1nrel);
			double vaezrel = -r1[1] / (pow(r1[2], 2.0) + pow(r1[1], 2.0)) * (1 + r1[0] / r1nrel);
			//後曳渦右
			double vebxrel = 0.0;
			double vebyrel = -r2[2] / (pow(r2[2], 2.0) + pow(r2[1], 2.0)) * (1 + r2[0] / r2nrel);
			double vebzrel = r2[1] / (pow(r2[2], 2.0) + pow(r2[1], 2.0)) * (1 + r2[0] / r2nrel);

			double vijx = (vabrel[0] + vaexrel + vebxrel) / (4.0 * M_PI);
			double vijy = (vabrel[1] + vaeyrel + vebyrel) / (4.0 * M_PI);
			double vijz = (vabrel[2] + vaezrel + vebzrel) / (4.0 * M_PI);

			//翼素ベクトルに垂直な成分のみを計算
			qijz[i][j] = -vijx*cos(_element[i].phi*M_PI / 180.0)*sin(_element[i].seta*M_PI / 180.0)
				- vijy*sin(_element[i].phi*M_PI / 180.0)*cos(_element[i].seta*M_PI / 180.0)
				+ vijz*cos(_element[i].phi*M_PI / 180.0)*cos(_element[i].seta*M_PI / 180.0);
		}
	}

	//連立方程式を作成
	double** c = new double*[_element.size()];
	double* d = new double[_element.size()];
	for (int i = 0; i < _element.size(); i++) {
		c[i] = new double[_element.size()];
		for (int j = 0; j < _element.size(); j++) {
			c[i][j] = qijz[i][j];
		}
		d[i] = sin((_alpha + _element[i].seta)*M_PI / 180.0 + _element[i].dzdx);
	}

	//ガウスの消去法で連立方程式を解く
	MSolve(_element.size(), c, d);

	//解を_elementに代入
	for (int i = 0; i < _element.size(); i++) {
		_element[i].ganma = d[i];
	}

	//吹き下ろしおよび圧力分布を計算
	for (int i = 0; i < _element.size(); i++) {
		_element[i].downwash = 0.0;
		for (int j = 0; j < _element.size(); j++) {
			_element[i].downwash += qijz[i][j] * _element[j].ganma * _element[j].dl;
		}
		_element[i].cp = -2.0*fabs(_element[i].ganma) / _element[i].dc;
	}

	//メモリの開放
	for (int i = 0; i < _element.size(); i++) {
		delete(qijz[i]);
		delete(c[i]);
	}
	delete(qijz);
	delete(c);
	delete(d);
}

//*************************************************************************************************
//キャンバーを取得する
//		_cpos		：コード長%位置
//		_foilname	：翼型名
//Foilsフォルダ内に翼型datファイルが必要
//*************************************************************************************************
double getCamber(double _cpos, string _foilname) {
	//----------ファイルから翼型の座標データを読み込む----------
	ifstream fin("./Foils/" + _foilname + ".dat");
	if (!fin) {														
		cout << "Foil File Open Error ./Foils/" + _foilname + ".dat\n";
		return 0.0;										//ファイル読み込みに失敗したときはキャンバー傾斜0で値を返す
	}
	stringstream strstream;									
	strstream << fin.rdbuf();
	fin.close();
	string data(strstream.str());
	vector<string> linedata = split(data, '\n');	
	int plotnum = linedata.size() - 1;

	double *x = new double[plotnum];
	double *y = new double[plotnum];
	for (int j = 0; j < plotnum; j++) {				
		istringstream ssi(linedata[j + 1]);
		ssi >> x[j] >> y[j];
	}

	//----------前縁探索----------
	int leadingedgei = 0;
	for (int i = 1; i < plotnum - 1; i++) {
		if (x[i] <= x[i - 1] && x[i] < x[i + 1]) {
			leadingedgei = i;
			break;
		}
	}

	//----------上面側両隣の点を取得----------
	int upperneari = 0;
	for (int i = 1; i < leadingedgei; i++) {
		if (x[i] > _cpos && _cpos >= x[i + 1]) {
			upperneari = i;
			break;
		}
	}

	//----------下面の点をもとにキャンバー座標を取得----------
	double canbery[4] = { 0, 0, 0, 0 };
	for (int k = 0; k < 4; k++) {
		int bottomi = leadingedgei;
		for (int j = leadingedgei; j < plotnum; j++) {
			if (x[j - 1] < x[upperneari + k - 1] && x[upperneari + k - 1] <= x[j]) {
				bottomi = j;
				break;
			}
		}
		double bottomy = (y[bottomi] - y[bottomi - 1])*(x[upperneari + k - 1] - x[bottomi - 1]) / (x[bottomi] - x[bottomi - 1]) + y[bottomi - 1];
		canbery[k] = 0.5*(y[upperneari + k - 1] + bottomy);
	}

	//----------キャンバーの値をもとに二次精度差分近似でキャンバー傾斜を取得----------
	double dcanber[2] = { 0,0 };
	for (int k = 0; k < 2; k++) {
		dcanber[k] = (canbery[k + 2] - canbery[k]) / (x[upperneari + k - 1] - x[upperneari + k + 1]);
	}

	return (dcanber[1] - dcanber[0]) * (_cpos - x[upperneari]) / (x[upperneari + 1] - x[upperneari]) + dcanber[0];
}

//*************************************************************************************************
//パネル生成
//		_element	：パネル格納配列
//		_fnname		：翼ファイル名
//		_sym		：trueの場合左右の複製を行う
//		_setpos		：翼の中心座標
//*************************************************************************************************
int getpanels(vector<WingElementClass> &_element, char *_fname, bool _sym, double *_setpos) {
	int elementnum = _element.size();

	//----------主翼データを取り込む----------
	ifstream fin(_fname);
	if (!fin) {
		std::cout << "File Open Error" << std::endl;
		return 1;
	}
	stringstream strstream;
	strstream << fin.rdbuf();
	fin.close();
	string data(strstream.str());
	vector<string> linedata = split(data, '\n');

	//----------パネルデータを生成----------
	double rf[3], rr[3], tf[3], tr[3];								//各セクションの角の座標

																	//翼根のデータを読み込む
	istringstream ssk(linedata[1]);
	double yk, ck, osk, dk, ak, xnumk, ynumk, dam0, dam1, dam2, ykp1, ckp1, oskp1, dkp1, akp1, xnumkp1, ynumkp1;
	string fnamek, fnamekp1;
	ssk >> yk >> ck >> osk >> dk >> ak >> xnumk >> ynumk >> dam0 >> dam1 >> dam2 >> fnamek;
	ak += _setpos[3];

	double sx = 0.0, sy = yk, sz = 0.0;

	rf[0] = osk*cos(ak*M_PI / 180.0) + _setpos[0];
	rf[1] = sy + _setpos[1];
	rf[2] = -osk*sin(ak*M_PI / 180.0) + _setpos[2];
	rr[0] = (ck + osk)*cos(ak*M_PI / 180.0) + _setpos[0];
	rr[1] = sy + _setpos[1];
	rr[2] = -(ck + osk)*sin(ak*M_PI / 180.0) + _setpos[2];

	//セクション毎にデータを読み込む
	for (int k = 2; k < linedata.size(); k++) {
		//セクションの翼端のデータを読み込む
		istringstream ssk(linedata[k]);
		ssk >> ykp1 >> ckp1 >> oskp1 >> dkp1 >> akp1 >> xnumkp1 >> ynumkp1 >> dam0 >> dam1 >> dam2 >> fnamekp1;
		akp1 += _setpos[3];

		sy += (ykp1 - yk) * cos(dk*M_PI / 180.0);
		sz += (ykp1 - yk) * sin(dk*M_PI / 180.0);

		tf[0] = oskp1*cos(akp1*M_PI / 180.0) + _setpos[0];
		tf[1] = sy - oskp1*sin(akp1*M_PI / 180.0)*sin(dk*M_PI / 180.0) + _setpos[1];
		tf[2] = sz - oskp1*sin(akp1*M_PI / 180.0)*cos(dk*M_PI / 180.0) + _setpos[2];
		tr[0] = (ckp1 + oskp1)*cos(akp1*M_PI / 180.0) + _setpos[0];
		tr[1] = sy - (ckp1 + oskp1)*sin(akp1*M_PI / 180.0)*sin(dk*M_PI / 180.0) + _setpos[1];
		tr[2] = sz - (ckp1 + oskp1)*sin(akp1*M_PI / 180.0)*cos(dk*M_PI / 180.0) + _setpos[2];

		//y方向の分割
		for (int j = 0; j < ynumk; j++) {
			double tmprf[3], tmprr[3], tmptf[3], tmptr[3];
			for (int n = 0; n < 3; n++) {
				tmprf[n] = (tf[n] - rf[n])*j / (double)ynumk + rf[n];
				tmptf[n] = (tf[n] - rf[n])*(j + 1) / (double)ynumk + rf[n];
				tmprr[n] = (tr[n] - rr[n])*j / (double)ynumk + rr[n];
				tmptr[n] = (tr[n] - rr[n])*(j + 1) / (double)ynumk + rr[n];
			}

			//x方向の分割
			for (int i = 0; i < xnumk; i++) {
				WingElementClass tmpelement;

				for (int n = 0; n < 3; n++) {
					//パネルの角の座標
					tmpelement.lf[n] = (tmprr[n] - tmprf[n])*i / (double)xnumk + tmprf[n];
					tmpelement.lr[n] = (tmprr[n] - tmprf[n])*(i + 1) / (double)xnumk + tmprf[n];
					tmpelement.rf[n] = (tmptr[n] - tmptf[n])*i / (double)xnumk + tmptf[n];
					tmpelement.rr[n] = (tmptr[n] - tmptf[n])*(i + 1) / (double)xnumk + tmptf[n];

					//束縛渦端の座標
					tmpelement.pl[n] = (tmpelement.lf[n] * 3.0 + tmpelement.lr[n] * 1.0) / 4.0;
					tmpelement.pr[n] = (tmpelement.rf[n] * 3.0 + tmpelement.rr[n] * 1.0) / 4.0;

					//パネルコントロールポイント座標を計算
					tmpelement.pcp[n] = ((tmpelement.lf[n] * 1.0 + tmpelement.lr[n] * 3.0) + (tmpelement.rf[n] * 1.0 + tmpelement.rr[n] * 3.0)) / 8.0;
				}

				tmpelement.phi = dk;
				tmpelement.seta = (akp1 - ak)*(j + 0.5) / (double)ynumk + ak;
				tmpelement.dl = (ykp1 - yk) / (double)ynumk;
				tmpelement.dc = ((ckp1 - ck)*(j + 0.5) / (double)ynumk + ck) / (double)xnumk;
				tmpelement.ds = tmpelement.dl * tmpelement.dc;
				double camberl = atan(getCamber((0.75 + i) / (double)xnumk, fnamek));
				double camberr = atan(getCamber((0.75 + i) / (double)xnumk, fnamekp1));
				tmpelement.dzdx = (camberr - camberl)*(j + 0.5) / (double)ynumk + camberl;

				_element.push_back(tmpelement);
			}
		}

		for (int n = 0; n < 3; n++) {
			rf[n] = tf[n];
			rr[n] = tr[n];
		}
		yk = ykp1;	ck = ckp1;	osk = oskp1;	dk = dkp1;	ak = akp1;	xnumk = xnumkp1;	ynumk = ynumkp1;	fnamek = fnamekp1;
	}

	//----------対称パネルの生成----------
	if (_sym) {
		for (int i = _element.size() - 1; i >= elementnum; i--) {
			WingElementClass tmpelement;

			for (int n = 0; n < 3; n++) {
				//パネルの角の座標
				tmpelement.lf[n] = _element[i].rf[n] * pow(-1.0, n);
				tmpelement.lr[n] = _element[i].rr[n] * pow(-1.0, n);
				tmpelement.rf[n] = _element[i].lf[n] * pow(-1.0, n);
				tmpelement.rr[n] = _element[i].lr[n] * pow(-1.0, n);

				//束縛渦端の座標
				tmpelement.pl[n] = _element[i].pr[n] * pow(-1.0, n);
				tmpelement.pr[n] = _element[i].pl[n] * pow(-1.0, n);

				//パネルコントロールポイント座標を計算
				tmpelement.pcp[n] = _element[i].pcp[n] * pow(-1.0, n);
			}

			tmpelement.phi = -_element[i].phi;
			tmpelement.seta = _element[i].seta;
			tmpelement.dl = _element[i].dl;
			tmpelement.dc = _element[i].dc;
			tmpelement.ds = _element[i].ds;
			tmpelement.dzdx = _element[i].dzdx;

			_element.push_back(tmpelement);
		}
	}

	return 0;
}

//*************************************************************************************************
//gnuplotの初期化
//*************************************************************************************************
void initgnuplot(FILE *&_fp, char *_title, double *_range) {
	_fp = _popen("gnuplot -persist", "w");									//パイプを開き、gnuplotの立ち上げ
	fprintf(_fp, "set title '%s'\n", _title);								//グラフタイトル
	fprintf(_fp, "set multiplot\n");										//マルチプロットモード
	fprintf(_fp, "set view equal xyz\n");
	fprintf(_fp, "set xrange [%lf:%lf]\n", _range[0], _range[1]);			//範囲の指定
	fprintf(_fp, "set yrange [%lf:%lf]\n", _range[2], _range[3]);
	fprintf(_fp, "set zrange [%lf:%lf]\n", _range[4], _range[5]);
	fprintf(_fp, "set xlabel \"x\"\n");										//ラベル表示
	fprintf(_fp, "set ylabel \"y\"\n");
	fprintf(_fp, "set zlabel \"z\"\n");
	fprintf(_fp, "set ticslevel 0\n");
	fprintf(_fp, "set style arrow 1 nohead linecolor rgb 'black'\n");
}

//*************************************************************************************************
//パネルの描画
//*************************************************************************************************
void showpanels(vector<WingElementClass> _element) {
	FILE *gp;
	double range[6] = { -1.0, 1.0, -13.0, 13.0, 0.0, 2.0 };
	initgnuplot(gp, "Panels", range);

	//パネルの描画
	for (int i = 0; i < _element.size(); i++) {
		fprintf(gp, "set arrow from %lf,%lf,%lf to %lf,%lf,%lf arrowstyle 1\n", _element[i].lf[0], _element[i].lf[1], _element[i].lf[2], _element[i].rf[0], _element[i].rf[1], _element[i].rf[2]);
		fprintf(gp, "set arrow from %lf,%lf,%lf to %lf,%lf,%lf arrowstyle 1\n", _element[i].rf[0], _element[i].rf[1], _element[i].rf[2], _element[i].rr[0], _element[i].rr[1], _element[i].rr[2]);
		fprintf(gp, "set arrow from %lf,%lf,%lf to %lf,%lf,%lf arrowstyle 1\n", _element[i].rr[0], _element[i].rr[1], _element[i].rr[2], _element[i].lr[0], _element[i].lr[1], _element[i].lr[2]);
		fprintf(gp, "set arrow from %lf,%lf,%lf to %lf,%lf,%lf arrowstyle 1\n", _element[i].lr[0], _element[i].lr[1], _element[i].lr[2], _element[i].lf[0], _element[i].lf[1], _element[i].lf[2]);
	}

	fprintf(gp, "splot '-' with points pointtype 1\n");
	fprintf(gp, "e\n");
	fprintf(gp, "set nomultiplot\n");						// マルチプロットモード終了
	fprintf(gp, "exit\n");									// gnuplotの終了
	fflush(gp);
}

//*************************************************************************************************
//Cp分布の描画
//*************************************************************************************************
void showcp(vector<WingElementClass> _element) {
	FILE *fp;
	double range[6] = { -1.0, 1.0, -13.0, 13.0, 0.0, 2.0 };
	initgnuplot(fp, "Cp Distribution", range);

	//Cp分布の描画
	for (int i = 0; i < _element.size(); i++) {
		fprintf(fp, "set arrow from %lf,%lf,%lf to %lf,%lf,%lf\n", _element[i].pcp[0], _element[i].pcp[1], 0.0, _element[i].pcp[0], _element[i].pcp[1], _element[i].cp);
	}

	fprintf(fp, "splot '-' with points pointtype 1\n");
	fprintf(fp, "e\n");
	fprintf(fp, "set nomultiplot\n");						// マルチプロットモード終了
	fprintf(fp, "exit\n");									// gnuplotの終了
	fflush(fp);
}

//*************************************************************************************************
//メイン処理
//*************************************************************************************************
int main()
{
	//----------VLMの設定----------
	double mainwingpos[4] = { 0.0, 0.0, 0.0, 0.0 };					//主翼配置位置と取付角
	double htailwingpos[4] = { 3.5, 0.0, 0.0, -2.0 };				//水平尾翼配置位置と取付角

	//----------VLM用のパネルを生成----------
	cout << "Make Panels\n";
	vector<WingElementClass> elements;								//翼要素配列	
	getpanels(elements, "./mainwing.xwimp", true, mainwingpos);
	getpanels(elements, "./elevetor.xwimp", true, htailwingpos);

	//----------VLMによる解析----------
	cout << "Culculate AeroDynamicCoefficients\n";
	vlm_solver(elements, 6.0);

	double lift = 0.0, idrag = 0.0, space = 0.0;
	for (int i = 0; i < elements.size(); i++) {
		//揚力の計算
		lift += fabs(elements[i].ganma * elements[i].dl);
		//誘導抗力の計算
		idrag += fabs(elements[i].ganma) * elements[i].downwash * elements[i].dl;
		//翼面積の計算
		space += elements[i].ds;
	}
	cout << "CL=" << lift / (0.5*space) << "\tCDi=" << idrag / (0.5*space) << "\n";

	//----------結果の表示----------
	//翼パネルの表示
	showpanels(elements);

	//圧力分布の表示
	showcp(elements);

	return 0;
}