//-------------------------------------------------------------------------------------------------
//翼素定義クラス
//主翼の離散要素の情報を管理する
//-------------------------------------------------------------------------------------------------
#pragma once

class WingElementClass
{
public:
	WingElementClass();
	~WingElementClass();

	//-------------------------------------------
	//パネル座標
	//-------------------------------------------
	double lf[3], lr[3], rf[3], rr[3];					//パネル角の座標
	double pl[3], pr[3];								//パネルの渦配置座標
	double pcp[3];										//コントロールポイント座標

	//-------------------------------------------
	//パネル定数
	//-------------------------------------------
	double phi;											//パネル上反角
	double seta;										//パネル迎角

	double ds;											//パネル面積
	double dl;											//パネル幅
	double dc;											//パネルコード長

	double dzdx;										//パネルキャンバー傾斜
	double ganma;										//パネル循環
	double downwash;									//パネル吹き下ろし
	double cp;											//パネル圧力
};

