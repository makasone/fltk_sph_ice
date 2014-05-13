#ifndef _QURYCOUTER_
#define _QURYCOUTER_

#include <windows.h>//QueryCounter、及びSleep()用
#include <sstream> //ods()用

//デバッグウインドウへ出力
//template < typename T, typename T2 > void ods( T tep, T2 tep2, T tep3)
//{
//std::wstringstream ss;
//std::wstring st;
//
//ss << tep << tep2<< tep3 << L"\n";
//st = ss.str();
//OutputDebugString(st.c_str());
//};

//ミリ秒単位でStart()からEnd()までの時間を取得するクラス
//表示される値の小数点下四桁くらいまでなら正確
class QueryCounter
{
public:
QueryCounter(void)
{
LARGE_INTEGER freq;
QueryPerformanceFrequency(&freq);
dFreq_ = (double)freq.QuadPart / 1000.0;
}

void Start()
//時間計測の開始時に呼んで下さい
{
QueryPerformanceCounter(&start_);
}

double End()
//時間計測の終了時に呼んで下さい
{
LARGE_INTEGER end;
QueryPerformanceCounter(&end);
return (double (end.QuadPart - start_.QuadPart) / dFreq_ );
}

private:
LARGE_INTEGER start_;
double dFreq_;
};

#endif