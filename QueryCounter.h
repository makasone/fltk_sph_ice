#ifndef _QURYCOUTER_
#define _QURYCOUTER_

#include <windows.h>//QueryCounter�A�y��Sleep()�p
#include <sstream> //ods()�p

//�f�o�b�O�E�C���h�E�֏o��
//template < typename T, typename T2 > void ods( T tep, T2 tep2, T tep3)
//{
//std::wstringstream ss;
//std::wstring st;
//
//ss << tep << tep2<< tep3 << L"\n";
//st = ss.str();
//OutputDebugString(st.c_str());
//};

//�~���b�P�ʂ�Start()����End()�܂ł̎��Ԃ��擾����N���X
//�\�������l�̏����_���l�����炢�܂łȂ琳�m
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
//���Ԍv���̊J�n���ɌĂ�ŉ�����
{
QueryPerformanceCounter(&start_);
}

double End()
//���Ԍv���̏I�����ɌĂ�ŉ�����
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