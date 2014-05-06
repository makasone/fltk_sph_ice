
#include <Surf_SM.h>

//Vec3 Surf_SM::GetPos()
//{
//	//Vec3 A_
//}
//
//rxMatrix3 Surf_SM::GetApq()
//{
//}
/*!
 * prefixSumを計算するパス作成
 * @param[in]
 * @param[in]
 * @param[in]
 */
void Surf_SM::MakePath(const float* pos, int particleNum, int pathSize)
{	cout << __FUNCTION__ << endl;
	//まずは，立方体の表面で試すために適当なパスを作ってみる．
	//なぜかだめとされている１本のパスでやってみる．
	//
	m_mk2DiPTHoPRT.SetSize(1, particleNum);

	//パス→粒子
	for(int i = 0; i < particleNum; i++)
	{
		m_mk2DiPTHoPRT(0, i) = i;
	}

	//粒子→パス
	m_viPRTtoPTH.resize(particleNum);

	for(int i = 0; i < particleNum; i++)
	{
		m_viPRTtoPTH[i] = 0;
	}

	//位置
	

	//行列

}

void Surf_SM::CalcPrefixSum()
{	cout << __FUNCTION__ << endl;

}