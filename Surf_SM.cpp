
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
 * prefixSum���v�Z����p�X�쐬
 * @param[in]
 * @param[in]
 * @param[in]
 */
void Surf_SM::MakePath(const float* pos, int particleNum, int pathSize)
{	cout << __FUNCTION__ << endl;
	//�܂��́C�����̂̕\�ʂŎ������߂ɓK���ȃp�X������Ă݂�D
	//�Ȃ������߂Ƃ���Ă���P�{�̃p�X�ł���Ă݂�D
	//
	m_mk2DiPTHoPRT.SetSize(1, particleNum);

	//�p�X�����q
	for(int i = 0; i < particleNum; i++)
	{
		m_mk2DiPTHoPRT(0, i) = i;
	}

	//���q���p�X
	m_viPRTtoPTH.resize(particleNum);

	for(int i = 0; i < particleNum; i++)
	{
		m_viPRTtoPTH[i] = 0;
	}

	//�ʒu
	

	//�s��

}

void Surf_SM::CalcPrefixSum()
{	cout << __FUNCTION__ << endl;

}