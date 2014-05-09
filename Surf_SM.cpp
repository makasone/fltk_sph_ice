
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
 * prefixSum���v�Z���邽�߂̃p�X�쐬
 * @param[in]
 * @param[in]
 * @param[in]
 */
void Surf_SM::MakePath(const float* pos, int prtNum, int pthSize)
{	cout << __FUNCTION__ << endl;
	//�܂��́C�����̂̕\�ʂŎ������߂ɓK���ȃp�X������Ă݂�D
	//�Ȃ������߂Ƃ���Ă���P�{�̃p�X�ł���Ă݂�D
	//

	//�p�X�����q
	int pthNum = 1;							//�Ƃ肠����1�{
	m_mk2DiPTHtoPRT.SetSize(pthNum, pthSize);

	for(int i = 0; i < pthNum; i++)
	{
		for(int j = 0; j < pthSize; j++)
		{
			m_mk2DiPTHtoPRT(i, j) = -1;		//-1�ŏ�����
		}
	}

	for(int i = 0; i < prtNum; i++)
	{
		m_mk2DiPTHtoPRT(0, i) = i;			//���ۂ̃p�X�̏����� 1�{�p
	}

	//���q���p�X
	m_mk2DiPRTtoPTH.SetSize(prtNum, 2);

	//�P�{�p
	for(int i = 0; i < prtNum; i++)
	{
		m_mk2DiPRTtoPTH(i, 0) = 0;
		m_mk2DiPRTtoPTH(i, 1) = i;
	}

	//prefixSum��������
	m_mk2Dvec3_PrfxPos.SetSize(pthNum, prtNum);	//�ʒu
	CalcPrefixSumPos(pos);
	
	m_mk2Dmat3_PrfxApq.SetSize(pthNum, prtNum);	//�ό`�s��
	InitOrgPos(pos, prtNum);					//�����ʒu
	CalcPrefixSumApq(pos);

	////�f�o�b�O
	//DebugPathDataPos();
	//DebugPathDataApq();
}

//�d�S��prefix sum
//�ȈՉ��̂��߂Ɏ��ʂ����Ƃ��Ă���
void Surf_SM::CalcPrefixSumPos(const float* p)
{	cout << __FUNCTION__ << endl;
	
	//�p�X�̐��~�p�X�Ɋ܂܂��ő嗱�q���@�����J��Ԃ�
	for(int indxX = 0; indxX < m_mk2DiPTHtoPRT.GetSizeX(); indxX++)
	{
		Vec3 preVec = Vec3(0.0f, 0.0f, 0.0f);
		float massSum = 0.0f;

		for(int indxY = 0; indxY < m_mk2DiPTHtoPRT.GetSizeY(); indxY++)
		{
			//TODO:�����͕K��-1�ł��邱�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			int pIndx = m_mk2DiPTHtoPRT(indxX, indxY);
			float mass = 1.0f;		//�Ƃ肠����1.0f�ŌŒ�

			if(pIndx == -1)
			{
				break;
			}

			massSum += mass;

			m_mk2Dvec3_PrfxPos(indxX, indxY) = Vec3(p[pIndx*4+0], p[pIndx*4+1], p[pIndx*4+2]) * mass + preVec;
			//m_mk2Dvec3_PrfxPos(indxX, indxY) = m_mk2Dvec3_PrfxPos(indxX, indxY) / massSum;
			preVec = m_mk2Dvec3_PrfxPos(indxX, indxY);
		}
	}
}

//�����ʒu�̏�����
void Surf_SM::InitOrgPos(const float* p, int pNum)
{
	m_vvec3OrgPos.resize(pNum);

	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		m_vvec3OrgPos[pIndx] = Vec3(p[pIndx*4+0], p[pIndx*4+1], p[pIndx*4+2]);
	}
}

//�ό`�s��imoment matrix�jApq��prefix sum
void Surf_SM::CalcPrefixSumApq(const float *pos)
{	cout << __FUNCTION__ << endl;

	//�p�X�̐��~�p�X�Ɋ܂܂��ő嗱�q���@�����J��Ԃ�
	for(int indxX = 0; indxX < m_mk2DiPTHtoPRT.GetSizeX(); indxX++)
	{
		rxMatrix3 preMat(0.0f);
		float massSum = 0.0f;

		for(int indxY = 0; indxY < m_mk2DiPTHtoPRT.GetSizeY(); indxY++)
		{
			//TODO:�����͕K��-1�ł��邱�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			int pIndx = m_mk2DiPTHtoPRT(indxX, indxY);
			float mass = 1.0f;		//�Ƃ肠����1.0f�ŌŒ�

			if(pIndx == -1)
			{
				break;
			}

			massSum += mass;

			Vec3 p = Vec3(pos[pIndx*4+0], pos[pIndx*4+1], pos[pIndx*4+2]);
			Vec3 q = m_vvec3OrgPos[pIndx];

			//���݂�Aij
			m_mk2Dmat3_PrfxApq(indxX, indxY)(0,0) = mass * p[0] * q[0];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(0,1) = mass * p[0] * q[1];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(0,2) = mass * p[0] * q[2];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(1,0) = mass * p[1] * q[0];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(1,1) = mass * p[1] * q[1];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(1,2) = mass * p[1] * q[2];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(2,0) = mass * p[2] * q[0];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(2,1) = mass * p[2] * q[1];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(2,2) = mass * p[2] * q[2];

			m_mk2Dmat3_PrfxApq(indxX, indxY) += preMat;		//����܂ł�Aij�����Z
			preMat = m_mk2Dmat3_PrfxApq(indxX, indxY);
		}
	}
}

void Surf_SM::DebugPathDataPos()
{	cout << __FUNCTION__ << endl;
	for(int indxX = 0; indxX < m_mk2Dvec3_PrfxPos.GetSizeX(); indxX++)
	{
		for(int indxY = 0; indxY < m_mk2Dvec3_PrfxPos.GetSizeY(); indxY++)
		{
			cout << "m_mk2Dvec3_PrfxPos(" << indxX << ", " << indxY << ") = " << m_mk2Dvec3_PrfxPos(indxX, indxY) << endl;
		}
	}
}

void Surf_SM::DebugPathDataApq()
{	cout << __FUNCTION__ << endl;
	for(int indxX = 0; indxX < m_mk2Dmat3_PrfxApq.GetSizeX(); indxX++)
	{
		for(int indxY = 0; indxY < m_mk2Dmat3_PrfxApq.GetSizeY(); indxY++)
		{
			cout << "m_mk2Dmat3_PrfxApq(" << indxX << ", " << indxY << ") = " << m_mk2Dmat3_PrfxApq(indxX, indxY) << endl;
		}
	}
}