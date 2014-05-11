
#include <Surf_SM.h>
#include <IceStructure.h>

//Vec3 Surf_SM::GetPos()
//{
//	//Vec3 A_
//}
//
//rxMatrix3 Surf_SM::GetApq()
//{
//}

//------------------------------------------------������----------------------------------------------
/*!
 * prefixSum���v�Z���邽�߂̃p�X������
 * @param[in]
 * @param[in]
 * @param[in]
 */
void Surf_SM::InitPath(const float* pos, const float* vel, const vector<Ice_SM*> iceSM, IceStructure* strct, int pthSize)
{	cout << __FUNCTION__ << endl;

	m_strct = strct;

	//�܂��́C�����̂̕\�ʂŎ������߂ɓK���ȃp�X������Ă݂�D
	//�Ȃ������߂Ƃ���Ă���P�{�̃p�X�ł���Ă݂�D
	//

	//�p�X�����q
	int pthNum = 1;							//�Ƃ肠����1�{
	int prtNum = strct->GetParticleNum();

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

	InitOrgPos(pos, prtNum);					//�����ʒu

	//prefixSum��������
	m_mk2Dvec3_PrfxPos.SetSize(pthNum, prtNum);	//�ʒu
	UpdatePrefixSumPos(pos, vel);
	
	m_mk2Dmat3_PrfxApq.SetSize(pthNum, prtNum);	//�ό`�s��
	UpdatePrefixSumApq(pos);

	//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z
	InitPathPrfxIndxSet(iceSM, strct);

	////�f�o�b�O
	//DebugPathDataPos();
	//DebugPathDataApq();
	//DebugPathPrfxIndxSet();

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

//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z
void Surf_SM::InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct)
{
	//�T�C�Y�́C�N���X�^���~�ߖT���q��
	//����Ȃɂ�������͂���Ȃ��͂������C�Ìł��l�����Ĕ{���炢�ɂ��Ă����H�D
	m_mk3DiPTHandPrfxSet.SetSize( strct->GetClusterNum(), strct->GetClusterNum()*2, 2 );

	for(unsigned i = 0; i < m_mk3DiPTHandPrfxSet.Get().size(); i++)
	{
		m_mk3DiPTHandPrfxSet[i] = -1;
	}

	//�T���������ǂ����̃t���O�쐬
	vector<bool> searchFlag;
	searchFlag.resize(strct->GetCtoPMax());

	//�S�ẴN���X�^�ŒT��
	for(unsigned iclstr = 0; iclstr < iceSM.size(); iclstr++)
	{
		//�t���O�̏�����
		for(unsigned j = 0; j < searchFlag.size(); j++)
		{
			searchFlag[j] = false;
		}

		int isetIndx = 0;	//�Z�b�g�̂��߂̃C���f�b�N�X

		//�N���X�^���̑S�Ă̗��q�ŒT��
		for(unsigned jprt = 0; jprt < iceSM[iclstr]->GetNumVertices(); jprt++)
		{
			if(searchFlag[jprt] == true){  continue;	}
			searchFlag[jprt] = true;									//�T���ς݂Ȃ̂Ńt���O�I��

			unsigned jprtIndx = iceSM[iclstr]->GetParticleIndx(jprt);	//�͂��߂̗��q�ԍ��擾
			unsigned jpthIndx = m_mk2DiPRTtoPTH(jprtIndx, 0);			//���q���ǂ̃p�X�ɏ������Ă��邩�擾
			unsigned jordIndx = m_mk2DiPRTtoPTH(jprtIndx, 1);			//���q���p�X���ŉ��Ԗڂ����擾

			if(jordIndx == -1)
			{
				cout << "jordIndx == -1" << endl;
			}

			//�������̓p�X���ŉ��Ԗڂ��̏��
			m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 0) = jordIndx;			//�n�_
			m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 1) = jordIndx;			//�I�_

			//�n�_�T��
			for(int kord = jordIndx-1; -1 < kord; kord--)
			{
				unsigned kprtIndx = m_mk2DiPTHtoPRT(jpthIndx, kord);	//�O�̗��q
				int ksearchIndx = iceSM[iclstr]->SearchIndx(kprtIndx);	//�N���X�^���ł̏��Ԃ��擾
																		//iceSM����擾�ł����������C���͂ł��Ȃ��̂ɒ��Ӂ@�T�����Ȃ��Ƃ����Ȃ�
				if(ksearchIndx != -1)
				{
					m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 0) = kord;	//���݂���Ȃ�n�_���X�V
					searchFlag[ksearchIndx] = true;						//�T���ς݂Ȃ̂Ńt���O�I��
				}
				else
				{
					break;		//�N���X�^�ɑ��݂��Ȃ��Ȃ�T���I��
				}
			}

			//�I�_�T��
			for(unsigned kord = jordIndx+1; kord < m_mk2DiPTHtoPRT.GetSizeY(); kord++)
			{
				unsigned kprtIndx = m_mk2DiPTHtoPRT(jpthIndx, kord);	//���̗��q
				int ksearchIndx = iceSM[iclstr]->SearchIndx(kprtIndx);	//�N���X�^���ł̏��Ԃ��擾

				if(ksearchIndx != -1)
				{
					m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 1) = kord;		//���݂���Ȃ�I�_���X�V
					searchFlag[ksearchIndx] = true;						//�T���ς݂Ȃ̂Ńt���O�I��
				}
				else
				{
					break;		//�N���X�^�ɑ��݂��Ȃ��Ȃ�T���I��
				}
			}

			isetIndx++;
		}
	}
}


void Surf_SM::UpdatePrefixSum(const float* p, const float* v)
{
	UpdatePrefixSumPos(p, v);
	UpdatePrefixSumApq(p);
}

//�d�S��prefix sum
//�ȈՉ��̂��߂Ɏ��ʂ����Ƃ��Ă���
void Surf_SM::UpdatePrefixSumPos(const float* p, const float* v)
{//	cout << __FUNCTION__ << endl;
	
	//�p�X�̐��~�p�X�Ɋ܂܂��ő嗱�q���@�����J��Ԃ�
	for(int indxX = 0; indxX < m_mk2DiPTHtoPRT.GetSizeX(); indxX++)
	{
		Vec3 preVec = Vec3(0.0, 0.0, 0.0);

		for(int indxY = 0; indxY < m_mk2DiPTHtoPRT.GetSizeY(); indxY++)
		{
			//TODO:�����͕K��-1�ł��邱�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			int pIndx = m_mk2DiPTHtoPRT(indxX, indxY);
			double mass = 1.0;		//�Ƃ肠����1.0f�ŌŒ�

			if(pIndx == -1)
			{
				break;
			}

			m_mk2Dvec3_PrfxPos(indxX, indxY) = Vec3(p[pIndx*4+0], p[pIndx*4+1], p[pIndx*4+2]) * mass
				+ (Vec3(v[pIndx*4+0], v[pIndx*4+1], v[pIndx*4+2])/* + Vec3(0.0, -9.81, 0.0)*0.01*/ ) * 0.01
				+ preVec;
			preVec = m_mk2Dvec3_PrfxPos(indxX, indxY);
		}
	}
}

//�ό`�s��imoment matrix�jApq��prefix sum
void Surf_SM::UpdatePrefixSumApq(const float *pos)
{//	cout << __FUNCTION__ << endl;

	//�p�X�̐��~�p�X�Ɋ܂܂��ő嗱�q���@�����J��Ԃ�
	for(int indxX = 0; indxX < m_mk2DiPTHtoPRT.GetSizeX(); indxX++)
	{
		rxMatrix3 preMat(0.0);
		double massSum = 0.0;

		for(int indxY = 0; indxY < m_mk2DiPTHtoPRT.GetSizeY(); indxY++)
		{
			//TODO:�����͕K��-1�ł��邱�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			int pIndx = m_mk2DiPTHtoPRT(indxX, indxY);
			double mass = 1.0;		//�Ƃ肠����1.0f�ŌŒ�

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
//�d�S�̑��a��Ԃ�
Vec3 Surf_SM::ClacCmSum(int cIndx, const float* p)
{	//cout << __FUNCTION__ << " start" << endl;
	Vec3 cmSum(0.0);

	//�e�N���X�^���ɗp�ӂ����f�[�^�Z�b�g���g���ďd�S�����߂�
	for(int iprt = 0; iprt < m_strct->GetClusterNum(); iprt++)
	{
		int prtIndx = m_strct->GetCtoP(cIndx, iprt, 0);
		int pthIndx = m_mk2DiPRTtoPTH(prtIndx, 0);
		int start	= m_mk3DiPTHandPrfxSet(cIndx, iprt, 0);
		int end		= m_mk3DiPTHandPrfxSet(cIndx, iprt, 1);

		if(start == -1 || end == -1){	break;	}

		cmSum += CalcCmFromPrfxSm(pthIndx, start, end)/* + Vec3(p[start*4+0], p[start*4+1], p[start*4+2])*/;
	}

	return cmSum;
}

//prfixSum����l��Ԃ�
Vec3 Surf_SM::CalcCmFromPrfxSm(int path, int start, int end)
{	//cout << __FUNCTION__ << endl;
	if(start == 0)
	{
		return m_mk2Dvec3_PrfxPos(path, end);
	}
	else
	{
		return m_mk2Dvec3_PrfxPos(path, end) - m_mk2Dvec3_PrfxPos(path, start-1);
	}
}





//--------------------------------------------------------�f�o�b�O-------------------------------------------
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

void Surf_SM::DebugPathPrfxIndxSet()
{	cout << __FUNCTION__ << endl;
	for(int iX = 0; iX < m_mk3DiPTHandPrfxSet.GetSizeX(); iX++)
	{
		for(int iY = 0; iY < m_mk3DiPTHandPrfxSet.GetSizeY(); iY++)
		{
			if( m_mk3DiPTHandPrfxSet(iX, iY, 0) == -1 || m_mk3DiPTHandPrfxSet(iX, iY, 1) == -1 ){	break;	}
			
			cout << "m_mk3DiPTHandPrfxSet(" << iX << ", " << iY << ", 0) = " << m_mk3DiPTHandPrfxSet(iX, iY, 0)
				 << ", (" << iX << ", " << iY << ", 1) = " << m_mk3DiPTHandPrfxSet(iX, iY, 1) << endl;
		}
	}
}