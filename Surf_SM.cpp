
#include <Surf_SM.h>
#include <IceStructure.h>

//------------------------------------------------������----------------------------------------------
/*!
 * prefixSum���v�Z���邽�߂̃p�X������
 * @param[in]
 * @param[in]
 * @param[in]
 */
void Surf_SM::InitPath(const float* pos, const float* vel, const vector<Ice_SM*> iceSM, IceStructure* strct, int pthSize)
{	cout << __FUNCTION__ << endl;

	m_iceSM = iceSM;
	m_strct = strct;	//�����̂����ł��Ƃ̏������|�C���^�o�R�ŃA�N�Z�X���邽�߁C�x���Ȃ��Ă���

	m_fPos = pos;
	m_fVel = vel;

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

	InitOrgPos(prtNum);					//�����ʒu

	//prefixSum��������
	m_mk2Dvec3_PrfxPos.SetSize(pthNum, prtNum);	//�ʒu
	UpdatePrefixSumPos();

	m_mk2Dmat3_PrfxApq.SetSize(pthNum, prtNum);	//�ό`�s��
	UpdatePrefixSumApq();
	cout << __FUNCTION__ << " check2" << endl;
	//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z
	InitPathPrfxIndxSet(iceSM, strct);
	cout << __FUNCTION__ << " check4" << endl;
	InitOrgCm();						//�����d�S
	cout << __FUNCTION__ << " check3" << endl;
	////�f�o�b�O
	//DebugPathDataPos();
	//DebugPathDataApq();
	//DebugPathPrfxIndxSet();

}

//�����ʒu�̏�����
void Surf_SM::InitOrgPos(int prtNum)
{
	m_vvec3OrgPos.resize(prtNum);

	for(int pIndx = 0; pIndx < prtNum; pIndx++)
	{
		m_vvec3OrgPos[pIndx] = Vec3(m_fPos[pIndx*4+0], m_fPos[pIndx*4+1], m_fPos[pIndx*4+2]);
	}
}

//�����d�S�̏�����
void Surf_SM::InitOrgCm()
{
	m_vvec3OrgCm.resize(m_strct->GetClusterNum());

	for(int iclstr = 0; iclstr < m_strct->GetClusterNum(); iclstr++)
	{
		Vec3 vec(0.0);
		 CalcCmSum(iclstr, vec);
		m_vvec3OrgCm[iclstr] = vec;
	}
}

//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z
void Surf_SM::InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct)
{
	//�T�C�Y�́C�N���X�^���~�ߖT���q��
	//����Ȃɂ�������͂���Ȃ��͂������C�Ìł��l�����Ĕ{���炢�ɂ��Ă����H
	//m_mk3DiPTHandPrfxSet.SetSize( strct->GetClusterNum(), strct->GetClusterNum()*2, 2 );
	//�������������ꍇ�́C�Œ����
	m_mk3DiPTHandPrfxSet.SetSize( strct->GetClusterNum(), strct->GetClusterNum()/20, 2 );

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
		for(int jprt = 0; jprt < iceSM[iclstr]->GetNumVertices(); jprt++)
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
			for(int kord = jordIndx+1; kord < m_mk2DiPTHtoPRT.GetSizeY(); kord++)
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


void Surf_SM::UpdatePrefixSum()
{
	//UpdatePrefixSumPos();
	//UpdatePrefixSumApq();

	//���Ɍ��ʂȂ�
	#pragma omp parallel
	#pragma omp sections
	{
	    #pragma omp section
	    {
	        UpdatePrefixSumPos();
	    }
	    #pragma omp section
	    {
	        UpdatePrefixSumApq();
	    }
	}
}

//�d�S��prefix sum
//�ȈՉ��̂��߂Ɏ��ʂ����Ƃ��Ă���
void Surf_SM::UpdatePrefixSumPos()
{//	cout << __FUNCTION__ << endl;

	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	//�p�X�̐��~�p�X�Ɋ܂܂��ő嗱�q���@�����J��Ԃ�
	for(int indxX = 0; indxX < sizeX; indxX++)
	{
		Vec3 preVec(0.0, 0.0, 0.0);
		int pIndx = 0;

		for(int indxY = 0; indxY < sizeY; indxY++)
		{
			//TODO:�����͕K��-1�ł��邱�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*4;
			if(pIndx < 0)
			{
				break;
			}

			double mass = 1.0;		//�Ƃ肠����1.0�ŌŒ�

			m_mk2Dvec3_PrfxPos(indxX, indxY) = Vec3(m_fPos[pIndx+0], m_fPos[pIndx+1], m_fPos[pIndx+2]) * mass
				+ (Vec3(m_fVel[pIndx+0], m_fVel[pIndx+1], m_fVel[pIndx+2])/* + Vec3(0.0, -9.81, 0.0)*0.01*/ ) * 0.02
				+ preVec;
			preVec = m_mk2Dvec3_PrfxPos(indxX, indxY);
		}
	}
}

//�ό`�s��imoment matrix�jApq��prefix sum
void Surf_SM::UpdatePrefixSumApq()
{//	cout << __FUNCTION__ << endl;
	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	//�p�X�̐��~�p�X�Ɋ܂܂��ő嗱�q���@�����J��Ԃ�
	for(int indxX = 0; indxX < sizeX; indxX++)
	{
		rxMatrix3 preMat(0.0);
		int pIndx = 0;
		Vec3 p;
		Vec3 q;

		for(int indxY = 0; indxY < sizeY; indxY++)
		{
			//TODO:�����͕K��-1�ł��邱�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*4;
			if(pIndx < 0)
			{
				break;
			}
			double mass = 1.0;		//�Ƃ肠����1.0�ŌŒ�

			/*Vec3 */p = Vec3(m_fPos[pIndx+0], m_fPos[pIndx+1], m_fPos[pIndx+2])
					+ (Vec3(m_fVel[pIndx+0], m_fVel[pIndx+1], m_fVel[pIndx+2])/* + Vec3(0.0, -9.81, 0.0)*0.01*/ ) * 0.02;
			/*Vec3 */q = m_vvec3OrgPos[pIndx*1/4];

			//���݂�Aij
			//m_mk2Dmat3_PrfxApq(indxX, indxY)(0,0) = mass * p[0] * q[0];
			//m_mk2Dmat3_PrfxApq(indxX, indxY)(0,1) = mass * p[0] * q[1];
			//m_mk2Dmat3_PrfxApq(indxX, indxY)(0,2) = mass * p[0] * q[2];
			//m_mk2Dmat3_PrfxApq(indxX, indxY)(1,0) = mass * p[1] * q[0];
			//m_mk2Dmat3_PrfxApq(indxX, indxY)(1,1) = mass * p[1] * q[1];
			//m_mk2Dmat3_PrfxApq(indxX, indxY)(1,2) = mass * p[1] * q[2];
			//m_mk2Dmat3_PrfxApq(indxX, indxY)(2,0) = mass * p[2] * q[0];
			//m_mk2Dmat3_PrfxApq(indxX, indxY)(2,1) = mass * p[2] * q[1];
			//m_mk2Dmat3_PrfxApq(indxX, indxY)(2,2) = mass * p[2] * q[2];

			m_mk2Dmat3_PrfxApq(indxX, indxY)(0,0) = p[0] * q[0];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(0,1) = p[0] * q[1];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(0,2) = p[0] * q[2];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(1,0) = p[1] * q[0];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(1,1) = p[1] * q[1];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(1,2) = p[1] * q[2];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(2,0) = p[2] * q[0];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(2,1) = p[2] * q[1];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(2,2) = p[2] * q[2];

			m_mk2Dmat3_PrfxApq(indxX, indxY) += preMat;		//����܂ł�Aij�����Z
			preMat = m_mk2Dmat3_PrfxApq(indxX, indxY);
		}
	}
}
//�d�S�̑��a��Ԃ�
void Surf_SM::CalcCmSum(const int& cIndx, Vec3& vec)
{	//cout << __FUNCTION__ << " start" << endl;
	Vec3 cmSum(0.0);

	int start = -1;
	int end = -1;

	int prtIndx = -1;
	int pthIndx = -1;

	//�N���X�^�ɗp�ӂ����f�[�^�Z�b�g���g���ďd�S�����߂�
	for(int iprt = 0, ctopNum = m_strct->GetCtoPNum(cIndx); iprt < ctopNum; iprt++)
	{
		start	= m_mk3DiPTHandPrfxSet(cIndx, iprt, 0);
		end		= m_mk3DiPTHandPrfxSet(cIndx, iprt, 1);

		if(start == -1 || end == -1){	break;	}
		
		prtIndx = m_strct->GetCtoP(cIndx, iprt, 0);
		pthIndx = m_mk2DiPRTtoPTH(prtIndx, 0);

		//cmSum += CalcCmFromPrfxSm(pthIndx, start, end);

		//�֐����g��Ȃ��o�[�W����
		if(start == 0)
		{
			cmSum += m_mk2Dvec3_PrfxPos(pthIndx, end);
		}
		else
		{
			cmSum += m_mk2Dvec3_PrfxPos(pthIndx, end) - m_mk2Dvec3_PrfxPos(pthIndx, start-1);
		}
	}

	vec = cmSum;
}

//prfixSum����l��Ԃ�
const Vec3 Surf_SM::CalcCmFromPrfxSm(const int& path, const int& start, const int& end)
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

void Surf_SM::CalcApqSum(const int& cIndx, rxMatrix3& matrix)
{
	rxMatrix3 ApqSum(0.0);
	rxMatrix3 mtt0T(0.0);			//M_i t_i (t^0_i)^T
	double mass = 1.0;
	Vec3 t(0.0);
	CalcCmSum(cIndx, t);		//�ꎞ�I�u�W�F�N�g������Ȃ��悤��

	t *= 1.0 / (m_strct->GetCtoPNum(cIndx));

	Vec3 t0T( m_iceSM[cIndx]->GetOrgCm() );

	mtt0T(0,0) = mass * t[0] * t0T[0];
	mtt0T(0,1) = mass * t[0] * t0T[1];
	mtt0T(0,2) = mass * t[0] * t0T[2];
	mtt0T(1,0) = mass * t[1] * t0T[0];
	mtt0T(1,1) = mass * t[1] * t0T[1];
	mtt0T(1,2) = mass * t[1] * t0T[2];
	mtt0T(2,0) = mass * t[2] * t0T[0];
	mtt0T(2,1) = mass * t[2] * t0T[1];
	mtt0T(2,2) = mass * t[2] * t0T[2];

	ApqSum -= (m_strct->GetCtoPNum(cIndx)) * mtt0T;

	int start = -1;
	int end = -1;

	int prtIndx = -1;
	int pthIndx = -1;

	for(int iprt = 0, ctopNum = m_strct->GetCtoPNum(cIndx); iprt < ctopNum; iprt++)
	{
		start	= m_mk3DiPTHandPrfxSet(cIndx, iprt, 0);
		end		= m_mk3DiPTHandPrfxSet(cIndx, iprt, 1);

		if(start == -1 || end == -1){	break;	}

		prtIndx = m_strct->GetCtoP(cIndx, iprt, 0);
		pthIndx = m_mk2DiPRTtoPTH(prtIndx, 0);

		//ApqSum += CalcApqFromPrfxSm(pthIndx, start, end);
		//�֐����g��Ȃ��o�[�W����
		if(start == 0)
		{
			ApqSum += m_mk2Dmat3_PrfxApq(pthIndx, end);
		}
		else
		{
			ApqSum += m_mk2Dmat3_PrfxApq(pthIndx, end) - m_mk2Dmat3_PrfxApq(pthIndx, start-1);
		}
	}

	matrix = ApqSum;
}

const rxMatrix3 Surf_SM::CalcApqFromPrfxSm(const int& path, const int& start, const int& end)
{
	if(start == 0)
	{
		return m_mk2Dmat3_PrfxApq(path, end);
	}
	else
	{
		return m_mk2Dmat3_PrfxApq(path, end) - m_mk2Dmat3_PrfxApq(path, start-1);
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
			if( m_mk3DiPTHandPrfxSet(iX, iY, 0) == -1 || m_mk3DiPTHandPrfxSet(iX, iY, 1) == -1 ){	cout << iY-1 << " ";	break;	}
			
			//cout << "m_mk3DiPTHandPrfxSet(" << iX << ", " << iY << ", 0) = " << m_mk3DiPTHandPrfxSet(iX, iY, 0)
			//	 << ", (" << iX << ", " << iY << ", 1) = " << m_mk3DiPTHandPrfxSet(iX, iY, 1) << endl;
		}
	}
}