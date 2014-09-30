
#include <Surf_SM.h>

//------------------------------------------------������----------------------------------------------
/*!
 * prefixSum���v�Z���邽�߂̃p�X������
 * @param[in]
 * @param[in]
 * @param[in]
 */
void Surf_SM::InitPath(const float* pos, const float* vel, const vector<Ice_SM*> iceSM, IceStructure* strct, int pthSize, int prtNum)
{	cout << __FUNCTION__ << endl;

	m_iceSM = iceSM;
	m_strct = strct;

	m_fPos = pos;
	m_fVel = vel;
	m_iPrtclNum = prtNum;

	//�܂��́C�����̂̕\�ʂŎ������߂ɓK���ȃp�X������Ă݂�D
	//�Ȃ������߂Ƃ���Ă���P�{�̃p�X�ł���Ă݂�D
	//�p�X�����q
	int pthNum = 1;								//�Ƃ肠����1�{
	cout << "check1" << endl;
	m_mk2DiPTHtoPRT.SetSize(pthNum, pthSize);

	for(int i = 0; i < pthNum; i++)
	{
		for(int j = 0; j < pthSize; j++)
		{
			m_mk2DiPTHtoPRT(i, j) = -1;		//-1�ŏ�����
		}
	}
	cout << "check2" << endl;
	for(int i = 0; i < prtNum; i++)
	{
		m_mk2DiPTHtoPRT(0, i) = i;			//���ۂ̃p�X�̏����� 1�{�p
	}
	cout << "check3" << endl;
	//���q���p�X
	m_mk2DiPRTtoPTH.SetSize(prtNum, 2);

	//�P�{�p
	for(int i = 0; i < prtNum; i++)
	{
		m_mk2DiPRTtoPTH(i, 0) = 0;
		m_mk2DiPRTtoPTH(i, 1) = i;
	}

	InitOrgPos(prtNum);							//�����ʒu

	//prefixSum��������
	m_mk2Dvec3_PrfxPos.SetSize(pthNum, prtNum);	//�ʒu
	UpdatePrefixSumPos();

	m_mk2Dmat3_PrfxApq.SetSize(pthNum, prtNum);	//�ό`�s��
	UpdatePrefixSumApq();

	InitPathPrfxIndxSet(iceSM, strct);			//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z

	InitOrgCm();								//�����d�S

//�f�o�b�O
	//DebugPathDataPos();
	//DebugPathDataApq();
	//DebugPathPrfxIndxSet();
}

//GPU
void Surf_SM::InitPathGPU()
{	cout << __FUNCTION__ << endl;
	
	//�������m��
	int pthNum = 1;
	int pthSize = m_iPrtclNum;

	cudaMalloc((void**)&md_2DiPTHtoPRT,	sizeof(int) * pthNum * pthSize);		//�p�X�����q
	cudaMalloc((void**)&md_2DiPRTtoPTH,	sizeof(int) * m_iPrtclNum * 2);			//���q���p�X

	md_f3OrgPos = Ice_SM::GetOrgPosPointer();									//�N���X�^���̗��q�̏����ʒu
	md_f3OrgCm	= Ice_SM::GetOrgCmPointer();									//�N���X�^�̏����d�S		

	cudaMalloc((void**)&md_2Df3PrfxPos,	sizeof(float) * pthNum * pthSize * 3);	//�ʒu��prefixSum
	cudaMalloc((void**)&md_2Df9PrfxApq,	sizeof(float) * pthNum * pthSize * 9);	//�ό`��prefixSum

	cudaMalloc((void**)&md_3DiPTHandPrfxSet,	sizeof(int) * m_iPrtclNum * m_iPrtclNum * 2);	//�ό`��prefixSum


	//������
	//vector���g���Ă���̂ŁC��������������Ƃ�₱���������邪�A�h���X��n���Ă��邾���D
	cudaMemcpy(md_2DiPTHtoPRT, &m_mk2DiPTHtoPRT.Get()[0], sizeof(int) * pthNum * pthSize, cudaMemcpyHostToDevice);
	cudaMemcpy(md_2DiPRTtoPTH, &m_mk2DiPRTtoPTH.Get()[0], sizeof(int) * m_iPrtclNum * 2 , cudaMemcpyHostToDevice);

	cudaMemcpy(md_2Df3PrfxPos, &m_mk2Dvec3_PrfxPos.Get()[0], sizeof(float) * pthNum * pthSize * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(md_2Df9PrfxApq, &m_mk2Dmat3_PrfxApq.Get()[0], sizeof(float) * pthNum * pthSize * 9, cudaMemcpyHostToDevice);

	cudaMemcpy(md_3DiPTHandPrfxSet, &m_mk3DiPTHandPrfxSet.Get()[0], sizeof(int) * m_iPrtclNum * m_iPrtclNum * 2, cudaMemcpyHostToDevice);


//�f�o�b�O
	////�f�[�^��]��
	//int* PTHtoPRT = new int[pthNum * pthSize];
	//int* PRTtoPTH = new int[m_iPrtclNum * 2];

	//cudaMemcpy(PTHtoPRT, md_2DiPTHtoPRT, sizeof(int) * pthNum * pthSize, cudaMemcpyDeviceToHost);
	//cudaMemcpy(PRTtoPTH, md_2DiPRTtoPTH, sizeof(int) * m_iPrtclNum * 2 , cudaMemcpyDeviceToHost);

	////�z�X�g���̃f�[�^��]���������ʂ��_���v
	//ofstream ofs( "Surf_SM.txt" );
	//ofs << __FUNCTION__ << endl;
	//
	////�f�o�C�X���̃f�[�^��]��
	//for(int i = 0; i < pthNum * pthSize; i++)
	//{
	//	ofs << i << " ";
	//	ofs << "md_2DiPTHtoPRT  " << PTHtoPRT[i] << "	" << "m_mk2DiPTHtoPRT " << m_mk2DiPTHtoPRT.Get()[i] << endl;
	//}

	//ofs << endl;

	//for(int i = 0; i < m_iPrtclNum * 2; i++)
	//{
	//	ofs << i << " ";
	//	ofs << "md_2DiPRTtoPTH  " << PRTtoPTH[i] << "	" << "m_mk2DiPRTtoPTH " << m_mk2DiPRTtoPTH.Get()[i] << endl;
	//}

	//delete[] PTHtoPRT;
	//delete[] PRTtoPTH;
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
	m_vvec3OrgCm.resize(m_iPrtclNum);

	for(int iclstr = 0; iclstr < m_iPrtclNum; iclstr++)
	{
		Vec3 vec(0.0);
		vec = CalcCmSum(iclstr);
		m_vvec3OrgCm[iclstr] = vec;
	}
}

//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z
void Surf_SM::InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct)
{
	//�T�C�Y�́C�N���X�^���~�ߖT���q��
	//����Ȃɂ�������͂���Ȃ��͂������C�Ìł��l�����Ĕ{���炢�ɂ��Ă����H
	//m_mk3DiPTHandPrfxSet.SetSize(m_iPrtclNum, m_iPrtclNum*2, 2);
	//�������������ꍇ�́C�Œ����
	m_mk3DiPTHandPrfxSet.SetSize(m_iPrtclNum, m_iPrtclNum, 2);

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
					m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 1) = kord;	//���݂���Ȃ�I�_���X�V
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

void Surf_SM::UpdatePrefixSumItr()
{
	//UpdatePrefixSumPosItr();
	//UpdatePrefixSumApqItr();

	//���Ɍ��ʂȂ�
	#pragma omp parallel
	#pragma omp sections
	{
	    #pragma omp section
	    {
	        UpdatePrefixSumPosItr();
	    }
	    #pragma omp section
	    {
	        UpdatePrefixSumApqItr();
	    }
	}
}

//GPU
void Surf_SM::UpdatePrefixSumGPU()
{
	int PosSizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int PosSizeY = m_mk2DiPTHtoPRT.GetSizeY();

	int ApqSizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int ApqSizeY = m_mk2DiPTHtoPRT.GetSizeY();

	LaunchUpdatePrefixSumGPU(
		m_iPrtclNum,
		PosSizeX,
		PosSizeY,
		ApqSizeX,
		ApqSizeY,
		md_2DiPTHtoPRT,
		md_2DiPRTtoPTH,
		md_2Df3PrfxPos,
		md_2Df9PrfxApq,
		md_3DiPTHandPrfxSet,
		md_f3OrgPos,
		md_f3OrgCm,
		md_fPos,
		md_fVel
		);
}

//�d�S��prefix sum
void Surf_SM::UpdatePrefixSumPos()
{//	cout << __FUNCTION__ << endl;

	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	//�p�X�̐��~�p�X�Ɋ܂܂��ő嗱�q���@�����J��Ԃ�
	for(int indxX = 0; indxX < sizeX; ++indxX)
	{
		Vec3 preVec(0.0, 0.0, 0.0);
		int pIndx = 0;

		for(int indxY = 0; indxY < sizeY; ++indxY)
		{
			//TODO:�����͕K��-1�ł��邱�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*4;
			if(pIndx < 0){	break;	}

			double mass = 1.0;		//�Ƃ肠����1.0�ŌŒ�

			m_mk2Dvec3_PrfxPos(indxX, indxY) = Vec3(m_fPos[pIndx+0], m_fPos[pIndx+1], m_fPos[pIndx+2]) * mass
											+ Vec3(m_fVel[pIndx+0], m_fVel[pIndx+1], m_fVel[pIndx+2]) * 0.02
											+ preVec;

			preVec = m_mk2Dvec3_PrfxPos(indxX, indxY);
		}
	}
}

//�d�S��prefix sum�@��������
void Surf_SM::UpdatePrefixSumPosItr()
{//	cout << __FUNCTION__ << endl;

	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	float* pos = Ice_SM::GetSldPosPointer();
	float* vel = Ice_SM::GetSldVelPointer();

	//�p�X�̐��~�p�X�Ɋ܂܂��ő嗱�q���@�����J��Ԃ�
	for(int indxX = 0; indxX < sizeX; ++indxX)
	{
		Vec3 preVec(0.0, 0.0, 0.0);
		int pIndx = 0;

		for(int indxY = 0; indxY < sizeY; ++indxY)
		{
			//TODO:�����͕K��-1�ł��邱�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*3;
			if(pIndx < 0){	break;	}

			double mass = 1.0;		//�Ƃ肠����1.0�ŌŒ�

			m_mk2Dvec3_PrfxPos(indxX, indxY) = Vec3(pos[pIndx+0], pos[pIndx+1], pos[pIndx+2]) * mass + preVec;

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
	for(int indxX = 0; indxX < sizeX; ++indxX)
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
			if(pIndx < 0){	break;	}
			double mass = 1.0;		//�Ƃ肠����1.0�ŌŒ�

			p = Vec3(m_fPos[pIndx+0], m_fPos[pIndx+1], m_fPos[pIndx+2])
				+ Vec3(m_fVel[pIndx+0], m_fVel[pIndx+1], m_fVel[pIndx+2]) * 0.02;

			q = m_vvec3OrgPos[pIndx/4];

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

//�ό`�s��imoment matrix�jApq��prefix sum ��������
void Surf_SM::UpdatePrefixSumApqItr()
{//	cout << __FUNCTION__ << endl;
	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	float* pos = Ice_SM::GetSldPosPointer();
	float* vel = Ice_SM::GetSldVelPointer();

	//�p�X�̐��~�p�X�Ɋ܂܂��ő嗱�q���@�����J��Ԃ�
	for(int indxX = 0; indxX < sizeX; ++indxX)
	{
		rxMatrix3 preMat(0.0);
		int pIndx = 0;
		Vec3 p;
		Vec3 q;

		for(int indxY = 0; indxY < sizeY; indxY++)
		{
			//TODO:�����͕K��-1�ł��邱�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*3;
			if(pIndx < 0){	break;	}
			double mass = 1.0;		//�Ƃ肠����1.0�ŌŒ�

			p = Vec3(pos[pIndx+0], pos[pIndx+1], pos[pIndx+2]);

			q = m_vvec3OrgPos[pIndx/3];

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
Vec3 Surf_SM::CalcCmSum(const int& cIndx)
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

		//�֐����g���o�[�W����
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

	return cmSum;
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

rxMatrix3 Surf_SM::CalcApqSum(const int& cIndx)
{
	rxMatrix3 ApqSum(0.0);
	rxMatrix3 mtt0T(0.0);			//M_i t_i (t^0_i)^T
	double mass = 1.0;
	Vec3 t(0.0);
	t = CalcCmSum(cIndx);			//TODO: ������v�Z�����ɓǂݍ��߂΂����Ƒ����Ȃ�

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

	return ApqSum;
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