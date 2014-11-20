
#include <Surf_SM.h>

//------------------------------------------------������----------------------------------------------
Surf_SM::Surf_SM(const float* pos, const float* vel, const vector<Ice_SM*>& iceSM, IceStructure* strct, int pthSize, int prtNum)
{
	InitPath(pos, vel, iceSM, strct, pthSize, prtNum);
}

/*!
 * prefixSum���v�Z���邽�߂̃p�X������
 * @param[in]
 * @param[in]
 * @param[in]
 */
void Surf_SM::InitPath(const float* pos, const float* vel, const vector<Ice_SM*>& iceSM, IceStructure* strct, int pthSize, int prtNum)
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

	InitOrgPos(prtNum, iceSM);							//�����ʒu

	//prefixSum��������
	m_mk2Dvec3_PrfxPos.SetSize(pthNum, prtNum);	//�ʒu
	UpdatePrefixSumPos();

	cout << "check4" << endl;
	m_mk2Dmat3_PrfxApq.SetSize(pthNum, prtNum);	//�ό`�s��
	UpdatePrefixSumApq();

	cout << "check5" << endl;
	InitPathPrfxIndxSet(iceSM, strct);			//�ǂ̃p�X�̂ǂ̕������K�v�Ȃ̂��C���N���X�^���ƂɌv�Z

	InitOrgCm();								//�����d�S
	cout << "end" << endl;

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

	md_f3ClusterOrgPos	= Ice_SM::GetOrgPosPointer();							//�N���X�^���̗��q�̏����ʒu
	md_f3ClusterOrgCm	= Ice_SM::GetOrgCmPointer();							//�N���X�^�̏����d�S		

	cudaMalloc((void**)&md_2DiPTHtoPRT,	sizeof(int) * pthNum * pthSize);		//�p�X�����q
	cudaMalloc((void**)&md_2DiPRTtoPTH,	sizeof(short int) * m_iPrtclNum * 2);			//���q���p�X

	cudaMalloc((void**)&md_f3OrgPos,	sizeof(float) * m_iPrtclNum * 3);		//���q�̏����ʒu

	cudaMalloc((void**)&md_2Df3PrfxPos,	sizeof(float) * pthNum * pthSize * 3);	//�ʒu��prefixSum
	cudaMalloc((void**)&md_2Df9PrfxApq,	sizeof(float) * pthNum * pthSize * 9);	//�ό`��prefixSum

	cudaMalloc((void**)&md_3DiPTHandPrfxSet, sizeof(short int) * m_iPrtclNum * m_iPrtclNum * 2);	//�e�N���X�^�ɂ�����C�p�X��prefixSum�̔Ԓn�Z�b�g

	//Group Pos
	unsigned* groupPos = new unsigned[m_iPrtclNum*3];
	for(int i = 0; i < m_iPrtclNum; i++)
	{
		groupPos[m_iPrtclNum*0+i] = 0;
		groupPos[m_iPrtclNum*1+i] = 1;
		groupPos[m_iPrtclNum*2+i] = 2;
	}

	cudaMalloc((void**)&md_uPosGroup,	sizeof(unsigned) * m_iPrtclNum * 3);
	cudaMemcpy(md_uPosGroup, groupPos, sizeof(unsigned) * m_iPrtclNum * 3, cudaMemcpyHostToDevice);

	//Group Apq
	unsigned* groupApq = new unsigned[m_iPrtclNum*9];
	for(int i = 0; i < m_iPrtclNum; i++)
	{
		for(int j = 0; j < 9; j++)
		{
			groupApq[m_iPrtclNum*j+i] = j;
		}
	}
	cudaMalloc((void**)&md_uApqGroup,	sizeof(unsigned) * m_iPrtclNum * 9);
	cudaMemcpy(md_uApqGroup, groupApq, sizeof(unsigned) * m_iPrtclNum * 9, cudaMemcpyHostToDevice);

	//vector���g���Ă���̂ŁC��������������Ƃ�₱���������邪�A�h���X��n���Ă��邾���D
	//�p�X�����q
	cudaMemcpy(md_2DiPTHtoPRT, &m_mk2DiPTHtoPRT.Get()[0], sizeof(int) * pthNum * pthSize, cudaMemcpyHostToDevice);
	//���q���p�X
	cudaMemcpy(md_2DiPRTtoPTH, &m_mk2DiPRTtoPTH.Get()[0], sizeof(short int) * m_iPrtclNum * 2 , cudaMemcpyHostToDevice);

	//���q�̏����ʒu
	float* orgPos = new float[m_iPrtclNum * 3];
	for(int i = 0; i < m_iPrtclNum; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			int pIndx  = i*4+j;
			int opIndx = i*3+j;
			orgPos[opIndx] = m_fPos[pIndx];
		}
	}

	cudaMemcpy(md_f3OrgPos, orgPos, sizeof(float) * m_iPrtclNum * 3, cudaMemcpyHostToDevice);

	//�ʒu��prefixSum
	//TODO: �ꎟ���Ɍ��肵�Ă���̂ɒ���
	//m_mk2Dvec2_PrfxPos���g���Ƃł��Ȃ��݂����Ȃ̂ŁC��xfloat*�ɕϊ�
	float* prfxPos = new float[pthNum * pthSize * 3];
	for(int i = 0; i < m_iPrtclNum; i++)
	{
		Vec3 preVec = m_mk2Dvec3_PrfxPos(0, i);
		prfxPos[i*3+0] = preVec[0];		prfxPos[i*3+1] = preVec[1];		prfxPos[i*3+2] = preVec[2];
	}

	//cudaMemcpy(md_2Df3PrfxPos, &m_mk2Dvec3_PrfxPos.Get()[0], sizeof(float) * pthNum * pthSize * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(md_2Df3PrfxPos, prfxPos, sizeof(float) * pthNum * pthSize * 3, cudaMemcpyHostToDevice);
	
	//�ό`��prefixSum
	float* prfxApq = new float[pthNum * pthSize * 9];
	for(int i = 0; i < m_iPrtclNum; i++)
	{
		rxMatrix3 preMatrix = m_mk2Dmat3_PrfxApq(0, i);
		prfxApq[i*9+0] = preMatrix(0, 0);	prfxApq[i*9+1] = preMatrix(0, 1);	prfxApq[i*9+2] = preMatrix(0, 2);
		prfxApq[i*9+3] = preMatrix(1, 0);	prfxApq[i*9+4] = preMatrix(1, 1);	prfxApq[i*9+5] = preMatrix(1, 2);
		prfxApq[i*9+6] = preMatrix(2, 0);	prfxApq[i*9+7] = preMatrix(2, 1);	prfxApq[i*9+8] = preMatrix(2, 2);
	}
	
	//cudaMemcpy(md_2Df9PrfxApq, &m_mk2Dmat3_PrfxApq.Get()[0], sizeof(float) * pthNum * pthSize * 9, cudaMemcpyHostToDevice);
	cudaMemcpy(md_2Df9PrfxApq, prfxApq, sizeof(float) * pthNum * pthSize * 9, cudaMemcpyHostToDevice);

	//�e�N���X�^�ɂ�����C�p�X��prefixSum�̔Ԓn�Z�b�g *2�͊J�n�ʒu�ƏI���ʒu
	cudaMemcpy(md_3DiPTHandPrfxSet, &m_mk3DiPTHandPrfxSet.Get()[0], sizeof(short int) * m_iPrtclNum * m_iPrtclNum * 2, cudaMemcpyHostToDevice);

	delete[] orgPos;
	delete[] prfxPos;
	delete[] prfxApq;

//�f�o�b�O
	//DebugInit();
}


//�����ʒu�̏�����
void Surf_SM::InitOrgPos(int prtNum, const vector<Ice_SM*>& iceSM)
{
	m_vvec3OrgPos.resize(prtNum);

	for(int pIndx = 0; pIndx < prtNum; pIndx++)
	{
		//�N���X�^���痱�q��T�����ď����ʒu��o�^
		int cIndx = iceSM[pIndx]->SearchIndx(pIndx);
		if(cIndx == MAXINT){	cout << __FUNCTION__ << " Error!" << endl;	}

		m_vvec3OrgPos[pIndx] = iceSM[pIndx]->GetOrgPos(cIndx);
	}
}

//�����d�S�̏�����
void Surf_SM::InitOrgCm()
{
	m_vvec3OrgCm.resize(m_iPrtclNum);

	for(int iclstr = 0; iclstr < m_iPrtclNum; iclstr++)
	{
		m_vvec3OrgCm[iclstr] = CalcCmSum(iclstr);
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
				if(ksearchIndx != MAXINT)
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

				if(ksearchIndx != MAXINT)
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

//GPU prefixSum�̍X�V
void Surf_SM::UpdatePrefixSumGPU()
{
	int PosSizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int PosSizeY = m_mk2DiPTHtoPRT.GetSizeY();

	int ApqSizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int ApqSizeY = m_mk2DiPTHtoPRT.GetSizeY();

	//GPU�Ńf�[�^�X�V
	LaunchUpdatePrefixSumGPU
	(
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
		md_f3ClusterOrgCm,
		md_uPosGroup,
		md_uApqGroup,
		md_fPos,
		md_fVel
	);
	
//�e�X�g
	//TestUpdatePrefixSum();

//�f�o�b�O
}

void Surf_SM::TestUpdatePrefixSum()
{
	QueryCounter counter;
	QueryCounter counter2;	
	QueryCounter counter3;

counter.Start();

	//int dataNum = (1 << 13);					//�f�[�^���͂��̐����ŕύX�@TESTDATANUM������������
	int dataNum = m_iPrtclNum;
	float* inData = new float[dataNum*9];
	float* outData;
	//outData = new float[dataNum*2];

cout << "�e�X�g�J�n�F�v�f��" << dataNum << "��" << endl;

	for(int indx = 0; indx < dataNum*9; indx++)
	{
		inData[indx] = (float)1.1;
	}

	float* dInData;
	float* dOutData;

	double makeTime = counter.End()/1000.0;
cout << "�f�[�^�쐬�I���F	" << makeTime << endl;

counter2.Start();

	int PosSizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int PosSizeY = m_mk2DiPTHtoPRT.GetSizeY();

	int ApqSizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int ApqSizeY = m_mk2DiPTHtoPRT.GetSizeY();

	//GPU�Ńf�[�^�X�V
	LaunchUpdatePrefixSumGPU
	(
		dataNum,			//�e�X�g�f�[�^
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
		md_f3ClusterOrgCm,
		md_uPosGroup,
		md_uApqGroup,
		md_fPos,
		md_fVel
	);

	double endTime = counter2.End()/1000.0;
cout << "�v�Z�I���F		" << endTime << endl;

counter3.Start();

	//�v�Z���ʂ̊m�F
	//UpdatePrefixSumPos();
	float* debugPos = new float[dataNum*3];
	float* debugApq = new float[dataNum*9];

	cudaMemcpy(debugPos, md_2Df3PrfxPos, sizeof(float) * dataNum * 3, cudaMemcpyDeviceToHost);
	cudaMemcpy(debugApq, md_2Df9PrfxApq, sizeof(float) * dataNum * 9, cudaMemcpyDeviceToHost);

	//Pos
	for(int indx = 0; indx < dataNum; indx++)
	{
		cout << "GPU[" << indx << "] =	(" 
			<< debugPos[dataNum*0+indx] << ", " << debugPos[dataNum*1+indx] << ", " << debugPos[dataNum*2+indx] << ")"
		<< endl;

		Vec3 pos = m_mk2Dvec3_PrfxPos(0, indx);
		cout << "CPU[" << indx << "] =	(" 
			<< pos[0] << ", " << pos[1] << ", " << pos[2] << ")"
		<< endl;
	}

	////Apq
	//for(int indx = 0; indx < dataNum; indx++)
	//{
	//	cout << "GPU:md_2Df9PrfxApq[" << indx << "] = " << endl;
	//	cout <<	"(" << debugApq[indx*9+0] << ", " << debugApq[indx*9+1] << ", " << debugApq[indx*9+2] << ")" << endl;
	//	cout <<	"(" << debugApq[indx*9+3] << ", " << debugApq[indx*9+4] << ", " << debugApq[indx*9+5] << ")" << endl;
	//	cout <<	"(" << debugApq[indx*9+6] << ", " << debugApq[indx*9+7] << ", " << debugApq[indx*9+8] << ")" << endl;
	//	cout << endl;

	//	rxMatrix3 matrix = m_mk2Dmat3_PrfxApq(0, indx);
	//	cout << "CPU:m_mk2Dmat3_PrfxApq[" << indx << "] = " << endl;
	//	cout <<	"(" << matrix(0, 0) << ", " << matrix(0, 1) << ", " << matrix(0, 2) << ")" << endl;
	//	cout <<	"(" << matrix(1, 0) << ", " << matrix(1, 1) << ", " << matrix(1, 2) << ")" << endl;
	//	cout <<	"(" << matrix(2, 0) << ", " << matrix(2, 1) << ", " << matrix(2, 2) << ")" << endl;
	//	cout << endl;
	//}

	delete[] debugPos;
	delete[] debugApq;
	delete[] inData;

	double finishTime = counter3.End()/1000.0;
cout << "��n���I���F		" << finishTime << endl;
cout << "���v�F		" << makeTime+endTime+finishTime << endl;
}


//�d�S��prefix sum
void Surf_SM::UpdatePrefixSumPos()
{//	cout << __FUNCTION__ << endl;

	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	//�p�X�̐�(sizeX)�~�p�X�Ɋ܂܂��ő嗱�q��(sizeY)�@�����J��Ԃ�
	for(int indxX = 0; indxX < sizeX; ++indxX)
	{
		Vec3 preVec(0.0, 0.0, 0.0);
		int pIndx = 0;

		for(int indxY = 0; indxY < sizeY; ++indxY)
		{
			//TODO:�����͕K��-1�ł���C�Ƃ�z�肵�Ă���
			//TODO:�������̏ꍇ������̂ŁC�p�X�ɏ������闱�q���Ƃ����ۑ������ق�����������
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*4;
			if(pIndx < 0){	break;	}

			double mass = 1.0;		//�Ƃ肠����1.0�ŌŒ�

			m_mk2Dvec3_PrfxPos(indxX, indxY) = Vec3(m_fPos[pIndx+0], m_fPos[pIndx+1], m_fPos[pIndx+2]) * mass
											+ Vec3(m_fVel[pIndx+0], m_fVel[pIndx+1], m_fVel[pIndx+2]) * 0.01
											+ preVec;

			preVec = m_mk2Dvec3_PrfxPos(indxX, indxY);
		}
	}

//�f�o�b�O
	//DebugPathDataPos();
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
				+ Vec3(m_fVel[pIndx+0], m_fVel[pIndx+1], m_fVel[pIndx+2]) * 0.01;

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
		//�N���X�^�̓p�X�̈ꕔ���𕡐��g�ݍ��킹�č����C�ƍl����
		//�Ȃ̂ŁC�p�X�̂ǂ̕�������ǂ̕����܂ł𔲂��o���ďd�S�����߂�C�𕡐���s��
		start	= m_mk3DiPTHandPrfxSet(cIndx, iprt, 0);
		end		= m_mk3DiPTHandPrfxSet(cIndx, iprt, 1);

		if(start == -1 || end == -1){	break;	}
		
		prtIndx = m_strct->GetCtoP(cIndx, iprt, 0);
		pthIndx = m_mk2DiPRTtoPTH(prtIndx, 0);

		//cmSum += Vec3(ctopNum, prtIndx , pthIndx);

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
void Surf_SM::DebugInit()
{
	////�f�[�^��]��
	//int* PTHtoPRT = new int[pthNum * pthSize];
	//int* PRTtoPTH = new int[m_iPrtclNum * 2];

	//float* prefixPos = new float[pthNum * pthSize * 3];
	//float* prefixApq = new float[pthNum * pthSize * 9];

	//cudaMemcpy(PTHtoPRT, md_2DiPTHtoPRT, sizeof(int) * pthNum * pthSize, cudaMemcpyDeviceToHost);
	//cudaMemcpy(PRTtoPTH, md_2DiPRTtoPTH, sizeof(int) * m_iPrtclNum * 2 , cudaMemcpyDeviceToHost);

	//cudaMemcpy(prefixPos, md_2Df3PrfxPos, sizeof(float) * pthNum * pthSize * 3, cudaMemcpyDeviceToHost);
	//cudaMemcpy(prefixApq, md_2Df9PrfxApq, sizeof(float) * pthNum * pthSize * 9, cudaMemcpyDeviceToHost);

	////�z�X�g���̃f�[�^��]���������ʂ��_���v
	//ofstream ofs( "Surf_SM.txt" );
	//ofs << __FUNCTION__ << endl;
	//
	//�f�o�C�X���̃f�[�^��]��

	//PTHtoPRT
	//for(int i = 0; i < pthNum * pthSize; i++)
	//{
	//	ofs << i << " ";
	//	ofs << "md_2DiPTHtoPRT  " << PTHtoPRT[i] << "	" << "m_mk2DiPTHtoPRT " << m_mk2DiPTHtoPRT.Get()[i] << endl;
	//}

	//ofs << endl;

	//PRTtoPTH
	//for(int i = 0; i < m_iPrtclNum * 2; i++)
	//{
	//	ofs << i << " ";
	//	ofs << "md_2DiPRTtoPTH  " << PRTtoPTH[i] << "	" << "m_mk2DiPRTtoPTH " << m_mk2DiPRTtoPTH.Get()[i] << endl;
	//}

	//ofs << endl;

	////PrfxPos
	//for(int i = 0; i < m_iPrtclNum; i++)
	//{
	//	ofs << i << " ";
	//	ofs << "GPU = (" << prefixPos[i*3+0] << ", " << prefixPos[i*3+1] << ", " << prefixPos[i*3+2] << ")		";

	//	Vec3 preVec = m_mk2Dvec3_PrfxPos(0, i);
	//	ofs << "CPU = (" << preVec[0] << ", " << preVec[1] << ", " << preVec[2] << ")" << endl;
	//}

	////PrfxApq
	//for(int i = 0; i < m_iPrtclNum; i++)
	//{
	//	ofs << i << " ";
	//	ofs << "GPU =		(" << prefixApq[i*9+0] << ", " << prefixApq[i*9+1] << ", " << prefixApq[i*9+2] << ")" << endl;
	//	ofs << "		(" << prefixApq[i*9+3] << ", " << prefixApq[i*9+4] << ", " << prefixApq[i*9+5] << ")" << endl;
	//	ofs << "		(" << prefixApq[i*9+6] << ", " << prefixApq[i*9+7] << ", " << prefixApq[i*9+8] << ")" << endl;

	//	rxMatrix3 preMatrix = m_mk2Dmat3_PrfxApq(0, i);
	//	ofs << i << " ";
	//	ofs << "CPU =		(" << preMatrix(0, 0) << ", " << preMatrix(0, 1) << ", " << preMatrix(0, 2) << ")" << endl;
	//	ofs << "		(" << preMatrix(1, 0) << ", " << preMatrix(1, 1) << ", " << preMatrix(1, 2) << ")" << endl;
	//	ofs << "		(" << preMatrix(2, 0) << ", " << preMatrix(2, 1) << ", " << preMatrix(2, 2) << ")" << endl;
	//}

	//delete[] PTHtoPRT;
	//delete[] PRTtoPTH;

	//delete[] prefixPos;
	//delete[] prefixApq;
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
			cout << "m_mk2Dmat3_PrfxApq(" << indxX << ", " << indxY << ") = \n" 
				<< m_mk2Dmat3_PrfxApq(indxX, indxY) << endl;
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