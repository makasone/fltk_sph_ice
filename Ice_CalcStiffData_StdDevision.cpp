#include "Ice_CalcStiffData_StdDevision.h"

typedef Ice_CalcStiffData_StdDevision StiffDivision;

StiffDivision::Ice_CalcStiffData_StdDevision(const vector<Ice_SM*>& iceSM, Ice_JudgeMove* judge)
{
	m_iceSM = iceSM;
	m_iceJudge = judge;
}

StiffDivision::~Ice_CalcStiffData_StdDevision()
{

}

//�e�N���X�^�̕W���΍��̑��a��Ԃ�
//�኱�������d��
float StiffDivision::StepCalcData()
{
	unsigned pNum = IceObject::GetParticleNum();
	vector<Vec3> finalPos(pNum);
	//vector<Vec3> finalPos2(pNum);
	vector<unsigned> addParticleNum(pNum, 0);

	//�e���q�̕��ψʒu
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		//�ŏI���ʎZ�o�ɗp����N���X�^�̔���
		if(m_iceJudge->JudgeMove(cIndx) == false){	continue;	}

		//�N���X�^���܂�ł���e���q�̈ʒu�E���x��static�ȍŏI�ʒu�ɑ���
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			Vec3 pos = m_iceSM[cIndx]->GetVertexPos(oIndx);

			finalPos[pIndx] += pos;
			//finalPos2[pIndx] += pos*pos;

			//���q���̃J�E���g
			addParticleNum[pIndx] += 1;
		}
	}

	#pragma omp parallel for
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		finalPos[pIndx] /= addParticleNum[pIndx];
	}

	//�e���q�̈ʒu�̕W���΍�
	vector<float> stdDives(pNum, 0.0f);
	float stdDivSum = 0.0f;

	//�e���q�̕W���΍��i�����܂ł��ƕ�����������Ă��Ȃ��̂ŕ��U�j
	for(int cIndx = 0; cIndx < pNum; ++cIndx)
	{
		//�^���v�Z�����N���X�^�݂̂�Ώ�
		if(m_iceJudge->JudgeMoveDebug(cIndx) == false){	continue;	}

		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); ++oIndx)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			Vec3 pos = m_iceSM[cIndx]->GetVertexPos(oIndx);
			
			//���U�����߂�@�e���ő������킹
			stdDives[pIndx] += pow(pos[0] - finalPos[pIndx][0], 2.0);
			stdDives[pIndx] += pow(pos[1] - finalPos[pIndx][1], 2.0);
			stdDives[pIndx] += pow(pos[2] - finalPos[pIndx][2], 2.0);

			////�������ɂ���Ə��������₷���Ȃ�
			//stdDives[pIndx] += abs(pos[0] - finalPos[pIndx][0]);
			//stdDives[pIndx] += abs(pos[1] - finalPos[pIndx][1]);
			//stdDives[pIndx] += abs(pos[2] - finalPos[pIndx][2]);
		}
	}

	//�W���΍�
	#pragma omp parallel for
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		stdDives[pIndx] = stdDives[pIndx]/(float)addParticleNum[pIndx];
	}

	//�W���΍��̑��a
	#pragma omp parallel for reduction(+:stdDivSum)
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		stdDivSum += sqrt(stdDives[pIndx]);
		//stdDivSum += stdDives[pIndx];
	}

	return stdDivSum;

//��������������������
//#include <stdio.h>
//#include <math.h>
//
//void main() {
//    int i;
//    double data[] = { 3.5, 3.0, 2.8, 4.2 };
//    int n = sizeof(data)/(sizeof double); 
//    double avg, dev, sum = 0.0, sum2 = 0.0;
//
//    for (i = 0; i < n; i++) {
//        sum  += data[i];
//        sum2 += data[i]*data[i];
//    }
//    avg = sum/n;
//    dev = sqrt(sum2/n - avg*avg);
//    printf("����=%.3f, �W���΍�=%.3f", avg, dev);
}

void StiffDivision::StepCalcDataDebug()
{
}