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

//各クラスタの標準偏差の総和を返す
//若干処理が重い
float StiffDivision::StepCalcData()
{
	unsigned pNum = IceObject::GetParticleNum();
	vector<Vec3> finalPos(pNum);
	//vector<Vec3> finalPos2(pNum);
	vector<unsigned> addParticleNum(pNum, 0);

	//各粒子の平均位置
	for(int cIndx = 0; cIndx < pNum; cIndx++)
	{
		//最終結果算出に用いるクラスタの判定
		if(m_iceJudge->JudgeMove(cIndx) == false){	continue;	}

		//クラスタが含んでいる各粒子の位置・速度をstaticな最終位置に足す
		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); oIndx++)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			Vec3 pos = m_iceSM[cIndx]->GetVertexPos(oIndx);

			finalPos[pIndx] += pos;
			//finalPos2[pIndx] += pos*pos;

			//粒子数のカウント
			addParticleNum[pIndx] += 1;
		}
	}

	#pragma omp parallel for
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		finalPos[pIndx] /= addParticleNum[pIndx];
	}

	//各粒子の位置の標準偏差
	vector<float> stdDives(pNum, 0.0f);
	float stdDivSum = 0.0f;

	//各粒子の標準偏差（ここまでだと平方根を取っていないので分散）
	for(int cIndx = 0; cIndx < pNum; ++cIndx)
	{
		//運動計算したクラスタのみを対象
		if(m_iceJudge->JudgeMoveDebug(cIndx) == false){	continue;	}

		for(int oIndx = 0; oIndx < m_iceSM[cIndx]->GetIndxNum(); ++oIndx)
		{
			int pIndx = m_iceSM[cIndx]->GetParticleIndx(oIndx);
			if(pIndx == MAXINT){	continue;	}

			Vec3 pos = m_iceSM[cIndx]->GetVertexPos(oIndx);
			
			//分散を求める　各軸で足し合わせ
			stdDives[pIndx] += pow(pos[0] - finalPos[pIndx][0], 2.0);
			stdDives[pIndx] += pow(pos[1] - finalPos[pIndx][1], 2.0);
			stdDives[pIndx] += pow(pos[2] - finalPos[pIndx][2], 2.0);

			////こっちにすると少し扱いやすくなる
			//stdDives[pIndx] += abs(pos[0] - finalPos[pIndx][0]);
			//stdDives[pIndx] += abs(pos[1] - finalPos[pIndx][1]);
			//stdDives[pIndx] += abs(pos[2] - finalPos[pIndx][2]);
		}
	}

	//標準偏差
	#pragma omp parallel for
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		stdDives[pIndx] = stdDives[pIndx]/(float)addParticleNum[pIndx];
	}

	//標準偏差の総和
	#pragma omp parallel for reduction(+:stdDivSum)
	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		stdDivSum += sqrt(stdDives[pIndx]);
		//stdDivSum += stdDives[pIndx];
	}

	return stdDivSum;

//こういう書き方もある
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
//    printf("平均=%.3f, 標準偏差=%.3f", avg, dev);
}

void StiffDivision::StepCalcDataDebug()
{
}