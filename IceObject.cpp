#include "IceObject.h"

float* IceObject::s_sphPrtPos;
float* IceObject::s_sphPrtVel;

float* IceObject::m_fInterPolationCoefficience;

//デバイスポインタ
float* IceObject::sd_sldPrtPos;	
float* IceObject::sd_sldPrtVel;

float* IceObject::sd_ObjPrtPos;
float* IceObject::sd_ObjPrtVel;

int IceObject::sm_particleNum;
int IceObject::sm_tetraNum;
int IceObject::sm_clusterNum;


IceObject::IceObject(float* pos, float* vel, int pMaxNum, int cMaxNum, int tMaxNum)
{
	s_sphPrtPos = pos;
	s_sphPrtVel = vel;

	InitIceObj(pMaxNum, cMaxNum, tMaxNum);

}

IceObject::~IceObject()
{
}

//それぞれのクラス・変数の初期化
void IceObject::InitIceObj(int pMaxNum, int cMaxNum, int tMaxNum)
{	cout << __FUNCTION__ << endl;
	//物体の構造の初期化
	m_iceStrct = new IceStructure(pMaxNum, cMaxNum, tMaxNum);

	//運動計算を行うSMクラスタの初期化

	//補間処理のためのパラメータの初期化
	InitInterPolation();
}

//補間処理のためのパラメータの初期化
void IceObject::InitInterPolation()
{
	m_fInterPolationCoefficience = new float[sm_particleNum];	//線形補間係数
	cout << __FUNCTION__ << " pnum = " << sm_particleNum << endl;

	for(int i = 0; i < sm_particleNum; ++i)
	{
		m_fInterPolationCoefficience[i] = 1.0f;
	}
}

//固体の運動計算
void IceObject::StepObjMove()
{
	//GPUを用いたクラスタの運動計算
	//Ice_SM::UpdateGPU();

	//TODO::ここをGPUで処理する．
	for(int i = 0; i < sm_particleNum; i++)
	{
		if(GetPtoCNum(i) == 0){	continue;	}
		m_iceMove[i]->CopyDeviceToInstance(i);
	}
}

//流体と固体の最終的な運動計算
void IceObject::StepInterPolation()
{
	Vec3 pos,vel;

	for(int i = 0; i < sm_particleNum; ++i)
	{
		if(GetParticleNum() <= i){	continue;	}	//融解のみの実験のときに必要になる．
		if(GetPtoCNum(i) <= 0){		continue;	}

		//固体運動の最終位置計算
		CalcAverageCPU(i, pos, vel);

		//液体と固体の補間
		LinerInterPolationCPU(i, pos, vel);
	}
}

void IceObject::CalcAverageCPU(const int pIndx, Vec3& pos, Vec3& vel)
{
		//それぞれのベクトルを合成し平均をとる
		pos = Vec3(0.0, 0.0, 0.0);
		vel = Vec3(0.0, 0.0, 0.0);
		double shapeNum = 0.0;		//クラスタの数

		for(int j = 0; j < GetPtoCIndx(pIndx); ++j)
		{
			int jcIndx = GetPtoC(pIndx, j, 0);
			int joIndx = GetPtoC(pIndx, j, 1);

			if(jcIndx == -1 || joIndx == -1){	continue;	}

			pos += m_iceMove[jcIndx]->GetVertexPos(joIndx);
			vel += m_iceMove[jcIndx]->GetVertexVel(joIndx);

			shapeNum += 1.0;
		}

		int jpIndx = pIndx*4;
		

		//クラスタの数で割る
		//どのクラスタにも含まれていない場合，運動はSPH法に従う
		if(shapeNum != 0.0)
		{
			pos /= shapeNum;
			vel /= shapeNum;
		}		
		else
		{
			pos = Vec3(s_sphPrtPos[jpIndx+0], s_sphPrtPos[jpIndx+1], s_sphPrtPos[jpIndx+2]);
			vel = Vec3(s_sphPrtVel[jpIndx+0], s_sphPrtVel[jpIndx+1], s_sphPrtVel[jpIndx+2]);
		}
}

//SPH法とSM法で求めた速度と位置を線形補間 CPU
void IceObject::LinerInterPolationCPU(const int pIndx, const Vec3& pos, const Vec3& vel)
{
	int sphIndx = pIndx*4;
	double intrps = 1.0-m_fInterPolationCoefficience[pIndx];	//補間係数

	s_sphPrtVel[sphIndx+0] = vel[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+0] * intrps;
	s_sphPrtVel[sphIndx+1] = vel[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+1] * intrps;
	s_sphPrtVel[sphIndx+2] = vel[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtVel[sphIndx+2] * intrps;

	s_sphPrtPos[sphIndx+0] = pos[0] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+0] * intrps;
	s_sphPrtPos[sphIndx+1] = pos[1] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+1] * intrps;
	s_sphPrtPos[sphIndx+2] = pos[2] * m_fInterPolationCoefficience[pIndx] + s_sphPrtPos[sphIndx+2] * intrps;
}