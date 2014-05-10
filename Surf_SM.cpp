
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

//------------------------------------------------初期化----------------------------------------------
/*!
 * prefixSumを計算するためのパス初期化
 * @param[in]
 * @param[in]
 * @param[in]
 */
void Surf_SM::InitPath(const float* pos, const vector<Ice_SM*> iceSM, IceStructure* strct, int pthSize)
{	cout << __FUNCTION__ << endl;
	//まずは，立方体の表面で試すために適当なパスを作ってみる．
	//なぜかだめとされている１本のパスでやってみる．
	//

	//パス→粒子
	int pthNum = 1;							//とりあえず1本
	int prtNum = strct->GetParticleNum();

	m_mk2DiPTHtoPRT.SetSize(pthNum, pthSize);

	for(int i = 0; i < pthNum; i++)
	{
		for(int j = 0; j < pthSize; j++)
		{
			m_mk2DiPTHtoPRT(i, j) = -1;		//-1で初期化
		}
	}

	for(int i = 0; i < prtNum; i++)
	{
		m_mk2DiPTHtoPRT(0, i) = i;			//実際のパスの初期化 1本用
	}

	//粒子→パス
	m_mk2DiPRTtoPTH.SetSize(prtNum, 2);

	//１本用
	for(int i = 0; i < prtNum; i++)
	{
		m_mk2DiPRTtoPTH(i, 0) = 0;
		m_mk2DiPRTtoPTH(i, 1) = i;
	}

	InitOrgPos(pos, prtNum);					//初期位置

	//prefixSumを初期化
	m_mk2Dvec3_PrfxPos.SetSize(pthNum, prtNum);	//位置
	CalcPrefixSumPos(pos);
	
	m_mk2Dmat3_PrfxApq.SetSize(pthNum, prtNum);	//変形行列
	CalcPrefixSumApq(pos);

	//どのパスのどの部分が必要なのか，をクラスタごとに計算
	InitPathPrfxIndxSet(iceSM, strct);

	////デバッグ
	//DebugPathDataPos();
	//DebugPathDataApq();
	//DebugPathPrfxIndxSet();

}

//初期位置の初期化
void Surf_SM::InitOrgPos(const float* p, int pNum)
{
	m_vvec3OrgPos.resize(pNum);

	for(int pIndx = 0; pIndx < pNum; pIndx++)
	{
		m_vvec3OrgPos[pIndx] = Vec3(p[pIndx*4+0], p[pIndx*4+1], p[pIndx*4+2]);
	}
}

//どのパスのどの部分が必要なのか，をクラスタごとに計算
void Surf_SM::InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct)
{
	//サイズは，クラスタ数×近傍粒子数
	//そんなにたくさんはいらないはずだが，凝固を考慮して倍ぐらいにしておく？．
	m_mk3DiPTHandPrfxSet.SetSize( strct->GetClusterNum(), strct->GetClusterNum()*2, 2 );

	for(unsigned i = 0; i < m_mk3DiPTHandPrfxSet.Get().size(); i++)
	{
		m_mk3DiPTHandPrfxSet[i] = -1;
	}

	//探索したかどうかのフラグ作成
	vector<bool> searchFlag;
	searchFlag.resize(strct->GetCtoPMax());

	//全てのクラスタで探索
	for(unsigned iclstr = 0; iclstr < iceSM.size(); iclstr++)
	{
		//フラグの初期化
		for(unsigned j = 0; j < searchFlag.size(); j++)
		{
			searchFlag[j] = false;
		}

		int isetIndx = 0;	//セットのためのインデックス

		//クラスタ内の全ての粒子で探索
		for(unsigned jprt = 0; jprt < iceSM[iclstr]->GetNumVertices(); jprt++)
		{
			if(searchFlag[jprt] == true){  continue;	}
			searchFlag[jprt] = true;									//探索済みなのでフラグオン

			unsigned jprtIndx = iceSM[iclstr]->GetParticleIndx(jprt);	//はじめの粒子番号取得
			unsigned jpthIndx = m_mk2DiPRTtoPTH(jprtIndx, 0);			//粒子がどのパスに所属しているか取得
			unsigned jordIndx = m_mk2DiPRTtoPTH(jprtIndx, 1);			//粒子がパス内で何番目かを取得

			if(jordIndx == -1)
			{
				cout << "jordIndx == -1" << endl;
			}

			//初期化はパス内で何番目かの情報
			m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 0) = jordIndx;			//始点
			m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 1) = jordIndx;			//終点

			//始点探索
			for(int kord = jordIndx-1; -1 < kord; kord--)
			{
				unsigned kprtIndx = m_mk2DiPTHtoPRT(jpthIndx, kord);	//前の粒子
				int ksearchIndx = iceSM[iclstr]->SearchIndx(kprtIndx);	//クラスタ内での順番を取得
																		//iceSMから取得できそうだが，実はできないのに注意　探索しないといけない
				if(ksearchIndx != -1)
				{
					m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 0) = kord;	//存在するなら始点を更新
					searchFlag[ksearchIndx] = true;						//探索済みなのでフラグオン
				}
				else
				{
					break;		//クラスタに存在しないなら探索終了
				}
			}

			//終点探索
			for(unsigned kord = jordIndx+1; kord < m_mk2DiPTHtoPRT.GetSizeY(); kord++)
			{
				unsigned kprtIndx = m_mk2DiPTHtoPRT(jpthIndx, kord);	//次の粒子
				int ksearchIndx = iceSM[iclstr]->SearchIndx(kprtIndx);	//クラスタ内での順番を取得

				if(ksearchIndx != -1)
				{
					m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 1) = kord;		//存在するなら終点を更新
					searchFlag[ksearchIndx] = true;						//探索済みなのでフラグオン
				}
				else
				{
					break;		//クラスタに存在しないなら探索終了
				}
			}

			isetIndx++;
		}
	}
}







//重心のprefix sum
//簡易化のために質量が一定としている
void Surf_SM::CalcPrefixSumPos(const float* p)
{	cout << __FUNCTION__ << endl;
	
	//パスの数×パスに含まれる最大粒子数　だけ繰り返す
	for(int indxX = 0; indxX < m_mk2DiPTHtoPRT.GetSizeX(); indxX++)
	{
		Vec3 preVec = Vec3(0.0f, 0.0f, 0.0f);
		float massSum = 0.0f;

		for(int indxY = 0; indxY < m_mk2DiPTHtoPRT.GetSizeY(); indxY++)
		{
			//TODO:末尾は必ず-1であることを想定している
			//TODO:穴あきの場合もあるので，パスに所属する粒子数とかも保存したほうがいいかも
			int pIndx = m_mk2DiPTHtoPRT(indxX, indxY);
			float mass = 1.0f;		//とりあえず1.0fで固定

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

//変形行列（moment matrix）Apqのprefix sum
void Surf_SM::CalcPrefixSumApq(const float *pos)
{	cout << __FUNCTION__ << endl;

	//パスの数×パスに含まれる最大粒子数　だけ繰り返す
	for(int indxX = 0; indxX < m_mk2DiPTHtoPRT.GetSizeX(); indxX++)
	{
		rxMatrix3 preMat(0.0f);
		float massSum = 0.0f;

		for(int indxY = 0; indxY < m_mk2DiPTHtoPRT.GetSizeY(); indxY++)
		{
			//TODO:末尾は必ず-1であることを想定している
			//TODO:穴あきの場合もあるので，パスに所属する粒子数とかも保存したほうがいいかも
			int pIndx = m_mk2DiPTHtoPRT(indxX, indxY);
			float mass = 1.0f;		//とりあえず1.0fで固定

			if(pIndx == -1)
			{
				break;
			}

			massSum += mass;

			Vec3 p = Vec3(pos[pIndx*4+0], pos[pIndx*4+1], pos[pIndx*4+2]);
			Vec3 q = m_vvec3OrgPos[pIndx];

			//現在のAij
			m_mk2Dmat3_PrfxApq(indxX, indxY)(0,0) = mass * p[0] * q[0];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(0,1) = mass * p[0] * q[1];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(0,2) = mass * p[0] * q[2];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(1,0) = mass * p[1] * q[0];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(1,1) = mass * p[1] * q[1];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(1,2) = mass * p[1] * q[2];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(2,0) = mass * p[2] * q[0];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(2,1) = mass * p[2] * q[1];
			m_mk2Dmat3_PrfxApq(indxX, indxY)(2,2) = mass * p[2] * q[2];

			m_mk2Dmat3_PrfxApq(indxX, indxY) += preMat;		//これまでのAijを加算
			preMat = m_mk2Dmat3_PrfxApq(indxX, indxY);
		}
	}
}


//--------------------------------------------------------デバッグ-------------------------------------------
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