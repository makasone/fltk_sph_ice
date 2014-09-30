
#include <Surf_SM.h>

//------------------------------------------------初期化----------------------------------------------
/*!
 * prefixSumを計算するためのパス初期化
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

	//まずは，立方体の表面で試すために適当なパスを作ってみる．
	//なぜかだめとされている１本のパスでやってみる．
	//パス→粒子
	int pthNum = 1;								//とりあえず1本
	cout << "check1" << endl;
	m_mk2DiPTHtoPRT.SetSize(pthNum, pthSize);

	for(int i = 0; i < pthNum; i++)
	{
		for(int j = 0; j < pthSize; j++)
		{
			m_mk2DiPTHtoPRT(i, j) = -1;		//-1で初期化
		}
	}
	cout << "check2" << endl;
	for(int i = 0; i < prtNum; i++)
	{
		m_mk2DiPTHtoPRT(0, i) = i;			//実際のパスの初期化 1本用
	}
	cout << "check3" << endl;
	//粒子→パス
	m_mk2DiPRTtoPTH.SetSize(prtNum, 2);

	//１本用
	for(int i = 0; i < prtNum; i++)
	{
		m_mk2DiPRTtoPTH(i, 0) = 0;
		m_mk2DiPRTtoPTH(i, 1) = i;
	}

	InitOrgPos(prtNum);							//初期位置

	//prefixSumを初期化
	m_mk2Dvec3_PrfxPos.SetSize(pthNum, prtNum);	//位置
	UpdatePrefixSumPos();

	m_mk2Dmat3_PrfxApq.SetSize(pthNum, prtNum);	//変形行列
	UpdatePrefixSumApq();

	InitPathPrfxIndxSet(iceSM, strct);			//どのパスのどの部分が必要なのか，をクラスタごとに計算

	InitOrgCm();								//初期重心

//デバッグ
	//DebugPathDataPos();
	//DebugPathDataApq();
	//DebugPathPrfxIndxSet();
}

//GPU
void Surf_SM::InitPathGPU()
{	cout << __FUNCTION__ << endl;
	
	//メモリ確保
	int pthNum = 1;
	int pthSize = m_iPrtclNum;

	cudaMalloc((void**)&md_2DiPTHtoPRT,	sizeof(int) * pthNum * pthSize);		//パス→粒子
	cudaMalloc((void**)&md_2DiPRTtoPTH,	sizeof(int) * m_iPrtclNum * 2);			//粒子→パス

	md_f3OrgPos = Ice_SM::GetOrgPosPointer();									//クラスタ内の粒子の初期位置
	md_f3OrgCm	= Ice_SM::GetOrgCmPointer();									//クラスタの初期重心		

	cudaMalloc((void**)&md_2Df3PrfxPos,	sizeof(float) * pthNum * pthSize * 3);	//位置のprefixSum
	cudaMalloc((void**)&md_2Df9PrfxApq,	sizeof(float) * pthNum * pthSize * 9);	//変形のprefixSum

	cudaMalloc((void**)&md_3DiPTHandPrfxSet,	sizeof(int) * m_iPrtclNum * m_iPrtclNum * 2);	//変形のprefixSum


	//初期化
	//vectorを使っているので，第二引数がちょっとややこしく見えるがアドレスを渡しているだけ．
	cudaMemcpy(md_2DiPTHtoPRT, &m_mk2DiPTHtoPRT.Get()[0], sizeof(int) * pthNum * pthSize, cudaMemcpyHostToDevice);
	cudaMemcpy(md_2DiPRTtoPTH, &m_mk2DiPRTtoPTH.Get()[0], sizeof(int) * m_iPrtclNum * 2 , cudaMemcpyHostToDevice);

	cudaMemcpy(md_2Df3PrfxPos, &m_mk2Dvec3_PrfxPos.Get()[0], sizeof(float) * pthNum * pthSize * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(md_2Df9PrfxApq, &m_mk2Dmat3_PrfxApq.Get()[0], sizeof(float) * pthNum * pthSize * 9, cudaMemcpyHostToDevice);

	cudaMemcpy(md_3DiPTHandPrfxSet, &m_mk3DiPTHandPrfxSet.Get()[0], sizeof(int) * m_iPrtclNum * m_iPrtclNum * 2, cudaMemcpyHostToDevice);


//デバッグ
	////データを転送
	//int* PTHtoPRT = new int[pthNum * pthSize];
	//int* PRTtoPTH = new int[m_iPrtclNum * 2];

	//cudaMemcpy(PTHtoPRT, md_2DiPTHtoPRT, sizeof(int) * pthNum * pthSize, cudaMemcpyDeviceToHost);
	//cudaMemcpy(PRTtoPTH, md_2DiPRTtoPTH, sizeof(int) * m_iPrtclNum * 2 , cudaMemcpyDeviceToHost);

	////ホスト側のデータを転送した結果をダンプ
	//ofstream ofs( "Surf_SM.txt" );
	//ofs << __FUNCTION__ << endl;
	//
	////デバイス側のデータを転送
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


//初期位置の初期化
void Surf_SM::InitOrgPos(int prtNum)
{
	m_vvec3OrgPos.resize(prtNum);

	for(int pIndx = 0; pIndx < prtNum; pIndx++)
	{
		m_vvec3OrgPos[pIndx] = Vec3(m_fPos[pIndx*4+0], m_fPos[pIndx*4+1], m_fPos[pIndx*4+2]);
	}
}

//初期重心の初期化
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

//どのパスのどの部分が必要なのか，をクラスタごとに計算
void Surf_SM::InitPathPrfxIndxSet(const vector<Ice_SM*> iceSM, IceStructure* strct)
{
	//サイズは，クラスタ数×近傍粒子数
	//そんなにたくさんはいらないはずだが，凝固を考慮して倍ぐらいにしておく？
	//m_mk3DiPTHandPrfxSet.SetSize(m_iPrtclNum, m_iPrtclNum*2, 2);
	//高速化を試す場合は，最低限に
	m_mk3DiPTHandPrfxSet.SetSize(m_iPrtclNum, m_iPrtclNum, 2);

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
		for(int jprt = 0; jprt < iceSM[iclstr]->GetNumVertices(); jprt++)
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
			for(int kord = jordIndx+1; kord < m_mk2DiPTHtoPRT.GetSizeY(); kord++)
			{
				unsigned kprtIndx = m_mk2DiPTHtoPRT(jpthIndx, kord);	//次の粒子
				int ksearchIndx = iceSM[iclstr]->SearchIndx(kprtIndx);	//クラスタ内での順番を取得

				if(ksearchIndx != -1)
				{
					m_mk3DiPTHandPrfxSet(iclstr, isetIndx, 1) = kord;	//存在するなら終点を更新
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


void Surf_SM::UpdatePrefixSum()
{
	//UpdatePrefixSumPos();
	//UpdatePrefixSumApq();

	//特に効果なし
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

	//特に効果なし
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

//重心のprefix sum
void Surf_SM::UpdatePrefixSumPos()
{//	cout << __FUNCTION__ << endl;

	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	//パスの数×パスに含まれる最大粒子数　だけ繰り返す
	for(int indxX = 0; indxX < sizeX; ++indxX)
	{
		Vec3 preVec(0.0, 0.0, 0.0);
		int pIndx = 0;

		for(int indxY = 0; indxY < sizeY; ++indxY)
		{
			//TODO:末尾は必ず-1であることを想定している
			//TODO:穴あきの場合もあるので，パスに所属する粒子数とかも保存したほうがいいかも
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*4;
			if(pIndx < 0){	break;	}

			double mass = 1.0;		//とりあえず1.0で固定

			m_mk2Dvec3_PrfxPos(indxX, indxY) = Vec3(m_fPos[pIndx+0], m_fPos[pIndx+1], m_fPos[pIndx+2]) * mass
											+ Vec3(m_fVel[pIndx+0], m_fVel[pIndx+1], m_fVel[pIndx+2]) * 0.02
											+ preVec;

			preVec = m_mk2Dvec3_PrfxPos(indxX, indxY);
		}
	}
}

//重心のprefix sum　反復処理
void Surf_SM::UpdatePrefixSumPosItr()
{//	cout << __FUNCTION__ << endl;

	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	float* pos = Ice_SM::GetSldPosPointer();
	float* vel = Ice_SM::GetSldVelPointer();

	//パスの数×パスに含まれる最大粒子数　だけ繰り返す
	for(int indxX = 0; indxX < sizeX; ++indxX)
	{
		Vec3 preVec(0.0, 0.0, 0.0);
		int pIndx = 0;

		for(int indxY = 0; indxY < sizeY; ++indxY)
		{
			//TODO:末尾は必ず-1であることを想定している
			//TODO:穴あきの場合もあるので，パスに所属する粒子数とかも保存したほうがいいかも
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*3;
			if(pIndx < 0){	break;	}

			double mass = 1.0;		//とりあえず1.0で固定

			m_mk2Dvec3_PrfxPos(indxX, indxY) = Vec3(pos[pIndx+0], pos[pIndx+1], pos[pIndx+2]) * mass + preVec;

			preVec = m_mk2Dvec3_PrfxPos(indxX, indxY);
		}
	}
}

//変形行列（moment matrix）Apqのprefix sum
void Surf_SM::UpdatePrefixSumApq()
{//	cout << __FUNCTION__ << endl;
	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	//パスの数×パスに含まれる最大粒子数　だけ繰り返す
	for(int indxX = 0; indxX < sizeX; ++indxX)
	{
		rxMatrix3 preMat(0.0);
		int pIndx = 0;
		Vec3 p;
		Vec3 q;

		for(int indxY = 0; indxY < sizeY; indxY++)
		{
			//TODO:末尾は必ず-1であることを想定している
			//TODO:穴あきの場合もあるので，パスに所属する粒子数とかも保存したほうがいいかも
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*4;
			if(pIndx < 0){	break;	}
			double mass = 1.0;		//とりあえず1.0で固定

			p = Vec3(m_fPos[pIndx+0], m_fPos[pIndx+1], m_fPos[pIndx+2])
				+ Vec3(m_fVel[pIndx+0], m_fVel[pIndx+1], m_fVel[pIndx+2]) * 0.02;

			q = m_vvec3OrgPos[pIndx/4];

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

//変形行列（moment matrix）Apqのprefix sum 反復処理
void Surf_SM::UpdatePrefixSumApqItr()
{//	cout << __FUNCTION__ << endl;
	int sizeX = m_mk2DiPTHtoPRT.GetSizeX();
	int sizeY = m_mk2DiPTHtoPRT.GetSizeY();

	float* pos = Ice_SM::GetSldPosPointer();
	float* vel = Ice_SM::GetSldVelPointer();

	//パスの数×パスに含まれる最大粒子数　だけ繰り返す
	for(int indxX = 0; indxX < sizeX; ++indxX)
	{
		rxMatrix3 preMat(0.0);
		int pIndx = 0;
		Vec3 p;
		Vec3 q;

		for(int indxY = 0; indxY < sizeY; indxY++)
		{
			//TODO:末尾は必ず-1であることを想定している
			//TODO:穴あきの場合もあるので，パスに所属する粒子数とかも保存したほうがいいかも
			pIndx = m_mk2DiPTHtoPRT(indxX, indxY)*3;
			if(pIndx < 0){	break;	}
			double mass = 1.0;		//とりあえず1.0で固定

			p = Vec3(pos[pIndx+0], pos[pIndx+1], pos[pIndx+2]);

			q = m_vvec3OrgPos[pIndx/3];

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

//重心の総和を返す
Vec3 Surf_SM::CalcCmSum(const int& cIndx)
{	//cout << __FUNCTION__ << " start" << endl;
	Vec3 cmSum(0.0);

	int start = -1;
	int end = -1;

	int prtIndx = -1;
	int pthIndx = -1;

	//クラスタに用意したデータセットを使って重心を求める
	for(int iprt = 0, ctopNum = m_strct->GetCtoPNum(cIndx); iprt < ctopNum; iprt++)
	{
		start	= m_mk3DiPTHandPrfxSet(cIndx, iprt, 0);
		end		= m_mk3DiPTHandPrfxSet(cIndx, iprt, 1);

		if(start == -1 || end == -1){	break;	}
		
		prtIndx = m_strct->GetCtoP(cIndx, iprt, 0);
		pthIndx = m_mk2DiPRTtoPTH(prtIndx, 0);

		//関数を使うバージョン
		//cmSum += CalcCmFromPrfxSm(pthIndx, start, end);

		//関数を使わないバージョン
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

//prfixSumから値を返す
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
	t = CalcCmSum(cIndx);			//TODO: これを計算せずに読み込めばもっと速くなる

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
		//関数を使わないバージョン
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
			if( m_mk3DiPTHandPrfxSet(iX, iY, 0) == -1 || m_mk3DiPTHandPrfxSet(iX, iY, 1) == -1 ){	cout << iY-1 << " ";	break;	}
			
			//cout << "m_mk3DiPTHandPrfxSet(" << iX << ", " << iY << ", 0) = " << m_mk3DiPTHandPrfxSet(iX, iY, 0)
			//	 << ", (" << iX << ", " << iY << ", 1) = " << m_mk3DiPTHandPrfxSet(iX, iY, 1) << endl;
		}
	}
}