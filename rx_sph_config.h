/*!
  @file rx_sph_config.h
	
  @brief SPH�̃V�[���ݒ���t�@�C������ǂݍ���
 
  @author Makoto Fujisawa
  @date 2012-08
*/
// FILE --rx_sph_config.h--

#ifndef _RX_SPH_CONFIG_H_
#define _RX_SPH_CONFIG_H_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
// STL
#include <vector>
#include <string>

// boost
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include <boost/ref.hpp>

// ���[�e�B���e�B
#include "rx_utility.h"

// �V�~�����[�V����
#include "rx_sph_commons.h"
#include "rx_sph.h"

// 3D���f��
#include "rx_model.h"
#include "rx_fltk_widgets.h"

#define RX_USE_BOOST

// �ݒ�t�@�C��
#include "rx_atom_ini.h"

using namespace std;


//-----------------------------------------------------------------------------
// �����񏈗��֐�
//-----------------------------------------------------------------------------
/*!
 * �����񂩂�l(Vec3)���擾
 * @param[out] val �l
 * @param[in] str ������
 * @param[in] rel true�ŃV�~�����[�V������Ԃ̑傫���ɑ΂���W���Ƃ��Čv�Z
 * @param[in] cen �V�~�����[�V������Ԃ̒��S���W
 * @param[in] ext �V�~�����[�V������Ԃ̑傫��(�e�ӂ̒�����1/2)
 */
inline void GetValueFromString(Vec3 &val, const string &str, bool rel = false, Vec3 cen = Vec3(0.0), Vec3 ext = Vec3(1.0))
{
	int n = StringToVec3(str, val);
	if(rel){
		val = cen+(n == 1 ? RXFunc::Min3(ext) : ext)*val;
	}
}

/*!
 * �����񂩂�l(double)���擾
 * @param[out] val �l
 * @param[in] str ������
 * @param[in] rel true�ŃV�~�����[�V������Ԃ̑傫���ɑ΂���W���Ƃ��Čv�Z
 * @param[in] cen �V�~�����[�V������Ԃ̒��S���W
 * @param[in] ext �V�~�����[�V������Ԃ̑傫��(�e�ӂ̒�����1/2)
 */
inline void GetValueFromString(double &val, const string &str, bool rel = false, Vec3 cen = Vec3(0.0), Vec3 ext = Vec3(1.0))
{
	val = atof(str.c_str());
	if(rel){
		val = RXFunc::Min3(ext)*val;
	}
}



/*!
 * OBJ�t�@�C���ǂݍ���
 * @param[in] filename wrl�t�@�C���̃p�X
 */
inline static void ReadOBJ(const string filename, rxPolygons &polys, Vec3 cen, Vec3 ext, Vec3 ang)
{
	if(!polys.vertices.empty()){
		polys.vertices.clear();
		polys.normals.clear();
		polys.faces.clear();
		polys.materials.clear();
	}
	rxOBJ obj;
	if(obj.Read(filename, polys.vertices, polys.normals, polys.faces, polys.materials, true)){
		RXCOUT << filename << " have been read." << endl;

		if(polys.normals.empty()){
			CalVertexNormals(polys);
		}

		RXCOUT << " the number of vertex   : " << polys.vertices.size() << endl;
		RXCOUT << " the number of normal   : " << polys.normals.size() << endl;
		RXCOUT << " the number of polygon  : " << polys.faces.size() << endl;
		RXCOUT << " the number of material : " << polys.materials.size() << endl;

		//FitVertices(Vec3(0.0), Vec3(1.0), polys.vertices);
		AffineVertices(polys, cen, ext, ang);

		// �e�N�X�`���ǂݍ���
		if(!polys.materials.empty()){
			rxMTL::iterator iter = polys.materials.begin();
			for(; iter != polys.materials.end(); ++iter){
				if(iter->second.tex_file.empty()) continue;

				RXCOUT << iter->first << " : " << iter->second.tex_file;
				LoadGLTexture(iter->second.tex_file, iter->second.tex_name, true, false);

				RXCOUT << " : " << iter->second.tex_name << endl;
			}
		}

		polys.open = 1;
	}
}

//-----------------------------------------------------------------------------
//! rxSPHConfig�N���X - SPH�̃V�[���ݒ���t�@�C������ǂݍ���
//-----------------------------------------------------------------------------
class rxSPHConfig
{
protected:
	rxParticleSystemBase *m_pPS;	//!< SPH
	rxSPHEnviroment m_SphEnv;		//!< SPH���ݒ�
	
	string m_strCurrentScene;			//!< ���݂̃V�[���̖��O
	vector<string> m_vSceneFiles;		//!< �V�[���t�@�C�����X�g
	int m_iSceneFileNum;				//!< �V�[���t�@�C���̐�

	vector<string> m_vSceneTitles;		//!< �V�[���t�@�C�����X�g
	int m_iCurrentSceneIdx;				//!< ���݂̃V�[���t�@�C��

	vector<rxPolygons> m_vSolidPoly;	//!< �ő̃��b�V��

public:
	//! �f�t�H���g�R���X�g���N�^
	rxSPHConfig() : m_pPS(0)
	{
		m_vSceneFiles.resize(12, "");	// �V�[���t�@�C�����X�g
		m_vSceneTitles.resize(12, "");	// �V�[���^�C�g�����X�g
		m_iCurrentSceneIdx = -1;		// ���݂̃V�[���t�@�C��
		m_strCurrentScene = "null";		// ���݂̃V�[���̖��O
		m_iSceneFileNum = 0;

		Clear();
	}

	//! �f�X�g���N�^
	~rxSPHConfig(){}
	
	//! �ݒ菉����
	void Clear(void)
	{
		m_pPS = 0;
		m_vSolidPoly.clear();
	}


	//! �V�[���^�C�g�����X�g
	vector<string> GetSceneTitles(void) const { return m_vSceneTitles; }

	//! ���݂̃V�[��
	int GetCurrentSceneIdx(void) const { return m_iCurrentSceneIdx; }

	//! SPH��
	rxSPHEnviroment GetSphEnv(void) const { return m_SphEnv; }

	//! SPH�N���X
	void SetPS(rxParticleSystemBase *ps){ m_pPS = ps; }

	//! �ő̃|���S��
	int GetSolidPolyNum(void) const { return (int)m_vSolidPoly.size(); }
	vector<rxPolygons>& GetSolidPolys(void){ return m_vSolidPoly; }

public:
	/*!
	 * �ݒ肩��p�����[�^�̓ǂݍ���
	 * @param[in] names ���ږ����X�g
	 * @param[in] values �l���X�g
	 * @param[in] n ���X�g�̃T�C�Y
	 * @param[in] header �w�b�_��
	 */
	void SetSphSpace(string *names, string *values, int n, string header)
	{
		rxSPHEnviroment sph_env;
		sph_env.use_inlet = 0;
		sph_env.et_cri = 1.0;
		int idr = 0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")							GetValueFromString(sph_env.boundary_cen, values[i], false);
			else if(names[i] == "ext")						GetValueFromString(sph_env.boundary_ext, values[i], false);
			else if(names[i] == "max_particle_num")			sph_env.max_particles = atoi(values[i].c_str());
			else if(names[i] == "density")					sph_env.dens = atof(values[i].c_str());
			else if(names[i] == "mass")						sph_env.mass = atof(values[i].c_str());
			else if(names[i] == "kernel_particles")			sph_env.kernel_particles = atof(values[i].c_str());
			else if(names[i] == "mesh_res_max")				sph_env.mesh_max_n = atoi(values[i].c_str());
			else if(names[i] == "inlet_boundary")			sph_env.use_inlet = atoi(values[i].c_str());
			else if(names[i] == "dt")						sph_env.dt = atof(values[i].c_str());
			else if(names[i] == "viscosity")				sph_env.viscosity = atof(values[i].c_str());
			else if(names[i] == "gas_stiffness")			sph_env.gas_k = atof(values[i].c_str());
			else if(names[i] == "init_vertex_store")		sph_env.mesh_vertex_store = atoi(values[i].c_str());
			else if(names[i] == "et_cri")					sph_env.et_cri = atof(values[i].c_str());
			else if(names[i] == "use_delete_region")		sph_env.use_delete_region = atoi(values[i].c_str());
			else if(names[i] == "delete_region_min"){		GetValueFromString(sph_env.delete_region[idr/2][0], values[i], false); idr++; }
			else if(names[i] == "delete_region_max"){		GetValueFromString(sph_env.delete_region[idr/2][1], values[i], false); idr++; }
			else if(names[i] == "epsilon")					sph_env.epsilon = atof(values[i].c_str());
			else if(names[i] == "dens_fluctuation")			sph_env.eta = atof(values[i].c_str());
			else if(names[i] == "min_iterations")			sph_env.min_iter = atoi(values[i].c_str());
			else if(names[i] == "max_iterations")			sph_env.max_iter = atoi(values[i].c_str());
			else if(names[i] == "use_artificial_pressure")	sph_env.use_ap = atoi(values[i].c_str());
			else if(names[i] == "ap_k")						sph_env.ap_k = atof(values[i].c_str());
			else if(names[i] == "ap_n")						sph_env.ap_n = atof(values[i].c_str());
			else if(names[i] == "ap_q")						sph_env.ap_q = atof(values[i].c_str());
			/*!
			 * �ǉ��F�F�M�����Ɋւ���p�����[�^�����ݒ���t�@�C������ǂݍ���
			 */
			else if(names[i] == "htTimeStep")				sph_env.htTimeStep = atof(values[i].c_str());
			else if(names[i] == "tempMax")					sph_env.tempMax = atof(values[i].c_str());
			else if(names[i] == "tempMin")					sph_env.tempMin = atof(values[i].c_str());
			else if(names[i] == "latentHeat")				sph_env.latentHeat = atof(values[i].c_str());
			else if(names[i] == "cffcntHt")					sph_env.cffcntHt = atof(values[i].c_str());
			else if(names[i] == "cffcntTd")					sph_env.cffcntTd = atof(values[i].c_str());
			/*!
			 * �ǉ��F�F�X�\���Ɋւ���p�����[�^�����ݒ���t�@�C������ǂݍ���
			 */
			else if(names[i] == "smTimeStep")				sph_env.smTimeStep = atof(values[i].c_str());
			else if(names[i] == "smItr")					sph_env.smItr = atof(values[i].c_str());
			
			/*!
			 * �ǉ��F�F�X�\���Ɋւ���p�����[�^�����ݒ���t�@�C������ǂݍ���
			 */
			else if(names[i] == "layer")					sph_env.layer = atof(values[i].c_str());

		}
		if(sph_env.mesh_vertex_store < 1) sph_env.mesh_vertex_store = 1;

		sph_env.mesh_boundary_ext = sph_env.boundary_ext;
		sph_env.mesh_boundary_cen = sph_env.boundary_cen;

		m_SphEnv = sph_env;

	}

	//! �t�� : ���`
	void SetLiquidBox(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		Vec3 cen(0.0), ext(0.0), vel(0.0);
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")      GetValueFromString(cen, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "ext") GetValueFromString(ext, values[i], rel, Vec3(0.0), 0.5*m_pPS->GetDim());
			else if(names[i] == "vel") GetValueFromString(vel, values[i], false);
		}
		//m_pPS->AddBox(-1, cen, ext, vel, -1);
		m_pPS->AddBox(-1, cen, ext, vel, 2*m_pPS->GetParticleRadius());
		RXCOUT << "set liquid box : " << cen << ", " << ext << ", " << vel << endl;
	}

	//! �ǉ��H�@������
	//! �t�� : ��
	void SetLiquidSphere(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		Vec3 cen(0.0), vel(0.0);
		double rad = 0.0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")      GetValueFromString(cen, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "rad") GetValueFromString(rad, values[i], rel, Vec3(0.0), 0.5*m_pPS->GetDim());
			else if(names[i] == "vel") GetValueFromString(vel, values[i], false);
		}
		//m_pPS->AddSphere(-1, cen, rad, vel, -1);
		RXCOUT << "set liquid sphere : " << cen << ", " << rad << endl;
	}

	//! �ǉ��F�t�́F���`�F�\�ʂ̂�
	void SetLiquidBoxSurface(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		Vec3 cen(0.0), ext(0.0), vel(0.0);
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")      GetValueFromString(cen, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "ext") GetValueFromString(ext, values[i], rel, Vec3(0.0), 0.5*m_pPS->GetDim());
			else if(names[i] == "vel") GetValueFromString(vel, values[i], false);
		}
		
		m_pPS->AddBoxSurface(-1, cen, ext, vel, 2*m_pPS->GetParticleRadius());
		RXCOUT << "set liquid box surface : " << cen << ", " << ext << ", " << vel << endl;
	}

	//! �t�̗��� : ����
	void SetInletLine(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		Vec3 pos1(0.0), pos2(0.0), vel(0.0), up(0.0, 1.0, 0.0);
		int  span = -1, accum = 1;
		double spacing = 1.0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "pos1")      GetValueFromString(pos1, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "pos2") GetValueFromString(pos2, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "vel")  GetValueFromString(vel,  values[i], false);
			else if(names[i] == "up")   GetValueFromString(up,   values[i], false);
			else if(names[i] == "span") span = atoi(values[i].c_str());
			else if(names[i] == "accum") accum = atoi(values[i].c_str());
			else if(names[i] == "spacing") spacing = atof(values[i].c_str());
		}

		rxInletLine inlet;
		inlet.pos1 = pos1;
		inlet.pos2 = pos2;
		inlet.vel = vel;
		inlet.span = span;
		inlet.up = up;
		inlet.accum = accum;
		inlet.spacing = spacing;

		int num_of_inlets = m_pPS->AddLine(inlet);
		//((RXSPH*)m_pPS)->AddSubParticles(g_iInletStart, count);
		RXCOUT << "set inlet boundary : " << pos1 << "-" << pos2 << ", " << vel << endl;
		RXCOUT << "                     span=" << span << ", up=" << up << ", accum=" << accum << ", spacing=" << spacing << endl;
	}

	//! �ő� : ���`
	void SetSolidBox(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		Vec3 cen(0.0), ext(0.0), ang(0.0), vel(0.0);
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")      GetValueFromString(cen, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "ext") GetValueFromString(ext, values[i], rel, Vec3(0.0), 0.5*m_pPS->GetDim());
			else if(names[i] == "ang") GetValueFromString(ang, values[i], false);
			else if(names[i] == "vel") GetValueFromString(vel, values[i], false);
		}
		m_pPS->SetBoxObstacle(cen, ext, ang, vel, 1);
		RXCOUT << "set solid box : " << cen << ", " << ext << ", " << ang << endl;
	}

	//! �ő� : ��
	void SetSolidSphere(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		Vec3 cen(0.0), move_pos1(0.0), move_pos2(0.0), vel(0.0);
		int  move = 0, move_start = -1;
		double rad = 0.0, move_max_vel = 0.0, lap = 1.0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")      GetValueFromString(cen, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "rad") GetValueFromString(rad, values[i], rel, Vec3(0.0), 0.5*m_pPS->GetDim());
			else if(names[i] == "vel")  GetValueFromString(vel, values[i], false);
			else if(names[i] == "move_pos1") GetValueFromString(move_pos1, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "move_pos2") GetValueFromString(move_pos2, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "move") move = atoi(values[i].c_str());
			else if(names[i] == "move_start") move_start = atoi(values[i].c_str());
			else if(names[i] == "move_max_vel") move_max_vel = atof(values[i].c_str());
			else if(names[i] == "lap") lap = atof(values[i].c_str());
		}
		m_pPS->SetSphereObstacle(cen, rad, vel, 1);
		//if(move){
		//	m_vMovePos[0] = move_pos1;
		//	m_vMovePos[1] = move_pos2;
		//	m_fMoveMaxVel = move_max_vel;

		//	//m_bMoveSolid = false;
		//	if(move_start >= 0){
		//		m_iMoveStart = move_start;
		//	}
		//}
		RXCOUT << "set solid sphere : " << cen << ", " << rad << endl;
	}

	//! �ő� : �|���S��
	void SetSolidPolygon(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		string fn_obj;
		Vec3 cen(0.0), ext(0.0), ang(0.0), vel(0.0);
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")       GetValueFromString(cen, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "ext")  GetValueFromString(ext, values[i], rel, Vec3(0.0), 0.5*m_pPS->GetDim());
			else if(names[i] == "ang")  GetValueFromString(ang, values[i], false);
			else if(names[i] == "vel")  GetValueFromString(vel, values[i], false);
			else if(names[i] == "file") fn_obj = values[i];
		}
		if(!fn_obj.empty()){
			rxPolygons poly;
			ReadOBJ(fn_obj, poly, cen, ext, ang);
			m_vSolidPoly.push_back(poly);
			vector< vector<int> > tris;
			int pn = poly.faces.size();
			tris.resize(pn);
			for(int i = 0; i < pn; ++i){
				tris[i].resize(3);
				for(int j = 0; j < 3; ++j){
					tris[i][j] = poly.faces[i][j];
				}
			}
			m_pPS->SetPolygonObstacle(poly.vertices, poly.normals, tris, vel);
			RXCOUT << "set solid polygon : " << fn_obj << endl;
			RXCOUT << "                  : " << cen << ", " << ext << ", " << ang << endl;
		}
	}

	//! ���b�V���������E�͈�
	void SetMeshBoundary(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		Vec3 cen(0.0), ext(0.0);
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")       GetValueFromString(cen, values[i], rel, m_pPS->GetCen(), 0.5*m_pPS->GetDim());
			else if(names[i] == "ext")  GetValueFromString(ext, values[i], rel, Vec3(0.0), 0.5*m_pPS->GetDim());
		}

		m_SphEnv.mesh_boundary_cen = cen;
		m_SphEnv.mesh_boundary_ext = ext;
		RXCOUT << "boundary for mash : " << cen << ", " << ext << endl;
	}

	/*!
	 * �V�~�����[�V������Ԃ̐ݒ�ǂݍ���
	 */
	bool LoadSpaceFromFile(void)
	{
		bool ok = true;
		rxINI *cfg = new rxINI();
		cfg->SetHeaderFunc("space", boost::bind(&rxSPHConfig::SetSphSpace, this, _1, _2, _3, _4));
		if(!(cfg->Load(m_strCurrentScene))){
			RXCOUT << "Failed to open the " << m_strCurrentScene << " file!" << endl;
			ok = false;
		}
		delete cfg;
		return ok;
	}
	
	/*!
	 * �p�[�e�B�N����ő̃I�u�W�F�N�g�̐ݒ�ǂݍ���
	 */
	bool LoadSceneFromFile(void)
	{
		bool ok = true;
		rxINI *cfg = new rxINI();
		cfg->SetHeaderFunc("liquid box",		boost::bind(&rxSPHConfig::SetLiquidBox, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("liquid box (r)",    boost::bind(&rxSPHConfig::SetLiquidBox, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("liquid sphere",		boost::bind(&rxSPHConfig::SetLiquidSphere, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("liquid sphere (r)", boost::bind(&rxSPHConfig::SetLiquidSphere, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("solid box",			boost::bind(&rxSPHConfig::SetSolidBox, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("solid box (r)",		boost::bind(&rxSPHConfig::SetSolidBox, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("solid sphere",		boost::bind(&rxSPHConfig::SetSolidSphere, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("solid sphere (r)",	boost::bind(&rxSPHConfig::SetSolidSphere, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("solid polygon",		boost::bind(&rxSPHConfig::SetSolidPolygon, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("solid polygon (r)", boost::bind(&rxSPHConfig::SetSolidPolygon, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("inlet line",		boost::bind(&rxSPHConfig::SetInletLine, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("inlet line (r)",	boost::bind(&rxSPHConfig::SetInletLine, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("mesh grid",			boost::bind(&rxSPHConfig::SetMeshBoundary, this, _1, _2, _3, _4));
		cfg->SetHeaderFunc("mesh grid (r)",		boost::bind(&rxSPHConfig::SetMeshBoundary, this, _1, _2, _3, _4));
		//�ǉ��F���̕\�ʂ݂̂ɗ��q��z�u
		cfg->SetHeaderFunc("liquid box surface(r)",		boost::bind(&rxSPHConfig::SetLiquidBoxSurface, this, _1, _2, _3, _4));
		if(!(cfg->Load(m_strCurrentScene))){
			RXCOUT << "Failed to open the " << m_strCurrentScene << " file!" << endl;
			ok = false;
		}
		delete cfg;
		return ok;
	}

	void LoadSpaceFromFile(const string input)
	{
		// SPH�ݒ���t�@�C������ǂݍ���
		ifstream fsin;
		fsin.open(input.c_str());

		Vec3 bmin, bmax;
		fsin >> m_SphEnv.max_particles;
		fsin >> bmin[0] >> bmin[1] >> bmin[2];
		fsin >> bmax[0] >> bmax[1] >> bmax[2];
		fsin >> m_SphEnv.dens;
		fsin >> m_SphEnv.mass;
		fsin >> m_SphEnv.kernel_particles;

		m_SphEnv.boundary_ext = 0.5*(bmax-bmin);
		m_SphEnv.boundary_cen = 0.5*(bmax+bmin);

		fsin.close();

		RXCOUT << "[SPH - " << input << "]" << endl;
		RXCOUT << " num. of particles : " << m_SphEnv.max_particles << endl;
		RXCOUT << " boundary min      : " << bmin << endl;
		RXCOUT << " boundary max      : " << bmax << endl;
		RXCOUT << " boundary cen      : " << m_SphEnv.boundary_cen << endl;
		RXCOUT << " boundary ext      : " << m_SphEnv.boundary_ext << endl;
		RXCOUT << " density           : " << m_SphEnv.dens << endl;
		RXCOUT << " mass              : " << m_SphEnv.mass << endl;
		RXCOUT << " kernel particles  : " << m_SphEnv.kernel_particles << endl;

		m_SphEnv.mesh_boundary_cen = m_SphEnv.boundary_cen;
		m_SphEnv.mesh_boundary_ext = m_SphEnv.boundary_ext;

	}

	/*!
	 * �w�肵���t�H���_�ɂ���ݒ�t�@�C���̐��ƃV�[���^�C�g����ǂݎ��
	 * @param[in] dir �ݒ�t�@�C��������t�H���_(�����w�肵�Ȃ���Ύ��s�t�H���_)
	 */
	void ReadSceneFiles(string dir = "")
	{
		m_vSceneFiles.resize(12, "");	// �V�[���t�@�C�����X�g
		m_vSceneTitles.resize(12, "");	// �V�[���^�C�g�����X�g

		ifstream scene_ifs;
		string scene_fn = "null";
		for(int i = 1; i <= 12; ++i){
			if(ExistFile((scene_fn = CreateFileName(dir+"sph_scene_", ".cfg", i, 1)))){
				RXCOUT << "scene " << i << " : " << scene_fn << endl;
				m_vSceneFiles[i-1] = scene_fn;
				m_vSceneTitles[i-1] = scene_fn.substr(0, 11);

				// �V�[���^�C�g���̓ǂݎ��
				scene_ifs.open(scene_fn.c_str(), ios::in);
				string title_buf;
				getline(scene_ifs, title_buf);
				if(!title_buf.empty() && title_buf[0] == '#'){
					m_vSceneTitles[i-1] = title_buf.substr(2, title_buf.size()-2);
				}
				scene_ifs.close();

				m_iSceneFileNum++;
			}
		}

		if(m_iSceneFileNum){
			SetCurrentScene(0);
		}
	}

	/*!
	 * �J�����g�̃V�[���ݒ�
	 * @param[in] �V�[���C���f�b�N�X
	 */
	bool SetCurrentScene(int idx)
	{
		if(idx < 0 || idx >= m_iSceneFileNum || m_vSceneFiles[idx] == ""){
			cout << "There is no scene files!" << endl;
			return false;
		}
		m_iCurrentSceneIdx = idx;
		m_strCurrentScene = m_vSceneFiles[m_iCurrentSceneIdx];
		return true;
	}

	/*!
	 * �^�C�g������J�����g�̃V�[���ݒ�
	 * @param[in] label �V�[���^�C�g��
	 */
	bool SetCurrentSceneFromTitle(const string label)
	{
		int scene = 0;
		while(label.find(m_vSceneTitles[scene]) == string::npos) scene++;

		if(m_iCurrentSceneIdx != -1 && m_vSceneFiles[scene] != ""){
			RXCOUT << "scene " << scene+1 << " : " << m_vSceneFiles[scene] << endl;
			m_iCurrentSceneIdx = scene;
			m_strCurrentScene = m_vSceneFiles[scene];
			return true;
		}

		return false;
	}


	vector<string> GetSceneTitleList(void){ return m_vSceneTitles; }
};


#endif // #ifndef _RX_SPH_CONFIG_H_