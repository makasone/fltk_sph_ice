/*! 
  @file rx_ps.cpp

  @brief �p�[�e�B�N���������V�~�����[�V�����̊��N���X�̎���
 
  @author Makoto Fujisawa
  @date 2011-06
*/
// FILE --rx_ps.cpp--


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "rx_ps.h"

//-----------------------------------------------------------------------------
// �O���[�o���ϐ�
//-----------------------------------------------------------------------------
double g_fSurfThr[2] = {0.25, 0.35};

//-----------------------------------------------------------------------------
// rxParticleSystemBase�N���X�̎���
//-----------------------------------------------------------------------------
/*!
 * �p�[�e�B�N���f�[�^�̃Z�b�g
 * @param[in] ppos �p�[�e�B�N�����W
 * @param[in] pvel �p�[�e�B�N�����x
 */
bool rxParticleSystemBase::Set(const vector<Vec3> &ppos, const vector<Vec3> &pvel)
{
	// MARK:Set
	if(ppos.empty() || (int)ppos.size() != m_uNumParticles){
		return false;
	}

	int p = 0;
	for(uint i = 0; i < m_uNumParticles; ++i){
		Vec3 p0 = ppos[i];
		Vec3 v0 = pvel[i];
		for(uint j = 0; j < 3; ++j){
			m_hPos[DIM*i+j] = p0[j];
			m_hVel[DIM*i+j] = v0[j];
		}
		if(DIM == 4){
			m_hPos[DIM*i+3] = 0.0f;
			m_hVel[DIM*i+3] = 0.0f;
		}
	}

	SetArrayVBO(RX_POSITION, m_hPos, 0, m_uNumParticles);
	SetArrayVBO(RX_VELOCITY, m_hVel, 0, m_uNumParticles);

	SetParticlesToCell();

	SetColorVBO(m_iColorType);

	return true;
}


/*!
 * �V�[���̃��Z�b�g
 * @param[in] config �V�[���̎��
 */
void rxParticleSystemBase::Reset(rxParticleConfig config)
{
	// MRK:Reset
	switch(config){
	default:
	case RX_CONFIG_RANDOM:
		{
			int p = 0;
			for(uint i = 0; i < m_uNumParticles; ++i){
				Vec3 p0 = 2*Vec3(RX_FRAND(), RX_FRAND(), RX_FRAND())-0.5;
				for(uint j = 0; j < 3; ++j){
					m_hPos[DIM*i+j] = p0[j];
					m_hVel[DIM*i+j] = 0.0f;
				}
				if(DIM == 4){
					m_hPos[DIM*i+3] = 0.0f;
					m_hVel[DIM*i+3] = 0.0f;
				}
				m_hAttr[i] = 0;
			}
		}
		break;

	case RX_CONFIG_GRID:
		{
			uint s = (int) ceilf(powf((RXREAL) m_uNumParticles, 1.0f / 3.0f));
			AddBox(0, GetMin(), Vec3(s*m_fParticleRadius), Vec3(0.0), m_fParticleRadius*2);
		}
		break;

	case RX_CONFIG_NONE:
		for(uint i = 0; i < m_uMaxParticles; ++i){
			for(int j = 0; j < DIM; ++j){
				m_hVel[DIM*i+j] = 0.0f;
			}
		}
		SetArrayVBO(RX_VELOCITY, m_hVel, 0, m_uMaxParticles);
		m_uNumParticles = 0;
		break;
	}

	if(m_uNumParticles){
		SetArrayVBO(RX_POSITION, m_hPos, 0, m_uNumParticles);
		SetArrayVBO(RX_VELOCITY, m_hVel, 0, m_uNumParticles);

		SetParticlesToCell();

		SetColorVBO(m_iColorType);
	}
}


/*!
 * ���`����ɕ��ׂ�ꂽ�p�[�e�B�N���̒ǉ�
 * @param[in] start �ǉ��J�n�C���f�b�N�X
 * @param[in] pos[DIM] �ǉ��ʒu
 * @param[in] vel[DIM] �������x
 * @param[in] r ���a
 * @param[in] spacing �p�[�e�B�N���Ԋu
 * @return 
 */
void rxParticleSystemBase::AddSphere(int start, RXREAL *pos, RXREAL *vel, int r, RXREAL spacing, int attr)
{
	uint index;
	if(start < 0){
		index = m_uNumParticles;
		start = m_uNumParticles;
	}
	else{
		index = start;
	}

	int count = 0;
	bool over = false;
	for(int z = -r; z <= r; ++z){
		for(int y = -r; y <= r; ++y){
			for(int x = -r; x <= r; ++x){
				if(index >= m_uMaxParticles){
					index = 0;
					over = true;
				}

				RXREAL dx[3];
				dx[0] = x*spacing;
				dx[1] = y*spacing;
				dx[2] = z*spacing;
				RXREAL l = sqrtf(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
				RXREAL jitter = spacing*0.01f;
				if((l <= spacing*r) && (index < m_uNumParticles)) {
					for(uint j = 0; j < 3; ++j){
						m_hPos[DIM*index+j] = pos[j]+dx[j]+(RX_FRAND()*2.0f-1.0f)*jitter;
						m_hVel[DIM*index+j] = vel[j];
					}
					if(DIM == 4){
						m_hPos[DIM*index+3] = 0.0f;
						m_hVel[DIM*index+3] = 0.0f;
					}
					m_hAttr[index] = attr;

					index++;
					count++;
				}
			}
		}
	}
	RXCOUT << "num : " << count << endl;

	if(over){
		m_uNumParticles = m_uMaxParticles;
	}
	else{
		m_uNumParticles += count;
	}

	SetArrayVBO(RX_POSITION, m_hPos, start, index);
	SetArrayVBO(RX_VELOCITY, m_hVel, start, index);

	SetParticlesToCell();

	SetColorVBO(m_iColorType);
}

/*!
 * �{�b�N�X�`����ɕ��ׂ�ꂽ�p�[�e�B�N���̒ǉ�
 * @param[in] 
 * @param[in] 
 * @param[in] 
 * @param[in] 
 */
void rxParticleSystemBase::AddBox(int start, Vec3 cen, Vec3 ext, Vec3 vel, RXREAL spacing, int attr)
{
	uint index;
	if(start < 0){
		index = m_uNumParticles;
		start = m_uNumParticles;
	}
	else{
		index = start;
	}

	int count = 0;
	bool over = false;

	srand((unsigned)time(NULL));

	int sx = (int)(ext[0]/spacing)-1;
	int sy = (int)(ext[1]/spacing)-1;
	int sz = (int)(ext[2]/spacing)-1;
	RXREAL jitter = spacing*0.001f;

	for(int z = -sz; z <= sz; ++z){
		for(int y = -sy; y <= sy; ++y){
			for(int x = -sx; x <= sx; ++x){
				if(index >= m_uMaxParticles){
					index = 0;
					over = true;
				}

				RXREAL dx[3];
				dx[0] = x*spacing;
				dx[1] = y*spacing;
				dx[2] = z*spacing;

				for(uint j = 0; j < 3; ++j){
					m_hPos[DIM*index+j] = cen[j]+dx[j]+(RX_FRAND()*2.0f-1.0f)*jitter;
					m_hVel[DIM*index+j] = vel[j];
				}
				if(DIM == 4){
					m_hPos[DIM*index+3] = 0.0f;
					m_hVel[DIM*index+3] = 0.0f;
				}
				m_hAttr[index] = attr;

				index++;
				count++;
			}
		}
	}
	RXCOUT << "num : " << count << endl;

	if(over){
		m_uNumParticles = m_uMaxParticles;
	}
	else{
		m_uNumParticles += count;
	}

	//SetArrayVBO(RX_POSITION, m_hPos, start, count);
	//SetArrayVBO(RX_VELOCITY, m_hVel, start, count);
	SetArrayVBO(RX_POSITION, m_hPos, 0, m_uNumParticles);
	SetArrayVBO(RX_VELOCITY, m_hVel, 0, m_uNumParticles);

	SetParticlesToCell();

	SetColorVBO(m_iColorType);
}
/*!
 * �t�̃p�[�e�B�N���̗������C�����V�[���ɒǉ�
 * @param[in] line  �������C��
 * @return �ݒ肳�ꂽ�������C����
 */
int rxParticleSystemBase::AddLine(rxInletLine line)
{
	m_vInletLines.push_back(line);
	return (int)m_vInletLines.size();
}

/*!
 * �t�̃p�[�e�B�N���̗������C�����V�[���ɒǉ�
 * @param[in] start �ǉ��J�n�C���f�b�N�X
 * @param[in] line  �������C��
 * @return �ǉ����ꂽ�p�[�e�B�N����
 */
int rxParticleSystemBase::addParticles(int &start, rxInletLine line, int attr)
{
	// Vec3 pos1, Vec3 pos2, Vec3 vel, Vec3 up, int accum, int n, double s
	RXREAL rel[DIM];	// �[�_�Ԃ̑��Έʒu�x�N�g��
	RXREAL l = 0.0;		// �[�_�Ԃ̋���
	for(int i = 0; i < 3; ++i){
		rel[i] = line.pos2[i]-line.pos1[i];
		l += rel[i]*rel[i];
	}
	l = sqrt(l);

	RXREAL pr = GetParticleRadius();	// �p�[�e�B�N�����a
	int n = l/(pr);	// ���ׂ�p�[�e�B�N���̐�
	if(!n) return 0;

	RXREAL spacing = l/n;
	RXREAL jitter = spacing*0.01f;

	uint index;
	if(start < 0){
		index = m_uNumParticles;
		start = m_uNumParticles;
	}
	else{
		index = start;
	}

	int count = 0;
	bool over = false;
	for(int j = 0; j < line.accum; ++j){
		for(int i = 0; i < n; ++i){
			if(index >= m_uMaxParticles){
				index = 0;
				over = true;
			}

			for(int k = 0; k < 3; ++k){
				m_hPos[DIM*index+k] = line.pos1[k]+rel[k]/l*(i+0.5)*spacing+(RX_FRAND()*2.0f-1.0f)*jitter;
				m_hVel[DIM*index+k] = line.vel[k];
			}
			if(DIM == 4){
				m_hPos[DIM*index+3] = 0.0f;
				m_hVel[DIM*index+3] = 0.0f;
			}
			m_hAttr[index] = attr;
			index++;
			count++;
		}

		line.pos1 += line.up*line.spacing*spacing;
	}

	if(over){
		m_uNumParticles = m_uMaxParticles;

		count = m_uMaxParticles-start;
		//SetArrayVBO(RX_POSITION, m_hPos, start, count);
		//SetArrayVBO(RX_VELOCITY, m_hVel, start, count);

		start = 0;
	}
	else{
		if(m_uNumParticles < m_uMaxParticles){
			m_uNumParticles += count;
		}

		//SetArrayVBO(RX_POSITION, m_hPos, start, count);
		//SetArrayVBO(RX_VELOCITY, m_hVel, start, count);

		start += count;
	}

	//SetParticlesToCell();

	SetColorVBO(m_iColorType);

	return count;
}




/*!
 * �p�[�e�B�N���̃J���[�o�b�t�@�l���v�Z
 * @param[in] hVal �z�X�g��������̔z��
 * @param[in] d �z��̃X�e�b�v
 * @param[in] n �v�f��(�z��̃T�C�Y��n*d)
 * @param[in] use_max �ő�l�Œl�𐳋K��
 * @param[in] vmax �蓮�ݒ�ő�l
 * @return �l�̍ő�l
 */
RXREAL rxParticleSystemBase::SetColorVBOFromArray(RXREAL *hVal, int d, bool use_max, RXREAL vmax)
{
	if(m_bUseOpenGL){
		RXREAL l = 1.0;

		RXREAL max_val = 0.0;
		for(uint i = 0; i < m_uNumParticles; ++i){
			if(hVal[d*i] > max_val) max_val = hVal[d*i];
		}

		l = max_val;
		if(!use_max) l = vmax;
		//RXCOUT << "max val : " << max_val << endl;

		// �p�[�e�B�N���J���[�o�b�t�@
		glBindBufferARB(GL_ARRAY_BUFFER, m_colorVBO);
		RXREAL *data = (RXREAL*)glMapBufferARB(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
		//RXREAL *ptr = data;
		for(uint i = 0; i < m_uNumParticles; ++i){
			RXREAL value = hVal[d*i]/l;

			double col[3];
			Gradation(col, value, 0.0, 1.0);
			data[4*i+0] = col[0];
			data[4*i+1] = col[1];
			data[4*i+2] = col[2];

			//value = RX_CLAMP(value, (RXREAL)0.0, (RXREAL)1.0);
			//data[4*i+0] = value;
			//data[4*i+1] = 0.0;
			//data[4*i+2] = 1.0-value;
			data[4*i+3] = 1.0;
		}
		glUnmapBufferARB(GL_ARRAY_BUFFER);

		return max_val;
	}

	return (RXREAL)0.0;
}


/*!
 * �p�[�e�B�N�����̏o��
 * @param[in] fn �o�̓t�@�C����
 */
int rxParticleSystemBase::OutputParticles(string fn)
{
	ofstream fout;
	fout.open(fn.c_str(), ios::out|ios::binary);
	if(!fout){
		RXCOUT << fn << " couldn't open." << endl;
		return 0;
	}

	fout.write((char*)&m_uNumParticles, sizeof(uint));
	//fout << m_uNumParticles << endl;
	for(uint i = 0; i < m_uNumParticles; ++i){
		for(int j = 0; j < 3; ++j){
			fout.write((char*)&m_hPos[DIM*i+j], sizeof(RXREAL));
		}
	}
	for(uint i = 0; i < m_uNumParticles; ++i){
		for(int j = 0; j < 3; ++j){
			fout.write((char*)&m_hVel[DIM*i+j], sizeof(RXREAL));
		}
	}

	fout.close();

	return 1;
}


/*!
 * �t�@�C������p�[�e�B�N������ǂݍ���
 * @param[in] stp �X�e�b�v��(sph_�X�e�b�v��.dat�̃t�@�C��������ǂݍ���)
 * @param[out] ppos �p�[�e�B�N�����W
 * @param[out] pvel �p�[�e�B�N�����x
 */
int rxParticleSystemBase::InputParticles(string fn)
{
	ifstream fin;
	fin.open(fn.c_str(), ios::in|ios::binary);
	if(!fin){
		RXCOUT << fn << " couldn't find." << endl;
		return 0;
	}

	uint n;
	fin.read((char*)&n, sizeof(uint));

	for(uint i = 0; i < n; ++i){
		for(int j = 0; j < 3; ++j){
			fin.read((char*)&m_hPos[DIM*i+j], sizeof(RXREAL));
		}
	}

	if(!fin.eof()){
		for(uint i = 0; i < n; ++i){
			for(int j = 0; j < 3; ++j){
				fin.read((char*)&m_hVel[DIM*i+j], sizeof(RXREAL));
			}
		}
	}

	fin.close();
	
	return 1;
}
