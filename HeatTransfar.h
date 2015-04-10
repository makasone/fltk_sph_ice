//---------------------------------------------------------------------------
//���̃t�@�C���̃N���X�E�֐������ȏ�R���p�C�������̂�h�����߂̏����i�C���N���[�h�K�[�h�j
#ifndef HeatTransfarH
#define HeatTransfarH
//---------------------------------------------------------------------------

#include <cuda_runtime.h>

#include <vector>
#include "Math2d\math2d.h"

//GPU����
extern void LaunchHeatTransferGPU
(
	float* heats,
	float* temps,
	float* dtemps,
	const int* surfParticles,
	float airTemp,
	float cffCnHt,
	int prtNum
	//const vector<vector<rxNeigh>>& neights,
	//const vector<int>& objNeight,
	//float floor,
	//float effRadius,
	//const float* pos,
	//const float* dens
);

class HeatTransfar
{
public:
	HeatTransfar(int num);
	~HeatTransfar(void);

	void InitGPU();

	void MeltParticle(int pIndx);
	void WarmParticle(int pIndx, float temp, float heat);

//----------------------------------------GPU----------------------------------------------
	float* getHostHeats(){	return sd_Heats;	}
	float* getHostTemps(){	return sd_Temps;	}
	float* getHostDTemps(){	return sd_DTemps;	}

	int* getHostPhase(){	return sd_Phase;	}
	int* getHostPhaseChangeFlag(){	return sd_PhaseChangeFlag;	}
//----------------------------------------GPU----------------------------------------------

	float* getTemps(){	return mTemps;	}								//�e���q�̉��x�z����擾
	float* getHeats(){	return mHeats;	}								//�e���q�̔M�ʔz����擾
	int* getSurfaceParticleNums(){	return mSurfaceParticleNums;	}	//�e�\�ʗ��q�̔ԍ����擾

	int getAirTemp() {	return mAirTemp;	}							//��C���x���擾

	float getTempMax(){	return mTempMax;	}							//��E�i�ő�j���x���擾
	float getTempMin(){ return mTempMin;	}							//���E�i�ŏ��j���x���擾

	float getCffCntHt(){ return mHT;	}
	float getCffCntTd(){ return mTD;	}

	float getLatentHeat(){ return mLatentHeat; }						//�Z����M���擾

	int getPhase(int i){		return mPhase[i];	}					//���݂̏�Ԃ�Ԃ��i���E�X�j
	int getPhaseChange(int i){	return mPhaseChange[i]; }

	void setNumVertices(int num){	mNumVertices = num;	};				//���q�̌���ݒ�
	
	void setTimeStep( float time ){ timeStep = time; };
	
	void setTempMax(float max){	mTempMax = max; }
	void setTempMin(float min){ mTempMin = min; }
	
	void setLatentHeat(float heat){	mLatentHeat = heat; }

	void setTemps(int nr, double tmp){	mTemps[nr] = tmp;	};			//����Y�����̉��x��ݒ�
	void setHeats(int nr, double heat){	mHeats[nr] = heat;	};									//����Y�����̔M�ʂ�ݒ�
	void setAirTemp(float temp){	mAirTemp = temp;	};				//��C���x��ݒ�

	void setCffCntHt(float h){ mHT = h; };
	void setCffCntTd(float h){ mTD = h;	};

	void setPhase(int i, int phase){ mPhase[i] = phase; }
	void setPhaseChange(int i, int phase){ mPhaseChange[i] = phase; }

	void setSurfaceParticleNums(int nr, int num){	mSurfaceParticleNums[nr] = num;	}	//����\�ʗ��q�̓Y�������i�[

	void resetNeighborhoodsId();										//�e�ߖT���q�̓Y������������
	void resetNeighborhoodsDist();										//�e�ߖT���q�̋�����������

	void AddNeighborhoodsId(std::vector<int> ids);						//�e�ߖT���q�̓Y������ݒ�
	void AddNeighborhoodsDist(std::vector<float> dists);				//�e�ߖT���q�̋�����ݒ�

	void AddParticle(int nowVerticesNum);								//sph�@�ŗ��q���ǉ����ꂽ�ۂ̏���

	void heatAirAndParticle();											//��C�ƃp�[�e�B�N���̔M����
	void heatObjAndParticle(const std::vector<int>& neight);					//�I�u�W�F�N�g�ƃp�[�e�B�N���̔M����
	void heatParticleAndParticle(const float* d, double h);				//�p�[�e�B�N���Ԃ̔M����

	void calcTempAndHeat();												//�M�ʂ��牷�x�����߂�
	void calcTempAndHeat(int pIndx);									//1�̗��q�����M�ʂ��牷�x�����߂�

	// spline�J�[�l���@�J�[�l���֐�
	inline double KernelSpline(const double &r, const double &h);
//	inline Vec2   KernelSplineG(const double &r, const Vec2 &rij, const double &h);
	inline double KernelSplineL(const double &r, const double &h);

	void setCarnelConstant(float radius);

private:
	void initState();						//������

//----------------------------------------GPU----------------------------------------------
	static float* sd_Heats;
	static float* sd_Temps;
	static float* sd_DTemps;

	static int* sd_Phase;
	static int* sd_PhaseChangeFlag;
//----------------------------------------GPU----------------------------------------------

	float timeStep;							//�^�C���X�e�b�v

	float mAirHeat;							//��C�̔M��	���g�p
	float mAirTemp;							//��C�̉��x

	float mTempMax;							//���x�̏���l
	float mTempMin;							//���x�̉����l

	float mLatentHeat;						//�Z����M	readOnly�ɂ���H

	double mHT;								//�M�`�B�W��	coefficient of heat transfer
	double mTD;								//�M�g�U�W��	coefficient of thermal diffusivity


	int mNumVertices;						//���_��

	float *mTemps;							//���x
	float *mTempsDelta;						//�ω����x
	float *mHeats;							//�M��

	int *mPhase;							//��ԁ@-2���X�@-1���X���ԁ@1�������ԁ@2����
	int *mPhaseChange;						//�O�t���[���̏�ԁ@�O�����ω��Ȃ��@�P�����ω�����

	int *mSurfaceParticleNums;				//�\�ʗ��q�̋ߖT���q��

	std::vector < std::vector<int> >	mNeighborhoodsId;				//�e���q�̋ߖT�ƂȂ闱�q��ID
	std::vector < std::vector<float> >	mNeighborhoodsDis;				//�e���q�̋ߖT�ƂȂ闱�q�̋���

	// �J�[�l���֐��̌v�Z�̍ۂɗp������萔�W��
	double m_fWpoly6;				//!< Pory6�J�[�l���̒萔�W��
	double m_fGWpoly6;				//!< Pory6�J�[�l���̌��z�̒萔�W��
	double m_fLWpoly6;				//!< Pory6�J�[�l���̃��v���V�A���̒萔�W��
	double m_fWspiky;				//!< Spiky�J�[�l���̒萔�W��
	double m_fGWspiky;				//!< Spiky�J�[�l���̌��z�̒萔�W��
	double m_fLWspiky;				//!< Spiky�J�[�l���̃��v���V�A���̒萔�W��
	double m_fWvisc;				//!< Viscosity�J�[�l���̒萔�W��
	double m_fGWvisc;				//!< Viscosity�J�[�l���̌��z�̒萔�W��
	double m_fLWvisc;				//!< Viscosity�J�[�l���̃��v���V�A���̒萔�W��

	double m_fWspline;				//!< Spline�J�[�l���̒萔�W��
	double m_fGWspline;				//!< Spline�J�[�l���̌��z�̒萔�W��
	double m_fLWspline;				//!< Spline�J�[�l���̃��v���V�A���̒萔�W��

};
#endif