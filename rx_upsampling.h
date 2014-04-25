/*! 
 @file rx_upsampling.h

 @brief パーティクルのアップサンプリング
		Solenthaler, B.; Zhang, Y. & Pajarola, R. 
		Efficient Refinement of Dynamic Point Data 
		Proceedings Eurographics/IEEE VGTC Symposium on Point-Based Graphics, 2007
 
 @author Go Mimura, Makoto Fujisawa
 @date 2010
*/


#ifndef _RX_UPSAMPLING_H_
#define _RX_UPSAMPLING_H_


#include "rx_sph_commons.h"

class rxPoint;

typedef list<rxPoint*> POINTLIST;
typedef vector<rxPoint*> POINTVECTOR;

class rxPoint
{
public:
	Vec3 Pos;
	Vec3 Normal;
	int UpsamplingStep;
	int Step;

	POINTVECTOR Neighbor;
	POINTLIST Children;
	POINTLIST Checked;

	void ClearChild(){
		Children.clear();
	}

	void ClearNeighbor(){
		Neighbor.clear();
	}

//private:
	double Wp;

public:
	rxPoint()
	{
		Wp = 1.0;
		//Pos = Vec3(0.0, 0.0, 0.0);
		//Normal = Vec3(0.0, 0.0, 1.0);
	}

	~rxPoint()
	{
		ClearChild();
		ClearNeighbor();
	}

	inline rxPoint(const Vec3 &x, const Vec3 &n, const int step)
	{
		Pos = x;
		Normal = n; 
		Step = step;
		Wp = 1.0;
	}

	void AddChild(rxPoint* point)
	{
		Children.push_back(point);
	}

	void AddNeighbor(rxPoint* point)
	{
		Neighbor.push_back(point);
	}

	void CollisionHandling(const rxPoint *child, const double &wc)
	{
		if(child != NULL){
			Pos		= ( Wp * Pos + wc * child->Pos) / (Wp+wc);
			Normal = ( Wp * Normal + wc * child->Normal);
			normalize(Normal);
			Wp += wc;
		}
	}
};



//-----------------------------------------------------------------------------
// アップサンプル点作成，格納
//-----------------------------------------------------------------------------
class rxUpsamplingPoint
{
public:
	int Step;
	double Ri;
	int IterationCount;
	int Num_p0;

//private:
	POINTVECTOR PointList;
private:
	double div_wp;

public:
	rxUpsamplingPoint(const double r, const int iteration)
	{
		Ri = r;
		div_wp = 5.0*r;
		Step = 0;
		Num_p0 = 0;
		IterationCount = iteration;
	}

	~rxUpsamplingPoint()
	{
		POINTVECTOR::iterator itend = PointList.end();
		for(POINTVECTOR::iterator it = PointList.begin(); it != PointList.end();){
			if( (*it)->UpsamplingStep != 0){
				delete *it;
				it = PointList.erase(it);
			}
			else ++it;
		}

		//clearAllChildrenAndChecked();
		clearAllChildren();
		clearAllNeighbor();

	}
	
	
public:
	void Init(void);
	void AddPoint(const Vec3 &pos, const Vec3 &normal, const int step);
	void AddInitialPoint(rxPoint* particle);
	void Draw(void);
	void Upsampling(void);

private:
	void updateNeighbor(const double r);
	void clearAllChildren();
	void clearAllNeighbor();
	void clearAllChecked();
	void clearAllChildrenAndChecked();
	void update();
	rxPoint* makeChild(const rxPoint *parent1,const rxPoint *parent2,const int step,const double &hmax, const double &div_wp);
	void calUI(POINTVECTOR &list1,POINTVECTOR &list2, POINTVECTOR &uList, POINTVECTOR &iList);
	double dis2Point(const rxPoint* a, const rxPoint* b);
	double disPoint(const rxPoint* a, const rxPoint* b);
	void addNeighbor(rxPoint* a, rxPoint* b);
};


#endif // #ifndef _RX_UPSAMPLING_H_