/*! 
 @file rx_upsampling.cpp

 @brief パーティクルのアップサンプリング
		Solenthaler, B.; Zhang, Y. & Pajarola, R. 
		Efficient Refinement of Dynamic Point Data 
		Proceedings Eurographics/IEEE VGTC Symposium on Point-Based Graphics, 2007
 
 @author Go Mimura, Makoto Fujisawa
 @date 2010
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "rx_upsampling.h"


//-----------------------------------------------------------------------------
// rxUpsamplingPointクラスの実装
//-----------------------------------------------------------------------------

void rxUpsamplingPoint::updateNeighbor(const double r)
{
	const double rr = r*r;
	POINTVECTOR::iterator pend = PointList.end();
	for(POINTVECTOR::iterator pit = PointList.begin(); pit != pend; ++pit){
		rxPoint* particle = *pit;
		for(POINTVECTOR::iterator npit = particle->Neighbor.begin(); npit != particle->Neighbor.end();){
		rxPoint* nparticle = *npit;
			double dis2 = dis2Point(particle,nparticle);;
			if(dis2 >= rr) npit = particle->Neighbor.erase(npit);
			else ++npit;
		}
	}
}


void rxUpsamplingPoint::clearAllChildren()
{
	POINTVECTOR::iterator pend = PointList.end();
	for(POINTVECTOR::iterator pit = PointList.begin(); pit != pend; ++pit){
		(*pit)->Children.clear();
	}
}

void rxUpsamplingPoint::clearAllNeighbor()
{
	POINTVECTOR::iterator pend = PointList.end();
	for(POINTVECTOR::iterator pit = PointList.begin(); pit != pend; ++pit){
		(*pit)->Neighbor.clear();
	}
}

void rxUpsamplingPoint::clearAllChecked()
{
	POINTVECTOR::iterator pend = PointList.end();
	for(POINTVECTOR::iterator pit = PointList.begin(); pit != pend; ++pit){
		(*pit)->Checked.clear();
	}
}

void rxUpsamplingPoint::clearAllChildrenAndChecked()
{
	POINTVECTOR::iterator pend = PointList.end();
	for(POINTVECTOR::iterator pit = PointList.begin(); pit != pend; ++pit){
		(*pit)->Children.clear();
		(*pit)->Checked.clear();
	}
}

void rxUpsamplingPoint::Init()
{
	Step = 0;
	Num_p0 = PointList.size();
}

void rxUpsamplingPoint::AddPoint(const Vec3 &pos, const Vec3 &normal, const int step)
{
	rxPoint* point = new rxPoint;
	point->Pos = pos;
	point->Normal = normal;
	point->UpsamplingStep = Step;
	PointList.push_back(point);
}

void rxUpsamplingPoint::AddInitialPoint(rxPoint* particle)
{
	if(Step == 0) PointList.push_back(particle);
}



void rxUpsamplingPoint::Draw(void)
{
	glPushMatrix();
	glPointSize(2.0);
	glBegin(GL_POINTS);
	
	POINTVECTOR::iterator pend = PointList.end();
	for(POINTVECTOR::iterator pit = PointList.begin(); pit != pend; ++pit){
		switch( (*pit)->UpsamplingStep ){
			case 0: glColor3d(1.0, 0.0, 0.0);
					break;
			case 1:	glColor3d(0.0, 1.0, 0.0);
					break;
			default:glColor3d(0.0, 0.0, 1.0);
					break;
		}
		glVertex3dv((*pit)->Pos.data);;
	}
	glEnd();
	glPopMatrix();

}

void rxUpsamplingPoint::Upsampling(void)
{
	POINTVECTOR::iterator pend = PointList.end();
	for(POINTVECTOR::iterator pit = PointList.begin(); pit != pend; ++pit){
		(*pit)->UpsamplingStep = 0;
	}
	
	while(Step < IterationCount){
		Step++;
		update();
	}
	
}

void rxUpsamplingPoint::update()
{
	const double wc = 0.1;
	const double rr = Ri*Ri;
	const double newri = Ri * 0.50 * sqrt(4.0 - 2.0 * sqrt(2.0) );
	const double newrr = newri * newri;

	const double beta = 0.20 * Ri;
	const double beta2 = beta * beta;
	const double hmax = Ri * 0.5 * ( sqrt(2.0) - 1.0);	

	for(int i=PointList.size()-1;i>=0;i--){
		rxPoint* particle = PointList[i];
	/*POINTLIST::iterator pend = PointList.end();
	for(POINTLIST::iterator pit = PointList.begin(); pit != PointList.end(); ++pit){
		rxPoint* particle = *pit;*/

		if(particle->UpsamplingStep != Step){

			for(int j=particle->Neighbor.size()-1;j>=0;j--){
				rxPoint* nparticle = particle->Neighbor[j];
				//for(POINTLIST::iterator npit = particle->Neighbor.begin(); npit != particle->Neighbor.end(); ++npit){
				//rxPoint* nparticle = *npit;

				if(particle != nparticle){
					bool NOTMAKEPOINT = false;

					nparticle->Checked.push_back(particle);//重複を避ける
					POINTLIST::iterator cend = particle->Checked.end();
					POINTLIST::iterator cit = find(particle->Checked.begin(), cend, nparticle);
					if(cit == cend ){

						const double disdis = dis2Point(particle,nparticle);

						//Collision Type 1
						if(disdis <= rr && disdis > newrr){
							POINTVECTOR iList,uList;
							rxPoint *newParticle = makeChild(particle,nparticle,Step,hmax,div_wp);
							calUI(particle->Neighbor, nparticle->Neighbor,uList,iList);

							POINTVECTOR::iterator uLbegin	= uList.begin();
							POINTVECTOR::iterator uLend	= uList.end();
							POINTVECTOR::iterator iLbegin	= iList.begin();
							POINTVECTOR::iterator iLend	= iList.end();
							
							//Collision Type 2
							for(POINTVECTOR::iterator iLit = iLbegin; iLit != iLend && NOTMAKEPOINT == false; ++iLit){
								const double dis2 = dis2Point(*iLit, newParticle);
								if(dis2 < beta2){
									(*iLit)->CollisionHandling(newParticle,wc);
									delete newParticle;
									NOTMAKEPOINT = true;
									break;
								}
							}

							//Collision Type 3
							for(POINTVECTOR::iterator uLit = uLbegin; uLit != uLend && NOTMAKEPOINT == false; ++uLit){
									POINTLIST::iterator chend = (*uLit)->Children.end();
									for(POINTLIST::iterator chit = (*uLit)->Children.begin(); chit != chend; ++chit){
									double dis2 = dis2Point(newParticle,*chit);
									if(dis2 < beta2){
										(*chit)->CollisionHandling(newParticle,wc );
										delete newParticle;
										NOTMAKEPOINT = true;
										break;
									}
								}
							}

							//近傍粒子追加
							if(NOTMAKEPOINT == false){
								PointList.push_back(newParticle);
								particle->Children.push_back(newParticle);
								nparticle->Children.push_back(newParticle);

								//Neighbor Type 2								
								for(POINTVECTOR::iterator uLit = uLbegin; uLit != uLend; ++uLit){
									double dis2 = dis2Point(*uLit,newParticle);
									if(dis2 <= newrr){
										addNeighbor(*uLit,newParticle);
									}
								}

								//Neighbor Type 3
								for(POINTVECTOR::iterator uLit = uLbegin; uLit != uLend; ++uLit){
									POINTLIST::iterator chend = (*uLit)->Children.end();
									for(POINTLIST::iterator chit = (*uLit)->Children.begin(); chit != chend; ++chit){
										double dis2 = dis2Point(*chit,newParticle);
										if(dis2 <= newrr){
											addNeighbor(*chit,newParticle);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
	}

	//Neighbor Type 1
	updateNeighbor(newri);
	Ri = newri;
	clearAllChildrenAndChecked();
}

rxPoint* rxUpsamplingPoint::makeChild(const rxPoint *parent1,const rxPoint *parent2,const int step,const double &hmax, const double &div_wp)
{
	rxPoint* child;
	Vec3 pos;

	double h;

	const Vec3 d = parent1->Pos - parent2->Pos;
	const double norm_d = norm(d);
	const Vec3 a = d / norm_d;
	Vec3 tmpVec = parent1->Normal + parent2->Normal;
	tmpVec = tmpVec - dot(tmpVec,a) * a;
	const Vec3 nc = normalize(tmpVec);
	
	const double s1 = abs(dot(parent1->Normal,a));
	const double s2 = abs(dot(parent2->Normal,a));

	const double t1 = dot(parent1->Normal,nc) + 1.0;//norm(nc);
	const double t2 = dot(parent2->Normal,nc) + 1.0;//norm(nc);

	const double h1 = s1 / t1 * norm_d * 0.5;
	const double h2 = s2 / t2 * norm_d * 0.5;

	if(abs(h1) > abs(h2)) h=h2;
	else h=h1;

	if(abs(h)>hmax){
		h = hmax * h / h;
	}

	pos = 0.5 * (parent1->Pos + parent2->Pos) + h * nc;

	child = new rxPoint();
	if(div_wp != 0.0) child->Wp = disPoint(child,parent1) / div_wp;
	child->Pos = pos;
	child->Normal = nc;
	child->UpsamplingStep = step;

	return child;
}

void rxUpsamplingPoint::calUI(POINTVECTOR &list1,POINTVECTOR &list2, POINTVECTOR &uList, POINTVECTOR &iList)
{
	uList.clear();
	iList.clear();
	uList = list1;

	for(POINTVECTOR::iterator l2it = list2.begin(); l2it != list2.end(); ++l2it){
		POINTVECTOR::iterator l1it = find( list1.begin(), list1.end(), *l2it );
		if( l1it == list1.end() ) uList.push_back(*l2it);
		else iList.push_back(*l2it);
	}


}

double rxUpsamplingPoint::dis2Point(const rxPoint* a, const rxPoint* b)
{
	return norm2(a->Pos - b->Pos);
}

double rxUpsamplingPoint::disPoint(const rxPoint* a, const rxPoint* b)
{
	return norm(a->Pos - b->Pos);
}

void rxUpsamplingPoint::addNeighbor(rxPoint* a,rxPoint* b)
{
	a->Neighbor.push_back(b);
	b->Neighbor.push_back(a);
}
