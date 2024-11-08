#pragma once

#include "Points.h"
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace std;

class Spring
{
public:
    Node *node1;
    Node *node2;
	double restLength;
    double Ks;
    double Kd;
    
	Spring(Node *n1, Node *n2, double k)
	{
        node1 = n1;
        node2 = n2;
		
        Vec3 currSp = node2->position - node1->position;
        restLength = currSp.length();   //rest length
        Ks = k;
        Kd = 5.0;
       

	}

	void applyInternalForce(double timeStep) // Compute spring internal force
	{
        double currLen = Vec3::dist(node1->position, node2->position);
        Vec3 fDir1 = (node2->position - node1->position)/currLen;
        Vec3 diffV1 = node2->velocity - node1->velocity;
        Vec3 f1 = fDir1 * ((currLen-restLength)*Ks + Vec3::dot(diffV1, fDir1)*Kd);
        node1->addForce(f1);
        node2->addForce(f1.minus());
	}

    void springForceDerivative(){
        Vec3 deltaP = node2->position - node1->position;
        double dist = deltaP.length();
        Vec3 dir = deltaP/dist;

        Mat3 dp1p1 =(Mat3(1.0) - Mat3::outerProduct(dir, dir))* Ks ;
        Mat3 dp1p2 = dp1p1*(-1.0);
        Mat3 dp2p1 = dp1p2;
        Mat3 dp2p2 = dp1p1;

        Mat3 p1_df_dx =dp1p1;        
        Mat3 p1_df_dv = dp1p1*Kd;
        
        Mat3 p2_df_dx = dp2p2;
        Mat3 p2_df_dv = dp2p2*Kd;
        node1->addForceDerivative(p1_df_dx,p1_df_dv);
        node2->addForceDerivative(p2_df_dx,p2_df_dv);
    }
};
