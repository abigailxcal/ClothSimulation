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
    double C;
    double c_dot;
    Vec3 deltaP2;
    Vec3 deltaP;
    Vec3 dc_dp;
    Mat3 df_dx,df_dv;

    
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
        Vec3 fDir1 = (node2->position - node1->position)/currLen;   // dc_dp
        Vec3 diffV1 = node2->velocity - node1->velocity;
        Vec3 f1 = fDir1 * ((currLen-restLength)*Ks + Vec3::dot(diffV1, fDir1)*Kd);
        node1->addForce(f1);
        node2->addForce(f1.minus());

        deltaP = node1->position - node2->position; 
        C = deltaP.length()-restLength;
        dc_dp = deltaP/(deltaP.length());
        // C_Dot[i] = glm::dot(v1, -dc_dp[i]) + glm::dot(v2, dc_dp[i]);
        c_dot = Vec3::dot(node1->velocity,(dc_dp*-1.0))+Vec3::dot(node2->velocity,(dc_dp)); //check negative signs
        deltaP2 = deltaP*deltaP;
        
    }

    void springForceDerivative(){

        float c1 = C;
        Vec3 deltaP2_c1 = deltaP2*c1;
        Mat3 d2C_dp2_00(0.0);
        Mat3 d2C_dp2_01(0.0);

        //clear derivatives in Points.h?

        Vec3 val = (deltaP2_c1*(-1.0)) + c1;
        d2C_dp2_00.x1 = val.x;
        d2C_dp2_00.y2 = val.y;
        d2C_dp2_00.z3 = val.z;
        d2C_dp2_01.x1 = -val.x;
        d2C_dp2_01.y2 = -val.y;
        d2C_dp2_01.z3 = -val.z;
        
        Mat3 dp1 = Mat3::outerProduct(dc_dp, dc_dp);
        Mat3 dp2 = Mat3::outerProduct(dc_dp, (dc_dp*-1.0));
        Mat3 dp3 = Mat3::outerProduct((dc_dp*-1.0),(dc_dp*-1.0));

        Mat3 term1 =  ( ((d2C_dp2_00 * C)+dp1) * (-Ks)) -(d2C_dp2_00 * c_dot)*Kd;
        Mat3 term2 =  ( ((d2C_dp2_01 * C)+dp2) * (-Ks)) -(d2C_dp2_01 * c_dot)*Kd;

        // are df_dx,df_dv really dim(3x3)?
        df_dx += term1 + ( term2*2); // term2 is used twice
        df_dv += (dp1 + dp2 + dp3)*(-Kd);

        node1->addForceDerivative(df_dx,df_dv);
        
        
        
    
        
    
    }
};
