#pragma once
#include "Points.h"
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cmath>

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
    double C_Dot;
    Vec3 deltaP2;
    // Vec3 deltaP;
    Vec3 dc_dp;
    Mat3 df_dx,df_dv;
    Mat3 J;
    Mat3 E;
    float inv_len;
    Mat3 I;
    Mat3 d2C_dp1p1;
    Mat3 d2C_dp1p2;
    Mat3 d2C_dp2p1;
    Mat3 d2C_dp2p2;
    const double DEFAULT_DAMPING = 5.0;


    
	Spring(Node *n1, Node *n2, double k)
	{
        node1 = n1;
        node2 = n2;
        Vec3 currSp = node2->position - node1->position;
        restLength = currSp.length();   //rest length
        Ks = k;
        Kd = 70.0;
        C = 0.0;
        
        
        // deltaP.setZeroVec();
        deltaP2.setZeroVec();
        dc_dp.setZeroVec();
        df_dx.setZeromat3();
        df_dv.setZeromat3();
        
        E.setIdentity();
       

	}

	// void applyInternalForce(double timeStep) // Compute spring internal force
	// {
    //     deltaP = node1->position - node2->position; 
    //     double currLen = Vec3::dist(node1->position, node2->position);
    //     Vec3 dc_dp = (deltaP)/currLen;   // dc_dp
    //     // Vec3 dc_dp = (node2->position - node1->position)/currLen;   // dc_dp
    //     Vec3 deltaV = node2->velocity - node1->velocity;
    //     Vec3 f1 = dc_dp * ((currLen-restLength)*Ks + Vec3::dot(deltaV, dc_dp)*Kd);
    //     node1->addForce(f1);
    //     node2->addForce(f1.minus());

    
    //     double dist = deltaP.length();

    //     double deltaP22 = deltaP.length()*deltaP.length();
    //     double lo_l = restLength/dist;

    //     J = (((E-(Mat3::outerProduct(deltaP,deltaP)*(1/deltaP22)))*lo_l)-E)*Ks;
    //     // node1->addForceDerivative(J);

    //     C = deltaP.length()-restLength; //dist2
    //     if (isinf(C)){
    //         printf("Printing node1 and node 2 pos:");
    //         node1->printPosition();
    //         node2->printPosition();
    //         printf("node1->position =  %f , restLength =  %f \n", node1->position.x,restLength); 
    //     }
    //     printf("C = %f \n", C); 
       

    //     // dc_dp = deltaP/(deltaP.length());
    //     // C_Dot[i] = glm::dot(v1, -dc_dp[i]) + glm::dot(v2, dc_dp[i]);
    //     c_dot = Vec3::dot(node1->velocity,(dc_dp*-1.0))+Vec3::dot(node2->velocity,(dc_dp)); //check negative signs
        
        
    // }

    void applyInternalForce(double timeStep) // Compute spring internal force
	{
        double currLen = Vec3::dist(node1->position, node2->position);
        // Vec3 dc_dp = (deltaP)/currLen;   // dc_dp
        Vec3 dc_dp = (node2->position - node1->position)/currLen;   // orig: p1-p2
        C = currLen-restLength; //see if restlength matches
        Vec3 deltaV = node2->velocity - node1->velocity;
        Vec3 deltaP = node2->position - node1->position; //originally p1-p1
        Vec3 v1 = node1->velocity;
        Vec3 v2 = node2->velocity;
        Vec3 f = dc_dp * ((C)*Ks + Vec3::dot(deltaV, dc_dp)*Kd);
        C_Dot = Vec3::dot(v1,dc_dp) + Vec3::dot(v2,dc_dp*(-1.0));
        deltaP2 = Vec3(deltaP.x * deltaP.x, deltaP.y * deltaP.y, deltaP.z * deltaP.z);
        node1->addForce(f);
        node2->addForce(f.minus());
        
        // todo: look at gravity and damping velocity if statement

        // Vec3 deltaP = node1->position - node2->position; 
        // Vec3 deltaV = node1->velocity - node2->velocity;
        // Vec3 v1 = node1->velocity;
        // Vec3 v2 = node2->velocity;
        // float dist = Vec3::dist(node1->position, node2->position);
        // double currLen = Vec3::dist(node1->position, node2->position);
        // inv_len = 1.0f/dist;
        // C = dist-restLength; //see if restlength matches
        // dc_dp = (deltaP)/currLen;   // dc_dp
        // C_Dot =  Vec3::dot(v1, dc_dp*(-1.0)) +  Vec3::dot(v2, dc_dp*(-1.0));
        // deltaP2 = Vec3(deltaP.x * deltaP.x, deltaP.y * deltaP.y, deltaP.z * deltaP.z);
        // float leftTerm = -Ks * (dist-restLength);
		// float rightTerm = Kd * (Vec3::dot(deltaV, deltaP)/dist);
        // deltaP.normalize();
		// Vec3 springForce = deltaP*(leftTerm + rightTerm);
        // // node1->addForce(springForce);
        // // node2->addForce(springForce.minus());
        // deltaP = node1->position - node2->position;

        

        // need 
       //C_Dot =  Vec3::dot(v1, dc_dp*(-1.0)) +  Vec3::dot(v2, dc_dp*(-1.0));
        //deltaP2 = Vec3(deltaP.x * deltaP.x, deltaP.y * deltaP.y, deltaP.z * deltaP.z);

        // Vec3 deltaP = node1->position - node2->position; 
        // Vec3 deltaV = node1->velocity - node2->velocity;
        // Vec3 v1 = node1->velocity;
        // Vec3 v2 = node2->velocity;
        // float dist = Vec3::dist(node1->position, node2->position);
        // double currLen = Vec3::dist(node1->position, node2->position);
        // inv_len = 1.0f/dist;
        // C = dist-restLength; //see if restlength matches
        // dc_dp = (deltaP)/currLen;   // dc_dp
        // C_Dot =  Vec3::dot(v1, dc_dp*(-1.0)) +  Vec3::dot(v2, dc_dp*(-1.0));
        // deltaP2 = Vec3(deltaP.x * deltaP.x, deltaP.y * deltaP.y, deltaP.z * deltaP.z);
        // float leftTerm = -Ks * (dist-restLength);
		// float rightTerm = Kd * (Vec3::dot(deltaV, deltaP)/dist);
        // //deltaP.normalize();
		// Vec3 f = deltaP*(leftTerm + rightTerm);
        // node1->addForce(springForce);
        // node2->addForce(springForce.minus());
        

    
        

        
        // double lo_l = restLength/dist;

      
        // // node1->addForceDerivative(J);

       
        // if (isinf(C)){
        //     printf("Printing node1 and node 2 pos:");
        //     node1->printPosition();
        //     node2->printPosition();
        //     printf("node1->position =  %f , restLength =  %f \n", node1->position.x,restLength); 
        // }

       

        
        
    }

    void springForceDerivative(){
        I.setIdentity();
        float c1 = C;
        d2C_dp1p1.identityMult ( (deltaP2 * (c1* (-1.0) ) ) + c1);
        d2C_dp1p2.identityMult((deltaP2 * c1) - c1);
        d2C_dp2p1 = d2C_dp1p2;
        d2C_dp2p2 = d2C_dp1p1;
        df_dx = Mat3::outerProduct(dc_dp, dc_dp) * -Ks+ Mat3(1e-8);
        df_dv = Mat3::outerProduct(dc_dp, dc_dp) * Kd+ Mat3(1e-8);


        //Mat3 df_dx = Mat3::outerProduct(dc_dp, dc_dp) * -Ks + Mat3(1e-8); // Adding small regularization
        Mat3 dp1 = Mat3::outerProduct(dc_dp, dc_dp);
        Mat3 dp2 = Mat3::outerProduct(dc_dp, dc_dp*(-1.0));
        Mat3 dp3 = Mat3::outerProduct(dc_dp*(-1.0),dc_dp*(-1.0));

        Mat3 term1 =  ( (dp1+(d2C_dp1p1 * c1)) * (-Ks)) - ((d2C_dp1p1 * C_Dot)*(Kd)); //replaced kd with -0.125
        Mat3 term2 =  ( (dp2+(d2C_dp1p2 * c1)) * (-Ks)) - ((d2C_dp1p2 * C_Dot)*(Kd));

        //df_dx = term1 + ( term2*2); // term2 is used twice
        //df_dv = (dp1 + dp2 + dp3)*(Kd);
        node1->addForceDerivative(df_dx,df_dv);
    }
};
