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
    Vec3 deltaV;
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
    


    
	Spring(Node *n1, Node *n2, double k)
	{
        node1 = n1;
        node2 = n2;
        Vec3 currSp = node2->position - node1->position;
        restLength = currSp.length();   //rest length
        Ks = k;
        Kd = k*0.02;
        C = 0.0;
        //printf("Spring initialized: Ks = %f, Kd = %f\n", Ks, Kd);
        // deltaP.setZeroVec();
        deltaP2.setZeroVec();
        dc_dp.setZeroVec();
        df_dx.setZeromat3();
        df_dv.setZeromat3();
        
        E.setIdentity();
       

	}

    Vec3 computeDcDp(Vec3 deltaP, float currentSpringLength) {
    const float epsilon = 1e-6; // Threshold for numerical stability
    if (currentSpringLength > epsilon) {
        Vec3 dc_dp = deltaP / currentSpringLength;
        // Clamp to remove numerical noise
        if (fabs(dc_dp.x) < epsilon) dc_dp.x = 0.0;
        if (fabs(dc_dp.y) < epsilon) dc_dp.y = 0.0;
        if (fabs(dc_dp.z) < epsilon) dc_dp.z = 0.0;
        return dc_dp;
    } else {
        printf("Warning: Degenerate spring detected with length: %f\n", currentSpringLength);
        return Vec3(0.0, 0.0, 0.0); // Handle edge case
    }
}


    void applyInternalForce(double timeStep) // Compute spring internal force
	{
        double currLen = Vec3::dist(node1->position, node2->position);

        //debugging
        // printf("Node1 Position: ");
        // node1->position.printVec3();
        // printf("Node2 Position: ");
        // node2->position.printVec3();
        // printf("Current Spring Length: %f\n", currLen);
        // printf("Computed dc_dp: ");
        // dc_dp.printVec3();
        // printf("Outer Product of dc_dp: ");
        // Mat3::outerProduct(dc_dp, dc_dp).printMat3();

       // dc_dp = (node2->position - node1->position)/currLen;   // orig: p1-p2
        // if (currLen < 1e-6) {
            // printf("Warning: Spring length is near zero, skipping dc_dp computation.\n");
            // dc_dp.setZeroVec();
        // }


        // printf("Computed dc_dp: ");
        // dc_dp.printVec3();

        //C = currLen-restLength; //see if restlength matches
        deltaV = node2->velocity - node1->velocity;
        Vec3 deltaP = node2->position - node1->position; //originally p1-p1
        //printf("deltaP: (%f, %f, %f)\n", deltaP.x, deltaP.y, deltaP.z);

        float currentSpringLength = deltaP.length(); // Compute magnitude
        dc_dp = computeDcDp(deltaP,currentSpringLength);
  
 
        C = currentSpringLength-restLength;
        Vec3 v1 = node1->velocity;
        Vec3 v2 = node2->velocity;
        Vec3 f = dc_dp * ((C)*Ks + Vec3::dot(deltaV, dc_dp)*Kd);
        C_Dot = Vec3::dot(v1,dc_dp) + Vec3::dot(v2,dc_dp*(-1.0));
        deltaP2 = Vec3(deltaP.x * deltaP.x, deltaP.y * deltaP.y, deltaP.z * deltaP.z);
        node1->addForce(f);
        node2->addForce(f.minus());
    }

    



    void springForceDerivative(){
        I.setIdentity();
        float c1 = C;
        d2C_dp1p1.identityMult ( (deltaP2 * (c1* (-1.0) ) ) + c1);
        d2C_dp1p2.identityMult((deltaP2 * c1) - c1);
        d2C_dp2p1 = d2C_dp1p2;
        d2C_dp2p2 = d2C_dp1p1;
        //df_dx = Mat3::outerProduct(dc_dp, dc_dp) * -Ks+ Mat3(1e-8);
       // df_dv = Mat3::outerProduct(dc_dp, dc_dp) * Kd+ Mat3(1e-8);
        //Mat3 df_dx = Mat3::outerProduct(dc_dp, dc_dp) * -Ks + Mat3(1e-8); // Adding small regularization
        
        //printf("Spring Deformation Derivative (dc_dp): ");
        //dc_dp.printVec3();
        //printf("Relative Velocity (deltaV): ");
        //deltaV.printVec3();

        Mat3 dp1 = Mat3::outerProduct(dc_dp, dc_dp);
        Mat3 dp2 = Mat3::outerProduct(dc_dp, dc_dp*(-1.0));
        Mat3 dp3 = Mat3::outerProduct(dc_dp*(-1.0),dc_dp*(-1.0));

        Mat3 term1 =  ( (dp1+(d2C_dp1p1 * c1)) * (-Ks)) - ((d2C_dp1p1 * C_Dot)*(Kd)); //replaced kd with -0.125
        Mat3 term2 =  ( (dp2+(d2C_dp1p2 * c1)) * (-Ks)) - ((d2C_dp1p2 * C_Dot)*(Kd));
        // Compute term1 and term2
        // printf("Spring Derivative Calculation:\n");
        // printf("Term1: ");
        // term1.printMat3();
        // printf("Term2: ");
        // term2.printMat3();
        
        df_dx = term1 + ( term2*2); // term2 is used twice
        df_dv = (dp1 + dp2 + dp3)*(Kd);

        // printf("deltaV: ");
        // deltaV.printVec3();
        // printf("Outer Product of dc_dp: ");
        // Mat3::outerProduct(dc_dp, dc_dp).printMat3();
        // printf("Kd: %f\n", Kd);
        // printf("Computed df_dv: ");
        // df_dv.printMat3();

        
        // printf("df_dx:\n");
        // df_dx.printMat3();
        // printf("df_dv:\n");
        // df_dv.printMat3();
        node1->addForceDerivative(df_dx,df_dv);
    }
};
