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

    // actual spring properties
    Node *node1;
    Node *node2;
	double restLength;
    double Ks;
    double Kd;

    // force and derivative properties
    double C;         // constraint 
    double C_Dot;     // time derivative of constraint 
    Vec3 deltaP; 
    Vec3 deltaP2;     
    Vec3 deltaV;      
    Vec3 dc_dp;       // gradient of constraint wrt position
    Mat3 df_dx;       // derivative of force wrt position
    Mat3 df_dv;       // derivative of force wrt velocitie
    Mat3 J;
    Mat3 E;
    float inv_len;
    Mat3 I;
    Mat3 d2C_dp1p1;
    Mat3 d2C_dp1p2;
    Mat3 d2C_dp2p1;
    Mat3 d2C_dp2p2;
    
	Spring(Node *n1, Node *n2, double k){
        node1 = n1;
        node2 = n2;
        Vec3 currSp = node2->position - node1->position;
        restLength = currSp.length();   
        Ks = k;
        Kd = 0.25;
        C = 0.0;
        deltaP2 = Vec3(0.0);
        dc_dp = Vec3(0.0);
        df_dx = Mat3(0.0);
        df_dv = Mat3(0.0);
        
       

	}
    // might remove this for simplicity's sake
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

    void applyInternalForce(double timeStep) 
	{
        Vec3 v1 = node1->velocity;
        Vec3 v2 = node2->velocity;

        deltaV = v2-v1;
        deltaP = node2->position - node1->position; //originally p1-p1
        float currentSpringLength = deltaP.length(); // redundant i think 
        dc_dp = computeDcDp(deltaP,currentSpringLength);
        C = currentSpringLength-restLength;
       
        Vec3 f = dc_dp * ((C)*Ks + Vec3::dot(deltaV, dc_dp)*Kd);
        C_Dot = Vec3::dot(v1,dc_dp) + Vec3::dot(v2,dc_dp*(-1.0));  //  might need to flip the negative signs
        deltaP2 = Vec3(deltaP.x * deltaP.x, deltaP.y * deltaP.y, deltaP.z * deltaP.z);

        node1->addForce(f);
        node2->addForce(f.minus());
    }

    void springForceDerivative(){
        float c1 = C;
        d2C_dp1p1.identityMult ( (deltaP2 * (c1* (-1.0) ) ) + c1);
        d2C_dp1p2.identityMult((deltaP2 * c1) - c1);

        // for symmetry, only useful for preconditioned CGM
        // uhm ok i guess not
        d2C_dp2p1 = d2C_dp1p2; 
        d2C_dp2p2 = d2C_dp1p1;

        Mat3 dp1 = Mat3::outerProduct(dc_dp, dc_dp);
        Mat3 dp2 = Mat3::outerProduct(dc_dp, dc_dp*(-1.0));
        Mat3 dp3 = Mat3::outerProduct(dc_dp*(-1.0),dc_dp*(-1.0));

        Mat3 term1 =  ( (dp1+(d2C_dp1p1 * c1)) * (-Ks)) - ((d2C_dp1p1 * C_Dot)*(Kd)); //replaced kd with -0.125
        Mat3 term2 =  ( (dp2+(d2C_dp1p2 * c1)) * (-Ks)) - ((d2C_dp1p2 * C_Dot)*(Kd));
        Mat3 term3 = ( (dp2+(d2C_dp1p1 * c1)) * (-Ks)) - ((d2C_dp1p1 * C_Dot)*(Kd));
        df_dx = term1 + term2 + term3; 
        df_dv = (dp1 + dp2 + dp3)*(-Kd);  // should it be -kd or +kd?
        node1->addForceDerivative(df_dx,df_dv);
    }

    void logDerivatives() {
        std::cout << "df_dx: \n";
        df_dx.printMat3();
        std::cout << "df_dv: \n";
        df_dv.printMat3();

    }
};
