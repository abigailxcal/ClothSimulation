#pragma once

#include "Vectors.h"
#include <vector>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


struct Vertex
{
public:
    Vec3 position;
    Vec3 normal;

    Vertex() {}
    Vertex(Vec3 pos)
    {
        position = pos;
    }
    ~Vertex() {}
};

class Node
{
public:
    double mass;   // In this project it will always be 1
    bool isFixed;  // Use to pin the cloth
    Vec2 texCoord; // Texture coord
    Vec3 normal;   // For smoothly shading
    Vec3 position;
    Vec3 velocity;
    Vec3 force;
    Vec3 acceleration;

    // For implicit integration + CGM :
    Mat3 df_dx;
    Mat3 df_dv;
    Vec3 deltaV;
    Mat3 j;
    Mat3 A;
    Vec3 b;
    Mat3 P;
    Mat3 P_inv;
    Vec3 C;
    Vec3 dc_dp;
    Vec3 C_dot;
    Mat3 M;

public:
    Node(void)
    {
        mass = 1.0;
        isFixed = false;
        velocity.setZeroVec();
        force.setZeroVec();
        acceleration.setZeroVec();
        M.setIdentity();
        df_dx.setZeromat3();
        df_dv.setZeromat3();
        deltaV.setZeroVec();
        j.setZeromat3();
    }
    Node(Vec3 pos)
    {
        mass = 1.0;
        isFixed = false;
        position = pos;
        velocity.setZeroVec();
        force.setZeroVec();
        acceleration.setZeroVec();
        M.setIdentity();
        df_dx.setZeromat3();
        df_dv.setZeromat3();
        deltaV.setZeroVec();
    }
    ~Node(void) {}

    void addForce(Vec3 f)
    {
        force += f;
    }
    void addForce(Vec3 f, Vec3 c, Vec3 c_dot)
    {
        force += f;
    }

    void addForceDerivative(Mat3 dx, Mat3 dv)
    {
        df_dx += dx;
        df_dv += dv;
    }
    void addForceDerivative(Mat3 jacobian)
    {
        j += jacobian;
    }

    void resetForces()
    {
        force.setZeroVec();
        df_dx.setZeromat3();
        df_dv.setZeromat3();
    }

    /*---------- CGM3 Implementation: solves dv individually for each point */
    void implicit_integration(double deltaTime){
        if (!isFixed){
            float h = deltaTime;
            float y = 0.0; //correciton term
            A =  M - ((df_dv + (df_dx * h))*h);
            b = ((force + (df_dx * ( velocity + y) * h)) * h); //A was modified by b hasn't been 
            
            //auto start = std::chrono::high_resolution_clock::now();
            solveConjugateGradient(A, deltaV, b, deltaTime);
            //auto end = std::chrono::high_resolution_clock::now();
            //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            //printf("CGM Convergence Time: %lld ms\n", duration.count());
            velocity += (deltaV*deltaTime);
            position += velocity*deltaTime;
        }
        resetForces();
    }

    void solveConjugateGradient(Mat3 A, Vec3 &x, Vec3 b,double timeStep)
    {
        int i_max = 10;
        float EPS = 0.001f;
        float EPS2 = EPS * EPS;
        float i = 0;
        Vec3 r = b - A * x; // residual
        Vec3 d = r;
        Vec3 q = Vec3(0);
        float beta = 0;
        float delta_old = 0;
        float delta_new = Vec3::dot(r, r);
        float delta0 = delta_new;
        float alpha = delta_new / (Vec3::dot(d, q) + 1e-8f);
        while (i < i_max && delta_new > EPS2 * delta0)
        {
            q = A * d;
            alpha = delta_new / Vec3::dot(d, q);
            x = x + (d * alpha); // updates position
            r = r - (q * alpha);
            delta_old = delta_new;
            delta_new = Vec3::dot(r, r);
            beta = delta_new / delta_old;
            d = r + (d * beta);
            i++;
        }
        //logIterations(timeStep, i);
        if (i >2)
        {
            printf("Converged at iteration %f \n", i);
        }
    }


     // --------------- for preconditioned conjugate gradient method ----------------------
    glm::mat3 calculateA(double deltaTime)
    {
        float h = deltaTime;
        float y = 0.0;                    
        A = M - ((df_dv + (df_dx * h)) * h); //  A = M - K*deltaTime
        return glm::mat3(A.x1, A.y1, A.z1,
                         A.x2, A.y2, A.z2,
                         A.x3, A.y3, A.z3);
    }

    glm::vec3 calculate_b(double deltaTime)
    {
        glm::vec3 b_glm;
        float h = deltaTime;
        float y = 0.0; // correciton term
        b = ((force + (df_dx * (velocity + y) * h)) * h);
        b_glm = glm::vec3(b.x, b.y, b.z);
        return b_glm;
    }

    glm::vec3 calculateP()
    {
        return glm::vec3(A.x1, A.y2, A.z3);
    }
    glm::vec3 calculateP_inv()
    {
        float x_inv = A.x1 != 0.0f ? 1.0f / A.x1 : 1e-8f;
        float y_inv = A.y2 != 0.0f ? 1.0f / A.y2 : 1e-8f;
        float z_inv = A.z3 != 0.0f ? 1.0f / A.z3 : 1e-8f;
        return glm::vec3(x_inv, y_inv, z_inv);
    }

    void apply_PCGM(double deltaTime, glm::vec3 dv)
    {
        if (glm::any(glm::isnan(dv)))
        {
            std::cerr << "Invalid dv detected!" << std::endl;
            return;
        }
        Vec3 dV = Vec3(dv[0], dv[1], dv[2]);
        velocity += (dV * deltaTime);
        position += velocity * deltaTime;
        resetForces();
    }

    void logIterations(double timeStep, int iterations) {
    std::ofstream logFile("iterations.log", std::ios::app); // Open file in append mode
    if (logFile.is_open()) {
        logFile << timeStep << " " << iterations << "\n";   // Log time step and iterations
        logFile.close();                                   // Close the file
    } else {
        std::cerr << "Error: Could not open iterations.log for writing." << std::endl;
    }
}
};
