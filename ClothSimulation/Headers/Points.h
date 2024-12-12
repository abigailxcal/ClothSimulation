#pragma once

#include "Vectors.h"
#include <vector>
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
        // M = Mat3(mass);
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

    void printPosition()
    {
        printf("printing position ");
        printf("(%f, %f, %f)\n", position.x, position.y, position.z);
    }
    void printForce()
    {
        printf("printing force ");
        printf("(%f, %f, %f)\n", force.x, force.y, force.z);
    }
    void printVelocity()
    {
        printf("printing velocity ");
        printf("(%f, %f, %f)\n", velocity.x, velocity.y, velocity.z);
    }
    void printForceDerivatives()
    {
        printf("printing df_dx and df_dv: ");
        // printf("printing jacobian: \n ");
        // j.printMat3();
        df_dx.printMat3();
        df_dv.printMat3();
    }

    void implicit_integration(double deltaTime)
    {
        if (!isFixed)
        {
            float h = deltaTime;
            float y = 0.0;                                    // correciton term
            A = M - ((df_dv + (df_dx * h)) * h);              //  A = M - K*deltaTime
            b = ((force + (df_dx * (velocity + y) * h)) * h); // A was modified but b hasn't been
            printf("Dimension of A");
            A.printMat3();
            printf("Dimension of deltaV");
            A.printMat3();
            printf("Dimension of b");
            b.printVec3();

            solveConjugateGradient(A, deltaV, b);
            velocity += (deltaV * deltaTime);
            position += velocity * deltaTime;
            // printf("new pos: ");
            // position.printVec3();
            // printf("new velocity: ");
            // velocity.printVec3();
        }
        force.setZeroVec();
        df_dx.setZeromat3();
        df_dv.setZeromat3();
    }

    // to solve CGM Preconditioned, i think i need to create functions in cloth that populate the large matrices/vectors
    // A, b, P_, and P_inv with corresp values from each node in the cloth
    // Then, i call the preconditioned CGM in Cloth and distribute the new velocity and position to each node
    // Also, if i were to do that, i would need to convert Vec3/Mat3 into datatypes compatible with LargeVector<glm>
    // After CGM is solved, convert LargeVector<glm> back into Vec3/Mat3
    void solveConjugateGradient(Mat3 A, Vec3 &x, Vec3 b)
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
        if (i > 3)
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
};
