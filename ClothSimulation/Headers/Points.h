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
    double  mass;           // In this project it will always be 1
    bool    isFixed;        // Use to pin the cloth
    Vec2    texCoord;       // Texture coord
    Vec3    normal;         // For smoothly shading
	Vec3	position;
    Vec3    velocity;
    Vec3    force;
	Vec3	acceleration;
    Mat3    df_dx;
    Mat3    df_dv;

    //For implicit integration:
    Mat3 A;
    Vec3 b;
    Mat3 P_matrix;
    Mat3 P_inv;
    // Vec3    df_dx1; //might need to change the vector type
    // Vec3    df_dx2;
    // Vec3    df_dx3;

    // Vec3    df_dv1; 
    // Vec3    df_dv2;
    // Vec3    df_dv3;




public:
    Node(void) {
        mass = 1.0;
        isFixed = false;
        velocity.setZeroVec();
        force.setZeroVec();
        acceleration.setZeroVec();
    }
	Node(Vec3 pos)
    {
        mass = 1.0;
        isFixed = false;
        position = pos;
        velocity.setZeroVec();
        force.setZeroVec();
        acceleration.setZeroVec();
    }
	~Node(void) {}

	void addForce(Vec3 f)
	{
        force += f;
	}

    void addForceDerivative(Mat3 dx,Mat3 dv){
        df_dx += dx;
        df_dv += dv;

    }
    void implicit_integration(double timeStep){

    }
           

	void integrate(double timeStep) // Only non-fixed nodes take integration
	{
		if (!isFixed) // Verlet integration
		{
            acceleration = force/mass;
            velocity += acceleration*timeStep;
            position += velocity*timeStep;
        }
        force.setZeroVec();
	}
};
