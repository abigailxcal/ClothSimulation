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
    Vec3    deltaV;
    Mat3 j;
    

    //For implicit integration:
    Mat3 A;
    Vec3 b;
    Mat3 P_matrix;
    Mat3 P_inv;
    Vec3 C;
    Vec3 dc_dp;
    Vec3 C_dot;
    Mat3 M;

    // for Conjugate Gradient
   
   
   




public:
    Node(void) {
        mass = 1.0;
        isFixed = false;
        velocity.setZeroVec();
        force.setZeroVec();
        acceleration.setZeroVec();
        M.setIdentity();
        //M = Mat3(mass);
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
        //M = Mat3(mass);
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

    void addForceDerivative(Mat3 dx,Mat3 dv){
        df_dx += dx;
        df_dv += dv;
        }
        void addForceDerivative(Mat3 jacobian){
            j += jacobian;
            
            //printPosition();
            //printForce();
            //printForceDerivatives();
            // Debug info
                // printf("df_dx: (%f, %f, %f)  \n", df_dx.x1, df_dx.y1, df_dx.z1);
        }
                //  std::cout << "Failed to create GLFW window." << std::endl;
        
    void printPosition(){
        printf("printing position ");
        printf("(%f, %f, %f)\n" ,position.x, position.y, position.z);
    }
    void printForce(){
        printf("printing force ");
        printf("(%f, %f, %f)\n" ,force.x, force.y, force.z);
    }
    void printForceDerivatives(){
        printf("printing df_dx and df_dv: ");
        // printf("printing jacobian: \n ");
        // j.printMat3();
        df_dx.printMat3();
        df_dv.printMat3();
    }

    void implicit_integration(double deltaTime){
        if (!isFixed){
            float h = deltaTime;
            float y = 0.0; //correciton term
            A =  M - ((df_dv + (df_dx * h))*h);
            b = ((force + (df_dx * ( velocity + y) * h)) * h); //A was modified by b hasn't been 
            //printf("old pos: ");
            //position.printVec3();
            //printf("old velocity: ");
            //velocity.printVec3();
            // Debugging prints
    // printf("Mass matrix M:\n");
    // M.printMat3();
    // printf("Force derivative df_dx:\n");
    // df_dx.printMat3();
    // printf("Damping derivative df_dv:\n");
    // df_dv.printMat3();
    // printf("Resulting matrix A:\n");
    // A.printMat3();
    //         printf("Matrix A:\n");
    //         A.printMat3();
    //         printf("Vector b: ");
    //         b.printVec3();
    
            solveConjugateGradient(A, deltaV, b);
            velocity += (deltaV*deltaTime);
            position += velocity*deltaTime;
            //printf("new pos: ");
            //position.printVec3();
            //printf("new velocity: ");
            //velocity.printVec3();
          
        }
       
        force.setZeroVec();
        df_dx.setZeromat3();
        df_dv.setZeromat3();


    }
    void solveConjugateGradient(Mat3 A, Vec3 &x, Vec3 b){
        int i_max = 10;
        float EPS=0.001f;
        float EPS2 = EPS*EPS;
    	float i =0;
	    Vec3 r = b - A*x; //residual
	    Vec3 d = r;
	    Vec3 q = Vec3(0);
	    float beta  = 0;
	    float delta_old = 0;
	    float delta_new = Vec3::dot(r,r);
	    float delta0    = delta_new;
        float alpha = delta_new / (Vec3::dot(d, q) + 1e-8f);
        // printf("Initial b: "); b.printVec3();
        // printf("Initial x: "); x.printVec3();
        // printf("Initial A:\n"); A.printMat3();
        // printf("old delta v: " ); 
        //x.printVec3();
	    while(i<i_max && delta_new> EPS2*delta0) {
	    	q = A*d;
	    	alpha = delta_new/Vec3::dot(d,q);
	    	x = x + (d*alpha); //updates position
	    	r = r - (q*alpha);
	    	delta_old = delta_new;
	    	delta_new = Vec3::dot(r,r);
	    	beta = delta_new/delta_old;
	    	d = r + (d*beta);
	    	i++;
            // printf("Iteration %d:\n", (int)i);
            // printf("Residual r: "); r.printVec3();
            // printf("Search direction d: "); d.printVec3();
            // printf("Alpha: %f\n", alpha);
            // printf("Delta new: %f\n", delta_new);
	    } 
        if (i>3){
        printf("Converged at iteration %f \n", i); 
        }
        // printf("new deltaV: " ); 
       // x.printVec3();
    }
	// }
};
    
	