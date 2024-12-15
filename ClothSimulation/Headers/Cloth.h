#pragma once
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <vector>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Spring.h"
#include "Rigid.h"
#include "large_vector.h"
#include <list>

 //Structural coeff: how well a cloth miantains its basic grid structure
 // high values = more rigid, low values = more stretch
 // rec: 200-500 N/m 
 //shearing coeff: how cloth miantains shape when stretched diagonally
 // low = unrealistic shearing
 // red: slightly lower than structural stiffness, 100-300 N/m
 // bending: resistance to out of plane bending
 // low = floppy, high = paper like 
 // 
 // damping coeff: controls the dissipation of energy to prevent endlesss oscillations
 // low damping = allow oscilllations to persist 

// The whole point of this program is to use implicit integration and CGM to solve the update equation
// for veloctiy at the current time step, which is then used to solve the update equation for position.
// the position vectors of each node at the current time step is what the cloth renderers use to 
// simulate the behavior. 
class Cloth
{
public:
    const int nodesDensity = 4; //4
    const int iterationFreq = 15; //25
    // const double structuralCoef = 1000;
    // const double shearCoef = 300;
    // const double bendingCoef = 20;
    // const double DEFAULT_DAMPING =  45.0;
    const double structuralCoef = 0.950;
    const double shearCoef = 0.15;
    const double bendingCoef = 0.02;
    const double DEFAULT_DAMPING =  0.7;    // maybe when time step exceeds 0.3, damping needs to be bigger than 0.5?
    LargeVector<glm::vec3> dV;
    LargeVector<glm::mat3> A;
    glm::mat3 M = glm::mat3(1.0f);
    LargeVector<glm::vec3> b;
    LargeVector<glm::vec3> P_;
    LargeVector<glm::vec3> P_inv;
    vector<float> inv_len;

// ------------- Logging --------------
    int targetNode = 10; 
    int targetSpring = 5;    // The spring to log data for
    bool logDerivatives = true; // Toggle for force derivatives
    std::list<int>targets;

    enum DrawModeEnum{
        DRAW_NODES,
        DRAW_LINES,
        DRAW_FACES
    };
    // DrawModeEnum drawMode = DRAW_FACES;
    DrawModeEnum drawMode = DRAW_LINES;
    
    Vec3 clothPos;
    
    int width, height;
    int nodesPerRow, nodesPerCol;
    
    std::vector<Node*> nodes;
	std::vector<Spring*> springs;
	std::vector<Node*> faces;
    
    Vec2 pin1;
    Vec2 pin2;
    
	Cloth(Vec3 pos, Vec2 size){
        clothPos = pos;
        width = size.x;
        height = size.y;
        init();
	}
	~Cloth(){ 
		for (int i = 0; i < nodes.size(); i++) { delete nodes[i]; }
		for (int i = 0; i < springs.size(); i++) { delete springs[i]; }
		nodes.clear();
		springs.clear();
		faces.clear();
	}
 
public:
    Node* getNode(int x, int y) { return nodes[y*nodesPerRow+x]; }
    Vec3 computeFaceNormal(Node* n1, Node* n2, Node* n3){
        return Vec3::cross(n2->position - n1->position, n3->position - n1->position);
    }
    
    void pin(Vec2 index, Vec3 offset){ // Pin cloth's (x, y) node with offset
        if (!(index.x < 0 || index.x >= nodesPerRow || index.y < 0 || index.y >= nodesPerCol)) {
            getNode(index.x, index.y)->position += offset;
            getNode(index.x, index.y)->isFixed = true;
        }
    }
    void unPin(Vec2 index){ // Unpin cloth's (x, y) node
        if (!(index.x < 0 || index.x >= nodesPerRow || index.y < 0 || index.y >= nodesPerCol)) {
            getNode(index.x, index.y)->isFixed = false;
        }
    }
    // Initial Cloth: add nodes, add springs to govern the behavior of those nodes, pin two corners/nodes of the cloth
	void init(){
        nodesPerRow = width * nodesDensity;
        nodesPerCol = height * nodesDensity;
        
        pin1 = Vec2(0, 0);
        pin2 = Vec2(nodesPerRow-1, 0);
        
        /** Add nodes **/
        printf("Init cloth with %d nodes\n", nodesPerRow*nodesPerCol);
        for (int i = 0; i < nodesPerRow; i ++) {
            for (int j = 0; j < nodesPerCol; j ++) {
                /** Create node by position **/
                Node* node = new Node(Vec3((double)j/nodesDensity, -((double)i/nodesDensity), 0));
                /** Set texture coordinates **/
                node->texCoord.x = (double)j/(nodesPerRow-1);
                node->texCoord.y = (double)i/(1-nodesPerCol);
                /** Add node to cloth **/
                nodes.push_back(node); 
            }
            std::cout << std::endl;
        }
        
        /** Add springs **/
        for (int i = 0; i < nodesPerRow; i ++) {
            for (int j = 0; j < nodesPerCol; j ++) {
                /** Structural **/
                if (i < nodesPerRow-1) springs.push_back(new Spring(getNode(i, j), getNode(i+1, j), structuralCoef));
                if (j < nodesPerCol-1) springs.push_back(new Spring(getNode(i, j), getNode(i, j+1), structuralCoef));
                /** Shear **/
                if (i < nodesPerRow-1 && j < nodesPerCol-1) {
                    springs.push_back(new Spring(getNode(i, j), getNode(i+1, j+1), shearCoef));
                    springs.push_back(new Spring(getNode(i+1, j), getNode(i, j+1), shearCoef));
                }
                /** Bending **/
                if (i < nodesPerRow-2) springs.push_back(new Spring(getNode(i, j), getNode(i+2, j), bendingCoef));
                if (j < nodesPerCol-2) springs.push_back(new Spring(getNode(i, j), getNode(i, j+2), bendingCoef));
            }
        }
        
        pin(pin1, Vec3(1.0, 0.0, 0.0));
        pin(pin2, Vec3(-1.0, 0.0, 0.0));
        
		/** Triangle faces **/
        for (int i = 0; i < nodesPerRow-1; i ++) {
            for (int j = 0; j < nodesPerCol-1; j ++) {
                // Left upper triangle
                faces.push_back(getNode(i+1, j));
                faces.push_back(getNode(i, j));
                faces.push_back(getNode(i, j+1));
                // Right bottom triangle
                faces.push_back(getNode(i+1, j+1));
                faces.push_back(getNode(i+1, j));
                faces.push_back(getNode(i, j+1));
            }
        }
        targets.push_back(5);
        targets.push_back(22);
        targets.push_back(193);
        
	}
    void createTargetNodes(){
        for (int idx: targets){
            nodes[idx]->istarget(idx);
        }
    }
	
	void computeNormal(){
        /** Reset nodes' normal **/
        Vec3 normal(0.0, 0.0, 0.0);
        for (int i = 0; i < nodes.size(); i ++) {
            nodes[i]->normal = normal;
        }
        /** Compute normal of each face **/
        for (int i = 0; i < faces.size()/3; i ++) { // 3 nodes in each face
            Node* n1 = faces[3*i+0];
            Node* n2 = faces[3*i+1];
            Node* n3 = faces[3*i+2];
            
            // Face normal
            normal = computeFaceNormal(n1, n2, n3);
            // Add all face normal
            n1->normal += normal;
            n2->normal += normal;
            n3->normal += normal;
        }
        
        for (int i = 0; i < nodes.size(); i ++) {
            nodes[i]->normal.normalize();
        }
	}
	//
	void addForce(Vec3 f){		 
		for (int i = 0; i < nodes.size(); i++){
			nodes[i]->addForce(f);
		}
	}

	void computeForce(double timeStep, Vec3 gravity){
        /** Nodes **/
		for (int i = 0; i < nodes.size(); i++){
            
            nodes[i]->addForce(gravity * nodes[i]->mass);
            nodes[i]->addForce (nodes[i]->velocity*(DEFAULT_DAMPING)*(-1.0));
            
            
		}
		/** Springs **/
		for (int i = 0; i < springs.size(); i++){
			springs[i]->applyInternalForce(timeStep);
		} 
        // for (int i = 0; i < nodes.size(); i++){
        //     if (i==targetNode){
        //         nodes[i]->istarget();
        //     //nodes[i]->resetForces();
        //     std::cout << "==================== Node" << i << " ====================\n";
        //     std::cout << "Position: (" << nodes[i]->position.x << ", " << nodes[i]->position.y << ", " << nodes[i]->position.z << ")\n";
        //     std::cout << "Internal Forces: (" << nodes[i]->force.x << ", " << nodes[i]->force.y << ", " << nodes[i]->force.z << ")\n";
        //     std::cout << "Velocity: (" << nodes[i]->velocity.x << ", " << nodes[i]->velocity.y << ", " << nodes[i]->velocity.z << ")\n";
        //     }
        // }  
	}
    // forces are calculated via springs, which get stored in the node
    // EDIT: to keep things consistent, cloth should call functions of nodes 
    // as much as possible. So change this so that node has a function to get it's 
    // force derivates from springs
    void computeForceDerivatives(double timeStep){
     for (int i = 0; i < springs.size(); i++){
            springs[i]->springForceDerivative(); //accesses df_dx,df_dv in Spring.h
            // if (i==targetSpring) {
            //     springs[i]->logDerivatives();

            // }
		}   
    }

    void implicit_integration_simple(double timeStep){
        for (int i = 0; i < nodes.size(); i++){
            
            nodes[i]->implicit_integration(timeStep);
        }     
    }

    
    Vec3 getWorldPos(Node* n) { return clothPos + n->position; }
    void setWorldPos(Node* n, Vec3 pos) { n->position = pos - clothPos; }
    
	void collisionResponse(Ground* ground)
    // void collisionResponse(Ground* ground, Ball* ball)
	{
        for (int i = 0; i < nodes.size(); i++)
        {
            /** Ground collision **/
            if (getWorldPos(nodes[i]).y < ground->position.y) {
                nodes[i]->position.y = ground->position.y - clothPos.y + 0.01;
                nodes[i]->velocity = nodes[i]->velocity * ground->friction;
            }
            
            /** Ball collision **/
            // Vec3 distVec = getWorldPos(nodes[i]) - ball->center;
            // double distLen = distVec.length();
            // double safeDist = ball->radius*1.05;
            // if (distLen < safeDist) {
            //     distVec.normalize();
            //     setWorldPos(nodes[i], distVec*safeDist+ball->center);
            //     nodes[i]->velocity = nodes[i]->velocity*ball->friction;
            // }
        }
	}
    // --------------- for preconditioned conjugate gradient method ----------------------
    void implicit_integration_cgm(double timeStep){
        A.resize(nodes.size());
        b.resize(nodes.size());
        dV.resize(nodes.size());
        P_.resize(nodes.size());
        P_inv.resize(nodes.size());
        if (A.size() != nodes.size()) {
    std::cerr << "Mismatch between A and nodes sizes!" << std::endl;
    return;
}
     for (int i = 0; i < nodes.size(); i++){
         if (nodes[i] && !nodes[i]->isFixed) {
        A[i] = nodes[i]->calculateA(timeStep);
        b[i] = nodes[i]->calculate_b(timeStep);
        P_[i] = nodes[i]->calculateP();
        P_inv[i] = nodes[i]->calculateP_inv();
        }
     }

     SolveConjugateGradientPreconditioned(A,dV,b,P_,P_inv);
     for(int i=0;i<nodes.size();i++) {
        if (nodes[i] && !nodes[i]->isFixed) {
		nodes[i]->apply_PCGM(timeStep,dV[i]);
        }
	}
 }
 void SolveConjugateGradientPreconditioned(LargeVector<glm::mat3> A, LargeVector<glm::vec3>& x, LargeVector<glm::vec3> b,LargeVector<glm::vec3> P, LargeVector<glm::vec3> P_inv) {
	float i =0;
    int i_max = 10;
    float EPS=0.001f;
    float EPS2 = EPS*EPS;
	LargeVector<glm::vec3> r =  (b - A*x);
	LargeVector<glm::vec3> d = P_inv*r;
	LargeVector<glm::vec3> q;
	float alpha_new = 0;
	float alpha = 0;
	float beta  = 0;
	float delta_old = 0;
	float delta_new = dot(r,P*r);
	float delta0    = delta_new;
	while(i<i_max && delta_new> EPS2*delta0) {
		q = A*d;
		alpha = delta_new/dot(d,q);
		x = x + alpha*d;
		r = r - alpha*q;
		delta_old = delta_new;
		delta_new = dot(r,r);
		beta = delta_new/delta_old;
		d = r + beta*d;
		i++;
	}
    if (i > 3)
        {
            printf("Converged at iteration %f \n", i);
        }

}

};
