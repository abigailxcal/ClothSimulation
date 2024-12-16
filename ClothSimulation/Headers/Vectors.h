#pragma once

#include <math.h>

struct Vec2
{
    double x;
    double y;
    
    Vec2(void)
    {
        x = 0.0;
        y = 0.0;
    }
    Vec2(double x0, double y0)
    {
        x = x0;
        y = y0;
    }
    ~Vec2() {}
    
    Vec2 operator+(Vec2 v)
    {
        return Vec2(x+v.x, y+v.y);
    }
    Vec2 operator-(Vec2 v)
    {
        return Vec2(x-v.x, y-v.y);
    }

    void operator+=(Vec2 v)
    {
        x += v.x;
        y += v.y;
    }
    void operator-=(Vec2 v)
    {
        x -= v.x;
        y -= v.y;
    }
    
};

struct Vec3
{
    double x;
    double y;
    double z;
    
    Vec3(void)
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }
    Vec3(double x0, double y0, double z0)
    {
		x = x0;
		y = y0;
		z = z0;
	}
    Vec3(double val)
    {
		x = val;
		y = val;
		z = val;
	}
	~Vec3(){}
    
    static Vec3 cross(Vec3 v1, Vec3 v2)
    {
        return Vec3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
    }
    static double dot(Vec3 v1, Vec3 v2)
    {
        return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
    }
    static double dist(Vec3 v1, Vec3 v2)
    {
        return sqrt(pow(v1.x - v2.x, 2) + pow(v1.y - v2.y, 2) + pow(v1.z - v2.z, 2));
    }

    
    Vec3 minus()
    {
        return Vec3(-x, -y, -z);
    }
    
	Vec3 operator+(Vec3 v)
	{
		return Vec3(x+v.x, y+v.y, z+v.z);
	}

    Vec3 operator+(double v)
	{
		return Vec3(x+v, y+v, z+v);
	}
	Vec3 operator-(Vec3 v)
	{
		return Vec3(x-v.x, y-v.y, z-v.z);
	}

	void operator+=(Vec3 v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}
	void operator-=(Vec3 v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}
    
	Vec3 operator*(double n)
	{
		return Vec3(x*n, y*n, z*n);
	}
    Vec3 operator*(Vec3 v2)
	{
		return Vec3(x*v2.x, y*v2.y, z*v2.z);
	}
	Vec3 operator/(double n)
	{
		return Vec3(x/n, y/n, z/n);
	}
    
    bool operator ==(Vec3 &v)
    {
        return x == v.x && y == v.y&&z == v.z;
    }
    bool operator!=(Vec3 &v)
    {
        return x != v.x || y != v.y || z != v.z;
    }

    double length()
    {
        return sqrt(x*x + y*y + z*z);
    }
    
    void normalize()
    {
        double w = length();
        if (w < 0.00001) return;
        
        x /= w;
        y /= w;
        z /= w;
    }
    void printVec3(){
        printf("(%f,%f,%f)\n", x,y,z);
    }
	
	void setZeroVec(){
        x = 0.0;
        y = 0.0;
        z = 0.0;
	}
};

    struct Mat3
{
    double x1,y1,z1;
    double x2,y2,z2;
    double x3,y3,z3;
    
    Mat3(void)
    {
        x1 = 0.0, y1= 0.0,z1 = 0.0;
        x2 = 0.0, y2= 0.0,z2 = 0.0;
        x3 = 0.0, y3= 0.0,z3 = 0.0;
    }
    Mat3(float n)
    {
    x1 = n, y1= n,z1 = n;
    x2 = n, y2= n,z2 = n;
    x3 = n, y3= n,z3 = n;
    }
    Mat3(double x1_0,double y1_0,double z1_0,
         double x2_0,double y2_0,double z2_0,
         double x3_0,double y3_0,double z3_0)
    {
		x1 = x1_0, y1= y1_0,z1 = z1_0;
		x2 = x2_0, y2= y2_0,z2 = z2_0;
		x3 = x3_0, y3= y3_0,z3 = z3_0;
	}

    Mat3(Vec3 v1, Vec3 v2, Vec3 v3) {
        x1 = v1.x, y1= v1.y,z1 = v1.z;
        x2 = v2.x, y2= v2.y,z2 = v2.z;
        x3 = v3.x, y3= v3.y,z3 = v3.z;
    }

	~Mat3(){}

    static Mat3 outerProduct(Vec3 v1, Vec3 v2){
        return Mat3(v1.x*v2.x, v1.x*v2.y,v1.x*v2.z,
            v1.y*v2.x, v1.y*v2.y,v1.y*v2.z,
            v1.z*v2.x, v1.z*v2.y,v1.z*v2.z);
    }

    Mat3 operator*(double val){
        return Mat3(x1*val,y1*val,z1*val,
        x2*val,y2*val,z2*val,
        x3*val,y3*val,z3*val );  

    }
    Mat3 operator*(Mat3 other){
        return Mat3(
            (x1*other.x1+y1*other.x2+z1*other.x3) ,(x1*other.y1+y1*other.y2+z1*other.y3),(x1*other.z1+y1*other.z2+z1*other.z3),
            (x2*other.x1+y2*other.x2+z2*other.x3) ,(x2*other.y1+y2*other.y2+z2*other.y3),(x2*other.z1+y2*other.z2+z2*other.z3),
            (x3*other.x1+y3*other.x2+z3*other.x3) ,(x3*other.y1+y3*other.y2+z3*other.y3),(x3*other.z1+y3*other.z2+z3*other.z3)

        );
    }
    //edited operator mat*vec to return mat instead of vec
    Vec3 operator*(Vec3 v){
        return Vec3(x1*v.x + y1*v.y + z1*v.z,
                    x2*v.x + y2*v.y + z2*v.z,
                    x3*v.x + y3*v.y + z3*v.z );  
    } 
    Mat3 identityMult(Vec3 v) {
        return Mat3(
        x1 = v.x,  y1= 0.0,  z1 = 0.0,
        x2 = 0.0,  y2 = v.y, z2 = 0.0,
        x3 = 0.0,  y3= 0.0,  z3 = v.z);
        
    }
    // Mat3 operator*(Vec3 v) {
        // return Mat3(
            // (x1*v.x+ y1*other.x2+z1*other.x3) ,(x1*other.y1+y1*other.y2+z1*other.y3),(x1*other.z1+y1*other.z2+z1*other.z3)
            // (x2*v.x+ y2*other.x2+z2*other.x3) ,(x2*other.y1+y2*other.y2+z2*other.y3),(x2*other.z1+y2*other.z2+z2*other.z3)
            // (x3*v.x+ y3*other.x2+z3*other.x3) ,(x3*other.y1+y3*other.y2+z3*other.y3),(x3*other.z1+y3*other.z2+z3*other.z3)
        // )
    // }
    

    Mat3 operator-(Mat3 m2){
    return Mat3(x1-m2.x1,y1-m2.y1,z1-m2.z1,
                x2-m2.x2,y2-m2.y2,z2-m2.z2,
                x3-m2.x3,y3-m2.y3,z3-m2.z3 );  
    }
    void operator+=(Mat3 m2){
        x1+=m2.x1,y1+=m2.y1,z1+=m2.z1,
        x2+=m2.x2,y2+=m2.y2,z2+=m2.z2,
        x3+=m2.x3,y3+=m2.y3,z3+=m2.z3 ;
    }
    Mat3 operator+(Mat3 m2){
        return Mat3(
        x1+m2.x1,y1+m2.y1,z1+m2.z1,
        x2+m2.x2,y2+m2.y2,z2+m2.z2,
        x3+m2.x3,y3+m2.y3,z3+m2.z3) ;
    }
    void setZeromat3(){
        x1 = 0.0, y1= 0.0,z1 = 0.0;
        x2 = 0.0, y2= 0.0,z2 = 0.0;
        x3 = 0.0, y3= 0.0,z3 = 0.0;
    }
    void printMat3(){
        printf("(%f, %f, %f)\n(%f, %f, %f)\n(%f, %f, %f)\n", x1,y1,z1,x2,y2,z2,x3,y3,z3);
        }
    void setIdentity(){
        x1 = 1.0; y1 = 0.0; z1 = 0.0;
        x2 = 0.0; y2 = 1.0; z2 = 0.0;
        x3 = 0.0; y3 = 0.0; z3 = 1.0;
    }

 
};
