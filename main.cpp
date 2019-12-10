#ifdef WIN32
#include <windows.h>
#endif

#if defined (__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
//#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#else //linux
#include <GL/gl.h>
#include <GL/glut.h>
#endif

//other includes
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <iomanip>
#include <ctime>

#define PI 3.141592654

using namespace std;

/*----------------------------- Global Variables -----------------------------*/
/*set in main()*/
//the number of pixels in the grid
int grid_width;
int grid_height;

//the size of pixels sets the inital window height and width
//don't make the pixels too large or the screen size will be larger than
//your display size
float pixel_size;

/*Window information*/
int win_height;
int win_width;

/* refraction coefficient */
//float refrCoef = 0.667;
float refrCoef = 1.5;
/* maximum bumber of recursion levels */
int recursionLevel = 3;

/*---------------------------------- Vector ----------------------------------*/
struct Vector {
    float x,y,z;
    Vector(float x1, float y1, float z1) {
        x = x1; y = y1; z = z1;
    }
    Vector(void) {
        x = 0; y = 0; z = 0;
    }
    inline Vector operator + (const Vector& v1)const {
        return Vector(x + v1.x, y + v1.y, z + v1.z);
    }
    inline Vector operator - (void)const {
        return Vector(-x, -y, -z);
    }
    inline Vector operator - (const Vector& v1)const {
        return Vector(x - v1.x, y - v1.y, z - v1.z);
    }
    inline float operator * (const Vector& v1)const {
        return x * v1.x + y * v1.y + z * v1.z;
    }
    inline Vector operator * (const float f)const {
        return Vector(x*f, y*f, z*f);
    }
    inline Vector operator / (const float f)const {
        return Vector(x/f, y/f, z/f);
    }
    inline Vector& operator += (const Vector& v1) {
        x += v1.x; y += v1.y; z += v1.z;
        return *this;
    };
//    inline Vector& operator -= (const Vector& v1) {
//        x -= v1.x; y -= v1.y; z -= v1.z;
//        return *this;
//    };
//    inline Vector& operator *= (float f) {
//        x *= f; y *= f; z *= f;
//        return *this;
//    };
//    inline Vector& operator /= (float f) {
//        x /= f; y /= f; z /= f;
//        return *this;
//    };
    inline float length(void)const {
        return sqrt(x*x + y*y + z*z);
    }
    inline Vector& norm(void) {
        float length = sqrt(x*x + y*y + z*z);
        x /= length; y /= length; z /= length;
        return *this;
    }
};

/*----------------------------------- Point ----------------------------------*/
struct Point {
    float x,y,z;
    Point(float x1, float y1, float z1) {
        x = x1; y = y1; z = z1;
    }
    Point(void) {
        x = 0; y = 0; z = 0;
    }
    inline Vector operator - (const Point& p1)const {
        return Vector(x - p1.x, y - p1.y, z - p1.z);
    }
    inline float operator * (const Point& p1)const {
        return x * p1.x + y * p1.y + z * p1.z;
    }
    inline Point operator * (const float f)const {
        return Point(x*f, y*f, z*f);
    }
    Vector toVec(void)const {
        return Vector(x,y,z);
    }
};

Vector toVec(Point& p) {
    return Vector(p.x, p.y, p.z);
}

Point toPnt(Vector& v) {
    return Point(v.x, v.y, v.z);
}

float dot(const Vector& v1, const Vector& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vector cross(const Vector& v1, const Vector& v2) {
    return Vector(v1.y*v2.z - v1.z*v2.y,v1.z*v2.x - v1.x*v2.z,v1.x*v2.y - v1.y*v2.x);
}

/*----------------------------------- Color ----------------------------------*/
struct Color {
    float r,g,b;
    Color(float r1, float g1, float b1) {
        r = r1; g = g1; b = b1;
    }
    Color(void) {
        r = 0; g = 0; b = 0;
    }
    inline Color operator + (const Color& c)const {
        return Color(r + c.r, g + c.g, b + c.b);
    }
//    inline Color operator - (const Color& c)const {
//        return Color(r - c.r, g - c.g, b - c.b);
//    }
    inline Color operator * (const Vector& v)const {
        return Color(r * v.x, g * v.y, b * v.z);
    }
    inline Color operator * (const Color& c)const {
        return Color(r * c.r, g * c.g, b * c.b);
    }
    inline Color operator * (const float f)const {
        return Color(r * f, g * f, b * f);
    }
    inline Color operator / (const float f)const {
        return Color(r / f, g / f, b / f);
    }
    inline Color &operator += (const Color& c) {
        r += c.r; g += c.g; b += c.b;
        return *this;
    };
//    inline Color &operator -= (const Color& c) {
//        r -= c.r; g -= c.g; b -= c.b;
//        return *this;
//    };
//    inline Color &operator *= (const Vector& c) {
//        r *= c.x; g *= c.y; b *= c.z;
//        return *this;
//    };
//    inline Color &operator *= (float f) {
//        r *= f; g *= f; b *= f;
//        return *this;
//    };
//    inline Color &operator /= (float f) {
//        r /= f; g /= f; b /= f;
//        return *this;
//    };
};
/* phong variables */
Point x = Point(0,5,0);     // lightSource
float K = 10;               // bigK
float IA = 0.2;
float IL = 1;
Color ka = Color(1.0,1.0,1.0);
Color kd = Color(1.0,1.0,1.0);
Color ks = Color(1.0,1.0,1.0);
int specularity = 1;        // bigN

/* global illumination */
Color kr = Color(0.5,0.5,0.5);
Color kt = Color(0.5,0.5,0.5);

/* store info of and pixels' color*/
vector<Color> pixels;

/*------------------------------------ Ray -----------------------------------*/
/* camera ray */
Point from = Point(0.0,-3.0,0.0);
Point at   = Point(0.0,0.0,0.0);
Vector up   = Vector(0.0,0.0,1.0);
Vector hori = Vector(1.0,0.0,0.0);
//float fov = 0.7853; // in radian
float fov = 45;

struct Ray {
    Point from;        // starting location
    Vector dir;         // unit direction vector
    bool inObj;         // if line is in an object or not
    Ray(const Point& from1, const Vector& dir1, bool inObj1) {
        from = from1;
        dir = dir1;
        inObj = inObj1;
    }
    Ray reflectRay(const Point& intersect, const Vector& N)const {
        Vector u = dir;
        Vector reflectVector = u - N * 2 * (u * N);
        return Ray(intersect, reflectVector.norm(), inObj);
    }
    Ray refractRay(const Point& intersect, const Vector& N)const {
        Vector u = dir;
        
//        float coef = refrCoef;
//        if (inObj == true) {
//            coef = 1 / coef;
//        }
        float coef = inObj ? refrCoef : 1 / refrCoef;
        
        float cosi = -N * u;
        float cosr = pow(1 - pow(coef,2) * (1 - pow(cosi,2)),0.5);
        Vector refractVector = u * coef - N * (cosr - coef * cosi);
        return Ray(intersect, refractVector.norm(), !inObj);
    }
};

Ray cameraRay(int x, int y) {
    Vector xv = hori * tan(fov / 180 * PI) * (2.0 * x / win_width - 1);
    Vector yv = up * tan(fov / 180 * PI) * (2.0 * y / win_height - 1);
    Vector fromToPix = at.toVec() + xv + yv + (at-from).norm();
    return Ray(from, fromToPix.norm(), false);
}

/*---------------------------------- Polygon ---------------------------------*/
struct Polygon {
    vector<Point> vertices;     // list of vertices
    Vector normal;              // normal vector of Polygon
    Color color;                // color of Polygon
    bool isTrans;               // is transparent
    Polygon(const vector<Point>& vertices1, const Vector& normal1,
          const Color& color1, bool isTrans1) {
        /* assign vertices */
        for (int i = 0; i < (int)vertices1.size(); i++) {
            vertices.push_back(vertices1.at(i));
        }
        
        /* calculate normal vector */
//        getNormal();
        normal = normal1;
        
        /* assign color and is Trans */
        color = color1;
        isTrans = isTrans1;
    }
    
    void getNormal(void) {
        /* calculate normal vector */
        Point v1 = vertices.at(0);
        Point v2 = vertices.at(1);
        Point v3 = vertices.at(2);
        
        Vector U = v2 - v1;
        Vector V = v3 - v1;
        
        float nx = (U.y * V.z) - (U.z * V.y);
        float ny = (U.z * V.x) - (U.x * V.z);
        float nz = (U.x * V.y) - (U.y * V.x);
        
        normal = Vector(nx,ny,nz);
    }
    
    bool intersect(const Ray& ray, Point& intersectPoint) {
        float D = -(normal * vertices.at(0).toVec());
        Vector P0 = ray.from.toVec();
        Vector u = ray.dir;
        /* distance between ray start point and intersect */
        float s = -(D + normal * P0) / (normal * u);
        /* if the ray is heading to the polygon's direction. */
        if (s >= 1e-4) {
            /* intersect between the plane of the polygon and the ray */
            Vector temp = P0 + u * s;
            Point P = toPnt(temp);
            for (int i = 0; i < (vertices.size()-1); i++) {
                /* get local cross product */
                Vector local = cross(vertices.at(i) - P, vertices.at(i+1) - P);
                if (local * normal <= 0) {
                    return false;
                }
            }
            if (cross(vertices.at((int)vertices.size()-1) - P, vertices.at(0) - P) * normal <= 0) {
                return false;
            }
            intersectPoint = P;
            return true;
        } else {
            return false;
        }
    }
};

/* related globals */
vector<Polygon> polygons;

/*---------------------------------- Sphere ----------------------------------*/
struct Sphere {
    Point center;
    float radius;
    Color color;
    bool isTrans;
    Sphere(const Point& center1, float radius1, const Color& color1, bool isTrans1) {
        center = center1;
        radius = radius1;
        color = color1;
        isTrans = isTrans1;
    }
    Vector getNormal(const Point& p) {
        return (p - center).norm();
    }
    bool intersect(const Ray& ray, Point& intersect) {
        Vector deltaP = center - ray.from;
        float udeltaP = ray.dir * deltaP;
        float sqrt1 = pow(udeltaP,2) - pow(deltaP.length(),2) + pow(radius,2);
        if (sqrt1 < 0) {
            return false;
        }
        sqrt1 = sqrt(sqrt1);

        Vector p1 = ray.from.toVec() + ray.dir * (udeltaP - sqrt1);
        Vector p2 = ray.from.toVec() + ray.dir * (udeltaP + sqrt1);

        bool v1 = ((p1 - ray.from.toVec()).length() >= 1e-4) && (udeltaP - sqrt1 >= 1e-4);
        bool v2 = ((p2 - ray.from.toVec()).length() >= 1e-4) && (udeltaP + sqrt1 >= 1e-4);
        if (v1 && v2) {
            intersect = (p1 - ray.from.toVec()).length() < (p2 - ray.from.toVec()).length() ? toPnt(p1) : toPnt(p2);
            return true;
        } else if (v1) {
            intersect = toPnt(p1);
            return true;
        } else if (v2) {
            intersect = toPnt(p2);
            return true;
        } else {
            return false;
        }
    }
};

/* related globals */
vector<Sphere> spheres;

/*--------------------------------- Ellipsoid --------------------------------*/
struct Ellipsoid {
    float A, B, C, D, E, F, G, H, I, J;
    Color color;
    bool isTrans;
    Ellipsoid(const Point& center, float radius, const Color& color1, bool isTrans1) {
        A = B = C = 1;
        D = E = F = 0;
        G = - 2 * center.x;
        H = - 2 * center.y;
        I = - 2 * center.z;
        J = center * center - radius * radius;
        
        color = color1;
        isTrans = isTrans1;
    }
    Ellipsoid(const Point& center, float a, float b, float c, const Color& color1, bool isTrans1) {
        A = b * b + c * c;
        B = a * a + c * c;
        C = a * a + b * b;
        D = E = F = 0;
        G = b * b * c * c * -2 * center.x;
        H = a * a * c * c * -2 * center.y;
        I = a * a * b * b * -2 * center.z;
        J = b * b * c * c * center.x * center.x
            + a * a * c * c * center.y * center.y
            + a * a * b * b * center.z * center.z
            - a * a * b * b * c * c;
        
        color = color1;
        isTrans = isTrans1;
    }
    Vector getNormal(const Point& p) {
        Vector ret;
        ret.x = 2*A*p.x + D*p.y + E*p.z + G;
        ret.y = 2*B*p.y + D*p.x + F*p.z + H;
        ret.z = 2*C*p.z + E*p.x + F*p.y + I;
        ret.norm();
        return ret;
    }
    bool intersect(const Ray& ray, Point& intersect) {
        float Aq = A * ray.dir.x * ray.dir.x + B * ray.dir.y * ray.dir.y + C * ray.dir.z * ray.dir.z
                    + D * ray.dir.x * ray.dir.y + E * ray.dir.x * ray.dir.z + F * ray.dir.y * ray.dir.z;
        float Bq = 2*A * ray.from.x * ray.dir.x + 2*B * ray.from.y * ray.dir.y + 2*C * ray.from.z * ray.dir.z
                    + D * (ray.from.x * ray.dir.y + ray.from.y * ray.dir.x)
                    + E * (ray.from.x * ray.dir.z + ray.from.z * ray.dir.x)
                    + F * (ray.from.y * ray.dir.z + ray.from.z * ray.dir.y)
                    + G * ray.dir.x + H * ray.dir.y + I * ray.dir.z;
        float Cq = A * ray.from.x * ray.from.x + B * ray.from.y * ray.from.y + C * ray.from.z * ray.from.z
                    + D * ray.from.x * ray.from.y + E * ray.from.x * ray.from.z + F * ray.from.y * ray.from.z
                    + G * ray.from.x + H * ray.from.y + I * ray.from.z + J;
        float delta = Bq * Bq - 4 * Aq * Cq;
        if (delta < 0) {
            return false;
        }
        delta = sqrt(delta);
        Vector p1, p2;
        if (Aq == 0) {
            p1 = p2 = ray.from.toVec() + ray.dir * (-Cq / Bq);
        }
        else {
            p1 = ray.from.toVec() + ray.dir * (- Bq - delta) / 2 / Aq;
            p2 = ray.from.toVec() + ray.dir * (- Bq + delta) / 2 / Aq;
        }
        bool v1 = ((p1 - ray.from.toVec()).length() >= 1e-4) && (- Bq - delta >= 1e-4);
        bool v2 = ((p2 - ray.from.toVec()).length() >= 1e-4) && (- Bq + delta >= 1e-4);
        if (v1 && v2)
            intersect = (p1 - ray.from.toVec()).length() < (p2 - ray.from.toVec()).length() ? toPnt(p1) : toPnt(p2);
        else if (v1)
            intersect = toPnt(p1);
        else if (v2)
            intersect = toPnt(p2);
        else return false;
        return true;
    }
};

/* related globals */
vector<Ellipsoid> ellipsoids;

bool intersectObject(const Ray& ray, Point& intersectPoint, int& id) {
//    Point resultPoint;
//    int resultID = -1;
//    float resultDist = 0;
    
    /* check if intersect with a polygon */
//    for (int i = 0; i < (int)polygons.size(); i++) {
//        Point curPoint;
//        if (polygons.at(i).intersect(ray, curPoint)) {
//            float curDist = (curPoint - ray.from).length();
//            if (resultID == -1) {
//                resultDist = curDist;
//                resultID = i;
//                resultPoint = curPoint;
//            } else if (resultID != -1 && resultDist > curDist) {
//                resultDist = curDist;
//                resultID = i;
//                resultPoint = curPoint;
//            }
//        }
//    }
    
    /* check if intersect with a sphere */
//    for (int i = 0; i < (int)spheres.size(); i++) {
//        Point curPoint;
//        if (spheres.at(i).intersect(ray, curPoint)) {
//            float curDist = (curPoint - ray.from).length();
//            if (resultID == -1) {
//                resultDist = curDist;
//                resultID = i + (int)polygons.size();
//                resultPoint = curPoint;
//            } else if (resultID != -1 && resultDist > curDist) {
//                resultDist = curDist;
//                resultID = i + (int)polygons.size();
//                resultPoint = curPoint;
//            }
//        }
//    }
    
    /* check if intersect with an ellipsoid */
//    for (int i = 0; i < (int)ellipsoids.size(); i++) {
//        Point curPoint;
//        if (ellipsoids.at(i).intersect(ray, curPoint)) {
//            float curDist = (curPoint - ray.from).length();
//            if (resultID == -1) {
//                resultDist = curDist;
//                resultID = i + (int)polygons.size();
//                resultPoint = curPoint;
//            } else if (resultID != -1 && resultDist > curDist) {
//                resultDist = curDist;
//                resultID = i + (int)polygons.size();
//                resultPoint = curPoint;
//            }
//        }
//    }

//    if (resultID != -1) {
//        id = resultID;
//        intersectPoint = resultPoint;
//        return true;
//    }
    Point best, now;
    int bestId = -1;
    float bestDist = 0, nowDist;
    bool ret = false;
    
    for (int i = 0; i < polygons.size(); i++) {
        if (polygons.at(i).intersect(ray, now)) {
            nowDist = (now - ray.from).length();
            if (!ret || (nowDist < bestDist)) {
                ret = true;
                best = now;
                bestId = i;
                bestDist = nowDist;
            }
        }
    }
    
    for (int i = 0; i < ellipsoids.size(); i++) {
        if (ellipsoids.at(i).intersect(ray, now)) {
            nowDist = (now - ray.from).length();
            if (!ret || (nowDist < bestDist)) {
                ret = true;
                best = now;
                bestId = i + (int)polygons.size();
                bestDist = nowDist;
            }
        }
    }
    if (ret) {
        intersectPoint = best;
        id = bestId;
    }
    return ret;
}

Color phong(const Point& point, const Point& f, const Vector& normal) {
    Color intensity = ka * IA;
    Vector l = (x - point).norm();
    Vector v = (f - point).norm();
    Vector r = (normal * 2.0 * (normal * l) - l).norm();
    Ray ray = Ray(point, l, false);
    
    Point intersectPoint;
    int id;
    if (intersectObject(ray, intersectPoint, id)) {
        return intensity;
    }
    
    Color diff;
    if (l * normal > 0) {
        diff = kd * (l * normal);
    }
    
    Color spec;
    if (l * normal > 0 && r * v > 0) {
        spec = ks * (r * v);
    }
    
    intensity += (diff + spec) * IL / ((f - point).length() + K);
    return intensity;
}

Color light(const Ray& ray, int recurlvl) {
    Color result;
    
    // base case
    if (recurlvl == 0) {
        return result;
    }
    
    // general case
    Point intersect;
    int id;
    if (intersectObject(ray, intersect, id)) {
        /* if intersecting object is a polygon */
        if (id < polygons.size()) {
            Polygon polygon = polygons.at(id);
            Vector normal = polygon.normal;
            result = polygon.color * phong(intersect, ray.from, normal);
            Ray refl = ray.reflectRay(intersect, normal);
            Color IR = light(refl, recurlvl-1);
//            result += IR * kr;
            result += IR * 0.5;
            if (polygon.isTrans) {
                Ray refr = ray.refractRay(intersect, normal);
                Color IT = light(refr, recurlvl-1);
//                result += IT * kt;
                result += IT * 0.5;
            }
        /* if intersecting object is a sphere */
//        } else {
//            Sphere sphere = spheres.at(id - polygons.size());
//            Vector normal = sphere.getNormal(intersect);
//            result = sphere.color * phong(intersect, ray.from, normal);
//            Ray refl = ray.reflectRay(intersect, normal);
//            Color IR = light(refl, recurlvl-1);
//            result += IR * kr;
//            if (sphere.isTrans) {
//                Ray refr = ray.refractRay(intersect, normal);
//                Color IT = light(refr, recurlvl-1);
//                result += IT * kt;
//            }
        /* if intersecting object is an ellipsoid */
        } else {
            Ellipsoid ellipsoid = ellipsoids.at(id - polygons.size());
            Vector normal = ellipsoid.getNormal(intersect);
            result = ellipsoid.color * phong(intersect, ray.from, normal);
            Ray refl = ray.reflectRay(intersect, normal);
            Color IR = light(refl, recurlvl-1);
//            result += IR * kr;
            result += IR * 0.5;
            if (ellipsoid.isTrans) {
                Ray refr = ray.refractRay(intersect, normal);
                Color IT = light(refr, recurlvl-1);
//                result += IT * kt;
                result += IT * 0.5;
            }
        }
    }
    return result;
}

void updateVector() {
    up.norm();
    hori = cross(at - from, up).norm();
}

void addObject() {
    vector<Point> vertices1 {Point(-1,2,2), Point(-1,2,1), Point(-1,1,2)};
//    vector<Point> vertices1 {Point(-1,0.5,0.5), Point(-1,0.5,-0.5), Point(-1,-0.5,0.5)};
    Polygon t1 = Polygon(vertices1, Vector(-1,0,0), Color(1,1,1), true);
    polygons.push_back(t1);
    
    vector<Point> vertices2 {Point(0,2,2), Point(-1,2,1), Point(-1,2,2)};
//    vector<Point> vertices2 {Point(0,0.5,0.5), Point(-1,0.5,-0.5), Point(-1,0.5,0.5)};
    Polygon t2 = Polygon(vertices2, Vector(0,1,0), Color(1,1,1), true);
    polygons.push_back(t2);
    
    vector<Point> vertices3 {Point(-1,2,2), Point(-1,1,2), Point(0,2,2)};
//    vector<Point> vertices3 {Point(-1,0.5,0.5), Point(-1,-0.5,0.5), Point(0,0.5,0.5)};
    Polygon t3 = Polygon(vertices3, Vector(0,0,1), Color(1,1,1), true);
    polygons.push_back(t3);
    
    vector<Point> vertices4 {Point(0,2,2), Point(-1,1,2), Point(-1,2,1)};
//    vector<Point> vertices4 {Point(0,0.5,0.5), Point(-1,-0.5,0.5), Point(-1,0.5,-0.5)};
    Polygon t4 = Polygon(vertices4, Vector(1,-1,-1), Color(1,1,1), true);
    polygons.push_back(t4);
    
//    Sphere s1 = Sphere(Point(0,0,0), 0.1, Color(1,1,1), true);
//    spheres.push_back(s1);
//    Sphere s2 = Sphere(Point(-0.5,1,-0.5), 0.4, Color(1,1,1), false);
//    spheres.push_back(s2);
//    Sphere s3 = Sphere(Point(1,1.5,1), 0.6, Color(1,1,1), true);
//    spheres.push_back(s3);
    
    Ellipsoid s1 = Ellipsoid(Point(0,0,0), 0.1, Color(1,1,1), true);
    ellipsoids.push_back(s1);
    Ellipsoid s2 = Ellipsoid(Point(-0.5,1,-0.5), 0.4, Color(1,1,1), false);
    ellipsoids.push_back(s2);
    Ellipsoid s3 = Ellipsoid(Point(1,1.5,1), 0.6, Color(1,1,1), true);
    ellipsoids.push_back(s3);
    
    updateVector();
}

void init();
void idle();
void display();
void draw_pix(int x, int y);
void key(unsigned char ch, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void check();

int main(int argc, char **argv)
{
    
    //the number of pixels in the grid
    grid_width = 500;
    grid_height = 500;
    
    //the size of pixels sets the inital window height and width
    //don't make the pixels too large or the screen size will be larger than
    //your display size
    pixel_size = 1;
    
    /*Window information*/
    win_height = grid_height*pixel_size;
    win_width = grid_width*pixel_size;
    
    /* add objects */
    addObject();
    
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    /*initialize variables, allocate memory, create buffers, etc. */
    //create window of size (win_width x win_height
    glutInitWindowSize(win_width,win_height);
    //windown title is "glut demo"
    glutCreateWindow("glut demo");
    
    glClearColor(0,0,0,0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,win_height,0,win_width);
    glClear(GL_COLOR_BUFFER_BIT);
    
    /*defined glut callback functions*/
    glutDisplayFunc(display); //rendering calls here
    glutMouseFunc(mouse);     //mouse button events
    glutMotionFunc(motion);   //mouse movement events
    glutKeyboardFunc(key);    //Keyboard events
    glutIdleFunc(idle);       //Function called while program is sitting "idle"
    
    //initialize opengl variables
    //init();
    //start glut event loop
    glutMainLoop();
    return 0;
}

/*initialize gl stufff*/
void init()
{
    //checks for OpenGL errors
    check();
}

//called repeatedly when glut isn't doing anything else
void idle()
{
    //redraw the scene over and over again
    glutPostRedisplay();
}

//this is where we render the screen
void display()
{
    //clears the screen
    glClear(GL_COLOR_BUFFER_BIT);
    
    /* update pixels */
    pixels.clear();
    for (int x = 0; x < win_width; x++) {
        for (int y = 0; y < win_height; y++) {
            Ray ray = cameraRay(x, y);
            Color color = light(ray, recursionLevel);
            pixels.push_back(color);
        }
    }
    /* normalizing */
    float brightest = 0.0;
    for (int i = 0; i < win_height*win_width; i++) {
        if (pixels.at(i).r > brightest) {
            brightest = pixels.at(i).r;
        }
        if (pixels.at(i).g > brightest) {
            brightest = pixels.at(i).g;
        }
        if (pixels.at(i).b > brightest) {
            brightest = pixels.at(i).b;
        }
    }
    if (brightest > 0.0) {
        for (int i = 0; i < win_height*win_width; i++) {
            pixels.at(i).r /= brightest;
            pixels.at(i).g /= brightest;
            pixels.at(i).b /= brightest;
        }
    }
    
    /* draw screen */
    for (int x = 0; x < win_width; x++) {
        for (int y = 0; y < win_height; y++) {
            Color cur = pixels.at(x * win_width + y);
            glBegin(GL_POINTS);
            glColor3f(cur.r,cur.g,cur.b);
            glVertex2f(x,y);
            glEnd();
        }
    }
    
    //blits the current opengl framebuffer on the screen
    glutSwapBuffers();
    //checks for opengl errors
    check();
}


//Draws a single "pixel" given the current grid size
//don't change anything in this for project 1
void draw_pix(int x, int y){
    glBegin(GL_POINTS);
    glColor3f(.2,.2,1.0);
    glVertex3f(x+.5,y+.5,0);
    glEnd();
}

//gets called when a key is pressed on the keyboard
void key(unsigned char ch, int x, int y)
{
    switch(ch)
    {
        default:
            //prints out which key the user hit
            printf("User hit the \"%c\" key\n",ch);
            break;
    }
    //redraw the scene after keyboard input
    glutPostRedisplay();
}


//gets called when a mouse button is pressed
void mouse(int button, int state, int x, int y)
{
    //print the pixel location, and the grid location
    printf ("MOUSE AT PIXEL: %d %d, GRID: %d %d\n",x,y,(int)(x/pixel_size),(int)((win_height-y)/pixel_size));
    switch(button)
    {
        case GLUT_LEFT_BUTTON: //left button
            printf("LEFT ");
            break;
        case GLUT_RIGHT_BUTTON: //right button
            printf("RIGHT ");
        default:
            printf("UNKNOWN "); //any other mouse button
            break;
    }
    if(state !=GLUT_DOWN)  //button released
        printf("BUTTON UP\n");
    else
        printf("BUTTON DOWN\n");  //button clicked
    
    //redraw the scene after mouse click
    glutPostRedisplay();
}

//gets called when the curser moves accross the scene
void motion(int x, int y)
{
    //redraw the scene after mouse movement
    glutPostRedisplay();
}

//checks for any opengl errors in the previous calls and
//outputs if they are present
void check()
{
    GLenum err = glGetError();
    if(err != GL_NO_ERROR)
    {
        printf("GLERROR: There was an error %s\n",gluErrorString(err) );
        exit(1);
    }
}
