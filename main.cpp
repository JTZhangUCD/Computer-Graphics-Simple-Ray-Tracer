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

using namespace std;

/*----------------------------- Global Variables -----------------------------*/
/*set in main()*/
//the number of pixels in the grid
int grid_width;
int grid_height;

//the size of pixels sets the inital window height and width
//don't make the pixels too large or the screen size will be larger than
//your display size
double pixel_size;

/*Window information*/
int win_height;
int win_width;

/* resolution */
int resol = 500;
/* refraction coefficient */
double refrCoef = 0.667;
/* maximum bumber of recursion levels */
int recursionLevel = 3;

/*---------------------------------- Vector ----------------------------------*/
struct Vector {
    double x,y,z;
    Vector(double x1, double y1, double z1) {
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
    inline double operator * (const Vector& v1)const {
        return x * v1.x + y * v1.y + z * v1.z;
    }
    inline Vector operator * (const double f)const {
        return Vector(x*f, y*f, z*f);
    }
//    inline Vector operator / (const double f)const {
//        return Vector(x/f, y/f, z/f);
//    }
//    inline Vector& operator += (const Vector& v1) {
//        x += v1.x; y += v1.y; z += v1.z;
//        return *this;
//    };
//    inline Vector& operator -= (const Vector& v1) {
//        x -= v1.x; y -= v1.y; z -= v1.z;
//        return *this;
//    };
//    inline Vector& operator *= (double f) {
//        x *= f; y *= f; z *= f;
//        return *this;
//    };
//    inline Vector& operator /= (double f) {
//        x /= f; y /= f; z /= f;
//        return *this;
//    };
    inline double length(void)const {
        return sqrt(x*x + y*y + z*z);
    }
    inline Vector& norm(void) {
        double length = sqrt(x*x + y*y + z*z);
        x /= length; y /= length; z /= length;
        return *this;
    }
};

/*----------------------------------- Point ----------------------------------*/
struct Point {
    double x,y,z;
    Point(double x1, double y1, double z1) {
        x = x1; y = y1; z = z1;
    }
    Point(void) {
        x = 0; y = 0; z = 0;
    }
    inline Vector operator - (const Point& p1)const {
        return Vector(x - p1.x, y - p1.y, z - p1.z);
    }
//    inline double operator * (const Point& p1)const {
//        return x * p1.x + y * p1.y + z * p1.z;
//    }
//    inline Point operator * (const double f)const {
//        return Point(x*f, y*f, z*f);
//    }
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

double dot(const Vector& v1, const Vector& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vector cross(const Vector& v1, const Vector& v2) {
    return Vector(v1.y*v2.z-v1.z*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x);
}

/*----------------------------------- Color ----------------------------------*/
struct Color {
    double r,g,b;
    Color(double r1, double g1, double b1) {
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
//    inline Color operator * (const Vector& v)const {
//        return Color(r * v.x, g * v.y, b * v.z);
//    }
    inline Color operator * (const Color& c)const {
        return Color(r * c.r, g * c.g, b * c.b);
    }
    inline Color operator * (const double f)const {
        return Color(r * f, g * f, b * f);
    }
//    inline Color operator / (const double f)const {
//        return Color(r / f, g / f, b / f);
//    }
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
//    inline Color &operator *= (double f) {
//        r *= f; g *= f; b *= f;
//        return *this;
//    };
//    inline Color &operator /= (double f) {
//        r /= f; g /= f; b /= f;
//        return *this;
//    };
};

/* phong variables */
Point x = Point(0,5,0);     // light source
double IA = 0.25;
double IL = 1;
double K = 2;
int specularity = 1;
Color ka = Color(1.0,1.0,1.0);
Color kd = Color(1.0,1.0,1.0);
Color ks = Color(1.0,1.0,1.0);

/* global illumination */
Color kr = Color(0.75,0.75,0.75);
Color kt = Color(0.75,0.75,0.75);

/* store info of and pixels' color*/
vector<Color> pixels;

/*------------------------------------ Ray -----------------------------------*/
/* camera ray */
Point from = Point(0.0,-2.5,0.0);
Point at   = Point(0.0,0.0,0.0);
Vector up   = Vector(1.0,0.0,0.0);
Vector hori = Vector(0.0,0.0,-1.0);
double fov = 0.7853; // in radian

struct Ray {
    Point p0;        // starting location
    Vector u;         // unit direction vector
    bool inObject;         // if line is in an object or not
    Ray(const Point& p01, const Vector& u1, bool inObject1) {
        p0 = p01;
        u = u1;
        inObject = inObject1;
    }
    Ray reflectRay(const Point& intersect, const Vector& N)const {
        Vector reflectVector = u - N * 2 * (u * N);
        return Ray(intersect, reflectVector.norm(), inObject);
    }
    Ray refractRay(const Point& intersect, const Vector& N)const {
        double coef = refrCoef;
        if (inObject == true) {
            coef = 1 / coef;
        }
        
        double cosi = -N * u;
        double cosr = pow(1 - pow(coef,2) * (1 - pow(cosi,2)),0.5);
        Vector refractVector = u * coef - N * (cosr - coef * cosi);
        return Ray(intersect, refractVector.norm(), !inObject);
    }
};

Ray cameraRay(int x, int y) {
    Vector xv = hori * tan(fov) * (2.0 * x / win_width - 1);
    Vector yv = up * tan(fov) * (2.0 * y / win_height - 1);
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
        //getNormal();
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
        
        double nx = (U.y * V.z) - (U.z * V.y);
        double ny = (U.z * V.x) - (U.x * V.z);
        double nz = (U.x * V.y) - (U.y * V.x);
        
        normal = Vector(nx,ny,nz);
    }
    
    bool intersect(const Ray& ray, Point& intersect) {
        double negD = normal * vertices.at(0).toVec();
        double s = (negD - (normal * ray.p0.toVec())) / (normal * ray.u);
        bool isInter = s >= 1e-4;
        Vector intersectVec = ray.p0.toVec() + ray.u * s;
        intersect = toPnt(intersectVec);
        
        for (int i = 0; i < vertices.size(); i++) {
            if (isInter) {
                Point p1, p2;
                if (i == vertices.size() - 1) {
                    p1 = vertices.at(i);
                    p2 = vertices.at(0);
                } else {
                    p1 = vertices.at(i);
                    p2 = vertices.at(i+1);
                }
                Vector side = cross(p1 - intersect, p2 - intersect);
                isInter = isInter && (side * normal > 0);
            }
        }
        return isInter;
    }
};

/* related globals */
vector<Polygon> polygons;

/*---------------------------------- Ellipsoid ----------------------------------*/
struct Ellipsoid {
    Point center;
    double radius;
    Color color;
    bool isTrans;
    Ellipsoid(const Point& center1, double radius1, const Color& color1, bool isTrans1) {
        center = center1;
        radius = radius1;
        color = color1;
        isTrans = isTrans1;
    }
    Vector getNormal(const Point& p) {
        return (p - center).norm();
    }
    bool intersect(const Ray& ray, Point& intersect) {
        Vector deltaP = center - ray.p0;
        double udeltaP = ray.u * deltaP;
        double sqrt1 = sqrt(pow(udeltaP,2) - pow(deltaP.length(),2) + pow(radius,2));
        if (sqrt1 < 0) return false;

        Vector intersect1 = ray.p0.toVec() + ray.u * (udeltaP - sqrt1);
        Vector intersect2 = ray.p0.toVec() + ray.u * (udeltaP + sqrt1);

        bool hasIntersect1 = udeltaP - sqrt1 >= 1e-4;
        bool hasIntersect2 = udeltaP + sqrt1 >= 1e-4;
        
        if (hasIntersect1 && hasIntersect2) {
            double length1 = (intersect1 - ray.p0.toVec()).length();
            double length2 = (intersect2 - ray.p0.toVec()).length();
            if (length1 >= length2) {
                intersect = toPnt(intersect2);
            } else {
                intersect = toPnt(intersect1);
            }
            return true;
        } else if (hasIntersect1) {
            intersect = toPnt(intersect1);
            return true;
        } else if (hasIntersect2) {
            intersect = toPnt(intersect2);
            return true;
        } else {
            return false;
        }
    }
};

/* related globals */
vector<Ellipsoid> ellipsoids;

/*------------------------------- Bounding Box -------------------------------*/
struct BoundingBox {
    vector<Polygon> faces;
    vector<Vector> normals;
    BoundingBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
        Point p1(xmin,ymin,zmin);
        Point p2(xmin,ymin,zmax);
        Point p3(xmin,ymax,zmin);
        Point p4(xmin,ymax,zmax);
        Point p5(xmax,ymin,zmin);
        Point p6(xmax,ymin,zmax);
        Point p7(xmax,ymax,zmin);
        Point p8(xmax,ymax,zmax);
        
        vector<Point> vertices1 {p1,p2,p4,p3};
        Polygon poly1 = Polygon(vertices1, Vector(1,0,0), Color(0,0,0), true);
        faces.push_back(poly1);
        vector<Point> vertices2 {p1,p5,p7,p3};
        Polygon poly2 = Polygon(vertices2, Vector(0,0,1), Color(0,0,0), true);
        faces.push_back(poly2);
        vector<Point> vertices3 {p1,p5,p6,p2};
        Polygon poly3 = Polygon(vertices3, Vector(0,1,0), Color(0,0,0), true);
        faces.push_back(poly3);
        vector<Point> vertices4 {p3,p4,p8,p7};
        Polygon poly4 = Polygon(vertices4, Vector(0,-1,0), Color(0,0,0), true);
        faces.push_back(poly4);
        vector<Point> vertices5 {p2,p4,p8,p6};
        Polygon poly5 = Polygon(vertices5, Vector(0,0,-1), Color(0,0,0), true);
        faces.push_back(poly5);
        vector<Point> vertices6 {p5,p6,p8,p7};
        Polygon poly6 = Polygon(vertices6, Vector(-1,0,0), Color(0,0,0), true);
        faces.push_back(poly6);
        
    }
    bool intersect(const Ray& ray, Point& intersect) {
        bool hitBox = false;
        for (int index = 0; index < faces.size(); index++) {
            vector<Point> vertices = faces.at(index).vertices;
            Vector normal = normals.at(index);
            double negD = normal * vertices.at(0).toVec();
            double s = (negD - (normal * ray.p0.toVec())) / (normal * ray.u);
            bool isInter = s >= 1e-4;
            Vector intersectVec = ray.p0.toVec() + ray.u * s;
            intersect = toPnt(intersectVec);
            for (int i = 0; i < vertices.size(); i++) {
                if (isInter) {
                    Point p1, p2;
                    if (i == vertices.size() - 1) {
                        p1 = vertices.at(i);
                        p2 = vertices.at(0);
                    } else {
                        p1 = vertices.at(i);
                        p2 = vertices.at(i+1);
                    }
                    Vector side = cross(p1 - intersect, p2 - intersect);
                    isInter = isInter && (side * normal > 0);
                }
            }
            hitBox = hitBox || isInter;
        }
        return hitBox;
    }
};

BoundingBox boundingBox(-10,10,-10,10,-10,10);

bool intersectObject(const Ray& ray, Point& intersectPoint, int& id) {
    Point resultPoint;
    int resultID = -1;
    double resultDist = 0;
    
    /* check if intersect with a polygon */
    for (int i = 0; i < (int)polygons.size(); i++) {
        Point curPoint;
        if (polygons.at(i).intersect(ray, curPoint)) {
            double curDist = (curPoint - ray.p0).length();
            if (resultID == -1) {
                resultDist = curDist;
                resultID = i;
                resultPoint = curPoint;
            } else if (resultID != -1 && resultDist > curDist) {
                resultDist = curDist;
                resultID = i;
                resultPoint = curPoint;
            }
        }
    }
    
    /* check if intersect with a ellipsoid */
    for (int i = 0; i < (int)ellipsoids.size(); i++) {
        Point curPoint;
        if (ellipsoids.at(i).intersect(ray, curPoint)) {
            double curDist = (curPoint - ray.p0).length();
            if (resultID == -1) {
                resultDist = curDist;
                resultID = i + (int)polygons.size();
                resultPoint = curPoint;
            } else if (resultID != -1 && resultDist > curDist) {
                resultDist = curDist;
                resultID = i + (int)polygons.size();
                resultPoint = curPoint;
            }
        }
    }

    if (resultID != -1) {
        id = resultID;
        intersectPoint = resultPoint;
        return true;
    }
    return false;
}

Color phong(const Point& p, const Point& f, const Vector& n) {
    Vector l = (x - p).norm();
    Vector v = (f - p).norm();
    Vector r = (n * 2.0 * (n * l) - l).norm();
    Ray ray = Ray(p, l, false);
    double phongCoef = IL / ((f - p).length() + K);
    
    Color ambient = ka * IA;
    Color diffuse = (l * n > 0) ? kd * (l * n) : Color();
    Color specular = (l * n > 0 && r * v > 0) ? ks * pow((r * v),specularity) : Color();
    
    Point p1; int id;
    if (intersectObject(ray, p1, id)) {
        return ambient;
    } else {
        return ambient + (diffuse + specular) * phongCoef;
    }
}

Color processRay(const Ray& ray, int recurlvl) {
    Color result;
    
    // base case
    if (recurlvl == 0) {
        return result;
    }
    
    // general case
    Point intersect; int objectID;
    if (intersectObject(ray, intersect, objectID)) {
        /* if intersecting object is a polygon */
        if (objectID < polygons.size()) {
            Polygon polygon = polygons.at(objectID);
            Vector n = polygon.normal;
            result = polygon.color * phong(intersect, ray.p0, n);
            Ray refl = ray.reflectRay(intersect, n);
            Color IR = processRay(refl, recurlvl-1);
            result += IR * kr;
            if (polygon.isTrans) {
                Ray refr = ray.refractRay(intersect, n);
                Color IT = processRay(refr, recurlvl-1);
                result += IT * kt;
            }
        /* if intersecting object is a ellipsoid */
        } else {
            Ellipsoid ellipsoid = ellipsoids.at(objectID - polygons.size());
            Vector n = ellipsoid.getNormal(intersect);
            result = ellipsoid.color * phong(intersect, ray.p0, n);
            Ray refl = ray.reflectRay(intersect, n);
            Color IR = processRay(refl, recurlvl-1);
            result += IR * kr;
            if (ellipsoid.isTrans) {
                Ray refr = ray.refractRay(intersect, n);
                Color IT = processRay(refr, recurlvl-1);
                result += IT * kt;
            }
        }
        return result;
    } else { // hits bounding box if the ray does not intersect any object.
        // return black color
        return result;
    }
}

void initObjects() {
    /* create a triangle */
    vector<Point> vertices1 {Point(-1,2,1.7), Point(-1,2,0.7), Point(-1,1,1.7)};
    Polygon t1 = Polygon(vertices1, Vector(-1,0,0), Color(1,1,0.5), true);
    polygons.push_back(t1);
    vector<Point> vertices2 {Point(0,2,1.7), Point(-1,2,0.7), Point(-1,2,1.7)};
    Polygon t2 = Polygon(vertices2, Vector(0,1,0), Color(1,1,0.5), true);
    polygons.push_back(t2);
    vector<Point> vertices3 {Point(-1,2,1.7), Point(-1,1,1.7), Point(0,2,1.7)};
    Polygon t3 = Polygon(vertices3, Vector(0,0,1), Color(1,1,0.5), true);
    polygons.push_back(t3);
    vector<Point> vertices4 {Point(0,2,1.7), Point(-1,1,1.7), Point(-1,2,0.7)};
    Polygon t4 = Polygon(vertices4, Vector(1,-1,-1), Color(1,1,0.5), true);
    polygons.push_back(t4);

    /* create a ellipsoid */
    Ellipsoid e1 = Ellipsoid(Point(-0.5,1,-0.5), 0.6, Color(1,0.5,1), false);
    ellipsoids.push_back(e1);
    
    /* create a ellipsoid */
    Ellipsoid e2 = Ellipsoid(Point(1,0.5,0.5), 0.4, Color(0.5,1,1), true);
    ellipsoids.push_back(e2);
    
    /* create a ellipsoid */
    Ellipsoid e3 = Ellipsoid(Point(1,0.5,-0.5), 0.3, Color(0.5,1,0.5), true);
    ellipsoids.push_back(e3);
    
    /* create a plane */
    vector<Point> vertices5 {Point(-2,-6,-4), Point(-2,6,-4), Point(-2,6,4), Point(-2,-6,4)};
    Polygon p1 = Polygon(vertices5, Vector(1,0,0), Color(1,1,1), false);
    polygons.push_back(p1);
}

void introduction(void) {
    cout << "******************** Welcome to Project 5! ********************" << endl;
    cout << "******** Note: press i on screen to enable interface! *********" << endl;
}

int userInput = -1;
void interface(void) {
    cout << "> To change camera variables    enter 1" << endl;
    cout << "> To change light position      enter 2" << endl;
    cout << "> To change resolution          enter 3" << endl;
    cout << "> To change phong variables     enter 4" << endl;
    cout << "> To change refractive coef     enter 5" << endl;
    cout << "> To change max recursion level enter 6" << endl;
    cout << "> To print current variables    enter 7" << endl;
    cout << "> To exit                       enter 0" << endl;
    
    cin >> userInput;
    while (!(0 <= userInput && userInput <= 7)) {
        cout << "> invalid input(" << userInput << "), try again." << endl;
        cin >> userInput;
    }
    
    if (userInput == 0) {
        cout << "exit success!" << endl;
        exit(0);
    }
    
    if (userInput == 1) {
        cout << "> Current from Vector: (";
        cout << from.x << " " << from.y << " " << from.z << ")" << endl;
        cout << "> Current at Vector:   (";
        cout << at.x << " " << at.y << " " << at.z << ")" << endl;
        cout << "> Current up Vector:   (";
        cout << up.x << " " << up.y << " " << up.z << ")" << endl;
        cout << "> Current viewing angle (in radian):";
        cout << fov << endl;
        
        Point tempFrom;
        Point tempAt;
        Vector tempUp;
        double tempfov;
        cout << "> Enter the new from, at, and up vectors in following format:" << endl;
        cout << "> from.x from.y from.z at.x at.y at.z up.x up.y up.z angle" << endl;
        cin >> tempFrom.x >> tempFrom.y >> tempFrom.z;
        cin >> tempAt.x >> tempAt.y >> tempAt.z;
        cin >> tempUp.x >> tempUp.y >> tempUp.z;
        cin >> tempfov;
        
        /* check if from - at && up are perpendicular */
        if (dot(tempFrom - tempAt, tempUp) != 0) {
            cout << "> invalid input: (from-at) is not perpendicular to up" << endl;
            return;
        }
        
        from = tempFrom;
        at = tempAt;
        up = tempUp;
        fov = tempfov;
        up.norm();
        hori = cross(at - from, up).norm();
    }
    
    if (userInput == 2) {
        cout << "> Current light position: (";
        cout << x.x << " " << x.y << " " << x.z << ")" << endl;
        
        cout << "> Enter the new light position in following format:" << endl;
        cout << "> light.x light.y light.z" << endl;
        cin >> x.x >> x.y >> x.z;
    }
    
    if (userInput == 3) {
        cout << "> Current resolution: " << resol << endl;
        cout << "> Enter the new resolution (in integer):" << endl;
        cin >> resol;
    }
    
    if (userInput == 4) {
        cout << "> Current ka RGB: (";
        cout << ka.r << " " << ka.g << " " << ka.b << ")" << endl;
        cout << "> Current kd RGB: (";
        cout << kd.r << " " << kd.g << " " << kd.b << ")" << endl;
        cout << "> Current ks RGB: (";
        cout << ks.r << " " << ks.g << " " << ks.b << ")" << endl;
        cout << "> Current kr RGB: (";
        cout << kr.r << " " << kr.g << " " << kr.b << ")" << endl;
        cout << "> Current kt RGB: (";
        cout << kt.r << " " << kt.g << " " << kt.b << ")" << endl;
        cout << "> Current IA: " << IA << endl;
        cout << "> Current IL: " << IL << endl;
        cout << "> Current K: " << K << endl;
        cout << "> Current Specularity: " << specularity << endl;
        cout << "> Current light position: (";
        cout << x.x << " " << x.y << " " << x.z << ")" << endl;
        
        cout << "> To change ka/kd/ks/kr/kt,      enter 1" << endl;
        cout << "> To change IA/IL/K/Specularity, enter 2" << endl;
        cout << "> To change light Position,      enter 3" << endl;
        
        cout << "> Enter ka/kd/ks/kr/kt followed by RGB values." << endl;
        cout << "> Example: ka 1 1 1" << endl;
        cout << "> or" << endl;
        cout << "> Enter IA/IL/K/Specularity followed by its value." << endl;
        cout << "> Example: IA 0.5" << endl;
        cout << "> Example: Specularity 1" << endl;
        cout << "> or" << endl;
        cout << "> Enter light followed by its xyz values." << endl;
        cout << "> Example: light 0 0 0" << endl;
        
        string name = "";
        cin >> name;
        if (name == "ka") {
            cin >> ka.r >> ka.g >> ka.b;
        } else if (name == "kd") {
            cin >> kd.r >> kd.g >> kd.b;
        } else if (name == "ks") {
            cin >> ks.r >> ks.g >> ks.b;
        } else if (name == "kr") {
            cin >> kr.r >> kr.g >> kr.b;
        } else if (name == "kt") {
            cin >> kt.r >> kt.g >> kt.b;
        } else if (name == "IA") {
            cin >> IA;
        } else if (name == "IL") {
            cin >> IL;
        } else if (name == "K") {
            cin >> K;
        } else if (name == "Specularity") {
            cin >> specularity;
        } else if (name == "light") {
            cin >> x.x >> x.y >> x.z;
        } else {
            cout << "invalid argument" << endl;
            return;
        }
    }
    
    if (userInput == 5) {
        cout << "> Current refractive coef (from air to object): " << refrCoef << endl;
        cout << "> Enter the new refractive coef:" << endl;
        cin >> refrCoef;
    }
    
    if (userInput == 6) {
        cout << "> Current maximum recursion level: " << recursionLevel << endl;
        cout << "> Enter the new maximum recursion level (as integer):" << endl;
        cin >> recursionLevel;
    }
    
    if (userInput == 7) {
        cout << "> from Vector: (";
        cout << from.x << " " << from.y << " " << from.z << ")" << endl;
        cout << "> at Vector:   (";
        cout << at.x << " " << at.y << " " << at.z << ")" << endl;
        cout << "> up Vector:   (";
        cout << up.x << " " << up.y << " " << up.z << ")" << endl;
        cout << "> viewing angle (in radian):";
        cout << fov << endl;
        
        cout << "> light position: (";
        cout << x.x << " " << x.y << " " << x.z << ")" << endl;
        
        cout << "> resolution: " << resol << endl;
        
        cout << "> ka RGB: (";
        cout << ka.r << " " << ka.g << " " << ka.b << ")" << endl;
        cout << "> kd RGB: (";
        cout << kd.r << " " << kd.g << " " << kd.b << ")" << endl;
        cout << "> ks RGB: (";
        cout << ks.r << " " << ks.g << " " << ks.b << ")" << endl;
        cout << "> kr RGB: (";
        cout << kr.r << " " << kr.g << " " << kr.b << ")" << endl;
        cout << "> kt RGB: (";
        cout << kt.r << " " << kt.g << " " << kt.b << ")" << endl;
        cout << "> IA: " << IA << endl;
        cout << "> IL: " << IL << endl;
        cout << "> K: " << K << endl;
        cout << "> Specularity: " << specularity << endl;
    }
    
    return;
}

void init();
void idle();
void display();
void draw_pix(int x, int y);
void reshape(int width, int height);
void key(unsigned char ch, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void check();

int main(int argc, char **argv)
{
    //the number of pixels in the grid
    grid_width  = resol;
    grid_height = resol;
    
    //the size of pixels sets the inital window height and width
    //don't make the pixels too large or the screen size will be larger than
    //your display size
    pixel_size = 1;
    
    /*Window information*/
    win_height = grid_height*pixel_size;
    win_width = grid_width*pixel_size;
    
    /* add objects */
    initObjects();
    
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    /*initialize variables, allocate memory, create buffers, etc. */
    //create window of size (win_width x win_height
    glutInitWindowSize(win_width,win_height);
    //windown title is "glut demo"
    glutCreateWindow("Project5");
    
    /* print introduction */
    introduction();
    
    /*defined glut callback functions*/
    glutDisplayFunc(display); //rendering calls here
    glutReshapeFunc(reshape); //update GL on window size change
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
            Color color = processRay(ray, recursionLevel);
            pixels.push_back(color);
        }
    }
    /* normalizing */
    double brightest = 0.0;
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

/*Gets called when display size changes, including initial craetion of the display*/
void reshape(int width, int height)
{
    /*set up projection matrix to define the view port*/
    //update the ne window width and height
    win_width = width;
    win_height = height;
    
    //creates a rendering area across the window
    glViewport(0,0,width,height);
    // up an orthogonal projection matrix so that
    // the pixel space is mapped to the grid space
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,grid_width,0,grid_height,-10,10);
    
    //clear the modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    //set pixel size based on width, if the aspect ratio
    //changes this hack won't work as well
    pixel_size = width/(double)grid_width;
    
    //set pixel size relative to the grid cell size
    glPointSize(pixel_size);
    //check for opengl errors
    check();
}

//gets called when a key is pressed on the keyboard
void key(unsigned char ch, int x, int y)
{
    switch(ch)
    {
        case 'i': {
            interface();
        }
        default: {
            //prints out which key the user hit
            //printf("User hit the \"%c\" key\n",ch);
            break;
        }
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

