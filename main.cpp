#include <GLUT/glut.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>

#define WINDOW_SIZE 600
#define PI 3.141593
#define ZERO(x) (-0.00001 <= x && x < 0.00001)
#define NEG(x) (x <= -0.00001)
#define POS(x) (0.00001 <= x)

using namespace std;

//============================================================================== POINT
struct Point {
    float x,y,z;
    Point(float x1, float y1, float z1) {
        x = x1; y = y1; z = z1;
    }
    Point(void) {
        x = 0; y = 0; z = 0;
    }
    inline Point operator + (const Point& v1)const {
        return Point(x+v1.x, y+v1.y, z+v1.z);
    }
    inline Point operator - (void)const {
        return Point(-x, -y, -z);
    }
    inline Point operator - (const Point& v1)const {
        return Point(x-v1.x, y-v1.y, z-v1.z);
    }
    inline double operator * (const Point& v1)const{
        return x*v1.x + y*v1.y + z*v1.z;
    }
    inline Point operator * (const float f)const{
        return Point(x*f, y*f, z*f);
    }
    inline Point operator / (const float f)const{
        return Point(x/f, y/f, z/f);
    }
    inline Point &operator += (const Point& v1) {
        x += v1.x; y += v1.y; z += v1.z;
        return *this;
    }
    inline Point &operator -= (const Point& v1) {
        x -= v1.x; y -= v1.y; z -= v1.z;
        return *this;
    }
    inline Point &operator *= (float f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    inline Point &operator /= (float f) {
        x /= f; y /= f; z /= f;
        return *this;
    }
    inline float length(void)const {
        return sqrt(x*x + y*y + z*z);
    }
    inline Point& norm(void) {
        float length = sqrt(x*x + y*y + z*z);
        x /= length; y /= length; z /= length;
        return *this;
    }
};

Point cross(const Point& v1, const Point& v2) {
    return Point(v1.y*v2.z - v1.z*v2.y,v1.z*v2.x - v1.x*v2.z,v1.x*v2.y - v1.y*v2.x);
}

float dot(const Point& v1, const Point& v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

istream& operator>>(istream& in, Point& p) {
    in >> p.x >> p.y >> p.z;
    return in;
}

ostream& operator<<(ostream& out, const Point& p) {
    out << p.x << ' ' << p.y << ' ' << p.z;
    return out;
}
    
//============================================================================== COLOR

class Color {
public:
    float r,g,b;
    Color(float r1, float g1, float b1) {
        r = r1; g = g1; b = b1;
    }
    Color(void) {
        r = 0; g = 0; b = 0;
    }
    inline Color operator + (const Color& v1)const {
        return Color(r+v1.r, g+v1.g, b+v1.b);
    }
    inline Color operator - (const Color& v1)const {
        return Color(r-v1.r, g-v1.g, b-v1.b);
    }
    inline Color operator * (const Point& v1)const {
        return Color(r-v1.x, g-v1.y, b-v1.z);
    }
    inline Color operator * (const float f)const {
        return Color(r*f, g*f, b*f);
    }
    inline Color operator / (const float f)const {
        return Color(r/f, g/f, b/f);
    }
    inline Color &operator += (const Color& v1) {
        r += v1.r; g += v1.g; b += v1.b;
        return *this;
    }
    inline Color &operator -= (const Color& v1) {
        r -= v1.r; g -= v1.g; b -= v1.b;
        return *this;
    }
    inline Color &operator *= (const Point& v1) {
        r *= v1.x; g *= v1.y; b *= v1.z;
        return *this;
    }
    inline Color &operator *= (float f) {
        r *= f; g *= f; b *= f;
        return *this;
    }
    inline Color &operator /= (float f) {
        r /= f; g /= f; b /= f;
        return *this;
    }
};

istream& operator>>(istream& in, Color& p) {
    in >> p.r >> p.g >> p.b;
    return in;
}
ostream& operator<<(ostream& out, const Color& p) {
    out << p.r << ' ' << p.g << ' ' << p.b;
    return out;
}

//============================================================================== OBJECT
float iA = 0.2, iL = 1;
Point kA(1, 1, 1), kD(1, 1, 1), kS(1, 1, 1);
int bigN = 5;
float bigK = 10;
float kR = 0.5, kT = 0.5;
float eta = 1.5;
Point lightSource(0, 1, 0);

class Line {
public:
    Point start;
    Point u;
    bool inside;
public:
    Line(const Point& start_, const Point& u_, bool inside_):
        start(start_), u(u_), inside(inside_){}
    Line():inside(false){}
};

class Primitive {
public:
    bool transparent;
    Color color;
public:
    virtual bool intersect(const Line& l, Point& result) = 0;
    virtual Point normal(const Point& p) = 0;
};

class Polygon : public Primitive {
public:
    vector<Point> v;
    Point n;
public:
    bool intersect(const Line& l, Point& result);
    Point normal(const Point& p);
};

class Sphere : public Primitive {
public:
    Point center;
    float radius;
public:
    bool intersect(const Line& l, Point& result);
    Point normal(const Point& p);
};

class Quadratic : public Primitive {
public:
    float A, B, C, D, E, F, G, H, I, J;
public:
    bool intersect(const Line& l, Point& result);
    Point normal(const Point& p);
};

inline void makePixel(int x, int y, const Color& c, GLfloat *pixels, int resolution) {
    if (0 <= x && x < resolution && 0 <= y && y < resolution) {
        pixels[(y * resolution + x) * 3] = c.r;
        pixels[(y * resolution + x) * 3 + 1] = c.g;
        pixels[(y * resolution + x) * 3 + 2] = c.b;
    }
}
inline void normalize(GLfloat *pixels, int resolution) {
    GLfloat best = 0.0;
    for (int i = 0; i < resolution * resolution * 3; i++) {
        if (pixels[i] > best)
            best = pixels[i];
    }
    if (best > 0.0) {
        for (int i = 0; i < resolution * resolution * 3; i++)
            pixels[i] /= best;
    }
}

bool Polygon::intersect(const Line& l, Point& result) {
    float D = -dot(n, v[0]);
    float s = -(D+dot(n, l.start)) / dot(n, l.u);
    bool ret = POS(s);
    result = l.start + l.u * s;
    for (int i = 0; i < v.size(); i++)
        if (ret)
            ret = ret && dot(cross(v[i % v.size()] - result, v[(i+1) % v.size()] - result), n) > 0;
    return ret;
}

Point Polygon::normal(const Point& p) {
    return n;
}

bool Sphere::intersect(const Line& l, Point& result) {
    Point diff = center - l.start;
    float base = dot(l.u, diff);
    float delta = base * base - dot(diff, diff) + radius * radius;
    if (delta < 0)
        return false;
    delta = sqrt(delta);
    Point p1 = l.start + l.u * (base - delta);
    Point p2 = l.start + l.u * (base + delta);
    bool v1 = POS((p1 - l.start).length()) && POS(base - delta);
    bool v2 = POS((p2 - l.start).length()) && POS(base + delta);
    if (v1 && v2)
        result = (p1 - l.start).length() < (p2 - l.start).length() ? p1 : p2;
    else if (v1)
        result = p1;
    else if (v2)
        result = p2;
    else return false;
    return true;
}

Point Sphere::normal(const Point& p) {
    Point ret = (p - center).norm();
    return ret;
}

bool Quadratic::intersect(const Line& l, Point& result) {
    float Aq = A * l.u.x * l.u.x + B * l.u.y * l.u.y + C * l.u.z * l.u.z
                + D * l.u.x * l.u.y + E * l.u.x * l.u.z + F * l.u.y * l.u.z;
    float Bq = 2*A * l.start.x * l.u.x + 2*B * l.start.y * l.u.y + 2*C * l.start.z * l.u.z
                + D * (l.start.x * l.u.y + l.start.y * l.u.x)
                + E * (l.start.x * l.u.z + l.start.z * l.u.x)
                + F * (l.start.y * l.u.z + l.start.z * l.u.y)
                + G * l.u.x + H * l.u.y + I * l.u.z;
    float Cq = A * l.start.x * l.start.x + B * l.start.y * l.start.y + C * l.start.z * l.start.z
                + D * l.start.x * l.start.y + E * l.start.x * l.start.z + F * l.start.y * l.start.z
                + G * l.start.x + H * l.start.y + I * l.start.z + J;
    float delta = Bq * Bq - 4 * Aq * Cq;
    if (delta < 0)
        return false;
    delta = sqrt(delta);
    Point p1, p2;
    if (Aq == 0)
        p1 = p2 = l.start + l.u * (-Cq / Bq);
    else {
        p1 = l.start + l.u * (- Bq - delta) / 2 / Aq;
        p2 = l.start + l.u * (- Bq + delta) / 2 / Aq;
    }
    bool v1 = POS((p1 - l.start).length()) && POS(- Bq - delta);
    bool v2 = POS((p2 - l.start).length()) && POS(- Bq + delta);
    if (v1 && v2)
        result = (p1 - l.start).length() < (p2 - l.start).length() ? p1 : p2;
    else if (v1)
        result = p1;
    else if (v2)
        result = p2;
    else return false;
    return true;
}

Point Quadratic::normal(const Point& p) {
    Point ret;
    ret.x = 2*A*p.x + D*p.y + E*p.z + G;
    ret.y = 2*B*p.y + D*p.x + F*p.z + H;
    ret.z = 2*C*p.z + E*p.x + F*p.y + I;
    ret.norm();
    return ret;
}

Polygon* newTriangle(const Point& v1, const Point& v2, const Point& v3, const Point& normal) {
    Polygon* ret = new Polygon();
    ret->v.push_back(v1);
    ret->v.push_back(v2);
    ret->v.push_back(v3);
    ret->n = normal;
    return ret;
}

Polygon* newSquare(const Point& v1, const Point& v2, const Point& v3, const Point& v4, const Point& normal) {
    Polygon* ret = new Polygon();
    ret->v.push_back(v1);
    ret->v.push_back(v2);
    ret->v.push_back(v3);
    ret->v.push_back(v4);
    ret->n = normal;
    return ret;
}

Quadratic* newSphere(const Point& center, float radius) {
    Quadratic* ret = new Quadratic();
    ret->A = ret->B = ret->C = 1;
    ret->D = ret->E = ret->F = 0;
    ret->G = - 2 * center.x;
    ret->H = - 2 * center.y;
    ret->I = - 2 * center.z;
    ret->J = dot(center, center) - radius * radius;
    return ret;
}

Quadratic* newEllipsoid(const Point& center, float a, float b, float c) {
    Quadratic* ret = new Quadratic();
    ret->A = b * b + c * c;
    ret->B = a * a + c * c;
    ret->C = a * a + b * b;
    ret->D = ret->E = ret->F = 0;
    ret->G = b * b * c * c * -2 * center.x;
    ret->H = a * a * c * c * -2 * center.y;
    ret->I = a * a * b * b * -2 * center.z;
    ret->J = b * b * c * c * center.x * center.x
        +   a * a * c * c * center.y * center.y
        +   a * a * b * b * center.z * center.z
        -   a * a * b * b * c * c;
    return ret;
}

Point refract(const Line& l, const Point& normal) {
    float r = l.inside ? eta : 1 / eta;
    float cosThetaI = dot(-l.u, normal);
    float cosThetaR = sqrt(1 - (r * r) * (1 - cosThetaI * cosThetaI));
    return l.u * r - normal * (cosThetaR - r * cosThetaI);
}

Point reflect(const Point& u, const Point& normal) {
    return u - normal * 2 * dot(normal, u);
}
    
bool intersect(const Line& line, const vector<Primitive*>& primitives, Point& result, int& resultId) {
    Point best, now;
    int bestId;
    float bestDist = 0, nowDist;
    bool ret = false;
    for (int i = 0; i < primitives.size(); i++) {
        if (primitives[i]->intersect(line, now)) {
            nowDist = (now - line.start).length();
            if (!ret || (nowDist < bestDist)) {
                ret = true;
                best = now;
                bestId = i;
                bestDist = nowDist;
            }
        }
    }
    if (ret) {
        result = best;
        resultId = bestId;
    }
    return ret;
}

Point phong(const Point& point, const Point& ref, const Point& normal, const vector<Primitive*>& primitives) {
    Point intensity = kA * iA;
    Point l = (lightSource - point).norm();
    Point v = (ref - point).norm();
    Point r = (normal * 2.0 * dot(normal, l) - l).norm();
    Point diff;
    Line line; line.start = point; line.u = (lightSource - point).norm();
    Point result;
    int resultId;
    if (intersect(line, primitives, result, resultId))
        return intensity;
    float productDiff = dot(l, normal);
    if (productDiff > 0) diff = kD * productDiff;
    Point spec;
    float productSpec = dot(r, v);
    if (productDiff > 0 && productSpec > 0) spec = kS * productSpec;
    intensity += (diff + spec) * iL / ((ref - point).length() + bigK);
    return intensity;
}

Color light(const Line& line, const vector<Primitive*>& primitives, int depth) {
    Color ret;
    Primitive* prim;
    int id;
    Point point;
    Point normal;
    if (depth == 0) return ret;
    if (intersect(line, primitives, point, id)) {
        prim = primitives[id];
        normal = prim->normal(point);
        Point ph = phong(point, line.start, normal, primitives);
        ret = prim->color * ph;
        Line refl;
        refl.start = point;
        refl.u = reflect(line.u, normal).norm();
        refl.inside = line.inside;
        Color rf = light(refl, primitives, depth - 1);
        ret += rf * kR;
        if (prim->transparent) {
            Line refr;
            refr.start = point;
            refr.u = refract(line, normal).norm();
            refr.inside = !line.inside;
            Color rt = light(refr, primitives, depth - 1);
            ret += rt * kT;
        }
    }
    return ret;
}

//============================================================================== MAIN

int window;

float canvas[WINDOW_SIZE * WINDOW_SIZE * 3];

vector<Primitive*> objects;

Point from(0, -3, 0);
Point at(0, 0, 0);
Point up(0, 0, 1);
Point rv;
float angDeg = 45;
int resolution = WINDOW_SIZE;
int depth = 5;

Line ray(int x, int y) {
    Line ret;
    ret.start = from;
    float xx = (float)x/(float)resolution*2 - 1;
    float yy = (float)y/(float)resolution*2 - 1;
    Point p = from;
    p += rv * tan(angDeg / 180 * PI) * xx;
    p += up * tan(angDeg / 180 * PI) * yy;
    p += (at - from).norm();
    ret.u = (p - from).norm();
    return ret;
}

void updatePoint() {
    up.norm();
    rv = cross(at - from, up).norm();
}

void initWindow() {
    glClearColor(0, 0, 0, 0);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    memset(canvas, 0.0f, sizeof(canvas));
    for (int x = 0; x < resolution; x++)
        for (int y = 0; y < resolution; y++) {
            Line r = ray(x, y);
            Color f = light(r, objects, depth);
            makePixel(x, y, f, canvas, resolution);
        }
    normalize(canvas, resolution);
    glDrawPixels(resolution, resolution, GL_RGB, GL_FLOAT, canvas);
    glutSwapBuffers();
}

void refreshFunc() {
    glutPostWindowRedisplay(window);
}

void get() {
    cout << "from=" << from << endl;
    cout << "at=" << at << endl;
    cout << "up=" << up << endl;
    cout << "angle=" << angDeg << endl;
    cout << "I_a=" << iA << " I_l=" << iL << endl;
    cout << "k_a=" << kA << " k_d=" << kD << " k_s=" << kS << endl;
    cout << "K=" << bigK << " n=" << bigN << endl;
    cout << "k_r=" << kR << " k_t=" << kT << endl;
    cout << "eta=" << eta << endl;
    cout << "depth=" << depth << endl;
    cout << "resolution=" << resolution << endl;
}

void set(const string& command) {
    if (!command.size()) return;
    int pos = command.find('=');
    string name = command.substr(0, pos);
    stringstream value(command.substr(pos + 1, command.size() - pos - 1));
    cout << name << " " << value.str() << endl;
    if (name == "from") value >> from;
    if (name == "at") value >> at;
    if (name == "up") value >> up;
    if (name == "angle") value >> angDeg;
    if (name == "I_a") value >> iA;
    if (name == "I_l") value >> iL;
    if (name == "k_a") value >> kA;
    if (name == "k_d") value >> kD;
    if (name == "k_s") value >> kS;
    if (name == "K") value >> bigK;
    if (name == "n") value >> bigN;
    if (name == "k_r") value >> kR;
    if (name == "k_t") value >> kT;
    if (name == "eta") value >> eta;
    if (name == "depth") value >> depth;
    if (name == "resolution") value >> resolution;
}

void set() {
    string command;
    cout << ">";
    getline(cin, command);
    set(command);
    updatePoint();
    cout << "OK!" << endl;
    refreshFunc();
}

void init(int argc, char** argv) {
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);

    Quadratic* org = newSphere(Point(0, 0, 0), 0.1);
    org->color = Color(1, 1, 1);
    org->transparent = true;
    Quadratic* s = newSphere(Point(-0.5, 1, -0.5), 0.4);
    s->color = Color(1, 1, 1);
    s->transparent = false;
    Quadratic* s2 = newSphere(Point(1, 1.5, 1), 0.6);
    s2->color = Color(1, 1, 1);
    s2->transparent = true;
    Polygon* t = newTriangle(Point(-1, 2, 2),
                             Point(-1, 2, 1),
                             Point(-1, 1, 2),
                             Point(-1, 0, 0));
    t->color = Color(1, 1, 1);
    t->transparent = true;
    Polygon* t2 = newTriangle(Point(0, 2, 2),
                              Point(-1, 2, 1),
                              Point(-1, 2, 2),
                              Point(0, 1, 0));
    t2->color = Color(1, 1, 1);
    t2->transparent = true;
    Polygon* t3 = newTriangle(Point(-1, 2, 2),
                              Point(-1, 1, 2),
                              Point(0, 2, 2),
                              Point(0, 0, 1));
    t3->color = Color(1, 1, 1);
    t3->transparent = true;
    Polygon* t4 = newTriangle(Point(0, 2, 2),
                              Point(-1, 1, 2),
                              Point(-1, 2, 1),
                              Point(1, -1, -1));
    t4->color = Color(1, 1, 1);
    t4->transparent = true;
    Polygon* p = newSquare(Point(5, 6, -4), Point(5, 6, 4),
                           Point(-5, 6, 4), Point(-5, 6, -4),
                           Point(0, -1, 0));
    p->color = Color(1, 0.85, 0.25);
    p->transparent = false;
    Polygon* p2 = newSquare(Point(-5, -6, -4), Point(-5, 6, -4),
                           Point(-5, 6, 4), Point(-5, -6, 4),
                           Point(1, 0, 0));
    p2->color = Color(0, 0.5, 1);
    p2->transparent = false;
    Polygon* p3 = newSquare(Point(5, 6, -4), Point(5, -6, -4),
                           Point(5, -6, 4), Point(5, 6, 4),
                           Point(-1, 0, 0));
    p3->color = Color(1, 0.2, 0);
    p3->transparent = false;
    Polygon* p4 = newSquare(Point(5, 6, 4), Point(5, -6, 4),
                           Point(-5, -6, 4), Point(-5, 6, 4),
                           Point(0, 0, -1));
    p4->color = Color(0.5, 1, 0);
    p4->transparent = false;
    objects.push_back(org);
    objects.push_back(s);
    objects.push_back(s2);
//    objects.push_back(t);
//    objects.push_back(t2);
//    objects.push_back(t3);
//    objects.push_back(t4);
//    objects.push_back(p);
//    objects.push_back(p2);
//    objects.push_back(p3);
//    objects.push_back(p4);
    updatePoint();
}

void keyFunc(unsigned char ch, int x, int y) {
    if (ch == 'p')
        get();
    if (ch == 's')
        set();
}

int main(int argc, char** argv) {
    init(argc, argv);
    window = glutCreateWindow("OpenGL Project 5");
    initWindow();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyFunc);
    glutMainLoop();
    return 0;
}
