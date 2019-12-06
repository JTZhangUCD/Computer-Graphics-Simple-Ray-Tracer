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
struct Vector {
    float x,y,z;
    Vector(double x1, double y1, double z1) {
        x = x1; y = y1; z = z1;
    }
    Vector(void) {
        x = 0; y = 0; z = 0;
    }
    inline Vector operator + (const Vector& v1)const {
        return Vector(x+v1.x, y+v1.y, z+v1.z);
    }
    inline Vector operator - (void)const {
        return Vector(-x, -y, -z);
    }
    inline Vector operator - (const Vector& v1)const {
        return Vector(x-v1.x, y-v1.y, z-v1.z);
    }
    inline double operator * (const Vector& v1)const{
        return x * v1.x + y * v1.y + z * v1.z;
    }
    inline Vector operator * (const float f)const{
        return Vector(x * f, y * f, z * f);
    }
    inline Vector operator / (const float f)const{
        return Vector(x / f, y / f, z / f);
    }
    inline Vector &operator += (const Vector& v1) {
        x += v1.x; y += v1.y; z += v1.z;
        return *this;
    }
    inline Vector &operator -= (const Vector& v1) {
        x -= v1.x; y -= v1.y; z -= v1.z;
        return *this;
    }
    inline Vector &operator *= (float f) {
        x *= f; y *= f; z *= f;
        return *this;
    }
    inline Vector &operator /= (float f) {
        x /= f; y /= f; z /= f;
        return *this;
    }
    inline float length(void)const {
        return sqrt(x*x + y*y + z*z);
    }
    inline Vector& norm(void) {
        float length = sqrt(x*x + y*y + z*z);
        x /= length; y /= length; z /= length;
        return *this;
    }
};

Vector cross(const Vector& v1, const Vector& v2) {
    return Vector(v1.y*v2.z - v1.z*v2.y,v1.z*v2.x - v1.x*v2.z,v1.x*v2.y - v1.y*v2.x);
}

float dot(const Vector& v1, const Vector& v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

istream& operator>>(istream& in, Vector& p) {
    in >> p.x >> p.y >> p.z;
    return in;
}

ostream& operator<<(ostream& out, const Vector& p) {
    out << p.x << ' ' << p.y << ' ' << p.z;
    return out;
}
    
//============================================================================== COLOR

class Color3f {
public:
    float r, g, b;
    Color3f():r(0.0), g(0.0), b(0.0){}
    Color3f(float vr, float vg, float vb):r(vr),g(vg), b(vb){}
    Color3f& operator+=(const Color3f& rhs);
    Color3f& operator-=(const Color3f& rhs);
    Color3f& operator*=(const Vector& rhs);
    Color3f& operator*=(float a);
    Color3f& operator/=(float a);
};
    
Color3f& Color3f::operator+=(const Color3f& rhs) {
    r += rhs.r;
    g += rhs.g;
    b += rhs.b;
    return *this;
}

Color3f& Color3f::operator-=(const Color3f& rhs) {
    r -= rhs.r;
    g -= rhs.g;
    b -= rhs.b;
    return *this;
}

Color3f& Color3f::operator*=(const Vector& rhs) {
    r *= rhs.x;
    g *= rhs.y;
    b *= rhs.z;
    return *this;
}

Color3f& Color3f::operator*=(float a) {
    r *= a;
    g *= a;
    b *= a;
    return *this;
}

Color3f& Color3f::operator/=(float a) {
    r /= a;
    g /= a;
    b /= a;
    return *this;
}

Color3f operator+(const Color3f& l, const Color3f& rhs) {
    Color3f ret = l;
    ret += rhs;
    return ret;
}

Color3f operator-(const Color3f& l, const Color3f& rhs) {
    Color3f ret = l;
    ret -= rhs;
    return ret;
}

Color3f operator*(const Color3f& l, const Vector& r) {
    Color3f ret = l;
    ret *= r;
    return ret;
}

Color3f operator*(const Color3f& l, float a) {
    Color3f ret = l;
    ret *= a;
    return ret;
}

Color3f operator/(const Color3f& l, float a){
    Color3f ret = l;
    ret /= a;
    return ret;
}

istream& operator>>(istream& in, Color3f& p) {
    in >> p.r >> p.g >> p.b;
    return in;
}
ostream& operator<<(ostream& out, const Color3f& p) {
    out << p.r << ' ' << p.g << ' ' << p.b;
    return out;
}

//============================================================================== OBJECT
float iA = 0.2, iL = 1;
Vector kA(1, 1, 1), kD(1, 1, 1), kS(1, 1, 1);
int bigN = 5;
float bigK = 10;
float kR = 0.5, kT = 0.5;
float eta = 1.5;
Vector lightSource(0, 1, 0);

class Line {
public:
    Vector start;
    Vector u;
    bool inside;
public:
    Line(const Vector& start_, const Vector& u_, bool inside_):
        start(start_), u(u_), inside(inside_){}
    Line():inside(false){}
};

class Primitive {
public:
    bool transparent;
    Color3f color;
public:
    virtual bool intersect(const Line& l, Vector& result) = 0;
    virtual Vector normal(const Vector& p) = 0;
};

class Polygon : public Primitive {
public:
    vector<Vector> v;
    Vector n;
public:
    bool intersect(const Line& l, Vector& result);
    Vector normal(const Vector& p);
};

class Sphere : public Primitive {
public:
    Vector center;
    float radius;
public:
    bool intersect(const Line& l, Vector& result);
    Vector normal(const Vector& p);
};

class Quadratic : public Primitive {
public:
    float A, B, C, D, E, F, G, H, I, J;
public:
    bool intersect(const Line& l, Vector& result);
    Vector normal(const Vector& p);
};

inline void makePixel(int x, int y, const Color3f& c, GLfloat *pixels, int resolution) {
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

bool Polygon::intersect(const Line& l, Vector& result) {
    float D = -dot(n, v[0]);
    float s = -(D+dot(n, l.start)) / dot(n, l.u);
    bool ret = POS(s);
    result = l.start + l.u * s;
    for (int i = 0; i < v.size(); i++)
        if (ret)
            ret = ret && dot(cross(v[i % v.size()] - result, v[(i+1) % v.size()] - result), n) > 0;
    return ret;
}

Vector Polygon::normal(const Vector& p) {
    return n;
}

bool Sphere::intersect(const Line& l, Vector& result) {
    Vector diff = center - l.start;
    float base = dot(l.u, diff);
    float delta = base * base - dot(diff, diff) + radius * radius;
    if (delta < 0)
        return false;
    delta = sqrt(delta);
    Vector p1 = l.start + l.u * (base - delta);
    Vector p2 = l.start + l.u * (base + delta);
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

Vector Sphere::normal(const Vector& p) {
    Vector ret = (p - center).norm();
    return ret;
}

bool Quadratic::intersect(const Line& l, Vector& result) {
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
    Vector p1, p2;
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

Vector Quadratic::normal(const Vector& p) {
    Vector ret;
    ret.x = 2*A*p.x + D*p.y + E*p.z + G;
    ret.y = 2*B*p.y + D*p.x + F*p.z + H;
    ret.z = 2*C*p.z + E*p.x + F*p.y + I;
    ret.norm();
    return ret;
}

Polygon* newTriangle(const Vector& v1, const Vector& v2, const Vector& v3, const Vector& normal) {
    Polygon* ret = new Polygon();
    ret->v.push_back(v1);
    ret->v.push_back(v2);
    ret->v.push_back(v3);
    ret->n = normal;
    return ret;
}

Polygon* newSquare(const Vector& v1, const Vector& v2, const Vector& v3, const Vector& v4, const Vector& normal) {
    Polygon* ret = new Polygon();
    ret->v.push_back(v1);
    ret->v.push_back(v2);
    ret->v.push_back(v3);
    ret->v.push_back(v4);
    ret->n = normal;
    return ret;
}

Quadratic* newSphere(const Vector& center, float radius) {
    Quadratic* ret = new Quadratic();
    ret->A = ret->B = ret->C = 1;
    ret->D = ret->E = ret->F = 0;
    ret->G = - 2 * center.x;
    ret->H = - 2 * center.y;
    ret->I = - 2 * center.z;
    ret->J = dot(center, center) - radius * radius;
    return ret;
}

Quadratic* newEllipsoid(const Vector& center, float a, float b, float c) {
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

Vector refract(const Line& l, const Vector& normal) {
    float r = l.inside ? eta : 1 / eta;
    float cosThetaI = dot(-l.u, normal);
    float cosThetaR = sqrt(1 - (r * r) * (1 - cosThetaI * cosThetaI));
    return l.u * r - normal * (cosThetaR - r * cosThetaI);
}

Vector reflect(const Vector& u, const Vector& normal) {
    return u - normal * 2 * dot(normal, u);
}
    
bool intersect(const Line& line, const vector<Primitive*>& primitives, Vector& result, int& resultId) {
    Vector best, now;
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

Vector phong(const Vector& point, const Vector& ref, const Vector& normal, const vector<Primitive*>& primitives) {
    Vector intensity = kA * iA;
    Vector l = (lightSource - point).norm();
    Vector v = (ref - point).norm();
    Vector r = (normal * 2.0 * dot(normal, l) - l).norm();
    Vector diff;
    Line line; line.start = point; line.u = (lightSource - point).norm();
    Vector result;
    int resultId;
    if (intersect(line, primitives, result, resultId))
        return intensity;
    float productDiff = dot(l, normal);
    if (productDiff > 0) diff = kD * productDiff;
    Vector spec;
    float productSpec = dot(r, v);
    if (productDiff > 0 && productSpec > 0) spec = kS * productSpec;
    intensity += (diff + spec) * iL / ((ref - point).length() + bigK);
    return intensity;
}

Color3f light(const Line& line, const vector<Primitive*>& primitives, int depth) {
    Color3f ret;
    Primitive* prim;
    int id;
    Vector point;
    Vector normal;
    if (depth == 0) return ret;
    if (intersect(line, primitives, point, id)) {
        prim = primitives[id];
        normal = prim->normal(point);
        Vector ph = phong(point, line.start, normal, primitives);
        ret = prim->color * ph;
        Line refl;
        refl.start = point;
        refl.u = reflect(line.u, normal).norm();
        refl.inside = line.inside;
        Color3f rf = light(refl, primitives, depth - 1);
        ret += rf * kR;
        if (prim->transparent) {
            Line refr;
            refr.start = point;
            refr.u = refract(line, normal).norm();
            refr.inside = !line.inside;
            Color3f rt = light(refr, primitives, depth - 1);
            ret += rt * kT;
        }
    }
    return ret;
}

//============================================================================== MAIN

int window;

float canvas[WINDOW_SIZE * WINDOW_SIZE * 3];

vector<Primitive*> objects;

Vector from(0, -3, 0);
Vector at(0, 0, 0);
Vector up(0, 0, 1);
Vector rv;
float angDeg = 45;
int resolution = WINDOW_SIZE;
int depth = 5;

Line ray(int x, int y) {
    Line ret;
    ret.start = from;
    float xx = (float)x/(float)resolution*2 - 1;
    float yy = (float)y/(float)resolution*2 - 1;
    Vector p = from;
    p += rv * tan(angDeg / 180 * PI) * xx;
    p += up * tan(angDeg / 180 * PI) * yy;
    p += (at - from).norm();
    ret.u = (p - from).norm();
    return ret;
}

void updateVector() {
    up.norm();
    rv = cross(at - from, up).norm();
}

void initWindow() {
    glClearColor(0, 0, 0, 0);
}

void call() {
    1;
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    memset(canvas, 0.0f, sizeof(canvas));
    for (int x = 0; x < resolution; x++)
        for (int y = 0; y < resolution; y++) {
            Line r = ray(x, y);
            Color3f f = light(r, objects, depth);
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
    updateVector();
    cout << "OK!" << endl;
    refreshFunc();
}


void exitFunc() {
}

void init(int argc, char** argv) {
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);

    Quadratic* org = newSphere(Vector(0, 0, 0), 0.1);
    org->color = Color3f(1, 1, 1);
    org->transparent = true;
    Quadratic* s = newSphere(Vector(-0.5, 1, -0.5), 0.4);
    s->color = Color3f(1, 1, 1);
    s->transparent = false;
    Quadratic* s2 = newSphere(Vector(1, 1.5, 1), 0.6);
    s2->color = Color3f(1, 1, 1);
    s2->transparent = true;
    Polygon* t = newTriangle(Vector(-1, 2, 2),
                             Vector(-1, 2, 1),
                             Vector(-1, 1, 2),
                             Vector(-1, 0, 0));
    t->color = Color3f(1, 1, 1);
    t->transparent = true;
    Polygon* t2 = newTriangle(Vector(0, 2, 2),
                              Vector(-1, 2, 1),
                              Vector(-1, 2, 2),
                              Vector(0, 1, 0));
    t2->color = Color3f(1, 1, 1);
    t2->transparent = true;
    Polygon* t3 = newTriangle(Vector(-1, 2, 2),
                              Vector(-1, 1, 2),
                              Vector(0, 2, 2),
                              Vector(0, 0, 1));
    t3->color = Color3f(1, 1, 1);
    t3->transparent = true;
    Polygon* t4 = newTriangle(Vector(0, 2, 2),
                              Vector(-1, 1, 2),
                              Vector(-1, 2, 1),
                              Vector(1, -1, -1));
    t4->color = Color3f(1, 1, 1);
    t4->transparent = true;
    Polygon* p = newSquare(Vector(5, 6, -4), Vector(5, 6, 4),
                           Vector(-5, 6, 4), Vector(-5, 6, -4),
                           Vector(0, -1, 0));
    p->color = Color3f(1, 0.85, 0.25);
    p->transparent = false;
    Polygon* p2 = newSquare(Vector(-5, -6, -4), Vector(-5, 6, -4),
                           Vector(-5, 6, 4), Vector(-5, -6, 4),
                           Vector(1, 0, 0));
    p2->color = Color3f(0, 0.5, 1);
    p2->transparent = false;
    Polygon* p3 = newSquare(Vector(5, 6, -4), Vector(5, -6, -4),
                           Vector(5, -6, 4), Vector(5, 6, 4),
                           Vector(-1, 0, 0));
    p3->color = Color3f(1, 0.2, 0);
    p3->transparent = false;
    Polygon* p4 = newSquare(Vector(5, 6, 4), Vector(5, -6, 4),
                           Vector(-5, -6, 4), Vector(-5, 6, 4),
                           Vector(0, 0, -1));
    p4->color = Color3f(0.5, 1, 0);
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
    updateVector();
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
    atexit(exitFunc);
    glutMainLoop();
    return 0;
}

