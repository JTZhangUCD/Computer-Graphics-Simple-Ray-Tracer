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
    Vector(float x1, float y1, float z1) {
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
        return x*v1.x + y*v1.y + z*v1.z;
    }
    inline Vector operator * (const float f)const{
        return Vector(x*f, y*f, z*f);
    }
    inline Vector operator / (const float f)const{
        return Vector(x/f, y/f, z/f);
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

struct Color {
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
    inline Color operator * (const Vector& v1)const {
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
    inline Color &operator *= (const Vector& v1) {
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

//============================================================================== OBJECT
float iA = 0.2, iL = 1;
Vector kA(1, 1, 1), kD(1, 1, 1), kS(1, 1, 1);
int bigN = 1;
float bigK = 10;
float kR = 0.5, kT = 0.5;
float eta = 1.5; // refraction coef
Vector lightSource(0, 1, 0);
//float iA = 0.1, iL = 1;
//Vector kA(1, 1, 1), kD(1, 1, 1), kS(.75, .75, .75);
//int bigN = 1;
//float bigK = 2;
//float kR = 0.5, kT = 0.5;
//float eta = 1.5; // refraction coef
//Vector lightSource(0, 0, 1);
    
struct Line { // ray
    Vector start; // starting point
    Vector u; // unit-direction
    bool inside; // is inside of an object
    Line(const Vector& start1, const Vector& u1, bool inside1) {
        start = start1; u = u1; inside = inside1;
    }
    Line(void) {
        inside = false;
    }
};

class Primitive {
public:
    bool transparent; // is transparent
    Color color; // color of the object
public:
    virtual bool intersect(const Line& l, Vector& result) = 0; // check if intersect a line
    virtual Vector normal(const Vector& p) = 0; // get normal vector of a point
};

class Polygon : public Primitive { //这是一个面
public:
    vector<Vector> v; // list of vertices
    Vector n;
public:
    bool intersect(const Line& l, Vector& result) {
        float D = -dot(n, v[0]);
        float s = -(D+dot(n, l.start)) / dot(n, l.u);
        bool ret = POS(s);
        result = l.start + l.u * s;
        for (int i = 0; i < v.size(); i++)
            if (ret)
                ret = ret && dot(cross(v[i % v.size()] - result, v[(i+1) % v.size()] - result), n) > 0;
        return ret;
    }
    Vector normal(const Vector& p) {
        return n;
    }
};

class Sphere : public Primitive {
public:
    Vector center;
    float radius;
public:
    bool intersect(const Line& l, Vector& result) {
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
    Vector normal(const Vector& p) {
        Vector ret = (p - center).norm();
        return ret;
    }
};

class Quadratic : public Primitive {
public:
    float A, B, C, D, E, F, G, H, I, J;
public:
    bool intersect(const Line& l, Vector& result) {
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
    Vector normal(const Vector& p) {
        Vector ret;
        ret.x = 2*A*p.x + D*p.y + E*p.z + G;
        ret.y = 2*B*p.y + D*p.x + F*p.z + H;
        ret.z = 2*C*p.z + E*p.x + F*p.y + I;
        ret.norm();
        return ret;
    }
};

vector<Primitive*> objects; //=================================================================== global
int window;

float buffer[WINDOW_SIZE * WINDOW_SIZE * 3];
//vector<Color> pixelBuffer;
    
Vector from(0, -3, 0);
Vector at(0, 0, 0);
Vector up(0, 0, 1);
Vector rv;
float angDeg = 45;
int resolution = WINDOW_SIZE;
int depth = 5;

inline void makePixel(int x, int y, const Color& c) {
    buffer[(y * resolution + x) * 3] = c.r;
    buffer[(y * resolution + x) * 3 + 1] = c.g;
    buffer[(y * resolution + x) * 3 + 2] = c.b;
}
inline void normalize(void) {
    // best: 最亮的点
    GLfloat best = 0.0;
    for (int i = 0; i < resolution * resolution * 3; i++) {
        if (buffer[i] > best)
            best = buffer[i];
    }
    if (best > 0.0) {
        for (int i = 0; i < resolution * resolution * 3; i++)
            buffer[i] /= best;
    }
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
    
bool intersect(const Line& line, Vector& result, int& resultId) {
    Vector best, now;
    int bestId;
    float bestDist = 0, nowDist;
    bool ret = false;
    for (int i = 0; i < objects.size(); i++) {
        if (objects[i]->intersect(line, now)) {
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

Vector phong(const Vector& point, const Vector& ref, const Vector& normal) {
//    Vector intensity = kA * iA;
//    Vector l = (lightSource - point).norm();
//    Vector v = (ref - point).norm();
//    Vector r = (normal * 2.0 * dot(normal, l) - l).norm();
//    Vector diff;
//    Line line; line.start = point; line.u = (lightSource - point).norm();
//    Vector result;
//    int resultId;
//    if (intersect(line, objects, result, resultId))
//        return intensity;
//    float productDiff = dot(l, normal);
//    if (productDiff > 0) diff = kD * productDiff;
//    Vector spec;
//    float productSpec = dot(r, v);
//    if (productDiff > 0 && productSpec > 0) spec = kS * productSpec;
//    intensity += (diff + spec) * iL / ((ref - point).length() + bigK);
//    return intensity;
    
    Vector intensity;
    intensity += kA * iA;

    Vector x = lightSource;
    Vector p = point;
    Vector f = ref;
    float K = bigK;
    float coef = iL / (K + (f - p).length());

    Vector l = (x - p).norm();
    Vector n = normal;
    Vector r = -l + n * (2.0 * (l * n));
    Vector v = (f - p).norm();

    Line line = Line(p,l,false);
    Vector result;
    int resultId;
    if (intersect(line, result, resultId)) {
        return intensity;
    }

    if (n * v < 0) {
        // do nothing
    } else if ((n * l > 0 && n * v < 0) || (n * l < 0 && n * v > 0)) {
        // do nothing
    } else if (r * v < 0) {
        intensity += kD * (l * n) * coef;
    } else {
        intensity += kD * (l * n) * coef + kS * pow(r * v, bigN);
    }

    return intensity;
}

Color light(const Line& line, int depth) {
    Color ret;
    Primitive* prim;
    int id;
    Vector point;
    Vector normal;
    if (depth == 0) return ret;
    if (intersect(line, point, id)) {
        prim = objects[id];
        normal = prim->normal(point);
        Vector ph = phong(point, line.start, normal);
        ret = prim->color * ph;
        Line refl;
        refl.start = point;
        refl.u = reflect(line.u, normal).norm();
        refl.inside = line.inside;
        Color rf = light(refl, depth - 1);
        ret += rf * kR;
        if (prim->transparent) {
            Line refr;
            refr.start = point;
            refr.u = refract(line, normal).norm();
            refr.inside = !line.inside;
            Color rt = light(refr, depth - 1);
            ret += rt * kT;
        }
    }
    return ret;
}

//============================================================================== MAIN

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
    
void display() {
    glClear(GL_COLOR_BUFFER_BIT);
//    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
//    glLoadIdentity();
    memset(buffer, 0.0f, sizeof(buffer));
//    pixelBuffer.clear();
//    for (int i = 0; i < resolution * resolution; i++) {
//        pixelBuffer.at(i).r = 0.0f;
//        pixelBuffer.at(i).g = 0.0f;
//        pixelBuffer.at(i).b = 0.0f;
//    }
    for (int x = 0; x < resolution; x++) {
        for (int y = 0; y < resolution; y++) {
            Line r = ray(x, y);
            Color f = light(r, depth);
            makePixel(x, y, f);
//            pixelBuffer.push_back(f);
//            pixelBuffer.at(x * resolution + y).r = f.r;
//            pixelBuffer.at(x * resolution + y).g = f.g;
//            pixelBuffer.at(x * resolution + y).b = f.b;
        }
    }
    normalize();
    GLfloat best = 0.0;
    for (int i = 0; i < resolution * resolution * 3; i++) {
        if (buffer[i] > best)
            best = buffer[i];
    }
    if (best > 0.0) {
        for (int i = 0; i < resolution * resolution * 3; i++)
            buffer[i] /= best;
    }
//    for (int i = 0; i < resolution*resolution; i++) {
//        if (pixelBuffer.at(i).r > best) {
//            best = pixelBuffer.at(i).r;
//        }
//        if (pixelBuffer.at(i).g > best) {
//            best = pixelBuffer.at(i).g;
//        }
//        if (pixelBuffer.at(i).b > best) {
//            best = pixelBuffer.at(i).b;
//        }
//    }
//    cout << "best:" << best << endl;

//    if (best > 0.0) {
//        for (int i = 0; i < resolution*resolution; i++) {
//            pixelBuffer.at(i).r /= best;
//            pixelBuffer.at(i).g /= best;
//            pixelBuffer.at(i).b /= best;
//        }
//    }
    glDrawPixels(resolution, resolution, GL_RGB, GL_FLOAT, buffer);
//    for (int x = 0; x < resolution; x++) {
//        for (int y = 0; y < resolution; y++) {
//            // draw pixel
//            glPointSize(1);
//            Color cur = pixelBuffer.at(x * resolution + y);
//            glBegin(GL_POINTS);
//            glColor3f(cur.r,cur.g,cur.b);
//            glVertex2f(x,y);
//            glEnd();
//        }
//    }
//    cout << "Swap" << endl;
    glutSwapBuffers();
}

void refreshFunc() {
//    glutPostWindowRedisplay(window);
    glutPostRedisplay();
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

void init(int argc, char** argv) {
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(600, 600);
//    glutInitWindowPosition(100, 100);
    glutInitWindowPosition(600, 600);

    Quadratic* org = newSphere(Vector(0, 0, 0), 0.1);
    org->color = Color(1, 1, 1);
    org->transparent = true;
    Quadratic* s = newSphere(Vector(-0.5, 1, -0.5), 0.4);
    s->color = Color(1, 1, 1);
    s->transparent = false;
    Quadratic* s2 = newSphere(Vector(1, 1.5, 1), 0.6);
    s2->color = Color(1, 1, 1);
    s2->transparent = true;
    Polygon* t = newTriangle(Vector(-1, 2, 2),
                             Vector(-1, 2, 1),
                             Vector(-1, 1, 2),
                             Vector(-1, 0, 0));
    t->color = Color(1, 1, 1);
    t->transparent = true;
    Polygon* t2 = newTriangle(Vector(0, 2, 2),
                              Vector(-1, 2, 1),
                              Vector(-1, 2, 2),
                              Vector(0, 1, 0));
    t2->color = Color(1, 1, 1);
    t2->transparent = true;
    Polygon* t3 = newTriangle(Vector(-1, 2, 2),
                              Vector(-1, 1, 2),
                              Vector(0, 2, 2),
                              Vector(0, 0, 1));
    t3->color = Color(1, 1, 1);
    t3->transparent = true;
    Polygon* t4 = newTriangle(Vector(0, 2, 2),
                              Vector(-1, 1, 2),
                              Vector(-1, 2, 1),
                              Vector(1, -1, -1));
    t4->color = Color(1, 1, 1);
    t4->transparent = true;
    Polygon* p = newSquare(Vector(5, 6, -4), Vector(5, 6, 4),
                           Vector(-5, 6, 4), Vector(-5, 6, -4),
                           Vector(0, -1, 0));
    p->color = Color(1, 0.85, 0.25);
    p->transparent = false;
    Polygon* p2 = newSquare(Vector(-5, -6, -4), Vector(-5, 6, -4),
                           Vector(-5, 6, 4), Vector(-5, -6, 4),
                           Vector(1, 0, 0));
    p2->color = Color(0, 0.5, 1);
    p2->transparent = false;
    Polygon* p3 = newSquare(Vector(5, 6, -4), Vector(5, -6, -4),
                           Vector(5, -6, 4), Vector(5, 6, 4),
                           Vector(-1, 0, 0));
    p3->color = Color(1, 0.2, 0);
    p3->transparent = false;
    Polygon* p4 = newSquare(Vector(5, 6, 4), Vector(5, -6, 4),
                           Vector(-5, -6, 4), Vector(-5, 6, 4),
                           Vector(0, 0, -1));
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
    updateVector();
    
//    for (int i = 0; i < resolution * resolution; i++) {
//        pixelBuffer.push_back(Color());
//    }
}

void keyFunc(unsigned char ch, int x, int y) {
    if (ch == 'p')
        get();
    if (ch == 's')
        set();
    glutPostRedisplay();
}

void idle() {
    glutPostRedisplay();
}
    
int main(int argc, char** argv) {
    init(argc, argv);
    window = glutCreateWindow("OpenGL Project 5");
    initWindow();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyFunc);
    glutIdleFunc(idle);
    glClearColor(0,0,0,0);
    glutMainLoop();
    return 0;
}
