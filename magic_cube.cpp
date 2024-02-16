#include <bits/stdc++.h>
#include <GL/glut.h>
#ifdef _WIN32
#include <windows.h>
#endif

#define _USE_MATH_DEFINES

using namespace std;

// Debug utility
// #define DEBUG 1
#ifdef DEBUG
#define dbg(args...) { debugPrint(#args, args); }
template <typename... Args>
void debugPrint(const std::string& arg_names, const Args&... args) {
    stringstream ss;
    ss << ">> " << arg_names << " = ";
    ((ss << args << ", "), ...);
    string s(ss.str());
    s.erase(s.end()-2);
    std::cerr << s << std::endl;
}
#else
#define dbg(args...) (void) "Leap of faith!"
#endif

inline double degToRad(double deg) {
    return deg * M_PI / 180.0;
}

inline double RadToDeg(double rad) {
    return rad * 180.0 / M_PI;
}

// Constants
const double PHI = RadToDeg(acos(-(1.0 / 3.0))); // Dihedral angle of an octahedron
const double EPSILON = 1e-10;
const double MAGIC_CUBE_SCALE_DELTA = 0.04;
const double ROTATION_ANGLE_DELTA = 10.0;
const double CAMERA_DELTA = 0.5;
const double TILT_ANGLE_DELTA = 1.0;

// 3D Point class
struct Point2D {
    double x, y, z;
};

// 3D Vector class
struct Vec3 {
    double x, y, z;

    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3() : x(0), y(0), z(0) {}

    void normalize() {
        double val = sqrt(x*x + y*y + z*z);
        x /= val, y /= val, z /= val;
    }

    Vec3 operator+(const Vec3& o) const {
        return {
            x + o.x,
            y + o.y,
            z + o.z
        };
    }

    Vec3 operator-(const Vec3& o) const {
        return {
            x - o.x,
            y - o.y,
            z - o.z
        };
    }

    Vec3 operator*(double scaler) const {
        return {
            x * scaler,
            y * scaler,
            z * scaler
        };
    }

    Vec3 cross(const Vec3& o) const {
        return {
            y * o.z - o.y * z,
            o.x * z - x * o.z,
            x * o.y - o.x * y
        };
    }

    double dot(const Vec3& o) const {
        return x * o.x + y * o.y + z * o.z;
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec3& vec);
};

std::ostream& operator<<(std::ostream& os, const Vec3& vec) {
    return os << "(" << vec.x << "," << vec.y << "," << vec.z << ")";
}

// Camera class for managing camera movements
struct Camera {
    Vec3 eye, at, up;
    Vec3 xc, yc, zc;


    Camera(Vec3 eye, Vec3 at, Vec3 up) : eye(eye), at(at), up(up) {
        updateCamera();
    }

    Camera() : eye(4, 4, 4), at(0, 0, 0), up(0, 1, 0) {
        updateCamera();
    }

    void updateCamera() {
        zc = eye - at;
        zc.normalize();

        xc = up.cross(zc);
        xc.normalize();

        yc = zc.cross(xc);

        dbg(eye, at, up);
        dbg(xc, yc, zc);
    }

    void setupCamera() {
        gluLookAt(
            eye.x, eye.y, eye.z,
            at.x, at.y, at.z,
            up.x, up.y, up.z
        );
    }

    void moveUp() {
        eye = eye + yc * CAMERA_DELTA;
        at = at + yc * CAMERA_DELTA;
        updateCamera();
    }

    void moveDown() {
        eye = eye - yc * CAMERA_DELTA;
        at = at - yc * CAMERA_DELTA;
        updateCamera();
    }

    void moveLeft() {
        eye = eye - xc * CAMERA_DELTA;
        at = at - xc * CAMERA_DELTA;
        updateCamera();
    }

    void moveRight() {
        eye = eye + xc * CAMERA_DELTA;
        at = at + xc * CAMERA_DELTA;
        updateCamera();
    }

    void moveForward() {
        eye = eye - zc * CAMERA_DELTA;
        updateCamera();
    }

    void moveBackward() {
        eye = eye + zc * CAMERA_DELTA;
        updateCamera();
    }

    void lookLeft() {
        at = at - xc * CAMERA_DELTA;
        updateCamera();
    }

    void lookRight() {
        at = at + xc * CAMERA_DELTA;
        updateCamera();
    }

    void lookUp() {
        at = at + yc * CAMERA_DELTA;
        updateCamera();
    }

    void lookDown() {
        at = at - yc * CAMERA_DELTA;
        updateCamera();
    }

    void moveUpWithoutChangingReference() {
        eye = eye + up * CAMERA_DELTA;
        updateCamera();
    }

    void moveDownWithoutChangingReference() {
        eye = eye - up * CAMERA_DELTA;
        updateCamera();
    }

    void rotateCCW() {
        double tiltAngle = -TILT_ANGLE_DELTA;
        up = up * cos(degToRad(tiltAngle)) + zc.cross(up) * sin(degToRad(tiltAngle)) + zc * zc.dot(up) * (1 - cos(degToRad(tiltAngle)));
        updateCamera();
    }


    void rotateCW() {
        double tiltAngle = TILT_ANGLE_DELTA;
        up = up * cos(degToRad(tiltAngle)) + zc.cross(up) * sin(degToRad(tiltAngle)) + zc * zc.dot(up) * (1 - cos(degToRad(tiltAngle)));
        updateCamera();
    }
};

// Global variables
double globalScale = 2.0;
double magicCubeScale = 1.0;
double rotationAngle = 0.0;
bool showAxes = false;
bool showMagicCube = true;
double axesLength = 1.0;
Camera camera;

// Draw axes
void drawAxes() {
    glBegin(GL_LINES); {
        // X axis
        glColor3ub(255, 0, 0); // Red
        glVertex3f(0, 0, 0);
        glVertex3f(axesLength, 0, 0);

        // Y axis
        glColor3ub(0, 255, 0); // Green
        glVertex3f(0, 0, 0);
        glVertex3f(0, axesLength, 0);

        // Z axis
        glColor3ub(0, 0, 255); // Blue
        glVertex3f(0, 0, 0);
        glVertex3f(0, 0, axesLength);
    } glEnd();
}

// Draw the triangular face of an octahedron 
// with vertices at (1, 0, 0), (0, 1, 0) and (0, 0, 1)
void drawOctahedronBaseFace() {
    glBegin(GL_TRIANGLES); {
        glVertex3d(1, 0, 0);
        glVertex3d(0, 1, 0);
        glVertex3d(0, 0, 1);
    } glEnd();
}

// Draw octahedron with 8 scaled faces
void drawOctahedron() {
    // (+1, +1, +1) face
    glPushMatrix(); {
        glColor3ub(0, 255, 255);
        glScaled(1.0, 1.0, 1.0);
        glTranslated((1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0);
        glScaled(magicCubeScale, magicCubeScale, magicCubeScale);
        drawOctahedronBaseFace();
    } glPopMatrix();

    // (+1, +1, -1) face
    glPushMatrix(); {
        glColor3ub(255, 0, 255);
        glScaled(1.0, 1.0, -1.0);
        glTranslated((1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0);
        glScaled(magicCubeScale, magicCubeScale, magicCubeScale);
        drawOctahedronBaseFace();
    } glPopMatrix();

    // (+1, -1, +1) face
    glPushMatrix(); {
        glColor3ub(255, 0, 255);
        glScaled(1.0, -1.0, 1.0);
        glTranslated((1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0);
        glScaled(magicCubeScale, magicCubeScale, magicCubeScale);
        drawOctahedronBaseFace();
    } glPopMatrix();

    // (+1, -1, -1) face
    glPushMatrix(); {
        glColor3ub(0, 255, 255);
        glScaled(1.0, -1.0, -1.0);
        glTranslated((1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0);
        glScaled(magicCubeScale, magicCubeScale, magicCubeScale);
        drawOctahedronBaseFace();
    } glPopMatrix();

    // (-1, +1, +1) face
    glPushMatrix(); {
        glColor3ub(255, 0, 255);
        glScaled(-1.0, 1.0, 1.0);
        glTranslated((1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0);
        glScaled(magicCubeScale, magicCubeScale, magicCubeScale);
        drawOctahedronBaseFace();
    } glPopMatrix();

    // (-1, +1, -1) face
    glPushMatrix(); {
        glColor3ub(0, 255, 255);
        glScaled(-1.0, 1.0, -1.0);
        glTranslated((1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0);
        glScaled(magicCubeScale, magicCubeScale, magicCubeScale);
        drawOctahedronBaseFace();
    } glPopMatrix();

    // (-1, -1, +1) face
    glPushMatrix(); {
        glColor3ub(0, 255, 255);
        glScaled(-1.0, -1.0, 1.0);
        glTranslated((1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0);
        glScaled(magicCubeScale, magicCubeScale, magicCubeScale);
        drawOctahedronBaseFace();
    } glPopMatrix();

    // (-1, -1, -1) face
    glPushMatrix(); {
        glColor3ub(255, 0, 255);
        glScaled(-1.0, -1.0, -1.0);
        glTranslated((1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0, (1 - magicCubeScale) / 3.0);
        glScaled(magicCubeScale, magicCubeScale, magicCubeScale);
        drawOctahedronBaseFace();
    } glPopMatrix();
}

// Draw the base cylindrical face in (+1, +1, 0) plane
void drawCylinderBaseFace() {
    double r = (1 - magicCubeScale) / sqrt(3);
    double h = sqrt(2) * magicCubeScale;
    double a1 = PHI / 2;
    double a2 = 180.0 - (PHI / 2);
    int seg = 100;
    double d = (a2 - a1) / seg;

    glPushMatrix(); {
        glTranslated(0, magicCubeScale, 0);
        glRotated(-45.0, 0, 0, 1);

        glBegin(GL_QUADS); {
            double py = r * sin(degToRad(a1));
            double pz = r * cos(degToRad(a1));
            for (int i = 1; i <= seg; ++i) {
                double a = a1 + d * i;
                double cy = r * sin(degToRad(a));
                double cz = r * cos(degToRad(a));

                glVertex3d(0, py, pz);
                glVertex3d(h, py, pz);
                glVertex3d(h, cy, cz);
                glVertex3d(0, cy, cz);

                py = cy;
                pz = cz;
            }
        } glEnd();
    } glPopMatrix();
}

// Draw all 12 cylindrical faces
void drawCylinders() {
    double h = sqrt(2) * magicCubeScale;

    glColor3ub(255, 255, 0);

    // (+1, +1, 0) face
    drawCylinderBaseFace();

    // (-1, +1, 0) face
    glPushMatrix(); {
        glScaled(-1, 1, 1);
        drawCylinderBaseFace();
    } glPopMatrix();


    // (-1, -1, 0) face
    glPushMatrix(); {
        glScaled(-1, -1, 1);
        drawCylinderBaseFace();
    } glPopMatrix();

    // (+1, -1, 0) face
    glPushMatrix(); {
        glScaled(1, -1, 1);
        drawCylinderBaseFace();
    } glPopMatrix();

    // (0, +1, +1) face
    glPushMatrix(); {
        glRotated(-90, 0, 1, 0);
        drawCylinderBaseFace();
    } glPopMatrix();

    // (0, +1, -1) face
    glPushMatrix(); {
        glScaled(1, 1, -1);
        glRotated(-90, 0, 1, 0);
        drawCylinderBaseFace();
    } glPopMatrix();

    // (0, -1, -1) face
    glPushMatrix(); {
        glScaled(1, -1, -1);
        glRotated(-90, 0, 1, 0);
        drawCylinderBaseFace();
    } glPopMatrix();

    // (0, -1, +1) face
    glPushMatrix(); {
        glScaled(1, -1, 1);
        glRotated(-90, 0, 1, 0);
        drawCylinderBaseFace();
    } glPopMatrix();

    // (+1, 0, +1) face
    glPushMatrix(); {
        glTranslated(0, 0, h);
        glRotated(45.0, 0, 1, 0);
        glRotated(90.0, 1, 0, 0);
        glRotated(45.0, 0, 0, 1);
        glTranslated(0, -h, 0);
        drawCylinderBaseFace();
    } glPopMatrix();

    // (+1, 0, -1) face
    glPushMatrix(); {
        glScaled(1, 1, -1);
        glTranslated(0, 0, h);
        glRotated(45.0, 0, 1, 0);
        glRotated(90.0, 1, 0, 0);
        glRotated(45.0, 0, 0, 1);
        glTranslated(0, -h, 0);
        drawCylinderBaseFace();
    } glPopMatrix();


    // (-1, 0, -1) face
    glPushMatrix(); {
        glScaled(-1, 1, -1);
        glTranslated(0, 0, h);
        glRotated(45.0, 0, 1, 0);
        glRotated(90.0, 1, 0, 0);
        glRotated(45.0, 0, 0, 1);
        glTranslated(0, -h, 0);
        drawCylinderBaseFace();
    } glPopMatrix();

    // (-1, 0, +1) face
    glPushMatrix(); {
        glScaled(-1, 1, 1);
        glTranslated(0, 0, h);
        glRotated(45.0, 0, 1, 0);
        glRotated(90.0, 1, 0, 0);
        glRotated(45.0, 0, 0, 1);
        glTranslated(0, -h, 0);
        drawCylinderBaseFace();
    } glPopMatrix();
}

// Generate vertices for +X face only by intersecting 2 circular planes
// (longitudinal and latitudinal) at the given longitude/latitude angles
// Reference: https://www.songho.ca/opengl/gl_sphere.html
std::vector<std::vector<Point2D>> buildUnitPositiveXSphereFace(int subdivision) {
    const double DEG2RAD = acos(-1) / 180.0f;

    std::vector<std::vector<Point2D>> vertices;
    double n1[3];        // normal of longitudinal plane rotating along Y-axis
    double n2[3];        // normal of latitudinal plane rotating along Z-axis
    double v[3];         // direction vector intersecting 2 planes, n1 x n2
    double a1;           // longitudinal angle along Y-axis
    double a2;           // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (1 << subdivision) + 1;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for(int i = 0; i < pointsPerRow; ++i) {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        vertices.push_back(std::vector<Point2D>());
        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for(int j = 0; j < pointsPerRow; ++j) {
            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // find direction vector of intersected line, n1 x n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            // normalize direction vector
            float scale = 1 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;

            // add a vertex into array
            vertices.back().push_back(Point2D{v[0], v[1], v[2]});
        }
    }

    return vertices;
}

// Draw base spherical face in (+1, 0, 0) direction with radius 1
void drawSphereBaseFace() {
    int sd = 5;
    int n = (1 << sd) + 1;
    std::vector<std::vector<Point2D>> v = buildUnitPositiveXSphereFace(sd);

    glBegin(GL_QUADS); {
        for (int i = 0; i < n-1; ++i) {
            for (int j = 0; j < n-1; ++j) {
                Point2D a = v[i][j], b = v[i+1][j], c = v[i+1][j+1], d = v[i][j+1];

                glVertex3d(a.x, a.y, a.z);
                glVertex3d(b.x, b.y, b.z);
                glVertex3d(c.x, c.y, c.z);
                glVertex3d(d.x, d.y, d.z);
            }
        }
    } glEnd();
}

// Draw all 6 sperical faces
void drawSphere() {
    double r = (1 - magicCubeScale) / sqrt(3);

    // (+1, 0, 0) face
    glPushMatrix(); {
        glColor3ub(255, 0, 0);
        glTranslated(magicCubeScale, 0, 0);
        glScaled(r, r, r);
        drawSphereBaseFace();
    } glPopMatrix();

    // (-1, 0, 0) face
    glPushMatrix(); {
        glColor3ub(255, 0, 0);
        glScaled(-1, 1, 1);
        glTranslated(magicCubeScale, 0, 0);
        glScaled(r, r, r);
        drawSphereBaseFace();
    } glPopMatrix();

    // (0, +1, 0) face
    glPushMatrix(); {
        glColor3ub(0, 255, 0);
        glRotated(90.0, 0, 0, 1);
        glTranslated(magicCubeScale, 0, 0);
        glScaled(r, r, r);
        drawSphereBaseFace();
    } glPopMatrix();

    // (0, -1, 0) face
    glPushMatrix(); {
        glColor3ub(0, 255, 0);
        glScaled(1, -1, 1);
        glRotated(90.0, 0, 0, 1);
        glTranslated(magicCubeScale, 0, 0);
        glScaled(r, r, r);
        drawSphereBaseFace();
    } glPopMatrix();

    // (0, 0, +1) face
    glPushMatrix(); {
        glColor3ub(0, 0, 255);
        glRotated(-90.0, 0, 1, 0);
        glTranslated(magicCubeScale, 0, 0);
        glScaled(r, r, r);
        drawSphereBaseFace();
    } glPopMatrix();

    // (0, 0, -1) face
    glPushMatrix(); {
        glColor3ub(0, 0, 255);
        glScaled(1, 1, -1);
        glRotated(-90.0, 0, 1, 0);
        glTranslated(magicCubeScale, 0, 0);
        glScaled(r, r, r);
        drawSphereBaseFace();
    } glPopMatrix();
}

// Draw magic cube
void drawMagicCube() {
    drawOctahedron();
    drawCylinders();
    drawSphere();
}

// Handler for window-repaint event.
// Call back when the window first appears
// and whenever the window needs to be re-painted.
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear

    glMatrixMode(GL_MODELVIEW); // To operate on Model-View matrix
    glLoadIdentity(); // Reset the model-view matrix

    camera.setupCamera(); // Setup camera position

    // Perform camera tilt
    // glTranslated(-camera.eye.x, -camera.eye.y, -camera.eye.z);
    // glRotated(tiltAngle, camera.zc.x, camera.zc.y, camera.zc.z);
    // glTranslated(camera.eye.x, camera.eye.y, camera.eye.z);

    // Perform object rotation
    glRotated(rotationAngle, 0, 1, 0);

    // Perform global scaling
    glScaled(globalScale, globalScale, globalScale);
    

    if (showAxes) {
        drawAxes();
    }
    if (showMagicCube) {
        drawMagicCube();
    }

    glutSwapBuffers(); // Render
}

// Handler for window re-size event.
// Called back when the window first appears 
// and whenever the window is re-sized with its new width and height.
void reshapeListener(int width, int height) {
    // Compute new aspect ratio
    if (height == 0) {
        height = 1; // To prevent divide by 0
    }
    double aspect = (double) width / (double) height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity(); // Reset the projection matrix

    // // For orthographic projection
    // if (width >= height) {
    //     // aspect >= 1, set the height from -1 to 1, with larger width
    //     gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    // } else {
    //     // aspect < 1, set the width to -1 to 1, with larger height
    //     gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    // }

    // Enable perspective projection with fovy, aspect, zNear and zFar
    gluPerspective(45.0, aspect, 0.1, 100.0);
}

// Callback handler for normal-key event
void keyboardListener(unsigned char key, [[maybe_unused]] int x, [[maybe_unused]] int y) {
    switch (key) {
    case 'a':
        std::cerr << "Rotating object clockwise." << std::endl;
        rotationAngle -= ROTATION_ANGLE_DELTA;
        break;

    case 'd':
        std::cerr << "Rotating object counter-clockwise." << std::endl;
        rotationAngle += ROTATION_ANGLE_DELTA;
        break;

    case 'w':
        std::cerr << "Moving up without changing reference point." << std::endl;
        camera.moveUpWithoutChangingReference();
        break;

    case 's':
        std::cerr << "Moving down without changing reference point." << std::endl;
        camera.moveDownWithoutChangingReference();
        break;

    case '1':
        std::cerr << "Looking left." << std::endl;
        camera.lookLeft();
        break;

    case '2':
        std::cerr << "Looking right." << std::endl;
        camera.lookRight();
        break;

    case '3':
        std::cerr << "Looking up." << std::endl;
        camera.lookUp();
        break;

    case '4':
        std::cerr << "Looking down." << std::endl;
        camera.lookDown();
        break;

    case '5':
        std::cerr << "Tilting counter-clockwise." << std::endl;
        camera.rotateCCW();
        break;

    case '6':
        std::cerr << "Tilting clockwise." << std::endl;
        camera.rotateCW();
        break;

    case ',': // Morph from octahedron to sphere
        magicCubeScale = std::max(0.0, magicCubeScale - MAGIC_CUBE_SCALE_DELTA);
        if (magicCubeScale < EPSILON) {
            magicCubeScale = 0.0;
        }
        std::cerr << "Magic Cube Scale = " << magicCubeScale << std::endl;
        break;

    case '.': // Morph from sphere to octahedron
        magicCubeScale = std::min(1.0, magicCubeScale + MAGIC_CUBE_SCALE_DELTA);
        if (abs(magicCubeScale - 1) < EPSILON) {
            magicCubeScale = 1.0;
        }
        std::cerr << "Magic Cube Scale = " << magicCubeScale << std::endl;
        break;

    case 'c':
        showMagicCube = !showMagicCube;
        break;

    case 'x':
        showAxes = !showAxes;
        break;

    case '[':
        globalScale = std::max(0.0, globalScale - 0.1);
        if (globalScale < EPSILON) {
            globalScale = 0.0;
        }
        std::cerr << "Global Scale = " << globalScale << std::endl;
        break;

    case ']':
        globalScale = globalScale + 0.1;
        std::cerr << "Global Scale = " << globalScale << std::endl;
        break;

    case ';':
        axesLength = std::max(1.0, axesLength - 0.1);
        std::cerr << "Axes length = " << axesLength << std::endl;
        break;

    case '\'':
        axesLength = axesLength + 0.1;
        std::cerr << "Axes length = " << axesLength << std::endl;
        break;

    case 27: // ESC key
        exit(0);
        break;

    default:
        return;
    }
    glutPostRedisplay(); // Post a paint request to activate display()
}

// Callback handler for special-key event
void specialKeyListener(int key, [[maybe_unused]] int x, [[maybe_unused]] int y) {
    switch (key) {
    case GLUT_KEY_PAGE_UP:
        std::cerr << "Moving up." << std::endl;
        camera.moveUp();
        break;

    case GLUT_KEY_PAGE_DOWN:
        std::cerr << "Moving down." << std::endl;
        camera.moveDown();
        break;

    case GLUT_KEY_LEFT:
        std::cerr << "Moving left." << std::endl;
        camera.moveLeft();
        break;

    case GLUT_KEY_RIGHT:
        std::cerr << "Moving right." << std::endl;
        camera.moveRight();
        break;

    case GLUT_KEY_UP:
        std::cerr << "Moving forward." << std::endl;
        camera.moveForward();
        break;

    case GLUT_KEY_DOWN:
        std::cerr << "Moving backward." << std::endl;
        camera.moveBackward();
        break;

    default:
        return;
    }
    glutPostRedisplay(); // Post a paint request to activate display()
}

// Initialize OpenGL
void initGL() {
    glClearColor(0, 0, 0, 1); // Set background color
    glEnable(GL_DEPTH_TEST); // Enable depth testing for z-culling
}

//  Main function: GLUT runs as a console application starting at main()
int main(int argc, char** argv) {
    glutInit(&argc, argv);                      // Initialize GLUT
    glutInitWindowSize(640, 640);               // Set the window's initial width & height
    glutInitWindowPosition(50, 50);             // Position the window's initial top-left corner
	glutInitDisplayMode(
        GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	// Depth, Double buffer, RGB color
    glutCreateWindow("Magic Cube");             // Create a window with the given title
    glutDisplayFunc(display);                   // Register display callback handler for window re-paint
    glutReshapeFunc(reshapeListener);           // Register callback handler for window re-shape
    glutKeyboardFunc(keyboardListener);         // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);        // Register callback handler for special-key event
    initGL();                                   // Our own OpenGL initialization
    glutMainLoop();                             // Enter the event-processing loop
    return 0;
}
