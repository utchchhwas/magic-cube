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

inline int getCurrentTime() {
    time_t now = time(0);
    struct tm* tm = localtime(&now);
    return tm->tm_sec + tm->tm_min * 60 + tm->tm_hour * 3600;
}

inline int64_t getCurrentTimeMilli() {
    auto time = std::chrono::system_clock::now();
    auto since_epoch = time.time_since_epoch();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(since_epoch);
    return millis.count();
}

inline double degToRad(double deg) {
    return deg * M_PI / 180.0f;
}

struct Color {
    int r, g, b;
};

inline void glSetColor(Color c) {
    glColor3f((double) c.r / 255, (double) c.g / 255, (double) c.b / 255);
}

struct Point2D {
    double x, y;
};

void drawPoint(Point2D p) {
    glBegin(GL_POINTS); {
        glVertex2d(p.x, p.y);
    } glEnd();
}

void drawEmptyCircle(Point2D c, double r, int seg) {
    glBegin(GL_LINE_LOOP); {
        for (int i = 0; i < seg; ++i) {
            double a = (double) i * 360.0 / (double) seg;
            double x = r * cosf(degToRad(a));
            double y = r * sinf(degToRad(a));

            glVertex2d(c.x + x, c.y + y);
        }
    } glEnd();
}

void drawFilledCircle(Point2D c, double r, int seg) {
    glBegin(GL_TRIANGLE_FAN); {
        glVertex2d(c.x, c.y);
        for (int i = 0; i <= seg; ++i) {
            double a = (double) i * 360.0 / (double) seg;
            double x = r * cosf(degToRad(a));
            double y = r * sinf(degToRad(a));

            glVertex2f(c.x + x, c.y + y);
        }
    } glEnd();
}

void drawClockMarks(Point2D c, double r) {
    int d = 6;
    for (int a = 0; a < 360; a += d) {
        double x = c.x + r * cosf(degToRad(a));
        double y = c.y + r * sinf(degToRad(a));

        glColor3ub(242, 255, 64);
        if (a % 30 == 0) {
            glPushMatrix(); {
                glTranslatef(x, y, 0);
                drawFilledCircle({0, 0}, 0.012, 1000);
            } glPopMatrix();
        }
        else {
            glPointSize(5.0);
            drawPoint({x, y});
        }
    }
}

void drawSecondsHand(Point2D c, int t) {
    glPushMatrix(); {
        double a = 90.0 - (t % 60) * 6.0;

        glColor3ub(68, 214, 44);
        glRotatef(a, c.x, c.y, 1.0);
        glBegin(GL_POLYGON); {
            glVertex2d(-0.025, +0.0025);
            glVertex2d(-0.025, -0.0025);
            glVertex2d(+0.440, -0.0025);
            glVertex2d(+0.440, +0.0025);
        } glEnd();
    } glPopMatrix();
}

void drawMinutesHand(Point2D c, int t) {
    glPushMatrix(); {
        double a = 90.0 - (t % 3600) * 0.1;

        glColor3ub(219, 64, 117);
        glRotatef(a, c.x, c.y, 1.0);
        glBegin(GL_POLYGON); {
            glVertex2d(-0.025, +0.008);
            glVertex2d(-0.025, -0.008);
            glVertex2d(+0.350, -0.008);
            glVertex2d(+0.400, +0.000);
            glVertex2d(+0.350, +0.008);
        } glEnd();
    } glPopMatrix();
}

void drawHoursHand(Point2D c, int t) {
    glPushMatrix(); {
        double a = 90.0 - (t % 43200) * (1.0 / 120.0);

        glColor3ub(25, 230, 255);
        glRotatef(a, c.x, c.y, 1.0);
        glBegin(GL_POLYGON); {
            glVertex2d(-0.025, +0.015);
            glVertex2d(-0.025, -0.015);
            glVertex2d(+0.250, -0.015);
            glVertex2d(+0.340, +0.000);
            glVertex2d(+0.250, +0.015);
        } glEnd();
    } glPopMatrix();
}

void initGL() {
    glClearColor((double) 30 / 255, (double) 30 / 255, (double) 30 / 255, 1.0);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glPushMatrix(); {
        glTranslatef(0.0, 0.5, 0.0);

        glColor3ub(77, 77, 255);
        drawFilledCircle({0, 0}, 0.45, 10000); // clock

        drawClockMarks({0.0f, 0.0f}, 0.42f);

        int t = getCurrentTime();
        drawHoursHand({0, 0}, t);
        drawMinutesHand({0, 0}, t);
        drawSecondsHand({0, 0}, t);

        glColor3ub(0, 0, 0);
        drawFilledCircle({0, 0}, 0.015, 1000); // clock center

        glPushMatrix(); {
            glTranslated(0, -0.47, 0);

            int mt = getCurrentTimeMilli() % 2000;
            double mxa = 18.0;
            double a = mxa * cosf(2 * M_PI * mt / 2000.0 + degToRad(90.0));
            double r = 0.90;
            double x = r * sinf(degToRad(a));
            double y = r * cosf(degToRad(a));

            glColor3ub(255, 255, 255);
            glBegin(GL_LINES); { // pendulum wire
                glVertex2d(0, 0);
                glVertex2d(x, -y);
            } glEnd();

            glColor3ub(77, 77, 255);
            drawFilledCircle({0.0f, 0.0f}, 0.030f, 1000); // bod holder

            // glColor3ub(136, 139, 141);
            glColor3ub(255, 255, 255);
            drawFilledCircle({x, -y}, 0.1f, 1000); // bob
        } glPopMatrix();
    } glPopMatrix();

    glutSwapBuffers();
}

void idle() {
    glutPostRedisplay();
}

void reshape(GLsizei width, GLsizei height) {
    if (height == 0) {
        height = 1;
    }
    double aspect = (double) width / (double) height;

    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (width >= height) {
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(800, 800);
    glutCreateWindow("Clock");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutReshapeFunc(reshape);
    initGL();
    glutMainLoop();
    return 0;
}
