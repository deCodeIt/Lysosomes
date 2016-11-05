#include <GL/gl.h>
#include <GL/glu.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <vector>
#include "sdlglutils.h"

#define STEPS 10
#define CONTROL_POINTS 4
#define PI 3.14159265359
#define CELL_SECTORS 20
#define TCELL_CHANGE_STEP 30.0
#define CELL_STANDARD_DEVIATION 0.2
#define CELL_RADIUS 6.0

using namespace std;

float tCell=0.0;

vector< vector<GLfloat> > *ctrlpoints = new vector< vector <GLfloat> >(CELL_SECTORS+3, vector< GLfloat >(3)); // +4 for duplicate ending points for bsplines; // changes with tCell value as in interpolation
vector< vector<GLfloat> > *pointsNow = new vector< vector <GLfloat> >(CELL_SECTORS+3, vector< GLfloat >(3)); // cell points in current configuration
vector< vector<GLfloat> > *pointsLater = new vector< vector <GLfloat> >(CELL_SECTORS+3, vector< GLfloat >(3)); // cell points for next configuration

void getCellPoints(vector< vector< GLfloat> > *,float, float, float);
float getInBetween(float, float);
void getInterpolatedPointsForCell(vector< vector < GLfloat > > *);
void copyCellConfiguration(vector< vector<GLfloat> > *,vector< vector<GLfloat> > *);
void drawCell(vector< vector<GLfloat> > *,int);


void init(void)
{
   srand(time(NULL));
   glClearColor(0.0, 0.0, 0.0, 0.0);
   glShadeModel(GL_FLAT);
   
   glEnable(GL_MAP1_VERTEX_3);
   getCellPoints(ctrlpoints,CELL_RADIUS,CELL_STANDARD_DEVIATION,CELL_SECTORS); // initial cell configuration
   copyCellConfiguration(pointsNow,ctrlpoints); // the destination configuration is our current configuration
   getCellPoints(pointsLater,CELL_RADIUS,CELL_STANDARD_DEVIATION,CELL_SECTORS); // make a new destination configuration
}

float getInBetween(float a, float b){
   // returns a floating points between a and b
   return (rand()/(1.0*RAND_MAX))*(b-a) + a;
}

void getInterpolatedPointsForCell(vector< vector < GLfloat > > *pointsCurrent){
   // pointsNow ==> pointsLater
   // interpolated by tCell as pointsCurrent
   float minusTCell = 1-tCell;
   for(int i=0;i<pointsNow->size();i++){
      pointsCurrent->at(i).at(0) = pointsNow->at(i).at(0)*minusTCell + pointsLater->at(i).at(0)*tCell;
      pointsCurrent->at(i).at(1) = pointsNow->at(i).at(1)*minusTCell + pointsLater->at(i).at(1)*tCell;
      pointsCurrent->at(i).at(2) = pointsNow->at(i).at(2)*minusTCell + pointsLater->at(i).at(2)*tCell;
   }
}

void getCellPoints(vector< vector< GLfloat> > *points,float r, float sd, float sectors){
   // generates cell boundary points with radius r and some standard deviation sd and with the circle divided in sectors
   float rdash,theta,sectorAngle = 2*PI/sectors; // the variance limit across the circle/curve boundary
   // float points[sectors+4][3]; // three points repeating at beginning and 3 points at th end with bsplines
   for(int i=0;i<sectors;i++){
      // generate those points
      rdash = r + r*sd*(rand()/(1.0*RAND_MAX) - 0.5);
      theta = getInBetween(sectorAngle*i,sectorAngle*(i+1));
      points->at(i).at(0) = rdash*cos(theta); // x
      points->at(i).at(1) = rdash*sin(theta); // y
      points->at(i).at(2) = 0.0; // z
   }

   // for smooth loop at intersection of start and end points
   points->at(sectors).at(0) = points->at(0).at(0);
   points->at(sectors).at(1) = points->at(0).at(1);
   points->at(sectors).at(2) = points->at(0).at(2); // n-2th points over
   points->at(sectors+1).at(0) = points->at(1).at(0);
   points->at(sectors+1).at(1) = points->at(1).at(1);
   points->at(sectors+1).at(2) = points->at(1).at(2); // n-1th points over
   points->at(sectors+2).at(0) = points->at(2).at(0);
   points->at(sectors+2).at(1) = points->at(2).at(1);
   points->at(sectors+2).at(2) = points->at(2).at(2); // end point over
}

void copyCellConfiguration(vector< vector<GLfloat> > *dest,vector< vector<GLfloat> > *source){
   // copies the source control points to dest as its control points
   for(int i=0;i<source->size();i++){
      dest->at(i).at(0) = source->at(i).at(0);
      dest->at(i).at(1) = source->at(i).at(1);
      dest->at(i).at(2) = source->at(i).at(2);
   }
}

void drawCell(vector< vector<GLfloat> > *ctrlpts,int len){
   glBegin(GL_POLYGON);
   // glColor3f(rand()/(1.0*RAND_MAX), rand()/(1.0*RAND_MAX), rand()/(1.0*RAND_MAX));
   glColor3f(250.0/256.0,232.0/256.0,192.0/256.0);
      GLfloat u,u_2,u_3,b[ CONTROL_POINTS ]; // 'b' stores the blending polynomials for bsplines
      // int len = sizeof(ctrlpts)/sizeof(ctrlpts[0]);
      // cout << len ;
      for(int point = 0; point<len-3; point++){
         GLfloat curveVertex[3];
      
         for(float iter=0.0;iter<STEPS;iter+=1.0){
            u = iter/STEPS;
            u_2 = u*u;
            u_3 = u_2*u;

            b[0] = (1 - u_3 + 3*u_2 - 3*u)/6.0;
            b[1] = (4 - 6*u_2 + 3*u_3)/6.0;
            b[2] = (1 + 3*u + 3*u_2 - 3*u_3 )/6.0;
            b[3] = u_3/6.0;

            for(int k=0;k<3;k++){
               curveVertex[k] = b[0]*ctrlpts->at(point).at(k) + b[1]*ctrlpts->at(point+1).at(k) + b[2]*ctrlpts->at(point+2).at(k) + b[3]*ctrlpts->at(point+3).at(k);
            }

            glVertex3fv(curveVertex);

         }
      }
   glEnd();

   /* The following code displays the control points as dots. */
   // glPointSize(5.0);
   // glColor3f(1.0, 1.0, 0.0);
   // glBegin(GL_POINTS);
   //    for (int i = 0; i < len; i++) 
   //       glVertex3fv(&ctrlpts->at(i).at(0));
   // glEnd();
   glutSwapBuffers();
}

void display(void)
{

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   // cout << ctrlpoints.size();
   
   drawCell(ctrlpoints,ctrlpoints->size());

   // glColor3f(1.0, 0.0, 0.0);
   // glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 3, &ctrlpoints[0][0]);
   // glMapGrid1f(STEPS,0.0,1.0);
   // glEvalMesh1(GL_LINE,0,STEPS);
   // glColor3f(0.0, 1.0, 0.0);
   // glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 3, &ctrlpoints[2][0]);
   // glMapGrid1f(STEPS,0.0,1.0);
   // glEvalMesh1(GL_LINE,0,STEPS);
   // glColor3f(0.0, 0.0, 1.0);
   // glMap1f(GL_MAP1_VERTEX_3, 0.0, 1.0, 3, 3, &ctrlpoints[4][0]);
   // glMapGrid1f(STEPS,0.0,1.0);
   // glEvalMesh1(GL_LINE,0,STEPS);


   // glBegin(GL_LINE_STRIP);
   //    for (i = 0; i <= STEPS; i++) 
   //       glEvalCoord1f((GLfloat) i/STEPS);
   // glEnd();
   
}

void timer(int id){
   // cout << "Timer called at t: "<< tCell << endl;
   if(tCell<1.0){
      tCell+=1.0/TCELL_CHANGE_STEP;
   } else {
      // new end points as the cell has reached the destination configuration
      copyCellConfiguration(pointsNow,pointsLater); // the destination configuration is our current configuration
      getCellPoints(pointsLater,CELL_RADIUS,CELL_STANDARD_DEVIATION,CELL_SECTORS); // make a new destination configuration
      tCell=1.0/TCELL_CHANGE_STEP;
   }
   getInterpolatedPointsForCell(ctrlpoints);
   glutPostRedisplay();
   glutTimerFunc(100,timer,5);
}

void reshape(int w, int h)
{
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if (w <= h)
      glOrtho(-10.0, 10.0, -10.0*(GLfloat)h/(GLfloat)w, 
               10.0*(GLfloat)h/(GLfloat)w, -5.0, 5.0);
   else
      glOrtho(-5.0*(GLfloat)w/(GLfloat)h, 
               5.0*(GLfloat)w/(GLfloat)h, -5.0, 5.0, -5.0, 5.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glEnable(GL_TEXTURE_2D);
   glEnable(GL_DEPTH_TEST);
}

int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
   glutInitWindowSize (500, 500);
   glutInitWindowPosition (100, 100);
   glutCreateWindow (argv[0]);
   init ();
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutTimerFunc(100,timer,5);
   glutMainLoop();
   return 0;
}