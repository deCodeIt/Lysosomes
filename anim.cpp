#include <GL/gl.h>
#include <GL/glu.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <vector>
#include "sdlglutils.h"

#define CELL_GROW_RATE 0.01
#define STEPS 5
#define CONTROL_POINTS 4
#define PI 3.14159265359
#define CELL_SECTORS 16
#define TCELL_CHANGE_STEP 30.0
#define CELL_STANDARD_DEVIATION 0.2
#define CELL_RADIUS 6.0
#define CAPSULE_WIDTH 1.0
#define CAPSULE_HEIGHT 2.0
#define CAPSULE_CROSSSECTION_RADIUS 0.5
#define TCAPSULE_CHANGE_STEP 50.0
#define CAPSULE_ROTATION_LAGRANGE 0.1
#define PRIMARY_LYSOSOME_SECTORS 10
#define PRIMARY_LYSOSOME_RADIUS 1.0
#define PRIMARY_LYSOSOME_STANDARD_DEVIATION 0.5
#define TPRIMARY_LYSOSOME_CHANGE_STEP 30.0

using namespace std;

struct _capsule{
   GLfloat x,y,z,theta,penumbra;
   vector<vector<GLfloat>> *capsuleCtrlPoints;
   vector <vector<GLfloat> > *capsulePathCtrlPoints;
   // theta is the angle w.r.t x axis
   // when theta == 0 ; capsule is standing upright
   // penumbra is the membrane of the capsule
} capsule;

struct _primaryLysosome{
   GLfloat x,y,z;
   vector<vector<GLfloat>> *lysosomeCtrlPoints,*pointsNow,*pointsLater;

} primaryLysosome;

struct _scale{
   GLfloat x,y,z;
} scale;

float tCell=0.0,tPrimaryLysosome=0.0;
float tCapsule=0.0;

vector< vector<GLfloat> > *ctrlpoints = new vector< vector <GLfloat> >(CELL_SECTORS+3, vector< GLfloat >(3)); // +4 for duplicate ending points for bsplines; // changes with tCell value as in interpolation
vector< vector<GLfloat> > *pointsNow = new vector< vector <GLfloat> >(CELL_SECTORS+3, vector< GLfloat >(3)); // cell points in current configuration
vector< vector<GLfloat> > *pointsLater = new vector< vector <GLfloat> >(CELL_SECTORS+3, vector< GLfloat >(3)); // cell points for next configuration

void getCurvePoints(vector< vector< GLfloat> > *,float, float, float);
float getInBetween(float, float);
void getInterpolatedPointsForCurve(vector< vector < GLfloat > > *,vector< vector < GLfloat > > *,vector< vector < GLfloat > > *);
void copyCurveConfiguration(vector< vector<GLfloat> > *,vector< vector<GLfloat> > *);
void drawCell(vector< vector<GLfloat> > *,int);
void getCapsulePoints(vector< vector< GLfloat> > *);
void drawCapsule(vector< vector<GLfloat> > *);
void getCapsulePathPoints(vector< vector< GLfloat> > *);
void drawPath(vector< vector<GLfloat> > *);
void drawPrimaryLysosome(vector< vector<GLfloat> > *);

GLuint loadTexture( const char *filename, int width, int height ) {
   GLuint texture;
   unsigned char * data;
   FILE * file;

   file = fopen( filename, "rb" );
   if ( file == NULL ) return 0;
   // printf("Reading File: %s\n",filename);
   data = (unsigned char *)malloc( width * height * 4 );
   fread( data, width * height * 4, 1, file );
   fclose( file );

   glGenTextures( 1, &texture ); //generate the texture with the loaded data
   glBindTexture( GL_TEXTURE_2D, texture ); //bind the texture to it’s array
   glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE ); //set texture environment parameters

   glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
   glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

   glTexImage2D(GL_TEXTURE_2D, 0, 4, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

   free( data ); //free the texture
   return texture; //return whether it was successfull
}

void init(void)
{
   srand(time(NULL));
   
   // for initial animation
   scale.x = 0.01;
   scale.y = scale.x;
   scale.z = 1.0;

   glClearColor(1.0, 1.0, 1.0, 0.0);
   glShadeModel(GL_FLAT);
   
   glEnable(GL_MAP1_VERTEX_3);
   getCurvePoints(ctrlpoints,CELL_RADIUS,CELL_STANDARD_DEVIATION,CELL_SECTORS); // initial cell configuration
   copyCurveConfiguration(pointsNow,ctrlpoints); // the destination configuration is our current configuration
   getCurvePoints(pointsLater,CELL_RADIUS,CELL_STANDARD_DEVIATION,CELL_SECTORS); // make a new destination configuration

   capsule.x = 0.0;
   capsule.y = 10.0;
   capsule.z = 0.5;
   capsule.theta = -180.0;
   capsule.capsuleCtrlPoints = new vector< vector <GLfloat> >(10+3, vector< GLfloat >(3));
   capsule.capsulePathCtrlPoints = new vector< vector <GLfloat> >(12, vector< GLfloat >(3));
   getCapsulePoints(capsule.capsuleCtrlPoints);
   getCapsulePathPoints(capsule.capsulePathCtrlPoints);


   primaryLysosome.x=1.0;
   primaryLysosome.y=0.0;
   primaryLysosome.z=0.25;
   primaryLysosome.lysosomeCtrlPoints = new vector< vector <GLfloat> >(PRIMARY_LYSOSOME_SECTORS + 3, vector< GLfloat >(3));
   primaryLysosome.pointsNow = new vector< vector <GLfloat> >(PRIMARY_LYSOSOME_SECTORS + 3, vector< GLfloat >(3));
   primaryLysosome.pointsLater = new vector< vector <GLfloat> >(PRIMARY_LYSOSOME_SECTORS + 3, vector< GLfloat >(3));
   getCurvePoints(primaryLysosome.lysosomeCtrlPoints,PRIMARY_LYSOSOME_RADIUS,PRIMARY_LYSOSOME_STANDARD_DEVIATION,PRIMARY_LYSOSOME_SECTORS); // initial primary lysosome configuration
   copyCurveConfiguration(primaryLysosome.pointsNow,primaryLysosome.lysosomeCtrlPoints); // the destination configuration is our current configuration
   getCurvePoints(primaryLysosome.pointsLater,PRIMARY_LYSOSOME_RADIUS,PRIMARY_LYSOSOME_STANDARD_DEVIATION,PRIMARY_LYSOSOME_SECTORS); // initial primary lysosome configuration

}

float getInBetween(float a, float b){
   // returns a floating points between a and b
   return (rand()/(1.0*RAND_MAX))*(b-a) + a;
}

void getInterpolatedPointsForCurve(vector< vector < GLfloat > > *pointsCurrent,vector< vector < GLfloat > > *pointsNow,vector< vector < GLfloat > > *pointsLater){
   // pointsNow ==> pointsLater
   // interpolated by tCell as pointsCurrent
   float minusTCell = 1-tCell;
   for(int i=0;i<pointsNow->size();i++){
      pointsCurrent->at(i).at(0) = pointsNow->at(i).at(0)*minusTCell + pointsLater->at(i).at(0)*tCell;
      pointsCurrent->at(i).at(1) = pointsNow->at(i).at(1)*minusTCell + pointsLater->at(i).at(1)*tCell;
      pointsCurrent->at(i).at(2) = pointsNow->at(i).at(2)*minusTCell + pointsLater->at(i).at(2)*tCell;
   }
}

void getCurvePoints(vector< vector< GLfloat> > *points,float r, float sd, float sectors){
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

void copyCurveConfiguration(vector< vector<GLfloat> > *dest,vector< vector<GLfloat> > *source){
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
   // glColor3f(250.0/256.0,232.0/256.0,192.0/256.0);
   glColor3f(247.0/256.0,212.0/256.0,158.0/256.0);
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

   // control points
   // glPointSize(7.0);
   // glColor3f(1.0, 0.0, 0.0);
   // glBegin(GL_POINTS);
   // for (int i = 0; i < len; i++) {
   //    glVertex3fv(&ctrlpts->at(i).at(0));
//    }
   // glEnd();
}

void drawPrimaryLysosome(vector< vector<GLfloat> > *ctrlpts){

   glPushMatrix();
   int len = ctrlpts->size();
   glTranslatef(primaryLysosome.x,primaryLysosome.y,primaryLysosome.z);
   glBegin(GL_POLYGON);
   
   // glColor3f(rand()/(1.0*RAND_MAX), rand()/(1.0*RAND_MAX), rand()/(1.0*RAND_MAX));
   // glColor3f(250.0/256.0,232.0/256.0,192.0/256.0);
   glColor3f(0.0,0.0,1.0);
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

   // control points
   // glPointSize(7.0);
   // glColor3f(1.0, 0.0, 0.0);
   // glBegin(GL_POINTS);
   // for (int i = 0; i < len; i++) {
   //    glVertex3fv(&ctrlpts->at(i).at(0));
//    }
   // glEnd();
   glPopMatrix();
}

void getCapsulePoints(vector< vector< GLfloat> > *points){

   // right side points
   points->at(0).at(0) = CAPSULE_WIDTH/2.0;
   points->at(0).at(1) = -CAPSULE_HEIGHT/2.0;
   points->at(0).at(2) = capsule.z;
   points->at(1).at(0) = CAPSULE_WIDTH/2.0;
   points->at(1).at(1) = 0.0;
   points->at(1).at(2) = capsule.z;
   points->at(2).at(0) = CAPSULE_WIDTH/2.0;
   points->at(2).at(1) = CAPSULE_HEIGHT/2.0;
   points->at(2).at(2) = capsule.z;

   // top points
   points->at(3).at(0) = CAPSULE_CROSSSECTION_RADIUS*cos(60.0*PI/180.0);
   points->at(3).at(1) = CAPSULE_HEIGHT/2.0 + CAPSULE_CROSSSECTION_RADIUS*sin(60.0*PI/180.0);
   points->at(3).at(2) = capsule.z;
   points->at(4).at(0) = CAPSULE_CROSSSECTION_RADIUS*cos(120.0*PI/180.0);
   points->at(4).at(1) = CAPSULE_HEIGHT/2.0 + CAPSULE_CROSSSECTION_RADIUS*sin(120.0*PI/180.0);
   points->at(4).at(2) = capsule.z;

   // left points
   points->at(5).at(0) = -CAPSULE_WIDTH/2.0;
   points->at(5).at(1) = CAPSULE_HEIGHT/2.0;
   points->at(5).at(2) = capsule.z;
   points->at(6).at(0) = -CAPSULE_WIDTH/2.0;
   points->at(6).at(1) = 0.0;
   points->at(6).at(2) = capsule.z;
   points->at(7).at(0) = -CAPSULE_WIDTH/2.0;
   points->at(7).at(1) = -CAPSULE_HEIGHT/2.0;
   points->at(7).at(2) = capsule.z;

   // bottom points
   points->at(8).at(0) = CAPSULE_CROSSSECTION_RADIUS*cos(240.0*PI/180.0);
   points->at(8).at(1) = -CAPSULE_HEIGHT/2.0 + CAPSULE_CROSSSECTION_RADIUS*sin(240.0*PI/180.0);
   points->at(8).at(2) = capsule.z;
   points->at(9).at(0) = CAPSULE_CROSSSECTION_RADIUS*cos(300.0*PI/180.0);
   points->at(9).at(1) = -CAPSULE_HEIGHT/2.0 + CAPSULE_CROSSSECTION_RADIUS*sin(300.0*PI/180.0);
   points->at(9).at(2) = capsule.z;

   // duplicate right side points
   points->at(10).at(0) = CAPSULE_WIDTH/2.0;
   points->at(10).at(1) = -CAPSULE_HEIGHT/2.0;
   points->at(10).at(2) = capsule.z;
   points->at(11).at(0) = CAPSULE_WIDTH/2.0;
   points->at(11).at(1) = 0.0;
   points->at(11).at(2) = capsule.z;
   points->at(12).at(0) = CAPSULE_WIDTH/2.0;
   points->at(12).at(1) = CAPSULE_HEIGHT/2.0;
   points->at(12).at(2) = capsule.z;

}

void getCapsulePathPoints(vector< vector< GLfloat> > *points){

   // straight down
   points->at(0).at(0) = -4.0;
   points->at(0).at(1) = 17.0;
   points->at(0).at(2) = capsule.z;
   points->at(1).at(0) = -3.7;
   points->at(1).at(1) = 9.0;
   points->at(1).at(2) = capsule.z;
   points->at(2).at(0) = -3.5;
   points->at(2).at(1) = 8.0;
   points->at(2).at(2) = capsule.z;
   points->at(3).at(0) = -3.0;
   points->at(3).at(1) = 7.0;
   points->at(3).at(2) = capsule.z;
   points->at(4).at(0) = -3.0;
   points->at(4).at(1) = 6.0;
   points->at(4).at(2) = capsule.z;
   points->at(5).at(0) = -2.0;
   points->at(5).at(1) = 5.0;
   points->at(5).at(2) = capsule.z;
   points->at(6).at(0) = -1.0;
   points->at(6).at(1) = 4.5;
   points->at(6).at(2) = capsule.z;

   // time to go slightly towards right
   points->at(7).at(0) = 0.0;
   points->at(7).at(1) = 3.2;
   points->at(7).at(2) = capsule.z;
   points->at(8).at(0) = 1.0;
   points->at(8).at(1) = 3.0;
   points->at(8).at(2) = capsule.z;
   points->at(9).at(0) = 1.5;
   points->at(9).at(1) = 3.0;
   points->at(9).at(2) = capsule.z;

   // time to stop
   points->at(10).at(0) = 1.5;
   points->at(10).at(1) = 3.0;
   points->at(10).at(2) = capsule.z;
   points->at(11).at(0) = 1.5;
   points->at(11).at(1) = 3.0;
   points->at(11).at(2) = capsule.z;

}

void drawCapsule(vector< vector<GLfloat> > *ctrlpts){
   glPushMatrix();

   glTranslatef(capsule.x,capsule.y,capsule.z);
   glRotatef(capsule.theta,0.0,0.0,1.0);
   int len = ctrlpts->size();
   glBegin(GL_POLYGON);
   // glColor3f(rand()/(1.0*RAND_MAX), rand()/(1.0*RAND_MAX), rand()/(1.0*RAND_MAX));
   // glColor3f(250.0/256.0,232.0/256.0,192.0/256.0);
   glColor3f(1.0,0.0,0.0);
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

   glBegin(GL_POLYGON);
   glColor3f(1.0,1.0,1.0);
   for(GLfloat i=0.0;i<2*PI;i+=PI/30.0){
      glVertex3f(CAPSULE_WIDTH*0.7*cos(i),(CAPSULE_HEIGHT*0.7+CAPSULE_CROSSSECTION_RADIUS*0.8)*sin(i),capsule.z);
   }
   glEnd();

   glPopMatrix();
}

void drawPath(vector< vector<GLfloat> > *ctrlpts){
   glPushMatrix();

   int len = ctrlpts->size();
   glBegin(GL_LINE_STRIP);
   // glColor3f(rand()/(1.0*RAND_MAX), rand()/(1.0*RAND_MAX), rand()/(1.0*RAND_MAX));
   // glColor3f(250.0/256.0,232.0/256.0,192.0/256.0);
   glColor3f(1.0,0.0,0.0);
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

   // control points
   glPointSize(7.0);
   glColor3f(0.0, 1.0, 0.0);
   glBegin(GL_POINTS);
   for (int i = 0; i < len; i++) {
      glVertex3fv(&ctrlpts->at(i).at(0));
   }
   glEnd();

   glPopMatrix();
}

void display(void)
{

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   // cout << ctrlpoints.size();
   glPushMatrix();
   glScalef(scale.x,scale.y,scale.z);
   drawCell(ctrlpoints,ctrlpoints->size());
   drawPrimaryLysosome(primaryLysosome.lysosomeCtrlPoints);
   if(scale.x == scale.y && scale.x == 1.0){
      // hide till the cell has expanded
      drawCapsule(capsule.capsuleCtrlPoints);
      // drawPath(capsule.capsulePathCtrlPoints);
   }
   glScalef(1.0/scale.x,1.0/scale.y,1.0/scale.z);
   glPopMatrix();

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
   glutSwapBuffers();
   
}

void timer(int id){
   // cout << "Timer called at t: "<< tCell << endl;
   
   if(tCell<1.0){
      tCell+=1.0/TCELL_CHANGE_STEP;
   } else {
      // new end points as the cell has reached the destination configuration
      copyCurveConfiguration(pointsNow,pointsLater); // the destination configuration is our current configuration
      getCurvePoints(pointsLater,CELL_RADIUS,CELL_STANDARD_DEVIATION,CELL_SECTORS); // make a new destination configuration
      tCell=1.0/TCELL_CHANGE_STEP;
   }
   getInterpolatedPointsForCurve(ctrlpoints,pointsNow,pointsLater);

   // for primary lysosome
   if(tPrimaryLysosome<1.0){
      tPrimaryLysosome+=1.0/TPRIMARY_LYSOSOME_CHANGE_STEP;
   } else {
      // new end points as the cell has reached the destination configuration
      copyCurveConfiguration(primaryLysosome.pointsNow,primaryLysosome.pointsLater); // the destination configuration is our current configuration
      getCurvePoints(primaryLysosome.pointsLater,PRIMARY_LYSOSOME_RADIUS,PRIMARY_LYSOSOME_STANDARD_DEVIATION,PRIMARY_LYSOSOME_SECTORS); // initial primary lysosome configuration
      tPrimaryLysosome=1.0/TPRIMARY_LYSOSOME_CHANGE_STEP;
   }
   getInterpolatedPointsForCurve(primaryLysosome.lysosomeCtrlPoints, primaryLysosome.pointsNow, primaryLysosome.pointsLater);

   
   if(scale.x == scale.y && abs(scale.x-1.0) <= CELL_GROW_RATE){
      scale.x = 1.0;
      scale.y = scale.x;
      // move capsule
      
      if(tCapsule < capsule.capsulePathCtrlPoints->size()-3){
         cout << tCapsule << " :: theta: " << capsule.theta << "(" << capsule.x <<"," << capsule.y <<")" << endl;
         int point = floor(tCapsule);

         GLfloat u,u_2,u_3,b[ CONTROL_POINTS ];
         u = tCapsule-floor(tCapsule);
         u_2 = u*u;
         u_3 = u_2*u;

         b[0] = (1 - u_3 + 3*u_2 - 3*u)/6.0;
         b[1] = (4 - 6*u_2 + 3*u_3)/6.0;
         b[2] = (1 + 3*u + 3*u_2 - 3*u_3 )/6.0;
         b[3] = u_3/6.0;
         GLfloat old_x=capsule.x,old_y=capsule.y,old_z=capsule.z;
         capsule.x = b[0]*capsule.capsulePathCtrlPoints->at(point).at(0) + b[1]*capsule.capsulePathCtrlPoints->at(point+1).at(0) + b[2]*capsule.capsulePathCtrlPoints->at(point+2).at(0) + b[3]*capsule.capsulePathCtrlPoints->at(point+3).at(0);
         capsule.y = b[0]*capsule.capsulePathCtrlPoints->at(point).at(1) + b[1]*capsule.capsulePathCtrlPoints->at(point+1).at(1) + b[2]*capsule.capsulePathCtrlPoints->at(point+2).at(1) + b[3]*capsule.capsulePathCtrlPoints->at(point+3).at(1);
         
         // smooth rotation
         capsule.theta = CAPSULE_ROTATION_LAGRANGE * (tanh((capsule.y - old_y)/(capsule.x - old_x))*180.0/PI - 90.0) + (1- CAPSULE_ROTATION_LAGRANGE )*capsule.theta;
         
         tCapsule+=1.0/TCAPSULE_CHANGE_STEP;
         
      }
   } else {
      scale.x += CELL_GROW_RATE;
      scale.y = scale.x; 
   }
   glutPostRedisplay();
   glutTimerFunc(100,timer,5);
}

void keyboard(unsigned char key, int x, int y){
   switch(key) {
      case 's':
         // starts the motion
         timer(0);
         break;
      default:
         break;
   }
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
   glutKeyboardFunc(keyboard);
   // glutTimerFunc(100,timer,5);
   glutMainLoop();
   return 0;
}