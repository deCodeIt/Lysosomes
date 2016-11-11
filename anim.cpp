#include <GL/gl.h>
#include <GL/glu.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <limits>

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
#define TCAPSULE_CHANGE_STEP 30.0
#define CAPSULE_ROTATION_LAGRANGE 0.1

#define PRIMARY_LYSOSOME_SECTORS 10
#define PRIMARY_LYSOSOME_RADIUS 1.0
#define PRIMARY_LYSOSOME_STANDARD_DEVIATION 0.5
#define TPRIMARY_LYSOSOME_CHANGE_STEP 30.0
#define TPRIMARY_LYSOSOME_MOVE_STEP 60.0
#define TPRIMARY_LYSOSOME_SIZE_STEP 0.1
#define TPRIMARY_LYSOSOME_TEXTURE_CHANGE_STEP 2.0
#define TPRIMARY_LYSOSOME_MAGNIFICATION 1.0

#define NUCLEUS_SECTORS 15
#define NUCLEUS_RADIUS 1.5
#define NUCLEUS_STANDARD_DEVIATION 0.2
#define TNUCLEUS_CHANGE_STEP 30.0

#define NUCLEOLUS_SECTORS 10
#define NUCLEOLUS_RADIUS 0.7
#define NUCLEOLUS_STANDARD_DEVIATION 0.1
#define TNUCLEOLUS_CHANGE_STEP 30.0

#define MITOCHONDRIA_POINTS 30.0
#define MITOCHONDRIA_RADIUS_X 2.0
#define MITOCHONDRIA_RADIUS_Y 0.6
#define TMITOCHONDRIA_CHANGE_STEP 30.0
#define MITOCHONDRIA_STANDARD_DEVIATION_X 0.0
#define MITOCHONDRIA_STANDARD_DEVIATION_Y 0.4

#define GOLGI_DIMENSION 3.0

#define CAPSULE_DISAPPEAR 2.5

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
   GLfloat x,y,z,r,g,b,xdash,ydash,textureXdash,textureYdash;
   vector<vector<GLfloat>> *lysosomeCtrlPoints,*pointsNow,*pointsLater;
   // parameters for a primary lysosome for its display and animation

};

struct _mitochondria{
   GLfloat x,y,z,theta,xdash,ydash;
   // theta w.r.t x-axis
   vector<vector<GLfloat>> *ctrlpoints,*pointsNow,*pointsLater;
};

vector <_mitochondria> mitochondria(2); // number of mitochondria

vector <_primaryLysosome> primaryLysosome(4); // number of lysosomes + 1 nucleus + 1 nucleolus

GLfloat textureMitochondria,texturePrimaryLysosome,textureLysosome,textureGolgi;

bool flag = false; // to start the animation
bool pause = true; // pause tha animation

struct _scale{
   GLfloat x,y,z;
   // gets the user from normal magnification to cell level magnification
} scale;

float tCell=0.0,tPrimaryLysosome=0.0;
float tCapsule=0.0,tLysosomeMove=0.0,tMitochondria=0.0,tGolgi=0.0; // parameters for movement of cell and its constituents
int tLysosomeChange=0;

vector< vector<GLfloat> > *ctrlpoints = new vector< vector <GLfloat> >(CELL_SECTORS+3, vector< GLfloat >(3)); // +4 for duplicate ending points for bsplines; // changes with tCell value as in interpolation
vector< vector<GLfloat> > *pointsNow = new vector< vector <GLfloat> >(CELL_SECTORS+3, vector< GLfloat >(3)); // cell points in current configuration of the curve
vector< vector<GLfloat> > *pointsLater = new vector< vector <GLfloat> >(CELL_SECTORS+3, vector< GLfloat >(3)); // cell points for next configuration of the curve

void getCurvePoints(vector< vector< GLfloat> > *,float, float, float);
void getMitochondriaPoints(vector< vector< GLfloat> > *,float, float, float, float);
float getInBetween(float, float);
void getInterpolatedPointsForCurve(vector< vector < GLfloat > > *,vector< vector < GLfloat > > *,vector< vector < GLfloat > > *);
void copyCurveConfiguration(vector< vector<GLfloat> > *,vector< vector<GLfloat> > *);
void drawCell(vector< vector<GLfloat> > *,int);
void getCapsulePoints(vector< vector< GLfloat> > *);
void drawCapsule(vector< vector<GLfloat> > *);
void getCapsulePathPoints(vector< vector< GLfloat> > *);
void drawPath(vector< vector<GLfloat> > *);
void drawPrimaryLysosome(vector< vector<GLfloat> > *, int);
void drawMitochondria(int);

GLuint loadTexture( const char *filename, int width, int height ) {
   // loads the given texture from the file with specified width and height
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
   srand(time(NULL)); // seef for random number generation
   
   // for initial animation
   scale.x = 0.01;
   scale.y = scale.x; // scaling x and y directions equally
   scale.z = 1.0;

   glClearColor(1.0, 1.0, 1.0, 0.0);
   glShadeModel(GL_SMOOTH);
   
   // glEnable(GL_MAP1_VERTEX_3);
   getCurvePoints(ctrlpoints,CELL_RADIUS,CELL_STANDARD_DEVIATION,CELL_SECTORS); // initial cell configuration
   copyCurveConfiguration(pointsNow,ctrlpoints); // the destination configuration is our current configuration
   getCurvePoints(pointsLater,CELL_RADIUS,CELL_STANDARD_DEVIATION,CELL_SECTORS); // make a new destination configuration

   // capsule is our food to be eaten
   capsule.x = 0.0;
   capsule.y = 10.0;
   capsule.z = 0.5;
   capsule.theta = -180.0;
   capsule.capsuleCtrlPoints = new vector< vector <GLfloat> >(10+3, vector< GLfloat >(3));
   capsule.capsulePathCtrlPoints = new vector< vector <GLfloat> >(12, vector< GLfloat >(3));
   getCapsulePoints(capsule.capsuleCtrlPoints);
   getCapsulePathPoints(capsule.capsulePathCtrlPoints);


   // Nucleus
   primaryLysosome[0].x=0.0;
   primaryLysosome[0].y=0.0;
   primaryLysosome[0].z=0.25;

   // nucleolus
   primaryLysosome[1].x=0.0;
   primaryLysosome[1].y=0.0;
   primaryLysosome[1].z=0.26;

   // first lysosome
   primaryLysosome[2].x=3.5;
   primaryLysosome[2].y=0.2;
   primaryLysosome[2].z=0.90;

   //second lysosome
   primaryLysosome[3].x=0.75;
   primaryLysosome[3].y=-3.0;
   primaryLysosome[3].z=0.25;

   //third lysosome
   // primaryLysosome[4].x=-2.0;
   // primaryLysosome[4].y=-2.0;
   // primaryLysosome[4].z=0.25;

   for(int i=2;i<primaryLysosome.size();i++){
      primaryLysosome.at(i).textureXdash = primaryLysosome.at(i).x;
      primaryLysosome.at(i).textureYdash = primaryLysosome.at(i).y;
      primaryLysosome.at(i).r = 0.0;
      primaryLysosome.at(i).g = 0.5;
      primaryLysosome.at(i).b = 1.0;
      primaryLysosome.at(i).lysosomeCtrlPoints = new vector< vector <GLfloat> >(PRIMARY_LYSOSOME_SECTORS + 3, vector< GLfloat >(3));
      primaryLysosome.at(i).pointsNow = new vector< vector <GLfloat> >(PRIMARY_LYSOSOME_SECTORS + 3, vector< GLfloat >(3));
      primaryLysosome.at(i).pointsLater = new vector< vector <GLfloat> >(PRIMARY_LYSOSOME_SECTORS + 3, vector< GLfloat >(3));
      getCurvePoints(primaryLysosome.at(i).lysosomeCtrlPoints,PRIMARY_LYSOSOME_RADIUS,PRIMARY_LYSOSOME_STANDARD_DEVIATION,PRIMARY_LYSOSOME_SECTORS); // initial primary lysosome configuration
      copyCurveConfiguration(primaryLysosome.at(i).pointsNow,primaryLysosome.at(i).lysosomeCtrlPoints); // the destination configuration is our current configuration
      getCurvePoints(primaryLysosome.at(i).pointsLater,PRIMARY_LYSOSOME_RADIUS,PRIMARY_LYSOSOME_STANDARD_DEVIATION,PRIMARY_LYSOSOME_SECTORS); // initial primary lysosome configuration
   }
   // nucleus
   primaryLysosome.at(0).r = 0.0;
   primaryLysosome.at(0).g = 0.0;
   primaryLysosome.at(0).b = 0.0;
   primaryLysosome.at(0).lysosomeCtrlPoints = new vector< vector <GLfloat> >(NUCLEUS_SECTORS + 3, vector< GLfloat >(3));
   primaryLysosome.at(0).pointsNow = new vector< vector <GLfloat> >(NUCLEUS_SECTORS + 3, vector< GLfloat >(3));
   primaryLysosome.at(0).pointsLater = new vector< vector <GLfloat> >(NUCLEUS_SECTORS + 3, vector< GLfloat >(3));
   getCurvePoints(primaryLysosome.at(0).lysosomeCtrlPoints,NUCLEUS_RADIUS,NUCLEUS_STANDARD_DEVIATION,NUCLEUS_SECTORS); // initial primary lysosome configuration
   copyCurveConfiguration(primaryLysosome.at(0).pointsNow,primaryLysosome.at(0).lysosomeCtrlPoints); // the destination configuration is our current configuration
   getCurvePoints(primaryLysosome.at(0).pointsLater,NUCLEUS_RADIUS,NUCLEUS_STANDARD_DEVIATION,NUCLEUS_SECTORS); // initial primary lysosome configuration

   // nucleolus
   primaryLysosome.at(1).r = 0.2;
   primaryLysosome.at(1).g = 0.2;
   primaryLysosome.at(1).b = 0.2;
   primaryLysosome.at(1).lysosomeCtrlPoints = new vector< vector <GLfloat> >(NUCLEOLUS_SECTORS + 3, vector< GLfloat >(3));
   primaryLysosome.at(1).pointsNow = new vector< vector <GLfloat> >(NUCLEOLUS_SECTORS + 3, vector< GLfloat >(3));
   primaryLysosome.at(1).pointsLater = new vector< vector <GLfloat> >(NUCLEOLUS_SECTORS + 3, vector< GLfloat >(3));
   getCurvePoints(primaryLysosome.at(1).lysosomeCtrlPoints,NUCLEOLUS_RADIUS,NUCLEOLUS_STANDARD_DEVIATION,NUCLEOLUS_SECTORS); // initial primary lysosome configuration
   copyCurveConfiguration(primaryLysosome.at(1).pointsNow,primaryLysosome.at(1).lysosomeCtrlPoints); // the destination configuration is our current configuration
   getCurvePoints(primaryLysosome.at(1).pointsLater,NUCLEOLUS_RADIUS,NUCLEOLUS_STANDARD_DEVIATION,NUCLEOLUS_SECTORS); // initial primary lysosome configuration

   // mitochondria
   
   mitochondria.at(0).x = -2.5;
   mitochondria.at(0).y = 1.5;
   mitochondria.at(0).z = 0.25;
   mitochondria.at(0).theta = 60;

   mitochondria.at(1).x = 2.5;
   mitochondria.at(1).y = -1.5;
   mitochondria.at(1).z = 0.25;
   mitochondria.at(1).theta = 135;

   for(int i=0;i<mitochondria.size();i++){
      mitochondria.at(i).xdash = mitochondria.at(i).x;
      mitochondria.at(i).ydash = mitochondria.at(i).y;
      mitochondria.at(i).ctrlpoints = new vector<vector<GLfloat>>(MITOCHONDRIA_POINTS+3,vector<GLfloat>(3));
      mitochondria.at(i).pointsNow = new vector<vector<GLfloat>>(MITOCHONDRIA_POINTS+3,vector<GLfloat>(3));
      mitochondria.at(i).pointsLater = new vector<vector<GLfloat>>(MITOCHONDRIA_POINTS+3,vector<GLfloat>(3));
   }

   getMitochondriaPoints(mitochondria.at(0).ctrlpoints,MITOCHONDRIA_RADIUS_X - getInBetween(0,MITOCHONDRIA_RADIUS_X*0.1),MITOCHONDRIA_RADIUS_Y - getInBetween(0,MITOCHONDRIA_RADIUS_Y*0.2),MITOCHONDRIA_STANDARD_DEVIATION_X,MITOCHONDRIA_STANDARD_DEVIATION_Y);
   getMitochondriaPoints(mitochondria.at(1).ctrlpoints,MITOCHONDRIA_RADIUS_X - getInBetween(MITOCHONDRIA_RADIUS_X*0.3,MITOCHONDRIA_RADIUS_X*0.5),MITOCHONDRIA_RADIUS_Y - getInBetween(MITOCHONDRIA_RADIUS_Y*0.3,MITOCHONDRIA_RADIUS_Y*0.5),MITOCHONDRIA_STANDARD_DEVIATION_X,MITOCHONDRIA_STANDARD_DEVIATION_Y);
   copyCurveConfiguration(mitochondria.at(0).pointsNow, mitochondria.at(0).ctrlpoints);
   copyCurveConfiguration(mitochondria.at(1).pointsNow, mitochondria.at(1).ctrlpoints);
   getMitochondriaPoints(mitochondria.at(0).pointsLater,MITOCHONDRIA_RADIUS_X - getInBetween(0,MITOCHONDRIA_RADIUS_X*0.1),MITOCHONDRIA_RADIUS_Y - getInBetween(0,MITOCHONDRIA_RADIUS_Y*0.2),MITOCHONDRIA_STANDARD_DEVIATION_X,MITOCHONDRIA_STANDARD_DEVIATION_Y);
   getMitochondriaPoints(mitochondria.at(1).pointsLater,MITOCHONDRIA_RADIUS_X - getInBetween(MITOCHONDRIA_RADIUS_X*0.3,MITOCHONDRIA_RADIUS_X*0.5),MITOCHONDRIA_RADIUS_Y - getInBetween(MITOCHONDRIA_RADIUS_Y*0.3,MITOCHONDRIA_RADIUS_Y*0.5),MITOCHONDRIA_STANDARD_DEVIATION_X,MITOCHONDRIA_STANDARD_DEVIATION_Y);

   // load the texture
   textureMitochondria = loadTexture("mitochondriaTexNormalized.raw",95,37);
   texturePrimaryLysosome = loadTexture("primaryLysosomeNormalized.raw",33,34);
   textureLysosome = loadTexture("lysosomeNormalized.raw",181,180);
   textureGolgi = loadTexture("golgiNormalized.raw",170,140);
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

void getMitochondriaPoints(vector< vector< GLfloat> > *points,float radiusX, float radiusY, float sdx, float sdy){
   // generates mitochondria boundary points
   
   float sectors = MITOCHONDRIA_POINTS;
   float rdashX,rdashY,theta,sectorAngle = 2*PI/sectors; // the variance limit across the circle/curve boundary
   // three points repeating at beginning and 3 points at the end with bsplines
   for(int i=0;i<sectors;i++){
      // generate those points
      rdashX = radiusX*(1.0 + sdx*(rand()/(1.0*RAND_MAX) - 0.5));
      rdashY = radiusY*(1.0 + sdy*(rand()/(1.0*RAND_MAX) - 0.5));
      theta = getInBetween(sectorAngle*i,sectorAngle*(i+1));
      points->at(i).at(0) = rdashX*cos(theta); // x
      points->at(i).at(1) = rdashY*sin(theta); // y
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

   glColor3f(247.0/256.0,212.0/256.0,158.0/256.0);
      GLfloat u,u_2,u_3,b[ CONTROL_POINTS ]; // 'b' stores the blending polynomials for bsplines
      
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
            // interpolate bsplines
            for(int k=0;k<3;k++){
               curveVertex[k] = b[0]*ctrlpts->at(point).at(k) + b[1]*ctrlpts->at(point+1).at(k) + b[2]*ctrlpts->at(point+2).at(k) + b[3]*ctrlpts->at(point+3).at(k);
            }

            glVertex3fv(curveVertex);

         }
      }
   glEnd();
}

void drawPrimaryLysosome(vector< vector<GLfloat> > *ctrlpts,int i){

   glPushMatrix();
   int len = ctrlpts->size();
   glTranslatef(primaryLysosome.at(i).x,primaryLysosome.at(i).y,primaryLysosome.at(i).z);
   if(tLysosomeChange>=TPRIMARY_LYSOSOME_TEXTURE_CHANGE_STEP){
      // shift txture to show movement
      primaryLysosome.at(i).textureXdash = PRIMARY_LYSOSOME_RADIUS*(1+PRIMARY_LYSOSOME_STANDARD_DEVIATION)-getInBetween(PRIMARY_LYSOSOME_RADIUS*(1+PRIMARY_LYSOSOME_STANDARD_DEVIATION)*0.1,PRIMARY_LYSOSOME_RADIUS*(1+PRIMARY_LYSOSOME_STANDARD_DEVIATION)*0.2);
      primaryLysosome.at(i).textureYdash = PRIMARY_LYSOSOME_RADIUS*(1+PRIMARY_LYSOSOME_STANDARD_DEVIATION)-getInBetween(PRIMARY_LYSOSOME_RADIUS*(1+PRIMARY_LYSOSOME_STANDARD_DEVIATION)*0.1,PRIMARY_LYSOSOME_RADIUS*(1+PRIMARY_LYSOSOME_STANDARD_DEVIATION)*0.2);
   }
   if(i>=2){
      glEnable(GL_TEXTURE_2D);
      if(i==2 && tLysosomeMove>=CAPSULE_DISAPPEAR/2.0){
         glBindTexture(GL_TEXTURE_2D,textureLysosome);
      } else {
         glBindTexture(GL_TEXTURE_2D,texturePrimaryLysosome);
      }
      glColor3f(1.0,1.0,1.0);
   } else {
      glDisable(GL_TEXTURE_2D);
      glColor3f(primaryLysosome.at(i).r,primaryLysosome.at(i).g,primaryLysosome.at(i).b);
   }
   glBegin(GL_POLYGON);
   
      GLfloat u,u_2,u_3,b[ CONTROL_POINTS ]; // 'b' stores the blending polynomials for bsplines

      GLfloat minX = numeric_limits<float>::max(),maxX = numeric_limits<float>::min(),minY = numeric_limits<float>::max(),maxY = numeric_limits<float>::min();
      vector<vector<GLfloat> > curvePoints;
      for(int point = 0; point<len-3; point++){
         // interpolate bsplines
         for(float iter=0.0;iter<STEPS;iter+=1.0){
            vector<GLfloat> curveVertex(3);
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
            curvePoints.push_back(curveVertex);
            if(curveVertex[0] < minX){
               minX = curveVertex[0];
            }
            if(curveVertex[0] > maxX){
               maxX = curveVertex[0];
            }
            if(curveVertex[1] < minY){
               minY = curveVertex[1];
            }
            if(curveVertex[1] > maxY){
               maxY = curveVertex[1];
            }
         }
      }
      for(int k=0;k<curvePoints.size();k++){
         // cout << (curvePoints.at(k).at(0) - minX)/(maxX-minX)<< "," << (curvePoints.at(k).at(1) - minY)/(maxY-minY) << endl;
         glTexCoord2f((curvePoints.at(k).at(0) - minX)/(maxX-minX),(curvePoints.at(k).at(1) - minY)/(maxY-minY));
         glVertex3f(curvePoints.at(k).at(0),curvePoints.at(k).at(1),curvePoints.at(k).at(2));
      }
   glEnd();

   glDisable(GL_TEXTURE_2D);
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

   if(tLysosomeMove>=1.0){
         // time to reduce the food size as the lysosome is now consuming it
         tLysosomeMove+=TPRIMARY_LYSOSOME_SIZE_STEP;
         glScalef(1.0/tLysosomeMove,1.0/tLysosomeMove,1.0);
      }

   int len = ctrlpts->size();
   glBegin(GL_POLYGON);
   glColor3f(1.0,0.0,0.0);
      GLfloat u,u_2,u_3,b[ CONTROL_POINTS ]; // 'b' stores the blending polynomials for bsplines
      for(int point = 0; point<len-3; point++){
         GLfloat curveVertex[3];
      // mitochondria.a
         for(float iter=0.0;iter<STEPS;iter+=1.0){
            u = iter/STEPS;
            u_2 = u*u;
            u_3 = u_2*u;

            b[0] = (1 - u_3 + 3*u_2 - 3*u)/6.0;
            b[1] = (4 - 6*u_2 + 3*u_3)/6.0;
            b[2] = (1 + 3*u + 3*u_2 - 3*u_3 )/6.0;
            b[3] = u_3/6.0;
            // interpolate b-splines
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
      // if(i<=1.5*PI)
      //    glColor3f(0.0,1.0,0.0);
      // else
      //    glColor3f(0.0,0.0,1.0);
      glVertex3f(CAPSULE_WIDTH*0.7*cos(i),(CAPSULE_HEIGHT*0.7+CAPSULE_CROSSSECTION_RADIUS*0.8)*sin(i),capsule.z);
   }
   glEnd();

   if(tLysosomeMove >= 1.0){
      glScalef(tLysosomeMove,tLysosomeMove,1.0);
   }

   glPopMatrix();
}

void drawMitochondria(int i){
   glPushMatrix();

   glTranslatef(mitochondria.at(i).xdash,mitochondria.at(i).ydash,mitochondria.at(i).z);
   glRotatef(mitochondria.at(i).theta,0.0,0.0,1.0);
   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D,textureMitochondria);
   glBegin(GL_POLYGON);
   glColor3f(1.0,1.0,1.0);

   // compute normalized texture
   GLfloat minX = numeric_limits<float>::max(),maxX = numeric_limits<float>::min(),minY = numeric_limits<float>::max(),maxY = numeric_limits<float>::min();

   for(GLfloat j=0.0;j<mitochondria.at(i).ctrlpoints->size();j++){
      // for texture
      if(mitochondria.at(i).ctrlpoints->at(j).at(0) < minX){
         minX = mitochondria.at(i).ctrlpoints->at(j).at(0);
      }
      if(mitochondria.at(i).ctrlpoints->at(j).at(0) > maxX){
         maxX = mitochondria.at(i).ctrlpoints->at(j).at(0);
      }
      if(mitochondria.at(i).ctrlpoints->at(j).at(1) < minY){
         minY = mitochondria.at(i).ctrlpoints->at(j).at(1);
      }
      if(mitochondria.at(i).ctrlpoints->at(j).at(1) > maxY){
         maxY = mitochondria.at(i).ctrlpoints->at(j).at(1);
      }
   }

   for(GLfloat j=0.0;j<mitochondria.at(i).ctrlpoints->size();j++){
       
      glTexCoord2f((mitochondria.at(i).ctrlpoints->at(j).at(0) - minX)/(maxX-minX),(mitochondria.at(i).ctrlpoints->at(j).at(1) - minX)/(maxX-minX));
      glVertex3f(mitochondria.at(i).ctrlpoints->at(j).at(0),mitochondria.at(i).ctrlpoints->at(j).at(1),mitochondria.at(i).ctrlpoints->at(j).at(2));
   }
   glEnd();
   glDisable(GL_TEXTURE_2D);
   glPopMatrix();
}

void drawPath(vector< vector<GLfloat> > *ctrlpts){
   // additional function that draws out the path taken by the food/capsule
   // for debugging Purposes Only

   glPushMatrix();

   int len = ctrlpts->size();
   glBegin(GL_LINE_STRIP);
   glColor3f(1.0,0.0,0.0);
      GLfloat u,u_2,u_3,b[ CONTROL_POINTS ]; // 'b' stores the blending polynomials for bsplines
      for(int point = 0; point<len-3; point++){
         GLfloat curveVertex[3];
         // interpolate bsplines
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
   for(int i = 0; i < len; i++){
      glVertex3fv(&ctrlpts->at(i).at(0));
   }
   glEnd();

   glPopMatrix();
}

void drawGolgi(){
   // golgi is fixed and made completely as a texture inside a quad enclosure
   glPushMatrix();

      glTranslatef(-2.0,-2.5,0.25);
      glRotatef(135.0,0.0,0.0,1.0);
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D,textureGolgi);

      glBegin(GL_POLYGON);

      glTexCoord2f(0.0,0.0);
      glVertex3f(-GOLGI_DIMENSION/2.0,-GOLGI_DIMENSION/2.0,0.0);
      glTexCoord2f(1.0,0.0);
      glVertex3f(GOLGI_DIMENSION/2.0,-GOLGI_DIMENSION/2.0,0.0);
      glTexCoord2f(1.0,1.0);
      glVertex3f(GOLGI_DIMENSION,GOLGI_DIMENSION/2.0,0.0);
      glTexCoord2f(0.0,1.0);
      glVertex3f(-GOLGI_DIMENSION/2.0,GOLGI_DIMENSION/2.0,0.0);

      glEnd();
      glDisable(GL_TEXTURE_2D);
      
   glPopMatrix();
}

void display(void)
{

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   // cout << ctrlpoints.size();
   
   glPushMatrix();

   glScalef(scale.x,scale.y,scale.z); // magnification from normal level to cell level
   drawCell(ctrlpoints,ctrlpoints->size());
   tLysosomeChange++;
   for(int i=0;i<primaryLysosome.size();i++){
      drawPrimaryLysosome(primaryLysosome.at(i).lysosomeCtrlPoints,i);
   }
   if(tLysosomeChange>=TPRIMARY_LYSOSOME_TEXTURE_CHANGE_STEP){
      tLysosomeChange=0;
   }
   if(scale.x == scale.y && scale.x == 1.0){
      // hide till magnified till cell level
      if(tLysosomeMove<=CAPSULE_DISAPPEAR){
         drawCapsule(capsule.capsuleCtrlPoints);
      }
      // for debugging purposes      
      // drawPath(capsule.capsulePathCtrlPoints);
   }
   for(int i=0;i<mitochondria.size();i++){
      drawMitochondria(i);
   }
   drawGolgi();
   glScalef(1.0/scale.x,1.0/scale.y,1.0/scale.z);
   glPopMatrix();

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
      copyCurveConfiguration(primaryLysosome.at(0).pointsNow,primaryLysosome.at(0).pointsLater); // the destination configuration is our current configuration
      getCurvePoints(primaryLysosome.at(0).pointsLater,NUCLEUS_RADIUS,NUCLEUS_STANDARD_DEVIATION,NUCLEUS_SECTORS); // initial primary lysosome configuration
      copyCurveConfiguration(primaryLysosome.at(1).pointsNow,primaryLysosome.at(1).pointsLater); // the destination configuration is our current configuration
      getCurvePoints(primaryLysosome.at(1).pointsLater,NUCLEOLUS_RADIUS,NUCLEOLUS_STANDARD_DEVIATION,NUCLEOLUS_SECTORS); // initial primary lysosome configuration
      for(int i=2;i<primaryLysosome.size();i++){
         copyCurveConfiguration(primaryLysosome.at(i).pointsNow,primaryLysosome.at(i).pointsLater); // the destination configuration is our current configuration
         if(i==2 && tLysosomeMove>=CAPSULE_DISAPPEAR/2.0){
            // now change texture of primary lysosome as it has eaten up the capsule
            getCurvePoints(primaryLysosome.at(i).pointsLater,PRIMARY_LYSOSOME_RADIUS*1.5,PRIMARY_LYSOSOME_STANDARD_DEVIATION,PRIMARY_LYSOSOME_SECTORS); // initial primary lysosome configuration
         } else {            
            getCurvePoints(primaryLysosome.at(i).pointsLater,PRIMARY_LYSOSOME_RADIUS,PRIMARY_LYSOSOME_STANDARD_DEVIATION,PRIMARY_LYSOSOME_SECTORS); // initial primary lysosome configuration
         }
      }
      primaryLysosome.at(2).z = 1.0;
      tPrimaryLysosome=1.0/TPRIMARY_LYSOSOME_CHANGE_STEP;
   }
   for(int i=0;i<primaryLysosome.size();i++){
      getInterpolatedPointsForCurve(primaryLysosome.at(i).lysosomeCtrlPoints, primaryLysosome.at(i).pointsNow, primaryLysosome.at(i).pointsLater);
   }

   // for mitochondria
   fflush(stdout);
   if(tMitochondria<1.0){
      tMitochondria+=1.0/TMITOCHONDRIA_CHANGE_STEP;
   } else {
      copyCurveConfiguration(mitochondria.at(0).pointsNow,mitochondria.at(0).pointsLater); // the destination configuration is our current configuration
      getMitochondriaPoints(mitochondria.at(0).pointsLater,MITOCHONDRIA_RADIUS_X - getInBetween(0,MITOCHONDRIA_RADIUS_X*0.1),MITOCHONDRIA_RADIUS_Y - getInBetween(0,MITOCHONDRIA_RADIUS_Y*0.2),MITOCHONDRIA_STANDARD_DEVIATION_X,MITOCHONDRIA_STANDARD_DEVIATION_Y);
      copyCurveConfiguration(mitochondria.at(1).pointsNow,mitochondria.at(1).pointsLater); // the destination configuration is our current configuration
      getMitochondriaPoints(mitochondria.at(1).pointsLater,MITOCHONDRIA_RADIUS_X - getInBetween(MITOCHONDRIA_RADIUS_X*0.3,MITOCHONDRIA_RADIUS_X*0.5),MITOCHONDRIA_RADIUS_Y - getInBetween(MITOCHONDRIA_RADIUS_Y*0.3,MITOCHONDRIA_RADIUS_Y*0.5),MITOCHONDRIA_STANDARD_DEVIATION_X,MITOCHONDRIA_STANDARD_DEVIATION_Y);
      tMitochondria = 1.0/TMITOCHONDRIA_CHANGE_STEP;
   }
   for(int i=0;i<mitochondria.size();i++){
      getInterpolatedPointsForCurve(mitochondria.at(i).ctrlpoints, mitochondria.at(i).pointsNow, mitochondria.at(i).pointsLater);
   }

   // if magnification has reached cell level
   if(scale.x == scale.y && abs(scale.x-1.0) <= CELL_GROW_RATE){
      scale.x = 1.0;
      scale.y = scale.x;
      // move capsule
      
      if(tCapsule < capsule.capsulePathCtrlPoints->size()-3 && tLysosomeMove<CAPSULE_DISAPPEAR){
         // cout << tCapsule << " :: theta: " << capsule.theta << "(" << capsule.x <<"," << capsule.y <<")" << endl;
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
         
         // smooth rotation of capsule/food
         capsule.theta = CAPSULE_ROTATION_LAGRANGE * (tanh((capsule.y - old_y)/(capsule.x - old_x))*180.0/PI - 90.0) + (1- CAPSULE_ROTATION_LAGRANGE )*capsule.theta;
         
         tCapsule+=1.0/TCAPSULE_CHANGE_STEP;
         
      } else {
         if(tLysosomeMove==0.0){
            // initialize points as current points
            primaryLysosome.at(2).xdash = primaryLysosome.at(2).x;
            primaryLysosome.at(2).ydash = primaryLysosome.at(2).y;
         }
         if(tLysosomeMove<1.0){
            primaryLysosome.at(2).x = primaryLysosome.at(2).xdash*(1- tLysosomeMove) + capsule.x*tLysosomeMove;
            primaryLysosome.at(2).y = primaryLysosome.at(2).ydash*(1- tLysosomeMove) + capsule.y*tLysosomeMove;
            tLysosomeMove+=1/TPRIMARY_LYSOSOME_MOVE_STEP;
         }
      }
   } else {
      scale.x += CELL_GROW_RATE;
      scale.y = scale.x; 
   }
   glutPostRedisplay();
   if(!pause){
      // continue with the animation
      glutTimerFunc(100,timer,5);
   }
}

void keyboard(unsigned char key, int x, int y){
   switch(key) {
      case 's':
         // starts the animation
         if(!flag){
            pause = false;
            flag = true;
            timer(0);            
         }
         break;
      case 'p':
         // pause the animation
         pause = true;
         flag = false;
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
   // glEnable(GL_TEXTURE_2D);
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
   // glutTimerFunc(100,timer,5); // enable this to start animation without pressing the key
   glutMainLoop();
   return 0;
}