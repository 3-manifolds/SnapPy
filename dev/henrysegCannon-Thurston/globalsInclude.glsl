#version 300 es

//--------------------------------------------
//Global Constants
//--------------------------------------------
// const int MAX_MARCHING_STEPS = 127;
// const float MIN_DIST = 0.0;
// const float EPSILON = 0.0001;
// const vec4 ORIGIN = vec4(0,0,0,1);
//--------------------------------------------
//Global Variables
//--------------------------------------------
out vec4 out_FragColor;
// vec4 sampleEndPoint = vec4(1, 1, 1, 1);
// vec4 sampleTangentVector = vec4(1, 1, 1, 1);
// vec4 N = ORIGIN; //normal vector
// vec4 globalLightColor = ORIGIN;
// int hitWhich = 0;
//-------------------------------------------
//Translation & Utility Variables
//--------------------------------------------
// uniform int isStereo;
uniform vec2 screenResolution;
uniform float fov;

uniform int tetNum;
uniform mat4 currentBoost;
uniform float currentWeight;
// uniform mat4 stereoBoosts[2];
// uniform mat4 cellBoost; 
// uniform mat4 invCellBoost;
uniform int maxSteps;
uniform float maxDist;
uniform int subpixelCount;
uniform float edgeThickness;
uniform float contrast;
uniform int perspectiveType;
uniform int viewMode;
uniform int multiScreenShot;
uniform vec2 tile;
uniform vec2 numTiles;
//
uniform vec4 planes[##arrayLength##];   // ##arrayLength## gets replaced when we loadShaders in Main.js
uniform int otherTetNums[##arrayLength##]; 
uniform int entering_face_nums[##arrayLength##]; 
uniform float weights[##arrayLength##]; 
uniform mat4 SO31tsfms[##arrayLength##];

uniform float gradientThreshholds[5];
uniform vec3 gradientColours[5];

//--------------------------------------------
//Lighting Variables & Global Object Variables
//--------------------------------------------
// uniform vec4 lightPositions[4];
// uniform vec4 lightIntensities[6]; //w component is the light's attenuation -- 6 since we need controllers
// uniform int attnModel;
// uniform bool renderShadows[2];
// uniform float shadSoft;
// uniform sampler2D tex;
// uniform int controllerCount; //Max is two
// uniform mat4 controllerBoosts[2];
// uniform mat4 globalObjectBoost;
// uniform float globalObjectRadius;
//--------------------------------------------
//Scene Dependent Variables
//--------------------------------------------
// uniform vec4 halfCubeDualPoints[3];
// uniform float halfCubeWidthKlein;
// uniform float tubeRad;
// uniform vec4 cellPosition;
// uniform float cellSurfaceOffset;
// uniform vec4 vertexPosition;
// uniform float vertexSurfaceOffset;

// These are the planar mirrors of the fundamental simplex in the Klein (or analagous) model.
// Order is mirrors opposite: vertex, edge, face, cell.
// The xyz components of a vector give the unit normal of the mirror. The sense will be that the normal points to the outside of the simplex.
// The w component is the offset from the origin.
// uniform bool useSimplex;
// uniform vec4 simplexMirrorsKlein[4];
// uniform vec4 simplexDualPoints[4];

// The type of cut (1=sphere, 2=horosphere, 3=plane) for the vertex opposite the fundamental simplex's 4th mirror.
// These integers match our values for the geometry of the honeycomb vertex figure.
// We'll need more of these later when we support more symmetry groups.
// uniform int cut1;
// uniform int cut4;

// //Raymarch Functions
// float unionSDF(float d1, float d2){
//   return min(d1, d2);
// }

//--------------------------------------------------------------------
// Generalized Functions
//--------------------------------------------------------------------

// vec4 geometryNormalize(vec4 v, bool toTangent);
// vec4 geometryDirection(vec4 u, vec4 v);
// vec4 geometryFixDirection(vec4 u, vec4 v, mat4 fixMatrix);
// float geometryDot(vec4 u, vec4 v);
// float geometryDistance(vec4 u, vec4 v);
// float geometryNorm(vec4 v){
//   return sqrt(abs(geometryDot(v,v)));
// }

// vec4 pointOnGeodesic(vec4 u, vec4 vPrime, float dist);
// bool isOutsideCell(vec4 samplePoint, out mat4 fixMatrix);

//--------------------------------------------------------------------
// Generalized SDFs
//--------------------------------------------------------------------

// float globalSceneSDF(vec4 samplePoint, mat4 globalTransMatrix, bool collideWithLights);
// float localSceneSDF(vec4 samplePoint);

// float sphereSDF(vec4 samplePoint, vec4 center, float radius){
//   return geometryDistance(samplePoint, center) - radius;
// }