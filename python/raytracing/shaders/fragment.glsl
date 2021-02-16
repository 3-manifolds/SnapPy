#version 150

// GLSL version 1.50 corresponds to OpenGL 3.2 which is the version we target.

// History:
//
// This code was mostly written by Matthias Goerner during the
// "Illustrating Mathematics" workshop at ICERM in Fall 2019.
//
// Parts of it are taken from Cohomology fractals at
// https://github.com/henryseg/cohomology_fractals/tree/master/shaders
// The authors of Cohomology fractals and its history are listed at
// https://github.com/henryseg/cohomology_fractals/blob/master/README.md
//

#define COLOR_SCHEME 1

//--------------------------------------------
//Global Variables
//--------------------------------------------
out vec4 out_FragColor;

//-------------------------------------------
//Translation & Utility Variables
//--------------------------------------------

uniform vec2 screenResolution;
uniform float fov;

uniform int currentTetIndex;
uniform mat4 currentBoost;
uniform float currentWeight;

uniform int maxSteps;
uniform float maxDist;
uniform int subpixelCount;
uniform float edgeThickness;
uniform float contrast;

const int perspectiveTypeMaterial   = 0;
const int perspectiveTypeIdeal      = 1;
const int perspectiveTypeHyperideal = 2;

uniform int perspectiveType;
uniform int viewMode;
uniform int multiScreenShot;
uniform vec2 tile;
uniform vec2 numTiles;
uniform bool noGradient;
uniform bool showElevation = false;

// Convention: ##NAME## names a compile time constant.
// The string ##NAME## is replaced by the code in __init__.py
// by looking up its value in the dictionary constants_dict.

uniform int otherTetNums[4 * ##num_tets##]; 
uniform int otherFaceNums[4 * ##num_tets##]; 
uniform float weights[4 * ##num_tets##]; 

layout (std140) uniform TetrahedraBasics
{
    vec4 R13Vertices[4 * ##num_tets##];
    vec4 planes[4 * ##num_tets##];
    mat4 SO13tsfms[4 * ##num_tets##];
};

uniform float edgeTubeRadiusParam;

// +1 or -1 depending on orientation of tetrahedron.
uniform int orientations[##num_tets##];

uniform int face_color_indices[4 * ##num_tets##];
uniform int edge_color_indices[6 * ##num_tets##];
uniform int vertex_color_indices[4 * ##num_tets##];

uniform float gradientThreshholds[5];
uniform vec3 gradientColours[5];

uniform float lightBias;
uniform float lightFalloff;
uniform float brightness;
const vec4 lightSourcePosition = vec4(1.0, 0.0, 0.7, 0.0);

uniform bool isNonGeometric;
uniform sampler2D nonGeometricTexture;

const int num_tets = ##num_tets##;
const int num_edges = ##num_edges##;
const int num_cusps = ##num_cusps##;

#if ##finiteTrig##
layout (std140) uniform TetrahedraEdges
{
    // Ends of edge 01 are  0 + 12 * tetNum and  1 + 12 * tetNum
    // Ends of edge 02 are  2 + 12 * tetNum and  3 + 12 * tetNum
    // Ends of edge 12 are  4 + 12 * tetNum and  5 + 12 * tetNum
    // Ends of edge 03 ...
    // Ends of edge 13 ...
    // Ends of edge 23 are 10 + 12 * tetNum and 11 + 12 * tetNum
    //
    // Also see edgeToVertices
    vec4 R13EdgeEnds[12 * ##num_tets##];
};

uniform float vertexSphereRadiusParam;

#else

// For an incomplete cusp, the Margulis tube is a cylinder about
// the geodesic the triangulation spins about. In other words,
// it is the cylinder fixed by the peripheral group of that cusp.
// We encode it by its two end points and cosh(radius/2)^2/2.
layout (std140) uniform MargulisTubes
{
    vec4 margulisTubeTails[4 * ##num_tets##];
    vec4 margulisTubeHeads[4 * ##num_tets##];
};
uniform float margulisTubeRadiusParams[4 * ##num_tets##];

uniform float horosphereScales[4 * ##num_tets##];

// Heights of the Euclidean triangle obtained when intersecting
// horosphere with tetrahedron.
// There are four horospheres about the four edges of the tetrahedron,
// giving four triangles but they are all similar, so storing only
// one per tet.
uniform vec3 horotriangleHeights[##num_tets##];

// cosh(r)^2 where r is the radius of the sphere
// about the center of the tetrahedron.
uniform float insphereRadiusParams[##num_tets##];

// Matrix to convert between coordinates where the cusp is at
// infinity and the space of the tetrahedron
layout (std140) uniform TetCuspMatrices
{
    mat4 tetToCuspMatrices[4 * ##num_tets##];
    mat4 cuspToTetMatrices[4 * ##num_tets##];
};

uniform vec2 logAdjustments[4 * ##num_tets##];
uniform mat2 matLogs[4 * ##num_tets##];

const float peripheralCurveThickness = 0.015;

#if COLOR_SCHEME == 1
const vec3 longitudeColor = vec3( 1.0 , 1.0 , 1.0 );
const vec3 meridianColor  = vec3( 0.5 , 0.5 , 0.5 );
#else
const vec3 longitudeColor = vec3( 1.0 , 0.2 , 0.2 );
const vec3 meridianColor  = vec3( 0.2 , 1.0 , 0.2 );
#endif

#endif

// Colouring function. All components are in the range [0...1], including hue.
// from http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl
vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

// Lorentz dot product with signature -+++
float
R13Dot(vec4 u, vec4 v)
{
    return - u.x*v.x + u.y*v.y + u.z*v.z + u.w*v.w;
}

vec4
R13Normalise(vec4 v)
{
    return v * inversesqrt(abs(R13Dot(v,v)));
}

// Given a direction vector and a point on the hyperbolic model, produces
// a unit tangent vector (in Minkowski space) at that point.
//
// It does by subtracting a suitable multiple of point from direction to
// make it orthogonal to point and then normalising.
vec4
makeUnitTangentVector(vec4 dir, vec4 point)
{
    // Make dir orthogonal to point.
    // Note that R13Dot(point, point) == -1, so we add instead of
    // subtract.
    vec4 t = dir + R13Dot(dir, point) * point;

    return R13Normalise(t);
}

float
geodesicParameterPlanes(vec4 samplePoint, vec4 dualPoint1, vec4 dualPoint2){
  // "distance" from a geodesic defined by two (assumed perpendicular) geodesic planes, this is not quite distance, need to asinh(sqrt( result ))

  float dot1 = -R13Dot(samplePoint, dualPoint1);
  vec4 dualPointPerp = R13Normalise(dualPoint2 - R13Dot(dualPoint1, dualPoint2) * dualPoint1); // should be precalculated if this is a main feature
  float dot2 = -R13Dot(samplePoint, dualPointPerp);

  return dot1*dot1 + dot2*dot2;
}

float
triangleBdryParam(vec4 samplePoint, int tetNum, int exit_face){
  vec4 exit_dual_point = planes[4*tetNum + exit_face];
  float smallest_p = 100000000.0;
  for(int face=0; face<4; face++){
      if(face != exit_face){  // find p when we hit that face
          int index = 4*tetNum + face;
          float new_p = geodesicParameterPlanes(samplePoint, exit_dual_point, planes[index]);
          if(new_p < smallest_p){
            smallest_p = new_p;
          }   
      }
  }
  return smallest_p;
}

/// --- Ray-trace code --- ///

// Kind of object a ray hit.

const int object_type_nothing             = 0;
const int object_type_face                = 1;
const int object_type_edge_cylinder_enter = 2;
const int object_type_edge_cylinder_exit  = 3;
const int object_type_horosphere          = 4;
const int object_type_edge_fan            = 5;
const int object_type_insphere            = 6;
const int object_type_vertex_sphere       = 7;
const int object_type_margulis_tube       = 8;
const int object_type_elevation_enter     = 9;
const int object_type_elevation_exit      = 10;

// A ray consists of a point in the hyperbolid model and a
// unit tangent vector dir orthogonal to the point with respect
// to the Lorentz product.
struct Ray
{
    vec4 point;
    vec4 dir;
};

// Ray with extra information what object was hit.
struct RayHit
{
    Ray ray;
    // What tetrahedron we are in.
    int tet_num;
    vec4 light_source;
    // Distance the ray traveled from eye so far
    float dist;
    float weight;
    float distWhenLeavingCusp;
    // Type of object hit
    int object_type;
    // Index of object (within the tetrahedron), e.g.,
    // 0-3 for horospheres, 0-5 for edges.
    int object_index;
};

// We distinguish between:
// - colored ray hits: the geometry is lit, e.g., the cylinder
//   about an edge reacts to light)
// - valued ray hits: we use a value such as the distance or cohomology
//   fractal weight to compute a color using a color gradient
//
// Note that for subsampling, we need to average the value before
// applying the color gradient.
//
bool isColored(RayHit ray_hit)
{
    return
        ray_hit.object_type == object_type_vertex_sphere ||
        ray_hit.object_type == object_type_insphere ||
        ray_hit.object_type == object_type_horosphere ||
        ray_hit.object_type == object_type_edge_cylinder_enter ||
        ray_hit.object_type == object_type_edge_cylinder_exit ||
        ray_hit.object_type == object_type_margulis_tube ||
        ray_hit.object_type == object_type_edge_fan ||
        ray_hit.object_type == object_type_elevation_enter ||
        ray_hit.object_type == object_type_elevation_exit;
}

// Advances ray by distance atanh(p).
//
// It is often more convenient to work with tanh(distance) rather than
// distance and we refer to tanh(distance) as "distParam".
void
advanceRayByDistParam(inout Ray ray, float p)
{
    ray.point = R13Normalise(ray.point + p * ray.dir );
    ray.dir = makeUnitTangentVector(ray.dir, ray.point);
}

// Returned as distParam by some methods to indicate that ray does not
// intersect object.
const float unreachableDistParam = 1000.0;

// Returns the real roots of the equation a * x^2 + b * x + c = 0, but
// each root is returned only if it is less than min_val.
// The existence of no such root is indicated by returning
// unreachableDistParam.
vec2
realRootsOfQuadratic(float a, float b, float c,
                     float min_val)
{
    float d = b * b - 4 * a * c;
    if (d < 0) {
        return vec2(unreachableDistParam, unreachableDistParam);;
    }
    float offset = sign(a) * sqrt(d);
    vec2 result = vec2(-b - offset, -b + offset) / (2 * a);
    if (result.y < min_val) {
        return vec2(unreachableDistParam, unreachableDistParam);;
    }
    if (result.x < min_val) {
        return vec2(unreachableDistParam, result.y);
    }
    return result;
}

//--------------------------------------------
// Ray intersections
//
// Functions to compute the distParam (tanh(distance)) for the intersection
// of a ray with some object. Only the intersection of the ray entering an
// object is returned. No intersection is indicated by unreachableDistParam.

// Intersection with plane.
// The plane is given by all points such that the inner product
// of the point and the planeEqn is zero. 
float
distParamForPlaneIntersection(Ray ray,
                              vec4 planeEqn)
{
    // solve:
    //   R13Dot(planeEqn, ray.point  + p *                   ray.dir) = 0
    //   R13Dot(planeEqn, ray.point) + p * R13Dot(planeEqn, ray.dir) = 0

    float denom = R13Dot(planeEqn, ray.dir);
    if(denom == 0.0) {
        return unreachableDistParam;
    }
    return -R13Dot(planeEqn, ray.point) / denom;
}

// Intersection with sphere about given center (in hyperboloid model).
// The last parameter is the cosh(radius)^2.
vec2
distParamsForSphereIntersection(Ray ray,
                               vec4 center, float sphereRadiusParam)
{
    float startDot = R13Dot(center, ray.point);
    float dirDot   = R13Dot(center, ray.dir);
    
    return realRootsOfQuadratic(
        dirDot * dirDot + sphereRadiusParam,
        2.0 * dirDot * startDot,
        startDot * startDot - sphereRadiusParam,
        0.0);
}

// Intersection with horosphere defined by given light-like vector.
// The horosphere consists of all points such that the inner product
// with the light-like vector is -1.
vec2
distParamsForHorosphereIntersection(Ray ray,
                                   vec4 horosphere)
{
    return distParamsForSphereIntersection(ray, horosphere, 1);
}

// Intersection with cylinder about the geodesic with the two given
// (light-like) endpoints.
// tubeRadiusParam is cosh(radius/2)^2/2.
// This function can detect intersections of the ray that happen
// some distance before the start point of the ray. To use this feature,
// set minDistParam to a negative value, namely to tanh(-distance),
// where distance is how far back we want to track the ray.
vec2
distParamsForTubeIntersection(Ray ray,
                             vec4[2] endpoints,
                             float tubeRadiusParam,
                             float minDistParam)
{
    float start0Dot = R13Dot(endpoints[0], ray.point);
    float dir0Dot   = R13Dot(endpoints[0], ray.dir);
    float start1Dot = R13Dot(endpoints[1], ray.point);
    float dir1Dot   = R13Dot(endpoints[1], ray.dir);
    float endDot    = R13Dot(endpoints[0], endpoints[1]);

    return realRootsOfQuadratic(
        dir0Dot * dir1Dot - endDot * tubeRadiusParam,
        start0Dot * dir1Dot + start1Dot * dir0Dot,
        start0Dot * start1Dot + endDot * tubeRadiusParam,
        minDistParam);
}

//--------------------------------------------
// Object normals
//
// Function to compute the normal of an object at the given point


// Normal for a sphere about the center.
vec4
normalForSphere(vec4 point, vec4 center)
{
    vec4 t = center - point;

    return makeUnitTangentVector(t, point);
}

// Normal for a cylinder between the two given (light-like) endpoints.
vec4
normalForTube(vec4 point, vec4[2] endpoints)
{
    vec4 t =   endpoints[0] * R13Dot(point, endpoints[1])
             + endpoints[1] * R13Dot(point, endpoints[0]);

    return - makeUnitTangentVector(t, point);
}

//--------------------------------------------
// Object helpers
//
// Various helpers to obtain object representations from
// the uniforms.

// Indices of vertices that are the end points of an edge.
//
// Needs to be consistent with SnapPy/t3m conventions
// (i.e., t3m's OneSubsimplices).
const ivec2[6] edgeToVertices = ivec2[](ivec2(0, 1),
                                        ivec2(0, 2),
                                        ivec2(1, 2),
                                        ivec2(0, 3),
                                        ivec2(1, 3),
                                        ivec2(2, 3));

// The two endpoints of an edge of a tetrahedron
vec4[2]
endpointsForEdge(int tet, int edge)
{
#if ##finiteTrig##
    return vec4[](R13EdgeEnds[12 * tet + 2 * edge    ],
                  R13EdgeEnds[12 * tet + 2 * edge + 1]);
#else
    return vec4[](R13Vertices[4 * tet + edgeToVertices[edge].x],
                  R13Vertices[4 * tet + edgeToVertices[edge].y]);
#endif
}

#if !##finiteTrig##
// The two endpoints of a Margulis tube.
vec4[2]
endpointsForMargulisTube(int index)
{
    return vec4[](margulisTubeTails[index],
                  margulisTubeHeads[index]);
}

// The equation for the horosphere about a vertex in a tetrahedron.
// index is 4 * tetIndex + vertex.
vec4
horosphereEqn(int index)
{
    return horosphereScales[index] * R13Vertices[index];
}
#endif

vec4
normalForRayHit(RayHit ray_hit)
{
#if ##finiteTrig##
    if(ray_hit.object_type == object_type_vertex_sphere) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        return normalForSphere(ray_hit.ray.point, R13Vertices[index]);
    }
#else
    if(ray_hit.object_type == object_type_insphere) {
        return normalForSphere(ray_hit.ray.point, vec4(1,0,0,0));
    }

    if(ray_hit.object_type == object_type_horosphere) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        return horosphereEqn(index) - ray_hit.ray.point;
    }

    if(ray_hit.object_type == object_type_margulis_tube) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        return normalForTube(
            ray_hit.ray.point,
            endpointsForMargulisTube(index));
    }
#endif    
    
    if(ray_hit.object_type == object_type_edge_fan ||
       ray_hit.object_type == object_type_elevation_enter ||
       ray_hit.object_type == object_type_elevation_exit) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        return planes[index];
    }

    if(ray_hit.object_type == object_type_edge_cylinder_enter) {
        return normalForTube(
            ray_hit.ray.point,
            endpointsForEdge(ray_hit.tet_num, ray_hit.object_index));
    }

    if(ray_hit.object_type == object_type_edge_cylinder_exit) {
        return - normalForTube(
            ray_hit.ray.point,
            endpointsForEdge(ray_hit.tet_num, ray_hit.object_index));
    }

    return vec4(0,1,0,0);
}

// Convert point in hyperboloid model to upper halfspace
// model.
// The vec3 result corresponds to result.x + result.y * i + result.z * j,
// i.e., the last component of the vec3 is the height.
vec3
hyperboloidToUpperHalfspace(vec4 h)
{
    vec3 klein = h.yzw / h.x;
    vec3 poincare = klein / (1.0 + sqrt(1.0 - dot(klein, klein)));
    vec3 denom_helper = vec3(poincare.x - 1.0, poincare.yz);
    float denom = dot(denom_helper, denom_helper);
    
    return vec3(2.0 * poincare.yz, 1.0 - dot(poincare, poincare)) / denom;   
}

#if !##finiteTrig##
// Compute the coordinates of a ray hit on a horosphere in the upper
// half space model such that the cusp is at infinity.
vec3
preferredUpperHalfspaceCoordinates(RayHit ray_hit)
{
    int index = 4 * ray_hit.tet_num + ray_hit.object_index;

    return hyperboloidToUpperHalfspace(
        ray_hit.ray.point * tetToCuspMatrices[index]);
}

vec2
complexLog(vec2 z)
{
    return vec2(log(length(z)), atan(z.y, z.x));
}

vec2
MLCoordinatesForRayHit(RayHit rayHit)
{
    int index = 4 * rayHit.tet_num + rayHit.object_index;
    
    vec3 pointUpperHalfspace = preferredUpperHalfspaceCoordinates(rayHit);
    vec2 z = pointUpperHalfspace.xy;

    if (rayHit.object_type == object_type_margulis_tube) {
        z = complexLog(z) + logAdjustments[index];
    }

    return z * matLogs[index];
}
#endif

// Compute the SO13 transform corresponding to the PSL(2,C)-matrix
// [[1, z], [0, 1]].
// Special case of the kernel's Moebius_to_O31 for upper unit triangular
// matrices.
mat4
parabolicSO13(vec2 z)
{
    float t = dot(z, z) / 2.0;
    
    return mat4( 1.0 + t,     - t,     z.x,     z.y,
                       t, 1.0 - t,     z.x,     z.y,
                     z.x,    -z.x,     1.0,     0.0,
                     z.y,    -z.y,     0.0,     1.0 );
}

// Compute the SO13 transform corresponding to the PGL(2,C)-matrix
// [[ exp(z), 0], [0, 1]].
// Special case of the kernel's Moebius_to_O31 for diagonal matrices.
mat4
loxodromicSO13(vec2 z)
{
    return mat4( cosh(z.x), sinh(z.x),       0.0,       0.0,
                 sinh(z.x), cosh(z.x),       0.0,       0.0,
                       0.0,       0.0,  cos(z.y), -sin(z.y),
                       0.0,       0.0,  sin(z.y),  cos(z.y) );
}

void
ray_trace_through_hyperboloid_tet(inout RayHit ray_hit)
{
    int entry_object_type = ray_hit.object_type;
    int entry_object_index = ray_hit.object_index;
 
    ///Given shape of a tet and a ray, find where the ray exits and through which face
    float smallest_p = 100000000.0;

    for(int face = 0; face < 4; face++) {
        if (entry_object_type != object_type_face || entry_object_index != face) {
            // find p when we hit that face
            int index = 4 * ray_hit.tet_num + face;
            if(R13Dot(ray_hit.ray.dir, planes[index]) > 0.0){ 
                float p = distParamForPlaneIntersection(ray_hit.ray, planes[index]);
                // if ((-10000.0 <= p) && (p < smallest_p)) {
                if (p < smallest_p) {  
                    /// negative values are ok if we have to go backwards a little to get through the face we are a little the wrong side of
                    /// Although this can apparently get caught in infinite loops in an edge

                    /// if we are on an edge then we don't in fact move as we go through this tet: t = 0.0
                    /// also allow tiny negative values, which will come up from floating point errors. 
                    /// surface normals check should ensure that even in this case we make progress through 
                    /// the triangles around an edge
                    smallest_p = p;
                    ray_hit.object_type = object_type_face;
                    ray_hit.object_index = face;
                }
            }
        }
    }

#if ##finiteTrig##
    if (vertexSphereRadiusParam > 1.0001) {
        for (int vertex = 0; vertex < 4; vertex++) {
            int index = 4 * ray_hit.tet_num + vertex;
            float p = distParamsForSphereIntersection(
                ray_hit.ray,
                R13Vertices[index],
                vertexSphereRadiusParam).x;
            if (p < smallest_p) {
                smallest_p = p;
                ray_hit.object_type = object_type_vertex_sphere;
                ray_hit.object_index = vertex;
            }
        }
    }
#else
    {
        float r = insphereRadiusParams[ray_hit.tet_num];

        if (r > 1.0001) {
            float p = distParamsForSphereIntersection(
                ray_hit.ray,
                vec4(1,0,0,0),
                r).x;
            if (p < smallest_p) {
                smallest_p = p;
                ray_hit.object_type = object_type_insphere;
                ray_hit.object_index = 0;
            }
        }
    }

    for (int vertex = 0; vertex < 4; vertex++) {
        int index = 4 * ray_hit.tet_num + vertex;
        if (horosphereScales[index] != 0.0) {
            vec2 params = distParamsForHorosphereIntersection(ray_hit.ray,
                                                              horosphereEqn(index));
            if (params.x < smallest_p) {
                smallest_p = params.x;
                ray_hit.object_type = object_type_horosphere;
                ray_hit.object_index = vertex;
            }
        }

        if (margulisTubeRadiusParams[index] > 0.50001) {
            vec2 params = distParamsForTubeIntersection(
                ray_hit.ray,
                endpointsForMargulisTube(index),
                margulisTubeRadiusParams[index],
                0.0);

            if (params.x < smallest_p) {
                smallest_p = params.x;
                ray_hit.object_type = object_type_margulis_tube;
                ray_hit.object_index = vertex;
            }
        }
    }

#endif

    if (edgeTubeRadiusParam > 0.50001) {

#if ##finiteTrig##
        float backDistParam = 0.0;
#else
        float backDistParam = tanh(ray_hit.distWhenLeavingCusp-ray_hit.dist);
#endif    

        for (int edge = 0; edge < 6; edge++) {

            vec2 params = distParamsForTubeIntersection(
                ray_hit.ray,
                endpointsForEdge(ray_hit.tet_num, edge),
                edgeTubeRadiusParam,
                backDistParam);

            if (params.x < smallest_p) {
                smallest_p = params.x;
                ray_hit.object_type = object_type_edge_cylinder_enter;
                ray_hit.object_index = edge;
            } else if (params.y < smallest_p) {
                smallest_p = params.y;
                ray_hit.object_type = object_type_edge_cylinder_exit;
                ray_hit.object_index = edge;
            }
        }
    }

    ray_hit.dist += atanh(smallest_p);
    advanceRayByDistParam(ray_hit.ray, smallest_p);

    if(edgeThickness > 0.00001) {
        if (ray_hit.object_type == object_type_face) {
            if(triangleBdryParam(ray_hit.ray.point, ray_hit.tet_num, ray_hit.object_index) < edgeThickness) {
                ray_hit.object_type = object_type_edge_fan;
            }
        }
    }
}

int
is_elevation_hit(float old_weight, float new_weight)
{
    const float eps = 1e-4;
    float o = old_weight - eps;
    float n = new_weight - eps;

    if (o * n < 0.0) {
        if (n < 0.0) {
            return object_type_elevation_enter;
        } else {
            return object_type_elevation_exit;
        }
    }

    return object_type_nothing;
}

void
ray_trace(inout RayHit ray_hit) {

    for(int i = 0; i < maxSteps; i++){
        ray_trace_through_hyperboloid_tet(ray_hit);

        if (ray_hit.object_type != object_type_face) {
            break;
        }

        if (ray_hit.dist > maxDist) {
            break;
        }

        // in fact pow(sinh(radius in hyperbolic units),2.0). However, sinh^2 is monotonic for 
        // positive values so we get correct behaviour by comparing without the sinh^2. 
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;

        float new_weight = ray_hit.weight + weights[ index ];        

        if (showElevation) {
            int elevation_object_type =
                is_elevation_hit(ray_hit.weight, new_weight);

            if (elevation_object_type != object_type_nothing) {
                ray_hit.object_type = elevation_object_type;
                break;
            }
        }

        ray_hit.weight = new_weight;

        ray_hit.object_index = otherFaceNums[ index ];
        mat4 tsfm = SO13tsfms[ index ];

        ray_hit.light_source = ray_hit.light_source * tsfm;
        ray_hit.ray.point = ray_hit.ray.point * tsfm;
        ray_hit.ray.dir = R13Normalise( ray_hit.ray.dir * tsfm ); 
        ray_hit.tet_num = otherTetNums[ index ];
    }
}



/// --- Colour gradient code --- ///

int find_band(float t, float threshholds[5]){
    for(int j = 1; j < 4; j++) {
        if(t < threshholds[j]) {
            return j;
        }
    }
    return 4;
}

vec3 general_gradient(float t, float threshholds[5], vec3 colours[5]){
    int i = find_band(t, threshholds);
    return mix(colours[i-1],
               colours[i],
               (t - threshholds[i-1])/(threshholds[i] - threshholds[i-1]));
}

// Given the (average) value for a ray hit, apply gradient to get
// color.
vec3 colorForValue(float value)
{
    if (noGradient) {
        return vec3(value);
    }

    value = contrast * value;
    value = 0.5 + 0.5 * value/ (abs(value) + 1.0);  //faster than atan, similar

    return general_gradient(value, gradientThreshholds, gradientColours);
}

struct MaterialParams
{
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

    float shininess;
};

MaterialParams
material_params(RayHit ray_hit)
{
    MaterialParams result;

    result.diffuse  = vec3(0.2, 0.6,  0.3);
    result.ambient  = 0.5 * result.diffuse;
    result.specular  = vec3(0.5, 0.5, 0.5);
    result.shininess = 20;

#if ##finiteTrig##
    if (ray_hit.object_type == object_type_vertex_sphere) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = vertex_color_indices[index];

        result.diffuse = hsv2rgb(vec3(float(color_index)/float(num_cusps), 0.25, 1.0));
        result.ambient = 0.5 * result.diffuse;
    }        
#else
    if (ray_hit.object_type == object_type_horosphere ||
        ray_hit.object_type == object_type_margulis_tube) {
        
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = vertex_color_indices[index];

#if COLOR_SCHEME == 1
        result.diffuse = hsv2rgb(vec3(float(color_index)/float(num_cusps), 0.25, 1.0));
#else
        result.diffuse =
            vec3(0.5, 0.5, 0.5)
            + sin(color_index) * vec3( 0.3, -0.3, 0.0)
            + cos(color_index) * vec3(0.15, 0.15, -0.3);
#endif
        result.ambient = 0.5 * result.diffuse;

        vec2 coords = fract(MLCoordinatesForRayHit(ray_hit));
        
        if (coords.x <       peripheralCurveThickness ||
            coords.x > 1.0 - peripheralCurveThickness) {
            result.diffuse = longitudeColor;
            result.ambient = result.diffuse;
        }
        if (coords.y <       peripheralCurveThickness ||
            coords.y > 1.0 - peripheralCurveThickness) {
            result.diffuse = meridianColor;
            result.ambient = result.diffuse;
        }
    }

    if (ray_hit.object_type == object_type_insphere) {
        result.diffuse = hsv2rgb(vec3(float(ray_hit.tet_num)/float(num_tets), 0.5, 1.0));

        result.diffuse *= 0.5;

        result.ambient = 0.5 * result.diffuse;
    }

#endif

    if (ray_hit.object_type == object_type_edge_fan) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = face_color_indices[index];
        result.diffuse = hsv2rgb(vec3(float(color_index)/float(2*num_tets), 0.75, 0.5));
        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_edge_cylinder_enter) {
        int index = 6 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = edge_color_indices[index];
        
        //using num_tets = num_edges

#if COLOR_SCHEME == 1
        result.diffuse = hsv2rgb(vec3(float(color_index)/float(num_edges), 1.0, 1.0));
#else
        result.diffuse =
            vec3(0.5, 0.5, 0.5)
            + sin(color_index) * vec3( 0.3,  -0.3,  0.0)
            + cos(color_index) * vec3(0.15, 0.15, -0.3);
#endif

        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_edge_cylinder_exit) {
        int index = 6 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = edge_color_indices[index];
        
        //using num_tets = num_edges
        result.diffuse = 0.3 * hsv2rgb(vec3(float(color_index)/float(num_tets), 1.0, 1.0));
        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_elevation_enter) {
        result.diffuse = vec3(0.3,0.7,0.3);
        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_elevation_exit) {
        result.diffuse = vec3(0.7,0.3,0.3);
        result.ambient = 0.5 * result.diffuse;
    }

    return result;
}

// Compute the value for a valued ray hit (see isColored for explanation)
float valueForRayHit(RayHit ray_hit)
{
    if (viewMode == 0) {
        return ray_hit.weight;
    } else if (viewMode == 1) {
        return 0.5 * ray_hit.dist;
    } else {
        return float(ray_hit.tet_num);
    }
}

// Apply lighting for a colored ray hit (see isColored for explanation)
vec3 colorForRayHit(RayHit ray_hit)
{
    MaterialParams material = material_params(ray_hit);

    vec4 normal = normalForRayHit(ray_hit);

    vec4 light_position = R13Normalise(ray_hit.light_source);

    // Distance of light source to origin where the eye ray started
    float light_dist_origin = acosh(-R13Dot(R13Normalise(vec4(1,0,0.7,0)), vec4(1, 0, 0, 0)));

    // Distance of light source to ray hit
    float unsafe_dist = acosh(-R13Dot(ray_hit.ray.point, light_position));

    // Use triangle inequality to limit distance
    float dist = clamp(unsafe_dist,
                       ray_hit.dist - light_dist_origin,
                       ray_hit.dist + light_dist_origin);


    vec4 light_dir_at_hit = makeUnitTangentVector(
        - light_position, ray_hit.ray.point);
    
    float normal_light = clamp(R13Dot(normal, light_dir_at_hit), 0, 1);

    vec4 half_angle = R13Normalise(light_dir_at_hit + ray_hit.ray.dir);

    float blinn_term =
        normal_light > 0.0
        ? pow(clamp(R13Dot(half_angle, normal), 0, 1), material.shininess)
        : 0.0;

    return  brightness * (material.ambient
           + material.diffuse * normal_light
           + material.specular * blinn_term ) / pow((dist + lightBias) / lightBias, lightFalloff);
}

/// --- Graph-trace code --- ///

int faceFurthest(vec4 v, int tet_num, int entry_face)
{
    int result = -1;
    float biggest_amount = 0.0000001;
    for (int face = 0; face < 4; face++) {
        if (entry_face != face) {
            float amount = R13Dot( v, planes[4 * tet_num + face] );
            if (amount > biggest_amount) {
                biggest_amount = amount;
                result = face;
            }
        }
    }
    return result;
}

void graph_trace(inout RayHit ray)
{
  int entry_face = -1;
  mat4 tsfm = mat4(1.0);

  for(int i = 0; i < maxSteps; i++) {
      int face = faceFurthest(ray.ray.point, ray.tet_num, entry_face);
      if (face == -1) {
          break;
      }

      int index = 4 * ray.tet_num + face;
      entry_face = otherFaceNums[ index ];
      ray.tet_num = otherTetNums[ index ];
      ray.weight += weights[ index ];
      ray.ray.point = ray.ray.point * SO13tsfms[ index ];
      ray.ray.dir = ray.ray.dir * SO13tsfms[ index ];
      ray.light_source = ray.light_source * SO13tsfms[ index ];
  }
}

/// --- Ray init pt and directions code --- ///

Ray get_ray_eye_space(vec2 xy)
{
    Ray result;

    if (perspectiveType == perspectiveTypeMaterial) {
        result.point = vec4(1.0,0.0,0.0,0.0);
        result.dir = R13Normalise(vec4(0.0, 2.0 * xy, -1.0));
    } else if (perspectiveType == perspectiveTypeIdeal) {
        // parabolic transformation magic by Saul
        float r2 = 0.5 * dot(xy, xy);
        result.point = vec4(r2 + 1.0, xy, r2);
        result.dir   = vec4(r2,       xy, r2 - 1.0);
    } else { // perspectiveTypeHyperIdeal
        result.point = R13Normalise(vec4(1.0, 2.0 * xy, 0.0));
        result.dir = vec4(0.0, 0.0, 0.0, -1.0);
    }
    
    return result;
}

#if ##finiteTrig##
bool
leaveVertexNeighborhood(inout RayHit rayHit)
{
    if (vertexSphereRadiusParam <= 1.0001) {
        return false;
    }

    float smallest_p = unreachableDistParam;

    // For all vertices
    for (int vertex = 0; vertex < 4; vertex++) {
        int index = 4 * rayHit.tet_num + vertex;
        vec2 params = distParamsForSphereIntersection(
            rayHit.ray,
            R13Vertices[index],
            vertexSphereRadiusParam);
        if (params.x == unreachableDistParam) {
            if (params.y < smallest_p) {
                smallest_p = params.y;
                rayHit.object_type = object_type_vertex_sphere;
                rayHit.object_index = vertex;
            }
        }
    }

    if (smallest_p < unreachableDistParam) {
        rayHit.dist += atanh(smallest_p);
        rayHit.distWhenLeavingCusp = rayHit.dist;
        advanceRayByDistParam(rayHit.ray, smallest_p);
        graph_trace(rayHit);
    }

    return false;
}

#else

// Determine whether the ray starts in a horosphere or
// Margulis tube.
// If yes, move the ray to the point where we exit the
// horosphere and set rayHit.object_type to horosphere.
//
// The result is true if we are inside a horosphere AND
// the ray is hitting a peripheral curve on the horosphere.
// 
// For optimization, leaveVertexNeighborhood will also apply a
// parabolic transformation to the ray trying to bring the
// where we exit the horosphere closer to the entry point.
bool
leaveVertexNeighborhood(inout RayHit rayHit)
{
    float smallest_p = unreachableDistParam;

    // For all vertices
    for (int vertex = 0; vertex < 4; vertex++) {
        int index = 4 * rayHit.tet_num + vertex;
        // corresponding to complete cusps
        if (horosphereScales[index] != 0.0) {
            vec2 params = distParamsForHorosphereIntersection(
                rayHit.ray, horosphereEqn(index));
            if (params.x == unreachableDistParam) {
                // We are in the horosphere
                if (params.y < smallest_p) {
                    // Remember this
                    smallest_p = params.y;
                    rayHit.object_type = object_type_horosphere;
                    rayHit.object_index = vertex;
                }
            }
        } else if (margulisTubeRadiusParams[index] > 0.50001) {
            vec2 params = distParamsForTubeIntersection(
                rayHit.ray,
                endpointsForMargulisTube(index),
                margulisTubeRadiusParams[index],
                0.0);
            if (params.x == unreachableDistParam) {
                if (params.y < smallest_p) {
                    smallest_p = params.y;
                    rayHit.object_type = object_type_margulis_tube;
                    rayHit.object_index = vertex;
                }
            }
        }
    }
    
    // We are in a horosphere.
    if (smallest_p < unreachableDistParam) {

        // Book-keeping and advancing the ray to the exit
        // point
        rayHit.dist += smallest_p < 1.0 ? atanh(smallest_p) : 20.0;
        rayHit.distWhenLeavingCusp = rayHit.dist;
        advanceRayByDistParam(rayHit.ray, smallest_p);

        int index = 4 * rayHit.tet_num + rayHit.object_index;

        vec2 ml = MLCoordinatesForRayHit(rayHit);

        // Compute the coordinates of exit point in upper half space
        // such that the cusp is at infinity
        // Use complex part ignoring height in upper halfspace

        // Map R^2->torus
        vec2 coords = fract(ml);

        // Check whether we hit peripheral curve
        if (coords.x <       peripheralCurveThickness ||
            coords.x > 1.0 - peripheralCurveThickness ||
            coords.y <       peripheralCurveThickness ||
            coords.y > 1.0 - peripheralCurveThickness) {
            // Hit peripheral curve
            return true;
        }

        // Compute suitable multiple of merdian and longitude translation
        // bringing the exit point into the fundamental parallelogram
        // near zero.
        vec2 c = -round(ml) * inverse(matLogs[index]);

        mat4 tsfmCuspSpace =
            (rayHit.object_type == object_type_horosphere)
            ? parabolicSO13(c)
            : loxodromicSO13(c);
        
        // Convert O13 matrix from space where cusp was at infinity
        // to space of tetrahedron
        mat4 tsfm =
            tetToCuspMatrices[index] *
            tsfmCuspSpace *
            cuspToTetMatrices[index];
        
        // And apply transformation to ray.
        rayHit.light_source = rayHit.light_source * tsfm;
        rayHit.ray.point = rayHit.ray.point * tsfm;
        rayHit.ray.dir = R13Normalise( rayHit.ray.dir * tsfm ); 

        // If we are inside a horosphere, leaveVertexNeighborhood has computed
        // the point where we leave the horosphere. But that point
        // might not be inside the current tetrahedron, so fix it.
        graph_trace(rayHit);
    }

    return false;
}

#endif

RayHit computeRayHit(vec2 xy){
    Ray ray_eye_space = get_ray_eye_space(xy);

    RayHit ray_tet_space;
    ray_tet_space.ray.point = ray_eye_space.point * currentBoost;
    ray_tet_space.ray.dir   = ray_eye_space.dir   * currentBoost;
    ray_tet_space.dist = 0.0;
    ray_tet_space.distWhenLeavingCusp = 0.0;
    ray_tet_space.weight = currentWeight;
    ray_tet_space.tet_num = currentTetIndex;
    ray_tet_space.light_source = R13Normalise(lightSourcePosition * currentBoost);
    ray_tet_space.object_type = object_type_nothing;
    ray_tet_space.object_index = -1;

    // If using a camera where the ray's do not
    // all start from a common point, transform ray first
    // to be inside a tetrahedron.
    if (perspectiveType != perspectiveTypeMaterial) {
        graph_trace(ray_tet_space);
    }

    // Check whether we are in a horosphere and if yes,
    // whether the ray hit a peripheral curve.
    bool hitPeripheral = leaveVertexNeighborhood(ray_tet_space);

    if (hitPeripheral) {
        // If we hit a peripheral curve, leaveVertexNeighborhood has given us
        // the intersection point and we can immediately shade.
    } else {
        // In all other cases, we need to raytrace before we shade.
        ray_trace(ray_tet_space);
    }
    
    return ray_tet_space;
}

vec3 sampleNonGeometricTexture(vec2 fragCoord)
{
    vec2 coord = gl_FragCoord.xy - 0.5 * screenResolution.xy;
    coord.x /=  320.0;
    coord.y /= -100.0;
    coord += vec2(0.5, 0.5);

    if (coord.x < 0.002 || coord.x > 0.99 || coord.y < 0.01 || coord.y > 0.99) {
        return vec3(0.0);
    }

    return texture(nonGeometricTexture, coord).xyz;
}


void main(){
     
    // Show text "Non-geoemtric"
    if (isNonGeometric) {
        out_FragColor = vec4(sampleNonGeometricTexture(gl_FragCoord.xy), 1.0);
        return;
    }

    vec2 xy = (gl_FragCoord.xy - 0.5*screenResolution.xy)/screenResolution.x;
    if(multiScreenShot == 1) {
        // Return multiple 4096x4096 screenshots that can be combined in, e.g. Photoshop.
        // Here screenResolution is really tileResolution;
        xy = (xy + tile - 0.5*(numTiles - vec2(1.0,1.0))) / numTiles.x;
    }

    float min_depth = 1;
    vec3 total_color = vec3(0);
    int num_valued_subpixels = 0;
    float total_value = 0;

    for(int i=0; i<subpixelCount; i++){
        for(int j=0; j<subpixelCount; j++){
            vec2 offset =
                ( vec2(float(1+2*i), float(1+2*j))/float(2*subpixelCount) - vec2(0.5,0.5) )
                / screenResolution.x / numTiles.x;
            vec2 scaled_xy = tan(radians(fov * 0.5)) * (xy + offset);
            
            bool outsideView =
                perspectiveType == perspectiveTypeHyperideal &&
                length(scaled_xy) >= 0.5;

            if (outsideView) {
                total_color += vec3(1.0, 1.0, 1.0);
            } else {            
                RayHit ray_hit = computeRayHit(scaled_xy);
                if (ray_hit.object_type != object_type_nothing) {
                    min_depth = min(min_depth, tanh(ray_hit.dist));
                    if (isColored(ray_hit)) {
                        // Accumulate color for colored subpixels.
                        total_color += colorForRayHit(ray_hit);
                    } else {
                        // Accumulate value for valued subpixels.
                        // Count number of valued subpixels so that we can
                        // compute average value.
                        num_valued_subpixels += 1;
                        total_value += valueForRayHit(ray_hit);
                    }
                }
            }
        }
    }

    if (num_valued_subpixels > 0) {
        // Average value of valued ray hits.
        float value = total_value / num_valued_subpixels;

        // Compute gradient and add to total color.
        //
        // The contribution of the gradient color is proportional
        // to the number of valued subsamples for this pixel.
        total_color += num_valued_subpixels * colorForValue(value);
    }

    // Divide by total number of subsamples.
    out_FragColor = vec4(total_color / float(subpixelCount * subpixelCount), 1);
    gl_FragDepth = min_depth;
}
