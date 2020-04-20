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
uniform bool perspectiveType;
uniform int viewMode;
uniform int multiScreenShot;
uniform vec2 tile;
uniform vec2 numTiles;

// Convention: ##NAME## names a compile time constant.
// The string ##NAME## is replaced by the code in __init__.py
// by looking up its value in the dictionary constants_dict.

uniform int otherTetNums[4 * ##num_tets##]; 
uniform int enteringFaceNums[4 * ##num_tets##]; 
uniform float weights[4 * ##num_tets##]; 

// +1 or -1 depending on orientation of tetrahedron.
uniform int orientations[##num_tets##];

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

layout (std140) uniform TetrahedraBasics
{
    vec4 R13Vertices[4 * ##num_tets##];
    vec4 planes[4 * ##num_tets##];
    mat4 SO13tsfms[4 * ##num_tets##];
};

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

uniform float edgeTubeRadiusParam;

uniform float gradientThreshholds[5];
uniform vec3 gradientColours[5];

uniform int face_color_indices[4 * ##num_tets##];
uniform int edge_color_indices[6 * ##num_tets##];
uniform int horosphere_color_indices[4 * ##num_tets##];

uniform float lightBias;
uniform float lightFalloff;
uniform float brightness;

uniform bool isNonGeometric;
uniform sampler2D nonGeometricTexture;

const int num_tets = ##num_tets##;
const int num_cusps = ##num_cusps##;

const float peripheralCurveThickness = 0.015;

const vec3 longitudeColor = vec3( 1.0 , 1.0 , 1.0 );
const vec3 meridianColor  = vec3( 0.5 , 0.5 , 0.5 );

const vec4 lightSourcePosition = vec4(1.0, 0.0, 0.7, 0.0);

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

const int object_type_nothing             = 0;
const int object_type_face                = 1;
const int object_type_edge_cylinder_enter = 2;
const int object_type_edge_cylinder_exit  = 3;
const int object_type_horosphere          = 4;
const int object_type_edge_fan            = 5;
const int object_type_sphere              = 6;
const int object_type_margulis_tube       = 7;

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


// Normal for a sphere about the origin.
vec4
normalForSphere(vec4 point)
{
    vec4 t = vec4(1,0,0,0) - point;

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

// The equation for the horosphere about a vertex in a tetrahedron.
// index is 4 * tetIndex + vertex.
vec4
horosphereEqn(int index)
{
    return horosphereScales[index] * R13Vertices[index];
}

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
    return vec4[](R13Vertices[4 * tet + edgeToVertices[edge].x],
                  R13Vertices[4 * tet + edgeToVertices[edge].y]);
}

// The two endpoints of a Margulis tube.
vec4[2]
endpointsForMargulisTube(int index)
{
    return vec4[](margulisTubeTails[index],
                  margulisTubeHeads[index]);
}

vec4
normalForRayHit(RayHit ray_hit)
{
    if(ray_hit.object_type == object_type_horosphere) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        return horosphereEqn(index) - ray_hit.ray.point;
    }

    if(ray_hit.object_type == object_type_sphere) {
        return normalForSphere(ray_hit.ray.point);
    }
    
    if(ray_hit.object_type == object_type_edge_fan) {
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

    if(ray_hit.object_type == object_type_margulis_tube) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        return normalForTube(
            ray_hit.ray.point,
            endpointsForMargulisTube(index));
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

// Compute the SO13 transform coresponding to the PGL(2,C)-matrix
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

    {
        float r = insphereRadiusParams[ray_hit.tet_num];

        if (r > 1.0001) {
            float p = distParamsForSphereIntersection(
                ray_hit.ray,
                vec4(1,0,0,0),
                r).x;
            if (p < smallest_p) {
                smallest_p = p;
                ray_hit.object_type = object_type_sphere;
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
    }
                
    float backDistParam = tanh(ray_hit.distWhenLeavingCusp-ray_hit.dist);

    if (edgeTubeRadiusParam > 0.50001) {
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

    for (int vertex = 0; vertex < 4; vertex++) {
        int index = 4 * ray_hit.tet_num + vertex;

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
        ray_hit.weight += weights[ index ];

        ray_hit.object_index = enteringFaceNums[ index ];
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

float value_for_gradient(RayHit ray_hit)
{
    if (viewMode == 0) {
        return ray_hit.weight;
    } else if (viewMode == 1) {
        return 0.5 * ray_hit.dist;
    } else {
        return float(ray_hit.tet_num);
    }
}

vec3 shade_by_gradient(RayHit ray_hit)
{
    float value = value_for_gradient(ray_hit);
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

    if (ray_hit.object_type == object_type_horosphere ||
        ray_hit.object_type == object_type_margulis_tube) {
        
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = horosphere_color_indices[index];

        result.diffuse = hsv2rgb(vec3(float(color_index)/float(num_cusps), 0.25, 1.0));
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
        result.diffuse = hsv2rgb(vec3(float(color_index)/float(num_tets), 1.0, 1.0));
        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_edge_cylinder_exit) {
        int index = 6 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = edge_color_indices[index];
        
        //using num_tets = num_edges
        result.diffuse = 0.3 * hsv2rgb(vec3(float(color_index)/float(num_tets), 1.0, 1.0));
        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_sphere) {
        result.diffuse = hsv2rgb(vec3(float(ray_hit.tet_num)/float(num_tets), 0.5, 1.0));
        result.ambient = 0.5 * result.diffuse;
    }

    return result;
}

vec3 shade_with_lighting(RayHit ray_hit)
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

vec4 shade(RayHit ray_hit)
{
    float depth = tanh(ray_hit.dist);

    if (ray_hit.object_type == object_type_sphere ||
        ray_hit.object_type == object_type_horosphere ||
        ray_hit.object_type == object_type_edge_cylinder_enter ||
        ray_hit.object_type == object_type_edge_cylinder_exit ||
        ray_hit.object_type == object_type_margulis_tube ||
        ray_hit.object_type == object_type_edge_fan) {
        
        return vec4(shade_with_lighting(ray_hit), depth);
    } else if (ray_hit.object_type == object_type_nothing) {
        return vec4(0,0,0,1);
    } else {
        return vec4(shade_by_gradient(ray_hit), depth);
    }
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
      entry_face = enteringFaceNums[ index ];
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

    float fov_scale = 2.0 * tan(radians(fov * 0.5));
    
    if (perspectiveType) {
        // ideal
        
        // factor 0.5 so that far away objects look the same
        // in ideal and material view
        vec2 scaled_xy = 0.5 * fov_scale * xy;

        float foo = 0.5 * dot(scaled_xy, scaled_xy);
        // parabolic transformation magic by Saul
        result.point = vec4(foo + 1.0, scaled_xy, foo);
        result.dir   = vec4(foo,       scaled_xy, foo - 1.0);
    } else {
        // material

        vec2 scaled_xy = fov_scale * xy;

        result.point = vec4(1.0,0.0,0.0,0.0);
        result.dir = R13Normalise(vec4(0.0, scaled_xy, -1.0));
    }
    
    return result;
}

// Determine whether the ray starts in a horosphere or
// Margulis tube.
// If yes, move the ray to the point where we exit the
// horosphere and set rayHit.object_type to horosphere.
//
// The result is true if we are inside a horosphere AND
// the ray is hitting a peripheral curve on the horosphere.
// 
// For optimization, leaveCusp will also apply a
// parabolic transformation to the ray trying to bring the
// where we exit the horosphere closer to teh entry point.
bool
leaveCusp(inout RayHit rayHit)
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

        // If we are inside a horosphere, leaveCusp has computed
        // the point where we leave the horosphere. But that point
        // might not be inside the current tetrahedron, so fix it.
        graph_trace(rayHit);
    }

    return false;
}

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

    // If using "parabolic" camera where the ray's do not
    // all start from a common point, transform ray first
    // to be inside a tetrahedron.
    if (perspectiveType) {
        graph_trace(ray_tet_space);
    }

#if 1

    // Check whether we are in a horosphere and if yes,
    // whether the ray hit a peripheral curve.
    bool hitPeripheral = leaveCusp(ray_tet_space);

    if (hitPeripheral) {
        // If we hit a peripheral curve, leaveCusp has given us
        // the intersection point and we can immeadiately shade.
    } else {
        // In all other cases, we need to raytrace before we shade.
        ray_trace(ray_tet_space);
    }
#else
    ray_trace(ray_tet_space);
#endif
    
    return ray_tet_space;
}

vec4 get_color_and_depth(vec2 xy)
{
    return shade(computeRayHit(xy));
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
    vec3 total_color = vec3(0);
    float depth;
    for(int i=0; i<subpixelCount; i++){
        for(int j=0; j<subpixelCount; j++){
            vec2 offset = ( (float(1+2*i), float(1+2*j))/float(2*subpixelCount) - vec2(0.5,0.5) ) / screenResolution.x;
            vec4 color_and_depth = get_color_and_depth(xy + offset);
            gl_FragDepth = color_and_depth.w;
            total_color += color_and_depth.rgb;
        }
    }

    vec3 color = total_color/float(subpixelCount*subpixelCount); // average over all subpixels

    out_FragColor = vec4(color, 1);
}
