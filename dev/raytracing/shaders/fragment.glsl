#version 150

// GLSL version 1.50 corresponds to OpenGL 3.2 which is the minimal version we
// require.

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

uniform int tetNum;
uniform mat4 currentBoost;
uniform float currentWeight;

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

// Convention: ##NAME## names a compile time constant.
// The string ##NAME## is replaced by the code in __init__.py
// by looking up its value in the dictionary constants_dict.

uniform vec4 planes[4 * ##num_tets##];
uniform int otherTetNums[4 * ##num_tets##]; 
uniform int enteringFaceNums[4 * ##num_tets##]; 
uniform float weights[4 * ##num_tets##]; 
uniform mat4 SO13tsfms[4 * ##num_tets##];
uniform mat4 SO13EdgeInvolutions[3 * ##num_tets##];

// For an incomplete cusp, the Margulis tube is a cylinder about
// the geodesic the triangulation spins about. In other words,
// it is the cylinder fixed by the peripheral group of that cusp.
// We encode this cylinder by the involution fixing it (i.e., the
// rotation about the geodesic by pi) and the cosh of its radius.
uniform mat4 margulisTubeSO13Involutions[4 * ##num_tets##];
uniform float margulisTubeCoshThickness[4 * ##num_tets##];

// Light vector corresponding to a cusp neighborhood.
uniform vec4 horospheres[4 * ##num_tets##];

// Heights of the Euclidean triangle obtained when intersecting
// horosphere with tetrahedron.
// There are four horospheres about the four edges of the tetrahedron,
// giving four triangles but they are all similar, so storing only
// one per tet.
uniform vec3 triangleHeightVectors[##num_tets##];

uniform float insphere_radii[##num_tets##];

uniform mat3x2 planeDistToComplexCoordinates[4 * ##num_tets##];
uniform vec2 logAdjustments[4 * ##num_tets##];
uniform mat2 matLogs[4 * ##num_tets##];

uniform float edgeThicknessCylinder;

uniform float gradientThreshholds[5];
uniform vec3 gradientColours[5];

uniform int edge_color_indices[6 * ##num_tets##];
uniform int horosphere_color_indices[4 * ##num_tets##];

uniform float lightBias;
uniform float lightFalloff;
uniform float brightness;

const vec4[6] edgeSym = vec4[](
    vec4( 1, 1, 1, 1),
    vec4( 1, 1, 1, 1),
    vec4( 1, 1, 1, 1),
    vec4( 1,-1,-1, 1),
    vec4( 1, 1,-1,-1),
    vec4( 1, 1,-1,-1));

int edgeIndex(int tet, int e)
{
    if (e < 3) {
        return 3 * tet + e;
    }
    return 3 * tet + 5 - e;
}

float R13_dot(vec4 u, vec4 v)
{
    return - u.x*v.x + u.y*v.y + u.z*v.z + u.w*v.w; // Lorentz Dot
}

float R13_norm_inv(vec4 v)
{
    return inversesqrt(abs(R13_dot(v,v)));
}
  
vec4 R13_normalise(vec4 v)
{
    return v * R13_norm_inv(v);
}

vec4 R13_ortho_decomposition_time(vec4 v, vec4 time_vector)
{
    return v + R13_dot(v, time_vector) * time_vector;
}

float geodesicParameterPlanes(vec4 samplePoint, vec4 dualPoint1, vec4 dualPoint2){
  // "distance" from a geodesic defined by two (assumed perpendicular) geodesic planes, this is not quite distance, need to asinh(sqrt( result ))

  float dot1 = -R13_dot(samplePoint, dualPoint1);
  vec4 dualPointPerp = R13_normalise(dualPoint2 - R13_dot(dualPoint1, dualPoint2) * dualPoint1); // should be precalculated if this is a main feature
  float dot2 = -R13_dot(samplePoint, dualPointPerp);

  return dot1*dot1 + dot2*dot2;
}

float triangleBdryParam(vec4 samplePoint, int tetNum, int exit_face){
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

float hyp_dist(vec4 u, vec4 v) {
    float bUV = -R13_dot(u,v);
    if (bUV < 1.0) {
        return 0.0;
    }
    else {
        return acosh(bUV);
    }
} 

const int object_type_nothing       = 0;
const int object_type_face          = 1;
const int object_type_edge_cylinder = 2;
const int object_type_horosphere    = 3;
const int object_type_edge_fan      = 4;
const int object_type_sphere        = 5;
const int object_type_margulis_tube = 6;

struct Ray
{
    vec4 point;
    vec4 dir;
};

struct RayHit
{
    Ray ray;
    int tet_num;
    mat4 eye_space_to_tet_space;
    float dist;
    float weight;
    int object_type;
    int object_index;
};

float param_to_isect_line_with_sphere(Ray ray, vec4 center, float cosh_sqr_radius)
{
    float start_dot = R13_dot(center, ray.point);
    float dir_dot   = R13_dot(center, ray.dir);
    
    float a = dir_dot * dir_dot + cosh_sqr_radius;
    float b = 2 * dir_dot * start_dot;
    float c = start_dot * start_dot - cosh_sqr_radius;

    float disc = b * b - 4 * a * c;
    if (disc < 0) {
        return 200000000.0;
    }
    
    float result = (-b - sign(a) * sqrt(disc)) / (2 * a);
    if (result < 0) {
        return 200000000.0;
    }

    return result;
}

float param_to_isect_line_with_horosphere(Ray ray, vec4 horosphere)
{
    return param_to_isect_line_with_sphere(ray, horosphere, 1);
}

float param_to_isect_line_with_edge_cylinder(Ray ray, mat4 involution, float radius, float back_p, vec4 sym)
{
    vec4 pt  = ray.point * sym;
    vec4 dir = ray.dir   * sym;

    vec4 image_start = pt  * involution;
    vec4 image_dir   = dir   * involution;

    float a = R13_dot(image_dir,   dir)   - radius;
    float b = R13_dot(image_start, dir)   + R13_dot(image_dir,   pt );
    float c = R13_dot(image_start, pt ) + radius;

    float disc = b * b - 4 * a * c;
    if (disc < 0) {
        return 200000000.0;
    }
    
    float result = (-b - sign(a) * sqrt(disc)) / (2 * a);
    if (result < back_p) {
        return 200000000.0;
    }

    return result;
}

float param_to_isect_line_with_plane(Ray ray, vec4 plane) {
    float denom = R13_dot(plane, ray.dir);
    if(denom == 0.0){
        return 200000000.0;  // bigger than the initial smallest_p value we will accept
    }
    /// solve: R13_dot(plane, line_start + p * line_dir) = 0
    ///        R13_dot(plane, line_start) + p * R13_dot(plane, line_dir) = 0
    return (-R13_dot(plane, ray.point)) / denom;
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
            if(R13_dot(ray_hit.ray.dir, planes[index]) > 0.0){ 
                float p = param_to_isect_line_with_plane(ray_hit.ray, planes[index]);
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
        float p = param_to_isect_line_with_sphere(
            ray_hit.ray,
            vec4(1,0,0,0),
            insphere_radii[ray_hit.tet_num]);
        if (p < smallest_p) {
            smallest_p = p;
            ray_hit.object_type = object_type_sphere;
            ray_hit.object_index = 0;
        }
    }

    for (int vertex = 0; vertex < 4; vertex++) {
        int index = 4 * ray_hit.tet_num + vertex;
        if (horospheres[index] != vec4(0)) {
            float p = param_to_isect_line_with_horosphere(ray_hit.ray, horospheres[index]);
            if (p < smallest_p) {
                smallest_p = p;
                ray_hit.object_type = object_type_horosphere;
                ray_hit.object_index = vertex;
            }
        }
    }
                
    float back_p = tanh(-ray_hit.dist);

    if (edgeThicknessCylinder > 1.00002) {
        for (int edge = 0; edge < 6; edge++) {
            float p = param_to_isect_line_with_edge_cylinder(
                ray_hit.ray,
                SO13EdgeInvolutions[edgeIndex(ray_hit.tet_num,edge)],
                edgeThicknessCylinder,
                back_p,
                edgeSym[edge]);
            if (p < smallest_p) {
                smallest_p = p;
                ray_hit.object_type = object_type_edge_cylinder;
                ray_hit.object_index = edge;
            }
        }
    }

    for (int vertex = 0; vertex < 4; vertex++) {
        int index = 4 * ray_hit.tet_num + vertex;

        if (margulisTubeCoshThickness[index] > 1.00002) {
            float p = param_to_isect_line_with_edge_cylinder(
                ray_hit.ray,
                margulisTubeSO13Involutions[index],
                margulisTubeCoshThickness[index],
                back_p,
                vec4(1));
            if (p < smallest_p) {
                smallest_p = p;
                ray_hit.object_type = object_type_margulis_tube;
                ray_hit.object_index = vertex;
            }
        }
    }

    ray_hit.dist += atanh(smallest_p);
    ray_hit.ray.point = R13_normalise(
        ray_hit.ray.point + smallest_p * ray_hit.ray.dir );
    ray_hit.ray.dir = R13_normalise(
        R13_ortho_decomposition_time(ray_hit.ray.dir, ray_hit.ray.point));

    if(edgeThickness > 0.00001) {
        if (ray_hit.object_type == object_type_face) {
            if(triangleBdryParam(ray_hit.ray.point, ray_hit.tet_num, ray_hit.object_index) < edgeThickness) {
                ray_hit.object_type = object_type_edge_fan;
            }
        }
    }
}

RayHit
ray_trace(RayHit ray_hit) {

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

        ray_hit.eye_space_to_tet_space = ray_hit.eye_space_to_tet_space * tsfm;
        ray_hit.ray.point = ray_hit.ray.point * tsfm;
        ray_hit.ray.dir = R13_normalise( ray_hit.ray.dir * tsfm ); 
        ray_hit.tet_num = otherTetNums[ index ];
    }

    return ray_hit;
}

float compute_cusp_triangle_coordinate(RayHit ray_hit, int v1)
{
    int face = (ray_hit.object_index + v1 + 1) % 4;
    int plane_index = 4 * ray_hit.tet_num + face;
    return R13_dot(ray_hit.ray.point, planes[plane_index]);
}

vec2 compute_cusp_triangle_position(RayHit ray_hit)
{
    vec3 coordinates =
        vec3(compute_cusp_triangle_coordinate(ray_hit, 0),
             compute_cusp_triangle_coordinate(ray_hit, 1),
             compute_cusp_triangle_coordinate(ray_hit, 2));

    vec3 heights = triangleHeightVectors[ray_hit.tet_num];

    coordinates /= (ray_hit.object_index % 2 == 0) ? heights : heights.zyx;
    coordinates /= coordinates.x + coordinates.y + coordinates.z;

    int index = 4 * ray_hit.tet_num + ray_hit.object_index;

    return planeDistToComplexCoordinates[index] * coordinates;
}

vec2 compute_ML_coordinates_for_horosphere(RayHit ray_hit)
{
    return fract(compute_cusp_triangle_position(ray_hit));
}

vec2 compute_ML_coordinates_for_banana(RayHit ray_hit)
{
    vec2 result = compute_cusp_triangle_position(ray_hit);

    int index = 4 * ray_hit.tet_num + ray_hit.object_index;

    vec2 log_result =
        vec2(log(length(result)),
             atan(result.y, result.x)) +
        logAdjustments[index];
    
    return fract(log_result * matLogs[index]);
}

vec4 compute_normal(RayHit ray_hit)
{
    if(ray_hit.object_type == object_type_horosphere) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        return horospheres[index] - ray_hit.ray.point;
    }

    if(ray_hit.object_type == object_type_sphere) {
        vec4 diff = vec4(1,0,0,0) - ray_hit.ray.point;
        return R13_normalise(
            R13_ortho_decomposition_time(
                diff, ray_hit.ray.point));
    }
    
    if(ray_hit.object_type == object_type_edge_fan) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        return planes[index];
    }

    if(ray_hit.object_type == object_type_edge_cylinder) {
        int index = edgeIndex(ray_hit.tet_num, ray_hit.object_index);
        vec4 image_pt =
            edgeSym[ray_hit.object_index] *
            ((edgeSym[ray_hit.object_index] * ray_hit.ray.point)
             * SO13EdgeInvolutions[index]);
        vec4 diff = image_pt - ray_hit.ray.point;
        return R13_normalise(
            R13_ortho_decomposition_time(
                diff, ray_hit.ray.point));
    }

    if(ray_hit.object_type == object_type_margulis_tube) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        vec4 image_pt = ray_hit.ray.point * margulisTubeSO13Involutions[index];
        vec4 diff = image_pt - ray_hit.ray.point;
        return R13_normalise(
            R13_ortho_decomposition_time(
                diff, ray_hit.ray.point));
    }

    return vec4(0,1,0,0);
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

    if (ray_hit.object_type == object_type_horosphere) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = horosphere_color_indices[index];

        result.diffuse =
            vec3(0.5, 0.5, 0.5)
            + sin(color_index) * vec3( 0.3,  -0.3,   0.0)
            + cos(color_index) * vec3(0.15,   0.15, -0.3);
        result.ambient = 0.5 * result.diffuse;

        vec2 coords = compute_ML_coordinates_for_horosphere(ray_hit);
        
        if (coords.x < 0.015 || coords.x > 0.985) {
            result.diffuse = vec3(1,0.2,0.2);
            result.ambient = result.diffuse;
        }
        if (coords.y < 0.015 || coords.y > 0.985) {
            result.diffuse = vec3(0.2,1.0,0.2);
            result.ambient = result.diffuse;
        }
    }

    if (ray_hit.object_type == object_type_margulis_tube) {
        int index = 4 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = horosphere_color_indices[index];
        
        result.diffuse =
            vec3(0.5, 0.5, 0.5)
            + sin(color_index+0.4) * vec3( 0.3,  -0.3,   0.0)
            + cos(color_index+0.4) * vec3(0.15,   0.15, -0.3);
        result.ambient = 0.5 * result.diffuse;

        vec2 coords = compute_ML_coordinates_for_banana(ray_hit);

        if (coords.x < 0.015 || coords.x > 0.985) {
            result.diffuse = vec3(1,0.2,0.2);
            result.ambient = result.diffuse;
        }
        if (coords.y < 0.015 || coords.y > 0.985) {
            result.diffuse = vec3(0.2,1.0,0.2);
            result.ambient = result.diffuse;
        }
    }

    if (ray_hit.object_type == object_type_edge_fan) {
        result.diffuse = vec3(0.6, 0.2, 0.2);
        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_edge_cylinder) {
        int index = 6 * ray_hit.tet_num + ray_hit.object_index;
        int color_index = edge_color_indices[index];
        
        result.diffuse =
            vec3(0.5, 0.5, 0.5)
            + sin(color_index) * vec3( 0.3,  -0.3,   0.0)
            + cos(color_index) * vec3(0.15,   0.15, -0.3);

        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_sphere) {
        result.diffuse =
            vec3(0.5, 0.5, 0.5)
            + sin(ray_hit.tet_num + 0.3) * vec3( 0.3,  -0.3,   0.0)
            + cos(ray_hit.tet_num + 0.3) * vec3(0.15,   0.15, -0.3);

        result.ambient = 0.5 * result.diffuse;
    }

    return result;
}

vec3 shade_with_lighting(RayHit ray_hit)
{
    MaterialParams material = material_params(ray_hit);

    vec4 normal = R13_normalise(compute_normal(ray_hit));

    vec4 light_position = R13_normalise(
        vec4(1,0,0.7,0) * ray_hit.eye_space_to_tet_space);

    float dist = acosh(-R13_dot(ray_hit.ray.point, light_position));

    vec4 light_dir_at_hit = R13_normalise(
            R13_ortho_decomposition_time(ray_hit.ray.point - light_position,
                                         ray_hit.ray.point));
    
    float normal_light = max(0, R13_dot(normal, light_dir_at_hit));

    vec4 half_angle = R13_normalise(light_dir_at_hit + ray_hit.ray.dir);

    float blinn_term =
        normal_light > 0.0
        ? pow(max(0, R13_dot(half_angle, normal)), material.shininess)
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
        ray_hit.object_type == object_type_edge_cylinder ||
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

float amountOutsideTetrahedron(vec4 v, int tet_num, out int biggest_face)
{
  float biggest_amount = -100000.0;
  float amount;
  for(int i = 0; i < 4; i++){
    amount = R13_dot( v, planes[4*tet_num + i] );
    if( amount > biggest_amount ){
      biggest_amount = amount;
      biggest_face = i;
    }
  }
  return biggest_amount; 
}

void graph_trace(inout RayHit ray)
{
  int entry_face = -1;
  mat4 tsfm = mat4(1.0);

  for(int i = 0; i < maxSteps; i++) {
      int biggest_face;
      if ( amountOutsideTetrahedron(ray.ray.point, ray.tet_num, biggest_face) > 0.0000001 && biggest_face != entry_face ){
          int index = 4 * ray.tet_num + biggest_face;
          entry_face = enteringFaceNums[ index ];
          ray.tet_num = otherTetNums[ index ];
          ray.weight += weights[ index ];
          ray.ray.point = ray.ray.point * SO13tsfms[ index ];
          tsfm = tsfm * SO13tsfms[ index ];
          // if (R13_dot(goal_pt, goal_pt) > -0.5){ return -1000.0; } // errors accumulate and we get junk!
      }
      else{
          break;
      }
  }

  ray.eye_space_to_tet_space = ray.eye_space_to_tet_space * tsfm;
  ray.ray.dir = ray.ray.dir * tsfm;  // move the direction back to here
}

/// --- Ray init pt and directions code --- ///

Ray get_ray_eye_space(vec2 xy)
{
    Ray result;

    if(perspectiveType == 0) {
        // material
        result.point = vec4(1.0,0.0,0.0,0.0);
        float z = 0.5 / tan(radians(fov * 0.5));
        result.dir = R13_normalise(vec4(0.0, xy, -z));
    } else {
        // ideal
        float foo = 0.5 * dot(xy, xy);
        // parabolic transformation magic by Saul
        result.point = vec4(foo + 1.0, xy, foo);
        result.dir   = vec4(foo,       xy, foo - 1.0);
    }

    return result;
}

vec4 get_color_and_depth(vec2 xy){
    Ray ray_eye_space = get_ray_eye_space(xy);

    RayHit ray_tet_space;
    ray_tet_space.ray.point = ray_eye_space.point * currentBoost;
    ray_tet_space.ray.dir   = ray_eye_space.dir   * currentBoost;
    ray_tet_space.dist = 0.0;
    ray_tet_space.weight = currentWeight;
    ray_tet_space.tet_num = tetNum;
    ray_tet_space.eye_space_to_tet_space = currentBoost;
    ray_tet_space.object_type = object_type_nothing;
    ray_tet_space.object_index = -1;

    if (perspectiveType == 1) {
        graph_trace(ray_tet_space);
    }
    
    return shade(ray_trace(ray_tet_space));
}

void main(){
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
