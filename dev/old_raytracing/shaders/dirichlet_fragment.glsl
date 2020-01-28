#version 150

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
//

uniform float gradientThreshholds[5];
uniform vec3 gradientColours[5];

uniform float edgeThicknessCylinder;
uniform float sphereRadius;

uniform vec4 planes[##num_planes##];
uniform float weights[##num_planes##];
uniform vec4 vertex_positions[##num_vertices##];
uniform mat4 edge_involutions[##num_edges##];
uniform mat4 face_pairings[##num_planes##];

uniform int edge_color_indices[##num_edges##];
uniform int vertex_color_indices[##num_vertices##];

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

const int object_type_face          = 1;
const int object_type_edge_cylinder = 2;
const int object_type_horosphere    = 3;
const int object_type_edge_fan      = 4;
const int object_type_sphere        = 5;
const int object_type_nothing       = 6;

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

float param_to_isect_line_with_edge_cylinder(Ray ray, mat4 involution)
{
    vec4 image_start = ray.point * involution;
    vec4 image_dir   = ray.dir   * involution;

    float a = R13_dot(image_dir,   ray.dir)   - edgeThicknessCylinder;
    float b = R13_dot(image_start, ray.dir)   + R13_dot(image_dir,   ray.point);
    float c = R13_dot(image_start, ray.point) + edgeThicknessCylinder;

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
ray_trace_through_domain(inout RayHit ray_hit)
{
    int entry_object_type = ray_hit.object_type;
    int entry_object_index = ray_hit.object_index;
 
    ray_hit.object_type = object_type_nothing;

    ///Given shape of a tet and a ray, find where the ray exits and through which face
    float smallest_p = 100000000.0;

    for (int face = 0; face < ##num_planes##; face++) {
        float p = param_to_isect_line_with_plane(ray_hit.ray, planes[face]);
         
        if(R13_dot(ray_hit.ray.dir, planes[face]) > 0.0) {
            if (p < smallest_p) {
                smallest_p = p;
                ray_hit.object_type = object_type_face;
                ray_hit.object_index = face;
            }
        }
    }

    for (int vertex = 0; vertex < ##num_vertices##; vertex++) {
        float p = param_to_isect_line_with_sphere(
            ray_hit.ray,
            vertex_positions[vertex],
            sphereRadius);
        if (p < smallest_p) {
            smallest_p = p;
            ray_hit.object_type = object_type_sphere;
            ray_hit.object_index = vertex;
        }
    }

    for (int edge = 0; edge < ##num_edges##; edge++) {
        float p = param_to_isect_line_with_edge_cylinder(
            ray_hit.ray,
            edge_involutions[edge]);
        if (p < smallest_p) {
            smallest_p = p;
            ray_hit.object_type = object_type_edge_cylinder;
            ray_hit.object_index = edge;
        }
    }

    if (ray_hit.object_type == object_type_nothing) {
        return;
    }
    
    ray_hit.dist += atanh(smallest_p);
    ray_hit.ray.point = R13_normalise(
        ray_hit.ray.point + smallest_p * ray_hit.ray.dir );
    ray_hit.ray.dir = R13_normalise(
        R13_ortho_decomposition_time(ray_hit.ray.dir, ray_hit.ray.point));
}

RayHit
ray_trace(RayHit ray_hit) {

    for(int i = 0; i < maxSteps; i++){
        ray_trace_through_domain(ray_hit);

        if (ray_hit.object_type != object_type_face) {
            break;
        }

        if (ray_hit.dist > maxDist) {
            break;
        }

        ray_hit.weight += weights[ray_hit.object_index];

        ray_hit.object_index = (ray_hit.object_index & 0xfe) | ((~ray_hit.object_index) & 0x01);
        mat4 tsfm = face_pairings[ray_hit.object_index];
        
        ray_hit.eye_space_to_tet_space = ray_hit.eye_space_to_tet_space * tsfm;
        ray_hit.ray.point = ray_hit.ray.point * tsfm;
        ray_hit.ray.dir = R13_normalise( ray_hit.ray.dir * tsfm ); 

    }

    return ray_hit;
}

vec4 compute_normal(RayHit ray_hit)
{
    if(ray_hit.object_type == object_type_sphere) {
        vec4 diff = vertex_positions[ray_hit.object_index] - ray_hit.ray.point;
        return R13_normalise(
            R13_ortho_decomposition_time(
                diff, ray_hit.ray.point));
    }

    if(ray_hit.object_type == object_type_edge_cylinder) {
        vec4 image_pt = ray_hit.ray.point * edge_involutions[ray_hit.object_index];
        vec4 diff = image_pt - ray_hit.ray.point;
        return R13_normalise(
            R13_ortho_decomposition_time(
                diff, ray_hit.ray.point));
    }
    
    if(ray_hit.object_type == object_type_face) {
        return planes[ray_hit.object_index];
    }

    return vec4(0,1,0,0);
}

/// --- Colour gradient code --- ///

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

    if (ray_hit.object_type == object_type_sphere) {
        int index = vertex_color_indices[ray_hit.object_index];
        result.diffuse =
            vec3(0.5, 0.5, 0.5)
            + sin(index + 0.3) * vec3( 0.3,  -0.3,   0.0)
            + cos(index + 0.3) * vec3(0.15,   0.15, -0.3);

        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_edge_cylinder) {
        int index = edge_color_indices[ray_hit.object_index];
        result.diffuse =
            vec3(0.5, 0.5, 0.5)
            + sin(index + 1.3) * vec3( 0.3,  -0.3,   0.0)
            + cos(index + 1.3) * vec3(0.15,   0.15, -0.3);

        result.ambient = 0.5 * result.diffuse;
    }

    if (ray_hit.object_type == object_type_face) {
        result.diffuse =
            vec3(0.3, 0.3, 0.3)
            + sin(ray_hit.object_index + 1.1) * vec3( 0.3,  -0.3,   0.0)
            + cos(ray_hit.object_index + 1.1) * vec3(0.15,   0.15, -0.3);

        result.ambient = 0.5 * result.diffuse;
    }

    return result;
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

vec3 shade_with_lighting(RayHit ray_hit)
{
    MaterialParams material = material_params(ray_hit);

    vec4 normal = R13_normalise(compute_normal(ray_hit));

    vec4 light_position = R13_normalise(
        vec4(1,0,0.7,0) * ray_hit.eye_space_to_tet_space);

    vec4 light_dir_at_hit = R13_normalise(
            R13_ortho_decomposition_time(ray_hit.ray.point - light_position,
                                         ray_hit.ray.point));
    
    float normal_light = max(0, R13_dot(normal, light_dir_at_hit));

    vec4 half_angle = R13_normalise(light_dir_at_hit + ray_hit.ray.dir);

    float blinn_term =
        normal_light > 0.0
        ? pow(max(0, R13_dot(half_angle, normal)), material.shininess)
        : 0.0;

    return  material.ambient
          + material.diffuse * normal_light
          + material.specular * blinn_term;
}

vec3 shade(RayHit ray_hit)
{
    if (ray_hit.object_type == object_type_face) {
        return shade_by_gradient(ray_hit);
    }

    return shade_with_lighting(ray_hit);
}

/// --- Ray init pt and directions code --- ///

Ray get_ray_eye_space(vec2 xy)
{
    Ray result;

    // material
    result.point = vec4(1.0,0.0,0.0,0.0);
    float z = 0.5 / tan(radians(fov * 0.5));
    result.dir = R13_normalise(vec4(0.0, xy, -z));

    return result;
}

vec3 get_color(vec2 xy){
    Ray ray_eye_space = get_ray_eye_space(xy);

    RayHit ray_tet_space;
    ray_tet_space.ray.point = ray_eye_space.point * currentBoost;
    ray_tet_space.ray.dir   = ray_eye_space.dir   * currentBoost;
    ray_tet_space.dist = 0.0;
    ray_tet_space.weight = currentWeight;
    ray_tet_space.tet_num = tetNum;
    ray_tet_space.eye_space_to_tet_space = currentBoost;
    ray_tet_space.object_type = -1;
    ray_tet_space.object_index = -1;

    return shade(ray_trace(ray_tet_space));
}

void main(){
    vec2 xy = (gl_FragCoord.xy - 0.5*screenResolution.xy)/screenResolution.x;
    if(multiScreenShot == 1){  // Return multiple 4096x4096 screenshots that can be combined in, e.g. Photoshop.
        // Here screenResolution is really tileResolution;
        xy = (xy + tile - 0.5*(numTiles - vec2(1.0,1.0))) / numTiles.x;
    }
    vec3 total_color = vec3(0);
    for(int i=0; i<subpixelCount; i++){
        for(int j=0; j<subpixelCount; j++){
            vec2 offset = ( (float(1+2*i), float(1+2*j))/float(2*subpixelCount) - vec2(0.5,0.5) ) / screenResolution.x;
            total_color += get_color(xy + offset);
        }
    }
    vec3 color = total_color/float(subpixelCount*subpixelCount); // average over all subpixels
    
    out_FragColor = vec4(color, 1);
}
