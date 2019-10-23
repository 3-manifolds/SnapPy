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
uniform vec4 planes[4 * ##num_tets##];   // ##num_tets## gets replaced when we loadShaders in Main.js
uniform int otherTetNums[4 * ##num_tets##]; 
uniform int enteringFaceNums[4 * ##num_tets##]; 
uniform float weights[4 * ##num_tets##]; 
uniform mat4 SO13tsfms[4 * ##num_tets##];

uniform vec2 barycentricToMLCoordinates[3 * 4 * ##num_tets##];

uniform float gradientThreshholds[5];
uniform vec3 gradientColours[5];

uniform vec4 horospheres[4 * ##num_tets##];

uniform float edgeThicknessCylinder;
uniform mat4 SO13EdgeInvolutions[6 * ##num_tets##];

vec4 general_gradient(float t, float threshholds[5], vec3 colours[5]);

// inf             1
//   v0 -------- v2
//    | `.    .' |
//    |   `. '   |
//    |   . `.   |
//    | .'    `. |
//   v3 -------- v1
// z               0

  float R13_dot(vec4 u, vec4 v){
    return - u.x*v.x + u.y*v.y + u.z*v.z + u.w*v.w; // Lorentz Dot
  }

  float R13_norm_inv(vec4 v){
    return inversesqrt(abs(R13_dot(v,v)));
  }
  
  vec4 R13_normalise(vec4 v){
    return v*R13_norm_inv(v);
  }

float geodesicParameterPlanes(vec4 samplePoint, vec4 dualPoint1, vec4 dualPoint2){
  // "distance" from a geodesic defined by two (assumed perpendicular) geodesic planes, this is not quite distance, need to asinh(sqrt( result ))
  float dot1 = -R13_dot(samplePoint, dualPoint1);
  vec4 dualPointPerp = R13_normalise(dualPoint2 - R13_dot(dualPoint1, dualPoint2) * dualPoint1); // should be precalculated if this is a main feature
  float dot2 = -R13_dot(samplePoint, dualPointPerp);
  // float dot3 = -R13_dot(dualPoint1, dualPoint2);
  return dot1*dot1 + dot2*dot2;
  // return dot1*dot1 * (1.0 + dot3*dot3) + dot2*dot2 + 2.0*dot1*dot2*dot3;
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

  float hyp_dist(vec4 u, vec4 v){
    float bUV = -R13_dot(u,v);
    if (bUV < 1.0) {return 0.0;}
    else {return acosh(bUV);}  
  } 

float param_to_isect_line_with_horosphere(vec4 line_start, vec4 line_dir, vec4 horosphere)
{
    float start_dot = R13_dot(horosphere, line_start);
    float dir_dot = R13_dot(horosphere, line_dir);
    
    float a = dir_dot * dir_dot + 1;
    float b = 2 * dir_dot * start_dot;
    float c = start_dot * start_dot - 1;
    
    float disc = b * b - 4 * a * c;
    if (disc < 0) {
        return 200000000.0;
    }
    
    float result = (-b - sqrt(disc)) / (2 * a);
    if (result < 0) {
        return 200000000.0;
    }

    return result;
}

float param_to_isect_line_with_edge_cylinder(vec4 line_start, vec4 line_dir, mat4 involution)
{
    vec4 t_line_start = line_start * involution;
    vec4 t_line_dir   = line_dir   * involution;

    float a = R13_dot(t_line_dir,   line_dir)
            - edgeThicknessCylinder;
    float b = R13_dot(t_line_start, line_dir)
            + R13_dot(t_line_dir,   line_start);
    float c = R13_dot(t_line_start, line_start)
            + edgeThicknessCylinder;

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

float param_to_isect_line_with_plane(vec4 line_start, vec4 line_dir, vec4 plane){
    float denom = R13_dot(plane, line_dir);
    if(denom == 0.0){
        return 200000000.0;  // bigger than the initial smallest_p value we will accept
    }
    /// solve: R13_dot(plane, line_start + p * line_dir) = 0
    ///        R13_dot(plane, line_start) + p * R13_dot(plane, line_dir) = 0
    return (-R13_dot(plane, line_start)) / denom;
}

vec3 override_color = vec3(0);

const int objectTypeFace         = 1;
const int objectTypeEdgeCylinder = 2;
const int objectTypeHorosphere   = 3;
const int objectTypeEdgeFan      = 4;

uniform int c;

struct RayHit
{
    int tetNum;
    int objectType;
    int objectIndex;

    vec4 hitPoint;

    float weight;
    float dist;

    mat4 trans;
};

RayHit
ray_trace_through_hyperboloid_tet(vec4 init_pos, vec4 init_dir, int tetNum, int entryObjectType, int entryObjectIndex)
{
    RayHit rayHit;
    rayHit.tetNum = tetNum;

    ///Given shape of a tet and a ray, find where the ray exits and through which face
    float smallest_p = 100000000.0;

    for(int face = 0; face < 4; face++) {
        if (entryObjectType != objectTypeFace || entryObjectIndex != face) {
            // find p when we hit that face
            int index = 4 * tetNum + face;
            if(R13_dot(init_dir, planes[index]) > 0.0){ 
                float p = param_to_isect_line_with_plane(init_pos, init_dir, planes[index]);
                // if ((-10000.0 <= p) && (p < smallest_p)) {
                if (p < smallest_p) {  
                    /// negative values are ok if we have to go backwards a little to get through the face we are a little the wrong side of
                    /// Although this can apparently get caught in infinite loops in an edge

                    /// if we are on an edge then we don't in fact move as we go through this tet: t = 0.0
                    /// also allow tiny negative values, which will come up from floating point errors. 
                    /// surface normals check should ensure that even in this case we make progress through 
                    /// the triangles around an edge
                    smallest_p = p;
                    rayHit.objectType = objectTypeFace;
                    rayHit.objectIndex = face;
                }
            }
        }
    }

    for (int vertex = 0; vertex < 4; vertex++) {
        if (entryObjectType != objectTypeHorosphere || entryObjectIndex != vertex) {
            int index = 4 * tetNum + vertex;
            float p = param_to_isect_line_with_horosphere(init_pos, init_dir, horospheres[index]);
            if (p < smallest_p) {
                smallest_p = p;
                rayHit.objectType = objectTypeHorosphere;
                rayHit.objectIndex = vertex;
            }
        }
    }
                
    for (int edge = 0; edge < 6; edge++) {
        if (entryObjectType != objectTypeEdgeCylinder || entryObjectIndex != edge) {
            float p = param_to_isect_line_with_edge_cylinder(
                init_pos, init_dir, SO13EdgeInvolutions[6 * tetNum + edge]);
            if (p < smallest_p) {
                smallest_p = p;
                rayHit.objectType = objectTypeEdgeCylinder;
                rayHit.objectIndex = edge;
            }
        }
    }
    
    rayHit.hitPoint = R13_normalise( init_pos + smallest_p * init_dir );

    if (rayHit.objectType == objectTypeFace) {
        if(edgeThickness > 0.00001) {
            if(triangleBdryParam(rayHit.hitPoint, tetNum, rayHit.objectIndex) < edgeThickness) {
                rayHit.objectType = objectTypeEdgeFan;
            }
        }

    }

    return rayHit;
}

RayHit
ray_trace(vec4 init_pt, vec4 init_dir, int tetNum, float weight, mat4 currentTransform){
    int entryObjectType  = -1;
    int entryObjectIndex = -1;
    int index;
    mat4 tsfm;
    vec4 new_dir;

    float dist = 0;

    RayHit rayHit;

    for(int i = 0; i < maxSteps; i++){
        rayHit = ray_trace_through_hyperboloid_tet(
            init_pt, init_dir, tetNum, entryObjectType, entryObjectIndex);

        rayHit.trans = currentTransform;
        dist += hyp_dist(init_pt, rayHit.hitPoint);
        rayHit.dist = dist;

        if (rayHit.objectType == objectTypeHorosphere) {
            break;
        }

        if (rayHit.objectType == objectTypeEdgeCylinder) {
            break;
        }

        if (rayHit.objectType == objectTypeEdgeFan) {
            break;
        }

        if (dist > maxDist) {
            break;
        }

        // in fact pow(sinh(radius in hyperbolic units),2.0). However, sinh^2 is monotonic for 
        // positive values so we get correct behaviour by comparing without the sinh^2. 
        index = 4 * tetNum + rayHit.objectIndex;
        weight += weights[ index ];

        rayHit.weight = weight;
        rayHit.dist = dist;

        entryObjectIndex = enteringFaceNums[ index ];
        entryObjectType = objectTypeFace;
        tsfm = SO13tsfms[ index ];
        tetNum = otherTetNums[ index ];

        currentTransform = currentTransform * tsfm;

        new_dir = init_dir + R13_dot(init_dir, rayHit.hitPoint) * rayHit.hitPoint; // orthonormal decomp, no normalisation yet
        init_pt = rayHit.hitPoint * tsfm;
        init_dir = R13_normalise( new_dir * tsfm ); 
    }

    return rayHit;
}

vec3 shade(vec4 init_pt, vec4 init_dir, RayHit rayHit)
{
    if (rayHit.objectType == objectTypeHorosphere) {
        vec2 coords = vec2(0);

        for (int v1 = 0; v1 < 3; v1++) {
            int face = (rayHit.objectIndex + v1 + 1) % 4;

            coords +=
                barycentricToMLCoordinates[12 * rayHit.tetNum + 3 * rayHit.objectIndex + v1]
                * abs(R13_dot(rayHit.hitPoint, planes[4 * rayHit.tetNum + face]));
        }
        
        coords = fract(coords);

        vec3 color = vec3(0.3, 0.5, 0.4);
        if (coords.x < 0.03) {
            color.x = 1.0;
        }
        if (coords.y < 0.03) {
            color.y = 1.0;
        }

        vec4 rayEnd = R13_normalise(
            cosh(rayHit.dist) * init_pt + sinh(rayHit.dist) * init_dir);

        vec4 normal =
//            horospheres[4 * rayHit.tetNum + rayHit.objectIndex] * inverse(rayHit.trans);
            - rayEnd;

        color = 0.5 + 0.1 * normal.yzw;
        
        return color;
    }

    if (rayHit.objectType == objectTypeEdgeCylinder) {
        return vec3(0.5, 0.3, 0.7);
    }

    float weight;

    if (viewMode == 0) {
        weight = rayHit.weight;
    } else if (viewMode == 1) {
        weight = 0.5 * rayHit.dist;
    } else {
        weight = float(rayHit.tetNum);
    }

    weight = contrast * weight;
    weight = 0.5 + 0.5*weight/(abs(weight) + 1.0);  //faster than atan, similar

    return general_gradient(weight, gradientThreshholds, gradientColours).xyz;
}

/// --- Graph-trace code --- ///

float amountOutsideTetrahedron(vec4 v, int tetNum, out int biggest_face) {
  float biggest_amount = -100000.0;
  float amount;
  for(int i = 0; i < 4; i++){
    amount = R13_dot( v, planes[4*tetNum + i] );
    if( amount > biggest_amount ){
      biggest_amount = amount;
      biggest_face = i;
    }
  }
  return biggest_amount; 
}

// Get point at distance dist on the geodesic from u in the direction vPrime
vec4 pointOnGeodesic(vec4 u, vec4 vPrime, float dist){
  return u*cosh(dist) + vPrime*sinh(dist);
}

float graph_trace(inout vec4 goal_pt, inout int tetNum, out mat4 tsfm){ // tsfm is matrix to send goal_pt to its image in the tetrahedron coordinates it is in
  // similar function to ray_trace, but different algorithm
  float total_face_weight = currentWeight;
  int entry_face = -1;
  int index;
  int biggest_face;
  tsfm = mat4(1.0);
  for(int i=0; i<maxSteps; i++){
      if ( amountOutsideTetrahedron(goal_pt, tetNum, biggest_face) > 0.0000001 && biggest_face != entry_face ){
        index = 4*tetNum + biggest_face;
        entry_face = enteringFaceNums[ index ];
        tetNum = otherTetNums[ index ];
        total_face_weight += weights[ index ];
        goal_pt *= SO13tsfms[ index ];
        tsfm *= SO13tsfms[ index ];
        // if (R13_dot(goal_pt, goal_pt) > -0.5){ return -1000.0; } // errors accumulate and we get junk!
      }
      else{ break; }
    }
    return total_face_weight;
  }

/// --- Colour gradient code --- ///

int find_band(float t, float threshholds[5]){
    for(int j=1;j<4;j++){
        if(t < threshholds[j]){return j;}
    }
    return 4;
}
vec4 general_gradient(float t, float threshholds[5], vec3 colours[5]){
    int i = find_band(t, threshholds);
    return vec4( mix(colours[i-1], colours[i],(t - threshholds[i-1])/(threshholds[i] - threshholds[i-1]) ), 1.0);
}

/// --- Ray init pt and directions code --- ///

vec4 get_ray_dir(vec2 xy){ 
    xy = 0.2 * xy;
    float z = 0.1/tan(radians(fov*0.5));
    vec4 p =  R13_normalise(vec4(0.0, xy,-z));
    return p;
}

vec3 get_signed_count(vec2 xy){
    vec4 camera_init_pt;
    vec4 camera_init_dir;
    vec4 init_pt;
    vec4 init_dir;
    mat4 init_transform;
    vec3 color = vec3(0);

    if(perspectiveType == 0){ // material
        init_transform = currentBoost;
        camera_init_pt = vec4(1.0,0.0,0.0,0.0);
        camera_init_dir = get_ray_dir(xy);
        init_pt = camera_init_pt * currentBoost;
        init_dir = camera_init_dir * currentBoost;
        color = shade(camera_init_pt, camera_init_dir, ray_trace(init_pt, init_dir, tetNum, currentWeight, init_transform));
    }
    else { // ideal
        float foo = 0.5*dot(xy, xy);
        init_transform = currentBoost;
        camera_init_pt = vec4(foo + 1.0, xy, foo);   // parabolic transformation magic by Saul
        camera_init_dir = vec4(foo, xy, foo - 1.0);
        init_pt = camera_init_pt * currentBoost;
        init_dir = camera_init_dir * currentBoost;
        mat4 tsfm = mat4(1.0);
        int currentTetNum = tetNum;  // gets modified inside graph_trace
        
        float weight = graph_trace(init_pt, currentTetNum, tsfm);  // get us to the tetrahedron containing init_pt
        // init_pt *= tsfm;  // the point gets moved back in graph_trace
        init_transform = init_transform * tsfm;
        init_dir *= tsfm;  // move the direction back to here
//    if(viewMode != 0){ weight = 0.0; } // the other modes don't have a cumulative count from moving to the correct tetrahedron
        color = shade(camera_init_pt, camera_init_dir, ray_trace(init_pt, init_dir, currentTetNum, weight, init_transform));
    }
    
  /*
  if(viewMode == 0){ // Cannon-Thurston
    weight += currentWeight;
  }
  */
  return color;
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
      total_color += get_signed_count(xy + offset);
    }
  }
  vec3 color = total_color/float(subpixelCount*subpixelCount); // average over all subpixels

  out_FragColor = vec4(color, 1);

  /*
  weight = contrast * weight;
  weight = 0.5 + 0.5*weight/(abs(weight) + 1.0);  //faster than atan, similar
  // weight = 0.5 + atan(0.3 * weight)/PI;  // between 0.0 and 1.0
  out_FragColor = general_gradient(weight, gradientThreshholds, gradientColours);

  if (override_color != vec3(0)) {
      out_FragColor = vec4(override_color,1);
  }
  */

}
