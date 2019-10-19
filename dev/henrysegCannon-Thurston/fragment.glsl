
// inf             1
//   v0 -------- v2
//    | `.    .' |
//    |   `. '   |
//    |   . `.   |
//    | .'    `. |
//   v3 -------- v1
// z               0

  float R31_dot(vec4 u, vec4 v){
    return u.x*v.x + u.y*v.y + u.z*v.z - u.w*v.w; // Lorentz Dot
  }

  float R31_norm_inv(vec4 v){
    return inversesqrt(abs(R31_dot(v,v)));
  }
  
  vec4 R31_normalise(vec4 v){
    return v*R31_norm_inv(v);
  }

float geodesicParameterPlanes(vec4 samplePoint, vec4 dualPoint1, vec4 dualPoint2){
  // "distance" from a geodesic defined by two (assumed perpendicular) geodesic planes, this is not quite distance, need to asinh(sqrt( result ))
  float dot1 = -R31_dot(samplePoint, dualPoint1);
  vec4 dualPointPerp = R31_normalise(dualPoint2 - R31_dot(dualPoint1, dualPoint2) * dualPoint1); // should be precalculated if this is a main feature
  float dot2 = -R31_dot(samplePoint, dualPointPerp);
  // float dot3 = -R31_dot(dualPoint1, dualPoint2);
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
    float bUV = -R31_dot(u,v);
    if (bUV < 1.0) {return 0.0;}
    else {return acosh(bUV);}  
  } 

float param_to_isect_line_with_plane(vec4 line_start, vec4 line_dir, vec4 plane){
    float denom = R31_dot(plane, line_dir);
    if(denom == 0.0){ return 200000000.0; }  // bigger than the initial smallest_p value we will accept
    /// solve: R31_dot(plane, line_start + p * line_dir) = 0
    ///        R31_dot(plane, line_start) + p * R31_dot(plane, line_dir) = 0
    return (-R31_dot(plane, line_start)) / denom;
  }

vec4 ray_trace_through_hyperboloid_tet(vec4 init_pos, vec4 init_dir, int tetNum, int entry_face, out int exit_face){
    ///Given shape of a tet and a ray, find where the ray exits and through which face
    float smallest_p = 100000000.0;
    for(int face=0; face<4; face++){
        if(face != entry_face){  // find p when we hit that face
            int index = 4*tetNum + face;
            if(R31_dot(init_dir, planes[index]) > 0.0){ 
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
                    exit_face = face;
                }
            }
        }
    }
    return R31_normalise( init_pos + smallest_p * init_dir );
}

float ray_trace(vec4 init_pt, vec4 init_dir, float dist_to_go, int tetNum){
    int entry_face = -1;   /// starts off with no entry face
    int exit_face = -1;
    float total_face_weight = 0.0;
    vec4 new_pt;
    int index;
    mat4 tsfm;
    vec4 new_dir;
    for(int i=0; i<maxSteps; i++){
      new_pt = ray_trace_through_hyperboloid_tet(init_pt, init_dir, tetNum, entry_face, exit_face);
      dist_to_go -= hyp_dist(init_pt, new_pt);
      if (dist_to_go <= 0.0){ break; }
      if(edgeThickness > 0.00001){
        if(triangleBdryParam(new_pt, tetNum, exit_face) < edgeThickness){ break; }}
        // in fact pow(sinh(radius in hyperbolic units),2.0). However, sinh^2 is monotonic for 
        // positive values so we get correct behaviour by comparing without the sinh^2. 
      index = 4*tetNum + exit_face;
      total_face_weight += weights[ index ];
      entry_face = entering_face_nums[ index ];
      tsfm = SO31tsfms[ index ];
      tetNum = otherTetNums[ index ];

      new_dir = init_dir + R31_dot(init_dir, new_pt) * new_pt; // orthonormal decomp, no normalisation yet
      init_pt = new_pt * tsfm;  
      init_dir = R31_normalise( new_dir * tsfm ); 
    }
    if(viewMode == 0){ return total_face_weight; } // Cannon-Thurston
    else if(viewMode == 1){ return 0.5*maxDist - dist_to_go; } // Distance
    else{ return float(tetNum);}
}

/// --- Graph-trace code --- ///

float amountOutsideTetrahedron(vec4 v, int tetNum, out int biggest_face) {
  float biggest_amount = -100000.0;
  float amount;
  for(int i = 0; i < 4; i++){
    amount = R31_dot( v, planes[4*tetNum + i] );
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
  float total_face_weight = 0.0;
  int entry_face = -1;
  int index;
  int biggest_face;
  tsfm = mat4(1.0);
  for(int i=0; i<maxSteps; i++){
      if ( amountOutsideTetrahedron(goal_pt, tetNum, biggest_face) > 0.0000001 && biggest_face != entry_face ){
        index = 4*tetNum + biggest_face;
        entry_face = entering_face_nums[ index ];
        tetNum = otherTetNums[ index ];
        total_face_weight += weights[ index ];
        goal_pt *= SO31tsfms[ index ];
        tsfm *= SO31tsfms[ index ];
        // if (R31_dot(goal_pt, goal_pt) > -0.5){ return -1000.0; } // errors accumulate and we get junk!
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
    vec4 p =  R31_normalise(vec4(xy,-z,0.0));
    return p;
}

float get_signed_count(vec2 xy){
  vec4 init_pt;
  vec4 init_dir;
  float weight = 0.0;
  if(perspectiveType == 0){ // material
    init_pt = vec4(0.0,0.0,0.0,1.0);
    init_dir = get_ray_dir(xy);
    init_pt *= currentBoost;
    init_dir *= currentBoost; 
    weight = ray_trace(init_pt, init_dir, maxDist, tetNum);
  }
  else{ // ideal
    float foo = 0.5*dot(xy, xy);
    init_pt = vec4(xy.x, xy.y, foo, foo + 1.0);   // parabolic transformation magic by Saul
    init_dir = vec4(xy.x, xy.y, foo - 1.0, foo);
    init_pt *= currentBoost;
    init_dir *= currentBoost; 
    mat4 tsfm = mat4(1.0);
    int currentTetNum = tetNum;  // gets modified inside graph_trace

    weight = graph_trace(init_pt, currentTetNum, tsfm);  // get us to the tetrahedron containing init_pt
    // init_pt *= tsfm;  // the point gets moved back in graph_trace
    init_dir *= tsfm;  // move the direction back to here
    if(viewMode != 0){ weight = 0.0; } // the other modes don't have a cumulative count from moving to the correct tetrahedron
    weight += ray_trace(init_pt, init_dir, maxDist, currentTetNum);
  }

  if(viewMode == 0){ // Cannon-Thurston
    weight += currentWeight;
  }
  return weight;
}

void main(){
  vec2 xy = (gl_FragCoord.xy - 0.5*screenResolution.xy)/screenResolution.x;
  if(multiScreenShot == 1){  // Return multiple 4096x4096 screenshots that can be combined in, e.g. Photoshop.
    // Here screenResolution is really tileResolution;
    xy = (xy + tile - 0.5*(numTiles - vec2(1.0,1.0))) / numTiles.x;
  }
  float total_weight = 0.0;
  for(int i=0; i<subpixelCount; i++){
    for(int j=0; j<subpixelCount; j++){
      vec2 offset = ( (float(1+2*i), float(1+2*j))/float(2*subpixelCount) - vec2(0.5,0.5) ) / screenResolution.x;
      total_weight += get_signed_count(xy + offset);
    }
  }
  float weight = total_weight/float(subpixelCount*subpixelCount); // average over all subpixels

  weight = contrast * weight;
  weight = 0.5 + 0.5*weight/(abs(weight) + 1.0);  //faster than atan, similar
  // weight = 0.5 + atan(0.3 * weight)/PI;  // between 0.0 and 1.0
  out_FragColor = general_gradient(weight, gradientThreshholds, gradientColours);
}
