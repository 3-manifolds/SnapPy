
#include "CLinkProjectionPanorama.h"
#include "string.h"

/* provides the function:
 * 
 * Triangulation* triangulate_link_complement_from_file(char* file_name, char* path)
 *
 */

/*  This function goes ahead and triangulates the complement because:
 *   A) that's the only kernel function that applies to KLPProjections
 *   B)  since there is not kernel function to delete KLPProjections, we
 *   would leak a little memory.
 *
 */

extern "C" Triangulation* triangulate_link_complement_from_file(char* file_name, char *path){
  FILE* fp;
  CLinkProjectionWrapper myWrapper;
  CLinkProjectionPanorama  myPanorama;
  KLPProjection* projection;
  Triangulation* manifold;
  char full_name[500] = "";

  strcat(full_name, path);
  if (strlen(path) > 0 && path[strlen(path) - 1] != '/')
    strcat(full_name, "/");
  strcat(full_name, file_name);
  fp = fopen(full_name, "r");
  if (fp){
    myWrapper.ReadFile(fp);
    fclose(fp);
    myPanorama.ContentsToWindow(&myWrapper);
    projection =  myPanorama.CreateKLPProjection();
    manifold = triangulate_link_complement(projection);
    myPanorama.FreeKLPProjection(projection);
    set_triangulation_name(manifold, file_name);
    return manifold;
  }
  return NULL;
}

KLPProjection* load_link_complement_from_file(char* file_name, char *path){
  FILE* fp;
  CLinkProjectionWrapper myWrapper;
  CLinkProjectionPanorama  myPanorama;
  KLPProjection* projection;
  Triangulation* manifold;
  char full_name[500] = "";

  strcat(full_name, path);
  if (strlen(path) > 0 && path[strlen(path) - 1] != '/')
    strcat(full_name, "/");
  strcat(full_name, file_name);
  fp = fopen(full_name, "r");
  if (fp){
    myWrapper.ReadFile(fp);
    fclose(fp);
    myPanorama.ContentsToWindow(&myWrapper);
    projection =  myPanorama.CreateKLPProjection();
    return projection;
  }
  return NULL;
}

void* update_link_format(char* orig_file, char* new_file){
  FILE* fp;
  CLinkProjectionWrapper myWrapper;
  fp = fopen(orig_file, "r");
  if (fp){
    myWrapper.ReadFile(fp);
    fclose(fp);
    fp = fopen(new_file, "w");
    myWrapper.WriteFile(fp);
    fclose(fp);
  }
  return NULL;
}
