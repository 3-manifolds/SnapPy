/* Prototypes for functions defined in gluing_equations.c */

extern int** get_gluing_equations(Triangulation* manifold, 
				  int* num_rows, 
				  int* num_cols);

extern void free_gluing_equations(int** equations, 
				  int num_rows);

extern int* get_cusp_equation(Triangulation* manifold,
			      int cusp_num,
			      int m,
			      int l,
			      int* num_rows);

extern void free_cusp_equation(int* equation);

/* Prototype for the function defined in load_link_proj.cc */

extern Triangulation*    triangulate_link_complement_from_file(char* file_name,
								 char* path);

/* Prototype for the function defined in braid.cc */

extern Triangulation* fibered_manifold_associated_to_braid(int numStrands,
							     int braidLength,
							     int* word);
