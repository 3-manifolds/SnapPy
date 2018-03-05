#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <string>
#include "global.h"
#include "parsing.h"
#include "twister.h"

void construct_manifold(manifold &M, std::string surface_file_contents, std::string gluing, std::string handles);

#endif
