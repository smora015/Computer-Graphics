#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include <vector>
#include "point2d.h"



// Determines max y value out of 3 points.
// @param: number   The starting point in our points array for the 3 vertices.
int ymax( int number );

// Determines minimum y value out of 3 points
// @param: number   The starting point in our points array for the 3 vertices
int ymin( int number );

// Renders a triangular polygon
// @param: number   The starting point in our points array for the 3 vertices
void renderTPolygon( int number );

// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(int x, int y, float r, float g, float b);

// Renders line based on incremental DDA
void renderLine( int x0, int y0, int x1, int y1);





#endif
