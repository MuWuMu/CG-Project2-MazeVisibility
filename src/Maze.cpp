/************************************************************************
	 File:        Maze.cpp

	 Author:
				  Stephen Chenney, schenney@cs.wisc.edu
	 Modifier
				  Yu-Chi Lai, yu-chi@cs.wisc.edu

	 Comment:
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for Maze class. Manages the maze.


	 Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "Maze.h"
#include "Edge.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <FL/Fl.h>
#include <FL/fl_draw.h>
#include <GL/GL.h>

#include <iostream>

const char Maze::X = 0;
const char Maze::Y = 1;
const char Maze::Z = 2;

const float Maze::BUFFER = 0.1f;

//**********************************************************************
//
// * Constructor for the maze exception
//======================================================================
MazeException::
	MazeException(const char *m)
//======================================================================
{
	message = new char[strlen(m) + 4];
	strcpy(message, m);
}

//**********************************************************************
//
// * Constructor to create the default maze
//======================================================================
Maze::
	Maze(const int nx, const int ny, const float sx, const float sy)
//======================================================================
{
	// Build the connectivity structure.
	Build_Connectivity(nx, ny, sx, sy);

	// Make edges transparent to create a maze.
	Build_Maze();

	// Set the extents of the maze
	Set_Extents();

	// Default values for the viewer.
	viewer_posn[X] = viewer_posn[Y] = viewer_posn[Z] = 0.0;
	viewer_dir = 0.0;
	viewer_fov = 45.0;

	// Always start on the 0th frame.
	frame_num = 0;
}

//**********************************************************************
//
// * Construtor to read in precreated maze
//======================================================================
Maze::
	Maze(const char *filename)
//======================================================================
{
	char err_string[128];
	FILE *f;
	int i;

	// Open the file
	if (!(f = fopen(filename, "r")))
		throw new MazeException("Maze: Couldn't open file");

	// Get the total number of vertices
	if (fscanf(f, "%d", &num_vertices) != 1)
		throw new MazeException("Maze: Couldn't read number of vertices");

	// Read in each vertices
	vertices = new Vertex *[num_vertices];
	for (i = 0; i < num_vertices; i++)
	{
		float x, y;
		if (fscanf(f, "%g %g", &x, &y) != 2)
		{
			sprintf(err_string, "Maze: Couldn't read vertex number %d", i);
			throw new MazeException(err_string);
		}
		vertices[i] = new Vertex(i, x, y);
	}

	// Get the number of edges
	if (fscanf(f, "%d", &num_edges) != 1)
		throw new MazeException("Maze: Couldn't read number of edges");

	// read in all edges
	edges = new Edge *[num_edges];
	for (i = 0; i < num_edges; i++)
	{
		int vs, ve, cl, cr, o;
		float r, g, b;
		if (fscanf(f, "%d %d %d %d %d %g %g %g",
				   &vs, &ve, &cl, &cr, &o, &r, &g, &b) != 8)
		{
			sprintf(err_string, "Maze: Couldn't read edge number %d", i);
			throw new MazeException(err_string);
		}
		edges[i] = new Edge(i, vertices[vs], vertices[ve], r, g, b);
		edges[i]->Add_Cell((Cell *)cl, Edge::LEFT);
		edges[i]->Add_Cell((Cell *)cr, Edge::RIGHT);
		edges[i]->opaque = o ? true : false;
	}

	// Read in the number of cells
	if (fscanf(f, "%d", &num_cells) != 1)
		throw new MazeException("Maze: Couldn't read number of cells");

	// Read in all cells
	cells = new Cell *[num_cells];
	for (i = 0; i < num_cells; i++)
	{
		int epx, epy, emx, emy;
		if (fscanf(f, "%d %d %d %d", &epx, &epy, &emx, &emy) != 4)
		{
			sprintf(err_string, "Maze: Couldn't read cell number %d", i);
			throw new MazeException(err_string);
		}
		cells[i] = new Cell(i, epx >= 0 ? edges[epx] : NULL,
							epy >= 0 ? edges[epy] : NULL,
							emx >= 0 ? edges[emx] : NULL,
							emy >= 0 ? edges[emy] : NULL);
		if (cells[i]->edges[0])
		{
			if (cells[i]->edges[0]->neighbors[0] == (Cell *)i)
				cells[i]->edges[0]->neighbors[0] = cells[i];
			else if (cells[i]->edges[0]->neighbors[1] == (Cell *)i)
				cells[i]->edges[0]->neighbors[1] = cells[i];
			else
			{
				sprintf(err_string,
						"Maze: Cell %d not one of edge %d's neighbors",
						i, cells[i]->edges[0]->index);
				throw new MazeException(err_string);
			}
		}

		if (cells[i]->edges[1])
		{
			if (cells[i]->edges[1]->neighbors[0] == (Cell *)i)
				cells[i]->edges[1]->neighbors[0] = cells[i];
			else if (cells[i]->edges[1]->neighbors[1] == (Cell *)i)
				cells[i]->edges[1]->neighbors[1] = cells[i];
			else
			{
				sprintf(err_string,
						"Maze: Cell %d not one of edge %d's neighbors",
						i, cells[i]->edges[1]->index);
				throw new MazeException(err_string);
			}
		}
		if (cells[i]->edges[2])
		{
			if (cells[i]->edges[2]->neighbors[0] == (Cell *)i)
				cells[i]->edges[2]->neighbors[0] = cells[i];
			else if (cells[i]->edges[2]->neighbors[1] == (Cell *)i)
				cells[i]->edges[2]->neighbors[1] = cells[i];
			else
			{
				sprintf(err_string,
						"Maze: Cell %d not one of edge %d's neighbors",
						i, cells[i]->edges[2]->index);
				throw new MazeException(err_string);
			}
		}
		if (cells[i]->edges[3])
		{
			if (cells[i]->edges[3]->neighbors[0] == (Cell *)i)
				cells[i]->edges[3]->neighbors[0] = cells[i];
			else if (cells[i]->edges[3]->neighbors[1] == (Cell *)i)
				cells[i]->edges[3]->neighbors[1] = cells[i];
			else
			{
				sprintf(err_string,
						"Maze: Cell %d not one of edge %d's neighbors",
						i, cells[i]->edges[3]->index);
				throw new MazeException(err_string);
			}
		}
	}

	if (fscanf(f, "%g %g %g %g %g",
			   &(viewer_posn[X]), &(viewer_posn[Y]), &(viewer_posn[Z]),
			   &(viewer_dir), &(viewer_fov)) != 5)
		throw new MazeException("Maze: Error reading view information.");

	// Some edges have no neighbor on one side, so be sure to set their
	// pointers to NULL. (They were set at -1 by the save/load process.)
	for (i = 0; i < num_edges; i++)
	{
		if (edges[i]->neighbors[0] == (Cell *)-1)
			edges[i]->neighbors[0] = NULL;
		if (edges[i]->neighbors[1] == (Cell *)-1)
			edges[i]->neighbors[1] = NULL;
	}

	fclose(f);

	Set_Extents();

	// Figure out which cell the viewer is in, starting off by guessing the
	// 0th cell.
	Find_View_Cell(cells[0]);

	frame_num = 0;
}

//**********************************************************************
//
// * Destructor must free all the memory allocated.
//======================================================================
Maze::
	~Maze(void)
//======================================================================
{
	int i;

	for (i = 0; i < num_vertices; i++)
		delete vertices[i];
	delete[] vertices;

	for (i = 0; i < num_edges; i++)
		delete edges[i];
	delete[] edges;

	for (i = 0; i < num_cells; i++)
		delete cells[i];
	delete[] cells;
}

//**********************************************************************
//
// * Randomly generate the edge's opaque and transparency for an empty maze
//======================================================================
void Maze::
	Build_Connectivity(const int num_x, const int num_y,
					   const float sx, const float sy)
//======================================================================
{
	int i, j, k;
	int edge_i;

	// Ugly code to allocate all the memory for a new maze and to associate
	// edges with vertices and faces with edges.

	// Allocate and position the vertices.
	num_vertices = (num_x + 1) * (num_y + 1);
	vertices = new Vertex *[num_vertices];
	k = 0;
	for (i = 0; i < num_y + 1; i++)
	{
		for (j = 0; j < num_x + 1; j++)
		{
			vertices[k] = new Vertex(k, j * sx, i * sy);
			k++;
		}
	}

	// Allocate the edges, and associate them with their vertices.
	// Edges in the x direction get the first num_x * ( num_y + 1 ) indices,
	// edges in the y direction get the rest.
	num_edges = (num_x + 1) * num_y + (num_y + 1) * num_x;
	edges = new Edge *[num_edges];
	k = 0;
	for (i = 0; i < num_y + 1; i++)
	{
		int row = i * (num_x + 1);
		for (j = 0; j < num_x; j++)
		{
			int vs = row + j;
			int ve = row + j + 1;
			edges[k] = new Edge(k, vertices[vs], vertices[ve],
								rand() / (float)RAND_MAX * 0.5f + 0.25f,
								rand() / (float)RAND_MAX * 0.5f + 0.25f,
								rand() / (float)RAND_MAX * 0.5f + 0.25f);
			k++;
		}
	}

	edge_i = k;
	for (i = 0; i < num_y; i++)
	{
		int row = i * (num_x + 1);
		for (j = 0; j < num_x + 1; j++)
		{
			int vs = row + j;
			int ve = row + j + num_x + 1;
			edges[k] = new Edge(k, vertices[vs], vertices[ve],
								rand() / (float)RAND_MAX * 0.5f + 0.25f,
								rand() / (float)RAND_MAX * 0.5f + 0.25f,
								rand() / (float)RAND_MAX * 0.5f + 0.25f);
			k++;
		}
	}

	// Allocate the cells and associate them with their edges.
	num_cells = num_x * num_y;
	cells = new Cell *[num_cells];
	k = 0;
	for (i = 0; i < num_y; i++)
	{
		int row_x = i * (num_x + 1);
		int row_y = i * num_x;
		for (j = 0; j < num_x; j++)
		{
			int px = edge_i + row_x + 1 + j;
			int py = row_y + j + num_x;
			int mx = edge_i + row_x + j;
			int my = row_y + j;
			cells[k] = new Cell(k, edges[px], edges[py], edges[mx], edges[my]);
			edges[px]->Add_Cell(cells[k], Edge::LEFT);
			edges[py]->Add_Cell(cells[k], Edge::RIGHT);
			edges[mx]->Add_Cell(cells[k], Edge::RIGHT);
			edges[my]->Add_Cell(cells[k], Edge::LEFT);
			k++;
		}
	}
}

//**********************************************************************
//
// * Add edges from cell to the set that are available for removal to
//   grow the maze.
//======================================================================
static void
Add_To_Available(Cell *cell, int *available, int &num_available)
//======================================================================
{
	int i, j;

	// Add edges from cell to the set that are available for removal to
	// grow the maze.

	for (i = 0; i < 4; i++)
	{
		Cell *neighbor = cell->edges[i]->Neighbor(cell);

		if (neighbor && !neighbor->counter)
		{
			int candidate = cell->edges[i]->index;
			for (j = 0; j < num_available; j++)
				if (candidate == available[j])
				{
					printf("Breaking early\n");
					break;
				}
			if (j == num_available)
			{
				available[num_available] = candidate;
				num_available++;
			}
		}
	}

	cell->counter = 1;
}

//**********************************************************************
//
// * Grow a maze by removing candidate edges until all the cells are
//   connected. The edges are not actually removed, they are just made
//   transparent.
//======================================================================
void Maze::
	Build_Maze()
//======================================================================
{
	Cell *to_expand;
	int index;
	int *available = new int[num_edges];
	int num_available = 0;
	int num_visited;
	int i;

	srand(time(NULL));

	// Choose a random starting cell.
	index = (int)floor((rand() / (float)RAND_MAX) * num_cells);
	to_expand = cells[index];
	Add_To_Available(to_expand, available, num_available);
	num_visited = 1;

	// Join cells up by making edges opaque.
	while (num_visited < num_cells && num_available > 0)
	{
		int ei;

		index = (int)floor((rand() / (float)RAND_MAX) * num_available);
		to_expand = NULL;

		ei = available[index];

		if (edges[ei]->neighbors[0] &&
			!edges[ei]->neighbors[0]->counter)
			to_expand = edges[ei]->neighbors[0];
		else if (edges[ei]->neighbors[1] &&
				 !edges[ei]->neighbors[1]->counter)
			to_expand = edges[ei]->neighbors[1];

		if (to_expand)
		{
			edges[ei]->opaque = false;
			Add_To_Available(to_expand, available, num_available);
			num_visited++;
		}

		available[index] = available[num_available - 1];
		num_available--;
	}

	for (i = 0; i < num_cells; i++)
		cells[i]->counter = 0;
}

//**********************************************************************
//
// * Go through all the vertices looking for the minimum and maximum
//   extents of the maze.
//======================================================================
void Maze::
	Set_Extents(void)
//======================================================================
{
	int i;

	min_xp = vertices[0]->posn[Vertex::X];
	max_xp = vertices[0]->posn[Vertex::X];
	min_yp = vertices[0]->posn[Vertex::Y];
	max_yp = vertices[0]->posn[Vertex::Y];
	for (i = 1; i < num_vertices; i++)
	{
		if (vertices[i]->posn[Vertex::X] > max_xp)
			max_xp = vertices[i]->posn[Vertex::X];
		if (vertices[i]->posn[Vertex::X] < min_xp)
			min_xp = vertices[i]->posn[Vertex::X];
		if (vertices[i]->posn[Vertex::Y] > max_yp)
			max_yp = vertices[i]->posn[Vertex::Y];
		if (vertices[i]->posn[Vertex::Y] < min_yp)
			min_yp = vertices[i]->posn[Vertex::Y];
	}
}

//**********************************************************************
//
// * Figure out which cell the view is in, using seed_cell as an
//   initial guess. This procedure works by repeatedly checking
//   whether the viewpoint is in the current cell. If it is, we're
//   done. If not, Point_In_Cell returns in new_cell the next cell
//   to test. The new cell is the one on the other side of an edge
//   that the point is "outside" (meaning that it might be inside the
//   new cell).
//======================================================================
void Maze::
	Find_View_Cell(Cell *seed_cell)
//======================================================================
{
	Cell *new_cell;

	//
	while (!(seed_cell->Point_In_Cell(viewer_posn[X], viewer_posn[Y],
									  viewer_posn[Z], new_cell)))
	{
		if (new_cell == 0)
		{
			// The viewer is outside the top or bottom of the maze.
			throw new MazeException("Maze: View not in maze\n");
		}

		seed_cell = new_cell;
	}

	view_cell = seed_cell;
}

float *Maze::LookAt(void)
{
	// Look at
	float viewer_pos[3] = {viewer_posn[Maze::Y], 0.0f, viewer_posn[Maze::X]}; // viewer(camera position)
	float center[3] = {														  // target position
					   viewer_pos[Maze::X] + sin(Maze::To_Radians(viewer_dir)),
					   viewer_pos[Maze::Y],
					   viewer_pos[Maze::Z] + cos(Maze::To_Radians(viewer_dir))};
	float up[3] = {0.0f, 1.0f, 0.0f}; // up vector

	float foward[3] = {
		-(center[0] - viewer_pos[0]),
		-(center[1] - viewer_pos[1]),
		-(center[2] - viewer_pos[2])};
	float foward_length = sqrt(foward[0] * foward[0] + foward[1] * foward[1] + foward[2] * foward[2]);
	float foward_normalized[3] = {foward[0] / foward_length, foward[1] / foward_length, foward[2] / foward_length};

	float left[3] = {
		(up[1] * foward_normalized[2] - up[2] * foward_normalized[1]),
		(up[2] * foward_normalized[0] - up[0] * foward_normalized[2]),
		(up[0] * foward_normalized[1] - up[1] * foward_normalized[0])};

	float lookat[16] = {
		left[0], up[0], foward_normalized[0], 0.0f,
		left[1], up[1], foward_normalized[1], 0.0f,
		left[2], up[2], foward_normalized[2], 0.0f,
		-(viewer_pos[0] * left[0] + viewer_pos[1] * left[1] + viewer_pos[2] * left[2]), -(viewer_pos[0] * up[0] + viewer_pos[1] * up[1] + viewer_pos[2] * up[2]), -(viewer_pos[0] * foward_normalized[0] + viewer_pos[1] * foward_normalized[1] + viewer_pos[2] * foward_normalized[2]), 1.0f};

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	return lookat;
}

float *Maze::Perspective(const float wOverh)
{
	// Perspective
	float fovy = viewer_fov;
	float aspect = wOverh;
	float zNear = 0.01;
	float zFar = 200;

	float f = 1.0f / tanf(fovy * 0.5f * (M_PI / 180.0f)); // calc cotangent

	float perspective[16] = {
		f / aspect, 0.0f, 0.0f, 0.0f,
		0.0f, f, 0.0f, 0.0f,
		0.0f, 0.0f, (zFar + zNear) / (zNear - zFar), -1.0f,
		0.0f, 0.0f, (2 * zFar * zNear) / (zNear - zFar), 0.0f};

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	return perspective;
}

void Maze::NDC(float edge[4][4], const float start[2], const float end[2], const float color[3], const float modelView[16], const float projection[16])
{
	float tmp[4][4] = {0.0f};

	for (int i = 0; i <= 3; i++)
	{
		for (int j = 0; j <= 3; j++)
		{
			for (int k = 0; k <= 3; k++)
			{
				tmp[k][i] += edge[k][j] * modelView[j * 4 + i];
			}
		}
	}
	for (int i = 0; i <= 3; i++)
	{
		for (int k = 0; k <= 3; k++)
		{
			edge[k][i] = tmp[k][i];
		}
	}

	for (int i = 0; i <= 3; i++)
	{
		for (int j = 0; j <= 3; j++)
		{
			for (int k = 0; k <= 3; k++)
			{
				tmp[k][i] += edge[k][j] * projection[j * 4 + i];
			}
		}
	}
	for (int i = 0; i <= 3; i++)
	{
		for (int k = 0; k <= 3; k++)
		{
			edge[k][i] = tmp[k][i];
		}
	}

	for (int i = 0; i <= 2; i++)
	{
		for (int k = 0; k <= 3; k++)
		{
			edge[k][i] /= edge[k][3];
		}
	}
}

void Maze::getInterSectionPoint(std::array<float, 2> &readyToPush, std::array<float, 2> start, std::array<float, 2> end, float frustStart[2], float frustEnd[2])
{
	float x1 = frustStart[X], y1 = frustStart[Y];
	float x2 = frustEnd[X], y2 = frustEnd[Y];
	float x3 = start[X], y3 = start[Y];
	float x4 = end[X], y4 = end[Y];

	readyToPush[X] = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) /
					 ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
	readyToPush[Y] = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) /
					 ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
}

std::vector<std::array<float, 2>> Maze::Clipping(float edge[4][4])
{
	// Frustum points
	const float POSITIVE = 0.5f;
	const float NEGATIVE = -0.5f;
	const float FRUSTUM_EDGES[4][2] = {{POSITIVE, POSITIVE}, {POSITIVE, NEGATIVE}, {NEGATIVE, NEGATIVE}, {NEGATIVE, POSITIVE}};

	std::vector<std::array<float, 2>> edgesClippedLastTime;
	for (int i = 0; i < 4; i++)
	{
		std::array<float, 2> tmp = {edge[i][X], edge[i][Y]};
		edgesClippedLastTime.push_back(tmp);
	}

	for (int frustIt = 0; frustIt < 4; frustIt++)
	{
		std::vector<std::array<float, 2>> clippedEdges; // store the clipped edges which clipped by this one frustum edge
		for (int wallEdgeIt = 0; wallEdgeIt < edgesClippedLastTime.size(); wallEdgeIt++)
		{
			float frustStart[2] = {FRUSTUM_EDGES[frustIt][X], FRUSTUM_EDGES[frustIt][Y]};
			float frustEnd[2] = {FRUSTUM_EDGES[(frustIt + 1) % 4][X], FRUSTUM_EDGES[(frustIt + 1) % 4][Y]};

			// discriminant, d = (x2 - x1)(y - y1) - (y2 - y1)(x - x1)
			// if d < 0, the point is on the right side of the line (in the frustum)
			// if d > 0, the point is on the left side of the line (out of the frustum)
			// if d = 0, the point is on the line
			float d1 =
				(frustEnd[X] - frustStart[X]) * (edgesClippedLastTime[wallEdgeIt][Y] - frustStart[Y]) - (frustEnd[Y] - frustStart[Y]) * ((edgesClippedLastTime[wallEdgeIt][X] - frustStart[X]));
			float d2 =
				(frustEnd[X] - frustStart[X]) * (edgesClippedLastTime[(wallEdgeIt + 1) % edgesClippedLastTime.size()][Y] - frustStart[Y]) - (frustEnd[Y] - frustStart[Y]) * ((edgesClippedLastTime[(wallEdgeIt + 1) % edgesClippedLastTime.size()][X] - frustStart[X]));

			// 4 cases
			if (d1 <= 0 && d2 <= 0) // in, in
			{						// return start point
				clippedEdges.push_back(edgesClippedLastTime[wallEdgeIt]);
			}
			else if (d1 > 0 && d2 > 0) // out, out
			{						   // return nothing
			}
			else if (d1 <= 0 && d2 > 0) // in, out
			{							// return start point and intersection point
				clippedEdges.push_back(edgesClippedLastTime[wallEdgeIt]);
				std::array<float, 2> readyToPush;
				getInterSectionPoint(readyToPush, edgesClippedLastTime[wallEdgeIt], edgesClippedLastTime[(wallEdgeIt + 1) % edgesClippedLastTime.size()], frustStart, frustEnd);
				clippedEdges.push_back(readyToPush);
			}
			else if (d1 > 0 && d2 <= 0) // out, in
			{							// return intersection
				std::array<float, 2> readyToPush;
				getInterSectionPoint(readyToPush, edgesClippedLastTime[wallEdgeIt], edgesClippedLastTime[(wallEdgeIt + 1) % edgesClippedLastTime.size()], frustStart, frustEnd);
				clippedEdges.push_back(readyToPush);
			}
		}
		edgesClippedLastTime = clippedEdges;
		clippedEdges.clear();
	}

	return edgesClippedLastTime;
}

void Maze::ClipIn2D(float wall_start[2], float wall_end[2], float frustum_edge[2][2], float color[3])
{
}

void Maze::Draw_Cell(Cell *targetCell, const float LPoint[2], const float RPoint[2])
{
}

//**********************************************************************
//
// * Move the viewer's position. This method will do collision detection
//   between the viewer's location and the walls of the maze and prevent
//   the viewer from passing through walls.
//======================================================================
void Maze::
	Move_View_Posn(const float dx, const float dy, const float dz)
//======================================================================
{
	Cell *new_cell;
	float xs, ys, zs, xe, ye, ze;

	// Move the viewer by the given amount. This does collision testing to
	// prevent walking through walls. It also keeps track of which cells the
	// viewer is in.

	// Set up a line segment from the start to end points of the motion.
	xs = viewer_posn[X];
	ys = viewer_posn[Y];
	zs = viewer_posn[Z];
	xe = xs + dx;
	ye = ys + dy;
	ze = zs + dz;

	// Fix the z to keep it in the maze.
	if (ze > 1.0f - BUFFER)
		ze = 1.0f - BUFFER;
	if (ze < BUFFER - 1.0f)
		ze = BUFFER - 1.0f;

	// Clip_To_Cell clips the motion segment to the view_cell if the
	// segment intersects an opaque edge. If the segment intersects
	// a transparent edge (through which it can pass), then it clips
	// the motion segment so that it _starts_ at the transparent edge,
	// and it returns the cell the viewer is entering. We keep going
	// until Clip_To_Cell returns NULL, meaning we've done as much of
	// the motion as is possible without passing through walls.
	while ((new_cell = view_cell->Clip_To_Cell(xs, ys, xe, ye, BUFFER)))
		view_cell = new_cell;

	// The viewer is at the end of the motion segment, which may have
	// been clipped.
	viewer_posn[X] = xe;
	viewer_posn[Y] = ye;
	viewer_posn[Z] = ze;
}

//**********************************************************************
//
// * Set the viewer's location
//======================================================================
void Maze::
	Set_View_Posn(float x, float y, float z)
//======================================================================
{
	// First make sure it's in some cell.
	// This assumes that the maze is rectangular.
	if (x < min_xp + BUFFER)
		x = min_xp + BUFFER;
	if (x > max_xp - BUFFER)
		x = max_xp - BUFFER;
	if (y < min_yp + BUFFER)
		y = min_yp + BUFFER;
	if (y > max_yp - BUFFER)
		y = max_yp - BUFFER;
	if (z < -1.0f + BUFFER)
		z = -1.0f + BUFFER;
	if (z > 1.0f - BUFFER)
		z = 1.0f - BUFFER;

	viewer_posn[X] = x;
	viewer_posn[Y] = y;
	viewer_posn[Z] = z;

	// Figure out which cell we're in.
	Find_View_Cell(cells[0]);
}

//**********************************************************************
//
// * Set the angle in which the viewer is looking.
//======================================================================
void Maze::
	Set_View_Dir(const float d)
//======================================================================
{
	viewer_dir = d;
}

//**********************************************************************
//
// * Set the horizontal field of view.
//======================================================================
void Maze::
	Set_View_FOV(const float f)
//======================================================================
{
	viewer_fov = f;
}

//**********************************************************************
//
// * Draws the map view of the maze. It is passed the minimum and maximum
//   corners of the window in which to draw.
//======================================================================
void Maze::
	Draw_Map(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int height;
	float scale_x, scale_y, scale;
	int i;

	// Figure out scaling factors and the effective height of the window.
	scale_x = (max_x - min_x - 10) / (max_xp - min_xp);
	scale_y = (max_y - min_y - 10) / (max_yp - min_yp);
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * (max_yp - min_yp));

	min_x += 5;
	min_y += 5;

	// Draw all the opaque edges.
	for (i = 0; i < num_edges; i++)
		if (edges[i]->opaque)
		{
			float x1, y1, x2, y2;

			x1 = edges[i]->endpoints[Edge::START]->posn[Vertex::X];
			y1 = edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
			x2 = edges[i]->endpoints[Edge::END]->posn[Vertex::X];
			y2 = edges[i]->endpoints[Edge::END]->posn[Vertex::Y];

			fl_color((unsigned char)floor(edges[i]->color[0] * 255.0),
					 (unsigned char)floor(edges[i]->color[1] * 255.0),
					 (unsigned char)floor(edges[i]->color[2] * 255.0));
			fl_line_style(FL_SOLID);
			fl_line(min_x + (int)floor((x1 - min_xp) * scale),
					min_y + height - (int)floor((y1 - min_yp) * scale),
					min_x + (int)floor((x2 - min_xp) * scale),
					min_y + height - (int)floor((y2 - min_yp) * scale));
		}
}

//**********************************************************************
//
// * Draws the first-person view of the maze. It is passed the focal distance.
//   THIS IS THE FUINCTION YOU SHOULD MODIFY.
//======================================================================
void Maze::
	Draw_View(const float aspect)
//======================================================================
{
	frame_num++;
	// std::cout << "Frame number: " << frame_num << std::endl;

	glClear(GL_DEPTH_BUFFER_BIT);

	float half_fov_rad = To_Radians(viewer_fov) * 0.5f;
	float right_dir = To_Radians(viewer_dir) - half_fov_rad;
	float left_dir = To_Radians(viewer_dir) + half_fov_rad;

	float right_dir_vector[2] = {cos(right_dir), sin(right_dir)};
	float left_dir_vector[2] = {cos(left_dir), sin(left_dir)};

	float right_line_end[2] = {viewer_posn[X], viewer_posn[Y]};
	float right_line_start[2] = {viewer_posn[X] + right_dir_vector[X], viewer_posn[Y] + right_dir_vector[Y]};
	float left_line_start[2] = {viewer_posn[X], viewer_posn[Y]};
	float left_line_end[2] = {viewer_posn[X] + left_dir_vector[X], viewer_posn[Y] + left_dir_vector[Y]};

	for (int i = 0; i < (int)this->num_edges; i++)
	{
		if (!this->edges[i]->opaque)
			continue;

		float edge_start[2] = {
			this->edges[i]->endpoints[Edge::START]->posn[Vertex::X],
			this->edges[i]->endpoints[Edge::START]->posn[Vertex::Y]};
		float edge_end[2] = {
			this->edges[i]->endpoints[Edge::END]->posn[Vertex::X],
			this->edges[i]->endpoints[Edge::END]->posn[Vertex::Y]};

		
		// Clipping the edge with the right line
		//======================================================================
		float d1 = (right_line_end[X] - right_line_start[X]) * (edge_start[Y] - right_line_start[Y]) -
				   (right_line_end[Y] - right_line_start[Y]) * (edge_start[X] - right_line_start[X]);
		float d2 = (right_line_end[X] - right_line_start[X]) * (edge_end[Y] - right_line_start[Y]) -
				   (right_line_end[Y] - right_line_start[Y]) * (edge_end[X] - right_line_start[X]);

		if (d1 <= 0 && d2 <= 0)
		{
			// Both points are on the right side, keep the entire edge
			// Do nothing
		}
		else if (d1 > 0 && d2 > 0)
		{
			// Both points are on the left side, skip the entire edge
			// Not drawing the edge
			continue;
		}
		else
		{
			// One point is on the right side, one point is on the left side
			// Calculate the intersection point
			float intersection[2] = {
				(((right_line_start[X] * right_line_end[Y] - right_line_start[Y] * right_line_end[X]) * (edge_start[X] - edge_end[X]) - 
				(right_line_start[X] - right_line_end[X]) * (edge_start[X] * edge_end[Y] - edge_start[Y] * edge_end[X])) /
				((right_line_start[X] - right_line_end[X]) * (edge_start[Y] - edge_end[Y]) -
				(right_line_start[Y] - right_line_end[Y]) * (edge_start[X] - edge_end[X]))),

				(((right_line_start[X] * right_line_end[Y] - right_line_start[Y] * right_line_end[X]) * (edge_start[Y] - edge_end[Y]) -
				(right_line_start[Y] - right_line_end[Y]) * (edge_start[X] * edge_end[Y] - edge_start[Y] * edge_end[X])) /
				((right_line_start[X] - right_line_end[X]) * (edge_start[Y] - edge_end[Y]) -
				(right_line_start[Y] - right_line_end[Y]) * (edge_start[X] - edge_end[X])))
				};

			if (d1 < 0)
			{
				// edge_start is on the right side, no need to update, in order, update edge_end
				edge_end[0] = intersection[0];
				edge_end[1] = intersection[1];
			}
			else
			{
				// edge_end is on the right side, no need to update, in order, update edge_start
				edge_start[0] = intersection[0];
				edge_start[1] = intersection[1];
			}
		}


		// Clipping the edge with the left line
		//======================================================================
		d1 = (left_line_end[X] - left_line_start[X]) * (edge_start[Y] - left_line_start[Y]) -
			 (left_line_end[Y] - left_line_start[Y]) * (edge_start[X] - left_line_start[X]);
		d2 = (left_line_end[X] - left_line_start[X]) * (edge_end[Y] - left_line_start[Y]) -
			 (left_line_end[Y] - left_line_start[Y]) * (edge_end[X] - left_line_start[X]);

		if (d1 <= 0 && d2 <= 0)
		{
			// Both points are on the right side, keep the entire edge
			// Do nothing
		}
		else if (d1 > 0 && d2 > 0)
		{
			// Both points are on the left side, skip the entire edge
			// Not drawing the edge
			continue;
		}
		else
		{
			// One point is on the right side, one point is on the left side
			// Calculate the intersection point
			float intersection[2] = {
				(((left_line_start[X] * left_line_end[Y] - left_line_start[Y] * left_line_end[X]) * (edge_start[X] - edge_end[X]) - 
				(left_line_start[X] - left_line_end[X]) * (edge_start[X] * edge_end[Y] - edge_start[Y] * edge_end[X])) /
				((left_line_start[X] - left_line_end[X]) * (edge_start[Y] - edge_end[Y]) -
				(left_line_start[Y] - left_line_end[Y]) * (edge_start[X] - edge_end[X]))),

				(((left_line_start[X] * left_line_end[Y] - left_line_start[Y] * left_line_end[X]) * (edge_start[Y] - edge_end[Y]) -
				(left_line_start[Y] - left_line_end[Y]) * (edge_start[X] * edge_end[Y] - edge_start[Y] * edge_end[X])) /
				((left_line_start[X] - left_line_end[X]) * (edge_start[Y] - edge_end[Y]) -
				(left_line_start[Y] - left_line_end[Y]) * (edge_start[X] - edge_end[X])))
				};

			if (d1 < 0)
			{
				// edge_start is on the right side, no need to update, in order, update edge_end
				edge_end[0] = intersection[0];
				edge_end[1] = intersection[1];
			}
			else
			{
				// edge_end is on the right side, no need to update, in order, update edge_start
				edge_start[0] = intersection[0];
				edge_start[1] = intersection[1];
			}
		}


		Draw_Wall(edge_start, edge_end, this->edges[i]->color, this->LookAt(), this->Perspective(aspect));
	}

	// // glEnable(GL_DEPTH_TEST);
	// for (int i = 0; i < (int)this->num_edges; i++)
	// {
	// 	float edge_start[2] = {
	// 		this->edges[i]->endpoints[Edge::START]->posn[Vertex::X],
	// 		this->edges[i]->endpoints[Edge::START]->posn[Vertex::Y]};
	// 	float edge_end[2] = {
	// 		this->edges[i]->endpoints[Edge::END]->posn[Vertex::X],
	// 		this->edges[i]->endpoints[Edge::END]->posn[Vertex::Y]};

	// 	float color[3] = {this->edges[i]->color[0], this->edges[i]->color[1], this->edges[i]->color[2]};

	// 	if (this->edges[i]->opaque)
	// 	{
	// 		Draw_Wall(edge_start, edge_end, color, this->LookAt(), this->Perspective(aspect));
	// 	}
	// }

	// remove the for loop above
	// do the clipping
	// and draw the wall with the clipped edges
}

void Maze::
	Draw_Wall(const float start[2], const float end[2], const float color[3], const float modelView[16], const float projection[16])
{
	float edge[4][4] = {
		{start[Y], 1.0f, start[X], 1.0f},
		{end[Y], 1.0f, end[X], 1.0f},
		{end[Y], -1.0f, end[X], 1.0f},
		{start[Y], -1.0f, start[X], 1.0f}};

	NDC(edge, start, end, color, modelView, projection);

	glBegin(GL_POLYGON);
	glColor3fv(color);
	for (int i = 0; i < 4; i++)
	{
		glVertex2f(edge[i][X], edge[i][Y]);
	}
	glEnd();


	// std::vector<std::array<float, 2>> clippedEdge = Clipping(edge);

	// glBegin(GL_POLYGON);
	// glColor3fv(color);
	// for (int i = 0; i < clippedEdge.size(); i++)
	// {
	// 	glVertex2f(clippedEdge[i][X], clippedEdge[i][Y]);
	// }
	// glEnd();
}

//**********************************************************************
//
// * Draws the frustum on the map view of the maze. It is passed the
//   minimum and maximum corners of the window in which to draw.
//======================================================================
void Maze::
	Draw_Frustum(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int height;
	float scale_x, scale_y, scale;
	float view_x, view_y;

	// Draws the view frustum in the map. Sets up all the same viewing
	// parameters as draw().
	scale_x = (max_x - min_x - 10) / (max_xp - min_xp);
	scale_y = (max_y - min_y - 10) / (max_yp - min_yp);
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * (max_yp - min_yp));

	min_x += 5;
	min_y += 5;

	view_x = (viewer_posn[X] - min_xp) * scale;
	view_y = (viewer_posn[Y] - min_yp) * scale;
	fl_line(min_x + (int)floor(view_x +
							   cos(To_Radians(viewer_dir + viewer_fov / 2.0)) * scale),
			min_y + height -
				(int)floor(view_y +
						   sin(To_Radians(viewer_dir + viewer_fov / 2.0)) *
							   scale),
			min_x + (int)floor(view_x),
			min_y + height - (int)floor(view_y));
	fl_line(min_x + (int)floor(view_x +
							   cos(To_Radians(viewer_dir - viewer_fov / 2.0)) *
								   scale),
			min_y + height -
				(int)floor(view_y + sin(To_Radians(viewer_dir - viewer_fov / 2.0)) *
										scale),
			min_x + (int)floor(view_x),
			min_y + height - (int)floor(view_y));
}

//**********************************************************************
//
// * Draws the viewer's cell and its neighbors in the map view of the maze.
//   It is passed the minimum and maximum corners of the window in which
//   to draw.
//======================================================================
void Maze::
	Draw_Neighbors(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int height;
	float scale_x, scale_y, scale;
	int i, j;

	// Draws the view cell and its neighbors in the map. This works
	// by drawing just the neighbor's edges if there is a neighbor,
	// otherwise drawing the edge. Every edge is shared, so drawing the
	// neighbors' edges also draws the view cell's edges.

	scale_x = (max_x - min_x - 10) / (max_xp - min_xp);
	scale_y = (max_y - min_y - 10) / (max_yp - min_yp);
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * (max_yp - min_yp));

	min_x += 5;
	min_y += 5;

	for (i = 0; i < 4; i++)
	{
		Cell *neighbor = view_cell->edges[i]->Neighbor(view_cell);

		if (neighbor)
		{
			for (j = 0; j < 4; j++)
			{
				Edge *e = neighbor->edges[j];

				if (e->opaque)
				{
					float x1, y1, x2, y2;

					x1 = e->endpoints[Edge::START]->posn[Vertex::X];
					y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
					x2 = e->endpoints[Edge::END]->posn[Vertex::X];
					y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

					fl_color((unsigned char)floor(e->color[0] * 255.0),
							 (unsigned char)floor(e->color[1] * 255.0),
							 (unsigned char)floor(e->color[2] * 255.0));
					fl_line_style(FL_SOLID);
					fl_line(min_x + (int)floor((x1 - min_xp) * scale),
							min_y + height - (int)floor((y1 - min_yp) * scale),
							min_x + (int)floor((x2 - min_xp) * scale),
							min_y + height - (int)floor((y2 - min_yp) * scale));
				}
			}
		}
		else
		{
			Edge *e = view_cell->edges[i];

			if (e->opaque)
			{
				float x1, y1, x2, y2;

				x1 = e->endpoints[Edge::START]->posn[Vertex::X];
				y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
				x2 = e->endpoints[Edge::END]->posn[Vertex::X];
				y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

				fl_color((unsigned char)floor(e->color[0] * 255.0),
						 (unsigned char)floor(e->color[1] * 255.0),
						 (unsigned char)floor(e->color[2] * 255.0));
				fl_line_style(FL_SOLID);
				fl_line(min_x + (int)floor((x1 - min_xp) * scale),
						min_y + height - (int)floor((y1 - min_yp) * scale),
						min_x + (int)floor((x2 - min_xp) * scale),
						min_y + height - (int)floor((y2 - min_yp) * scale));
			}
		}
	}
}

//**********************************************************************
//
// * Save the maze to a file of the given name.
//======================================================================
bool Maze::
	Save(const char *filename)
//======================================================================
{
	FILE *f = fopen(filename, "w");
	int i;

	// Dump everything to a file of the given name. Returns false if it
	// couldn't open the file. True otherwise.

	if (!f)
	{
		return false;
	}

	fprintf(f, "%d\n", num_vertices);
	for (i = 0; i < num_vertices; i++)
		fprintf(f, "%g %g\n", vertices[i]->posn[Vertex::X],
				vertices[i]->posn[Vertex::Y]);

	fprintf(f, "%d\n", num_edges);
	for (i = 0; i < num_edges; i++)
		fprintf(f, "%d %d %d %d %d %g %g %g\n",
				edges[i]->endpoints[Edge::START]->index,
				edges[i]->endpoints[Edge::END]->index,
				edges[i]->neighbors[Edge::LEFT] ? edges[i]->neighbors[Edge::LEFT]->index : -1,
				edges[i]->neighbors[Edge::RIGHT] ? edges[i]->neighbors[Edge::RIGHT]->index : -1,
				edges[i]->opaque ? 1 : 0,
				edges[i]->color[0], edges[i]->color[1], edges[i]->color[2]);

	fprintf(f, "%d\n", num_cells);
	for (i = 0; i < num_cells; i++)
		fprintf(f, "%d %d %d %d\n",
				cells[i]->edges[0] ? cells[i]->edges[0]->index : -1,
				cells[i]->edges[1] ? cells[i]->edges[1]->index : -1,
				cells[i]->edges[2] ? cells[i]->edges[2]->index : -1,
				cells[i]->edges[3] ? cells[i]->edges[3]->index : -1);

	fprintf(f, "%g %g %g %g %g\n",
			viewer_posn[X], viewer_posn[Y], viewer_posn[Z],
			viewer_dir, viewer_fov);

	fclose(f);

	return true;
}
