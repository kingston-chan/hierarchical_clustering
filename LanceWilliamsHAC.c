// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

// Struct that holds the index of the 
// two clusters that is going to be 
// merged
typedef struct Cluster {
	int i;		
	int j;
} Cluster;

typedef struct DendogramArray *DendA;
struct DendogramArray {
	Dendrogram node;
};

static void free_matrix(double **dist, int size);
static double **init_matrix(Graph g);
static Cluster find_cluster(double **dist, int size, Cluster c);
static Dendrogram new_node();
static void free_dendA(DendA dendA, int size);
static void update_dist(double **dist, int size, DendA dendA, Graph g, 
						int method);
static void closest_cluster(Dendrogram node, double **dist, int c_index, 
							Graph g, DendA dendA, int size, int method);		
static void in_cluster(Dendrogram node, int vertex, int *found);
static int find_vertex(DendA dendA, int size, int vertex);
static void insert_dist_cell(double **dist, int i, int j, double distance, 
							 int method);
static void find_dist(AdjList v_list, DendA dendA, double **dist, int c_index, 
					  int size, int method);

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {
	int num_v = GraphNumVertices(g);
	DendA dendA = malloc(num_v * sizeof(DendA));

	// Initialise the dendogram array to store pointers
	// to each vertex
	for (int i = 0; i < num_v; i++) {
		dendA[i].node = new_node();
		dendA[i].node->vertex = i;
	}

	// Create a dist grid to store distance between clusters
	double **dist = init_matrix(g);

	// Insert distances between clusters corresponding to 
	// method
	update_dist(dist, num_v, dendA, g, method);

	int i = num_v;
	while (i > 2) {
		// Find the two closest clusters
		Cluster c = find_cluster(dist, i, c);
		Dendrogram ci = dendA[c.i].node;
		Dendrogram cj = dendA[c.j].node;

		int k = 0;
		// Update the dendogram array
		for (int j = 0; j < i; j++) {
			if (j != c.i && j != c.j) {
				dendA[k].node = dendA[j].node;
				k++;
			}
		}

		// Merge two clusters
		dendA[k].node = new_node();
		dendA[k].node->right = cj;
		dendA[k].node->left = ci;
		dendA[k].node->vertex = -1;
		// Number of clusters is decreased by 1
		i--;

		// Update dist grid with merged clusters
		update_dist(dist, i, dendA, g, method);
	}

	// Complete the dendogram
	Dendrogram d = new_node();
	d->left = dendA[0].node;
	d->right = dendA[1].node;

	free_dendA(dendA, i);
	free_matrix(dist, num_v);

	return d;
}

/**
 * Frees all memory associated with the given Dendrogram structure.
 */
void freeDendrogram(Dendrogram d) {
	if (d == NULL) return;
	freeDendrogram(d->left);
	freeDendrogram(d->right);
	free(d);
}

			/////////////////////////////
			///// HELPER FUNCTIONS //////
			/////////////////////////////

// Helper function to free memory allocated to matrix
static void free_matrix(double **dist, int size) {
	for (int i = 0; i < size; i++) {
		free(dist[i]);
	}
	free(dist);
}

// Helper function to initialise the distance values for
// the matrix
static double **init_matrix(Graph g) {
	int size = GraphNumVertices(g);
	double **dist = malloc(size * sizeof(*dist));

	for (int i = 0; i < size; i++) {
		dist[i] = malloc(size * sizeof(double));
		for (int j = 0; j < size; j++) {
			dist[i][j] = DBL_MAX;
		}
	}

	return dist;
}

// Find the two clusters with the smallest distance between each
// other
static Cluster find_cluster(double **dist, int size, Cluster c) {
	double smallest = dist[0][1];
	c.i = 0;
	c.j = 1;

	for (int row = 1; row < size; row++) {
		for (int col = 0; col < row; col++) {
			if (dist[row][col] < smallest) {
				c.i = row;
				c.j = col;
				smallest = dist[row][col];
			}
		}
	}

	return c;
}

// Create a new dendogram node
static Dendrogram new_node() {
	Dendrogram dnode = malloc(sizeof(*dnode));
	dnode->right = NULL;
	dnode->left = NULL;
	return dnode;
}

// Free the dendogram array
static void free_dendA(DendA dendA, int size) {
	for (int i = 0; i < size; i++) {
		dendA[i].node = NULL;
	}
	free(dendA);
}

// Update the dist grid
static void update_dist(double **dist, int size, DendA dendA, Graph g, 
						int method) {
	// Set all dist grid values to infinity 
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			dist[i][j] = DBL_MAX;
		}
	}
	// Determine the values for dist grid corresponding to clustering
	// method
	for (int i = 0; i < size; i++) {
		closest_cluster(dendA[i].node, dist, i, g, dendA, size, method);
	}
}

// Determine the distance values for each cluster and input the distances 
// corresponding to the clustering method
static void closest_cluster(Dendrogram node, double **dist, int c_index, 
							Graph g, DendA dendA, int size, int method) {
	if (node->right == NULL && node->left == NULL) {
		// Found leaf node
		AdjList v_out = GraphOutIncident(g, node->vertex);
		AdjList v_in = GraphInIncident(g, node->vertex);

		find_dist(v_out, dendA, dist, c_index, size, method);
		find_dist(v_in, dendA, dist, c_index, size, method);
	} else {
		closest_cluster(node->left, dist, c_index, g, dendA, size, method);
		closest_cluster(node->right, dist, c_index, g, dendA, size, method);
	}
}

// Check if vertex is in current cluster
static void in_cluster(Dendrogram node, int vertex, int *found) {
	// Vertices will always be leaf nodes
	if (node->right == NULL && node->left == NULL) {
		if (node->vertex == vertex) {
			*found = 1;
		}
	} else {
		in_cluster(node->left, vertex, found);
		in_cluster(node->right, vertex, found);
	}
}

// Find the given vertex in the dendA array
static int find_vertex(DendA dendA, int size, int vertex) {
	int i = 0;
	for (; i < size; i++) {
		int found_v = 0;
		in_cluster(dendA[i].node, vertex, &found_v);
		if (found_v) break;
	}
	
	return i;
}

// Insert distance into dist grid corresponding to clustering method
// If method does not correspond to clustering method insert given 
// distance into grid
static void insert_dist_cell(double **dist, int i, int j, double distance, 
							 int method) {
	if (method == SINGLE_LINKAGE) {
		dist[i][j] = min(dist[i][j], distance);
		dist[j][i] = min(dist[j][i], distance);
	} else {
		if (dist[i][j] == DBL_MAX && dist[j][i] == DBL_MAX) {
			dist[i][j] = distance;
			dist[j][i] = distance;
		} else {
			dist[i][j] = max(dist[i][j], distance);
			dist[j][i] = max(dist[j][i], distance);
		}
	}
}

// Finds the distance between clusters and determines the distance
// to input into the dist grid by calling the insert_dist_cell function
static void find_dist(AdjList v_list, DendA dendA, double **dist, 
					  int c_index, int size, int method) {
	while (v_list != NULL) {
		int found_vtx = 0;
		// Checks if adjacent vertex is in given vertex's cluster
		in_cluster(dendA[c_index].node, v_list->v, &found_vtx);
		if (!found_vtx) {
			// If adjacent vertex is not in given vertex's cluster,
			// find the adjacent vertex's index in the dendA array
			int vtx_i = find_vertex(dendA, size, v_list->v);
			double distance = 1/(double)v_list->weight;
			// Insert to dist grid corresponding to clustering method
			insert_dist_cell(dist, c_index, vtx_i, distance, method);
		}
		v_list = v_list->next;
	}
}