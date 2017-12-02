/* Walter Martin, Ryan Wallace
 * CS 124 Programming Assignment 1
 * Feb 19, 2016 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>

#define HIGH 10.0

struct mst_node;

typedef struct vertex {
    int id;
    float* position;
    float p_len;
    float e_len;
    struct vertex* v_prev;
    struct mst_node* n_prev;
    bool S;
} vertex;

typedef struct {
    int id;
    int num_v;
    int dim;
    vertex** verts;
    float** edges;
} graph;

typedef struct mst_node {
    vertex* v_ptr;
    struct mst_node* parent;
    int num_children;
    struct mst_node** children;
} mst_node;

typedef struct {
    //int num_v;
    vertex** q;
    int last;
} min_heap;

int max(int a, int b);
int min(int a, int b);

float** getHalfMatrix (int rows, int columns)
{
    // allocate an array of row pointers to pointers to floats
    float** matrixMemory = calloc (rows, sizeof(float*));
    if (matrixMemory == NULL)
    {
        return NULL;
    }
    
    // for each pointer in the array
    for (int i = 0; i < rows; i++)
    {
        // allocate a pointer to a float
        matrixMemory[i] = calloc (i + 1, sizeof(float*));
        if (matrixMemory[i] == NULL)
        {
            return NULL;
        }
    }
    
    // return the allocated for the matrix
    return (matrixMemory);
}

graph* init_rand(int n, int id, float max_weight) 
{
    graph* gr = calloc (1, sizeof(graph));
    (*gr).id = id;
    (*gr).num_v = n;

    // get space for adjacency matrix
    (*gr).edges = getHalfMatrix(n, n);
    if ((*gr).edges == NULL)
    {
        return NULL;
    }

    // get space for array of vertex pointers
    (*gr).verts = calloc (n, sizeof (vertex*));
    if ((*gr).verts == NULL)
    {
        return NULL;
    }
    
    // initialize lower half of matrix randomly, upper half to complement
    for (int i = 0; i < n; i++) 
    {
        // calloc space for each vertex
        (*gr).verts[i] = calloc(1, sizeof(vertex));
        if ((*gr).verts[i] == NULL)
        {
            return NULL;
        }

        // initialize the vertex
        (*gr).verts[i]->S = false; // not in MST yet
        (*gr).verts[i]->id = i;
        (*gr).verts[i]->v_prev = NULL;
        (*gr).verts[i]->n_prev = NULL;
        (*gr).verts[i]->position = NULL;

        for (int j = 0; j < i; j++) 
        {
            (*gr).edges[i][j] = (float)(rand()%100000)/100000;
            if ((*gr).edges[i][j] > max_weight)
            {
                (*gr).edges[i][j] = HIGH;
            }
        }
    }

    // initialize diagonal with 0s
    for (int i = 0; i < n; i++)
        gr->edges[i][i] = 0;

    return gr;
}

graph* init_rand_dim(int n, int d, int id, float max_weight)
{
    graph* gr = calloc (1, sizeof(graph));
    (*gr).id = id;
    (*gr).num_v = n;

    // get space for adjacency matrix
    (*gr).edges = getHalfMatrix(n, n);
    if ((*gr).edges == NULL)
    {
        return NULL;
    }

    // get space for array of vertex pointers
    (*gr).verts = calloc (n, sizeof (vertex*));
    if ((*gr).verts == NULL)
    {
        return NULL;
    }
    
    for (int i = 0; i < n; i++) 
    {
        (*gr).verts[i] = calloc (1, sizeof(vertex));
        if ((*gr).verts[i] == NULL)
        {
            return NULL;
        }

        (*gr).verts[i]->S = false; // not in MST yet
        (*gr).verts[i]->id = i;
        (*gr).verts[i]->v_prev = NULL;
        (*gr).verts[i]->n_prev = NULL; 
        (*gr).verts[i]->position = NULL;         

        // get space for array of float pointers
        (*gr).verts[i]->position = calloc (d, sizeof (float));
        if ((*gr).verts[i]->position == NULL)
        {
            return NULL;
        }

        // assign random float to each coordinate
        for (int dim = 0; dim < d; dim++)
        {   
            (*gr).verts[i]->position[dim] = (float)(rand()%100000000)/100000000;
        }
    }

    for(int i = 0; i < n; i++) 
    {
        for(int j = 0; j < i + 1; j++) 
        {
            float sum_sq_diffs = 0;
            for(int dim = 0; dim < d; dim++) 
            {
                sum_sq_diffs += pow ((((*gr).verts[i]->position[dim]) - ((*gr).verts[j]->position[dim])), 2);
            }
            (*gr).edges[i][j] = sqrt (sum_sq_diffs);
            if ((*gr).edges[i][j] > max_weight)
            {
                (*gr).edges[i][j] = HIGH;
            }
        }
    }
    return gr;
}

vertex deletemin(min_heap* heap, int nodes) {
    
    vertex root = *heap->q[0];
    
    //if there's only one element in the heap
    if(heap->q[1] == NULL) {
        heap->q[0] = NULL;
        heap->last--; 
        return root;
    }
    
    //moves last element to the root
    heap->q[0] = heap->q[heap->last];
    heap->q[heap->last] = NULL;
    heap->last--;

    int index = 0;
    
    while(index < nodes && heap->q[2*(index+1) - 1] != NULL) {
        //path lengths of the vertices on hand
        float rt1 = heap->q[index]->p_len;
        float c1 = heap->q[2*(index+1)-1]->p_len;
            //Where second child doesn't exist        
            if(heap->q[2*(index+1)] != NULL) {
                float c2 = heap->q[2*(index+1)]->p_len;
                if(c1 <= c2) {
                    //child 1 has shortest path
                    if(c1 < rt1) {
                        vertex* temp = heap->q[index];
                        heap->q[index] = heap->q[2*(index+1)-1];
                        heap->q[2*(index+1)-1] = temp;
                        index = 2*(index+1)-1;
                    }
                    //root has shortest path
                    else {
                        return root;
                    }
                }
                else if (c2 < rt1) {
                    vertex* temp = heap->q[index];
                    heap->q[index] = heap->q[2*(index+1)];
                    heap->q[2*(index+1)] = temp;
                    index = 2*(index+1);
                }
                else { 
                    return root;
                }
            }
            else if (heap->q[2*(index+1)] == NULL && c1 < rt1) {
                vertex* temp = heap->q[index];
                heap->q[index] = heap->q[2*(index+1)-1];
                heap->q[2*(index+1)-1] = temp;
                index = 2*(index+1)-1;
            }
            else {
                return root;
            }
    } 
    
    return root;
}

void insert(min_heap* heap, vertex* v) 
{
    heap->last++;
    int index = heap->last;
    heap->q[index] = v;

    //path lengths for current vertex and parent vertex
    float cur = v->p_len;
    float par = 0;
    if(index == 0) {
        par = 0;
    }
    else {
        par = heap->q[(index+1)/2 - 1]->p_len;
    }
    while(cur < par && index != 0) 
    {
        //swap parent and child, update index to current position of v
        heap->q[index] = heap->q[(index+1)/2 - 1];
        heap->q[(index+1)/2 - 1] = v;
        index = (index+1)/2 - 1;
        if(index != 0) {
            par = heap->q[(index+1)/2 - 1]->p_len;
        }
        else {
            par = 0;
        }
    }
}

mst_node* prim(int size, graph* g) 
{
    // MST set and member index
    mst_node *top = calloc(1, sizeof(mst_node));
    mst_node *cur = top;
    top->num_children = 0;
    
    // priority heap of vertices in G
    min_heap* heap = calloc (1, sizeof(min_heap));
    heap->last = -1; 

    int nodes;
    if (size < 32768) {
        nodes = size*size;
    }
    else if (size < 65536) {
        nodes = size*366;
    }
    else {
        nodes = size*500;
    }

    heap->q = calloc(nodes, sizeof(vertex*));
    if (heap->q == NULL)
    {
        return NULL;
    }
    insert(heap, (g->verts[0]));
    
    // distances, initialized to inf
    g->verts[0]->p_len = 0;
    for(int i = 1; i < size; i++) 
    {
        g->verts[i]->p_len = FLT_MAX;
        g->verts[i]->e_len = 0;
    }

    // while the heap is non-empty
    while (heap->q[0] != NULL)
    {
        vertex v = deletemin(heap, nodes);
        g->verts[v.id]->S = true;
        if(v.id == 0) {
            top->v_ptr = &v;
        }
        else {
            mst_node *n_temp = calloc(1, sizeof(mst_node));
            n_temp->v_ptr = &v;
            cur = n_temp;
            cur->v_ptr->v_prev = v.v_prev;
            cur->parent = v.n_prev;
            cur->num_children = 0;
        }

        // point the vertex pointer to the given vertex in g
        cur->v_ptr = g->verts[v.id];
        
		// update the nodes parent to have the node as a child
        cur->children = NULL;
        if(cur->parent != NULL) {
            if (cur->parent->children == NULL) {
                cur->parent->children = calloc(1, sizeof(mst_node*));
            }
            else {
                mst_node** buffer = cur->parent->children;
                cur->parent->children = calloc(cur->parent->num_children + 1, sizeof(mst_node*));
                memcpy(cur->parent->children, buffer, cur->parent->num_children * sizeof(mst_node*));
                free(buffer);

            }
            if (cur->parent->children == NULL)
            {
                return NULL;
            }
            cur->parent->children[cur->parent->num_children] = cur;
            cur->parent->num_children++;
        }
        
        if(cur->v_ptr->id != 0) {
            cur->v_ptr->p_len = v.v_prev->p_len 
                + g->edges[max(cur->v_ptr->id,v.v_prev->id)]
                    [min(cur->v_ptr->id,v.v_prev->id)];
        }
        else {
            cur->v_ptr->p_len = 0;
        }

        for (int i = 0; i < size; i++)
        {
            // if attached vertex is not already in S and path length is short enough, insert
            if (!g->verts[i]->S && g->verts[i]->p_len > g->edges[max(v.id, i)][min(v.id, i)] && 
                g->edges[max(v.id, i)][min(v.id, i)] < HIGH)
            {
                g->verts[i]->p_len = g->edges[max(v.id, i)][min(v.id, i)]; 
                g->verts[i]->e_len = g->edges[max(v.id, i)][min(v.id, i)];
                g->verts[i]->v_prev = g->verts[v.id];
                g->verts[i]->n_prev = cur;
                insert (heap, (g->verts[i]));
            }
        }
    }  
    
    // free heap things
    free(heap->q);
    free(heap);

    return top;
}
 
// calculate total weight of the mst generated in prim
float calc_weight(graph* gr)
{
    float weight = 0;
    for(int i = 0; i < gr->num_v; i++) {
        weight += gr->verts[i]->e_len;
    }
    return weight;
}

void free_mst(mst_node* mst)
{
    if(mst != NULL)
    {
        for(int i = 0; i < mst->num_children; i++) {
            free_mst(mst->children[i]);
        }
        free(mst->children);
        
        free(mst);
    }
}

int main(int argc, char *argv[]) {
    // input validation
    int numpoints = atoi(argv[2]);
    int numtrials = atoi(argv[3]);
    int dimension = atoi(argv[4]);
    float max_weight = FLT_MAX;

    if (argc != 5 || numpoints < 1 || numtrials < 1 || dimension < 0 || dimension > 4 || dimension== 1)
    {
        printf("Use: randmst 0 numpoints numtrials dimension\n");
        return 1;
    }

    if (numpoints > 60000)
    {
        if (dimension == 0)
            max_weight = 0.01;

        if (dimension == 2)
            max_weight = 1.0;

        if (dimension == 3)
            max_weight = 1.3;

        if (dimension == 4)
            max_weight = 1.6;
    }
    
    // seed randomness
    int seed = time(NULL);
    srand(seed);

    graph* gr = NULL; // pointer to the graph
    mst_node* mst = NULL; //the number of the trial/id of the graph
    float weights[numtrials]; // sum of wieghts in each trial
    float average = 0; // average weight 

    for (int trial = 0; trial < numtrials; trial++)
    {
        if (dimension == 0)
        {
            gr = init_rand (numpoints, trial, max_weight);
            if (gr == NULL)
            {
                printf("Heap size insufficient.\n");
                return 1;
            }

            mst = prim (numpoints, gr);
            if (mst == NULL)
            {
                printf("Heap size insufficient.\n");
                return 1;
            }
        }

        else 
        {
            gr = init_rand_dim (numpoints, dimension, trial, max_weight);
            if (gr == NULL)
            {
                printf("Heap size insufficient.\n");
                return 1;
            }

            mst = prim (numpoints, gr);
            if (mst == NULL)
            {
                printf("Heap size insufficient.\n");
                return 1;
            }
        }

        // calc weights
        weights[trial] = calc_weight (gr); 
        
        // free everything

        // mst
        free_mst(mst);

        // edges matrix
        for (int i = 0; i < numpoints; i++)
        {
            free (gr->edges[i]);
        }
        free (gr->edges);

        // vertex matrix
        for (int i = 0; i < numpoints; i++) 
        {
            free(gr->verts[i]->position); 
            free(gr->verts[i]);
        }
        free(gr->verts);
        
        // graph
        free (gr);

        sleep (1); // for random number generator to re-seed
    }

    // calculate average weight of MST
    for (int trial = 0; trial < numtrials; trial++)
    {
        average += weights[trial];
    }
    average = average / numtrials;

    // output
    printf("%.5f %d %d %d\n", average, numpoints, numtrials, dimension);

    return 0;
}

int max(int a, int b) {
    if (a >= b) return a; else return b;
}

int min(int a, int b) {
    if (a <= b) return a; else return b;
}
