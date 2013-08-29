/********************************************************
 C code to handle directed graphs in the context of
 modelling processes modulating trait evolution along
 phylogeny.
 The code is a work in progress and may therefore contain
 more functionalities than those involved in the functions
 functions called from R.
 Guillaume Guenard - Universite de Montreal - 2010-2012
 C Header
 Notes:
 1. location of the R headers on Ubuntu 12.04:
 /usr/share/R/include
 2. LAPACK functions are in R_ext/Lapack.R
    dgesvd() for svd.
*********************************************************/

// Defines
#define length 0x0
#define meanchar 0x0
#define charvalue 0x1
#define missing 0x0       // ! check for a more suitable missing value marker.
// #define with_testing

#ifndef M_LN_2PI
#define M_LN_2PI 1.837877066409345     /* log(2*pi) */
#endif

// Type declarations

// Two mutually contained structures to represent directed graphs.
struct dedge {
  int id;                 // Edge identification number (should equal array indices).
  unsigned int nv;        // Number of values at a particular edge.
  double* v;              // Values at edge (most commonly its length).
  struct dvertex* u;      // upward-connected directed vertex.
  struct dvertex* d;      // downward-connected directed vertex.
};

struct dvertex {
  int id;                 // vertex identification number (should equal array indices).
  unsigned int nv;        // Number of values at a particular vertex.
  double* v;              // Values at vertex.
  unsigned int nu;        // Number of upward-connected directed edges.
  struct dedge** u;       // upward-connected directed edges.
  unsigned int nd;        // Number of downward-connected directed edges.
  struct dedge** d;       // downward-connected directed edges.
};

// Structure holding directed graphs.
struct dgraph {
  char *id;               // Graph identification string.
  unsigned int ne;        // Number of edge involved in the graph.
  struct dedge* de;       // The list of edges.
  char **elabels;         // Edge labels.
  unsigned int nn;        // Number of vertices involved in the graph.
  struct dvertex* dn;     // The list of vertices.
  char **nlabels;         // vertex labels.
};

// Structure to wrap arrays into matrices.
struct matrix {
  char* id;               // Matrix idenfification string.
  unsigned int nr;        // Number of row(s) of the matrix.
  unsigned int nc;        // Number of column(s) of the matrix.
  double* v;              // Values attached to the matrix, ordered by column(s).
};

// C functions declarations.
// Edge functions:
struct dedge* allocdedge(unsigned int ne);
struct dedge* reallocdedge(struct dedge* de, unsigned int ne);
struct dedge* initdedge(struct dedge* de, unsigned int start, unsigned int ne);
void assigndedgevalues(struct dedge* de, unsigned int ne, double* ev, unsigned int nev);
struct dedge* freededge(struct dedge* de);

// Vertex functions:
struct dvertex* allocdvertex(unsigned int nn);
struct dvertex* reallocdvertex(struct dvertex *dn, unsigned int nn);
struct dvertex* initdvertex(struct dvertex* dn, unsigned int start, unsigned int nn);
void assigndvertexvalues(struct dvertex* dn, unsigned int nn, double* nv, unsigned int nnv);
struct dvertex* evalallocdvertexres(struct dvertex* dn, unsigned int nn, int* a, int* b, unsigned int nr);
struct dvertex freedvertexres(struct dvertex dn);  // Called by freedvertex.
struct dvertex* freedvertex(struct dvertex* dn, unsigned int nn);

// Graph functions.
struct dgraph initdgraph(char* id, unsigned int ne, char** elabels, unsigned int nn, char** nlabels);
void assigndgraphvalues(struct dgraph* dgr, double* ev, unsigned int nev, double* nv, unsigned int nnv);
void makedgraph(int* a, int* b, struct dgraph* dgr);
void freedgraph(struct dgraph* dgr);

// Matrix functions.
struct matrix initmatrix(char* id, unsigned int nr, unsigned int nc);
struct matrix assignmatrix(char* id, unsigned int nr, unsigned int nc, double* v);
void freematrix(struct matrix* mat);
void deassignmatrix(struct matrix* mat);
struct matrix copymatrix(struct matrix *a);
void rowsums(struct matrix *a, double *s);
void colsums(struct matrix *a, double *s);
void rowcentering(struct matrix *a, struct matrix *b, double *c);
void colcentering(struct matrix *a, struct matrix *b, double *c);
void rowcentermeans(struct matrix *a, struct matrix *b, double *m);
void colcentermeans(struct matrix *a, struct matrix *b, double *m);
void rowweighting(struct matrix *a, struct matrix *b, double *w);
void colweighting(struct matrix *a, struct matrix *b, double *w);
void addmatrix(struct matrix *a, struct matrix *b, struct matrix *c);
void subtractmatrix(struct matrix *a, struct matrix *b, struct matrix *c);
void matrixscalar(struct matrix *a, double b, struct matrix *c);
void matrixdotproduct(struct matrix *a, struct matrix *b, struct matrix *c);
void matrixproduct(struct matrix *a, struct matrix *b, struct matrix *c);
void matrixweightedproduct(struct matrix *a, double*d, struct matrix *b, struct matrix *c);
void matrixtransproduct(struct matrix *a, struct matrix *b, struct matrix *c);
void matrixproducttrans(struct matrix *a, struct matrix *b, struct matrix *c);
void matrixproductweightedtrans(struct matrix *a, double *d, struct matrix *b, struct matrix *c);
void getdiagonal(struct matrix *mat, double *a);
void getrow(struct matrix *mat, unsigned int i, double *a);
void getcolumn(struct matrix *mat, unsigned int j, double *a);

// Auxiliary functions.
void InfluenceRD(struct dgraph* dgr, unsigned int e, int* out);
unsigned int rselect(double* prob, unsigned int n);
void evolveqcalongtree(struct dgraph* dgr, double* tw, unsigned int ntw, unsigned int sr, unsigned int nnv);
void OUdedgecoefs(double* ev, double* lg, unsigned int ne, double alpha, double sigma);
void simOUprocess(struct dgraph* dgr, unsigned int sr, unsigned int n, double* out);
void PEMvar(double* d, int* nd, double* a, double* psi, double* res);
void PEMweight(double* d, int* nd, double* a, double* psi, double* res);
void Psquared(double* p, double* o, int* n, double* res);

// Testing functions.
#ifdef with_testing
void checkdedge(struct dedge* de, unsigned int ne);
void checkdedgevalues(struct dedge* de, unsigned int ne);
void checkdvertex(struct dvertex* dn, unsigned int nn);
void checkdvertexvalues(struct dvertex* dn, unsigned int nn);
void checkdgraph(struct dgraph* dgr);
void checkdgraphvalues(struct dgraph* dgr);
void checkmatrix(struct matrix* mat);
// void test_function(double *mat1, double *mat2, double *res, int *rmat1, int *cmat2, int *p);
// void test_function2(double *mat1, double *mat2, double *res, int *rmat1, int *rmat2, int *cols);
// void test_function3(double *mat1, double *d, double *mat2, double *res, int *rmat1, int *rmat2, int *cols);
// void test_function4(double *mat1, int* dimmat1, double *m);
// void test_function5(double *mat1, int* dimmat1, double *m);
// void test_function6(double *mat1, double *d, double *mat2, double *res, int *rmat1, int *cmat2, int *crmats);
#endif

// R functions:
void PEMInfMat(int* from, int* to, int* ne, int* nn, int* out);
void EvolveQC(int* from, int* to, int* ne, int* nn, double* nv, double* tw, int* ntw, int* anc, int* n, int* sr);
void OUsim(int* from, int* to, int* ne, int* nn, double* lg, double* alpha, double* sigma, double* opt, int* n, int* sr, double* out);
void PEMbuildC(int* ne, int* nsp, double* Bc, double* m, double* d, double* a, double* psi, double* w, double* BcW);
void PEMupdateC(int* ne, int* nsp, double* Bc, double* d, double* a, double* psi, double* w, double* BcW);
void PEMLoc2Scores(int* ne, double* mw, int* ntgt, double* loc, double* a, double* psi, int* nd, double* d, double* vt, double* sc);
