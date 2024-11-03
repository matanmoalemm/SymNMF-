# ifndef SYMNMF_H_
# define SYMNMF_H_


int d,n;

typedef signed long ssize_t;

typedef struct Node {
    double *data;
    struct Node* next;
} Node;

typedef struct LinkedList{
    Node* head;
    Node* tail;
    int size;
}LinkedList;

double EDistance_squared(double *point1, double *point2,int d);
LinkedList** assign(LinkedList* data,double **currClusters,int d,int K);
void compute_d(char* path);
double** readData(char* path);
void freeArray(double** array, int K);
LinkedList* createLinkedList();
Node* createNode(double *data,int d);
void append(LinkedList* list, double *data,int d);
void printlinkedList(LinkedList* list, int d);
void freeList(LinkedList* list);
int handleMemoryFail();
void printArray(double** list,int s, int m);
double** matrix_mul(double **firstMatrix,double **secondMatrix, int row1,int col1,int row2, int col2);
double** sym(double **data);
double** ddg(double **A);
double** norm(double** D,double** A);
int convergence(double** prevH,double** newH,int k);
double** transpose(double** matrix,int m , int k );
double** updateH(double** H,double** W,int k);
double** finalH(double** H,double** W,int k);
void printArray(double** list,int s, int m);
double** copyH(double** H,int row,int col);
double** createArray(int row, int col);
# endif