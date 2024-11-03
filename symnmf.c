#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"
#define epsilon 1e-4
#define max_iter 300

/*d is a global variable represent the dimenstion of R^d*/
void compute_d(char* path) {
    FILE *file = NULL;
    int cnt = 0;
    ssize_t i;
    size_t len = 0;
    char *line = NULL;
    char *x_i; 
    file  = fopen(path,"r");
    i = getline(&line, &len, file);
    if (i == -1) handleMemoryFail();
    d = 1; /*handels edge case of d = 1*/
    x_i = strtok(line, ",");
    while (x_i != NULL) {
        cnt++; 
        x_i = strtok(NULL, ","); 
    }
    d = cnt;
    free(line);   
    fclose(file);
}


double** readData(char* path) {
    FILE *file = NULL; ssize_t i;
    size_t len = 0;
    LinkedList* data = createLinkedList();
    double** list;
    double *point; int cnt; char *x_i =NULL ;
    char *line = NULL; int j;
    Node* node;
    compute_d(path);
    file  = fopen( path,"r");

    while ((i = getline(&line, &len, file)) != -1) {
        if (i == -1) handleMemoryFail();
        point = (double*)calloc(d,sizeof(double));    
        if (point == NULL) handleMemoryFail(); 
        cnt = 0 ;
        x_i = strtok(line, ",");
        while (x_i != NULL) {
            point[cnt] = atof(x_i);
            x_i = strtok(NULL, ",");
            cnt++;
        }  
        append(data,point,d);
        free(point);
    }
    free(line);
    fclose(file);
    n = data->size; /*n is a global variable represents the number of points*/
    list = createArray(n,d);
    node = data->head; 
    for(j = 0; j < n ; j++){
        memcpy(list[j],node->data,d*sizeof(double));
        node = node->next;
    }
    freeList(data);
    return list;
}

LinkedList* createLinkedList() {
    LinkedList* list = (LinkedList*)malloc(sizeof(LinkedList));
    if (list == NULL) handleMemoryFail();
    list->head = NULL;
    list->tail = NULL;
    list->size = 0;
    return list;
}

Node* createNode(double *data,int d) {
    double *point;
    Node* newNode = (Node*)malloc(sizeof(Node));
    if (newNode == NULL) handleMemoryFail();
    point = (double*)calloc(d,sizeof(double));
    memcpy(point,data,d*sizeof(double));
    newNode->data = point;
    newNode->next = NULL;
    return newNode;
}

void append(LinkedList* list, double *data,int d) {
    Node* newNode = createNode(data,d);
    if (list->tail == NULL) {
        list->head = newNode;
        list->tail = newNode;
        
    } else {
        list->tail->next = newNode;
        list->tail = newNode;
    }
    list->size++;
}

void freeList(LinkedList* list) {
    Node* current = list->head;
    Node* tmp;
    while (current != NULL) {
        tmp = current;
        current = current->next;
        free(tmp->data);
        free(tmp);
    }
    free(list);
}
void freeArray(double** array, int K) {
    int i;
    for ( i = 0; i < K; i++) {
        free(array[i]);
    }
    free(array);
}

int handleMemoryFail() {
    fprintf(stderr, "An Error Has Occurred\n");
    exit(1);
}

double EDistance_squared(double *point1, double *point2,int d){
    int i;
    double sum = 0.0;
    for (i=0; i<d; i++){
    sum += pow(point1[i]-point2[i],2);
    }
    return sum;

}

double** matrix_mul(double **firstMatrix,double **secondMatrix, int row1,int col1,int row2, int col2){
    int i,j,k;
    double** result = createArray(row1,col2);
    if (col1 != row2) {
        fprintf(stderr, "col1 must be equal to row2\n");
        exit(1);
    }
    for (i = 0; i < row1; i++) {
        for (j = 0; j < col2; j++) {
            for (k = 0; k < col1; k++) {
                result[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
            }
        }
    }
    return result;
}

double** sym(double **data){
    int i,j;
    double dist;
    double** A = createArray(n,n);
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            dist = -1*EDistance_squared(data[i],data[j],d);
            if (i !=j) A[i][j] = exp(dist*0.5);
            else A[i][j] = 0;
        }
    }
    return A;  
}

double** ddg(double **A){
    int i,j;
    double sum;
    double** D = createArray(n,n);
    for (i = 0; i < n; i++){
        sum = 0.0;
        for (j = 0; j < n; j++){
            sum = sum + A[i][j];
        }
        D[i][i] = sum;
        }
    return D;
}

double** norm(double** D,double** A){
    int i;
    double **W,**mul;
    for (i = 0; i < n; i++){
        D[i][i] = 1 / sqrt(D[i][i]);
    }
    mul = matrix_mul(D,A,n,n,n,n);
    W = matrix_mul(mul,D,n,n,n,n);
    freeArray(mul,n);
    return W;
}

int convergence(double** prevH,double** newH,int k){
    int i,j;
    double x = 0;
    for (i = 0 ; i < n; i++){
        for (j = 0; j < k; j++){
            x += pow(newH[i][j]-prevH[i][j],2);
        }
    }
    if (x < epsilon) return 1;
    else return 0;
}

double** transpose(double** matrix,int m , int k ){
    int i,j;
    double** Tmatrix = createArray(k,m);
    for (i = 0 ; i < k ; i++){   
        for (j = 0; j < m; j++){
            Tmatrix[i][j] = matrix[j][i];
        }
    }
    return Tmatrix;
}

/*preforms one iteration of H update*/
double** updateH(double** H,double** W,int k){
    double **mat1, **mat2, **mat3 ,**T;
    double** newH;
    int i,j;
    mat1 = matrix_mul(W,H,n,n,n,k); /*W times H*/
    T = transpose(H,n,k); /* H transpose*/
    mat3 = matrix_mul(H,T,n,k,k,n); /*H times H transposed*/     
    mat2 = matrix_mul(mat3,H,n,n,n,k);/*H times H transposed times H*/
    newH = createArray(n,n);
    for (i = 0; i <n ; i++){
        for (j = 0 ; j < k ; j++){
            newH[i][j] = H[i][j]*(1-0.5 + 0.5*(mat1[i][j]/mat2[i][j]));
        }
    }
    freeArray(mat1,n);
    freeArray(mat2,n);
    freeArray(mat3,n);
    freeArray(T,k);
    return newH;
}


double** finalH(double** H,double** W,int k){
    int iter = 1;
    double **newH ,**res;
    int is_converged = 0;
    res = copyH(H,n,k);
    while(is_converged == 0 && iter < 200){     
        newH = updateH(res,W,k);
        is_converged = convergence(res,newH,k);
        freeArray(res,n);
        res = copyH(newH,n,k);
        freeArray(newH,n);
        iter++;
    }
    return res;
}


void printArray(double** list,int s, int m){    
    int i ,j;
    for (i = 0; i <s; i++){
        for(j =0; j<m;j++){
            printf("%.4f",list[i][j]);
            if(j!=m-1) printf("%s",",");
        }
        printf("\n");
        }
}



void printlinkedList(LinkedList* list, int d){
    Node* curr = list->head;
    int i;
    while (curr != NULL){ 
        for (i = 0 ; i < d ; i++){
            printf("%.4f",curr->data[i]);
            if(i!=d-1) printf("%s",",");
        }
        printf("\n");
        curr = curr->next;
    }
    
}

double** copyH(double** H,int row,int col){
    int i,j;
    double** array = createArray(row,col);
    for (i = 0; i <row ; i++){
        for (j = 0; j < col; j++){
            array[i][j] = H[i][j];
        }
    }
    return array;
}

double** createArray(int row, int col){
    int i;
    double** array = (double**)calloc(row , sizeof(double*));
    if (array == NULL) handleMemoryFail();
    for (i = 0; i <row ; i++){
        array[i] = (double*)calloc(col , sizeof(double));
        if (array[i] == NULL) handleMemoryFail();
    }
    return array;
}



int main(int argc, char *argv[]){
    double **data,**A, **D,**W;
    if (argc != 3){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    data = readData(argv[2]);
    A = sym(data);   
    if(strcmp(argv[1],"sym") == 0) printArray(A,n,n);
    else {
        D = ddg(A);
        if(strcmp(argv[1],"ddg") == 0) printArray(D,n,n);
        else {
            if(strcmp(argv[1],"norm") == 0) {
            W = norm(D,A);
            printArray(W,n,n);
            freeArray(W,n);
            }
            else {
                printf("An Error Has Occurred\n");
                exit(1);
            }
         }
    freeArray(D,n);
    }
    freeArray(data,n);
    freeArray(A,n);
    return 1;
}

