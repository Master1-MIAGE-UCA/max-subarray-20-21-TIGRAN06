#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
struct tablo {
    int * tab;
    int size;
};
#define MINGLOBAL -2147483640
int max(int a, int b){
    if (a>b){
        return a;
    }else{
        return b;
    }
}
void printArray(struct tablo * tmp) {
    printf("---- Array of size %i ---- \n", tmp->size);
    int size = tmp->size;
    int i;
    for (i = 0; i < size; ++i) {
        printf("%i ", tmp->tab[i]);
    }
    printf("\n");
}
void montee(struct tablo * source, struct tablo * destination) {
    destination->tab[0] = 0;
    for (int p = 0; p < source->size; p++) {
        destination->tab[p + source->size] = source->tab[p];
    }
    for (int l = log2(source->size) - 1; l >= 0; l--) {
        int inf = pow(2, l);
        int sup = pow(2, l + 1);
        #pragma omp parallel for
        for (int j = inf; j < sup; j++) {
            destination->tab[j] = destination->tab[2 * j] + destination->tab[2 * j + 1];
        }
    }
}
void monteeMAX(struct tablo * source, struct tablo * destination) {
    destination->tab[0] = MINGLOBAL;
    for (int p = 0; p < source->size; p++) {
        destination->tab[p + source->size] = source->tab[p];
    }

    for (int l = log2(source->size) - 1; l >= 0; l--) {
        int inf = pow(2, l);
        int sup = pow(2, l + 1);
        #pragma omp parallel for
        for (int j = inf; j < sup; j++) {
            destination->tab[j] = max(destination->tab[2 * j] , destination->tab[2 * j + 1]);
        }
    }

}
void descenteSuffixe(struct tablo * a, struct tablo * b) {
    b->tab[1]=0;
    for(int l=1;l<=log2(a->size/2);l++){
        int inf = pow(2,l);
        int sup = pow(2,l+1);
        #pragma omp parallel for
        for(int j = sup-1;j>=inf;j--){
            if(j%2==1){
                b->tab[j]=b->tab[(j-1)/2];
            }
            else{b->tab[j]=b->tab[j/2]+a->tab[j+1];}
        }

    }
}
void descenteSuffixeMAX(struct tablo * a, struct tablo * b) {
    b->tab[1]=MINGLOBAL;
    for(int l=1;l<=log2(a->size/2);l++){
        int inf = pow(2,l);
        int sup = pow(2,l+1);
        #pragma omp parallel for
        for(int j = sup-1;j>=inf;j--){
            if(j%2==1){
                b->tab[j]=b->tab[(j-1)/2];
            }
            else{b->tab[j]=max(b->tab[j/2],a->tab[j+1]);}
        }

    }
}
void descente(struct tablo * a, struct tablo * c) {
    c->tab[1]=0;
    for(int l=1;l<=log2(a->size/2);l++){
        int inf = pow(2,l);
        int sup = pow(2,l+1);
        #pragma omp parallel for
        for(int j = inf;j<sup;j++){
            if(j%2==0){
                c->tab[j]=c->tab[j / 2];
            }
            else{ c->tab[j]= c->tab[(j - 1) / 2] + a->tab[j - 1];}
        }

    }

}
void descenteMAX(struct tablo * a, struct tablo * c) {
    c->tab[1]=MINGLOBAL;
    for(int l=1;l<=log2(a->size/2);l++){
        int inf = pow(2,l);
        int sup = pow(2,l+1);
        #pragma omp parallel for
        for(int j = inf;j<sup;j++){
            if(j%2==0){
                c->tab[j]=c->tab[j / 2];
            }
            else{ c->tab[j]= max(c->tab[(j - 1) / 2] , a->tab[j - 1]);}
        }

    }

}
void final(struct tablo * a, struct tablo *b) {
    int inf = pow(2,log2(a->size/2));
    int sup = pow(2,log2(a->size/2)+1);

    #pragma omp parallel for
    for(int j=inf;j<sup;j++){
        b->tab[j]=b->tab[j]+a->tab[j];
    }
}
void finalMAX(struct tablo * a, struct tablo *b) {
    int inf = pow(2,log2(a->size/2));
    int sup = pow(2,log2(a->size/2)+1);

    #pragma omp parallel for
    for(int j=inf;j<sup;j++){
        b->tab[j]=max(b->tab[j],a->tab[j]);
    }
}
void generateArray(struct tablo * s,char * threads,int size){
    s->size=size;
    s->tab=malloc(s->size*sizeof(int));
    FILE *myFile;
    myFile = fopen(threads, "r");
    int i;
    for (i = 0; i < size; i++)
    {
        fscanf(myFile, "%d", &s->tab[i]);
    }
}
void initializeValgrind(struct tablo * dest){
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = 0;
    }

}
void inititalizeValgrindMax(struct tablo * dest){
    #pragma omp parallel for
    for (int i = 0; i < dest->size; i++) {
        dest->tab[i] = MINGLOBAL;
    }

}
void MakeArray(struct tablo * source, struct tablo * dest){
    dest->size =source->size;
    dest->tab=malloc(dest->size * sizeof(int));
    initializeValgrind(dest);
}
void MakeArrayMAX(struct tablo * source, struct tablo * dest){
    dest->size =source->size;
    dest->tab=malloc(dest->size * sizeof(int));
    inititalizeValgrindMax(dest);
}
void scan_prefixe(struct tablo *source, struct tablo *dest){
    struct tablo * a = malloc(sizeof(struct tablo));
    a->tab = malloc(source->size*2*sizeof(int));
    a->size =source->size*2;
    initializeValgrind(a);
    montee(source, a);
    struct tablo * b = malloc(sizeof(struct tablo));
    b->tab= malloc(source->size*2*sizeof(int));
    b->size=source->size*2;
    initializeValgrind(b);
    descente(a, b);
    final(a,b);
    #pragma omp parallel for
    for (int i=0; i < dest->size; i++){
        dest->tab[i]=b->tab[i + source->size];
    }
    free(a->tab);
    free(b->tab);
    free(a);
    free(b);
}
void scan_suffixe(struct tablo *source, struct tablo *dest){
    struct tablo * a = malloc(sizeof(struct tablo));
    a->tab = malloc(source->size*2*sizeof(int));
    a->size =source->size*2;
    initializeValgrind(a);
    montee(source, a);
    struct tablo * b = malloc(sizeof(struct tablo));
    b->tab= malloc(source->size*2*sizeof(int));
    b->size=source->size*2;
    initializeValgrind(b);
    descenteSuffixe(a, b);
    final(a,b);

    #pragma omp parallel for
    for (int i=0; i < dest->size; i++){
        dest->tab[i]=b->tab[i + source->size];
    }
    free(a->tab);
    free(b->tab);
    free(a);
    free(b);
}
void scan_prefixeMAX(struct tablo *source, struct tablo *dest){
    struct tablo * a = malloc(sizeof(struct tablo));
    a->tab = malloc(source->size*2*sizeof(int));
    a->size =source->size*2;
    inititalizeValgrindMax(a);
    monteeMAX(source, a);

    struct tablo * b = malloc(sizeof(struct tablo));
    b->tab= malloc(source->size*2*sizeof(int));
    b->size=source->size*2;
    inititalizeValgrindMax(b);
    descenteMAX(a, b);
    finalMAX(a,b);
    #pragma omp parallel for
    for (int i=0; i < dest->size; i++){
        dest->tab[i]=b->tab[i +  source->size];
    }
    free(a->tab);
    free(b->tab);
    free(a);
    free(b);
}
void scan_suffixeMAX(struct tablo *source, struct tablo *dest){
    struct tablo * a = malloc(sizeof(struct tablo));
    a->tab = malloc(source->size*2*sizeof(int));
    a->size =source->size*2;
    inititalizeValgrindMax(a);
    monteeMAX(source, a);

    struct tablo * b = malloc(sizeof(struct tablo));
    b->tab= malloc(source->size*2*sizeof(int));
    b->size=source->size*2;
    inititalizeValgrindMax(b);
    descenteSuffixeMAX(a, b);

    finalMAX(a,b);
    #pragma omp parallel for
    for (int i=0; i < dest->size; i++){
        dest->tab[i]=b->tab[i +  source->size];
    }
    free(a->tab);
    free(b->tab);
    free(a);
    free(b);
}
int countInt(char * threads){
    FILE* file;
    int array[100000];
    int count = 0;
    file = fopen(threads,"r");

    while(fscanf(file, "%d", &array[count]) == 1)
        count++;

    return count;
}

int main(int argc, char **argv){
    char * threads ="";
    if (argc>1) {
        threads= argv[1];}
    else {printf("WARNING No input parameter exit status 1");
            exit(1);}

    struct tablo Q;
    FILE* fileTest = fopen(threads, "r");
    fseek(fileTest, 0L, SEEK_END);
    Q.size = ftell(fileTest);
    fseek(fileTest, 0L, SEEK_SET);
    int count = 0;
    Q.tab = malloc(Q.size * sizeof(int));
    while (fscanf(fileTest, "%d", &Q.tab[count]) == 1) {
        count++;
    }
    fclose(fileTest);
    Q.size = count;
    //printArray(&Q);
    struct tablo * PSUM = malloc(sizeof(struct tablo));
    MakeArray(&Q, PSUM);
    scan_prefixe(&Q, PSUM);
    //printArray(PSUM);

    struct tablo * SSUM = malloc(sizeof(struct tablo));
    MakeArray(&Q, SSUM);
    scan_suffixe(&Q, SSUM);
    //printArray(SSUM);

    struct tablo * SMAX = malloc(sizeof(struct tablo));
    MakeArrayMAX(&Q, SMAX);
    scan_suffixeMAX(PSUM, SMAX);
    //printArray(SMAX);

    struct tablo * PMAX = malloc(sizeof(struct tablo));
    MakeArrayMAX(&Q, PMAX);
    scan_prefixeMAX(SSUM, PMAX);
    //printArray(PMAX);

    struct tablo * M = malloc(sizeof(struct tablo));
    MakeArray(&Q, M);

    #pragma omp parallel for
    for (int i=0; i < Q.size; i++) {
        M->tab[i] = SMAX->tab[i] - PSUM->tab[i] + Q.tab[i] + PMAX->tab[i] - SSUM->tab[i];
    }
    //printArray(M);
    int MaxGlobalValue = 0;
    int MaxGlobalIndex =0;
    #pragma omp parallel
    {
        int MaxValue = 0;
        int index = 0;
        #pragma omp critical
        for (int indexIterator = 0; indexIterator < M->size; ++indexIterator) {
            if (M->tab[indexIterator] > MaxValue) {
                MaxValue = M->tab[indexIterator];
                index = indexIterator;
            }
        }
        {
            if(MaxValue >= MaxGlobalValue) {
                if((MaxValue==MaxGlobalValue) & (index<MaxGlobalIndex)){
                    MaxGlobalIndex = index;
                }
                else{
                    MaxGlobalValue = MaxValue;
                    MaxGlobalIndex = index;}
            }

        }
    }

    int pointeur= MaxGlobalIndex;
    int pointeur2= MaxGlobalIndex;
    while(pointeur>0){
        if (!(M->tab[pointeur] == MaxGlobalValue)){
            pointeur++;
            break;
        }
        else {
            if (pointeur==-1){
                pointeur++;
                break;
            }
            else pointeur--;
        }
    }

    while(pointeur2 < M->size){
        if (!(M->tab[pointeur2] == MaxGlobalValue)){
            pointeur2--;
            break;
        }
        else {
            if (pointeur2 == M->size-1){
                break;
            }
            else pointeur2++;
        }
    }
    int MinIterator = pointeur;
    int MaxIterator = pointeur2;
    printf("\n%d ",MaxGlobalValue);
    for (int i=MinIterator;i<MaxIterator;i++){
        printf("%d ", Q.tab[i]);
    }
    printf("%d", Q.tab[MaxIterator]);
    free(Q.tab);
    free(M->tab);
    free(PSUM->tab);
    free(SSUM->tab);
    free(PMAX->tab);
    free(SMAX->tab);
    free(M);
    free(PSUM);
    free(SSUM);
    free(PMAX);
    free(SMAX);
}
