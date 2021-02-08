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
//Descente suffix de la meme maniere que prefix sauf qu'on utilise la symétrie au lieu de faire un reverse de façon naif
//
//Du coup on descends dans la arbre si on est le fils droit on prends la valeur du pere sinon on prends la valeur du pere + la valeur de son frere dans l'autre arbre
//
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
}//Descente suffix de la meme maniere que prefix sauf qu'on utilise la symétrie au lieu de faire un reverse de façon naif
//
//Du coup on descends dans la arbre si on est le fils droit on prends la valeur du pere sinon on prends le max de la valeur du pere et la valeur de son frere dans l'autre arbre
//
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

//Descente prefix
//
//Du coup on descends dans la arbre si on est le fils gauche on prends la valeur du pere sinon on prends la valeur du pere + la valeur de son frere dans l'autre arbre
//
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
//Descente prefix
//
//Du coup on descends dans la arbre si on est le fils gauche on prends la valeur du pere sinon on prends le mac de la valeur du pere et la valeur de son frere dans l'autre arbre
//
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
    for (int j=inf;j<sup;j++){
        b->tab[j]=b->tab[j]+a->tab[j];
    }
}
void finalMAX(struct tablo * a, struct tablo *b) {
    int inf = pow(2,log2(a->size/2));
    int sup = pow(2,log2(a->size/2)+1);

    #pragma omp parallel for
    for (int j=inf;j<sup;j++){
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
    free(myFile);
}

//l'algorithm Perumalla et al est vraiment rapide en parallel
void scan_prefixe(struct tablo *source, struct tablo *dest){
    struct tablo * a = malloc(sizeof(struct tablo));
    a->tab = malloc(source->size*2*sizeof(int));
    a->size =source->size*2;
    #pragma omp parallel for
    for  (int i = 0; i < a->size; i++) {
        a->tab[i] = 0;
    }
    montee(source, a);
    struct tablo * b = malloc(sizeof(struct tablo));
    b->tab= malloc(source->size*2*sizeof(int));
    b->size=source->size*2;
    #pragma omp parallel for
    for  (int i = 0; i < b->size; i++) {
        b->tab[i] = 0;
    }
    descente(a, b);
    final(a,b);
    #pragma omp parallel for
    for  (int i=0; i < dest->size; i++){
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
    #pragma omp parallel for
    for  (int i = 0; i < a->size; i++) {
        a->tab[i] = 0;
    }
    montee(source, a);
    struct tablo * b = malloc(sizeof(struct tablo));
    b->tab= malloc(source->size*2*sizeof(int));
    b->size=source->size*2;
    #pragma omp parallel for
    for  (int i=0; i < dest->size; i++){
        dest->tab[i]=b->tab[i + source->size];
    }
    descenteSuffixe(a, b);
    final(a,b);

    #pragma omp parallel for
    for  (int i=0; i < dest->size; i++){
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
    #pragma omp parallel for
    for  (int i = 0; i < a->size; i++) {
        a->tab[i] = MINGLOBAL;
    }
    monteeMAX(source, a);

    struct tablo * b = malloc(sizeof(struct tablo));
    b->tab= malloc(source->size*2*sizeof(int));
    b->size=source->size*2;
    #pragma omp parallel for
    for  (int i = 0; i < b->size; i++) {
        b->tab[i] = MINGLOBAL;
    }
    descenteMAX(a, b);
    finalMAX(a,b);
    #pragma omp parallel for
    for  (int i=0; i < dest->size; i++){
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
    #pragma omp parallel for
    for  (int i = 0; i < a->size; i++) {
        a->tab[i] = MINGLOBAL;
    }
    monteeMAX(source, a);

    struct tablo * b = malloc(sizeof(struct tablo));
    b->tab= malloc(source->size*2*sizeof(int));
    b->size=source->size*2;
    #pragma omp parallel for
    for  (int i=0; i < dest->size; i++){
        dest->tab[i]=b->tab[i +  source->size];
    }
    descenteSuffixeMAX(a, b);

    finalMAX(a,b);
    #pragma omp parallel for
    for  (int i=0; i < dest->size; i++){
        dest->tab[i]=b->tab[i +  source->size];
    }
    free(a->tab);
    free(b->tab);
    free(a);
    free(b);
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
    PSUM->size =Q.size;
    PSUM->tab=malloc(PSUM->size * sizeof(int));
    #pragma omp parallel for
    for  (int i = 0; i < PSUM->size; i++) {
        PSUM->tab[i] = 0;
    }
    scan_prefixe(&Q, PSUM);
    //printArray(PSUM);

    struct tablo * SSUM = malloc(sizeof(struct tablo));
    SSUM->size =Q.size;
    SSUM->tab=malloc(SSUM->size * sizeof(int));
    #pragma omp parallel for
    for  (int i = 0; i < SSUM->size; i++) {
        SSUM->tab[i] = 0;
    }
    scan_suffixe(&Q, SSUM);
    //printArray(SSUM);

    struct tablo * SMAX = malloc(sizeof(struct tablo));
    SMAX->size =Q.size;
    SMAX->tab=malloc(SMAX->size * sizeof(int));
    #pragma omp parallel for
    for  (int i = 0; i < SMAX->size; i++) {
        SMAX->tab[i] = 0;
    }
    scan_suffixeMAX(PSUM, SMAX);
    //printArray(SMAX);

    struct tablo * PMAX = malloc(sizeof(struct tablo));
    PMAX->size =Q.size;
    PMAX->tab=malloc(PMAX->size * sizeof(int));
    #pragma omp parallel for
    for  (int i = 0; i < SMAX->size; i++) {
        PMAX->tab[i] = 0;
    }
    scan_prefixeMAX(SSUM, PMAX);
    //printArray(PMAX);

    struct tablo * M = malloc(sizeof(struct tablo));
    M->size =Q.size;
    M->tab=malloc(M->size * sizeof(int));
    #pragma omp parallel for
    for  (int i = 0; i < M->size; i++) {
        M->tab[i] = 0;
    }

    #pragma omp parallel for
    for  (int i=0; i < Q.size; i++) {
        M->tab[i] = SMAX->tab[i] - PSUM->tab[i] + Q.tab[i] + PMAX->tab[i] - SSUM->tab[i];
    }
    //printArray(M);
    int MaxValue = MINGLOBAL;
    int index = Q.size;

    // STEP 6 with reduction
    // On cherche d'abord la valeur maximal en utilisant reduction
    #pragma omp parallel for reduction(max:MaxValue)
    for (int indexIterator = 0; indexIterator < M->size; ++indexIterator) {
        //printf("%d",omp_get_thread_num());
        if (M->tab[indexIterator] > MaxValue) {
            MaxValue = M->tab[indexIterator];
        }
    }
    //on cherche l'element dans le table ayant la meme valeur que MaxValue mais ayant l'indice minimal
    #pragma omp parallel for reduction(min:index)
    for (int indexIterator = 0; indexIterator < M->size; ++indexIterator) {
        //printf("   %d  ",index);
        if (M->tab[indexIterator] == MaxValue) {
            if(index>indexIterator){
                index=indexIterator;
            }
        }
    }
    // on se position à l'indice minimal dans le tableau M puis avec un pointeur on parcours jusqu'à que la valeur de l'element est differente de max value
    int pointeur= index;
    printf("%d",MaxValue);
    while(pointeur < M->size){
        if (!(M->tab[pointeur] == MaxValue)){
            pointeur--;
            break;
        }
        else {
            if (pointeur == M->size-1){
                printf(" %d", Q.tab[pointeur]);
                break;
            }
            else {
                if (pointeur <= M->size-1){
                    printf(" %d", Q.tab[pointeur]);
                }pointeur++;}
        }
    }


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
