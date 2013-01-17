#include<stdlib.h>
#include<stdio.h>
#include<omp.h>
#include<stdint.h>

#define DATAT uint64_t
#define SIZET long int 
#define SIZE (1 << 27)
#define SEGMENTS 8

DATAT* get_data(SIZET size, int which, int kind);
int check(DATAT* data, SIZET size);
DATAT* seq_merge(DATAT* d1, SIZET s1, DATAT* d2, SIZET s2, DATAT* res);
void segment(SIZET s, SIZET* res, int segments, int threads) ;
void segment_back(DATAT* d1, SIZET s1, DATAT* d2, SIZET s2, SIZET* seg1, SIZET* seg1b, int segments, int threads) ;
SIZET search(DATAT value, DATAT* data, SIZET size, SIZET offset) ;
void merge(DATAT* d1, SIZET s1, DATAT* d2, SIZET s2, DATAT* res) ;
DATAT* par_merge(DATAT * d1, SIZET s1, DATAT* d2, SIZET s2, DATAT* res, int threads) ;
void merge2(SIZET* d1, SIZET s1, SIZET* d2, SIZET s2, SIZET* res) ;
void make_testdata(DATAT* d1, SIZET s1, DATAT* d2, SIZET s2, int kind) ;

int main(int argc, char *argv[]) {
    int threads, myid;
    int i; threads = 1;
    int parallel = 1;
    int allthreads = 0;
    int testdata = 1;
    int passes = 1;
    SIZET size = SIZE;
    for (i=1; i<argc&&argv[i][0]=='-'; i++) {
        if (argv[i][1]=='t') sscanf(argv[++i],"%d",&threads);
        if (argv[i][1]=='k') sscanf(argv[++i],"%d",&testdata);
        if (argv[i][1]=='d') sscanf(argv[++i],"%ld",&size);
        if (argv[i][1]=='s') parallel=0;
        if (argv[i][1]=='a') allthreads=1;
        if (argv[i][1]=='p') sscanf(argv[++i],"%d",&passes);
    }
#ifdef DEBUG
    printf("Maximum number of threads possible is %d\n",
            omp_get_max_threads());
#endif

    if (threads<omp_get_max_threads()) {
        if (threads<1) threads = 1;
        omp_set_num_threads(threads);
    } else {
        threads = omp_get_max_threads();
    }
#ifdef DEBUG
#pragma omp parallel num_threads(threads)
    {
        myid = omp_get_thread_num();
        printf("Thread %d of %d active\n",myid,threads);
        //for(SIZET i = 0; i < SIZE; i++);
    }
#endif

#ifdef DEBUG
    printf("pre data generation\n");
#endif

    DATAT* data1 = calloc(size, sizeof(DATAT));//get_data(size, 1, testdata);
    DATAT* data2 = calloc(size, sizeof(DATAT));//get_data(size, 2, testdata);
    DATAT* res = calloc(2*size,sizeof(DATAT));

    make_testdata(data1, size, data2, size, testdata);

#ifdef DEBUG
    printf("post data generation\n");
#endif

    double t1;
    double t2;
    res = seq_merge(data1, size, data2, size, res);
    if(! check(res, size*2)) {
        printf("something is WRONG\n");
    } else {
    #ifdef DEBUG
        printf("it works!\n");
    #endif
    }
    res = par_merge(data1, size, data2, size, res, threads);
    if(! check(res, size*2)) {
        printf("something is WRONG\n");
    } else {
    #ifdef DEBUG
        printf("it works!\n");
    #endif
    }

    if(! parallel || allthreads) {
        printf("sequential, ");
        fflush(stdout);
        for(int i = passes; i > 0; i--) {
            if(testdata == 3) {
                make_testdata(data1, size, data2, size, testdata);
            }
            t1 = - omp_get_wtime();
            res = seq_merge(data1, size, data2, size, res);
            t1 += omp_get_wtime();
            printf("%lf", t1);
            fflush(stdout);
            if(i > 1) {
                printf(", ");
            }
        }
        printf("\n");
    }

    if(parallel && ! allthreads) {
        printf("%d threads, ", threads);
        fflush(stdout);
        for(int i = passes; i > 0; i--) {
            t2 = - omp_get_wtime();
            res = par_merge(data1, size, data2, size, res, threads);
            t2 += omp_get_wtime();
            printf("%lf", t2);
            fflush(stdout);
            if(i > 1) {
                printf(", ");
            }
        }
        printf("\n");
    }

    if(allthreads) {
        for(int i = 1; i <= threads; i++) {
            printf("%d threads, ", i);
            fflush(stdout);
            for(int j = passes; j > 0; j--) {
                if(testdata == 3) {
                    make_testdata(data1, size, data2, size, testdata);
                }
                t2 = - omp_get_wtime();
                res = par_merge(data1, size, data2, size, res, i);
                t2 += omp_get_wtime();
                printf("%lf", t2);
                fflush(stdout);
                if(j > 1) {
                    printf(", ");
                }
            }
            printf("\n");
        }
        
    }

#ifdef DEBUG
    if(data1[search(1336,data1,size,0)] == 1336)
        printf("woo\n");
    else
        printf("boo %ld\n",data1[search(1336,data1,size,0)]);
    printf("bla %ld %ld\n", data1[search(253522375,data1,size,0)], data1[size-1]);
#endif

    //printf("%lf\n%lf\n", t1, t2);
    return 0;
}

void make_testdata(DATAT* d1, SIZET s1, DATAT* d2, SIZET s2, int kind) {
    SIZET i = 0;
    SIZET j = 0;
    switch(kind) {
    case 1: //interleaved
        for(i = 0; i < s1; i++) {
            d1[i] = i + j;
            if(j < s2) {
                d2[j] = i + j + 1;
                j++;
            }
        }
        while(j < s2) {
            d2[j] = i + j;
            j++;
        }
    break;
    case 2: //1 before 2
        for(i = 0; i < s1; i++) {
            d1[i] = i;
        }
        for(j = 0; j < s2; j++) {
            d2[j] = i + j;
        }
    break;
    case 3: //RANDOM!!!11
        while(i < s1 && j < s2) {
            if (rand() <= RAND_MAX / 2) {
                d1[i] = i+j;
                i++;
            } else {
                d2[j] = i+j;
                j++;
            }
        }
        while(i < s1) {
            d1[i] = i+j;
            i++;
        }
        while(j < s2) {
            d2[j] = i+j;
            j++;
        }
    break;

    }
}

DATAT* get_data(SIZET size, int which, int kind) {
    DATAT* ret = malloc(size*sizeof(DATAT));

    switch(kind) {
    case 1: //interleaved
        for(SIZET i = 0; i < size; i++) {
            ret[i] = 2*i + which;
        }
        break;
    case 2: //1 before 2
        for(SIZET i = 0; i < size; i++) {
            ret[i] = i + (which-1)*size;
        }
        break;
    }
    return ret;
}

int check(DATAT* data, SIZET size) {
    int ret = 1;

    for(SIZET i = 0; ret && i < size - 1; i++) {
        //ret = data[i] <= data[i+1];
        ret = data[i] == i;
        if(! ret)
            printf("oh noes! %ld doesn't fit\n", data[i+1]);
    }

    return ret;
}

DATAT* seq_merge(DATAT* d1, SIZET s1, DATAT* d2, SIZET s2, DATAT* res) {

    //printf("sequential merge\n");

    merge(d1,s1,d2,s2,res);

    return res;
}

void merge(DATAT* d1, SIZET s1, DATAT* d2, SIZET s2, DATAT* res) {
    SIZET i1 = 0;
    SIZET i2 = 0;

#ifdef DEBUG
    printf("merging (%ld, %ld)[%ld] and (%ld, %ld)[%ld]\n", d1[0], d1[s1-1], s1, d2[0], d2[s2-1], s2);
#endif

    while(i1 < s1 && i2 < s2) {
        if(d1[i1] < d2[i2]) {
            res[i1 + i2] = d1[i1]; i1++;
        } else {
            res[i1 + i2] = d2[i2]; i2++;
        }
    }

    while(i1 < s1) {
        res[i1 + i2] = d1[i1]; i1++;
    }
    while(i2 < s2) {
        res[i1 + i2] = d2[i2]; i2++;
    }
}
void merge2(SIZET* d1, SIZET s1, SIZET* d2, SIZET s2, SIZET* res) {
    SIZET i1 = 0;
    SIZET i2 = 0;

    while(i1 < s1 && i2 < s2) {
        if(d1[i1] < d2[i2]) {
            res[i1 + i2] = d1[i1]; i1++;
        } else {
            res[i1 + i2] = d2[i2]; i2++;
        }
    }

    while(i1 < s1) {
        res[i1 + i2] = d1[i1]; i1++;
    }
    while(i2 < s2) {
        res[i1 + i2] = d2[i2]; i2++;
    }
}

DATAT* par_merge(DATAT * d1, SIZET s1, DATAT* d2, SIZET s2, DATAT* res, int threads) {
    SIZET segments = threads*2;

    SIZET seg1[segments];
    SIZET seg1b[segments];
    SIZET seg2[segments];
    SIZET seg2b[segments];
    SIZET segc1[segments*2+2];
    SIZET segc2[segments*2+2];

    //printf("parallel merge\n");

    segc1[0] = 0;
    segc2[0] = 0;
    segc1[segments*2+1] = s1;
    segc2[segments*2+1] = s2;

    segment(s1, seg1, segments, threads);
    segment(s2, seg2, segments, threads);
#ifdef DEBUG
    for(int i = 0; i < segments; i++) {
        printf("%ld ", seg1[i]);
    }
    printf("\n");
    for(int i = 0; i < segments; i++) {
        printf("%ld ", seg2[i]);
    }
    printf("\n");
#endif
    segment_back(d1, s1, d2, s2, seg2, seg1b, segments, threads);
    segment_back(d2, s2, d1, s1, seg1, seg2b, segments, threads);
    
#ifdef DEBUG
    for(int i = 0; i < segments; i++) {
        printf("%ld ", seg1[i]);
    }
    printf("\n");
    for(int i = 0; i < segments; i++) {
        printf("%ld ", seg2[i]);
    }
    printf("\n");
    for(int i = 0; i < segments; i++) {
        printf("%ld ", seg1b[i]);
    }
    printf("\n");
    for(int i = 0; i < segments; i++) {
        printf("%ld ", seg2b[i]);
    }
    printf("\n");
#endif

    merge2(seg1, segments, seg1b, segments, segc1+1);
    merge2(seg2, segments, seg2b, segments, segc2+1);

#ifdef DEBUG
    for(int i = 0; i < segments*2+2; i++) {
        printf("%ld ", segc1[i]);
    }
    printf("\n");
    for(int i = 0; i < segments*2+2; i++) {
        printf("%ld ", segc2[i]);
    }
    printf("\n");

    for(int i = 0; i < segments*2+2; i++) {
        printf("%ld ", d1[segc1[i]]);
    }
    printf("\n");
    for(int i = 0; i < segments*2+2; i++) {
        printf("%ld ", d2[segc2[i]]);
    }
    printf("\n");
#endif

    #pragma omp parallel for num_threads(threads)
    for(SIZET i = 0; i < segments*2+1; i++) {
       merge(d1+segc1[i], segc1[i+1] - segc1[i], d2+segc2[i], segc2[i+1] - segc2[i], res + segc1[i]+segc2[i]); 
    }
    
    return res;
}

void segment(SIZET s, SIZET* res, int segments, int threads) {
    #pragma omp parallel for num_threads(threads)
    for(int i = 0; i < segments; i++) {
        res[i] = ((i+1)*s)/(segments+1);
    }
}

void segment_back(DATAT* d1, SIZET s1, DATAT* d2, SIZET s2, SIZET* seg2, SIZET* seg1b, int segments, int threads) {
    #pragma omp parallel for num_threads(threads)
    for(int i = 0; i < segments; i++) {
        int tosearch = d2[seg2[i]];
        seg1b[i] = search(tosearch, d1, s1, 0);
    }
}

SIZET search(DATAT value, DATAT* data, SIZET size, SIZET offset) {
    if(size==1)
    {
        if(data[offset] >= value)
            return offset;
        else
            return offset+1;
    }
    if(size==0)
    {
        return offset;
    }
    if(data[offset + size/2] == value) {
        return offset + size/2;
    } else if(data[offset + size/2] < value) {
        return search(value, data, (size-1)/2, offset + size/2 + 1);
    } else {
        return search(value, data, size/2, offset);
    }
    return offset;
}
