#define CONCAT2(A, B) A ## B
#define CONCAT(A, B) CONCAT2(A, B)

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

#define LIST_TYPE CUDA_LIST
#define CONTAINED_TYPE CUDA_LIST_CONTAINED_TYPE
#define SUBLIST_SIZE CUDA_LIST_SUBLIST_SIZE

#define SUBLIST_TYPE CONCAT(CUDA_LIST, _sublist)
#define ITERATOR_TYPE CONCAT(CUDA_LIST, _iterator)

#define BACKING_STORAGE_BASE CONCAT(CUDA_LIST, _backing_storage_base)
#define BACKING_STORAGE_SIZE CONCAT(CUDA_LIST, _backing_storage_size)
#define NEXT_BACKING_BLOCK CONCAT(CUDA_LIST, _next_backing_block)

__managed__ char *BACKING_STORAGE_BASE;
__managed__ size_t BACKING_STORAGE_SIZE;
__managed__ char *NEXT_BACKING_BLOCK;

struct SUBLIST_TYPE;
struct ITERATOR_TYPE;

struct SUBLIST_TYPE {
    CONTAINED_TYPE content[SUBLIST_SIZE];
    SUBLIST_TYPE* __restrict__ next;
};

struct LIST_TYPE {
    SUBLIST_TYPE* __restrict__ first_sublist;
    SUBLIST_TYPE* __restrict__ current_sublist;
    int num_sublists;
    int index_in_current_sublist;

    static void backingStorageInit(size_t size);
    static void backingStorageDestruct();
    static void backingStoragePrefetch(int dstDevice);
    __device__ static SUBLIST_TYPE* claim_block();
    __device__ void init();
    __device__ void pushBack(CONTAINED_TYPE c);
    __device__ void append(CONTAINED_TYPE *cs, int num_c);
    size_t size();
    ITERATOR_TYPE begin();
    ITERATOR_TYPE end();
};

struct ITERATOR_TYPE{
    LIST_TYPE* list;
    SUBLIST_TYPE* current_sublist;
    int index_in_current_sublist;

    CONTAINED_TYPE& operator* () const{
        return current_sublist->content[index_in_current_sublist];
    };

    CONTAINED_TYPE* operator-> () const{
        return &current_sublist->content[index_in_current_sublist];
    };

    ITERATOR_TYPE operator++ (){
        index_in_current_sublist++;
        if(index_in_current_sublist == SUBLIST_SIZE){
            current_sublist = current_sublist->next;
            index_in_current_sublist = 0;
        }

        return *this;
    };

    ITERATOR_TYPE operator++ (int){
        ITERATOR_TYPE tmp = *this;
        ++*this;
        return tmp;
    };

    bool operator== (const ITERATOR_TYPE &other) const{
        bool equal = true;
        equal &= list==other.list;
        equal &= current_sublist==other.current_sublist;
        equal &= index_in_current_sublist==other.index_in_current_sublist;
        return equal;
    }

    bool operator!= (const ITERATOR_TYPE &other) const{
        return !operator==(other);
    }
};

#include <iostream>
void LIST_TYPE::backingStorageInit(size_t sublists){
    #ifdef DEBUG
        if(BACKING_STORAGE_BASE != NULL) throw logic_error("tried initializing backingStorage while it was already initialized");
    #endif

    BACKING_STORAGE_SIZE = sublists*sizeof(SUBLIST_TYPE);
    CUDACHK(cudaMallocManaged(&BACKING_STORAGE_BASE, BACKING_STORAGE_SIZE));
    NEXT_BACKING_BLOCK = BACKING_STORAGE_BASE;

};

void LIST_TYPE::backingStorageDestruct(){
    CUDACHK(cudaFree(BACKING_STORAGE_BASE));
    BACKING_STORAGE_BASE = NULL;
};

void LIST_TYPE::backingStoragePrefetch(int dstDevice){
    CUDACHK(cudaMemPrefetchAsync(BACKING_STORAGE_BASE, BACKING_STORAGE_SIZE, dstDevice));
};

__device__ SUBLIST_TYPE* LIST_TYPE::claim_block(){
    size_t addr = atomicAdd((unsigned long long *)&NEXT_BACKING_BLOCK, sizeof(SUBLIST_TYPE));
    #ifdef DEBUG
        if(addr + sizeof(SUBLIST_TYPE) > (size_t)BACKING_STORAGE_BASE + BACKING_STORAGE_SIZE){
            #ifdef DEBUG_OUTPUT
                printf(STRINGIFY(LIST_TYPE) " backing storage ran out of space %llx\n", addr);
                printf("Base:   %016llx\n", BACKING_STORAGE_BASE);
                printf("Size:   %016llx\n", BACKING_STORAGE_SIZE);
                printf("Limit:  %016llx\n", (size_t)BACKING_STORAGE_BASE + BACKING_STORAGE_SIZE);
                printf("AddrLo: %016llx\n", addr);
                printf("AddrHi: %016llx\n", addr + sizeof(SUBLIST_TYPE));
            #endif
            assert(false);
        }
    #endif
    return (SUBLIST_TYPE*)addr;
}

__device__ void LIST_TYPE::init(){
    first_sublist = claim_block();
    current_sublist = first_sublist;
    index_in_current_sublist = 0;
    num_sublists = 1;
}

__device__ void LIST_TYPE::pushBack(CONTAINED_TYPE c){
    #ifdef DEBUG
        if(index_in_current_sublist < 0 || index_in_current_sublist > SUBLIST_SIZE){
            #ifdef DEBUG_OUTPUT
                printf(STRINGIFY(LIST_TYPE) "index_in_current_sublist has invalid value %d\n", index_in_current_sublist);
            #endif
            assert(false);
        }
    #endif

    if(index_in_current_sublist == SUBLIST_SIZE){
        SUBLIST_TYPE *old_sublist = current_sublist;
        SUBLIST_TYPE *new_sublist = claim_block();
        new_sublist->next = NULL;

        old_sublist->next = new_sublist;
        current_sublist = new_sublist;
        index_in_current_sublist = 0;
        num_sublists++;
    }

    current_sublist->content[index_in_current_sublist++] = c;
}

__device__ void LIST_TYPE::append(CONTAINED_TYPE *cs, int num_c){
    for(int i=0; i<num_c; i++){
        pushBack(cs[i]);
    }
}

size_t LIST_TYPE::size(){
    size_t sum = 0;
    sum += (num_sublists-1)*SUBLIST_SIZE;
    sum += index_in_current_sublist;
    return sum;
}

ITERATOR_TYPE LIST_TYPE::begin(){
    ITERATOR_TYPE res;
    res.list = this;
    res.current_sublist = first_sublist;
    res.index_in_current_sublist = 0;
    return res;
};

ITERATOR_TYPE LIST_TYPE::end(){
    ITERATOR_TYPE res;
    res.list = this;

    if(index_in_current_sublist == SUBLIST_SIZE){
        res.current_sublist = current_sublist->next;
        res.index_in_current_sublist = 0;
    }
    else{
        res.current_sublist = current_sublist;
        res.index_in_current_sublist = index_in_current_sublist;
    }
    return res;
};

#undef CONCAT2
#undef CONCAT
#undef STRINGIFY2
#undef STRINGIFY
#undef LIST_TYPE
#undef CONTAINED_TYPE
#undef SUBLIST_SIZE
#undef SUBLIST_TYPE
#undef ITERATOR_TYPE
#undef BACKING_STORAGE_BASE
#undef BACKING_STORAGE_SIZE
#undef NEXT_BACKING_BLOCK
#undef CUDA_LIST
#undef CUDA_LIST_SUBLIST_SIZE
#undef CUDA_LIST_CONTAINED_TYPE
