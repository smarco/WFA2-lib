#pragma once

#define CUDACHK(code) { gpuAssert((code), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

bool canPrefetch(int dstDevice){
   if(dstDevice==cudaCpuDeviceId) dstDevice=0;
   int res;
   CUDACHK(cudaDeviceGetAttribute(&res, cudaDevAttrConcurrentManagedAccess, dstDevice));
   return res;
}

int smCount(int device){
   cudaDeviceProp prop;
   CUDACHK(cudaGetDeviceProperties(&prop, device));
   return prop.multiProcessorCount;
}

void smemCarveout(int percent, void* func){
   CUDACHK(cudaFuncSetAttribute(func, cudaFuncAttributePreferredSharedMemoryCarveout, percent));
}

static bool alreadySetMallocHeapLimit = false;
void setMallocHeapLimit(size_t bytes){
   if(alreadySetMallocHeapLimit) return;
   CUDACHK(cudaDeviceSetLimit(cudaLimitMallocHeapSize, bytes));
   alreadySetMallocHeapLimit = true;
}

int maximizeDynamicSmem(void* func, int device){
   int max_smem;
   CUDACHK(cudaDeviceGetAttribute(&max_smem, cudaDevAttrMaxSharedMemoryPerBlockOptin, device));
   cudaFuncAttributes attr;
   CUDACHK(cudaFuncGetAttributes(&attr, func));
   int static_smem = attr.sharedSizeBytes;
   int max_dynamic_smem = max_smem - static_smem;
   CUDACHK(cudaFuncSetAttribute(func, cudaFuncAttributeMaxDynamicSharedMemorySize, max_dynamic_smem));
   return max_dynamic_smem;
}
