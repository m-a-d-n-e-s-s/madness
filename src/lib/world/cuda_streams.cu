#include <cuda.h>
//#include <world/cuda_streams.h>
#include <cuda_runtime.h>
//#define NUM_STREAMS 16
//cudaStream_t streams[NUM_STREAMS];

extern "C" void ** streams_initialize(unsigned int streams){
    unsigned int i;
    void ** cast_streams = new void*[streams];
    for (i = 0; i < streams; i++){
        cudaStream_t * stream = new cudaStream_t;
        cudaStreamCreate(stream);
        cast_streams[i] = (void *)stream;
    }
    return cast_streams;
}


extern "C" void streams_destroy(void ** cast_streams, unsigned int streams){
    unsigned int i;
    for (i = 0; i < streams; i++){
        cudaStream_t * stream = (cudaStream_t *)cast_streams[i];
        cudaStreamDestroy(*stream);
        delete stream;
    }
}

extern "C" void streams_synchronize(void ** cast_streams, unsigned int streams){
    unsigned int i;
    for (i = 0; i < streams; i++){
        cudaStream_t * stream = (cudaStream_t *)cast_streams[i];
        cudaStreamSynchronize(*stream);
    }

}
