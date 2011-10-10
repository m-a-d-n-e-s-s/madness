#define NUM_STREAMS 16
cudaStream_t streams[NUM_STREAMS];

void streams_initialize(){
    unsigned int i;
    for (i = 0; i < NUM_STREAMS; i++){
        cudaStreamCreate(&streams[i]);
    }
}

void streams_destroy(){
    unsigned int i;
    for (i = 0; i < NUM_STREAMS; i++){
        cudaStreamDestroy(streams[i]);
    }
}

void streams_synchronize(){
    unsigned int i;
    for (i = 0; i < NUM_STREAMS; i++){
        cudaStreamSynchronize(streams[i]);
    }
}
