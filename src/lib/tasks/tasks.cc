#include <tasks/tasks.h>

namespace madness {
    TaskQueue taskq;

    void task_add_am(am_handlerT op, ProcessID src, const AMArg& arg) {
        taskq.add_local_hp(new TaskAM(src,op,arg));       
    };
    
    /// Might be beneficial to have this done directly in the AM handler.
    void task_generic_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
        long size = arg.arg0;
        long tag = arg.arg1;
        long handle = arg.arg2;
        //print(comm.rank(),"in TGH",size,tag,handle);
        GenericOpT op = taskq.get_generic_pointer(handle);
        // For now, remotely created tasks are always high priority
        taskq.add_local_hp(new TaskGeneric(size,src,tag,op));
        //print(comm.rank(),"added task in TGH");
    }     
    
    /// Handler for generic broadcast task creation
    void task_generic_broadcast_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
        long size = arg.arg0;
        long tag = arg.arg1;
        long handle = arg.arg2;
        long root = arg.arg3;
        std::vector<unsigned char> v(size);
        MPI::Request req;
        //print(comm.rank(),"in BGH",size,tag,handle,src,root);
        if (size>0) madness::comm_default->Recv(&v[0],size,src,tag);
        //print((int)v[0],(int)v[1],(int)v[2],(int)v[3],(int)v[4],(int)v[5],(int)v[6],(int)v[7]);
        GenericOpT op = taskq.get_generic_pointer(handle);
        taskq.add_generic_broadcast(handle,v,root);
        // For now, remotely created tasks are always high priority
        taskq.add_local_hp(new TaskGeneric(root,op,v));
        //print(comm.rank(),"added task in BGH");
    }     
    
    
    
}
