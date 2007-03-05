#ifndef WORLD_DEP_H
#define WORLD_DEP_H

/// \file worlddep.h
/// \brief Defines DependencyInterface and CallbackInterface

namespace madness {

    /// This class used for callbacks (e.g., for dependency tracking)
    class CallbackInterface {
    public:
        virtual void notify() = 0;

	virtual ~CallbackInterface(){};
    };
        

    /// Provides interface for tracking dependencies
    class DependencyInterface : public CallbackInterface {
    private:
        int ndepend;        //< Counts dependencies
        
        mutable std::vector<CallbackInterface*> callbacks; //< Called ONCE by dec() ndepend==0

        inline void do_callbacks() const {
            for (int i=0; i<(int)callbacks.size(); i++)
                callbacks[i]->notify();
            callbacks.clear();
        };

    public:
        DependencyInterface(int ndepend = 0) : 
            ndepend(ndepend), callbacks() {};

        /// Returns the number of unsatisfied dependencies
        inline int ndep() const {
            return ndepend;
        };

        /// Returns true if ndepend == 0
        inline bool probe() const {return ndep() == 0;};


        /// Invoked by callbacks to notifiy of dependencies being satisfied
        void notify() {dec();};


        /// Registers a callback for when ndepend=0

        /// If ndepend == 0, the callback is immediately invoked.
        inline void register_callback(CallbackInterface* callback) {
            this->callbacks.push_back(callback);
            if (ndepend == 0) do_callbacks();
        };
        
        /// Increment the number of dependencies
        inline void inc() {
            MADNESS_ASSERT(ndepend>=0);
            ndepend++;
        };
        
        /// Decrement the number of dependencies and invoke callback if ndepend=0
        inline void dec() {
            ndepend--;
            MADNESS_ASSERT(ndepend>=0);
            if (ndepend == 0) do_callbacks();
        };

        virtual ~DependencyInterface() {
            if (ndepend) {
                print("DependencyInterface: destructor with ndepend =",ndepend,"?");
            }
        };
    };
}
#endif
