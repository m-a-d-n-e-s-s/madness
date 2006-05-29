#ifndef SAV_H_
#define SAV_H_

/// \file sav.h
/// \brief Basic functionality for single assignment variables

/// How to think about using the task pool?  You can a insert custom
/// task which derives from TaskInterface but does not use the single
/// assignment variables.  Or you can use one of the existing
/// task templates which require you to use SAV as the arguments.
/// In addition to providing the necessary interface, the SAV has
/// the correct and efficient copy semantics (since it is actually
/// just a shared pointer).

#include <misc/shared_ptr.h>


template <class T>
class SAVInterface {
public:
	virtual void set(const T& t) = 0;
	virtual const T& get() const = 0;
	virtual bool probe() const = 0;
};

template <typename T>
class SAVImpl : public SAVInterface<T> {
private:
    mutable bool assigned;
    T t;  ///< The actual value

protected:
    /// Friends may need to assign via a pointer
    inline T* get_ptr() const {
        return const_cast<T*>(&t);
    };

    /// Friends assigning via a pointer will need to mark assignment

    /// This also calls assign_action.  Logically const since may be
    /// called from probe_action and it easier for us to deconst here.
    inline void set_assigned() const {
        assigned = true;
        assign_action();
    };

    /// Reset ... is necessary in some scenarios?
    inline void reset() {
        assigned = false;
    };

public:
    /// Default constructor is an unassigned variable
    SAVImpl()
    : assigned(false)
    , t() {};

    /// Constructor assigning from value
    SAVImpl(const T& t)
    : assigned(true)
    , t(t) {};

    /// Probe returns true if assigned & must not block
    virtual bool probe() const {    // Logically const ...
        if (!assigned) probe_action();
        return assigned;
    };

    /// Get const reference to the value ... throws exception if not assigned
    inline const T& get() const {
            if (!probe()) throw "SAV: trying to get unassigned SAV";
            return t;
        };

    /// Set the value ... throws exception if assigned
    inline void set(const T& value) {
        if (probe()) throw "SAV: trying to set assigned SAV";
        t = value;
        set_assigned();
    };

    /// Override this to have an action performed at assignment
    virtual void assign_action() const {};

    /// Override this to have an action performed at probe
    virtual void probe_action() const {};

    virtual ~SAVImpl() {};
};

/// Public wrapper for single assignment variables

/// In addition to easy memory management and privacy,
/// the wrapping of the implementation in a shared pointer
/// makes this quick to copy.  It is also safe since behind the
/// scenes there is only one instance of the variable
/// which ensures the correct semantics.  
template <class T, class Q = SAVImpl<T> >
class SAV : public SAVInterface<T> {
public:
	typedef madness::SharedPtr<Q> ptrT;

private:
    ptrT p;

public:
    SAV() : p(new Q()) {};
    
    SAV(const T& t) : p(new Q(t)) {};

    inline void set(const T& value) {
        p->set(value);
    };

    inline const T& get() const {
        return p->get();
    };

    inline bool probe() const {
        return p->probe();
    };
    
    virtual ~SAV() {};
};


#endif /*SAV_H_*/
