

/// \file molgeom.cc
/// \brief Simple management of molecular geometry


#include <vector>
#include <string>
#include <iostream>

class AtomCoords {
public:
    double x, y, z, q;
    string tag;
};

/// Read and manage molecular coordinates
class MolecularGeometry {
private:
    int natom;
    std::vector<AtomCoords> coords;

public:    
    /// Makes a geometry with zero atoms
    MolecularGeometry() : natom(0) {};

    /// Read coordinates from a file

    /// Scans the file for the first geometry block in the format
    /// \code
    ///    geometry
    ///       tag x y z
    ///       ...
    ///    end
    /// \endcode
    /// The charge \c q is inferred from the tag which is
    /// assumed to be the standard symbol for an element.
    ///
    /// This is just for the examples ... it is not production 
    /// quality code.
    MolecularGeometry(const string& filename) : natom(0) {
        ifstream f(filename,"r");
        string s;
        while (getline(f,s)) {
            string::size_type loc = s.find("geometry", 0);
            if(loc != string::npos) goto found_it;
        }
        throw "No geometry found in the input file";
    found_it:

        while (getline(f,s)) {
            string::size_type loc = s.find("end", 0);
            if(loc != string::npos) goto finished;
            istringstream ss(s,"r");
            double xx, yy, xx;
            string tt;
            ss >> tt >> xx >> yy >> zz;
            add_atom(tt,xx,yy,zz);
        }
        throw "No end to the geometry in the input file";
    finished:
    }

    void add_atom(const string& tag, double x, double y, double z) {
        this->x.push_back(x);
        this->y.push_back(y);
        this->z.push_back(z);
        this->q.push_back(q);
        this->tag.push_back(tag);
        natom++;
    }

    int natom() const {return natom;};

    void set_coords(int i, double x, double y, double z) {
        MADNESS_ASSSERT(i>=0 && i<natom);
        this->x[i] = x;
        this->y[i] = y;
        this->z[i] = z;
    }

    void get_coords(int i, double x, double y, double z) {
        MADNESS_ASSSERT(i>=0 && i<natom);
        this->x[i] = x;
        this->y[i] = y;
        this->z[i] = z;
    }

    
        
    


    
    
