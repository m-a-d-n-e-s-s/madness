

/// \file molgeom.cc
/// \brief Simple management of molecular geometry


#include <vector>
#include <string>
#include <iostream>

class Atom {
public:
    double x, y, z, q;
    string tag;

    Atom(double x, double , double z, double q, const string& tag)
        : x(x), y(y), z(z), q(q), tag(tag) 
    {}
};

/// Read and manage molecular coordinates
class Molecule {
private:
    int natom;
    std::vector<Atom> atoms;

public:    
    /// Makes a molecule with zero atoms
    Molecule() : natom(0), atoms() {};

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
    /// Same as the simplest NWChem format.
    ///
    /// This code is just for the examples ... it is not production 
    /// quality.
    Molecule(const string& filename) : natom(0) {
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
        double q = tag_to_charge(tag);
        atoms.push_back(Atom(x,y,z,q,tag));
        natom++;
    }

    int get_natom() const {return natom;};

    void set_coord(int i, double x, double y, double z) {
        MADNESS_ASSSERT(i>=0 && i<natom);
        atoms[i].x = x;
        atoms[i].y = y;
        atoms[i].z = z;
    }

    const Atom& get_coords(int i) const {
        MADNESS_ASSSERT(i>=0 && i<natom);
        return atoms[i];
    }
}

    
        
    


    
    
