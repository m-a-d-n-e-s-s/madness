/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/

/*!
  \file examples/znemo.cc
  \brief solve the HF equations using numerical exponential MOs

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/nemo.cc>here</a>.

*/

#include <madness/chem/zcis.h>

using namespace madness;


int main(int argc, char** argv) {

    World& world=initialize(argc, argv,false);
    if (world.rank() == 0) {
        print_header1(" ZCIS -- excited states in the CIS approximation using complex orbitals");
    }

    startup(world,argc,argv,true);
    std::cout.precision(6);
    if (world.rank()==0) print(info::print_revision_information());


    commandlineparser parser(argc,argv);
    if (parser.key_exists("help")) {
        Zcis::help();

    } else if (parser.key_exists("print_parameters")) {
        Zcis::print_parameters();

    } else {

        try {
            std::shared_ptr<Znemo> znemo(new Znemo(world, parser));
            znemo->value();
            Zcis zcis(world, parser, znemo);
            zcis.value();
        } catch (const SafeMPI::Exception& e) {
            print(e);
            error("caught an MPI exception");
        } catch (const madness::MadnessException& e) {
            print(e);
            error("caught a MADNESS exception");
        } catch (const madness::TensorException& e) {
            print(e);
            error("caught a Tensor exception");
        } catch (const char *s) {
            print(s);
            error("caught a string exception");
        } catch (const std::string& s) {
            print(s);
            error("caught a string (class) exception");
        } catch (const std::exception& e) {
            print(e.what());
            error("caught an STL exception");
        } catch (...) {
            error("caught unhandled exception");
        }
    }


    finalize();
    return 0;
}
