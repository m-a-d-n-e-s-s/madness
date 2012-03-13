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


  $Id $
*/

#ifndef MADNESS_WORLD_UNIQUEID_H__INCLUDED
#define MADNESS_WORLD_UNIQUEID_H__INCLUDED

#include <cstddef>
#include <iosfwd>

namespace madness {

    class World;

    class uniqueidT {
        friend class World;
    private:
        unsigned long worldid;
        unsigned long objid;

        uniqueidT(unsigned long worldid, unsigned long objid)
                : worldid(worldid), objid(objid) {};

    public:
        uniqueidT()
                : worldid(0), objid(0) {};

        bool operator==(const uniqueidT& other) const {
            return  objid==other.objid && worldid==other.worldid;
        }

        std::size_t operator()(const uniqueidT& id) const { // for GNU hash
            return id.objid;
        }

        operator bool() const {
            return objid!=0;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & worldid & objid;
        }

        unsigned long get_world_id() const {
            return worldid;
        }

        unsigned long get_obj_id() const {
            return objid;
        }

        friend std::ostream& operator<<(std::ostream& s, const uniqueidT& id) {
            s << "{" << id.get_world_id() << "," << id.get_obj_id() << "}";
            return s;
        }
    }; // class uniqueidT


}  // namespace madness


#endif // MADNESS_WORLD_UNIQUEID_H__INCLUDED
