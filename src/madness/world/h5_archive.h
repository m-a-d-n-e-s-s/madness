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
*/

#ifndef MADNESS_WORLD_H5_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_H5_ARCHIVE_H__INCLUDED

/**
 \file h5_archive.h
 \brief Implements an archive wrapping HDF5 filestream.
 \ingroup serialization
*/

#warning 'HEELLO'

#ifdef HAVE_HDF5
#  if ! __has_include("hdf5.h")
#    error "HAVE_HDF5 is on, but hdf5.h was not found."
#  endif

#include <type_traits>
#include <fstream>
#include <cstring>
#include <madness/world/archive.h>
#include <madness/world/print.h>  // this injects operator<<(std::ostream,T) for common Ts
#include <madness/world/vector_archive.h>
#include "hdf5.h"



//#include "hdf5.h"

namespace madness {
    namespace archive {

        /// \addtogroup serialization
        /// @{

        /// Wraps an archive around a text filestream for output.
        class H5OutputArchive : public BaseOutputArchive {
            mutable std::ofstream os; ///< The filestream.
	    hid_t   file_id; /* identifiers */
	    hid_t   dataset_id; /* identifiers */
	    hid_t   dataspace_id; /* identifiers */
	    herr_t  status;
	    mutable std::vector<unsigned char> v;
            madness::archive::VectorOutputArchive var;

        public:
            /// Default constructor.

            /// The filename and open modes are optional here; they can be
            /// specified later by calling \c open().
            /// \param[in] filename Name of the file to write to.
            /// \param[in] mode I/O attributes for opening the file.
            H5OutputArchive(const char* filename = nullptr
                    ): v(), var(v)
            {
                if (filename)
                    open(filename);
            };
            void open(const char* filename) {
		    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	    } 


            /// Store data to the filestream.

            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T The type of data to be written.
            /// \param[in] t Location of the data to be written.
            /// \param[in] n The number of data items to be written.
            template <class T>
            typename std::enable_if< madness::is_ostreammable_v<T> >::type
            store(const T* t, long n) const {
		if (n > 0) {
			var.store(t, n);
			//here flush, maybe
		}
            }

            void close(){
		   // flush();
		    status = H5Dclose(dataset_id);
		    status = H5Sclose(dataspace_id);
		    status = H5Fclose(file_id);
	    }

            /// Flush the filestream.
            void flush()  {
		   int dims = 1 ;
		   hsize_t dimens_1d = v.size();

		   dataspace_id = H5Screate_simple(dims, &dimens_1d, NULL);

//		    /* Create the dataset. */
                    dataset_id =
			    H5Dcreate2(file_id, "/dset", H5T_NATIVE_UCHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		    status = H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
            }

            /// Destructor.
//            ~H5OutputArchive() {
//                close();
//            }
        }; // class H5OutputArchive


        /// Wraps an archive around a text filestream for input.
        class H5InputArchive : public BaseInputArchive {
	    hid_t   file_id; /* identifiers */
	    hid_t   dataset_id; /* identifiers */
	    hid_t   dataspace_id; /* identifiers */
	    herr_t  status;
	    mutable std::vector<unsigned char> v;
	    madness::archive::VectorInputArchive var;
	    mutable std::size_t i; ///< Current input location.

        public:
            /// Default constructor.

            /// The filename and open modes are optional here; they can be
            /// specified later by calling \c open().
            /// \param[in] filename Name of the file to read from.
            /// \param[in] mode I/O attributes for opening the file.
            H5InputArchive(const char* filename = nullptr
                    ): v(), var(v)
            {
                if (filename)
                    open(filename);
            }
            void open(const char* filename) {

		    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
//
		    dataset_id = H5Dopen2(file_id, "/dset", H5P_DEFAULT);
//
		    hid_t datatype = H5Dget_type(dataset_id);
		    size_t size = H5Tget_size(datatype);
		    H5T_class_t t_class  = H5Tget_class(datatype);

		    hid_t dataspace_id = H5Dget_space(dataset_id); /* dataspace handle */
		    int rank = H5Sget_simple_extent_ndims(dataspace_id);
		    hsize_t dims_out[1];
		    status  = H5Sget_simple_extent_dims(dataspace_id, dims_out, NULL);
		    v.resize((unsigned long)(dims_out[0]));
		    status = H5Dread(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
	    } 
            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T The type of data to be read.
            /// \param[out] t Where to put the loaded data.
            /// \param[in] n The number of data items to be loaded.
            template <class T>
            typename std::enable_if< madness::is_istreammable_v<T> >::type
            load(T* t, long n) const {
		var.load(t,n);
            }
            /// Close the filestream.
            void close() {
		    status = H5Dclose(dataset_id);
//		    status = H5Sclose(dataspace_id);
		    status = H5Fclose(file_id);
            }
            /// Destructor.
//            ~H5InputArchive() {
//                close();
//            }
        }; // class H5InputArchive

        /// Implement pre/postamble storage routines for a \c TextFstreamOutputArchive.

        /// \tparam T The type to be stored.
        template <class T>
        struct ArchivePrePostImpl<H5OutputArchive, T> {
            /// Write the preamble to the archive.

            /// \param[in] ar The archive.
            static void preamble_store(const H5OutputArchive& ar) {
            }

            /// Write the postamble to the archive.

            /// \param[in] ar The archive.
            static inline void postamble_store(const H5OutputArchive& ar) {
            }
        }; // struct ArchivePrePostImpl<TextFstreamOutputArchive,T>

        /// Implement pre/postamble load routines for a \c TextFstreamInputArchive.

        /// \tparam T The expected type to be loaded.
        template <class T>
        struct ArchivePrePostImpl<H5InputArchive,T> {
            /// Load the preamble and perform type-checking.

            /// \param[in] ar The archive.
            static inline void preamble_load(const H5InputArchive& ar) {
            }

            /// Load the postamble and perform type-checking.

            /// \param[in] ar The archive.
            static inline void postamble_load(const H5InputArchive& ar) {
            }
        }; // struct ArchivePrePostImpl<TextFstreamInputArchive,T>

        /// @}
    }
}

#endif // HAVE_HDF5

#endif // MADNESS_WORLD_H5_ARCHIVE_H__INCLUDED
