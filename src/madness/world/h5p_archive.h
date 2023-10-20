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

#ifdef H5_HAVE_PARALLEL

#ifndef MADNESS_WORLD_H5P_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_H5P_ARCHIVE_H__INCLUDED

/**
 \file text_fstream_archive.h
 \brief Implements an archive wrapping text filestream.
 \ingroup serialization
*/

#include <type_traits>
#include <fstream>
#include <cstring>
#include <madness/world/archive.h>
#include <madness/world/print.h>  // this injects operator<<(std::ostream,T) for common Ts
#include <madness/world/vector_archive.h>

#include "hdf5.h"

namespace madness {
    namespace archive {

        /// \addtogroup serialization
        /// @{

        /// Wraps an archive around a text filestream for output.
        class H5POutputArchive : public BaseOutputArchive {
            mutable std::ofstream os; ///< The filestream.
	    hid_t   file_id;                                         /* HDF5 file IDs  */
            hid_t    acc_tpl1;                                        /* File access templates */
	    hid_t   dataset_id;                                      /* Dataset transfer properties list */
	    hid_t   dataspace_id;                                    /* Dataspace ID */
	    hid_t   xfer_plist;                                      /* Dataset transfer properties list */
	    hid_t   file_dataspace;                                  /* File dataspace ID */
	    hid_t   mem_dataspace;                                   /* memory dataspace ID */
	    hid_t   dataset1, dataset2;                              /* Dataset ID */

	    herr_t  status;
	    mutable std::vector<unsigned char> v;
            madness::archive::VectorOutputArchive var;
            mutable World* world; ///< The world.

        public:
            /// Default constructor.

            /// The filename and open modes are optional here; they can be
            /// specified later by calling \c open().
            /// \param[in] filename Name of the file to write to.
            /// \param[in] mode I/O attributes for opening the file.
            H5POutputArchive(World& world, const char* filename = nullptr
                    ):world(&world), v(), var(v)
            {
		    if (filename)
			    open( world, filename);
            };
            void open(World& world,const char* filename) {
		    std::cout << "h5p this open " << filename << std::endl;
		    auto iocomm = world.mpi.comm().Get_mpi_comm();
		    std::cout << "h5p this open comm " << iocomm << std::endl;
                    auto ioinfo = world.mpi.comm().Get_info();

                     herr_t ret;
//                    MPI_Info info = MPI_INFO_NULL;
//		    /* -------------------
//		     * START AN HDF5 FILE
//		     * -------------------*/
//		    /* setup file access template with parallel IO access. */
		    acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
		    std::cout << "H5Pcreate access succeed" << std::endl;
//
//		    /* set Parallel access with communicator */
		    ret = H5Pset_fapl_mpio(acc_tpl1, iocomm, ioinfo);
		    std::cout << "H5Pset_fapl_mpio succeed"  << std::endl;
//
//		    /* create the file collectively */
		    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);
		    std::cout << "H5Fcreate succeed" << std::endl;
//
		    /* Release file-access template */
		    ret = H5Pclose(acc_tpl1);

//		    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

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
                //using madness::operators::operator<<;
	        std::cout << "h5p this store 1 v " << v.size() << std::endl;
	        std::cout << "h5p this store 1 n " << n << std::endl;
		if (n > 0) {
			var.store(t, n);
//			//here flush
		}
            }

            void close(){
		   // flush();
		    std::cout << "h5p this close " <<  std::endl;
//		    status = H5Dclose(dataset_id);
//		    status = H5Sclose(dataspace_id);
		    status = H5Fclose(file_id);
	    }

            /// Flush the filestream.
            void flush()  {
		   std::cout << "h5p this flush " <<  std::endl;
		   int dims = 1 ;
		   hsize_t dimens_1d = v.size();
	           std::cout << "h5p this flush 1 v " << v.size() << std::endl;

                     herr_t ret;
                   std::cout <<"h5p flush my rank " <<world->rank() << std::endl;
                   int my =2 ;
                   int total =2 ;
                     
                   
//
////		   dataspace_id = H5Screate_simple(dims, &dimens_1d, NULL);
//                  dataspace_id = H5Screate_simple(dims, &dimens_1d, NULL);
//		   std::cout << "h5p this dataspace created " <<  std::endl;
//
////		    /* Create the dataset. */
//
//                    dataset_id =
//                        H5Dcreate2(file_id, "/dset", H5T_NATIVE_UCHAR,  dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
////			    H5Dcreate2(file_id, "/dset", H5T_NATIVE_UCHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
////		    status = H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
//                  file_dataspace = H5Dget_space (dataset_id);
//                  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
//                 ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
//		   std::cout << "h5p this flush H5Pcreate xfer succeed created " <<  std::endl;
//
//
//
//                 hsize_t stride = 1;
//                 hsize_t count  = dimens_1d / world.mpi.comm().Get_rank();
//                 hsize_t start  = mpi_rank * count[0];

//                 ret = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, &file_offset, &stride, &dsize, NULL); //stride=NULL?
//                 ret = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
//
//                 ret = H5Dwrite(dataset1, H5T_NATIVE_INT, mem_dataspace, file_dataspace, xfer_plist, data_array1);
            }

            /// Destructor.
//            ~H5POutputArchive() {
//                close();
//            }
        }; // class H5POutputArchive


        /// Wraps an archive around a text filestream for input.
        class H5PInputArchive : public BaseInputArchive {
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
            H5PInputArchive(const char* filename = nullptr
                    ): v(), var(v)
            {
                if (filename)
                    open(filename);
            }
            void open(const char* filename) {
		    std::cout << "h5p input this open " << filename << std::endl;
		    //file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
//
		    dataset_id = H5Dopen2(file_id, "/dset", H5P_DEFAULT);
//
		    hid_t datatype = H5Dget_type(dataset_id);
		    size_t size = H5Tget_size(datatype);
		    H5T_class_t t_class  = H5Tget_class(datatype);
		    if(t_class == H5T_INTEGER)
		    std::cout << "h5p input ajua"  << (int)t_class<< std::endl;
		    std::cout << "h5p input  size "<< (int)size << std::endl;
		    hid_t dataspace_id = H5Dget_space(dataset_id); /* dataspace handle */
		    int rank = H5Sget_simple_extent_ndims(dataspace_id);
		    std::cout << "h5p input  rank "<< (int)rank << std::endl;
		    hsize_t dims_out[1];
		    status  = H5Sget_simple_extent_dims(dataspace_id, dims_out, NULL);
		    printf("h5p input rank %d, dimensions %lu x %lu \n", rank, (unsigned long)(dims_out[0]),
           (unsigned long)(dims_out[0]));
		v.resize((unsigned long)(dims_out[0]));
		    status = H5Dread(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
//  	    status = H5Dread(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
//
//		    int dims = 1 ;
//		    hsize_t dimens_1d = v->size();
//
//		    dataspace_id = H5Screate_simple(dims, &dimens_1d, NULL);
//
//		    /* Create the dataset. */
//                    dataset_id =
//		      H5Dcreate2(file_id, "/dset", H5T_NATIVE_UCHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    } 
            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T The type of data to be read.
            /// \param[out] t Where to put the loaded data.
            /// \param[in] n The number of data items to be loaded.
            template <class T>
            typename std::enable_if< madness::is_istreammable_v<T> >::type
            load(T* t, long n) const {
		    std::cout << "h5p input this load " << std::endl;
	        std::cout << "h5p input this load 1 n " << n << std::endl;
                std::cout << "h5p this store 1 v " << v.size() << std::endl;
		var.load(t,n);
            }
            /// Close the filestream.
            void close() {
		    std::cout << "h5p inp this close " <<  std::endl;
		    status = H5Dclose(dataset_id);
//		    status = H5Sclose(dataspace_id);
		    status = H5Fclose(file_id);
            }
            /// Destructor.
//            ~H5PInputArchive() {
//                close();
//            }
        }; // class H5PInputArchive

        /// Implement pre/postamble storage routines for a \c TextFstreamOutputArchive.

        /// \tparam T The type to be stored.
        template <class T>
        struct ArchivePrePostImpl<H5POutputArchive, T> {
            /// Write the preamble to the archive.

            /// \param[in] ar The archive.
            static void preamble_store(const H5POutputArchive& ar) {
            }

            /// Write the postamble to the archive.

            /// \param[in] ar The archive.
            static inline void postamble_store(const H5POutputArchive& ar) {
            }
        }; // struct ArchivePrePostImpl<TextFstreamOutputArchive,T>

        /// Implement pre/postamble load routines for a \c TextFstreamInputArchive.

        /// \tparam T The expected type to be loaded.
        template <class T>
        struct ArchivePrePostImpl<H5PInputArchive,T> {
            /// Load the preamble and perform type-checking.

            /// \param[in] ar The archive.
            static inline void preamble_load(const H5PInputArchive& ar) {
            }

            /// Load the postamble and perform type-checking.

            /// \param[in] ar The archive.
            static inline void postamble_load(const H5PInputArchive& ar) {
            }
        }; // struct ArchivePrePostImpl<TextFstreamInputArchive,T>

        /// @}
    }

    template <>
    struct is_archive<archive::H5POutputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::H5PInputArchive> : std::true_type {};
    template <>
    struct is_output_archive<archive::H5POutputArchive> : std::true_type {};
    template <>
    struct is_input_archive<archive::H5PInputArchive> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::H5POutputArchive, T, std::enable_if_t<is_iostreammable_v<T> || std::is_function_v<T> || is_any_function_pointer_v<T>>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::H5PInputArchive, T, std::enable_if_t<is_iostreammable_v<T> || is_any_function_pointer_v<T>>> : std::true_type {};

} //namespace

#endif // MADNESS_WORLD_H5P_ARCHIVE_H__INCLUDED
#endif // H5_HAVE_PARALLEL
