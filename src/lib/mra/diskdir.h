#ifndef DISKDIR_H
#define DISKDIR_H
#include <vector>
using namespace std;
#include <misc/communicator.h>
#include <mra/mra.h>
using madness::Function;

#include <serialize/textfsar.h>
using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;

#include <serialize/binfsar.h>
using madness::archive::BinaryFstreamInputArchive;
using madness::archive::BinaryFstreamOutputArchive;

namespace madness {

    /// This class can serialize many datas of Function class into same file.
    /// rw_status     = "read"   : read mode
    ///                 "write"  : write mode
    /// delete_status = "keep"   : keep file mode(not finished yet)
    ///                 "delete" : delete file mode(not finished yet)
    template <typename T>
    class DiskDir {
        private:
	    char* rw_status, delete_status;
            Communicator DDcomm;
            TextFstreamInputArchive DDiar;
            TextFstreamOutputArchive DDoar;
            //BinaryFstreamInputArchive DDiar;
            //BinaryFstreamOutputArchive DDoar;

        public:
            DiskDir() {};

            DiskDir(const char* filename, const char* delete_status, 
                    const char* rw_status, Communicator& DDcomm) {
              if(DDcomm.rank() == 0) {
                if(rw_status == "read") {
                  DDiar.open(filename, ios::in);
                  DDoar.close();
                }
                else if (rw_status == "write") {
                  DDoar.open(filename, ios::out | ios::trunc);
                  DDiar.close();
                }
              }
              else {
                  DDiar.close();
                  DDoar.close();
              }
            };

            ~DiskDir() {
              classClose();
            };

            void operator<<(Function<T>& f) {
              if (DDcomm.rank() == 0) {
                f.saveMain(DDoar, DDcomm);
              }
              else {
                f.saveLoadWorker4DD(DDcomm, true);
              }
            };

            void operator>>(Function<T>& f) {
              if (DDcomm.rank() == 0) {
                DDiar & FunctionDefaults::k;
                DDiar & FunctionDefaults::thresh;
                DDiar & FunctionDefaults::initial_level;
                DDiar & FunctionDefaults::max_refine_level;
                DDiar & FunctionDefaults::truncate_method;
                DDiar & FunctionDefaults::autorefine;
              }
              DDcomm.Bcast(FunctionDefaults::k, 0);
              DDcomm.Bcast(FunctionDefaults::thresh, 0);
              DDcomm.Bcast(FunctionDefaults::initial_level, 0);
              DDcomm.Bcast(FunctionDefaults::max_refine_level, 0);
              DDcomm.Bcast(FunctionDefaults::truncate_method, 0);
              DDcomm.Bcast(FunctionDefaults::autorefine, 0);
              if (DDcomm.rank() == 0) {
                bool active_flag;
                DDiar & active_flag;
                if(active_flag) f.loadManager4DD(DDiar, DDcomm, active_flag);
                for ( int i = 1; i < DDcomm.size(); i++) {
                  archive::MPIOutputArchive arout(DDcomm, i);
                  arout & -1;
                }
              }
              else {
                f.saveLoadWorker4DD(DDcomm, false);
              }
            };
/*
            DiskDir<T> operator<<(const vector< Function<T>& > f) {
            };
            DiskDir<T> operator>>(const vector< Function<T> > f) {
            };
*/
            void classClose() {
               if(DDcomm.rank() == 0) {
                 if(rw_status == "read") {
                   DDiar.close();
	         }
	         else if (rw_status == "write") {
                   DDoar.close();
	         }
               }
            };
    };
};

#endif
