#define MAD_ARCHIVE_TYPE_NAMES_CC
#define ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE
#include <serialize/archive.h>

namespace madness {
  namespace archive {

    template <typename T> const unsigned char archive_typeinfo<T>::cookie;
    

    // Forces initialization of type names at startup
    // (breaks on shared libs ?)
    static BaseArchive fred_and_mary_sitting_under_a_tree;

    void archive_initialize_type_names() {
      static  bool initialized = false;
      if (initialized) return;
      initialized = true;

      
      for (int i=0; i<255; i++) archive_type_names[i] = "invalid";
      archive_type_names[255] = "unknown/user-defined";
      
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned char);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned short);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned int);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned long);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(unsigned long long);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed char);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed short);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed int);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed long);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(signed long long);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(bool);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(float);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(double);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(long double);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::complex<float>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::complex<double>);
      
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<char>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<unsigned char>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<short>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<unsigned short>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<int>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<unsigned int>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<long>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<unsigned long>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<bool>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<float>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::vector<double>);
      
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(std::string);
      
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor<int>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor<long>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor<float>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor<double>);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor< std::complex<float> >);
      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(Tensor< std::complex<double> >);

      ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(OctTree<double>);
    }
  }
}
