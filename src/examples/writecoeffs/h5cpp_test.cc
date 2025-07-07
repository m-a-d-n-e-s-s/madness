#include <h5cpp/h5cpp.hpp>
#include <iostream>
#include <vector>

using namespace hdf5;

//
// h5cpp has native support for std::vector IO. This example shows how to do it.
//

using DataType = std::vector<int>;

int main() {
  file::File f = file::create("std_vector_io.h5", file::AccessFlags::Truncate);
  node::Group root_group = f.root();

  //
  // writing data to a dataset in the file
  //
  DataType write_data{1, 2, 3, 4};
  node::Dataset dataset(root_group, "data", datatype::create<DataType>(),
                        dataspace::create(write_data));
  dataset.write(write_data);

  //
  // retrieving the data back from disk
  //
  DataType read_data(static_cast<size_t>(
      dataset.dataspace().size())); // allocate enough memory
  dataset.read(read_data);

  //
  // print the data
  //
  std::for_each(read_data.begin(), read_data.end(),
                [](int value) { std::cout << value << "\t"; });
  std::cout << std::endl;

  return 0;
}
