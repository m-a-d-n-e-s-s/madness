#include <h5cpp/hdf5.hpp>
#include <string>

using namespace hdf5;

void createFile() {
  file::File File = file::create("Tensor.hdf5", file::AccessFlags::Truncate);
  node::Group RootGroup = File.root();

  std::vector<int> Data{1, 2, 3, 4, 5, 6};
  Dimensions Shape{2, 3};
  Dimensions MaxShape{dataspace::Simple::unlimited, 3};
  Dimensions ChunkSize{512, 3};
  dataspace::Simple Dataspace{Shape, MaxShape};
  datatype::Datatype Datatype = datatype::create<std::int32_t>();
  auto Dataset = node::ChunkedDataset(RootGroup, "test_data", Datatype,
                                      Dataspace, ChunkSize);
  Dataset.write(Data);
}

hdf5::file::File create_file(const std::string &filename);
