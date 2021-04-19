#include "TDA2.h"
#include "TDA_Basic_Operators.h"

double my_zero(const coord_3d & r)
{
   return 0.0;
}
// Prints norms of the given vector
void print_norms(World & world,
                 std::vector<std::vector<real_function_3d>> f)
{
   // Container
   Tensor<double> norms(f.size(),f[0].size());

   // Calc the norms
   for(int i = 0; i < f.size(); i++)
   {
      for(int j = 0; j < f[0].size(); j++)
      {
         norms(i,j) = f[i][j].norm2();
      }
   }

   // Print em in a smart way
   if(world.rank() == 0) print(norms);

}

int main(int argc, char** argv)
{
   initialize(argc, argv);
   World world(SafeMPI::COMM_WORLD);
   startup(world, argc, argv);

   // Set box size to be reasonable
   FunctionDefaults<3>::set_cubic_cell(-30,30);

   // Create a zero function
   real_function_3d zero = real_factory_3d(world).f(my_zero);

   // Create gaussians far apart from each other
   real_function_3d gaussian1 = real_factory_3d(world).functor(GaussianGuess<3>(Vector<double,3>{0.0, 0.0, 10.0}, 0.5, std::vector<int>{0,0,0}));
   real_function_3d gaussian2 = real_factory_3d(world).functor(GaussianGuess<3>(Vector<double,3>{0.0, 10.0, 0.0}, 0.5, std::vector<int>{0,0,0}));
   real_function_3d gaussian3 = real_factory_3d(world).functor(GaussianGuess<3>(Vector<double,3>{10.0, 0.0, 0.0}, 0.5, std::vector<int>{0,0,0}));

   // Normalize
   gaussian1.scale(1.0/sqrt(inner(gaussian1, gaussian1)));
   gaussian2.scale(1.0/sqrt(inner(gaussian2, gaussian2)));
   gaussian3.scale(1.0/sqrt(inner(gaussian3, gaussian3)));
 
   // Make a vector of functions
   std::vector<real_function_3d> vec1;
   std::vector<real_function_3d> vec2;
   std::vector<real_function_3d> vec3;

   // Add in functions to each vector
   vec1.push_back(gaussian1); vec1.push_back(zero); vec1.push_back(zero);
   vec2.push_back(zero); vec2.push_back(gaussian2); vec2.push_back(zero);
   vec3.push_back(zero); vec3.push_back(zero); vec3.push_back(gaussian3);
  
   // Construct a vector of vectors of functions
   std::vector<std::vector<real_function_3d>> mat1; 
   std::vector<std::vector<real_function_3d>> mat2; 
   std::vector<std::vector<real_function_3d>> mat3; 

   // Add in functions to each vector
   for(int i = 0; i < 3; i++) mat1.push_back(vec1); 
   for(int i = 0; i < 3; i++) mat2.push_back(vec2); 
   for(int i = 0; i < 3; i++) mat3.push_back(vec3); 

   // Create tensor for multiplication
   Tensor<double> ten1(3);
   for(int i = 0; i < 3; i++) ten1(i) = i+1;

   // Create a zero matrix for inplace addition
   std::vector<std::vector<real_function_3d>> temp;
   for(int i = 0; i < 3; i++) temp.push_back(std::vector<real_function_3d>{zero,zero,zero});
 
   // Establish base case 
   if(world.rank() == 0) print("\n\nmat1 norms:");
   print_norms(world, mat1);

   if(world.rank() == 0) print("\nmat2 norms:");
   print_norms(world, mat2);

   if(world.rank() == 0) print("\nmat3 norms:");
   print_norms(world, mat3);

   // Test addition
   if(world.rank() == 0) print("\n\nmat1 + mat2 norms:");
   print_norms(world, mat1 + mat2);

   if(world.rank() == 0) print("\nmat2 + mat3 norms:");
   print_norms(world, mat2 + mat3);

   if(world.rank() == 0) print("\nmat1 + mat3 norms:");
   print_norms(world, mat1 + mat3);

   if(world.rank() == 0) print("\nmat1 + mat2 + mat3 norms:");
   print_norms(world, mat1 + mat2 + mat3);

   // Test inplace addition
   if(world.rank() == 0) print("\n\ntemp norms:");
   print_norms(world, temp);

   temp += mat1;
   if(world.rank() == 0) print("\ntemp += mat1 norms:");
   print_norms(world, temp);

   temp += mat2;
   if(world.rank() == 0) print("\ntemp += mat2 norms:");
   print_norms(world, temp);

   temp += mat3;
   if(world.rank() == 0) print("\ntemp += mat3 norms:");
   print_norms(world, temp);

   // Test addition of a constant   
   if(world.rank() == 0) print("\n\nmat1[0][0](0,0,20):", mat1[0][0](coord_3d{0.0,0.0,20.0}));
   mat1 = mat1 + 3;   
   if(world.rank() == 0) print("\nmat1[0][0](0,0,20) + 3:", mat1[0][0](coord_3d{0.0,0.0,20.0}));


   // Test subtraction of a constant
   if(world.rank() == 0) print("\n\nmat1[0][0](0,0,20) + 3:", mat1[0][0](coord_3d{0.0,0.0,20.0}));
   mat1 = mat1 - 3;   
   if(world.rank() == 0) print("\nmat1[0][0](0,0,20) - 3:", mat1[0][0](coord_3d{0.0,0.0,20.0}));

   // Test scaling by constant
   if(world.rank() == 0) print("\n\n3.0 * mat1 norms:");
   print_norms(world, scale(mat1, 3.0));
  
   // Test scaling by a tensor 
   if(world.rank() == 0) print("\n\nten1 * mat1 norms:");
   print_norms(world, scale(mat1, ten1)); 

   // Test inner

   finalize();
   return 0;
}
