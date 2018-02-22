/*
 * Some basic operators for std::vector<std::vector<real_function_3d>> objects
 */

#ifndef MADNESS_APPS_TDHF_OPS_H_INCLUDED
#define MADNESS_APPS_TDHF_OPS_H_INCLUDED

using namespace madness;

// Addition of two vectors of vectors g[i][j] = a[i][j] + b[i][j]
std::vector<std::vector<real_function_3d>> operator+(std::vector<std::vector<real_function_3d>> a,
                                                     std::vector<std::vector<real_function_3d>> b)
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a.size() == b.size()); 
 
   std::vector<std::vector<real_function_3d>> result;

   for(unsigned int i = 0; i < a.size(); i++)
   {
      // Using vmra.h definitions
      result.push_back(a[i] + b[i]); 
   }

   return result;
}


// Multiplication of a vector of vectors by a function g[i][j] = a[i][j] * b 
std::vector<std::vector<real_function_3d>> multiply(std::vector<std::vector<real_function_3d>> a,
                                                    real_function_3d b)
{
   MADNESS_ASSERT(a.size() > 0);   
   MADNESS_ASSERT(a[0].size() > 0); 

   std::vector<std::vector<real_function_3d>> result;

   for(unsigned int i = 0; i < a.size(); i++)
   {
      // Using vmra.h definitions
      // TESTING!!!!!!!!!!! (the false)
      result.push_back(mul(a[i][0].world(), b, a[i], false));
   }

   // TESTING!!!!!!!!!!
   // Not needed ??????	
   a[0][0].world().gop.fence();

   return result;
}

// Multiplication of a vector of vectors by a scalar g[i][j] = a[i][j] * b(j)
std::vector<std::vector<real_function_3d>> scale(std::vector<std::vector<real_function_3d>> a,
                                                 Tensor<double> b)
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a[0].size() > 0);
   MADNESS_ASSERT(a[0].size() == (unsigned int)b.size());    

   std::vector<std::vector<real_function_3d>> result;

   for(unsigned int i = 0; i < a.size(); i++)
   {   
      std::vector<real_function_3d> temp;

      for(unsigned int j = 0; j < a[0].size(); j++)
      { 
         temp.push_back(a[i][j] * b[j]);
      }
     
      result.push_back(temp);
   }

   return result;
}

// Returns a shallow copy of the transpose of a vector of vector of functions
std::vector<std::vector<real_function_3d>> transpose(std::vector<std::vector<real_function_3d>> f)
{
   MADNESS_ASSERT(f.size() > 0);
   MADNESS_ASSERT(f[0].size() > 0);

   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // Return container
   std::vector<std::vector<real_function_3d>> g(n, std::vector<real_function_3d>(m));

   // Now do shallow copies
   for(int i = 0; i < m; i++)
      for(int j = 0; j < n; j++)
         g[j][i] = f[i][j];
      
   // Done
   return g;
}

// Multiplication of a vector of vectors by a matrix, g[i][k] = \sum_{j} a[i][j] * b(j,k)
// NOTE: NO BOUNDS CHECKING ON THE TENSOR b!!!!
std::vector<std::vector<real_function_3d>> scale_2d(World & world,
                                                    std::vector<std::vector<real_function_3d>>& a,
                                                    Tensor<double>& b)
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a[0].size() > 0);

   std::vector<std::vector<real_function_3d>> result;

   for(unsigned int i = 0; i < a.size(); i++)
   {
      // Using vmra.h definitions
      result.push_back(transform(world, a[i], b, false));
   }

   // Does this need to be here?
   world.gop.fence(); 

   return result;
}

// Multiplication of a vector of vectors by a scalar g[i][j] = a[i][j] * b
std::vector<std::vector<real_function_3d>> scale(std::vector<std::vector<real_function_3d>> a,
                                                 double b)
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a[0].size() > 0);   
 
   std::vector<std::vector<real_function_3d>> result;

   for(unsigned int i = 0; i < a.size(); i++)
   {
      // Using vmra.h definitions
      result.push_back(a[i] * b);
   }

   return result;
}

// Truncate a vector of vector of functions
void truncate(World & world,
              std::vector<std::vector<real_function_3d>> & v,
              double tol = 0.0,
              bool fence = true)
{
   MADNESS_ASSERT(v.size() > 0);
   MADNESS_ASSERT(v[0].size() > 0);

   for(unsigned int i = 0; i < v.size(); i++)
   {
      truncate(world, v[i], tol, fence);
   }
}

// Returns a deep copy of a vector of vector of functions
std::vector<std::vector<real_function_3d>> copy(World & world,
                                                std::vector<std::vector<real_function_3d>> f)
{
   MADNESS_ASSERT(f.size() > 0);

   int m = f.size();
   std::vector<std::vector<real_function_3d>> answer;

   for(int i = 0; i < m; i++) answer.push_back(copy(world,f[i]));

   return answer;
}

// Apply a vector of vector of operators to a vector of vector of functions g[i][j] = op[i][j](f[i][j])
std::vector<std::vector<real_function_3d>> apply(World & world,
                                                 std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> & op,
                                                 std::vector<std::vector<real_function_3d>> & f)
{
   MADNESS_ASSERT(f.size() > 0);
   MADNESS_ASSERT(f.size() == op.size());
   MADNESS_ASSERT(f[0].size() == op[0].size());

   std::vector<std::vector<real_function_3d>> result(f.size());
   
   for(unsigned int i = 0; i < f.size(); i++)
   {
      // Using vmra.h function, line 889
      result[i] = apply(world, op[i], f[i]);
   }

   return result;
}

// Apply the derivative operator to a vector of vector of functions
std::vector<std::vector<real_function_3d>> apply(World & world,
                                                 real_derivative_3d & op,
                                                 std::vector<std::vector<real_function_3d>> f)
{
   MADNESS_ASSERT(f.size() > 0);

   std::vector<std::vector<real_function_3d>> result;

   for(unsigned int i = 0; i < f.size(); i++)
   {
      result.push_back(apply(world, op, f[i]));
   }

   return result;
}

/*
 *
 *  These functions are here so that the KAIN solver is happy
 *
 */

// Subtraction of two vectors of vectors g[i][j] = a[i][j] - b[i][j]
std::vector<std::vector<real_function_3d>> operator-(const std::vector<std::vector<real_function_3d>> & a,
                                                     const std::vector<std::vector<real_function_3d>> & b) 
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a.size() == b.size());
 
   std::vector<std::vector<real_function_3d>> result;

   for(unsigned int i = 0; i < a.size(); i++)
   {
      // Using vmra.h definitions
      result.push_back(a[i] - b[i]);
   }

   return result;
}

// Multiplication of a vector of vectors by a scalar g[i][j] = a[i][j] * b
std::vector<std::vector<real_function_3d>> operator*(std::vector<std::vector<real_function_3d>> & a,
                                                     double b)
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a[0].size() > 0);   
 
   std::vector<std::vector<real_function_3d>> result;

   for(unsigned int i = 0; i < a.size(); i++)
   {
      // Using vmra.h definitions
      result.push_back(a[i] * b);
   }

   return result;
}

// Addition in place of a vector of vector of functions
void operator+=(std::vector<std::vector<real_function_3d>> & a,
                const std::vector<std::vector<real_function_3d>> & b)
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a.size() == b.size()); 
 
   for(unsigned int i = 0; i < a.size(); i++)
   {
      // Using vmra.h definitions
      a[i] += b[i]; 
   }
}

// Inner product for std::vector<std::vector<real_function_3d>> and std::vector<std::vector<real_function_3d>>
double inner(std::vector<std::vector<real_function_3d>> a,
             std::vector<std::vector<real_function_3d>> b)
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a.size() == b.size());
   MADNESS_ASSERT(a[0].size() > 0);
   MADNESS_ASSERT(a[0].size() == b[0].size());

   double value = 0.0;

   for(unsigned int i = 0; i < a.size(); i++)
   {
      // vmra.h function, line 580 
      value += inner(a[i][0].world(), a[i], b[i]).sum();
   }

   return value;
}

#endif

// Deuces
