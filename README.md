# jacobi-eigenvalue-algorithm-cpp
[Jacobi eigenvalue algorithm](https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm) realization in C++

# Example
Compile sources:
```
g++ -std=c++11 main.cpp
```
Launch:
```
./a.out
```
Enter your symmetric matrix dimension and matrix itself (Note that jacobi eigenvalue algorithm works only for symmetric matrices) and get solution:
```
5
-51  21 -84 -60 -16 
 21 -81  75  86 -38 
-84  75 -31  24 -49 
-60  86  24  30 -38 
-16 -38 -49 -38 -67 
Spent time: 0 milliseconds
Eigenvalues: 
λ1: -29.1238
λ2: -201.963
λ3: -9.693
λ4: 149.671
λ5: -108.892
Eigenvectors:
X1: 0.398252 0.296203 0.451736 -0.537932 -0.51012 
X2: -0.513184 0.582347 -0.511987 -0.324183 -0.174034 
X3: -0.647166 -0.294936 0.508216 -0.42961 0.226584 
X4: -0.330274 0.411707 0.483278 0.64414 -0.270079 
X5: 0.223901 0.562712 0.204848 -0.0778015 0.764988 
Number of iterations: 34
```
# Some notes
* This program will output eigenvectors to output.txt file if your matrix dimension is more than 10
* On every iteration of algorithm we need to choose `i` and `j` indices to zero corresponding element in matrix. There are 3 strategies implemented: 
  1. off-diagonal element with maximum absolute value
  2. off-diagonal element with maximum absolute value on the row having maximum sum of squares of off-diagonal elements
  3. off-diagonal cycle element choice. Just choose off-diagonal elements by cycle: 
     ![note](https://latex.codecogs.com/gif.latex?a_%7B12%7D%2C%20a_%7B13%7D%2C%20a_%7B14%7D%2C%20...%2C%20a_%7B21%7D%2C%20a_%7B23%7D%2C%20a_%7B24%7D%2C%20...)  
  you can specify strategy via command line arguments: `[strategy]`
* You can change f function in sources:
  ```cpp
  calc_type f(int i, int j, int n) {
      return i == j;
  }
  ```
  And use `[use_f]` command line argument in order to use you `f` function to construct matrix elements.
  
# Command line arguments
```
Usage: ./a.out [strategy] [input_file] [use_f]
where options include:
    [strategy] can be 1 2 or 3. If not specified 1 is used by default
    [use_f] can be 1 or 0. If not specified 0 is used by default
```
