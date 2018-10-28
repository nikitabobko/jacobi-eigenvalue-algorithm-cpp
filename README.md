# jacobi-eigenvalue-algorithm-cpp
[Jacobi eigenvalue algorithm](https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm) realization on C++

# Example
Compile sources:
```
g++ -std=c++11 main.cpp
```
Launch:
```
./a.out
```
Enter your symmetric matrix dimenstion and matrix itself (Note that jacobi eigenvalue algorithm works only for symmetric matricies) and get solution:
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
X1: -0.780703 -0.580654 -0.88555 1.05452 1 
X2: 2.94876 -3.34618 2.94188 1.86276 1 
X3: -2.85619 -1.30166 2.24295 -1.89603 1 
X4: 1.22288 -1.5244 -1.78939 -2.385 1 
X5: 0.292686 0.735583 0.267779 -0.101703 1 
Number of iterations: 34
```
# Some notes
* This programm will output eigenvectors to output.txt file if your matrix dimension is more than 10
* On every iteration of algorithm we need to choose `i` and `j` indicies to zero corresponding element in matrix. There are 3 strategies implemented: 
  1. off-diagonal element with maximum absolute value
  2. off-diagonal element with maximum absolute value on the row having maximum sum of squares of off-diagonal elements
  3. off-diagonal cycle element choice. Just choose off-diagoanl elements by cycle: 
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
