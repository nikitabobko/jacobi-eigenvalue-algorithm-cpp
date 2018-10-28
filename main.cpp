#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>

using namespace std;

typedef double calc_type;

template<typename T>
class Matrix {
public:
    T *matrix_pointer = nullptr;
    int n = 0;
    int m = 0;

    Matrix() = default;

    Matrix(int n, int m) : matrix_pointer(new T[n * m]), n(n), m(m) {
        fill_with_zeros();
    }

    void fill_with_zeros() {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                at(i, j) = T();
            }
        }
    }

    bool is_diagonal(double eps) const {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                if (i != j && std::fabs(at(i, j)) >= eps) {
                    return false;
                }
            }
        }
        return true;
    }

    std::vector<calc_type> get_line(int line_index) const {
        std::vector<calc_type> ret(m, calc_type());
        for (int j = 0; j < m; ++j) {
            ret[j] = at(line_index, j);
        }
        return ret;
    }

    std::vector<calc_type> get_column(int column_index) const {
        std::vector<calc_type> ret(n, calc_type());
        for (int i = 0; i < n; ++i) {
            ret[i] = at(i, column_index);
        }
        return ret;
    }

    T &at(int i, int j) {
        return matrix_pointer[i * m + j];
    }

    T at(int i, int j) const {
        return matrix_pointer[i * m + j];
    }

    std::string to_string() {
        std::string ret;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                ret += std::to_string(at(i, j)) + " ";
            }
            ret += "\n";
        }
        return ret;
    }

    virtual ~Matrix() {
        delete[] matrix_pointer;
    }
};

void print_eigenvectors(const Matrix<calc_type> &eigenvectors, int n) {
    ostream *output_temp = &cout;
    if (n > 10) {
        output_temp = new ofstream("output.txt");
    }
    ostream &output = *output_temp;

    output << "Eigenvectors:" << std::endl;
    for (int j = 0; j < n; ++j) {
        output << "X" << j + 1 << ": ";
        for (int i = 0; i < n; ++i) {
            output << eigenvectors.at(i, j) << " ";
        }
        output << std::endl;
    }
    if (n > 10) {
        delete output_temp;
    }
}

std::pair<int, int> chose_i_and_j(const Matrix<calc_type> &matrix, int n, int strategy, int previous_i, int previous_j,
                                  double eps) {
    int chosen_j = 1;
    int chosen_i = 0;
    if (strategy == 1) {
        calc_type max = -1.;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j && fabs(matrix.at(i, j)) > max) {
                    max = fabs(matrix.at(i, j));
                    chosen_i = i;
                    chosen_j = j;
                }
            }
        }
    } else if (strategy == 2) {
        calc_type max_line_squared = 0;
        for (int i = 0; i < n; ++i) {
            calc_type max_on_line = -1.;
            calc_type line_squared = 0;
            int local_chosen_j = 1;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    line_squared += matrix.at(i, j) * matrix.at(i, j);
                    if (fabs(matrix.at(i, j)) > max_on_line) {
                        max_on_line = fabs(matrix.at(i, j));
                        local_chosen_j = j;
                    }
                }
            }
            if (line_squared > max_line_squared) {
                max_line_squared = line_squared;
                chosen_i = i;
                chosen_j = local_chosen_j;
            }
        }
    } else if (strategy == 3) {
        if (previous_i == -1 || previous_j == -1) {
            previous_i = 0;
            previous_j = 0;
        }
        chosen_i = previous_i;
        chosen_j = previous_j;

        do {
            chosen_j++;
            chosen_i += (chosen_j >= n);
            chosen_j %= n;
            chosen_i %= n;
        } while(chosen_i == chosen_j || fabs(matrix.at(chosen_i, chosen_j)) < eps);
    }
    return std::pair<int, int>(chosen_i, chosen_j);
}

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


std::vector<calc_type > jacobi_eigenvalue_algorithm(Matrix<calc_type> &matrix, Matrix<calc_type> &eigenvectors, int n,
                                                    int strategy, double eps, int *counter) {
    *counter = 0;
    int prev_i = -1, prev_j = -1;
    while (!matrix.is_diagonal(eps)) {
        std::pair<int, int> temp = chose_i_and_j(matrix, n, strategy, prev_i, prev_j, eps);
        int chosen_i = temp.first;
        int chosen_j = temp.second;

        calc_type cos, sin;
        calc_type x = -2*matrix.at(chosen_i, chosen_j);
        calc_type y = matrix.at(chosen_i, chosen_i) - matrix.at(chosen_j, chosen_j);
        if (fabs(y) < eps) {
            cos = sin = 1./sqrt(2);
        } else {
            cos = sqrt(1./2 + fabs(y)/(2*sqrt(x*x + y*y)));
            sin = sgn(x*y)*fabs(x)/(2*cos*sqrt(x*x + y*y));
        }

        std::vector<calc_type> line_i = matrix.get_line(chosen_i);
        std::vector<calc_type> line_j = matrix.get_line(chosen_j);
        for (int j = 0; j < n; ++j) {
            matrix.at(chosen_i, j) = line_i[j]*cos - line_j[j]*sin;
            matrix.at(chosen_j, j) = line_i[j]*sin + line_j[j]*cos;
        }

        std::vector<calc_type> column_i = matrix.get_column(chosen_i);
        std::vector<calc_type> column_j = matrix.get_column(chosen_j);
        for (int i = 0; i < n; ++i) {
            matrix.at(i, chosen_i) = column_i[i]*cos - column_j[i]*sin;
            matrix.at(i, chosen_j) = column_i[i]*sin + column_j[i]*cos;
        }

        std::vector<calc_type> eigenvector_column_i = eigenvectors.get_column(chosen_i);
        std::vector<calc_type> eigenvector_column_j = eigenvectors.get_column(chosen_j);
        for (int i = 0; i < n; ++i) {
            eigenvectors.at(i, chosen_i) = eigenvector_column_i[i]*cos - eigenvector_column_j[i]*sin;
            eigenvectors.at(i, chosen_j) = eigenvector_column_i[i]*sin + eigenvector_column_j[i]*cos;
        }

        prev_i = chosen_i;
        prev_j = chosen_j;
        ++(*counter);
    }
    std::vector<calc_type> ret(n, calc_type(0));
    int i = 0;
    for (auto &it : ret) {
        it = matrix.at(i, i);
        i++;
    }
    return ret;
}

calc_type f(int i, int j, int n) {
    return i == j;
}

void print_eigenvalues(const vector<calc_type> &eigenvalues) {
    std::cout << "Eigenvalues: " << std::endl;
    int i = 1;
    for (auto &it : eigenvalues) {
        std::cout << "λ" << i++ << ": " << it << std::endl;
    }
}

int main(int argc, const char **argv) {
    bool error = false;
    int strategy = 1;
    if (argc - 1 >= 1) {
        strategy = atoi(argv[1]);
        if (strategy < 1 || strategy > 3) {
            error = true;
        }
    }
    if (argc - 1 >= 2) {
        freopen(argv[2], "r", stdin);
    }

    bool use_f = false;
    if (argc - 1 >= 3) {
        int temp = atoi(argv[3]);
        if (temp != 0 && temp != 1) {
            error = true;
        } else {
            use_f = bool(temp);
        }
    }

    if (argc - 1 > 3 || error) {
        printf("Usage: %s [strategy] [input_file] [use_f]\n"
               "where options include:\n"
               "    [strategy] can be 1 2 or 3. If not specified 1 is used by default\n"
               "    [use_f] can be 1 or 0. If not specified 0 is used by default", argv[0]);
        exit(1);
    }

    int n;
    std::cin >> n;
    Matrix<calc_type> matrix(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (use_f) {
                matrix.at(i, j) = f(i, j, n);
            } else {
                std::cin >> matrix.at(i, j);
            }
        }
    }
    auto a = matrix.to_string();
    Matrix<calc_type> eigenvectors(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            eigenvectors.at(i, j) = (i == j);
        }
    }

    int counter = 0;

    std::chrono::milliseconds start = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );
    std::vector<calc_type> eigenvalues = jacobi_eigenvalue_algorithm(matrix, eigenvectors, n, strategy, 1e-12, &counter);
    auto spent = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    ) - start;
    std::cout << "Spent time: " << spent.count() << " milliseconds" << std::endl;

    print_eigenvalues(eigenvalues);

    print_eigenvectors(eigenvectors, n);

    std::cout << "Number of iterations: " << counter << std::endl;
}
