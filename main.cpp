#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>


std::vector<double> split_string_into_array(std::string);
template <typename T> int sgn(T);
void recount(std::vector<std::vector<double>>&, int, int);
std::vector<std::vector<double>>  main_nums(int, const std::vector<std::vector<double>>&, double);
std::vector<std::vector<double>> main_vectors(int, const std::vector<std::vector<double>>&, double);
double norm_second(const std::vector<double>&);
double aposterior_error_estimate(const std::vector <std::vector <double >>&, std::vector <double>, double);
int normalizing_vector(std::vector<double>&);
double exponent_method(std::vector<std::vector<double>>, const std::vector<double>&, double);
double scalar_method(std::vector<std::vector<double>>, const std::vector<double>&, double);
double spectr_board(double, std::vector<std::vector<double>>, double);
std::vector<std::vector<double>> inverse(std::vector<std::vector<double>>);
double wilandt_method(double, std::vector<std::vector<double>>, double);



std::vector<double> split_string_into_array(std::string str) {
    std::vector<double> arr;
    while(str.length() != 0){
        int ind_of_space = str.find(' ');
        arr.push_back(stod(str.substr(0, ind_of_space)));
        str = str.substr(ind_of_space + 1, str.length() - ind_of_space);
        if(str.find(' ') == -1){
            if(str[0] != ' ')
                arr.push_back(stod(str));
            break;
        }
        else if(str.length() == 2 && str[0] == '-'){
            arr.push_back(stod(str));
            break;
        }
    }
    return arr;
}

template<typename T>
int sgn(T val) {
    return(T(0) < val) - (T(0) > val);
}

void recount(std::vector<std::vector<double>>& a_matrix, int ind_i, int ind_j) {
    int size = a_matrix[0].size();
    double cos_value, sin_value;
    std::vector<std::vector<double>> result_matrix = a_matrix;
    cos_value = sqrt(0.5 * (1 + fabs(result_matrix[ind_i][ind_i] - result_matrix[ind_j][ind_j]) /
                                sqrt(pow(result_matrix[ind_i][ind_i] - result_matrix[ind_j][ind_j], 2) + 4 * pow(result_matrix[ind_i][ind_j], 2))));
    sin_value = sgn(result_matrix[ind_i][ind_j] * (result_matrix[ind_i][ind_i] - result_matrix[ind_j][ind_j])) *
                sqrt(0.5 * (1 - fabs(result_matrix[ind_i][ind_i] - result_matrix[ind_j][ind_j]) /
                                sqrt(pow(result_matrix[ind_i][ind_i] - result_matrix[ind_j][ind_j], 2) + 4 * pow(result_matrix[ind_i][ind_j], 2))));
    for (int i = 0; i < size; ++i) {
        for (int j = i; j < size; ++j) {
            if (i != ind_i && i != ind_j && j == ind_i) {
                result_matrix[i][j] = cos_value * a_matrix[i][ind_i] + sin_value * a_matrix[i][ind_j];
                result_matrix[j][i] = result_matrix[i][j];
            } else if (i != ind_i && i != ind_j && j == ind_j) {
                result_matrix[i][j] = -sin_value * a_matrix[i][ind_i] + cos_value * a_matrix[i][ind_j];
                result_matrix[j][i] = result_matrix[i][j];
            } else if (i == ind_i && j == i) {
                result_matrix[i][j] = pow(cos_value, 2) * a_matrix[i][i] + 2 * cos_value * sin_value * a_matrix[ind_i][ind_j] +
                                 pow(sin_value, 2) * a_matrix[ind_j][ind_j];
            } else if (i == ind_j && j == i) {
                result_matrix[i][j] = pow(sin_value, 2) * a_matrix[i][i] - 2 * cos_value * sin_value * a_matrix[ind_i][ind_j] +
                                 pow(cos_value, 2) * a_matrix[ind_j][ind_j];
            } else if (i == ind_i && j == ind_j) {
                result_matrix[i][j] = 0;
                result_matrix[j][i] = 0;
            }
        }
    }

    a_matrix = result_matrix;
}

std::vector<std::vector<double>> main_nums(int size, const std::vector<std::vector<double>>& matrix, double eps){
    //returns main nums of matrix as matrix
    std::vector<std::vector<double>> main_nums_matrix = matrix;
    int iterations = 0;
    bool accuracy = true;
    while(accuracy) {
        iterations++;
        double max = std::abs(main_nums_matrix[0][1]);
        int ind_i = 0, ind_j = 1;
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                if (std::abs(main_nums_matrix[i][j]) > max) {
                    max = std::abs(main_nums_matrix[i][j]);
                    ind_i = i;
                    ind_j = j;
                }
            }
        }
        if (max < eps) accuracy = false;
        else recount(main_nums_matrix, ind_i, ind_j);      //evaluating main-vectors first
    }
    std::cout << "\nNumber of iterations for main numbers: " << iterations << std::endl;
    return main_nums_matrix;
}


std::vector<std::vector<double>> main_vectors(int size, const std::vector<std::vector<double>>& matrix, double eps) {
    //returns main vectors
    std::vector<std::vector<double>> x_matrix(size, std::vector<double>(size, 1));
    int iterations = 0;
    bool accuracy = true;
    while(accuracy) {
        iterations++;
        double max = std::abs(x_matrix[0][1]);
        int ind_i = 0, ind_j = 1;
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                if (fabs(x_matrix[i][j]) > max) {
                    max = std::abs(x_matrix[i][j]);
                    ind_i = i;
                    ind_j = j;
                }
            }
        }
        if (max < eps) accuracy = false;
        else recount(x_matrix, ind_i, ind_j);      //evaluating main-vectors
    }
    std::cout << "\nNumber of iterations for main vectors: " << iterations << std::endl;
    return x_matrix;
}

double norm_second(const std::vector<double>& vect){
    double sum = 0;
    for(auto el : vect)
        sum += pow(el, 2);
    return sqrt(sum);
}

double aposterior_error_estimate(const std::vector <std::vector <double >>& a_matrix, std::vector <double> main_vect, double lambda){
    int n = main_vect.size();
    std::vector<double> matr_to_vect(n), num_to_vect(n);       //for formula ||A*Y - lambda*Y||_2, where norm is second Gelder's norm, A - matrix, Y - approximation
                                                                 //to main vector , lambda - approximation to main number
    for(int i = 0; i < n; ++i){ //i - cols, j - rows
        double sum_temp = 0;
        for(int j = 0; j < n; ++j)
            sum_temp += a_matrix[i][j] * main_vect[j];
        matr_to_vect[i] = sum_temp;
        num_to_vect[i] = lambda * main_vect[i];
    }
    for(int i = 0; i < n; ++i)
        matr_to_vect[i] -= num_to_vect[i];
    return norm_second(matr_to_vect)/norm_second(main_vect);
}


int normalizing_vector(std::vector<double>& main_vector){
    int ind = 0;
    double y_p = main_vector[0];
    for(int i = 1 ; i < main_vector.size(); ++i)
        if (std::abs(main_vector[i]) > std::abs(y_p)) {
            ind = i;
            y_p = main_vector[i];
        }
    for(auto& el: main_vector)
        el /= y_p;
    return ind;
}


double exponent_method(std::vector<std::vector<double>> a_matrix, const std::vector<double>& main_vector, double eps){
    std::vector<double> new_vector = main_vector;
    int p_ind, iterations = 0;
    double lambda;
    do {
        p_ind = normalizing_vector(new_vector);
        iterations++;
        std::vector<double> temp_vector(main_vector.size(), 0);
        for (int i = 0; i < main_vector.size(); ++i) {
            for (int j = 0; j < main_vector.size(); ++j) {
                temp_vector[i] += a_matrix[i][j] * new_vector[j];
            }
        }
        new_vector = temp_vector;
        lambda = new_vector[p_ind];
    }while(aposterior_error_estimate(a_matrix, new_vector, lambda) >= eps);
    std::cout << "Number of iterations for exponent method is: " << iterations << std::endl;
    std::cout << "Main vector for founded so far main number is: \n";
    for(auto el: new_vector){
        std::cout << el << " ";
    }
    std::cout << std::endl;
    return lambda;
}


double scalar_method(std::vector<std::vector<double>> a_matrix, const std::vector<double>& main_vector, double eps){
    std::vector<double> new_vector = main_vector;
    int p_ind, iterations = 0;
    double lambda;
    int n = new_vector.size();
    do {
        p_ind = normalizing_vector(new_vector);
        iterations++;
        std::vector<double> temp_vector(n, 0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                temp_vector[i] += a_matrix[i][j] * new_vector[j];
            }
        }
        double denominator = 0, numerator = 0;
        for (int i = 0; i < n; ++i) {
            numerator += temp_vector[i]*new_vector[i];
        }
        for (int i = 0; i < n; ++i) {
            denominator += pow(new_vector[i], 2);
        }
        new_vector = temp_vector;
        lambda = numerator/denominator;
    }while(aposterior_error_estimate(a_matrix, new_vector, lambda) >= eps);
    std::cout << "Number of iterations for scalar method is: " << iterations << std::endl;
    return lambda;
}


double spectr_board(double lambda, std::vector<std::vector<double> > a_matrix, double eps){
    int n = a_matrix[0].size();     //finding opposite board of a spectre - in our case it's minimal number of matrix A
    std::vector<std::vector<double> > b_matrix = a_matrix;
    double board;
    for (int i = 0; i < n; ++i) {
            b_matrix[i][i] -= lambda;
    }
    board = exponent_method(b_matrix, a_matrix[0], eps);
    return board + lambda;
}

std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> a) {
    int n = a.size();
    std::vector < std::vector <double> > ans(n, std::vector <double> (n, 0));
    for (int i = 0; i < n; i++){
        ans[i][i] = 1.0;
    }
    for (int i = 0; i < n; i++){
        int row = i;
        double mx = a[i][i];
        for(int k = i + 1; k < n; k++){
            if (std::abs(a[k][i]) > mx){
                row = k;
                mx = std::abs(a[k][i]);
            }
        }
        if (row != i) {
            swap(a[row], a[i]);
            swap(ans[row], ans[i]);
        }
        for (int j = i+1; j < n; j++){
            double e = a[j][i]/a[i][i];
            for (int k = 0; k < n; k++){
                a[j][k] -= e*a[i][k];
                ans[j][k] -= e*ans[i][k];
            }
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--){
            double e = a[j][i]/a[i][i];
            for (int k = 0; k < n; k++){
                a[j][k] -= e*a[i][k];
                ans[j][k] -= e*ans[i][k];
            }
        }
        for (int j = 0; j < n; j++) {
            ans[i][j] /= a[i][i];
        }
    }
    return ans;
}
double wilandt_method(double lambda, std::vector<std::vector<double>> a_matrix, double eps){
    double lambda_approx = lambda, lambda_approx_new = lambda_approx;
    int n = a_matrix[0].size();
    do {
        lambda_approx = lambda_approx_new;
        std::vector<std::vector<double> > w_matrix = a_matrix;
        for (int i = 0; i < n; ++i) {
            w_matrix[i][i] -= lambda_approx;
        }
        double mu = scalar_method(inverse(w_matrix), a_matrix[0], eps * pow(10, -3));
        lambda_approx_new += 1/mu;
    } while(std::abs(lambda_approx_new - lambda_approx) >= eps);
    return lambda_approx_new;
}

int main() {
    std::cout << "Computational workshop\n Kizeev Danil\n 321 group\n 2019\n";
    std::ifstream file("C:\\Games\\Jacobi method\\data.txt");
    std::cout << "Opening file with start data: \n";
    if(!file){
        std::cout << "File not found";
        exit(0);
    }
    else std::cout << "File has opened successfully\n";
    int n;
    std::cout << "Enter n (matrix size) number:";
    std::cin >> n;
    std::vector<std::vector<double>> x_matrix, main_nums_matrix, a_matrix(n, std::vector<double>(n));
    std::string string_curr;
    for(int i = 0; i < n; ++i){
        getline(file, string_curr);
        std::vector<double> vect_temp = split_string_into_array(string_curr);
        for(int j = 0; j <  n; ++j){
            a_matrix[i][j] = vect_temp[j];
        }
    }
    double eps;
    std::cout << "Enter accuracy criteria(epsilon): "; std::cin >> eps;
    //output
    std::cout << "\nMain matrix: " << std::endl;
    for(int i = 0; i < n; ++i){         //main Matrix
        for(int j = 0; j < n; ++j){
            std::cout << a_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    main_nums_matrix = main_nums(n, a_matrix, eps);
    std::cout << "Main numbers of matrix: ";
    for(int i = 0; i < n; ++i) std::cout << main_nums_matrix[i][i] << " ";
    x_matrix = main_vectors(n, a_matrix, eps);
    std::cout << "Main vectors of matrix: " << std::endl;
    for(int i = 0; i < n; ++i){         //main-vectors Matrix
        for(int j = 0; j < n; ++j){
            std::cout << x_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::vector<double> begin_vector(n);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j) {
            begin_vector[i] += (n - j) * x_matrix[i][j];
        }
    }
    std::cout << "Begin vector is: ";
    for(auto el : begin_vector)
        std::cout << el << " ";
    std::cout << std::endl;
    std::cout << "Enter new epsilon for finding max main number by exponent method: "; std::cin >> eps;
    double lambda = exponent_method(a_matrix, begin_vector, eps);
    std::cout << "Maximum for abs in main numbers founded by exp method is: " << lambda << std::endl;
    std::cout << "Enter new epsilon for finding max main number by scalar method: "; std::cin >> eps;
    lambda = scalar_method(a_matrix, begin_vector, eps);
    std::cout << "Maximum for abs in main numbers founded by scalar method is: " << lambda << std::endl;
    std::cout << "Enter new epsilon for finding opposite spectre board using exponent method: "; std::cin >> eps;
    lambda = spectr_board(lambda, a_matrix, eps);
    std::cout << "Opposite for main numbers founded by opposite board spectre is: " << lambda << std::endl;
    std::cout << "Enter new epsilon for main number clarification using Wilandt method: "; std::cin >> eps;
    lambda = wilandt_method(main_nums_matrix[2][2], a_matrix, eps);
    std::cout << "Your clarified lambda_3 value by Wilandt's method: " << lambda << std::endl;
    file.close();
    return 0;
}