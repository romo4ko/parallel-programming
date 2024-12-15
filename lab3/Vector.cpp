#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <thread>
#include <iomanip>
#include <future>

template <typename T>
class Vector {
private:
    size_t length;
    T* elements;
    bool initialized;

    // Helper function to check index bounds
    void validate_index(size_t index) const {
        if (index >= length) {
            throw std::out_of_range("Index out of range");
        }
    }

    template <typename Func>
    T parallel_reduce(Func func, size_t thread_count) const {
        ensure_initialized();
        if (length == 0) {
            return 0;
        }
        std::vector<std::future<T>> futures;
        size_t chunk_size = length / thread_count;

        for (size_t i = 0; i < thread_count; ++i) {
            size_t start = i * chunk_size;
            size_t end = (i == thread_count - 1) ? length : (i + 1) * chunk_size;
            futures.push_back(std::async(std::launch::async, func, start, end));
        }

        T result = 0;
        for (auto& future : futures) {
            result += future.get();
        }
        return result;
    }

    template <typename Func>
    auto parallel_find_min_max(Func func, size_t thread_count) const {
        ensure_initialized();
        if (length == 0) {
            return std::make_pair(std::make_pair(static_cast<T>(0), static_cast<size_t>(0)), std::make_pair(static_cast<T>(0), static_cast<size_t>(0)));
        }
        std::vector<std::future<std::pair<T, size_t>>> futures;
        size_t chunk_size = length / thread_count;

        for (size_t i = 0; i < thread_count; ++i) {
            size_t start = i * chunk_size;
            size_t end = (i == thread_count - 1) ? length : (i + 1) * chunk_size;
            futures.push_back(std::async(std::launch::async, func, start, end));
        }

        std::pair<T, size_t> min_result = futures[0].get();
        std::pair<T, size_t> max_result = futures[0].get();
        for (size_t i = 1; i < futures.size(); ++i) {
            auto result = futures[i].get();
            if (result.first < min_result.first) {
                min_result = result;
            }
            if (result.second > max_result.second) {
                max_result = result;
            }
        }

        return std::make_pair(min_result, max_result);
    }

public:
    // Constructor
    Vector(size_t size) : length(size), elements(nullptr), initialized(false) {
        if (size > 0) {
            try {
                elements = new T[length];
            } catch (const std::bad_alloc& e) {
                std::cerr << "Memory allocation failed: " << e.what() << std::endl;
                throw;
            }
        } else {
            throw std::invalid_argument("Vector size must be positive");
        }
    }

    // Destructor
    ~Vector() {
        delete[] elements;
    }

    // Ensure vector is initialized
    void ensure_initialized() const {
        if (!initialized) {
            throw std::runtime_error("Vector is not initialized");
        }
    }

    // Initialize with a constant value
    std::chrono::duration<double> initialize(const T& value) {
        auto start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < length; ++i) {
            elements[i] = value;
        }
        initialized = true;
        auto end = std::chrono::high_resolution_clock::now();
        return end - start;
    }

    // Access operator with bounds checking
    T& operator[](size_t index) {
        ensure_initialized();
        validate_index(index);
        return elements[index];
    }

    const T& operator[](size_t index) const {
        ensure_initialized();
        validate_index(index);
        return elements[index];
    }

    // Initialize with random values
    std::chrono::duration<double> initialize_random(T min, T max) {
        auto start = std::chrono::high_resolution_clock::now();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<T> dist(min, max);

        for (size_t i = 0; i < length; ++i) {
            elements[i] = dist(gen);
        }
        initialized = true;
        auto end = std::chrono::high_resolution_clock::now();
        return end - start;
    }

    // Export to file
    std::chrono::duration<double> export_to_file(const std::string& filename) {
        ensure_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        std::ofstream file(filename, std::ios::binary);
        if (file.is_open()) {
            file.write(reinterpret_cast<char*>(elements), sizeof(T) * length);
            file.close();
        } else {
            throw std::runtime_error("Unable to open file");
        }
        auto end = std::chrono::high_resolution_clock::now();
        return end - start;
    }

    // Import from file
    std::chrono::duration<double> import_from_file(const std::string& filename) {
        if (elements == nullptr) {
            elements = new T[length];
        }
        auto start = std::chrono::high_resolution_clock::now();
        std::ifstream file(filename, std::ios::binary);
        if (file.is_open()) {
            file.read(reinterpret_cast<char*>(elements), sizeof(T) * length);
            file.close();
            initialized = true;
        } else {
            throw std::runtime_error("Unable to open file");
        }
        auto end = std::chrono::high_resolution_clock::now();
        return end - start;
    }

    // Methods for finding min/max, average, sum, norms, and dot product
    std::pair<std::pair<T, size_t>, std::pair<T, size_t>> parallel_find_min_max(size_t thread_count) const {
        return parallel_find_min_max([this](size_t start, size_t end) {
            if (start >= end) {
                return std::make_pair(std::numeric_limits<T>::max(), size_t(-1));
            }

            T min_val = elements[start];
            T max_val = elements[start];
            size_t min_index = start;
            size_t max_index = start;

            for (size_t i = start + 1; i < end; ++i) {
                min_val = std::min(min_val, elements[i]);
                max_val = std::max(max_val, elements[i]);
                min_index = (elements[i] < min_val) ? i : min_index;
                max_index = (elements[i] > max_val) ? i : max_index;
            }

            return std::make_pair(std::make_pair(min_val, min_index), std::make_pair(max_val, max_index));
        }, thread_count);
    }

    // Parallel Euclidean norm
    double parallel_euclidean_norm(size_t thread_count) const {
        return std::sqrt(parallel_reduce([this](size_t start, size_t end) {
            double local_sum_of_squares = 0;
            for (size_t i = start; i < end; ++i) {
                local_sum_of_squares += std::pow(elements[i], 2);
            }
            return local_sum_of_squares;
        }, thread_count));
    }

    // Average value
    T average() const {
        ensure_initialized();
        if (length == 0) {
            throw std::runtime_error("Vector is empty, can't calculate average");
        }
        return std::accumulate(elements, elements + length, static_cast<T>(0)) / length;
    }

    // Sum of elements
    T sum() const {
        ensure_initialized();
        return std::accumulate(elements, elements + length, static_cast<T>(0));
    }

    T parallel_sum(size_t thread_count) const {
        return parallel_reduce([this](size_t start, size_t end) {
            T local_sum = 0;
            for (size_t i = start; i < end; ++i) {
                local_sum += elements[i];
            }
            return local_sum;
        }, thread_count);
    }

    // Parallel average
    T parallel_average(size_t thread_count) const {
        if (length == 0) {
            return 0;
        }
        return parallel_sum(thread_count) / length;
    }

    // Euclidean norm
    double euclidean_norm() const {
        ensure_initialized();
        double sum_of_squares = 0.0;
        for (size_t i = 0; i < length; ++i) {
            sum_of_squares += std::pow(elements[i], 2);
        }
        return std::sqrt(sum_of_squares);
    }

    // Manhattan norm
    T manhattan_norm() const {
        ensure_initialized();
        T sum = 0;
        for (size_t i = 0; i < length; ++i) {
            sum += std::abs(elements[i]);
        }
        return sum;
    }

    // Dot product
    T dot_product(const Vector<T>& other) const {
        ensure_initialized();
        other.ensure_initialized();
        if (length != other.length) {
            throw std::invalid_argument("Vectors must have the same size for dot product");
        }

        T result = 0;
        for (size_t i = 0; i < length; ++i) {
            result += elements[i] * other.elements[i];
        }
        return result;
    }

    // Parallel Manhattan norm
    T parallel_manhattan_norm(size_t thread_count) const {
        return parallel_reduce([this](size_t start, size_t end) {
            T local_sum = 0;
            for (size_t i = start; i < end; ++i) {
                local_sum += std::abs(elements[i]);
            }
            return local_sum;
        }, thread_count);
    }

    // Parallel dot product
    T parallel_dot_product(const Vector<T>& other, size_t thread_count) const {
        ensure_initialized();
        other.ensure_initialized();
        if (length != other.length) {
            throw std::invalid_argument("Vectors must have the same size for dot product");
        }
        return parallel_reduce([this, &other](size_t start, size_t end) {
            T local_result = 0;
            for (size_t i = start; i < end; ++i) {
                local_result += elements[i] * other.elements[i];
            }
            return local_result;
        }, thread_count);
    }

    size_t size() const { return length; } // Return vector size
};

int main() {
    try {
        size_t vector_size = 10000000;
        Vector<double> vec(vector_size);
        vec.initialize_random(-10.0, 10.0);

        auto measure_time = [&](const auto& func, const std::string& name) {
            auto start = std::chrono::high_resolution_clock::now();
            auto result = func();
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            std::cout << name << " time: " << duration << "ms\n";
            return std::make_pair(duration, result);
        };

        unsigned int thread_count = std::thread::hardware_concurrency();
        if (thread_count == 0) {
            thread_count = 2;
        }

        std::cout << "Number of threads: " << thread_count << std::endl;

        std::ofstream output_file("results.txt");
        if (!output_file.is_open()) {
            throw std::runtime_error("Unable to open output file");
        }

        output_file << std::fixed << std::setprecision(6);

        int iterations = 5;
        for (int i = 0; i < iterations; ++i) {
            output_file << "Iteration " << i + 1 << ":\n";

            output_file << "Sequential tests:\n";
            output_file << "Sum: " << measure_time([&]() { return vec.sum(); }, "Sequential sum").first << "\n";
            output_file << "Average: " << measure_time([&]() { return vec.average(); }, "Sequential average").first << "\n";
            Vector<double> vec2 = vec;
            output_file << "Dot product: " << measure_time([&]() { return vec.dot_product(vec2); }, "Sequential dot product").first << "\n";

            output_file << "\nParallel tests with " << thread_count << " threads:\n";
            output_file << "Sum: " << measure_time([&]() { return vec.parallel_sum(thread_count); }, "Parallel sum").first << "\n";
            output_file << "Average: " << measure_time([&]() { return vec.parallel_average(thread_count); }, "Parallel average").first << "\n";
            output_file << "Dot product: " << measure_time([&]() { return vec.parallel_dot_product(vec2, thread_count); }, "Parallel dot product").first << "\n";
            output_file << "\n";
        }
        output_file.close();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}