#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

struct Atom {
    char element;
    double x, y, z;
};

// Function to read atomic coordinates from a file
std::vector<Atom> readCoordinates(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<Atom> atoms;
    char element;
    double x, y, z;
    while (file >> element >> x >> y >> z) {
        Atom atom;
        atom.element = element;
        atom.x = x;
        atom.y = y;
        atom.z = z;
        atoms.push_back(atom);
    }

    file.close();
    return atoms;
}

// Function to solve Poisson equation using finite difference method in 2D
void solvePoisson2D(std::vector<std::vector<double> >& phi, const std::vector<Atom>& atoms, double dx, double dy, int max_iter, double tolerance) {
    int Nx = phi.size();
    int Ny = phi[0].size();

    // Set up charge density
    for (const auto& atom : atoms) {
        int i = static_cast<int>(round(atom.x / dx)); // Convert x-coordinate to grid index (using round)
        int j = static_cast<int>(round(atom.y / dy)); // Convert y-coordinate to grid index (using round)
        if (i >= 0 && i < Nx && j >= 0 && j < Ny) {
            phi[i][j] = atom.element == 'O' ? -8.0 : 1.0; // Adjust charge for Oxygen
        }
    }

    // Main iterative loop
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_error = 0.0;

        // Start timing the iteration
        auto start_time = std::chrono::steady_clock::now();

        // Iterate over interior points
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double new_phi = 0.25 * (phi[i-1][j] + phi[i+1][j] + phi[i][j-1] + phi[i][j+1]);
                double error = std::abs(new_phi - phi[i][j]);
                if (error > max_error) {
                    max_error = error;
                }
                phi[i][j] = new_phi;
            }
        }

        // End timing the iteration
        auto end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        // Output the iteration time
        std::cout << "Iteration " << iter + 1 << " time: " << duration.count() << " microseconds" << std::endl;

        // Check convergence
        if (max_error < tolerance) {
            std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;
            return;
        }
    }

    std::cerr << "Warning: Did not converge after " << max_iter << " iterations." << std::endl;
}

int main() {
    // Parameters
    const int Nx = 100;          // Number of grid points in x-direction
    const int Ny = 100;          // Number of grid points in y-direction
    const double Lx = 1.0;       // Length of the domain in x-direction
    const double Ly = 1.0;       // Length of the domain in y-direction
    const double dx = Lx / (Nx - 1); // Grid spacing in x-direction
    const double dy = Ly / (Ny - 1); // Grid spacing in y-direction
    const int max_iter = 10000;   // Maximum number of iterations
    const double tolerance = 1e-6; // Convergence criterion

    // Read atomic coordinates from file
    std::vector<Atom> atoms = readCoordinates("H2O.txt");

    // Initialize arrays
    std::vector<std::vector<double> > phi(Nx, std::vector<double>(Ny, 0.0)); // Electrostatic potential

    // Open a file for writing
    std::ofstream outputFile("poission_2D_finitediff.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        return EXIT_FAILURE;
    }

    // Redirect std::cout to the file stream
    std::streambuf* originalCoutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(outputFile.rdbuf());

    // Solve Poisson equation
    solvePoisson2D(phi, atoms, dx, dy, max_iter, tolerance);

    // Restore std::cout
    std::cout.rdbuf(originalCoutBuffer);

    // Close the output file
    outputFile.close();

    std::cout << "Output written to poission_2D_finitediff.txt" << std::endl;

    return 0;

}

