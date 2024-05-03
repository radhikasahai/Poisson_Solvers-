#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdlib>

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

// Function to solve Poisson equation using Jacobi iteration in 2D
void solvePoisson2D_Jacobi(std::vector<std::vector<double> >& phi, const std::vector<Atom>& atoms, double dx, double dy, int max_iter, double tolerance) {
    int Nx = phi.size();
    int Ny = phi[0].size();

    // Open the output file
    std::ofstream outputFile("poisson_jacobi.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output file." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Main iterative loop
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_error = 0.0;

        // Start timing the iteration
        auto start_time = std::chrono::steady_clock::now();

        // Temporary array to store the updated values
        std::vector<std::vector<double> > phi_new(Nx, std::vector<double>(Ny, 0.0));

        // Iterate over interior points
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double charge_density = 0.0;
                for (const auto& atom : atoms) {
                    double dx_atom = i * dx - atom.x;
                    double dy_atom = j * dy - atom.y;
                    double dz_atom = -atom.z; // Assuming the plane of the molecule is in the xy-plane
                    double r_squared = dx_atom * dx_atom + dy_atom * dy_atom + dz_atom * dz_atom;
                    charge_density += atom.element == 'O' ? 8.0 / sqrt(r_squared) : 1.0 / sqrt(r_squared);
                }
                phi_new[i][j] = 0.25 * (charge_density * dx * dy + phi[i - 1][j] + phi[i + 1][j] + phi[i][j - 1] + phi[i][j + 1]);
                double error = std::abs(phi_new[i][j] - phi[i][j]);
                if (error > max_error) {
                    max_error = error;
                }
            }
        }

        // Update potential with new values
        phi = phi_new;

        // End timing the iteration
        auto end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        // Output the iteration time to the file
        outputFile << "Iteration " << iter + 1 << " time: " << duration.count() << " microseconds" << std::endl;

        // Check convergence
        if (max_error < tolerance) {
            outputFile << "Converged after " << iter + 1 << " iterations." << std::endl;
            outputFile.close();
            return;
        }
    }

    outputFile << "Warning: Did not converge after " << max_iter << " iterations." << std::endl;
    outputFile.close();
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
    const double tolerance = 1e-5; // Convergence criterion

    // Read atomic coordinates from file
    std::vector<Atom> atoms = readCoordinates("H2O.txt");

    // Initialize arrays
    std::vector<std::vector<double> > phi(Nx, std::vector<double>(Ny, 0.0)); // Electrostatic potential

    // Solve Poisson equation
    solvePoisson2D_Jacobi(phi, atoms, dx, dy, max_iter, tolerance);

    std::cout << "Output written to poisson_jacobi.txt" << std::endl;

    return 0;
}
