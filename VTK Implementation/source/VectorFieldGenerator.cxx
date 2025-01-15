#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkMath.h>
#include <random>
#include <cmath>

int main() {
    // Grid dimensions and bounds
    const int nx = 21; // Points along x
    const int ny = 21; // Points along y
    const int nt = 20; // Time steps (ensemble members)

    const double x_min = -1.0, x_max = 1.0;
    const double y_min = -1.0, y_max = 1.0;
    const double t_min = 0.0, t_max = 0.75 * vtkMath::Pi();

    const double dx = (x_max - x_min) / (nx - 1);
    const double dy = (y_max - y_min) / (ny - 1);
    const double dt = (t_max - t_min) / (nt - 1);

    // Random noise generator
    std::default_random_engine generator;
    std::uniform_real_distribution<double> noise(-0.1, 0.1);

    // Create structured grid
    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> vectors = vtkSmartPointer<vtkDoubleArray>::New();

    vectors->SetNumberOfComponents(3); // Vector (u, v, w)
    vectors->SetName("VectorField");

    // Populate the structured grid
    for (int t = 0; t < nt; ++t) {
        double currentTime = t_min + t * dt;

        for (int j = 0; j < ny; ++j) {
            double y = y_min + j * dy;

            for (int i = 0; i < nx; ++i) {
                double x = x_min + i * dx;

                // Base vector components
                double u0 = sin(x);
                double v0 = sin(y);
                double w0 = 0.5;

                // Time-evolved vector field
                double ut = u0 * cos(currentTime) - v0 * sin(currentTime);
                double vt = u0 * sin(currentTime) + v0 * cos(currentTime);
                double wt = w0;

                // Add noise
                ut += noise(generator);
                vt += noise(generator);
                wt += noise(generator);

                // Add point and vector
                points->InsertNextPoint(x, y, currentTime);
                double vector[3] = { ut, vt, wt };
                vectors->InsertNextTuple(vector);
            }
        }
    }

    // Set up the structured grid
    grid->SetDimensions(nx, ny, nt);
    grid->SetPoints(points);
    grid->GetPointData()->SetVectors(vectors);

    // Write to a .vtk file
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
    writer->SetFileName("vector_field.vtk");
    writer->SetInputData(grid);
    writer->Write();

    std::cout << "Generated vector_field.vtk successfully." << std::endl;
    return EXIT_SUCCESS;
}
