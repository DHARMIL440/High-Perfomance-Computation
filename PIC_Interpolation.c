#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <mpi.h>

// Structure to store coordinates of points
typedef struct {
    double x, y;
} Point;

// Global variables
int gridWidth, gridHeight, numX, numY;
int totalPoints, numThreads, maxIterations;
double stepX, stepY;

// Function to initialize points with random values
void initializePoints(Point *points, int numPoints) {
    for (int i = 0; i < numPoints; i++) {
        points[i].x = (double) rand() / RAND_MAX;
        points[i].y = (double) rand() / RAND_MAX;
    }
}

// Function to read points from file
void loadPointsFromFile(FILE *file, Point *points, int numPoints) {
    for (int i = 0; i < numPoints; i++) {
        fread(&points[i].x, sizeof(double), 1, file);
        fread(&points[i].y, sizeof(double), 1, file);
    }
}

// Function to output the mesh to a file
void writeMeshToFile(double *mesh) {
    FILE *file = fopen("Mesh.out", "w");
    if (file == NULL) {
        printf("Error: Could not create Mesh.out\n");
        exit(1);
    }

    for (int i = 0; i < gridHeight; i++) {
        for (int j = 0; j < gridWidth - 1; j++) {
            fprintf(file, "%lf ", mesh[i * gridWidth + j]);
        }
        fprintf(file, "%lf\n", mesh[i * gridWidth + gridWidth - 1]);
    }
    fclose(file);
}

// Cloud-in-a-Cell Interpolation function
void cloudInACellInterpolation(double *mesh, Point *points, int numPoints) {
    memset(mesh, 0, gridWidth * gridHeight * sizeof(double));
    double *threadPrivateMesh = (double *) calloc(numThreads * gridWidth * gridHeight, sizeof(double));

    #pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < numPoints; i++) {
        int threadID = omp_get_thread_num();
        double px = points[i].x;
        double py = points[i].y;

        double weight = 1.0;

        int gridX = (int)(px / stepX);
        int gridY = (int)(py / stepY);

        double localX = px - gridX * stepX;
        double localY = py - gridY * stepY;

        int p1 = gridY * gridWidth + gridX;
        int p2 = gridY * gridWidth + (gridX + 1);
        int p3 = (gridY + 1) * gridWidth + gridX;
        int p4 = (gridY + 1) * gridWidth + (gridX + 1);

        double area1 = (stepX - localX) * (stepY - localY);
        double area2 = localX * (stepY - localY);
        double area3 = (stepX - localX) * localY;
        double area4 = localX * localY;

        int offset = threadID * gridWidth * gridHeight;
        threadPrivateMesh[offset + p1] += area1 * weight;
        threadPrivateMesh[offset + p2] += area2 * weight;
        threadPrivateMesh[offset + p3] += area3 * weight;
        threadPrivateMesh[offset + p4] += area4 * weight;
    }

    #pragma omp parallel for
    for (int i = 0; i < gridWidth * gridHeight; i++) {
        for (int t = 0; t < numThreads; t++) {
            mesh[i] += threadPrivateMesh[t * gridWidth * gridHeight + i];
        }
    }

    free(threadPrivateMesh);
}

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char inputFile[50];
    double start, end, totalTime = 0.0;

    if (argc != 3) {
        if (rank == 0)
            printf("Usage: %s <input_filename> <num_threads>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    strcpy(inputFile, argv[1]);
    numThreads = atoi(argv[2]);

    FILE *file = NULL;
    if (rank == 0) {
        file = fopen(inputFile, "rb");
        if (!file) {
            printf("Error: Unable to open file %s\n", inputFile);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fread(&numX, sizeof(int), 1, file);
        fread(&numY, sizeof(int), 1, file);
        fread(&totalPoints, sizeof(int), 1, file);
        fread(&maxIterations, sizeof(int), 1, file);
    }

    MPI_Bcast(&numX, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numY, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&totalPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxIterations, 1, MPI_INT, 0, MPI_COMM_WORLD);

    gridWidth = numX + 1;
    gridHeight = numY + 1;
    stepX = 1.0 / numX;
    stepY = 1.0 / numY;

    double *localMesh = (double *) calloc(gridWidth * gridHeight, sizeof(double));
    double *globalMesh = NULL;
    if (rank == 0) globalMesh = (double *) calloc(gridWidth * gridHeight, sizeof(double));

    int pointsPerProc = totalPoints / size;
    int remainder = totalPoints % size;
    int pointsForMe = pointsPerProc + (rank < remainder ? 1 : 0);
    Point *myPoints = (Point *) calloc(pointsForMe, sizeof(Point));

    for (int iter = 0; iter < maxIterations; iter++) {
        Point *allPoints = NULL;
        if (rank == 0) {
            allPoints = (Point *) calloc(totalPoints, sizeof(Point));
            loadPointsFromFile(file, allPoints, totalPoints);
        }

        int *counts = (int *) malloc(size * sizeof(int));
        int *displs = (int *) malloc(size * sizeof(int));
        int offset = 0;
        for (int i = 0; i < size; i++) {
            counts[i] = (i < remainder ? pointsPerProc + 1 : pointsPerProc) * 2;
            displs[i] = offset;
            offset += counts[i];
        }

        MPI_Scatterv(allPoints, counts, displs, MPI_DOUBLE,
                     myPoints, pointsForMe * 2, MPI_DOUBLE,
                     0, MPI_COMM_WORLD);

        free(allPoints);
        free(counts);
        free(displs);

        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        cloudInACellInterpolation(localMesh, myPoints, pointsForMe);

        end = MPI_Wtime();
        totalTime += (end - start);
    }

    MPI_Reduce(localMesh, globalMesh, gridWidth * gridHeight, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        writeMeshToFile(globalMesh);
        printf("Interpolation execution time = %lf seconds\n", totalTime);
        free(globalMesh);
        fclose(file);
    }

    free(localMesh);
    free(myPoints);

    MPI_Finalize();
    return 0;
}
