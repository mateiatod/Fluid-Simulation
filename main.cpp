#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <random>
#include <ctime>
#include <iomanip>
#include <cmath>
#include "lbfgs.h"

class Vector {
public:
    double x, y;

    Vector(double x = 0, double y = 0) : x(x), y(y) {}

    double& operator[](size_t idx) {
        return (idx == 0) ? x : y;
    }

    const double& operator[](size_t idx) const {
        return (idx == 0) ? x : y;
    }

    Vector operator-(const Vector& other) const {
        return Vector(x - other.x, y - other.y);
    }

    Vector operator+(const Vector& other) const {
        return Vector(x + other.x, y + other.y);
    }

    Vector operator*(double t) const {
        return Vector(x * t, y * t);
    }

    double dot(const Vector& other) const {
        return x * other.x + y * other.y;
    }

    double normSquared() const {
        return x * x + y * y;
    }

    bool operator==(const Vector& other) const {
        return (x == other.x) && (y == other.y);
    }

    void print() const {
        std::cout << "(" << x << ", " << y << ")";
    }
};

class Polygon {
public:
    std::vector<Vector> vertices;

    void addVertex(const Vector& vertex) {
        vertices.push_back(vertex);
    }
};

void normalizePolygon(Polygon& polygon) {
    const double minCoord = 0.0;
    const double maxCoord = 500.0;

    for (Vector& vertex : polygon.vertices) {
        vertex.x = std::min(std::max(vertex.x, minCoord), maxCoord);
        vertex.y = std::min(std::max(vertex.y, minCoord), maxCoord);
    }
}

void save_svg(const std::vector<Polygon>& polygons, const std::string& filename, const std::string& fillcol = "none") {
    std::ofstream f(filename);
    if (!f.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\">\n";

    for (const Polygon& polygon : polygons) {
        f << "<g>\n";
        f << "<polygon points=\"";
        for (const Vector& vertex : polygon.vertices) {
            f << vertex.x << "," << vertex.y << " ";
        }
        f << "\"\nfill=\"" << fillcol << "\" stroke=\"black\"/>\n";
        f << "</g>\n";
    }

    f << "</svg>\n";
    f.close();

    std::cout << "SVG file saved: " << filename << std::endl;
}

bool inside(const Vector& X, const Polygon& polygon) {
    size_t n = polygon.vertices.size();
    bool inside = false;
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        if (((polygon.vertices[i].y > X.y) != (polygon.vertices[j].y > X.y)) &&
            (X.x < (polygon.vertices[j].x - polygon.vertices[i].x) * (X.y - polygon.vertices[i].y) / (polygon.vertices[j].y - polygon.vertices[i].y))) {
            inside = !inside;
        }
    }
    return inside;
}

bool inside_power(const Vector& X, const Vector& M, const Vector& Pi, const Vector& Pj, double wi, double wj) {
    Vector M_prime = M + (Pj - Pi) * ((wi - wj) / (2 * (Pi - Pj).normSquared())) ;
    return ((X - M_prime).dot(Pj - Pi) < 0);
}

Vector intersect(const Vector& A, const Vector& B, const Vector& M, const Vector& Pi, const Vector& Pj) {
    Vector intersection;

    Vector AB = B - A;

    double dot1 = (M - A).dot(Pi - Pj);
    double dot2 = AB.dot(Pi - Pj);

    double t = dot1 / dot2;

    intersection = A + AB * t;

    return intersection;
}

Vector intersect_power(const Vector& A, const Vector& B, const Vector& M, const Vector& Pi, const Vector& Pj, double wi, double wj) {
    Vector intersection;

    Vector AB = B - A;
    Vector M_prime = M +  (Pj - Pi) * ((wi - wj) / (2 * (Pi - Pj).normSquared()));

    double dot1 = (M_prime - A).dot(Pi - Pj);
    double dot2 = AB.dot(Pi - Pj);

    double t = dot1 / dot2;

    intersection = A + AB * t;

    return intersection;
}

Polygon ClipPolygon_power(const Polygon& subjectPolygon, const Vector& M, const Vector& Pi, const Vector& Pj, double wi, double wj) {
    Polygon outPolygon;

    size_t numVertices = subjectPolygon.vertices.size();
    for (size_t i = 0; i < numVertices; ++i) {
        Vector curVertex = subjectPolygon.vertices[i];
        Vector prevVertex = subjectPolygon.vertices[(i > 0) ? (i - 1) : (numVertices - 1)];

        bool curInside = inside_power(curVertex, M, Pi, Pj, wi, wj);
        bool prevInside = inside_power(prevVertex, M, Pi, Pj, wi, wj);

        if (curInside) {
            if (!prevInside) {
                Vector intersection = intersect_power(prevVertex, curVertex, M, Pi, Pj, wi, wj);
                outPolygon.vertices.push_back(intersection);
            }
            outPolygon.vertices.push_back(curVertex);
        } else if (prevInside) {
            Vector intersection = intersect_power(prevVertex, curVertex, M, Pi, Pj, wi, wj);
            outPolygon.vertices.push_back(intersection);
        }
    }

    return outPolygon;
}

Polygon gen_cellpower(size_t i, const std::vector<Vector>& P, const std::vector<double>& weights) {
    Polygon initialPolygon;
    Vector Pi = P[i];
    initialPolygon.addVertex(Vector(0, 0));
    initialPolygon.addVertex(Vector(10000, 0));
    initialPolygon.addVertex(Vector(10000, 10000));
    initialPolygon.addVertex(Vector(0, 10000));

    Polygon clippedPolygon = initialPolygon;

    for (size_t j = 0; j < P.size(); ++j) {
        if (P[j] == Pi) {
            continue;
        }

        Vector M = (Pi + P[j]) * 0.5;
        Polygon tempPolygon = ClipPolygon_power(clippedPolygon, M, Pi, P[j], weights[i], weights[j]);
        clippedPolygon = tempPolygon;
    }

    normalizePolygon(clippedPolygon);

    return clippedPolygon;
}

std::vector<Vector> generateRandomPoints(int N) {
    std::vector<Vector> points;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 500.0);

    for (int i = 0; i < N; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        points.push_back(Vector(x, y));
    }

    return points;
}

void save_svg_points(const std::vector<Polygon>& polygons, const std::vector<Vector>& points, const std::string& filename, const std::string& fillcol = "none") {
    std::ofstream f(filename);
    if (!f.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    f << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\">\n";

    f << "<g>\n";
    f << "<circle cx=\"50\" cy=\"50\" r=\"3\" fill=\"red\" />\n";
    for (const auto& point : points) {
        f << "<circle cx=\"" << point.x << "\" cy=\"" << point.y << "\" r=\"3\" fill=\"red\" />\n";
    }
    f << "</g>\n";

    for (const Polygon& polygon : polygons) {
        f << "<g>\n";
        f << "<polygon points=\"";
        for (const Vector& vertex : polygon.vertices) {
            f << vertex.x << "," << vertex.y << " ";
        }
        f << "\"\nfill=\"" << fillcol << "\" stroke=\"black\"/>\n";
        f << "</g>\n";
    }

    f << "</svg>\n";
    f.close();

    std::cout << "SVG file saved: " << filename << std::endl;
}

class PowerDiagram {
private:
    std::vector<Vector> points;
    std::vector<double> weights;
    size_t N;

public:
    PowerDiagram(const std::vector<Vector>& points, const std::vector<double>& weights) : points(points), weights(weights), N(points.size()) {}

    void updateWeights(const std::vector<double>& newWeights) {
        if (newWeights.size() != N) {
            std::cerr << "Error: Number of weights does not match the number of points." << std::endl;
            return;
        }
        weights = newWeights;
    }

    static lbfgsfloatval_t _evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
        PowerDiagram* pd = static_cast<PowerDiagram*>(instance);
        lbfgsfloatval_t fx = 0.0;

        return fx;
    }

    static int _progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls) {
        return 0;
    }

    lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
        return _evaluate(this, x, g, n, step);
    }

    int progress(const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls) {
        return _progress(this, x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    void optimizeLBFGS() {
        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);
        param.max_iterations = 250;

        int ret;
        lbfgsfloatval_t *x = lbfgs_malloc(N);
        if (!x) {
            std::cerr << "ERROR: Failed to allocate a memory block for variables." << std::endl;
            return;
        }

        for (size_t i = 0; i < N; ++i) {
            x[i] = points[i].x;
        }

        ret = lbfgs(N, x, nullptr, &_evaluate, &_progress, this, &param);

        if (ret == LBFGS_SUCCESS) {
            std::cout << "LBFGS optimization terminated with status: SUCCESS." << std::endl;
        } else {
            std::cerr << "LBFGS optimization terminated with status: " << ret << std::endl;
        }

        lbfgs_free(x);
    }
};

double calculateWeight(const Vector& Pi) {
    Vector C(250.0, 250.0);
    double distSquared = (Pi - C).normSquared();
    double weight = exp(-distSquared / 0.02);
    return weight;
}

std::vector<double> calculateWeights(const std::vector<Vector>& points) {
    std::vector<double> weights;
    weights.reserve(points.size());

    for (const auto& point : points) {
        double weight = calculateWeight(point);
        weights.push_back(weight);
    }

    return weights;
}

int main() {
    std::vector<Vector> points = generateRandomPoints(500);
    std::vector<double> weights(points.size(), 1.0);
    
    weights = calculateWeights(points);
    // for (int i = 0; i <= points.size(); ++i){
    //     std::cout<<weights[i]<<" ";
    // }
    // std::cout<<std::endl;

    PowerDiagram pd(points, weights);
    pd.optimizeLBFGS();

    std::vector<Polygon> polygons;
    for (size_t i = 0; i < points.size(); ++i) {
        polygons.push_back(gen_cellpower(i, points, weights));
    }

    save_svg_points(polygons, points, "output.svg");

    return 0;
}
