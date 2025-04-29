#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <vtkm/Math.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <Windows.h>
#include <commdlg.h>

struct GlyphParams {
    vtkm::FloatDefault base_width = 0.8f;
    vtkm::FloatDefault base_depth = 0.4f;
    vtkm::FloatDefault base_height = 0.5f;
    vtkm::FloatDefault middle_height = 1.0f;
    vtkm::FloatDefault middle_radius = 0.2f;
    vtkm::FloatDefault arrow_height = 0.5f;
    vtkm::FloatDefault arrow_radius = 0.4f;
    vtkm::Id resolution = 24;
    vtkm::FloatDefault superellipse_n = 4.0f;
    vtkm::FloatDefault insertion_depth = 0.3f;
    vtkm::FloatDefault scale_factor = 0.1f;
};

struct VectorData {
    vtkm::Vec3f position;
    vtkm::Vec3f direction;
    float depth;
    vtkm::Vec3f principalDir0;
    vtkm::Vec3f principalDir1;
    float theta0;
    float theta1;
    float timestamp;
};

// File dialog for CSV selection
std::string OpenFileDialog() {
    OPENFILENAME ofn;
    char szFile[260] = { 0 };
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = nullptr;
    ofn.lpstrFile = szFile;
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "CSV Files\0*.csv\0All Files\0*.*\0";
    ofn.nFilterIndex = 1;
    ofn.lpstrFileTitle = nullptr;
    ofn.nMaxFileTitle = 0;
    ofn.lpstrInitialDir = nullptr;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileName(&ofn)) {
        return std::string(szFile);
    }
    return "";
}

// Read CSV with zero-vector handling
std::vector<VectorData> ReadCSV(const std::string& filename, float selected_timestamp) {
    std::vector<VectorData> data;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return data;
    }

    std::getline(file, line); // Skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<float> values;
        std::string value;

        while (std::getline(ss, value, ',')) {
            try {
                values.push_back(std::stof(value));
            } catch (...) {
                values.push_back(0.0f);
            }
        }

        if (values.size() >= 7) {
            VectorData vec;
            vec.position = vtkm::Vec3f(values[0], values[1], values[2]);
            vec.direction = vtkm::Vec3f(values[3], values[4], values[5]);
            float mag = vtkm::Magnitude(vec.direction);
            vec.direction = mag > 1e-6 ? vec.direction / mag : vtkm::Vec3f(0.0f, 0.0f, 1.0f);
            vec.timestamp = values[6];
            if (std::abs(vec.timestamp - selected_timestamp) < 1e-6) {
                data.push_back(vec);
            }
        }
    }
    file.close();
    return data;
}

// Generate ensemble: 9 vectors (7 normal, 2 elongated)
void GenerateEnsemble(std::vector<vtkm::Vec3f>& ensemble, const vtkm::Vec3f& base_dir, std::mt19937& gen) {
    std::normal_distribution<float> dist(0.0f, 0.1f);
    ensemble.clear();
    ensemble.push_back(base_dir);

    for (int i = 0; i < 8; ++i) {
        vtkm::Vec3f perturb(dist(gen), dist(gen), dist(gen));
        vtkm::Vec3f dir = base_dir + perturb;
        float mag = vtkm::Magnitude(dir);
        dir = mag > 1e-6 ? dir / mag : vtkm::Vec3f(0.0f, 0.0f, 1.0f);
        if (i >= 6) {
            dir = dir * 2.0f; // Elongate last two
        }
        ensemble.push_back(dir);
    }
}

// Compute PCA with numerical stability
void ComputePCA(const std::vector<vtkm::Vec3f>& directions, vtkm::Vec3f& principalDir0, vtkm::Vec3f& principalDir1, float& theta0, float& theta1) {
    vtkm::Vec3f meanDir(0.0f);
    for (const auto& dir : directions) {
        meanDir += dir;
    }
    meanDir = meanDir / static_cast<float>(directions.size());
    float mag = vtkm::Magnitude(meanDir);
    meanDir = mag > 1e-6 ? meanDir / mag : vtkm::Vec3f(0.0f, 0.0f, 1.0f);

    principalDir0 = meanDir;
    vtkm::Vec3f zAxis(0.0f, 0.0f, 1.0f);
    principalDir1 = vtkm::Normal(vtkm::Cross(meanDir, zAxis));
    if (vtkm::Magnitude(principalDir1) < 1e-6) {
        vtkm::Vec3f xAxis(1.0f, 0.0f, 0.0f);
        principalDir1 = vtkm::Normal(vtkm::Cross(meanDir, xAxis));
    }

    theta0 = 0.0f;
    theta1 = 0.0f;
    for (const auto& dir : directions) {
        float dot0 = std::min(1.0f, std::max(-1.0f, vtkm::Dot(dir, principalDir0)));
        float dot1 = std::min(1.0f, std::max(-1.0f, vtkm::Dot(dir, principalDir1)));
        float angle0 = std::acos(dot0);
        float angle1 = std::acos(dot1);
        theta0 = std::max(theta0, angle0);
        theta1 = std::max(theta1, angle1);
    }
}

// Compute vector depth
float ComputeVectorDepth(const vtkm::Vec3f& vec, const std::vector<vtkm::Vec3f>& ensemble) {
    vtkm::Vec3f mean(0.0f);
    for (const auto& v : ensemble) {
        mean += v;
    }
    mean = mean / static_cast<float>(ensemble.size());
    float dist = vtkm::Magnitude(vec - mean);
    float maxDist = 0.0f;
    for (const auto& v : ensemble) {
        maxDist = std::max(maxDist, vtkm::Magnitude(v - mean));
    }
    return maxDist > 0.0f ? (1.0f - dist / maxDist) : 1.0f;
}

// Rotate point around axis
vtkm::Vec3f RotatePoint(const vtkm::Vec3f& point, const vtkm::Vec3f& axis, float angle) {
    float c = std::cos(angle);
    float s = std::sin(angle);
    float t = 1.0f - c;
    vtkm::Vec3f normalizedAxis = vtkm::Normal(axis);

    float x = point[0], y = point[1], z = point[2];
    float u = normalizedAxis[0], v = normalizedAxis[1], w = normalizedAxis[2];

    return vtkm::Vec3f(
        (t * u * u + c) * x + (t * u * v - s * w) * y + (t * u * w + s * v) * z,
        (t * u * v + s * w) * x + (t * v * v + c) * y + (t * v * w - s * u) * z,
        (t * u * w - s * v) * x + (t * v * w + s * u) * y + (t * w * w + c) * z
    );
}

// Create glyph (unchanged from reference)
void CreateGlyph(const VectorData& vecData, const GlyphParams& params,
                std::vector<vtkm::Vec3f>& points,
                std::vector<vtkm::UInt8>& shapes,
                std::vector<vtkm::IdComponent>& numIndices,
                std::vector<vtkm::Id>& connectivity) {
    float magnitude = vtkm::Magnitude(vecData.direction);
    float scale = params.scale_factor * (1.0f + vecData.depth);

    float h = params.base_height;
    float deltaH = 0.1f;
    float r0 = (h + deltaH) * std::tan(vecData.theta0 / 2.0f);
    float r1 = r0 * (vecData.theta1 / vecData.theta0);
    if (vecData.theta0 == 0.0f) r1 = r0;

    float scaleX = r0 / (params.base_width / 2.0f) * scale;
    float scaleY = r1 / (params.base_depth / 2.0f) * scale;

    vtkm::Vec3f zAxis(0.0f, 0.0f, 1.0f);
    vtkm::Vec3f axis = vtkm::Cross(zAxis, vecData.principalDir0);
    float angle = std::acos(vtkm::Dot(zAxis, vecData.principalDir0));
    if (vtkm::Magnitude(axis) < 1e-6) {
        axis = vtkm::Vec3f(1.0f, 0.0f, 0.0f);
    }
    axis = vtkm::Normal(axis);

    vtkm::Id base_start = static_cast<vtkm::Id>(points.size());
    vtkm::FloatDefault a = params.base_width / 2;
    vtkm::FloatDefault b = params.base_depth / 2;
    vtkm::Id segments = params.resolution;

    // Base of the glyph
    for (vtkm::Id i = 0; i < segments; i++) {
        vtkm::FloatDefault theta = 2.0f * vtkm::Pi() * i / segments;
        vtkm::FloatDefault cosTheta = std::cos(theta);
        vtkm::FloatDefault sinTheta = std::sin(theta);

        vtkm::FloatDefault x = a * std::pow(std::abs(cosTheta), 2.0f / params.superellipse_n) * std::copysign(1.0f, cosTheta);
        vtkm::FloatDefault y = b * std::pow(std::abs(sinTheta), 2.0f / params.superellipse_n) * std::copysign(1.0f, sinTheta);

        x *= scaleX;
        y *= scaleY;

        vtkm::Vec3f point(x, y, 0.0f);
        point = RotatePoint(point, axis, angle);
        points.push_back(vecData.position + point);

        point = vtkm::Vec3f(x, y, params.base_height * scale);
        point = RotatePoint(point, axis, angle);
        points.push_back(vecData.position + point);
    }

    // Connect base points
    for (vtkm::Id i = 0; i < segments; i++) {
        vtkm::Id next = (i + 1) % segments;
        shapes.push_back(vtkm::CELL_SHAPE_QUAD);
        numIndices.push_back(4);
        connectivity.push_back(base_start + 2 * i);
        connectivity.push_back(base_start + 2 * next);
        connectivity.push_back(base_start + 2 * next + 1);
        connectivity.push_back(base_start + 2 * i + 1);
    }

    shapes.push_back(vtkm::CELL_SHAPE_POLYGON);
    numIndices.push_back(segments);
    for (vtkm::Id i = 0; i < segments; i++) {
        connectivity.push_back(base_start + 2 * i);
    }

    shapes.push_back(vtkm::CELL_SHAPE_POLYGON);
    numIndices.push_back(segments);
    for (vtkm::Id i = 0; i < segments; i++) {
        connectivity.push_back(base_start + 2 * i + 1);
    }

    // Middle part (cone, scaled)
    vtkm::FloatDefault cone_base_width = params.middle_radius * 2;
    vtkm::FloatDefault cone_height = params.middle_height;
    vtkm::Id cone_base_start = static_cast<vtkm::Id>(points.size());

    for (vtkm::Id i = 0; i < segments; i++) {
        vtkm::FloatDefault theta = 2.0f * vtkm::Pi() * i / segments;
        vtkm::FloatDefault cosTheta = std::cos(theta);
        vtkm::FloatDefault sinTheta = std::sin(theta);
        vtkm::FloatDefault x = cone_base_width / 2 * cosTheta;
        vtkm::FloatDefault y = cone_base_width * params.base_depth / 2 * sinTheta;

        x *= scaleX;
        y *= scaleY;

        vtkm::Vec3f point(x, y, params.base_height * scale);
        point = RotatePoint(point, axis, angle);
        points.push_back(vecData.position + point);
    }

    for (vtkm::Id i = 0; i < segments; i++) {
        vtkm::Id next = (i + 1) % segments;
        shapes.push_back(vtkm::CELL_SHAPE_QUAD);
        numIndices.push_back(4);
        connectivity.push_back(cone_base_start + i);
        connectivity.push_back(cone_base_start + next);
        connectivity.push_back(cone_base_start + next + 1);
        connectivity.push_back(cone_base_start + i + 1);
    }

    vtkm::Id cone_tip = static_cast<vtkm::Id>(points.size());
    vtkm::Vec3f tip(0.0f, 0.0f, (params.base_height + cone_height + params.insertion_depth) * scale);
    tip = RotatePoint(tip, axis, angle);
    points.push_back(vecData.position + tip);

    for (vtkm::Id i = 0; i < segments; i++) {
        vtkm::Id next = (i + 1) % segments;
        shapes.push_back(vtkm::CELL_SHAPE_TRIANGLE);
        numIndices.push_back(3);
        connectivity.push_back(cone_base_start + i);
        connectivity.push_back(cone_base_start + next);
        connectivity.push_back(cone_tip);
    }

    // Arrow part (scaled)
    vtkm::Id arrow_cone_base_start = static_cast<vtkm::Id>(points.size());
    vtkm::FloatDefault arrow_cone_base_width = params.arrow_radius * 2;
    vtkm::FloatDefault arrow_cone_height = params.arrow_height;

    for (vtkm::Id i = 0; i < segments; i++) {
        vtkm::FloatDefault theta = 2.0f * vtkm::Pi() * i / segments;
        vtkm::FloatDefault cosTheta = std::cos(theta);
        vtkm::FloatDefault sinTheta = std::sin(theta);
        vtkm::FloatDefault x = (arrow_cone_base_width / 2) *
            std::pow(std::abs(cosTheta), 2.0f / params.superellipse_n) *
            std::copysign(1.0f, cosTheta);
        vtkm::FloatDefault y = (arrow_cone_base_width * params.base_depth / 2) *
            std::pow(std::abs(sinTheta), 2.0f / params.superellipse_n) *
            std::copysign(1.0f, sinTheta);

        x *= scaleX;
        y *= scaleY;

        vtkm::Vec3f point(x, y, (params.base_height + params.middle_height) * scale);
        point = RotatePoint(point, axis, angle);
        points.push_back(vecData.position + point);
    }

    shapes.push_back(vtkm::CELL_SHAPE_POLYGON);
    numIndices.push_back(segments);
    for (vtkm::Id i = 0; i < segments; i++) {
        connectivity.push_back(arrow_cone_base_start + i);
    }

    vtkm::Id arrow_cone_tip = static_cast<vtkm::Id>(points.size());
    vtkm::Vec3f arrowTip(0.0f, 0.0f, (params.base_height + params.middle_height + arrow_cone_height) * scale);
    arrowTip = RotatePoint(arrowTip, axis, angle);
    points.push_back(vecData.position + arrowTip);

    for (vtkm::Id i = 0; i < segments; i++) {
        vtkm::Id next = (i + 1) % segments;
        shapes.push_back(vtkm::CELL_SHAPE_TRIANGLE);
        numIndices.push_back(3);
        connectivity.push_back(arrow_cone_base_start + i);
        connectivity.push_back(arrow_cone_base_start + next);
        connectivity.push_back(arrow_cone_tip);
    }
}

// Create raw vectors dataset
vtkm::cont::DataSet CreateVectorDataSet(const std::vector<VectorData>& vectorData, const std::vector<std::vector<vtkm::Vec3f>>& ensembles) {
    std::vector<vtkm::Vec3f> points;
    std::vector<vtkm::UInt8> shapes;
    std::vector<vtkm::IdComponent> numIndices;
    std::vector<vtkm::Id> connectivity;

    for (size_t i = 0; i < vectorData.size(); ++i) {
        vtkm::Id start = static_cast<vtkm::Id>(points.size());
        const auto& ensemble = ensembles[i];

        for (const auto& dir : ensemble) {
            points.push_back(vectorData[i].position);
            points.push_back(vectorData[i].position + dir * 0.5f);
            shapes.push_back(vtkm::CELL_SHAPE_LINE);
            numIndices.push_back(2);
            connectivity.push_back(start);
            connectivity.push_back(start + 1);
            start += 2;
        }
    }

    vtkm::cont::DataSetBuilderExplicit dsBuilder;
    return dsBuilder.Create(points, shapes, numIndices, connectivity);
}

int main() {
    GlyphParams params;
    std::random_device rd;
    std::mt19937 gen(rd());

    // Prompt for timestamp
    float selected_timestamp = 0.0f;
    std::cout << "Enter timestamp to visualize (default 0): ";
    std::string input;
    std::getline(std::cin, input);
    if (!input.empty()) {
        try {
            selected_timestamp = std::stof(input);
        } catch (...) {
            std::cerr << "Invalid timestamp, using 0.0" << std::endl;
        }
    }

    // Load dataset
    std::cout << "Please select a CSV file..." << std::endl;
    std::string filename = OpenFileDialog();
    if (filename.empty()) {
        std::cerr << "No file selected!" << std::endl;
        return 1;
    }
    std::cout << "Selected file: " << filename << std::endl;
    std::vector<VectorData> vectorData = ReadCSV(filename, selected_timestamp);
    if (vectorData.empty()) {
        std::cerr << "No valid data found for timestamp " << selected_timestamp << "!" << std::endl;
        return 1;
    }
    std::cout << "Found " << vectorData.size() << " vectors for timestamp " << selected_timestamp << "." << std::endl;

    // Generate ensembles and compute metrics
    std::vector<std::vector<vtkm::Vec3f>> ensembles(vectorData.size());
    for (size_t i = 0; i < vectorData.size(); ++i) {
        GenerateEnsemble(ensembles[i], vectorData[i].direction, gen);
        ComputePCA(ensembles[i], vectorData[i].principalDir0, vectorData[i].principalDir1, vectorData[i].theta0, vectorData[i].theta1);
        vectorData[i].depth = ComputeVectorDepth(vectorData[i].direction, ensembles[i]);
    }

    // Save raw vectors
    vtkm::cont::DataSet vectorDataSet = CreateVectorDataSet(vectorData, ensembles);
    vtkm::io::VTKDataSetWriter vectorWriter("raw_vectors.vtk");
    vectorWriter.WriteDataSet(vectorDataSet);
    std::cout << "Raw vectors saved to raw_vectors.vtk" << std::endl;

    // Generate glyphs
    std::vector<vtkm::Vec3f> allPoints;
    std::vector<vtkm::UInt8> allShapes;
    std::vector<vtkm::IdComponent> allNumIndices;
    std::vector<vtkm::Id> allConnectivity;

    for (size_t i = 0; i < vectorData.size(); i++) {
        CreateGlyph(vectorData[i], params, allPoints, allShapes, allNumIndices, allConnectivity);
        std::cout << "Added glyph at position " << i + 1 << " of " << vectorData.size() << std::endl;
    }

    // Save glyphs
    vtkm::cont::DataSetBuilderExplicit dsBuilder;
    vtkm::cont::DataSet combinedDataset = dsBuilder.Create(allPoints, allShapes, allNumIndices, allConnectivity);
    std::string outputFilename = "glyphs_timestamp_" + std::to_string(static_cast<int>(selected_timestamp)) + ".vtk";
    vtkm::io::VTKDataSetWriter writer(outputFilename);
    writer.WriteDataSet(combinedDataset);
    
    std::cout << "Glyphs saved to: " << outputFilename << std::endl;
    return 0;
}