#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

// GPU Characteristics class
class GpuCharacteristics {
public:
    string Model;
    int MemoryMB;
    int CoreCount;

    // Default constructor
    GpuCharacteristics() : Model("Unknown"), MemoryMB(0), CoreCount(0) {}

    GpuCharacteristics(const string& model, int memoryMB, int coreCount)
        : Model(model), MemoryMB(memoryMB), CoreCount(coreCount) {}

    // Method to display data
    void Display() const {
        cout << "GPU Model: " << Model << ", Memory: " << MemoryMB << " MB, Cores: " << CoreCount << endl;
    }

    // Method to save data to file
    void Save(ofstream& out) const {
        out << Model << " " << MemoryMB << " " << CoreCount << endl;
    }

    // Method to load data from file
    void Load(ifstream& in) {
        in >> Model >> MemoryMB >> CoreCount;
    }
};

// CPU Characteristics class
class CpuCharacteristics {
public:
    string Model;
    int CoreCount;
    double FrequencyGHz;

    CpuCharacteristics() : Model("Unknown"), CoreCount(0), FrequencyGHz(0.0) {}

    CpuCharacteristics(const string& model, int coreCount, double frequencyGHz)
        : Model(model), CoreCount(coreCount), FrequencyGHz(frequencyGHz) {}

    void Display() const {
        cout << "CPU Model: " << Model << ", Cores: " << CoreCount << ", Frequency: " << FrequencyGHz << " GHz" << endl;
    }

    void Save(ofstream& out) const {
        out << Model << " " << CoreCount << " " << FrequencyGHz << endl;
    }

    void Load(ifstream& in) {
        in >> Model >> CoreCount >> FrequencyGHz;
    }
};

// RAM Characteristics class
class RamCharacteristics {
public:
    string Type;
    int SizeMB;

    RamCharacteristics() : Type("Unknown"), SizeMB(0) {}

    RamCharacteristics(const string& type, int sizeMB)
        : Type(type), SizeMB(sizeMB) {}

    void Display() const {
        cout << "RAM Type: " << Type << ", Size: " << SizeMB << " MB" << endl;
    }

    void Save(ofstream& out) const {
        out << Type << " " << SizeMB << endl;
    }

    void Load(ifstream& in) {
        in >> Type >> SizeMB;
    }
};

// LAN Characteristics class
class LanCharacteristics {
public:
    string AdapterModel;
    double SpeedGbps;

    LanCharacteristics() : AdapterModel("Unknown"), SpeedGbps(0.0) {}

    LanCharacteristics(const string& adapterModel, double speedGbps)
        : AdapterModel(adapterModel), SpeedGbps(speedGbps) {}

    void Display() const {
        cout << "LAN Adapter: " << AdapterModel << ", Speed: " << SpeedGbps << " Gbps" << endl;
    }

    void Save(ofstream& out) const {
        out << AdapterModel << " " << SpeedGbps << endl;
    }

    void Load(ifstream& in) {
        in >> AdapterModel >> SpeedGbps;
    }
};

// Cluster Node class
class ClusterNode {
public:
    CpuCharacteristics Cpu;
    GpuCharacteristics Gpu;
    RamCharacteristics Ram;
    LanCharacteristics Lan;

    void Display() const {
        cout << "Cluster Node Details:" << endl;
        Cpu.Display();
        Gpu.Display();
        Ram.Display();
        Lan.Display();
    }

    void Save(ofstream& out) const {
        Cpu.Save(out);
        Gpu.Save(out);
        Ram.Save(out);
        Lan.Save(out);
    }

    void Load(ifstream& in) {
        Cpu.Load(in);
        Gpu.Load(in);
        Ram.Load(in);
        Lan.Load(in);
    }
};

// Cluster class
class Cluster {
public:
    vector<ClusterNode> Nodes;

    void AddNode(const ClusterNode& node) {
        Nodes.push_back(node);
    }

    void Display() const {
        cout << "Cluster contains " << Nodes.size() << " nodes:" << endl;
        for (const auto& node : Nodes) {
            node.Display();
            cout << endl;
        }
    }

    void Save(const string& filename) const {
        ofstream out(filename);
        if (out.is_open()) {
            out << Nodes.size() << endl;
            for (const auto& node : Nodes) {
                node.Save(out);
            }
            out.close();
        } else {
            cerr << "Error opening file for saving!" << endl;
        }
    }

    void Load(const string& filename) {
        ifstream in(filename);
        if (in.is_open()) {
            size_t nodeCount;
            if (!(in >> nodeCount)) {
                cerr << "Error reading node count" << endl;
                return;
            }

            Nodes.resize(nodeCount);
            for (auto& node : Nodes) {
                node.Load(in);
                if (in.fail()) {
                    cerr << "Error reading node data" << endl;
                    break;
                }
            }
            in.close();
        } else {
            cerr << "Error opening file for loading!" << endl;
        }
    }
};

int main() {
    Cluster cluster;

    // Create a cluster node
    ClusterNode node1;
    node1.Cpu = {"Intel-Xeon", 8, 3.6};
    node1.Gpu = {"NvidiaRtx3080", 10240, 8704};
    node1.Ram = {"DDR4", 32768};
    node1.Lan = {"IntelEthernet", 1.0};

    // Add node to cluster
    cluster.AddNode(node1);

    // Display cluster data
    cluster.Display();

    // Save cluster data to file
    cluster.Save("cluster_data.txt");

    // Load cluster data from file
    Cluster importedCluster;
    importedCluster.Load("cluster_data.txt");

    cout << "Loaded cluster data from file:\n";
    importedCluster.Display();

    return 0;
}
