#include <iostream>
#include <algorithm>
#include <filesystem> 
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <sstream>

using namespace std;
int latencyConstrain;

class Node {
public:
    int id;
    int duration = 0;
    int asap = 0;
    int alap = 0;
    int timeFrame = 0;
    int scheduledTime = -10000;
    float force;
    string operation;
    vector<float> pdf;
    vector<Node*> predecessors;
    vector<Node*> successors;
    
    Node(const int& id, const string&operation):id(id),operation(operation) {
       
        if (operation == "+") {
            this->duration = 1;
        }
        else if (operation == "*") {
            this->duration = 3;
        }
    }
};
class Graph {
public:
    // Variables
    vector<Node> vertices;
    vector<vector<Node*>> adjacencyList; // adjacency list


    // init
    Graph(vector<Node> nodes):vertices(nodes) {
        adjacencyList.resize(this->vertices.size());
    }

    // 將edge轉化為adjacency list
    void addEdge(vector<vector<int>> edges) {
        for (vector<int> edge : edges) {
            this->adjacencyList[edge[0]].push_back(&this->vertices[edge[1]]);
        }
    }
};
Graph create_graph(ifstream& file) {
    string line;
    // get Latency constrain
    while (getline(file, line)){
        if (line.find("Latency constrain") != std::string::npos) {
            size_t pos = line.find(": ");
            latencyConstrain = stoi(line.substr(pos + 1));
            cout << "Latency constrain: " << latencyConstrain << endl;
            getline(file, line);
            getline(file, line);
            break;
        }
    }
    vector<Node> nodes = {Node(0, "null")}; // intit
    vector<vector<int>> edges;
    // 逐行讀取testcase
    while(getline(file, line)){
        // 初始化變數
        std::istringstream iss(line);
        int id;
        string operation;
        vector<int> children;
        int child;
        iss >> id >> operation;
        while (iss >> child) {
            children.push_back(child);
        }

        // create Node
        nodes.push_back(Node(id, operation));

        // 把此Node的successors都加進去edges
        for (int child : children) {
            edges.push_back({id, child});
        }
        
    }
    Graph graph = Graph(nodes);
    // Node建立完成，把edge加到adjacency list
    graph.addEdge(edges);

    return graph;
};
vector<Node*> calculateASAP(Graph& graph) {
    // 回傳topological order
    // 建立一個Vector來儲存每個Node的Indegree
    vector<int> in_degree (size(graph.vertices), 0);
    vector<Node*> topologicalOrder;

    for (int i = 1; i < size(graph.vertices); i++) {
        for (Node* node : graph.adjacencyList[i]) {
            in_degree[node->id]++;
            node->predecessors.push_back(&graph.vertices[i]);
        }
    }

    // Kahn's algorithm求topological order
    queue<Node*> q;
    for (int i = 1; i < size(in_degree); i++) {
        if (in_degree[i] == 0) {
            graph.vertices[i].asap = 0;
            q.push(&graph.vertices[i]);
        }
    }

    while (!q.empty()) {
        Node* u = q.front();
        topologicalOrder.push_back(u);
        q.pop();
        for (Node* node : graph.adjacencyList[u->id]) { // 對所有u的succussors
            node->asap = max(node->asap, u->asap + u->duration); // 更新asap
            if (--in_degree[node->id] == 0) {
                q.push(node);//將indegree=0的加入queue
            }
        }
    }
    return topologicalOrder;
};
void calculateALAP(Graph& graph) {
    // 儲存所有Node的Outdegree
    vector<int> out_degree (size(graph.vertices), 0);

    // 根據每個Node的child數量增加Outdegree
    for (Node& node : graph.vertices) {
        node.alap = latencyConstrain;
        out_degree[node.id] = size(graph.adjacencyList[node.id]);
        for (Node* succ : graph.adjacencyList[node.id]) {
            node.successors.push_back(succ);
        }
    }

    // Kahn's algorithm求reversed topological order
    queue<Node*> q;
    for (int i = 1; i <= size(out_degree); i++) {
        if (out_degree[i] == 0) {
            q.push(&graph.vertices[i]);
        }
    }

    while (!q.empty()) {
        Node* u = (q.front());
        q.pop();

        for (Node& node : graph.vertices) {
            for (Node* succ : graph.adjacencyList[node.id]) {
                if (succ == u) {
                    node.alap = min(node.alap, u->alap - u->duration);
                    
                    if (--out_degree[node.id] == 0) {
                        q.push(&node);
                    }
                }
            }
        }
        
    }
    //記得要-原本的執行時間
    for (Node& node : graph.vertices) {
        node.alap = node.alap - node.duration;
    }
}
float calculateForce(Node& node, int time, vector<vector<float>> resource_loading) {
    // Self-force 公式:（1-node.pdf[time])* resource_loading[type][time] + sum(resource在時間i的loading * (0-node.pdf[time]))
    int multi = (node.operation == "*" ? 1 : 0);
    float self_force = resource_loading[multi][time] * (1-node.pdf[time]);
    for (int i = 0; i < latencyConstrain; i++) {
        self_force -= resource_loading[multi][i] * node.pdf[i];
    }
    return self_force;
}
namespace fs = std::filesystem;
void doTestcase(string filename){
    ifstream infile(filename);
    if (!infile.is_open()) {
        cout << "Can not open the file.";
        exit(1);
    }

    Graph graph = create_graph(infile);
    vector<Node*> topologicalOrder = calculateASAP(graph);
    calculateALAP(graph);

    // 計算每個Node的timeFrame
    for (Node& node : graph.vertices) {
        node.timeFrame = node.alap - node.asap + 1;
    }
    
    // 計算每個Node的PDF
    for (Node& node : graph.vertices) {
        node.pdf.resize(latencyConstrain + 1);
        float probability = (double)1.0 / node.timeFrame;
        for (int t = node.asap; t <= node.alap + node.duration - 1 && t < latencyConstrain; t++) {
            node.pdf[t] = probability;
        }
    }


    // 計算每種Resource的Loading distribution
    vector<vector<float>> loading_distribution(2, vector<float>(latencyConstrain + 1));
    for (Node& node : graph.vertices) {
        if (node.operation == "+") {
            for (int i = 0; i < size(node.pdf); i++){
                loading_distribution[0][i] += node.pdf[i];
            }
        }
        else if (node.operation == "*") {
            for (int i = 0; i < size(node.pdf); i++){
                loading_distribution[1][i] += node.pdf[i];
            }
        }
    }
    for (Node& node : graph.vertices) {
        if (node.operation != "+" && node.operation != "*") {
            continue;
        }
        // cout << "Operation " << node.id << " - ASAP: " << node.asap << ", ALAP: " << node.alap << ", pdf: ";
        for (auto& p : node.pdf) {
            // cout << p << ", ";
        }
        // cout << endl;
    }
    for (vector<float> resource : loading_distribution) {
        for (float loading : resource) {
            // cout << loading << ", ";
        }
        // cout << endl;
    }

    // 計算每個時間用的resources
    vector<vector<int>> resources(2, vector<int> (latencyConstrain + 1, 0));

    // 照拓樸排序來schedule就不用計算predecessor force

    for (Node* node : topologicalOrder) {
        // 只計算*和+
        if (node->operation != "*" && node->operation != "+"){
            continue;
        }
        vector<float> node_forces (latencyConstrain, 100000.0); // 用來儲存此node在timeFrame內每個time的force
        for (int j = node->asap; j < node->alap + 1; j++){
            // 計算此node在時間j的force
            node_forces[j] = calculateForce(*node, j, loading_distribution);
            
            // 如果時間j + 此operation的duration大於後續operation的asap，就要計算其successor force
            for (Node* succ : graph.adjacencyList[node->id]) {
                if (j + node->duration > succ->asap) {
                    float succ_force = 0;
                    int type = (succ->operation == "+") ? 0 : 1;
                    for (int k = j + node->duration; k < succ->alap; k++) {
                        succ_force += loading_distribution[type][k]; // 縮減後的resource loading * pdf(k)
                    }
                    succ_force /= succ->alap - (j + node->duration) + 1;
                    for (int k = succ->asap; k < succ->alap; k++) {
                        succ_force -= loading_distribution[type][k] / (succ->alap - succ->asap); // 縮減前的resource loading * pdf(k)
                    }
                    // cout << "Operation " << succ->id << "'s P/S force to operation " << node->id << " is " << succ_force << endl;
                    node_forces[j] += succ_force;
                }
            }
            // cout << "Operation " << node->id << "'s force at time " << j << " is " << node_forces[j] << ", loading_distribution at time " << j << " is " << loading_distribution[0][j]<< endl;
        }
        // 將node在timeFrame中所有時間的force算出來後，分派該node到最小force的時間
        auto minForce = min_element(node_forces.begin(), node_forces.end());
        int time = distance(node_forces.begin(), minForce);
        // cout << "Assign operation " << node->id << " to cycle " << min_index << "." << endl;
        node->scheduledTime = time;
        if (node->operation == "+") {
            resources[0][time]++;
        }
        else {
            if (node->operation == "*") {
                resources[1][time]++;
                resources[1][time + 1]++;
                resources[1][time + 2]++;
            }
        }

        fill(node->pdf.begin(), node->pdf.end(), 0);
        node->pdf[time] = 1;

        // 更動sucessors的timeFrame與PDF
        for (Node* succ : graph.adjacencyList[node->id]) {
            succ->asap = max(succ->asap, node->scheduledTime + node->duration);
            succ->timeFrame = succ->alap - succ->asap + 1;
            float prob = (double)1.0/succ->timeFrame;
            for (int i = succ->asap; i <= succ->alap + succ->duration - 1; i++) {
                succ->pdf[i] = prob;
            }
        }

        // 重新計算Resource distribution
        for (vector<float> type : loading_distribution) {
            type.clear();
        }
        for (Node& node : graph.vertices) {
            if (node.operation == "+") {
                for (int i = 0; i < size(node.pdf); i++){
                    loading_distribution[0][i] += node.pdf[i];
                }
            }
            else if (node.operation == "*") {
                for (int i = 0; i < size(node.pdf); i++){
                    loading_distribution[1][i] += node.pdf[i];
                }
            }
        }
    }
    vector<vector<int>> scheduled;
    for (int i = 0; i < latencyConstrain; i++) {
        vector<int> temp;
        for (Node& node : graph.vertices) {
            if (node.scheduledTime == i || (node.scheduledTime < i && node.scheduledTime + node.duration - 1 >= i) ) {
                // cout<<"Time: "<<i<<" Node iD:"<<node.id <<" ScheduledTime: "<<node.scheduledTime<<endl;
                temp.push_back(node.id);
            }
        }
        // cout<<"next iteration"<<endl;
        scheduled.push_back(temp);
    }

    auto adder = max_element(resources[0].begin(), resources[0].end());
    auto multiplier = max_element(resources[1].begin(), resources[1].end());
    ofstream outfile(filename + ".out");
    outfile << *adder << endl;
    outfile << *multiplier << endl;
    for (vector<int> cycle : scheduled) {
        for (int i = 0; i < size(cycle); i++) {
            outfile << cycle[i];
            if (i != size(cycle) - 1) {
                outfile << " ";
            }
        }
        outfile << endl;
    }
}
int main() {
    string directory_path = fs::current_path(); // 將此處目前的path
    bool found_testcase = false;
    for (const auto& entry : fs::directory_iterator(directory_path)) {
        if (entry.is_regular_file() && entry.path().filename().string().substr(0, 8)==("testcase")) {
            found_testcase = true;
            cout<<"Find: "<<entry.path().filename().string()<<endl;
            doTestcase(entry.path().filename().string());
            cout<<"Finish: "<<entry.path().filename().string()<<endl;
        }
    }
    return 0;
}