#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <random>
#include <stack>
#include <filesystem>
#include <chrono>

#include <cstdlib> // for rand
#include <ctime>

using namespace std;
using namespace std::filesystem;

const long double PI = 3.14159265359, INF = 1e8;
int W;
const unsigned int amount_of_mutations_in_percent = 20;
const unsigned int amount_of_mutations_of_genom_percent = 10;

double global_ratio;

std::chrono::steady_clock::time_point begin0;
const int time_duration = 360;


struct Genom
{
    vector <bool> B;
    int weight;
    int value;

    Genom () : weight(0), value(0) {};
    Genom (vector <pair <int, int> > &a, unsigned int n, unsigned int i)
    {
        B.resize(n, 0);
        B[i] = 1;
        weight = a[i].second;
        value = a[i].first;
    }

    bool operator [](unsigned int i)
    {
        return B[i];
    }

    unsigned int size()
    {
        return B.size();
    }

    double get_ratio()
    {
        double ans = value;
        ans /= weight;
        return ans;
    }

    Genom* sum(vector <pair <int, int> > &a, Genom* g)
    {
        Genom* new_Genom = new Genom();
        new_Genom->B.resize(g->size());
        for (int i = 0; i < g->size(); i++)
        {
            int flag = (B[i] > 0 || g->B[i] > 0);
            new_Genom->B[i] = flag;
            new_Genom->weight += a[i].second * flag;
            new_Genom->value += a[i].first * flag;
        }
        return new_Genom;
    }

    void mutate(unsigned int n) // n = amount_of_mutations_of_genom
    {
        unsigned int k = size();
        for (int i = 0; i < n; i++)
        {
            int j = abs(rand()) % k;
            B[j] = !B[j];
        }
    }
};

bool cmp(Genom* a, Genom* b)
{
//    return a->get_ratio() > b->get_ratio();
    return a->value > b->value;
}

struct Population
{
    vector <Genom*> G;

    Population (){};
    Population (vector <pair <int, int> > &a, unsigned int n, unsigned int max_P)
    {
        G.reserve(max_P);
        #pragma omp parallel for
        for (int i = 0; i < max_P; i++)
        {
            Genom* new_genom = new Genom(a, n, i%n);
            if (new_genom->weight <= W)
            {
                #pragma omp critical
                G.push_back(new_genom);
            }
            else
            {
                delete new_genom;
            }
        }
    }
    ~Population ()
    {
        for (int i = 0; i < G.size(); i++)
        {
            delete G[i];
        }
    }

    Genom* operator [](unsigned int i)
    {
        return G[i];
    }

    void add(Genom* a)
    {
        #pragma omp critical
        G.push_back(a);
    }

    unsigned int size()
    {
        return G.size();
    }

    void gen_sort()
    {
        sort(G.begin(), G.end(), cmp);
    }

    void selection(unsigned int max_P)
    {
        vector <Genom*> new_G(max_P);
        for (int i = 0; i < max_P; i++)
        {
            new_G[i] = G[i];
        }
        G = new_G;
    }

    void mutation(unsigned int n) // n = amount_of_mutations
    {
        unsigned int max_P = G.size();
        #pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            int j = abs(rand())%max_P;
            unsigned int amount_of_mutations_of_genom = amount_of_mutations_of_genom_percent * G[j]->size() / 100;
            //#pragma omp critical
            G[j]->mutate(amount_of_mutations_of_genom);
        }
    }
};

int epoch(vector <pair <int, int> > &a, Population &P, unsigned int max_P)
{
//    cout << "start" << endl;
    for (int i = 0; i < max_P; i++)
    {
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int j = i+1; j < max_P; j++)
        {
            Genom* new_Genom = P[i]->sum(a, P[j]);
            if (new_Genom->weight <= W)
            {
                P.add(new_Genom);
            }
            else
            {
                delete new_Genom;
            }

///--------------------------------------------------------------------------------------------------------------time
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            //cout << std::chrono::duration_cast<std::chrono::seconds>(end - begin0).count() << endl;
            if (std::chrono::duration_cast<std::chrono::seconds>(end - begin0).count() > time_duration) return -1;
///--------------------------------------------------------------------------------------------------------------time
        }
    }
//    for (int i = 0; i < P.size(); i++)
//    {
//        for (int j = 0; j < P[i]->size(); j++)
//        {
//            cout << (*P[i])[j] << " ";
//        }
//        cout << P[i]->value << " ' " << P[i]->weight << endl;
//    }

    P.gen_sort();

    P.selection(max_P);
    unsigned int amount_of_mutations = amount_of_mutations_in_percent * max_P / 100;
    P.mutation(amount_of_mutations);

    int value = 0;
    #pragma omp parallel for reduction(max:value)
    for (int i = 0; i < max_P; i++)
    {
        if (P[i]->value > value)
        {
            value = P[i]->value;
        }
    }
    return value;
}

int code(basic_ifstream<char, char_traits<char>> &in)
{
    std::srand(std::time({})); // use current time as seed for random generator
    const int random_value = std::rand();

//    vector <int> indexes(30);
//    for (int i = 0; i < 30; i++)
//    {
//        indexes[i] = abs(rand()%4);
//        cout << indexes[i] << endl;
//    }
//    return 0;

    W = 0;

    int n, w;
    vector <pair <int, int> > a;

    double max_V = 0, max_W = 0;

    if (in.is_open())
    {
        in >> n >> w;
        W = w;
//        cout << n << " " << w << endl;

        #pragma omp parallel for reduction(+:max_V, max_W)
        for (int i = 0; i < n; i++)
        {
            int vi, wi;
            #pragma omp critical
            {
                in >> vi >> wi;
            }
            #pragma omp critical
            {
                a.emplace_back(vi, wi);
            }
            //a.push_back(make_pair(vi, wi));
            max_V += vi;
            max_W += wi;
//            cout << vi << " " << wi << endl;
        }
//        cout << endl;
    }
//    cout << endl;
    in.close();

///--------------------------

    global_ratio = max_V / max_W;

    int prev_value = 0;
    int new_value = -1;
    int flag = 0;

    unsigned int max_P = 4*n;
    Population P(a, n, max_P);

//    int value = 0;
//    for (int i = 0; i < max_P; i++)
//    {
//        value = max(P[i]->value, value);
//    }
//    cout << value << endl;

    while(flag < 2)
    {
///--------------------------------------------------------------------------------------------------------------time
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        //cout << std::chrono::duration_cast<std::chrono::seconds>(end - begin0).count() << endl;
        if (std::chrono::duration_cast<std::chrono::seconds>(end - begin0).count() > time_duration) return new_value;
///--------------------------------------------------------------------------------------------------------------time
        new_value = epoch(a, P, max_P);
//        cout << new_value << endl;
//        for (int i = 0; i < max_P; i++)
//        {
//            for (int j = 0; j < P[i]->size(); j++)
//            {
//                cout << (*P[i])[j] << " ";
//            }
//            cout << endl;
//        }

        if (prev_value - new_value == 0) flag++;
        else flag = 0;

        prev_value = new_value;
    }

//    cout << new_value << endl;

    return new_value;
}

int main()
{
//    ifstream in("data/ks_4_0");

    vector <string> data;

    std::string path = "data/"; // ѕуть к папке
    for (auto& entry : directory_iterator(path))
    {
        std::string path_string{entry.path().u8string()};
        cout << path_string.size() << " ";
        //if (path_string.size() > 13) continue;
        data.push_back(path_string);
    }
    cout << endl;

    std::ofstream out("table.txt", std::ios::app);

    for (int i = 0; i < data.size(); i++)
    {
        cout << data[i] << " -- ";
        ifstream in(data[i]);
        begin0 = std::chrono::steady_clock::now();
        int value = code(in);
        cout << "Answer: " << value << endl;
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        out << data[i] << "   " << value << "       " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin0).count() << endl;
    }
}
