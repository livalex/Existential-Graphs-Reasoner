// Copyright 2019 Luca Istrate, Danut Matei
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <cassert>
#include "./aegraph.h"

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    // Assert checks if the expression in the brackets is real.
    // If not,it returns a fail error.
    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; i++) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }

    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; i++) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

// nu mergem pe atomi
std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
    // 10p
    std::vector<std::vector<int>> road;
    int len_subgraphs = num_subgraphs();
    for (int i = 0; i < len_subgraphs; i++) {
            if (subgraphs[i].num_subgraphs() == 1 &&
            subgraphs[i].num_atoms() == 0) {
                road.push_back({i});
            }
            auto r = subgraphs[i].possible_double_cuts();
            for (auto& v : r) {
                v.insert(v.begin(), i);
            }
            copy(r.begin(), r.end(), back_inserter(road));
    }
    return road;
}

void AEGraph::double_cut_helper(std::vector<int> where, AEGraph &node) const {
    if (where.capacity() != 1) {
        unsigned int index;
        index = where[0];
        where.erase(where.begin());
        node.double_cut_helper(where, node.subgraphs[index]);
    } else if (where.capacity() == 1) {
        AEGraph aux = node.subgraphs[where[0]].subgraphs[0];
        for (int i = 0; i < aux.num_subgraphs(); ++i) {
            node.subgraphs.push_back(aux.subgraphs[i]);
        }
        for (int i = 0; i < aux.num_atoms(); ++i) {
            node.atoms.push_back(aux.atoms[i]);
        }
        node.subgraphs.erase(node.subgraphs.begin() + where[0]);
    }
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
    // 10p
    AEGraph auxiliar = *this;
    double_cut_helper(where, auxiliar);
    return auxiliar;
}

// mergem pe atomi
std::vector<std::vector<int>> AEGraph::possible_erasures(int level) const {
    // 10p
    std::vector<std::vector<int>> road;
    int len_subgraphs = num_subgraphs();
    int len_atoms = num_atoms();
    for (int i = 0; i < (len_subgraphs + len_atoms); i++) {
        if (level % 2 != 0) {
            road.push_back({i});
            if (level != -1 && (!num_subgraphs() &&
            num_atoms() == 1) && road.size()) {
                road.pop_back();
            }
            if (level != -1 && (num_subgraphs() == 1 &&
            !num_atoms()) && road.size()) {
                road.pop_back();                }
            }
         	if (i < num_subgraphs()) {
            	auto r = subgraphs[i].possible_erasures(level + 1);
                for (auto& v : r) {
                    v.insert(v.begin(), i);
                }
                copy(r.begin(), r.end(), back_inserter(road));
            }
    }
    return road;
}

void AEGraph::erase_helper(std::vector<int> where, AEGraph &node) const {
    if (where.capacity() != 1) {
        unsigned int index;
        index = where[0];
        where.erase(where.begin());
        node.erase_helper(where, node.subgraphs[index]);
    } else if (where.capacity() == 1) {
        if (node.num_subgraphs() - where[0] <= 0) {
            node.atoms.erase(node.atoms.begin() -
            node.num_subgraphs() + where[0]);
        } else if (node.num_subgraphs() - where[0] > 0) {
            node.subgraphs.erase(node.subgraphs.begin()
            + where[0]);
        }
    }
}

AEGraph AEGraph::erase(std::vector<int> where) const {
    // 10p
    AEGraph auxiliar = *this;
    erase_helper(where, auxiliar);
    return auxiliar;
}


std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
    // 20p
    std::vector<std::vector<int>> road;
    for (int i = 0; i < num_subgraphs(); ++i) {
        for (int j = 0; j < num_subgraphs(); ++j) {
            if (i != j) {
                auto r = subgraphs[j].get_paths_to(subgraphs[i]);
                for (auto& v : r) {
                    v.insert(v.begin(), j);
                }
                copy(r.begin(), r.end(), back_inserter(road));
            }
        }
    }
    for (int i = 0; i < num_atoms(); ++i) {
        for (int j = 0; j < num_subgraphs(); ++j) {
            auto r = subgraphs[j].get_paths_to(atoms[i]);
            for (auto& v : r) {
                v.insert(v.begin(), j);
            }
            copy(r.begin(), r.end(), back_inserter(road));
        }
    }
    return road;
}

void AEGraph::deiterate_helper(std::vector<int> where, AEGraph &node) const {
    if (where.capacity() != 1) {
        unsigned int index;
        index = where[0];
        where.erase(where.begin());
        node.deiterate_helper(where, node.subgraphs[index]);
    } else if (where.capacity() == 1) {
        if (node.num_subgraphs() - where[0] <= 0) {
            node.atoms.erase(node.atoms.begin() -
            node.num_subgraphs() + where[0]);
        } else if (node.num_subgraphs() - where[0] > 0) {
            node.subgraphs.erase(node.subgraphs.begin() + where[0]);
        }
    }
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
    // 10p
    AEGraph auxiliar = *this;
    deiterate_helper(where, auxiliar);
    return auxiliar;
}

