#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <optional>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <cctype>

using namespace std;
using uint64 = unsigned long long;

const uint64 MOD = 10000000000007ULL;

// Helper: Trim whitespace
string trim(const string &s) {
    auto start = s.begin();
    while (start != s.end() && isspace(*start)) start++;
    auto end = s.end();
    while (end != start && isspace(*(end-1))) end--;
    return string(start, end);
}

// Helper: Convert to uppercase
void to_upper(string &s) {
    for (auto &c : s) c = toupper(c);
}

string reverse_dna(const string &dna) {
    string rev;
    rev.reserve(dna.size());
    for (auto it = dna.rbegin(); it != dna.rend(); ++it) {
        char c = *it;
        switch(c) {
            case 'A': rev.push_back('T'); break;
            case 'T': rev.push_back('A'); break;
            case 'C': rev.push_back('G'); break;
            case 'G': rev.push_back('C'); break;
            default:
                throw runtime_error("Invalid DNA character: '" + string(1, c) + "'");
        }
    }
    return rev;
}

uint64 dna_to_num(char dna) {
    switch(dna) {
        case 'A': return 1;
        case 'T': return 2;
        case 'C': return 3;
        case 'G': return 4;
        default:
            throw runtime_error("Invalid DNA character: '" + string(1, dna) + "'");
    }
}

struct RefSeq {
    uint64 start;
    uint64 end;
    bool reverse;
};

struct Trace {
    RefSeq ref_seq;
    uint64 next;
    uint64 query_start;
    uint64 query_end;
};

void build_reference_hash(const string &dna, unordered_map<uint64, RefSeq> &map, bool reverse) {
    const size_t dna_len = dna.size();
    const string seq = reverse ? reverse_dna(dna) : dna;
    
    for (size_t start = 0; start < dna_len; ++start) {
        uint64 hash = 0;
        for (size_t end = start; end < dna_len; ++end) {
            hash = (hash * 5 + dna_to_num(seq[end])) % MOD;
            if (map.find(hash) == map.end()) {
                RefSeq ref;
                if (reverse) {
                    ref.start = dna_len - end - 1;
                    ref.end = dna_len - start - 1;
                } else {
                    ref.start = start;
                    ref.end = end;
                }
                ref.reverse = reverse;
                map[hash] = ref;
            }
        }
    }
}

vector<optional<Trace>> find_optimal_path(const string &query, const unordered_map<uint64, RefSeq> &ref_map) {
    const size_t query_len = query.size();
    vector<uint64> dp(query_len + 1, numeric_limits<uint64_t>::max() - 20);
    dp[query_len] = 0;
    vector<optional<Trace>> trace(query_len + 1, nullopt);

    for (int start = query_len - 1; start >= 0; --start) {
        uint64 hash = 0;
        for (size_t end = start; end < query_len; ++end) {
            hash = (hash * 5 + dna_to_num(query[end])) % MOD;
            if (const auto it = ref_map.find(hash); it != ref_map.end()) {
                const uint64 new_cost = dp[end + 1] + 1;
                if (new_cost < dp[start] || (new_cost == dp[start] && !it->second.reverse)) {
                    dp[start] = new_cost;
                    trace[start] = Trace{it->second, static_cast<uint64>(end + 1), 
                                       static_cast<uint64>(start), static_cast<uint64>(end)};
                }
            }
        }
    }
    return trace;
}

struct MatchSegment {
    RefSeq ref_info;
    uint64 query_start;
    uint64 query_end;
};

vector<MatchSegment> reconstruct_path(const vector<optional<Trace>> &trace, size_t query_len) {
    vector<MatchSegment> result;
    size_t pos = 0;
    while (pos < query_len) {
        if (!trace[pos].has_value()) {
            throw runtime_error("Alignment break: No match found at position " + to_string(pos));
        }
        const Trace &t = trace[pos].value();
        result.push_back({t.ref_seq, t.query_start, t.query_end});
        pos = t.next;
    }
    return result;
}

void validate_dna(const string &dna, const string &name) {
    for (char c : dna) {
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
            string msg = name + " contains invalid character: '" + string(1, c) + "'. Only A/T/C/G allowed";
            if (islower(c)) msg += " (detected lowercase, auto-converted to uppercase)";
            throw runtime_error(msg);
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // UI Initialization
    cout << "\033[1;34m\n======== DNA Sequence Alignment Tool ========\033[0m\n";
    
    // Reference sequence input
    cout << "\n\033[1;32m>>> Step 1/2: Enter Reference Sequence (long)\033[0m\n";
    cout << "Enter reference sequence (A/T/C/G only): \033[36m" << flush;
    string ref_seq;
    if (!getline(cin, ref_seq)) {
        cerr << "\033[31m\nError: Failed to read input\033[0m\n";
        return 1;
    }
    ref_seq = trim(ref_seq);
    to_upper(ref_seq);
    if (ref_seq.empty()) {
        cerr << "\033[31m\nError: Reference sequence cannot be empty\033[0m\n";
        return 1;
    }
    cout << "\033[0m"; // Reset color

    // Query sequence input
    cout << "\n\033[1;32m>>> Step 2/2: Enter Query Sequence (short)\033[0m\n";
    cout << "Enter query sequence (A/T/C/G only): \033[36m" << flush;
    string query_seq;
    if (!getline(cin, query_seq)) {
        cerr << "\033[31m\nError: Failed to read input\033[0m\n";
        return 1;
    }
    query_seq = trim(query_seq);
    to_upper(query_seq);
    if (query_seq.empty()) {
        cerr << "\033[31m\nError: Query sequence cannot be empty\033[0m\n";
        return 1;
    }
    cout << "\033[0m"; // Reset color

    try {
        // Validation
        validate_dna(ref_seq, "Reference sequence");
        validate_dna(query_seq, "Query sequence");

        // Build hash map
        unordered_map<uint64, RefSeq> ref_map;
        build_reference_hash(ref_seq, ref_map, false);
        build_reference_hash(ref_seq, ref_map, true);

        // Find optimal path
        auto trace = find_optimal_path(query_seq, ref_map);
        auto result = reconstruct_path(trace, query_seq.size());

        // Output results
        cout << "\n\033[1;34m======== Alignment Results ========\033[0m\n";
        cout << "Reference length: \033[33m" << ref_seq.size() << " bp\033[0m\n";
        cout << "Query length: \033[33m" << query_seq.size() << " bp\033[0m\n";
        cout << "\033[1;36mMatched segments: " << result.size() << "\033[0m\n\n";
        
        for (size_t i = 0; i < result.size(); ++i) {
            const auto& seg = result[i];
            const string seq = ref_seq.substr(seg.ref_info.start, 
                                            seg.ref_info.end - seg.ref_info.start + 1);
            cout << "\033[1;95mSegment " << i + 1 << ":\033[0m\n"
                 << "  \033[90mRef position:\033[0m [\033[35m" << seg.ref_info.start 
                 << "\033[0m-\033[35m" << seg.ref_info.end << "\033[0m]\n"
                 << "  \033[90mQuery position:\033[0m [\033[35m" << seg.query_start 
                 << "\033[0m-\033[35m" << seg.query_end << "\033[0m]\n"
                 << "  \033[90mStrand:\033[0m " 
                 << (seg.ref_info.reverse ? "\033[33mReverse complement\033[0m" : "\033[33mForward\033[0m") << "\n"
                 << "  \033[90mMatched sequence:\033[0m \033[36m" << seq << "\033[0m\n"
                 << "  \033[90mLength:\033[0m \033[32m" << seq.size() << " bp\033[0m\n\n";
        }
        cout << "\033[1;34m==========================\033[0m\n";

    } catch (const exception &e) {
        cerr << "\n\033[31mError: " << e.what() << "\033[0m\n";
        return 1;
    }

    return 0;
}
