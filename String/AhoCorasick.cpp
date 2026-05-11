/*
 * Authors:
 *  - Mohsen Dehbashi
 *  - Ghazal Rabiei Faradonbeh
 *
 * Aho-Corasick multiple pattern matching algorithm
 *
 * Finds all occurrences of multiple patterns in a text in O(|text| + total patterns length + number of matches) time
 */

#include <vector>
#include <string>
#include <queue>
#include <unordered_map>

class AhoCorasick {
    struct Node {
        std::unordered_map<char, Node*> next;
        Node* link = nullptr;
        std::vector<size_t> output;
    };

    Node* root = new Node();

    /*
     * Inserts a pattern to the trie
     */
    void insert(const std::string& pattern, const size_t index) const {
        Node* node = root;
        for (char c : pattern) {
            if (!node->next.count(c))
                node->next[c] = new Node();
            node = node->next[c];
        }
        node->output.push_back(index);
    }

    /*
     * Builds failure links using BFS
     */
    void build() const {
        std::queue<Node*> q;
        root->link = root;

        // Initialize root children
        for (auto& [c, node] : root->next) {
            node->link = root;
            q.push(node);
        }

        while (!q.empty()) {
            Node* u = q.front();
            q.pop();

            for (auto& [c, v] : u->next) {
                Node* fallback = u->link;

                // Find best fallback
                while (fallback != root && !fallback->next.count(c)) {
                    fallback = fallback->link;
                }

                if (fallback->next.count(c) && fallback->next[c] != v)
                    v->link = fallback->next[c];
                else
                    v->link = root;

                // Merge outputs
                for (const size_t x : v->link->output)
                    v->output.push_back(x);

                q.push(v);
            }
        }
    }

public:
    /*
     * Finds all occurrences of patterns in text
     *
     * @param patterns list of patterns
     * @param text The text to search in
     * @return result[i] = list of starting indices of pattern i
     *
     * Time complexity: O(|text| + total patterns length + number of matches)
     * Space complexity: O(|text| + total patterns length + number of matches)
     */
    std::vector<std::vector<size_t>> findAll(
        const std::vector<std::string>& patterns,
        const std::string& text
    ) const {
        // Build trie
        for (size_t i = 0; i < patterns.size(); i++) {
            insert(patterns[i], i);
        }

        // Build links
        build();

        std::vector<std::vector<size_t>> result(patterns.size());

        Node* node = root;

        // Traverse on text and trie simultaneously
        for (size_t i = 0; i < text.size(); i++) {
            char c = text[i];

            while (node != root && !node->next.count(c)) {
                node = node->link;
            }

            if (node->next.count(c))
                node = node->next[c];

            for (const size_t patternIndex : node->output) {
                const size_t len = patterns[patternIndex].size();
                result[patternIndex].push_back(i - len + 1);
            }
        }

        return result;
    }
};