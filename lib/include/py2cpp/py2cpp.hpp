#ifndef _HOME_UBUNTU_GITHUB_PY2CPP_PY2CPP_HPP
#define _HOME_UBUNTU_GITHUB_PY2CPP_PY2CPP_HPP 1

#include <unordered_map>
#include <unordered_set>
#include <initializer_list>
#include <utility>

namespace py {

template <typename Key>
struct set : std::unordered_set<Key> {
    using _Self = set<Key>;

    explicit set() : 
        std::unordered_set<Key>{} {}

    explicit set(std::initializer_list<Key> init) : 
        std::unordered_set<Key>{init} {}

    bool contains(const Key& key) const {
        return this->find(key) != this->end();
    }

    _Self copy() const { return *this; }

    _Self& operator=(const _Self& ) = delete;

    _Self& operator=(_Self&& ) = delete;

private:
    set(const _Self& ) = default;
};

template <typename Key>
inline bool operator<(const Key& key, const set<Key>& m) {
    return m.contains(key);
}

template <typename Key>
inline size_t len(const set<Key>& m) {
    return m.size();
}

template <typename Key>
set(std::initializer_list<Key> ) -> set<Key>;

// template <typename Key>
// set(std::initializer_list<const char*> ) -> set<std::string>;


template <typename Key, typename T>
struct dict : std::unordered_map<Key, T> {
    using value_type = std::pair<const Key, T>;
    using _Self = dict<Key, T>;

    explicit dict() : 
        std::unordered_map<Key, T>{} {}

    explicit dict(std::initializer_list<value_type> init) : 
        std::unordered_map<Key, T>{init} {}

    bool contains(const Key& key) const {
        return this->find(key) != this->end();
    }

    T get(const Key& key, const T& default_value) {
        if (!contains(key)) return default_value;
        return (*this)[key];
    }

    _Self copy() const { return *this; }

    _Self& operator=(const _Self& ) = delete;

    _Self& operator=(_Self&& ) = delete;

private:
    dict(const _Self& ) = default;
};

template <typename Key, typename T>
inline bool operator<(const Key& key, const dict<Key, T>& m) {
    return m.contains(key);
}


template <typename Key, typename T>
inline size_t len(const dict<Key, T>& m) {
    return m.size();
}

template <typename Key, typename T>
dict(std::initializer_list<std::pair<const Key, T>> ) -> dict<Key, T>;


}

#endif
