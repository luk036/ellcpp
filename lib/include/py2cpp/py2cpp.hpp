#ifndef _HOME_UBUNTU_GITHUB_PY2CPP_PY2CPP_HPP
#define _HOME_UBUNTU_GITHUB_PY2CPP_PY2CPP_HPP 1

#include <unordered_map>
#include <unordered_set>
#include <initializer_list>
#include <utility>

namespace py {

/**
 * @brief 
 * 
 * @tparam Key 
 */
template <typename Key>
struct set : std::unordered_set<Key> {
    using _Self = set<Key>;

    /**
     * @brief Construct a new set object
     * 
     */
    explicit set() : 
        std::unordered_set<Key>{} {}

    /**
     * @brief Construct a new set object
     * 
     * @param init 
     */
    explicit set(std::initializer_list<Key> init) : 
        std::unordered_set<Key>{init} {}

    /**
     * @brief 
     * 
     * @param key 
     * @return true 
     * @return false 
     */
    bool contains(const Key& key) const {
        return this->find(key) != this->end();
    }

    /**
     * @brief 
     * 
     * @return _Self 
     */
    _Self copy() const { return *this; }

    /**
     * @brief 
     * 
     * @return _Self& 
     */
    _Self& operator=(const _Self& ) = delete;

    /**
     * @brief 
     * 
     * @return _Self& 
     */
    _Self& operator=(_Self&& ) = delete;

private:
    /**
     * @brief Construct a new set object
     * 
     */
    set(const _Self& ) = default;
};

/**
 * @brief 
 * 
 * @tparam Key 
 * @param key 
 * @param m 
 * @return true 
 * @return false 
 */
template <typename Key>
inline bool operator<(const Key& key, const set<Key>& m) {
    return m.contains(key);
}

/**
 * @brief 
 * 
 * @tparam Key 
 * @param m 
 * @return size_t 
 */
template <typename Key>
inline size_t len(const set<Key>& m) {
    return m.size();
}

/**
 * @brief Template Deduction Guide
 * 
 * @tparam Key 
 */
template <typename Key>
set(std::initializer_list<Key> ) -> set<Key>;

// template <typename Key>
// set(std::initializer_list<const char*> ) -> set<std::string>;

/**
 * @brief 
 * 
 * @tparam Key 
 * @tparam T 
 */
template <typename Key, typename T>
struct dict : std::unordered_map<Key, T> {
    using value_type = std::pair<const Key, T>;
    using _Self = dict<Key, T>;

    /**
     * @brief Construct a new dict object
     * 
     */
    explicit dict() : 
        std::unordered_map<Key, T>{} {}

    /**
     * @brief Construct a new dict object
     * 
     * @param init 
     */
    explicit dict(std::initializer_list<value_type> init) : 
        std::unordered_map<Key, T>{init} {}

    /**
     * @brief 
     * 
     * @param key 
     * @return true 
     * @return false 
     */
    bool contains(const Key& key) const {
        return this->find(key) != this->end();
    }

    /**
     * @brief 
     * 
     * @param key 
     * @param default_value 
     * @return T 
     */
    T get(const Key& key, const T& default_value) {
        if (!contains(key)) return default_value;
        return (*this)[key];
    }

    /**
     * @brief 
     * 
     * @return _Self 
     */
    _Self copy() const { return *this; }

    /**
     * @brief 
     * 
     * @return _Self& 
     */
    _Self& operator=(const _Self& ) = delete;

    /**
     * @brief 
     * 
     * @return _Self& 
     */
    _Self& operator=(_Self&& ) = delete;

private:
    /**
     * @brief Construct a new dict object
     * 
     */
    dict(const _Self& ) = default;
};

/**
 * @brief 
 * 
 * @tparam Key 
 * @tparam T 
 * @param key 
 * @param m 
 * @return true 
 * @return false 
 */
template <typename Key, typename T>
inline bool operator<(const Key& key, const dict<Key, T>& m) {
    return m.contains(key);
}

/**
 * @brief 
 * 
 * @tparam Key 
 * @tparam T 
 * @param m 
 * @return size_t 
 */
template <typename Key, typename T>
inline size_t len(const dict<Key, T>& m) {
    return m.size();
}

/**
 * @brief Template Deduction Guide
 * 
 * @tparam Key 
 * @tparam T 
 */
template <typename Key, typename T>
dict(std::initializer_list<std::pair<const Key, T>> ) -> dict<Key, T>;


}

#endif
