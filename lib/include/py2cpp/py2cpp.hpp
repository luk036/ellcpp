#pragma once

#include <initializer_list>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

template <typename T>
using Value_type = typename T::value_type;

namespace py
{

/*!
 * @brief
 *
 * @tparam T
 * @tparam decltype(std::begin(std::declval<T>()))
 * @tparam decltype(std::end(std::declval<T>()))
 * @param[in] iterable
 * @return constexpr auto
 */
template <typename T, typename TIter = decltype(std::begin(std::declval<T>())),
    typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T&& iterable)
{
    struct iterator
    {
        size_t i;
        TIter iter;
        bool operator!=(const iterator& other) const
        {
            return iter != other.iter;
        }
        void operator++()
        {
            ++i;
            ++iter;
        }
        auto operator*() const
        {
            return std::tie(i, *iter);
        }
        auto operator*()
        {
            return std::tie(i, *iter);
        }
    };
    struct iterable_wrapper
    {
        T iterable;
        auto begin()
        {
            return iterator {0, std::begin(iterable)};
        }
        auto end()
        {
            return iterator {0, std::end(iterable)};
        }
    };
    return iterable_wrapper {std::forward<T>(iterable)};
}

// template <typename T>
// constexpr auto range(T stop) {
//     struct iterator {
//         T i;
//         constexpr bool operator!=(const iterator &other) const { return i !=
//         other.i; } constexpr bool operator==(const iterator &other) const {
//         return i == other.i; } constexpr T operator*() const { return i; }
//         constexpr iterator &operator++() {
//             ++i;
//             return *this;
//         }
//     };

//     struct iterable_wrapper {
//         using value_type = T; // luk
//         T stop;
//         constexpr auto begin() const { return iterator{0}; }
//         constexpr auto end() const { return iterator{stop}; }
//         constexpr auto empty() const -> bool { return stop == 0; }
//         constexpr auto size() const -> size_t { return stop; }
//         constexpr auto operator[](size_t n) const -> T { return n; } // no
//         bounds checking constexpr auto contains(T n) const -> bool { return n
//         < stop; }
//     };

//     if (stop < 0) stop = 0;
//     return iterable_wrapper{stop};
// }

template <typename T>
inline constexpr auto range(T start, T stop)
{
    struct __iterator
    {
        T i;
        constexpr bool operator!=(const __iterator& other) const
        {
            return this->i != other.i;
        }
        constexpr bool operator==(const __iterator& other) const
        {
            return this->i == other.i;
        }
        constexpr T operator*() const
        {
            return this->i;
        }
        constexpr __iterator& operator++()
        {
            ++this->i;
            return *this;
        }
    };

    struct iterable_wrapper
    {
      public:
        using value_type [[maybe_unused]] = T; // luk:
        using iterator = __iterator;           // luk
        T start;
        T stop;
        [[nodiscard]] constexpr auto begin() const
        {
            return iterator {this->start};
        }
        [[nodiscard]] constexpr auto end() const
        {
            return iterator {this->stop};
        }
        [[nodiscard]] constexpr auto empty() const -> bool
        {
            return this->stop == this->start;
        }
        [[nodiscard]] constexpr auto size() const -> size_t
        {
            return this->stop - this->start;
        }
        constexpr auto operator[](size_t n) const -> T
        {
            return T(this->start + n);
        } // no bounds checking
        [[nodiscard]] constexpr auto contains(T n) const -> bool
        {
            return !(n < this->start) && n < this->stop;
        }
    };

    stop = std::max(stop, start);
    // if (stop < start) {
    //     stop = start;
    // }
    return iterable_wrapper {start, stop};
}

template <typename T>
inline constexpr auto range(T stop)
{
    return range(T(0), stop);
}

/*!
 * @brief
 *
 * @tparam Key
 */
template <typename Key>
class set : public std::unordered_set<Key>
{
    using _Self = set<Key>;

  public:
    /*!
     * @brief Construct a new set object
     *
     */
    set()
        : std::unordered_set<Key> {}
    {
    }

    /*!
     * @brief Construct a new set object
     *
     */
    template <typename FwdIter>
    set(const FwdIter& start, const FwdIter& stop)
        : std::unordered_set<Key>(start, stop)
    {
    }

    /*!
     * @brief Construct a new set object
     *
     * @param[in] init
     */
    set(std::initializer_list<Key> init)
        : std::unordered_set<Key> {init}
    {
    }

    /*!
     * @brief
     *
     * @param[in] key
     * @return true
     * @return false
     */
    bool contains(const Key& key) const
    {
        return this->find(key) != this->end();
    }

    /*!
     * @brief
     *
     * @return _Self
     */
    _Self copy() const
    {
        return *this;
    }

    /*!
     * @brief
     *
     * @return _Self&
     */
    _Self& operator=(const _Self&) = delete;

    /*!
     * @brief
     *
     * @return _Self&
     */
    _Self& operator=(_Self&&) noexcept = default;

    /*!
     * @brief Move Constructor (default)
     *
     */
    set(set<Key>&&) noexcept = default;

  private:
    /*!
     * @brief Copy Constructor (deleted)
     *
     * Copy through explicitly the public copy() function!!!
     */
    set(const set<Key>&) = default;
};

/*!
 * @brief
 *
 * @tparam Key
 * @param[in] key
 * @param[in] m
 * @return true
 * @return false
 */
template <typename Key>
inline bool operator<(const Key& key, const set<Key>& m)
{
    return m.contains(key);
}

/*!
 * @brief
 *
 * @tparam Key
 * @param[in] m
 * @return size_t
 */
template <typename Key>
inline size_t len(const set<Key>& m)
{
    return m.size();
}

/*!
 * @brief Template Deduction Guide
 *
 * @tparam Key
 */
// template <typename Key>
// set(std::initializer_list<Key>) -> set<Key>;

// template <typename Key>
// set(std::initializer_list<const char*> ) -> set<std::string>;

template <typename Iter>
struct key_iterator : Iter
{
    explicit key_iterator(Iter it)
        : Iter(it)
    {
    }
    auto operator*() const
    {
        return Iter::operator*().first;
    }
    key_iterator& operator++()
    {
        Iter::operator++();
        return *this;
    }
};

/*!
 * @brief
 *
 * @tparam Key
 * @tparam T
 */
template <typename Key, typename T>
class dict : public std::unordered_map<Key, T>
{
    using _Self = dict<Key, T>;
    using _Base = std::unordered_map<Key, T>;

  public:
    using value_type = std::pair<const Key, T>;

    /*!
     * @brief Construct a new dict object
     *
     */
    dict()
        : std::unordered_map<Key, T> {}
    {
    }

    /*!
     * @brief Construct a new dict object
     *
     * @param[in] init
     */
    dict(std::initializer_list<value_type> init)
        : std::unordered_map<Key, T> {init}
    {
    }

    /*!
     * @brief Construct a new dict object
     *
     * @tparam Sequence
     * @param[in] S
     */
    // template <class Sequence>
    // explicit dict(const Sequence &S) {
    //     this->reserve(S.size());
    //     for (auto&& [i_v, v] : py::enumerate(S)) {
    //         (*this)[v] = i_v;
    //     }
    // }

    /*!
     * @brief
     *
     * @param[in] key
     * @return true
     * @return false
     */
    bool contains(const Key& key) const
    {
        return this->find(key) != this->end();
    }

    /*!
     * @brief
     *
     * @param[in] key
     * @param[in] default_value
     * @return T
     */
    T get(const Key& key, const T& default_value)
    {
        if (!contains(key))
        {
            return default_value;
        }
        return (*this)[key];
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto begin() const
    {
        using Iter = decltype(std::unordered_map<Key, T>::begin());
        return key_iterator<Iter> {std::unordered_map<Key, T>::begin()};
    }

    /*!
     * @brief
     *
     * @return auto
     */
    auto end() const
    {
        using Iter = decltype(std::unordered_map<Key, T>::end());
        return key_iterator<Iter> {std::unordered_map<Key, T>::end()};
    }

    /*!
     * @brief
     *
     * @return std::unordered_map<Key, T>&
     */
    std::unordered_map<Key, T>& items()
    {
        return *this;
    }

    /*!
     * @brief
     *
     * @return const std::unordered_map<Key, T>&
     */
    const std::unordered_map<Key, T>& items() const
    {
        return *this;
    }

    /*!
     * @brief
     *
     * @return _Self
     */
    _Self copy() const
    {
        return *this;
    }

    /*!
     * @brief
     *
     * @return _Self&
     */
    _Self& operator=(const _Self&) = delete;

    /*!
     * @brief
     *
     * @return _Self&
     */
    _Self& operator=(_Self&&) noexcept = default;

    /*!
     * @brief Move Constructor (default)
     *
     */
    dict(dict<Key, T>&&) noexcept = default;

  private:
    /*!
     * @brief Construct a new dict object
     *
     * Copy through explicitly the public copy() function!!!
     */
    dict(const dict<Key, T>&) = default;
};

/*!
 * @brief
 *
 * @tparam Key
 * @tparam T
 * @param[in] key
 * @param[in] m
 * @return true
 * @return false
 */
template <typename Key, typename T>
inline bool operator<(const Key& key, const dict<Key, T>& m)
{
    return m.contains(key);
}

/*!
 * @brief
 *
 * @tparam Key
 * @tparam T
 * @param[in] m
 * @return size_t
 */
template <typename Key, typename T>
inline size_t len(const dict<Key, T>& m)
{
    return m.size();
}

/*!
 * @brief Template Deduction Guide
 *
 * @tparam Key
 * @tparam T
 */
// template <typename Key, typename T>
// dict(std::initializer_list<std::pair<const Key, T>>) -> dict<Key, T>;

// template <class Sequence>
// dict(const Sequence& S)
//     -> dict<std::remove_cv_t<decltype(*std::begin(S))>, size_t>;

} // namespace py
