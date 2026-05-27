#ifndef SEQTRIE_RADIX_CHILD_MAP_H
#define SEQTRIE_RADIX_CHILD_MAP_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace seqtrie {

// Small insertion-ordered child map specialized for RadixMap nodes.
// Keys are concrete char values and the inline capacity is fixed at four,
// which covers the common DNA/RNA alphabet fanout without heap allocation.
template <typename Value>
class radix_child_map {
public:
  using key_type    = char;
  using mapped_type = Value;
  using size_type   = std::uint8_t;

  static constexpr size_type inline_capacity = 4;
  static constexpr size_type max_capacity = std::numeric_limits<size_type>::max();

  static_assert(std::is_default_constructible<mapped_type>::value,
                "radix_child_map value must be default constructible");
  static_assert(std::is_move_assignable<mapped_type>::value,
                "radix_child_map value must be move assignable");
  static_assert(std::is_nothrow_move_assignable<mapped_type>::value,
                "radix_child_map value must be nothrow move assignable");

  // Reference proxy returned by iterator dereference.
  // Range-for should use `auto ch : c`; members hold real references.
  struct value_proxy {
    const key_type& first;
    mapped_type&    second;
  };
  struct const_value_proxy {
    const key_type&    first;
    const mapped_type& second;
  };

private:
  struct inline_storage {
    key_type keys[inline_capacity];
    mapped_type vals[inline_capacity];

    inline_storage() : keys{} {}
  };

  struct heap_storage {
    key_type* keys;
    mapped_type* vals;
  };

  union storage_u {
    inline_storage inl;
    heap_storage heap;

    storage_u() {}
    ~storage_u() {}
  } u_;

  size_type size_;
  size_type cap_; // inline_capacity => inline; > inline_capacity => heap

  bool is_inline() const noexcept { return cap_ == inline_capacity; }

  void construct_inline() {
    ::new (&u_.inl) inline_storage();
  }

  void construct_heap(key_type* keys, mapped_type* vals) noexcept {
    ::new (&u_.heap) heap_storage{keys, vals};
  }

  key_type*       keys_ptr()       noexcept { return is_inline() ? u_.inl.keys : u_.heap.keys; }
  const key_type* keys_ptr() const noexcept { return is_inline() ? u_.inl.keys : u_.heap.keys; }
  mapped_type*       vals_ptr()       noexcept { return is_inline() ? u_.inl.vals : u_.heap.vals; }
  const mapped_type* vals_ptr() const noexcept { return is_inline() ? u_.inl.vals : u_.heap.vals; }

  static void alloc_heap(const size_type n, key_type*& out_keys, mapped_type*& out_vals) {
    auto keys = std::make_unique<key_type[]>(n);
    auto vals = std::make_unique<mapped_type[]>(n);
    out_keys = keys.release();
    out_vals = vals.release();
  }

  void destroy_storage() noexcept {
    if(is_inline()) {
      u_.inl.~inline_storage();
    } else {
      delete[] u_.heap.keys;
      delete[] u_.heap.vals;
      u_.heap.~heap_storage();
    }
  }

  static size_type next_capacity(const size_type current, const std::size_t need) {
    if(need > max_capacity) {
      throw std::length_error("radix_child_map cannot hold more than 255 children");
    }

    size_type nc = current;
    while(static_cast<std::size_t>(nc) < need) {
      if(nc > static_cast<size_type>(max_capacity / 2)) {
        nc = max_capacity;
      } else {
        nc = static_cast<size_type>(nc * 2);
      }
    }
    return nc;
  }

  void grow(const size_type new_cap) {
    key_type* nk = nullptr;
    mapped_type* nv = nullptr;
    alloc_heap(new_cap, nk, nv);

    key_type* ok = keys_ptr();
    mapped_type* ov = vals_ptr();

    std::memcpy(nk, ok, static_cast<std::size_t>(size_) * sizeof(key_type));
    for(size_type i = 0; i < size_; ++i) {
      nv[i] = std::move(ov[i]);
    }

    destroy_storage();
    construct_heap(nk, nv);
    cap_ = new_cap;
  }

  void ensure_capacity(const std::size_t need) {
    if(need <= cap_) return;
    grow(next_capacity(cap_, need));
  }

public:
  class iterator {
    friend class radix_child_map;
    radix_child_map* c_;
    size_type i_;
    iterator(radix_child_map* c, size_type i) : c_(c), i_(i) {}
  public:
    iterator() : c_(nullptr), i_(0) {}
    value_proxy operator*() const { return value_proxy{ c_->keys_ptr()[i_], c_->vals_ptr()[i_] }; }

    struct arrow {
      value_proxy v;
      value_proxy* operator->() { return &v; }
    };
    arrow operator->() const { return arrow{ value_proxy{ c_->keys_ptr()[i_], c_->vals_ptr()[i_] } }; }

    iterator& operator++() { ++i_; return *this; }
    iterator operator++(int) { iterator t = *this; ++i_; return t; }
    bool operator==(const iterator& o) const { return c_ == o.c_ && i_ == o.i_; }
    bool operator!=(const iterator& o) const { return !(*this == o); }
    size_type index() const { return i_; }
  };

  class const_iterator {
    friend class radix_child_map;
    const radix_child_map* c_;
    size_type i_;
    const_iterator(const radix_child_map* c, size_type i) : c_(c), i_(i) {}
  public:
    const_iterator() : c_(nullptr), i_(0) {}
    const_value_proxy operator*() const { return const_value_proxy{ c_->keys_ptr()[i_], c_->vals_ptr()[i_] }; }

    struct arrow {
      const_value_proxy v;
      const_value_proxy* operator->() { return &v; }
    };
    arrow operator->() const { return arrow{ const_value_proxy{ c_->keys_ptr()[i_], c_->vals_ptr()[i_] } }; }

    const_iterator& operator++() { ++i_; return *this; }
    const_iterator operator++(int) { const_iterator t = *this; ++i_; return t; }
    bool operator==(const const_iterator& o) const { return c_ == o.c_ && i_ == o.i_; }
    bool operator!=(const const_iterator& o) const { return !(*this == o); }
    size_type index() const { return i_; }
  };

  radix_child_map() : size_(0), cap_(inline_capacity) {
    construct_inline();
  }

  ~radix_child_map() {
    destroy_storage();
  }

  radix_child_map(const radix_child_map&) = delete;
  radix_child_map& operator=(const radix_child_map&) = delete;

  radix_child_map(radix_child_map&& o) noexcept : size_(0), cap_(inline_capacity) {
    construct_inline();
    move_from(o);
  }

  radix_child_map& operator=(radix_child_map&& o) noexcept {
    if(this != &o) {
      move_from(o);
    }
    return *this;
  }

  void swap(radix_child_map& o) noexcept {
    radix_child_map tmp;
    tmp.move_from(*this);
    move_from(o);
    o.move_from(tmp);
  }

private:
  void move_from(radix_child_map& o) noexcept {
    destroy_storage();

    size_ = o.size_;
    cap_ = o.cap_;

    if(o.is_inline()) {
      construct_inline();
      for(size_type i = 0; i < o.size_; ++i) {
        u_.inl.keys[i] = o.u_.inl.keys[i];
        u_.inl.vals[i] = std::move(o.u_.inl.vals[i]);
      }
      o.size_ = 0;
    } else {
      construct_heap(o.u_.heap.keys, o.u_.heap.vals);
      o.u_.heap.keys = nullptr;
      o.u_.heap.vals = nullptr;
      o.u_.heap.~heap_storage();
      o.size_ = 0;
      o.cap_ = inline_capacity;
      o.construct_inline();
    }
  }

public:
  size_type size() const noexcept { return size_; }
  bool empty() const noexcept { return size_ == 0; }

  void reserve(const std::size_t need) {
    ensure_capacity(need);
  }

  void reserve_exact(const std::size_t need) {
    if(need > max_capacity) {
      throw std::length_error("radix_child_map cannot hold more than 255 children");
    }
    if(need <= cap_) return;
    grow(static_cast<size_type>(need));
  }

  void append_unchecked(const key_type k, mapped_type v) noexcept {
    const size_type i = size_++;
    keys_ptr()[i] = k;
    vals_ptr()[i] = std::move(v);
  }

  iterator begin() noexcept { return iterator(this, 0); }
  iterator end() noexcept { return iterator(this, size_); }
  const_iterator begin() const noexcept { return const_iterator(this, 0); }
  const_iterator end() const noexcept { return const_iterator(this, size_); }
  const_iterator cbegin() const noexcept { return begin(); }
  const_iterator cend() const noexcept { return end(); }

  iterator find(const key_type k) noexcept {
    const key_type* ks = keys_ptr();
    for(size_type i = 0; i < size_; ++i) {
      if(ks[i] == k) return iterator(this, i);
    }
    return end();
  }

  const_iterator find(const key_type k) const noexcept {
    const key_type* ks = keys_ptr();
    for(size_type i = 0; i < size_; ++i) {
      if(ks[i] == k) return const_iterator(this, i);
    }
    return end();
  }

  mapped_type& at(const key_type k) {
    iterator it = find(k);
    if(it == end()) throw std::out_of_range("radix_child_map::at");
    return vals_ptr()[it.index()];
  }

  const mapped_type& at(const key_type k) const {
    const_iterator it = find(k);
    if(it == end()) throw std::out_of_range("radix_child_map::at");
    return vals_ptr()[it.index()];
  }

  mapped_type& operator[](const key_type k) {
    iterator it = find(k);
    if(it != end()) return vals_ptr()[it.index()];
    if(size_ == max_capacity) {
      throw std::length_error("radix_child_map cannot hold more than 255 children");
    }
    ensure_capacity(static_cast<std::size_t>(size_) + 1);
    const size_type i = size_++;
    keys_ptr()[i] = k;
    vals_ptr()[i] = mapped_type{};
    return vals_ptr()[i];
  }

  std::pair<iterator, bool> emplace(const key_type k, mapped_type v) {
    iterator it = find(k);
    if(it != end()) return { it, false };
    if(size_ == max_capacity) {
      throw std::length_error("radix_child_map cannot hold more than 255 children");
    }
    ensure_capacity(static_cast<std::size_t>(size_) + 1);
    const size_type i = size_++;
    keys_ptr()[i] = k;
    vals_ptr()[i] = std::move(v);
    return { iterator(this, i), true };
  }

  iterator erase(iterator it) {
    const size_type i = it.index();
    const size_type last = static_cast<size_type>(size_ - 1);
    if(i != last) {
      keys_ptr()[i] = keys_ptr()[last];
      vals_ptr()[i] = std::move(vals_ptr()[last]);
    }
    vals_ptr()[last] = mapped_type{};
    --size_;
    return iterator(this, i);
  }
};

} // namespace seqtrie

#endif // SEQTRIE_RADIX_CHILD_MAP_H
