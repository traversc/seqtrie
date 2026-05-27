#ifndef seqtrie_NODE_POOL_H
#define seqtrie_NODE_POOL_H

#include <cstddef>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace seqtrie {

template <typename T>
class NodePool {
public:
  explicit NodePool(size_t block_size = 4096) : block_size_(block_size) {}

  ~NodePool() { clear(); }

  NodePool(const NodePool &) = delete;
  NodePool & operator=(const NodePool &) = delete;
  NodePool(NodePool &&) = delete;
  NodePool & operator=(NodePool &&) = delete;

  void reserve_exact(size_t capacity) {
    if(capacity == 0) return;
    if(!blocks_.empty()) {
      throw std::logic_error("NodePool::reserve_exact requires an empty pool");
    }
    add_block(capacity);
  }

  template <typename... Args>
  T * create(Args&&... args) {
    if(blocks_.empty() || blocks_.back().used == blocks_.back().capacity) {
      add_block();
    }
    Block & block = blocks_.back();
    auto * data = static_cast<T *>(block.storage);
    T * ptr = data + block.used;
    ::new (static_cast<void *>(ptr)) T(std::forward<Args>(args)...);
    ++block.used;
    return ptr;
  }

  size_t allocated_nodes() const noexcept {
    size_t n = 0;
    for(const auto & block : blocks_) n += block.used;
    return n;
  }

  template <typename Func>
  void for_each_node(Func && func) const {
    for(const auto & block : blocks_) {
      const auto * data = static_cast<const T *>(block.storage);
      for(size_t i = 0; i < block.used; ++i) {
        func(data[i]);
      }
    }
  }

  void swap(NodePool & other) noexcept {
    using std::swap;
    swap(block_size_, other.block_size_);
    swap(blocks_, other.blocks_);
  }

private:
  struct Block {
    void * storage = nullptr;
    size_t used = 0;
    size_t capacity = 0;
  };

  void add_block() {
    add_block(block_size_);
  }

  void add_block(size_t capacity) {
    if(capacity == 0) {
      throw std::logic_error("NodePool block capacity must be greater than zero");
    }
    Block block;
    block.capacity = capacity;
    block.storage = ::operator new(sizeof(T) * block.capacity);
    try {
      blocks_.push_back(block);
    } catch(...) {
      ::operator delete(block.storage);
      throw;
    }
  }

  void clear() noexcept {
    for(auto & block : blocks_) {
      auto * data = static_cast<T *>(block.storage);
      for(size_t i = 0; i < block.used; ++i) {
        data[i].~T();
      }
      ::operator delete(block.storage);
    }
  }

  size_t block_size_;
  std::vector<Block> blocks_;
};

} // namespace seqtrie

#endif // seqtrie_NODE_POOL_H
