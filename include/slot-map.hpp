#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

/// An index into a SlotMap.
/**
 * Serves a double purpose: holds an index into the indirection list,
 * which itself holds indices to the data array. Thus, every array access
 * goes through two indirections.
 *
 * The type parameter T is not actually necessary, as a handle only holds
 * integers. It's present to strengthten the connection between a SlotMap<T> and
 * its handles.
 */
template <typename T>
struct SlotMapHandle
{
  std::uint32_t id = 0, version = std::numeric_limits<uint32_t>::max();
};

/// An associative container of SlotMapHandles to items of type T.
template <typename T>
class SlotMap
{
 public:
  // Internal storage is a packed std::vector<T>, so iterator
  // types are just aliases to the vector's.
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  void clear()
  {
    data.clear();
    free_list.clear();
    indirection.clear();
  }

  /// Inserts item into the map, returning its identifier handle.
  SlotMapHandle<T> add(T item)
  {
    // Internal data is always packed and contiguous,
    // so we do a regular push_back into it.
    std::uint32_t item_idx = size();
    data.push_back(std::move(item));

    // Invariant: free list is non-empty on insertion.
    if (free_list.empty())
    {
      free_list.push_back(item_idx);
      indirection.push_back(SlotMapHandle<std::uint32_t>{item_idx, 0});
    }

    // Consume most recent spot in free list.
    auto free = free_list.back();
    free_list.pop_back();

    assert(free < indirection.size());
    // Update selected element in indirection list
    // to point to the newly added item in the data array.
    auto& handle = indirection[free];
    handle.id = item_idx;

    return SlotMapHandle<T>{free, handle.version};
  }

  /// Locates the item indexed by the given handle.
  /**
   * If the item has been removed, end() is returned.
   * Thus, to check if an element is present in the map:
   *
   * ```
   * auto it = map.find(handle);
   * if (it != map.end()) {
   *  Element exists!
   *	auto& elem = *it;
   * }
   * ```
   */
  const_iterator find(SlotMapHandle<T> handle) const
  {
    if (handle.id < indirection.size())
    {
      auto obj_handle = indirection[handle.id];
      if (obj_handle.version == handle.version)
      {
        return begin() + obj_handle.id;
      }
    }
    return end();
  }

  /// Locates the item indexed by the given handle.
  /**
   * If the item has been removed, end() is returned.
   * Thus, to check if an element is present in the map:
   *
   * ```
   * auto it = map.find(handle);
   * if (it != map.end()) {
   *  Element exists!
   *	auto& elem = *it;
   * }
   * ```
   */
  iterator find(SlotMapHandle<T> handle)
  {
    if (handle.id < indirection.size())
    {
      auto obj_handle = indirection[handle.id];
      if (obj_handle.version == handle.version)
      {
        return begin() + obj_handle.id;
      }
    }
    return end();
  }

  /// Returns a reference to the item indexed by the given handle.
  /**
   * If the handle doesn't point to a valid element, this function
   * will invoke undefined behavior: make sure to always check contains()
   * before calling this function!
   */
  const T& operator[](SlotMapHandle<T> handle) const
  {
    auto iter = find(handle);
    assert(iter != end());
    return *iter;
  }

  /// Returns a reference to the item indexed by the given handle.
  /**
   * If the handle doesn't point to a valid element, this function
   * will invoke undefined behavior: make sure to always check contains()
   * before calling this function!
   */
  T& operator[](SlotMapHandle<T> handle)
  {
    auto iter = find(handle);
    assert(iter != end());
    return *iter;
  }

  /// Removes the item indexed by the given handle from the map.
  /**
   * Future calls to find() with this handle will thus return end().
   */
  void remove(SlotMapHandle<T> handle)
  {
    auto& obj_handle = indirection.at(handle.id);
    if (obj_handle.version == handle.version)
    {
      // Increment the version number on deletion, so access
      // from copied handles will fail.
      obj_handle.version++;
      // Mark spot taken up by removed item as free to use.
      free_list.push_back(handle.id);

      // Use the swap-and-pop trick to remove the item
      // at the given index.
      std::swap(data.at(obj_handle.id), data.back());
      data.pop_back();

      auto previous_idx = data.size();
      // Find which slot pointed to the moved item and update it.
      // Slower than previous implementation, but that ignored the free_list.
      auto prev_slot = indirection.begin();
      while (prev_slot != indirection.end() &&
             (prev_slot->id != previous_idx ||
              std::find(free_list.begin(), free_list.end(),
                        prev_slot - indirection.begin()) != free_list.end()))
      {
        ++prev_slot;
      }
      if (prev_slot != indirection.end())
      {
        prev_slot->id = obj_handle.id;
      }
    }
  }

  /**
   * Returns whether or not the given handle points to a valid element
   * in the map (that hasn't been removed).
   */
  bool contains(SlotMapHandle<T> handle) const
  {
    return handle.id < indirection.size() &&
           indirection[handle.id].version == handle.version;
  }

  std::vector<SlotMapHandle<std::uint32_t>> get_indirection()
  {
    return indirection;
  }

  /// Returns the number of items currently contained in the map.
  std::uint32_t size() const { return data.size(); }

  /// Begin iterator to access the data contained in the map sequentially.
  iterator begin() { return data.begin(); }
  /// Begin iterator to access the data contained in the map sequentially.
  const_iterator begin() const { return data.begin(); }
  /// End iterator to access the data contained in the map sequentially.
  iterator end() { return data.end(); }
  /// End iterator to access the data contained in the map sequentially.
  const_iterator end() const { return data.end(); }

 private:
  std::vector<T> data;
  std::vector<SlotMapHandle<std::uint32_t>> indirection;
  std::vector<std::uint32_t> free_list;
};

/// Returns true if both handles have an equal ID and version.
template <typename T>
inline bool operator==(const SlotMapHandle<T>& a, const SlotMapHandle<T>& b)
{
  return a.id == b.id && a.version == b.version;
}

/// Returns true if both handles have an equal ID and version.
template <typename T>
inline bool operator<(const SlotMapHandle<T>& a, const SlotMapHandle<T>& b)
{
  return a.id < b.id;
}

/// Returns false if both handles have an equal ID and version.
template <typename T>
inline bool operator!=(const SlotMapHandle<T>& a, const SlotMapHandle<T>& b)
{
  return !(a == b);
}
