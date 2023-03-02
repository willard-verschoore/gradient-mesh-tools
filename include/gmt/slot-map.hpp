#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

namespace gmt
{

/// An index into a SlotMap.
/**
 * Serves a double purpose: holds an index into the indirection list, which
 * itself holds indices to the data array. Thus, every array access goes through
 * two indirections.
 *
 * The type parameter T is not actually necessary, as a handle only holds
 * integers. It's present to strengthen the connection between a SlotMap and its
 * handles.
 */
template <typename T>
struct SlotMapHandle
{
  std::uint32_t id = 0, version = std::numeric_limits<uint32_t>::max();
};

/// An associative container of SlotMapHandle indices to items of type T.
/**
 * The items are guaranteed to be stored consecutively so that iterating over
 * them is fast. This is achieved by maintaining an indirection list which is
 * indexed using a SlotMapHandle. Accessing an item using its handle therefore
 * goes through an additional layer of indirection, but it is still O(1).
 */
template <typename T>
class SlotMap
{
 public:
  // Internal storage is a packed std::vector<T>, so iterator types are just
  // aliases for the vector's.
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  // Erases all elements from the map.
  void clear()
  {
    data.clear();
    indirection.clear();
    free_list.clear();
  }

  /// Inserts \c item into the map, returning its identifier handle.
  /**
   * @param item The item to add to the map.
   * @return A SlotMapHandle identifying the newly added item in map.
   */
  SlotMapHandle<T> add(T const &item)
  {
    // Internal data is always packed and contiguous, so we do a regular
    // push_back into it.
    std::uint32_t index = size();
    data.push_back(item);

    // Invariant: free list is non-empty on insertion.
    if (free_list.empty())
    {
      free_list.push_back(index);
      indirection.push_back(SlotMapHandle<std::uint32_t>{index, 0});
    }

    // Consume most recent spot in free list.
    auto free = free_list.back();
    free_list.pop_back();

    assert(free < indirection.size());

    // Update selected element in indirection list to point to the newly added
    // item in the data array.
    auto &handle = indirection[free];
    handle.id = index;

    return SlotMapHandle<T>{free, handle.version};
  }

  /// Locates the item identified by \c handle.
  /**
   * If the item has been removed, end() is returned. Thus, to check if an
   * element is present in the map:
   *
   * ```
   * auto it = map.find(handle);
   * if (it != map.end()) // Element exists.
   *   auto &elem = *it;
   * ```
   *
   * @param handle A SlotMapHandle identifying the item to be found.
   * @return An iterator to the item in the map's data array.
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

  /// Locates the item identified by \c handle.
  /**
   * If the item has been removed, end() is returned. Thus, to check if an
   * element is present in the map:
   *
   * ```
   * auto it = map.find(handle);
   * if (it != map.end()) // Element exists.
   *   auto &elem = *it;
   * ```
   *
   * @param handle A SlotMapHandle identifying the item to be found.
   * @return An iterator to the item in the map's data array.
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

  /// Returns a reference to the item identified by \c handle.
  /**
   * If the handle doesn't point to a valid element, this function will invoke
   * undefined behavior: make sure to always check contains() before calling
   * this function!
   *
   * @param handle A SlotMapHandle identifying the item to be found.
   * @return A reference to the requested item.
   */
  T const &operator[](SlotMapHandle<T> handle) const
  {
    auto iter = find(handle);
    assert(iter != end());
    return *iter;
  }

  /// Returns a reference to the item identified by \c handle.
  /**
   * If the handle doesn't point to a valid element, this function will invoke
   * undefined behavior: make sure to always check contains() before calling
   * this function!
   *
   * @param handle A SlotMapHandle identifying the item to be found.
   * @return A reference to the requested item.
   */
  T &operator[](SlotMapHandle<T> handle)
  {
    auto iter = find(handle);
    assert(iter != end());
    return *iter;
  }

  /// Removes the item identified by \c handle from the map.
  /**
   * Future calls to find() with this handle will thus return end().

   * @param handle A SlotMapHandle identifying the item to be removed.
   */
  void remove(SlotMapHandle<T> handle)
  {
    auto &obj_handle = indirection.at(handle.id);
    if (obj_handle.version == handle.version)
    {
      // Increment the version number on deletion, so access from copied handles
      // will fail.
      obj_handle.version++;

      // Mark spot taken up by removed item as free to use.
      free_list.push_back(handle.id);

      // Use the swap-and-pop trick to remove the item at the given index.
      std::swap(data.at(obj_handle.id), data.back());
      data.pop_back();

      // Find which slot pointed to the moved item and update it.
      std::uint32_t moved_index = data.size();
      for (std::uint32_t slot = 0; slot < indirection.size(); ++slot)
      {
        // The slot's id must match the moved index and it should be valid (i.e.
        // not in the free list).
        if (indirection[slot].id == moved_index && is_valid(slot))
        {
          indirection[slot].id = obj_handle.id;
          break;
        }
      }
    }
  }

  /// Checks whether \c handle identifies a valid item in the map.
  /**
   * This will return \c false if remove() has been called with \c handle
   * previously.
   *
   * @param handle The SlotMapHandle to be checked.
   * @return \c true if \c handle identifies a valid item in the map, \c false
   * otherwise.
   */
  bool contains(SlotMapHandle<T> handle) const
  {
    return handle.id < indirection.size() &&
           indirection[handle.id].version == handle.version;
  }

  /// Gets the handle identifying the item at \c position.
  /**
   * If \c position is not a valid iterator this will return a default
   * constructed SlotMapHandle.
   *
   * @param position An iterator into the map's data pointing to the item for
   * which to get the handle.
   * @return A SlotMapHandle identifying the item at \c position if found.
   */
  SlotMapHandle<T> get_handle(const_iterator position)
  {
    std::uint32_t index = position - begin();
    for (std::uint32_t slot = 0; slot < indirection.size(); ++slot)
    {
      if (indirection[slot].id == index && is_valid(slot))
        return SlotMapHandle<T>{slot, indirection[slot].version};
    }

    return SlotMapHandle<T>{};
  }

  /// Returns the number of items currently contained in the map.
  /**
   * @return The number of items currently contained in the map.
   */
  std::uint32_t size() const { return data.size(); }

  /// Gets the begin iterator to access the items in the map sequentially.
  /**
   * No guarantees are made about the order of the items in the map, especially
   * if remove() has been called many times before.
   *
   * @return An iterator to the first item in the map.
   */
  iterator begin() { return data.begin(); }

  /// Gets the begin iterator to access the items in the map sequentially.
  /**
   * No guarantees are made about the order of the items in the map, especially
   * if remove() has been called many times before.
   *
   * @return An iterator to the first item in the map.
   */
  const_iterator begin() const { return data.begin(); }

  /// Gets the end iterator to access the items in the map sequentially.
  /**
   * No guarantees are made about the order of the items in the map, especially
   * if remove() has been called many times before.
   *
   * @return An iterator to the item following the last item in the map.
   */
  iterator end() { return data.end(); }

  /// Gets the end iterator to access the items in the map sequentially.
  /**
   * No guarantees are made about the order of the items in the map, especially
   * if remove() has been called many times before.
   *
   * @return An iterator to the item following the last item in the map.
   */
  const_iterator end() const { return data.end(); }

 private:
  std::vector<T> data; // Underlying data storage.
  std::vector<SlotMapHandle<std::uint32_t>> indirection; // Indices into data.
  std::vector<std::uint32_t> free_list; // Free slots in indirection.

  /// Checks whether a slot in the indirection list is valid.
  /**
   * A slot is valid if it is within the bounds of the indirection list and the
   * slot is not contained in the free list.
   *
   * @param slot The index into the indirection list to check.
   * @return Whether \c slot is valid.
   */
  bool is_valid(std::uint32_t slot)
  {
    return slot < indirection.size() &&
           std::find(free_list.begin(), free_list.end(), slot) ==
               free_list.end();
  }
};

/// Checks whether two handles are equal.
/**
 * Two handles are equal if they have identical IDs and versions.
 *
 * @param lhs,rhs The two SlotMapHandle objects to compare.
 * @return \c true if the handles are equal, \c false otherwise.
 */
template <typename T>
inline bool operator==(SlotMapHandle<T> const &lhs, SlotMapHandle<T> const &rhs)
{
  return lhs.id == rhs.id && lhs.version == rhs.version;
}

/// Checks whether two handles are unequal.
/**
 * Two handles are unequal if they have different IDs and/or versions.
 *
 * @param lhs,rhs The two SlotMapHandle objects to compare.
 * @return \c true if the handles are unequal, \c false otherwise.
 */
template <typename T>
inline bool operator!=(SlotMapHandle<T> const &lhs, SlotMapHandle<T> const &rhs)
{
  return !(lhs == rhs);
}

/// Checks whether one handles precedes another.
/**
 * One handle precedes another if its ID is smaller.
 *
 * @param lhs,rhs The two SlotMapHandle objects to compare.
 * @return \c true if \c lhs precedes \c rhs, \c false otherwise.
 */
template <typename T>
inline bool operator<(SlotMapHandle<T> const &lhs, SlotMapHandle<T> const &rhs)
{
  return lhs.id < rhs.id;
}

} // namespace gmt
