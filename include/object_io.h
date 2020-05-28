/*
 * MIT License
 * Copyright (c) 2019 Alon Rashelbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <basic_types.h>

/**
 * @brief Used for packing information to buffer
 */
class ObjectPacker {
private:

	unsigned int _size;
	void* _buffer;

public:

	/**
	 * @brief Initialize an empty packer
	 */
	ObjectPacker();
	~ObjectPacker();

	ObjectPacker(const ObjectPacker&) = delete;
	ObjectPacker(ObjectPacker&&);

	ObjectPacker& operator=(const ObjectPacker&) = delete;
	ObjectPacker& operator=(ObjectPacker&&);

	/**
	 * @brief Insert any object to the packer at the end
	 */
	template <typename T>
	ObjectPacker& operator<<(const T& obj);

	/**
	 * @brief Adds byte array at the end of this
	 */
	ObjectPacker& push(void* buffer, unsigned int size);

	/**
	 * @brief Inserts data to the beginning of the packer
	 * @note This method requires two copies
	 */
	template <typename T>
	ObjectPacker& insert(const T& obj);

	/**
	 * @brief Packs this into byte array.
	 * @param[out] dst The destination buffer
	 * @param[out] dst_size The destination size
	 */
	void pack(unsigned char** dst, unsigned int* dst_size) const;

	/**
	 * @brief Returns the size of this
	 */
	unsigned int size() const { return _size; }
};

/**
 * @brief Used for reading information from buffer
 */
class ObjectReader {
private:
	unsigned char* _buffer;
	unsigned int _size;

	// Used for mmap
	bool _mmap_used;
	unsigned char* _org_buffer;
	unsigned int _org_size;

public:

	/**
	 * @brief Initialize empty reader
	 */
	ObjectReader();

	/**
	 * @brief Initialize from buffer and size
	 */
	ObjectReader(void* buffer, unsigned int size);

	/**
	 * @brief Initialize an ObjectReader from an ObjectPacker
	 * @note The raw binary data is copied from rhs to this
	 */
	ObjectReader(const ObjectPacker& rhs);

	/**
	 * @brief Initialize from filename
	 * @throws In case of file read error
	 */
	ObjectReader(const char* filename);

	virtual ~ObjectReader();

	/**
	 * @brief Read any object from the packer
	 * @throws runtime_error in case of buffer overflow
	 */
	template <typename T>
	ObjectReader& operator>>(T& obj);

	/**
	 * @brief Read any object from the packer
     * @throws runtime_error in case of buffer overflow
	 */
	template <typename T>
	T read();

	/**
	 * @brief Extract an object-reader object from buffer
	 */
	ObjectReader extract();

	/**
	 * @brief Returns the size of this
	 */
	unsigned int size() { return _size; }

	/**
	 * @brief Returns a pointer to the buffer of this
	 */
	void* buffer() { return _buffer; }
};
