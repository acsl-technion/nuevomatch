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

#include <algorithm> // copy
#include <sstream>
#include <string>

#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include <object_io.h>
#include <logging.h>

/**
 * @brief Initialize an empty packer, used for packing data
 */
ObjectPacker::ObjectPacker() : _size(0) {
	_buffer = new std::stringstream();
}

ObjectPacker::~ObjectPacker() {
	delete static_cast<std::stringstream*>(_buffer);
}

ObjectPacker::ObjectPacker(ObjectPacker&& other) : _size(other._size), _buffer(other._buffer) {
	other._buffer = nullptr;
}

ObjectPacker& ObjectPacker::operator=(ObjectPacker&& other) {
	if (&other == this) {
		return *this;
	}
	delete static_cast<std::stringstream*>(_buffer);
	_buffer = other._buffer;
	_size = other._size;
	other._buffer = nullptr;
	return *this;
}

/**
 * @brief Insert any object to the packer at the end
 */
template <typename T>
ObjectPacker& ObjectPacker::operator<<(const T& obj) {
	static_cast<std::stringstream*>(_buffer)->write( (const char*)&obj, sizeof(obj) );
	_size+=sizeof(obj);
	return *this;
}

/**
 * @brief Adds another packer at the end of this
 */
template <>
ObjectPacker& ObjectPacker::operator<<(const ObjectPacker& obj) {
	unsigned int size;
	unsigned char* data;
	obj.pack(&data, &size);

	*this << size;
	static_cast<std::stringstream*>(_buffer)->write( (const char*)data, size );
	_size+=size;

	delete data;
	return *this;
}

/**
 * @brief Adds byte array at the end of this
 */
ObjectPacker& ObjectPacker::push(void* buffer, unsigned int size) {
	static_cast<std::stringstream*>(_buffer)->write( (const char*)buffer, size );
	_size += size;
	return *this;
}

/**
 * @brief Inserts data to the beginning of the packer
 */
template <typename T>
ObjectPacker& ObjectPacker::insert(const T& obj) {
	ObjectPacker tmp;
	tmp << obj;

	unsigned int size;
	unsigned char* data;
	this->pack(&data, &size);
	tmp.push(data, size);
	delete data;

	*this = std::move(tmp);
	return *this;
}

/**
 * @brief Inserts packer to the beginning of this
 */
template <>
ObjectPacker& ObjectPacker::insert(const ObjectPacker& obj) {

	unsigned int my_size;
	unsigned char* my_data;
	this->pack(&my_data, &my_size);

	unsigned int lhs_size;
	unsigned char* lhs_data;
	obj.pack(&lhs_data, &lhs_size);

	ObjectPacker tmp;
	tmp.push(lhs_data, lhs_size);
	tmp << my_size;
	tmp.push(my_data, my_size);

	delete my_data;
	delete lhs_data;

	*this = std::move(tmp);
	return *this;
}

/**
 * @brief Packs this into byte array.
 * @param[out] dst The destination buffer
 * @param[out] dst_size The destination size
 */
void ObjectPacker::pack(unsigned char** dst, unsigned int* dst_size) const {
	*dst = new unsigned char[_size];
	std::string my_data = static_cast<std::stringstream*>(_buffer)->str();
	const char* char_data = my_data.c_str();
	std::copy(char_data, char_data + _size, *dst);
	*dst_size = _size;
}

/**
 * @brief Initialize empty reader
 */
ObjectReader::ObjectReader()
	: _buffer(nullptr), _size(0), _mmap_used(false), _org_buffer(nullptr), _org_size(0) {}

/**
 * @brief Initialize from buffer and size
 */
ObjectReader::ObjectReader(void* buffer, unsigned int size)
	: _buffer((unsigned char*)buffer), _size(size),
	  _mmap_used(false), _org_buffer(nullptr), _org_size(0) { }


/**
 * @brief Initialize an ObjectReader from an ObjectPacker
 * @note The raw binary data is copied from rhs to this
 */
ObjectReader::ObjectReader(const ObjectPacker& rhs)
	: _buffer(nullptr), _size(0),
	  _mmap_used(false), _org_buffer(nullptr), _org_size(0)
{
	rhs.pack(&_buffer, &_size);
}

/**
 * @brief Initialize from filename
 * @throws In case of file read error
 */
ObjectReader::ObjectReader(const char* filename) {

	// Check that file file exist
	if (access(filename, F_OK) == -1) {
		throw error("Cannot read from file: file '" << filename << "' does not exist");
	}

	int fd = open(filename, O_RDONLY);

	// Check the total length of the file
	struct stat file_stat;
	fstat(fd, &file_stat);
	this->_size = file_stat.st_size;
	this->_org_size = file_stat.st_size;

	// Open the file as memory map
	// info("Reading from file " << filename << " as memory map (length " << this->_size << " bytes)");
	void* shmem = mmap(NULL, this->_size, PROT_READ, MAP_SHARED, fd, 0);

	// Check errors
	if (shmem == MAP_FAILED) {
		throw error("Cannot read file: cannot open memory map: " << strerror(errno));
	}

	// No need the file descriptor to be opened
	if (close(fd) != 0 ) {
		// Error
		throw error("Cannot close file descriptor: " << strerror(errno));
	}

	this->_buffer = (unsigned char*)shmem;
	this->_org_buffer = (unsigned char*)shmem;
	this->_mmap_used = true;
}

ObjectReader::~ObjectReader() {
	if (this->_mmap_used) {
		if (munmap(this->_org_buffer, this->_mmap_used) != 0){
			warning("Cannot detach memory map: " << strerror(errno));
		}
	}
}

/**
 * @brief Read object reader from the reader
 */
template <>
ObjectReader& ObjectReader::operator>>(ObjectReader& obj) {
	unsigned int new_size = 0;
	*this >> new_size;
	if (_size < new_size) {
		throw error("Cannot extract object from reader: buffer overflow");
	}
	obj._buffer = _buffer;
	obj._size = new_size;
	_buffer += +new_size;
	_size -= new_size;
	return *this;
}

/**
 * @brief Read any object from the packer
 * @throws runtime_error in case of buffer overflow
 */
template <typename T>
ObjectReader& ObjectReader::operator>>(T& obj) {
	if (_size < sizeof(T)) {
		throw error("Cannot read object from reader: buffer overflow");
	}
	obj = *(T*)_buffer;
	_size -= sizeof(T);
	_buffer=(unsigned char*)((T*)_buffer+1);
	return *this;
}

/**
 * @brief Read any object from the packer
 * @throws runtime_error in case of buffer overflow
 */
template <typename T>
T ObjectReader::read() {
	T output;
	(*this) >> output;
	return output;
}

/**
 * @brief Extract an object-reader object from buffer
 */
ObjectReader ObjectReader::extract() {
	ObjectReader output;
	(*this) >> output;
	return output;
}

// Initiates templates
#define INIT_TEMPLATE(X) \
	template ObjectReader& ObjectReader::operator>>(X&);		\
	template X ObjectReader::read();							\
	template ObjectPacker& ObjectPacker::operator<<(const X&);	\
	template ObjectPacker& ObjectPacker::insert(const X&);

INIT_TEMPLATE(unsigned long);
INIT_TEMPLATE(unsigned int);
INIT_TEMPLATE(unsigned char);
INIT_TEMPLATE(unsigned int volatile);
INIT_TEMPLATE(long);
INIT_TEMPLATE(int);
INIT_TEMPLATE(char);
INIT_TEMPLATE(bool);
INIT_TEMPLATE(float);
INIT_TEMPLATE(double);
