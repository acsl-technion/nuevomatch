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

#include <iostream>
#include <sstream>
#include <cstdarg>
#include <stdexcept>
#include <signal.h>
#include <semaphore.h>
#include <pthread.h>
#include <string.h>
#include <pthread.h>

#ifndef NDEBUG
 // Print a debug message to stderr using std::stream convention
 #define info(...) SimpleLogger::lock(); SimpleLogger::get() << \
 "Info: (" << __func__ <<") " << __VA_ARGS__ << \
 " (" << __FILE__ << ":" << __LINE__ << ")" << SimpleLogger::endl(); SimpleLogger::release()
 // Print a debug message to stderr using printf convention
 #define infof(...) SimpleLogger::lock(); SimpleLogger::get() << \
 "Info: (" << __func__ <<") " << SimpleLogger::format(__VA_ARGS__) << \
 " (" << __FILE__ << ":" << __LINE__ << ")" << SimpleLogger::endl(); SimpleLogger::release()
#else // NDEBUG
 #define info(...)
 #define infof(...)
#endif // NDEBUG

// Print a log message to stderr using std::stream convention
 #define logger(...) SimpleLogger::lock(); SimpleLogger::get() << \
 "Log: (" << __func__ <<") " << __VA_ARGS__ << \
 " (" << __FILE__ << ":" << __LINE__ << ")" << SimpleLogger::endl(); SimpleLogger::release()
 // Print a log message to stderr using printf convention
 #define loggerf(...) SimpleLogger::lock(); SimpleLogger::get() << \
 "Log: (" << __func__ <<") " << SimpleLogger::format(__VA_ARGS__) << \
 " (" << __FILE__ << ":" << __LINE__ << ")" << SimpleLogger::endl(); SimpleLogger::release()

// Print a message to stdout using std::stream convention
#define message_s(...) SimpleLogger::get(false) << SimpleLogger::lock() << __VA_ARGS__ << \
	SimpleLogger::endl(); SimpleLogger::release()

// Print a message to stdout using printf convention
#define messagef(...) SimpleLogger::get(false) << SimpleLogger::lock() << \
	SimpleLogger::format(__VA_ARGS__) << \
	SimpleLogger::endl(); SimpleLogger::release()

// Print a warning message to stderr using std::stream convention
#define warning(...) SimpleLogger::lock(); SimpleLogger::get() << \
	"Warning: (" << __func__ <<") " << __VA_ARGS__ << \
	" (" << __FILE__ << ":" << __LINE__ << ")" << SimpleLogger::endl(); SimpleLogger::release()
// Print a warning message to stderr using printf convention
#define warningf(...) SimpleLogger::lock(); SimpleLogger::get() << \
	"Warning: (" << __func__ <<") " << SimpleLogger::format(__VA_ARGS__) << \
	" (" << __FILE__ << ":" << __LINE__ << ")" << SimpleLogger::endl(); SimpleLogger::release()

// Create an exception with an arbitrary message using std::stream convention
#define error(...) SimpleException::create() <<  "Exception: (" <<  \
	__func__ << "@" << __FILE__ << ":" << __LINE__ << ") " << __VA_ARGS__
// Create an exception with an arbitrary message using printf convention
#define errorf(...) SimpleException::create() <<  "Exception: (" <<  \
	__func__ << "@" << __FILE__ << ":" << __LINE__ << ") " << \
	SimpleException::format(__VA_ARGS__)

/**
 * @brief Handles messages and errors
 */
class SimpleLogger {

	// Message buffer maximum size
	static constexpr int buffer_size = 4096;

	// Holds a string buffer
	char _buffer[buffer_size];
	int _buffer_cursor;

	// Use stderr for forceful printing
	bool _use_stderr;

	// A semaphore for syncing messages between threads
	sem_t _semaphore;

	// Sticky force
	bool _sticky_force;

	// Holds the singleton instance
	static SimpleLogger* _instance;

	SimpleLogger() {
		// Initiate the semaphore only once
		if (sem_init(&_semaphore, 0, 1)) {
			std::stringstream ss;
			ss << "Fatal error: cannot initialize logging semaphore: " << strerror(errno);
			throw std::runtime_error(ss.str());
		}
		for (int i=0; i< buffer_size; ++i) {
			_buffer[i] = 0;
		}
		_buffer_cursor = 0;
		_sticky_force = true;
		_use_stderr = true;
	}

public:

	/**
	 * @brief Get the singleton logger object
	 */
	static SimpleLogger& get(bool use_stderr = true) {
		if (_instance == nullptr) {
			_instance = new SimpleLogger();
		}
		_instance->_use_stderr = use_stderr;
		return *_instance;
	}

	/**
	 * @brief Adds new message to buffer
	 * @param message The message after formatting
	 * @param flush True iff the message should be flushed to stderr
	 */
	void add(const char* message, bool flush = false) {
		// Should I print message to screen
		if ((_sticky_force || flush) && _use_stderr){
			std::cerr << message << std::flush;
		}
		else if (!_use_stderr){
			std::cout << message;
		}

		// Copy message to buffer
		int msg_size = strlen(message);
		if ((_buffer_cursor+msg_size+1)>=buffer_size) {
			_buffer_cursor=0;
		}
		strcpy(&_buffer[_buffer_cursor], message);
	}

	/**
	 * @brief A logger command. Close the current line.
	 */
	static SimpleLogger& endl() {
		SimpleLogger::get().add("\n");
		return *_instance;
	}

	/**
	 * @brief A logger command. Acquire logging lock.
	 */
	static SimpleLogger& lock() {
		// Acquire log semaphore
		int s;
		while ( (s = sem_wait(&SimpleLogger::get()._semaphore) == -1) && (errno==EINTR) ) continue;
		if (s==-1) {
			std::stringstream ss;
			ss <<  "Fatal error while waiting for semaphore: " << strerror(errno);
			throw std::runtime_error(ss.str());
		}
		return *_instance;
	}

	/**
	 * @brief A logger command. Release logging lock.
	 */
	static SimpleLogger& release() {
		sem_post(&SimpleLogger::get()._semaphore);
		return *_instance;
	}

	/**
	 * @brief Adds nothing to the buffer.
	 *        Used to execute log commands.
	 */
	SimpleLogger& operator<<(const SimpleLogger& rhs) {
		return *this;
	}

	/**
	 * @brief Log arbitrary items
	 */
	template <typename T>
	SimpleLogger& operator<<(const T& rhs) {
		std::stringstream ss;
		ss << rhs;
		add(ss.str().c_str(), 0);
		return *this;
	}

	/**
	 * @brief A logger command. Adds message to log as a formatted string.
	 */
	static std::string format(const char* fmt, ...) {
		std::va_list args;
		va_start(args, fmt);
		size_t size = vsnprintf( nullptr, 0, fmt, args) + 1;
		char buffer[size];
		va_start(args, fmt);
		vsnprintf(buffer, size, fmt, args);
		return std::string(buffer);
	}

	/**
	 * @brief Returns the buffer as string
	 */
	const char* get_buffer() const { return this->_buffer; }

	/**
	 * @brief In case sticky force is set, every logged
	 *        message will be printed to stderr.
	 *        True by default
	 */
	void set_sticky_force(bool value) {
		_sticky_force = value;
	}

};

/**
 * @brief A general error class for constructing informative exceptions
 *        while updating log. Usage should be via MACRO error()
 */
class SimpleException : public std::exception {

	std::stringstream _buffer;
	std::string _message;
	SimpleException() {}

public:

	/**
	 * @brief Used by the throw mechanism
	 */
	SimpleException(const SimpleException& rhs) {
		this->_message = rhs._message;
	}

	/**
	 * @brief Creates new exception class
	 */
	static SimpleException create() {
		return SimpleException();
	}

	/**
	 * @brief Adds nothing to the message.
	 *        Used to execute log commands.
	 */
	SimpleException& operator<<(const SimpleLogger& rhs) {
		return *this;
	}

	/**
	 * @brief Append to message arbitrary info
	 */
	template <typename T>
	SimpleException& operator<<(const T& rhs) {
		_buffer << rhs;
		_message = _buffer.str();
		return *this;
	}

	/**
	 * @brief A log command. Adds formatted message.
	 */
	static std::string format(const char* fmt, ...) {
		std::va_list args;
		va_start(args, fmt);
		size_t size = vsnprintf( nullptr, 0, fmt, args) + 1;
		char buffer[size];
		va_start(args, fmt);
		vsnprintf(buffer, size, fmt, args);
		return std::string(buffer);
	}

	virtual const char* what() const noexcept {
		return _message.c_str();
	}
};
