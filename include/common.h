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

typedef unsigned long uint64_t;
typedef unsigned int  uint32_t;
typedef unsigned char uint8_t;

// Used for micro-debugging
#define MICRO_DEBUG fprintf(stderr, "Checkpoint at %s:%d \n", __FILE__, __LINE__);

/**
 * @brief Macro for fast loading from memory and perform size checks
 * @param T Typename to read from buffer
 * @param buff The buffer to read from
 * @param size Size argument for buffer overflow check.
 * @param label An error label to jump in case of an error
 */
#define safe_buffer_read(T, buf, size, label)                \
    *(T*)buf;                                                \
    if ((size-=sizeof(T))<0) {                               \
        error("error while reading input");                  \
        goto label;                                          \
    }                                                        \
    buf=(uint8_t*)((T*)buf+1);


/**
 * @brief Macro for fast writing tomemory and perform size checks
 * @param T Typename to read from buffer
 * @param buff The buffer to write to
 * @param size Size argument for buffer overflow check.
 * @param label An error label to jump in case of an error
 */
#define safe_buffer_write(T, buf, size, label)             \
    *(T*)buf;                                              \
    if ((size-=sizeof(T))<0) {                             \
        error("error while writing buffer");               \
        goto label;                                        \
    }                                                      \
    buf=(T*)buf+1;


/**
 * @brief Print messages to stderr
 */
#ifndef buffer_print
#    define buffer_print(type, ...) {                                                \
        char _msg[2048];                                                             \
        int _end=0;                                                                  \
        _end+=snprintf(&_msg[_end], 2048-_end, "%s: (%s) ", type, __func__);         \
        _end+=snprintf(&_msg[_end], 2048-_end, __VA_ARGS__);                         \
        _end+=snprintf(&_msg[_end], 2048-_end, " (%s:%d)\n", __FILE__, __LINE__);    \
        fprintf(stderr, "%s", _msg);                                                 \
    }
#endif

/**
 * @brief Pints error message to stderr and fail
 */
#ifndef fail
#    define fail(...) buffer_print("Error", __VA_ARGS__); exit(1);
#endif

/**
 * @brief Pints error message to stderr
 */
#ifndef error
#    define error(...) buffer_print("Error", __VA_ARGS__)
#endif

/**
 * @brief Pints warning message to stderr
 */
#ifndef warning
#    define warning(...) buffer_print("Warning", __VA_ARGS__)
#endif

/**
 * @brief Pints regular message to stdout
 */
#ifndef message
#    define message(...) fprintf(stdout, __VA_ARGS__); fprintf(stdout, "\n");
#endif

/**
 * @brief Logs in case DEBUG option is on
 */
#ifndef info
#    ifndef NDEBUG
#        define info(...) buffer_print("Log", __VA_ARGS__)
#    else /* NDEBUG */
#        define info(...)
#    endif
#endif

