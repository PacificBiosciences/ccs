// Author: David Seifert

// Reduce the number of exposed symbols in order to speed up
// DSO load times
// https://gcc.gnu.org/wiki/Visibility

#pragma once

#if defined _WIN32 || defined __CYGWIN__
#ifdef UNANIMITY_BUILDING_LIBRARY
#ifdef __GNUC__
#define UNANIMITY_PUBLIC_API __attribute__((dllexport))
#else
#define UNANIMITY_PUBLIC_API __declspec(dllexport)  // Note: gcc seems to also supports this syntax
#endif
#else
#ifdef __GNUC__
#define UNANIMITY_PUBLIC_API __attribute__((dllimport))
#else
#define UNANIMITY_PUBLIC_API __declspec(dllimport)  // Note: gcc seems to also supports this syntax
#endif
#endif
#define UNANIMITY_PRIVATE_API
#else
#if __GNUC__ >= 4
#define UNANIMITY_PUBLIC_API __attribute__((visibility("default")))
#define UNANIMITY_PRIVATE_API __attribute__((visibility("hidden")))
#else
#define UNANIMITY_PUBLIC_API
#define UNANIMITY_PRIVATE_API
#endif
#endif

// Compatibility guard to prevent old compilers
// choking on relaxed C++14 constexpr handling
// https://github.com/PacificBiosciences/pitchfork/issues/417

#if (__GNUC__ >= 6) || defined(__clang__)
#define UNANIMITY_CONSTEXPR constexpr
#else
#define UNANIMITY_CONSTEXPR
#endif
