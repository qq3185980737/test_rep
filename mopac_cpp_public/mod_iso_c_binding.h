// mod_iso_c_binding.h — C++ translation.
#pragma once
#include <cstdint>
namespace iso_c_binding {
    using c_intptr_t = intptr_t;
    constexpr int c_char = 1;
    constexpr int c_bool = 1;
    constexpr int c_int = 4;
    constexpr int c_float = 4;
    constexpr int c_long = 4;
    constexpr int c_double = 8;
    using c_size_t = size_t;
    struct c_ptr { c_intptr_t ptr; };
}
