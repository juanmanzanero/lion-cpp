#ifndef LION_FOUNDATION_TUPLE_HELPERS_H
#define LION_FOUNDATION_TUPLE_HELPERS_H


#include <tuple>
#include <array>
#include <algorithm>


//
// Defines handy functions and type traits to operate on
// std::tuples (most should also work on std::arrays).
//


#ifndef _MSC_VER

namespace lioncpp::detail::tuple_helpers {

template<typename IndexType, IndexType Offset, typename IntegerSequence>
struct offset_integer_sequence {};

template<typename IndexType, IndexType Offset, IndexType... Is>
struct offset_integer_sequence<IndexType, Offset, std::integer_sequence<IndexType, Is...> >
{
    using type = std::integer_sequence<IndexType, Is + Offset...>;
};

template<typename IndexType, IndexType Offset, typename IntegerSequence>
using offset_integer_sequence_t =
    typename offset_integer_sequence<IndexType, Offset, IntegerSequence>::type;


template<typename Fun, typename IndexType, IndexType... Is>
constexpr void index_for_impl(Fun &&f, std::integer_sequence<IndexType, Is ...>)
{
    (f(std::integral_constant<IndexType, Is>{}), ...);
}


template<typename Tuple, typename Fun, typename IndexType, IndexType... Is>
constexpr void tuple_for_impl(Tuple &&tuple, Fun &&f, std::integer_sequence<IndexType, Is...>)
{
    (f(std::get<Is>(std::forward<Tuple>(tuple))), ...);
}


template<typename Tuple, typename Fun, typename IndexType, IndexType... Is>
constexpr void tuple_index_for_impl(Tuple &&tuple, Fun &&f, std::integer_sequence<IndexType, Is...>)
{
    (f(std::get<Is>(std::forward<Tuple>(tuple)), std::integral_constant<IndexType, Is>{}), ...);
}

} // end namespace lioncpp::detail::tuple_helpers


//
// Loops the range "[Begin, End)", applying a function with
// signature "void f(auto index)" at each iteration.
//

template<std::size_t Begin, std::size_t End,
    typename Fun,
    typename std::enable_if_t<Begin < End>* = nullptr>
constexpr void index_for(Fun &&f)
{
    lioncpp::detail::tuple_helpers::index_for_impl(
        std::forward<Fun>(f),
        lioncpp::detail::tuple_helpers::offset_integer_sequence_t<std::size_t, Begin,
            std::make_integer_sequence<std::size_t, End - Begin> >{});
}

template<std::size_t Begin, std::size_t End,
    typename Fun,
    typename std::enable_if_t<Begin >= End>* = nullptr>
constexpr void index_for(Fun &&) {}


//
// Loops the range "[0, End)", applying a function with
// signature "void f(auto index)" at each iteration.
//

template<std::size_t End,
    typename Fun>
constexpr void index_for(Fun &&f)
{
    lioncpp::detail::tuple_helpers::index_for_impl(
        std::forward<Fun>(f),
        std::make_integer_sequence<std::size_t, End>{});
}


//
// Traverses an std::tuple or an std::array, applying a
// function with signature "void f(auto &&tuple_element)"
// on each of its elements.
//

template<typename Tuple, typename Fun>
constexpr void tuple_for(Tuple &&tuple, Fun &&f)
{
    lioncpp::detail::tuple_helpers::tuple_for_impl(
        std::forward<Tuple>(tuple),
        std::forward<Fun>(f),
        std::make_integer_sequence<std::size_t, std::tuple_size_v<std::decay_t<Tuple> > >{});
}


//
// Traverses an std::tuple or an std::array, applying a function
// with signature "void f(auto &&tuple_element, auto element_index)"
// on each of its elements.
//

template<typename Tuple, typename Fun>
constexpr void tuple_index_for(Tuple &&tuple, Fun &&f)
{
    lioncpp::detail::tuple_helpers::tuple_index_for_impl(
        std::forward<Tuple>(tuple),
        std::forward<Fun>(f),
        std::make_integer_sequence<std::size_t, std::tuple_size_v<std::decay_t<Tuple> > >{});
}

#else
// VS2017 hasn't got fold expressions... we'll implement
// the loops using compile-time recursion


//
// Loops the range "[Begin, End)", applying a function with
// signature "void f(auto index)" at each iteration.
//

template<std::size_t Begin, std::size_t End,
    typename Fun,
    typename std::enable_if_t<Begin >= End>* = nullptr>
constexpr void index_for(Fun &&) {}

template<std::size_t Begin, std::size_t End,
    typename Fun,
    typename std::enable_if_t<Begin < End>* = nullptr>
constexpr void index_for(Fun &&f)
{
    f(std::integral_constant<std::size_t, Begin>{});
    index_for<Begin + 1, End>(std::forward<Fun>(f));
}


//
// Loops the range "[0, End)", applying a function with
// signature "void f(auto index)" at each iteration.
//

template<std::size_t End,
         typename Fun>
constexpr void index_for(Fun &&f)
{
    index_for<0u, End>(std::forward<Fun>(f));
}


//
// Traverses an std::tuple or an std::array, applying a
// function with signature "void f(auto &&tuple_element)"
// on each of its elements.
//

template<typename Tuple, typename Fun,
    std::size_t Begin = 0u,
    typename std::enable_if_t<Begin >= std::tuple_size_v<std::decay_t<Tuple> > >* = nullptr>
constexpr void tuple_for(Tuple &&, Fun &&) {}

template<typename Tuple, typename Fun,
    std::size_t Begin = 0u,
    typename std::enable_if_t<Begin < std::tuple_size_v<std::decay_t<Tuple> > >* = nullptr>
constexpr void tuple_for(Tuple &&tuple, Fun &&f)
{
    f(std::get<Begin>(tuple));
    tuple_for<Tuple, Fun, Begin + 1u>(std::forward<Tuple>(tuple),
        std::forward<Fun>(f));
}


//
// Traverses an std::tuple or an std::array, applying a function
// with signature "void f(auto &&tuple_element, auto element_index)"
// on each of its elements.
//

template<typename Tuple, typename Fun,
    std::size_t Begin = 0u,
    typename std::enable_if_t<Begin >= std::tuple_size_v<std::decay_t<Tuple> > >* = nullptr>
constexpr void tuple_index_for(Tuple &&, Fun &&) {}

template<typename Tuple, typename Fun,
    std::size_t Begin = 0u,
    typename std::enable_if_t<Begin < std::tuple_size_v<std::decay_t<Tuple> > >* = nullptr>
constexpr void tuple_index_for(Tuple &&tuple, Fun &&f)
{
    f(std::get<Begin>(tuple), std::integral_constant<std::size_t, Begin>{});
    tuple_index_for<Tuple, Fun, Begin + 1>(
        std::forward<Tuple>(tuple), std::forward<Fun>(f));
}

#endif


//
// Performs std::decay on every member type of an std::tuple.
//

template<typename ... Ts>
constexpr auto decay_types_of_a_tuple(const std::tuple<Ts...> &) ->
    std::tuple<std::remove_cv_t<std::remove_reference_t<Ts> >... >;

template<typename Tuple>
using decay_tuple = decltype(decay_types_of_a_tuple(std::declval<Tuple>()));


//
// A type_trait to identify if a tuple type contains a given type.
//

template<typename T, typename Tuple>
struct tuple_has_type : std::false_type {};

template<typename T>
struct tuple_has_type<T, std::tuple<> > : std::false_type {};

template<typename T, typename... Ts>
struct tuple_has_type<T, std::tuple<T, Ts...> > : std::true_type {};

template<typename T, typename U, typename... Ts>
struct tuple_has_type<T, std::tuple<U, Ts...> > : tuple_has_type<T, std::tuple<Ts...> > {};

template<typename T, typename Tuple>
constexpr bool tuple_has_type_v = tuple_has_type<T, Tuple>::type::value;


//
// Concatenation.
//

#ifndef _MSC_VER
template<typename Type, std::size_t... sizes>
constexpr auto array_cat(const std::array<Type, sizes>&... arrays)
{
    //
    // Concatenates the input std::arrays into a single one.
    //

    std::array<Type, (sizes + ...)> result;
    std::size_t index = 0u;

    ((std::copy_n(arrays.cbegin(), sizes, std::next(result.begin(), index)), index += sizes), ...);

    return result;
}

#else

namespace lioncpp::detail::tuple_helpers {

template<typename... Args>
struct add_array_sizes {};

template<typename T, std::size_t size>
struct add_array_sizes<std::array<T, size> >
{
    static constexpr auto value = size;
};

template<typename T, std::size_t size, typename... RestOfArrays>
struct add_array_sizes<std::array<T, size>, RestOfArrays...>
{
    static constexpr std::size_t value = size + add_array_sizes<RestOfArrays...>::value;
};

} // end namespace lioncpp::detail::tuple_helpers

template<typename Type, std::size_t... sizes>
constexpr auto array_cat(const std::array<Type, sizes>&... arrays)
{
    //
    // Concatenates the input std::arrays into a single one
    // (implementation that works in VS2017).
    //

    std::array<Type, lioncpp::detail::tuple_helpers::add_array_sizes<std::array<Type, sizes>...>::value> result;
    std::size_t index = 0u;

    tuple_for(std::forward_as_tuple(arrays...),
        [&index, &result](const auto &arr)
        {
            constexpr auto size = std::tuple_size_v<std::decay_t<decltype(arr)> >;
            std::copy_n(arr.cbegin(), size, std::next(result.begin(), index));
            index += size;
        });

    return result;
}

#endif

#endif
