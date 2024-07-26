#ifndef LION_IO_DOCUMENT_ELEMENT_H
#define LION_IO_DOCUMENT_ELEMENT_H

#include "lion/foundation/utils.hpp"
#include "lion/math/vector3d.hpp"
#include "lion/math/matrix3x3.h"
#include "lion/math/matrix_extensions.h"


class Document_element
{
public:

    class value_ptr : public std::shared_ptr<Document_element>
    {
    public:
        using base_type = std::shared_ptr<Document_element>;

        template<typename ...Args>
        value_ptr(Args&& ... args) : base_type(std::forward<Args>(args)...) {}


        template<typename T>
        inline T& cast()
        {
            auto cast_ptr = dynamic_cast<T*>(&base_type::operator*());

            if (cast_ptr != nullptr) {
                return *cast_ptr;
            }
            else {
                throw lion_exception("[ERROR] Invalid cast: cannot cast 'Document_element' to requested type");
            }
        }

        template<typename T>
        inline bool can_cast()
        {
            return dynamic_cast<T*>(&base_type::operator*()) != nullptr ? true : false;
        }

    private:

        Document_element& operator*() 
        {
            return base_type::operator*();
        }
    };

    Document_element(void* e, const std::string& name) : _name(name), _e(e)
    {
        if (e == nullptr) {
            throw lion_exception("Attempted to construct Document_element from nullptr");
        }
    }

    const auto& get_name() const { return _name; }

    const void* e() const { return _e; }
    void* e() { return _e; }


    virtual std::string get_value() const = 0;

    std::string get_value(std::string&&) const
    {
        return get_value();
    }

    virtual int get_value(int&&) const = 0;

    virtual double get_value(double&&) const = 0;

    virtual std::vector<double> get_value(std::vector<double>&&) const = 0;

    virtual std::vector<float> get_value(std::vector<float>&&) const = 0;

    virtual sVector3d get_value(sVector3d&&) const = 0;

    virtual sMatrix3x3 get_value(sMatrix3x3&&) const = 0;

    virtual void set_value(const std::string& val) = 0;

    virtual void set_value(const char* val) = 0;

    virtual void set_value(const double& val) = 0;

    virtual void set_value(const int& val) = 0;

    virtual void set_value(const float& val) = 0;

    virtual void set_value(const std::vector<double>& val) = 0;

    virtual void set_value(const std::vector<float>& val) = 0;

    virtual void set_value(const std::vector<int>& val) = 0;

    virtual void set_value(const sVector3d& val) = 0;

    virtual void set_value(const sMatrix3x3& val) = 0;

    virtual std::vector<value_ptr> get_children() = 0;

    virtual value_ptr add_child(const std::string& name) = 0;

    virtual bool has_children() const = 0;

    virtual value_ptr get_child(const std::string& name) = 0;

    virtual void print(std::ostream& os) const = 0;

    virtual bool has_child(const std::string& name) = 0;

    virtual void copy_contents(Document_element::value_ptr other) = 0;

    // ! Get children & their values and emplace_back them in an "std::vector<std::pair<std::string, number_or_vector_or_array> >"
    //! @param[inout] vp: the vector of pairs holding numbers or std::vectors or std::arrays, whose keys are std::strings
    template<typename T>
    void emplace_back_children_and_values_in_vector_of_pairs(std::vector<std::pair<std::string, T> > &vp);

    // ! Set children & their values from an "std::(unordered)map<std::string, number_or_vector_or_array>"
    //! @param[in] m: the map holding numbers or std::vectors or std::arrays, whose keys are std::strings
    template<typename MapType>
    void add_children_and_values_from_map(const MapType &m);

    // ! Get children & their values and emplace them in an "std::(unordered)map<std::string, number_or_vector_or_array>"
    //! @param[inout] m: the map holding numbers or std::vectors or std::arrays, whose keys are std::strings
    template<typename MapType>
    void emplace_children_and_values_in_map(MapType &m);


    //! Check if this element and its childs have any attribute
    virtual bool this_and_childs_have_attributes() = 0;

    virtual bool this_and_childs_have_both_value_and_children() = 0;


private:

    std::string _name;
    void* _e;
};

using Document_element_ptr = Document_element::value_ptr;


inline std::ostream& operator<<(std::ostream& os, Document_element& doc)
{
    doc.print(os);
    return os;
}


template<typename T>
inline void Document_element::emplace_back_children_and_values_in_vector_of_pairs(std::vector<std::pair<std::string, T> > &vp)
{
    //
    // Emplaces the children and values of *this Xml_element
    // in the given std::vector of std::pairs. These pairs must be
    // "std::pair<std::string, T>" in which "T" can be either
    // "std::vector", "std::array" or a number.
    //

    constexpr auto pairs_of_vectors = is_specialization_v<T, std::vector>;
    constexpr auto pairs_of_arrays = is_std_array_v<T>;
    constexpr auto pairs_of_scalars = std::is_scalar_v<T>;

    static_assert(pairs_of_vectors || pairs_of_arrays || pairs_of_scalars,
                  "Xml_element::emplace_back_children_and_values_in_vector_of_pairs:: unsupported value type.");

    for (auto& c : get_children()) {
        vp.emplace_back(c->get_name(), c->get_value(T{}));
    }
}


template<typename MapType>
inline void Document_element::emplace_children_and_values_in_map(MapType &m)
{
    //
    // Emplaces the children and values of *this Xml_element
    // in the given map. This map must be an
    // "std::(unordered_)map<std::string, mapped_type>" in which
    // "mapped_type" can be either "std::vector", "std::array" or a number.
    // If a child cannot be emplaced, the function throws an
    // std::runtime_error.
    //

    using key_type = typename MapType::key_type;
    using mapped_type = typename MapType::mapped_type;

    constexpr auto map_key_is_string = std::is_same_v<key_type, std::string>;
    constexpr auto map_of_vectors = is_specialization_v<mapped_type, std::vector>;
    constexpr auto map_of_arrays = is_std_array_v<mapped_type>;
    constexpr auto map_of_scalars = std::is_scalar_v<mapped_type>;

    static_assert(map_key_is_string && (map_of_vectors || map_of_arrays || map_of_scalars),
                  "Xml_element::emplace_children_and_values_in_map:: unsupported map type.");

    for (auto& c : get_children()) {
        const auto [_, did_emplace] = m.emplace(c->get_name(), c->get_value(mapped_type{}));

        if (!did_emplace) {
            throw std::runtime_error(
                std::string{ "Xml_element::emplace_children_and_values_in_map:"
                             " could not emplace child \"" } +
                c->get_name() + "\" in the given map.");
        }
    }
}



template<typename MapType>
inline void Document_element::add_children_and_values_from_map(const MapType &m)
{
    //
    // Saves an "std::(unordered_)map<std::string, mapped_type>" into the Xml_element,
    // in which "mapped_type" can be either "std::vector", "std::array" or a number.
    //

    using key_type = typename MapType::key_type;
    using mapped_type = typename MapType::mapped_type;

    constexpr auto map_key_is_string = std::is_same_v<key_type, std::string>;
    constexpr auto map_of_vectors = is_specialization_v<mapped_type, std::vector>;
    constexpr auto map_of_arrays = is_std_array_v<mapped_type>;
    constexpr auto map_of_scalars = std::is_scalar_v<mapped_type>;

    static_assert(map_key_is_string && (map_of_vectors || map_of_arrays || map_of_scalars),
                  "Xml_element::add_children_and_values_from_map:: unsupported map type.");

    for (const auto &mi : m) {
        add_child(mi.first)->set_value(mi.second);
    }
}



#endif
