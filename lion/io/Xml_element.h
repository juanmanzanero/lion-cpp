#ifndef __XML_ELEMENT_H__
#define __XML_ELEMENT_H__

#include "tinyxml2.h"

#include "lion/thirdparty/include/logger.hpp"
#include "lion/foundation/utils.hpp"
#include "lion/math/vector3d.hpp"
#include "lion/math/matrix3x3.h"
#include "lion/math/matrix_extensions.h"

class Xml_element
{
 public:
    //! Constructor
    //! @param[in] e: pointer to a tinyxml2 XMLElement
    Xml_element(tinyxml2::XMLElement* e); 

    //! Get the element name
    std::string get_name() const { return _e->Name(); }

    //! Set the element name
    //! @param[in] name: new name for the element
    void set_name(const std::string& name) { _e->SetName(name.c_str()); }

    //! Get value as string
    std::string get_value()             const { return (_e->GetText() == nullptr) ? "" : _e->GetText(); }
    std::string get_value(std::string&&) const { return get_value(); }

    //! Get value as double
    //! @param[in] simply pass int() to overload this version
    int get_value(int&&) { return std::stoi(get_value()); }

    //! Get value as double
    //! @param[in] simply pass double() to overload this version
    double get_value(double&&) { return std::stod(get_value()); }

    //! Get value as double vector
    //! @param[in] simply pass std::vector<double>() to overload this version
    std::vector<double> get_value(std::vector<double>&&) { return string_to_double_vector<double>(get_value()); }

    //! Get value as float vector
    //! @param[in] simply pass std::vector<float>() to overload this version
    std::vector<float> get_value(std::vector<float>&&) { return string_to_double_vector<float>(get_value()); }

    //! Get value as vector3d
    //! @param[in] simply pass Vector3d() to overload this function
    sVector3d get_value(sVector3d&&) { return sVector3d(string_to_double_vector<double>(get_value())); }

    //! Get value as matrix3x3
    //! @param[in] simply pass Matrix3x3() to overload this function
    sMatrix3x3 get_value(sMatrix3x3&&) { return transpose(sMatrix3x3(string_to_double_vector<double>(get_value()))); }

    //! Get value as bool
    //! @param[in] simply pass bool() to overload this function
    bool get_value(bool&&) { return to_bool(get_value()); }

    //! Get value as an std::array<double, N>
    //! @param[in] simply pass and std::array<double, N>{} to overload this function
    template<std::size_t N>
    std::array<double, N> get_value(std::array<double, N>&&)
    {
        const auto vec = get_value(std::vector<double>{});
        if (vec.size() != N) {
            throw std::runtime_error("Xml_element::get_value(std::array<double, N>): incorrect size.");
        }
        std::array<double, N> arr{};
        std::copy_n(vec.cbegin(), N, arr.begin());
        return arr;
    }

    //! Set the value from string
    //! @param[in] val: new string value
    Xml_element& set_value(const std::string& val) { _e->SetText(val.c_str()); return *this; }


    //! Set the value from numbers (scalars, std::vector, std::array)
    template<typename T,
             typename std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
    Xml_element& set_value(const T &num) { return set_value(num2str(num)); }

    template<typename T>
    Xml_element& set_value(const std::vector<T> &vec) { return set_value(vec2str(vec)); }

    template<typename T, std::size_t N>
    Xml_element& set_value(const std::array<T, N> &arr) { return set_value(vec2str(arr)); }


    // ! Set children & their values from an "std::(unordered)map<std::string, number_or_vector_or_array>"
    //! @param[in] m: the map holding numbers or std::vectors or std::arrays, whose keys are std::strings
    template<typename MapType>
    void add_children_and_values_from_map(const MapType &m);

    // ! Get children & their values and emplace them in an "std::(unordered)map<std::string, number_or_vector_or_array>"
    //! @param[inout] m: the map holding numbers or std::vectors or std::arrays, whose keys are std::strings
    template<typename MapType>
    void emplace_children_and_values_in_map(MapType &m);

    // ! Get children & their values and emplace_back them in an "std::vector<std::pair<std::string, number_or_vector_or_array> >"
    //! @param[inout] vp: the vector of pairs holding numbers or std::vectors or std::arrays, whose keys are std::strings
    template<typename T>
    void emplace_back_children_and_values_in_vector_of_pairs(std::vector<std::pair<std::string, T> > &vp);


    //! Get all children to a std::vector
    std::vector<Xml_element> get_children() const;

    //! Add children
    Xml_element add_child(const std::string& name); 

    //! Add comment
    void add_comment(const std::string& text) { _e->InsertNewComment(text.c_str()); }

    //! See if the element has children
    bool has_children() const {return (_e->FirstChildElement() ? true : false);}

    //! Get the first child
    Xml_element get_first_child() const {return _e->FirstChildElement();}

    //! Get child by name
    Xml_element get_child(const std::string& name) const;

    //! See if a child exists
    bool has_child(const std::string& name) const;


    //! Try to get the value of a child, or return a default one
    template<typename ValueType>
    ValueType try_get_child_value(const std::string &childname, ValueType default_value,
                                  bool warn_if_returning_default_value = true) const
    {
        if (has_child(childname)) {
            return get_child(childname).get_value(ValueType{});
        }
        else {
            if (warn_if_returning_default_value) {
                std::cerr << "Xml_element::try_get_child_value: warning, xml element \""
                    << get_name()
                    << "\" does not contain child \""
                    << childname << "\", returning a default value of \""
                    << default_value
                    << "\"."
                    << std::endl;
            }
            return default_value;
        }
    }


    //! Try to get the value of an attribute, or return a default one
    template<typename ValueType>
    ValueType try_get_attribute(const std::string &attributename, ValueType default_value,
                                  bool warn_if_returning_default_value = true) const
    {
        if (has_attribute(attributename)) {
            return get_attribute(attributename, ValueType{});
        }
        else {
            if (warn_if_returning_default_value) {
                std::cerr << "Xml_element::try_get_attribute_value: warning, xml element \""
                    << get_name()
                    << "\" does not contain attribute \""
                    << attributename << "\", returning a default value of \""
                    << default_value
                    << "\"."
                    << std::endl;
            }
            return default_value;
        }
    }


    //! See if parent exists
    bool has_parent() const { return (_e->Parent()->ToElement() ? true : false); }

    //! Get parent
    Xml_element get_parent() const { return _e->Parent()->ToElement(); }

    //! See if sibling exists
    bool has_sibling() const { return (_e->NextSiblingElement() ? true : false); }

    //! Get Next Sibling
    Xml_element get_sibling() const { return _e->NextSiblingElement()->ToElement(); }

    //! See if an attribute exists
    //! @param[in] attribute: name of the attribute
    bool has_attribute(const std::string& attribute) const {return (_e->FindAttribute(attribute.c_str()) ? true : false);}

    //! Get attribute by name
    //! @param[in] attribute: name of the attribute
    std::string get_attribute(const std::string& attribute, std::string&& = "") const;

    //! Get attribute as double
    //! @param[in] attribute: name of the attribute
    //! @param[in] simply pass double() to overload this version
    double get_attribute(const std::string& attribute, double&&) const { return std::stod(_e->Attribute(attribute.c_str())); }

    //! Sets a new attribute or modifies the value of an existing one
    void set_attribute(const std::string& attrib_name, const std::string& attrib_value) { _e->SetAttribute(attrib_name.c_str(), attrib_value.c_str()); }

    template<typename T>
    void set_attribute(const std::string& attrib_name, const T& attrib_value) { _e->SetAttribute(attrib_name.c_str(), attrib_value); }

    //! Print
    void print(std::ostream& os) const;
    
    //! Copy the contents of other into this node
    void copy_contents(Xml_element other);

    //! Delete attribute
    void delete_attribute(const std::string& attribute_name)
    {

        if (has_attribute(attribute_name))
            _e->DeleteAttribute(attribute_name.c_str());
        else
            throw lion_exception("[ERROR] delete_attribute -> attribute does not exist");
    } 


    //! Check if this element and its childs have any attribute
    bool this_and_childs_have_attributes() const;


    //! Check if this element and its childs have value and children
    bool this_and_childs_have_both_value_and_children() const;


 private:
    tinyxml2::XMLElement* _e;
};


inline void Xml_element::print(std::ostream& os) const
{
    tinyxml2::XMLPrinter printer;
    _e->Accept(&printer);
    os << printer.CStr();
}


inline std::ostream& operator<<(std::ostream& os, const Xml_element& e)
{
    e.print(os);
    return os;
}


inline void Xml_element::copy_contents(Xml_element other)
{
    tinyxml2::XMLElement* element = other._e->FirstChildElement(); 
    tinyxml2::XMLNode* previous_sibling_ptr = _e->LastChild();
    while (element != nullptr)
    {
        auto* element_copy = element->DeepClone(_e->GetDocument()); 

        if (previous_sibling_ptr != nullptr)
            _e->InsertAfterChild(previous_sibling_ptr, element_copy);
        else
            _e->InsertFirstChild(element_copy);

        element = element->NextSiblingElement();
        previous_sibling_ptr = element_copy;
    }
}


inline Xml_element::Xml_element(tinyxml2::XMLElement* e) : _e(e) 
{
    if ( e == nullptr )
        throw lion_exception("Attempted to construct Xml_element from nullptr"); 
} 


inline std::vector<Xml_element> Xml_element::get_children() const
{
    std::vector<Xml_element> output;

    tinyxml2::XMLElement* element = _e->FirstChildElement();

    while (element != nullptr)
    {
        output.push_back({element});
        element = element->NextSiblingElement();
    }

    return output;
}


inline Xml_element Xml_element::get_child(const std::string& name) const
{
    // Divide the string in child_name/the_rest
    std::string::size_type pos = name.find('/');
    std::string child_name;
    std::string the_rest;

    if (pos != std::string::npos)
    {
        child_name = name.substr(0, pos);
        the_rest = name.substr(pos+1);
    }
    else
    {
        child_name = name;
    }

    // Look for child_name, and make sure it is the only occurrence
    Xml_element child(_e);

    if ( child_name == "." )
        child = _e;

    else if ( child_name == ".." )
        child = get_parent();

    else
    {
        if ( _e->FirstChildElement(child_name.c_str()) == nullptr )
        {
            std::ostringstream s_out;
            s_out << "[ERROR] Xml_element::get_child -> No child found with name \"" + child_name + "\"" << std::endl;
            print(s_out);
            throw lion_exception(s_out.str());
        }

        child = _e->FirstChildElement(child_name.c_str());
    }

    if (child._e->NextSiblingElement(child_name.c_str()) != nullptr ) 
    {
        std::ostringstream s_out;
        s_out << "[ERROR] Xml_element::get_child -> There are several children with name \"" + child_name + "\"" << std::endl;
        print(s_out);
        throw lion_exception("There are several children with name \"" + child_name + "\"");
    }

    if ( the_rest.size() == 0 )
        return child;

    else
        return child.get_child(the_rest);
}


inline Xml_element Xml_element::add_child(const std::string& name) 
{
    // Divide the string in child_name/the_rest
    std::string::size_type pos = name.find('/');
    std::string child_name;
    std::string the_rest;

    if (pos != std::string::npos)
    {
        child_name = name.substr(0, pos);
        the_rest = name.substr(pos+1);
    }
    else
    {
        child_name = name;
    }

    // Look for child_name: if exists, make sure there's only one occurrence. If it does not, create it
    Xml_element child(_e);

    if ( the_rest.size() > 0 )
    {
        // The child exists and its me
        if ( child_name == "." )
            child = _e;

        // The child exists and its my parent
        else if ( child_name == ".." )
            child = get_parent();

        else
        {
            // The children does not exist -> create it
            if ( _e->FirstChildElement(child_name.c_str()) == nullptr )
            {
                child = _e->InsertNewChildElement(child_name.c_str());
            }
            else
            {
                // The children does exist -> point to it
                child = _e->FirstChildElement(child_name.c_str());
        
                // Check uniqueness
                if (child._e->NextSiblingElement(child_name.c_str()) != nullptr ) 
                {
                    std::ostringstream s_out;
                    s_out << "[ERROR] Xml_element::add_child -> There are several children with name \"" + child_name + "\"" << std::endl;
                    print(s_out);
                    throw lion_exception(s_out.str());
                }
            }
        }

        // Call add child
        return child.add_child(the_rest);
    } 
    else
    {
        // Make sure that the node does not already exist
        if ( _e->FirstChildElement(child_name.c_str()) != nullptr )
        {
            std::ostringstream s_out;
            s_out << "[ERROR] Xml_element::add_child -> Node \"" + child_name + "\" already exists" << std::endl;
            print(s_out);
            throw lion_exception("Node \"" + child_name + "\" already exists");
        }

        return {_e->InsertNewChildElement(name.c_str())};
    }
}


template<typename MapType>
void Xml_element::add_children_and_values_from_map(const MapType &m)
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
        add_child(mi.first).set_value(mi.second);
    }
}


template<typename MapType>
void Xml_element::emplace_children_and_values_in_map(MapType &m)
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

    for (auto c : get_children()) {
        const auto [_, did_emplace] = m.emplace(c.get_name(), c.get_value(mapped_type{}));

        if (!did_emplace) {
            throw std::runtime_error(
                std::string{ "Xml_element::emplace_children_and_values_in_map:"
                             " could not emplace child \"" } +
                c.get_name() + "\" in the given map.");
        }
    }
}


template<typename T>
void Xml_element::emplace_back_children_and_values_in_vector_of_pairs(std::vector<std::pair<std::string, T> > &vp)
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

    for (auto c : get_children()) {
        vp.emplace_back(c.get_name(), c.get_value(T{}));
    }
}


inline bool Xml_element::has_child(const std::string& name) const
{
    // Divide the string in child_name/the_rest
    std::string::size_type pos = name.find('/');
    std::string child_name;
    std::string the_rest;

    if (pos != std::string::npos)
    {
        child_name = name.substr(0, pos);
        the_rest = name.substr(pos+1);
    }
    else
    {
        child_name = name;
    }

    // Look for child_name, and make sure it is the only occurrence
    Xml_element child(_e);

    if ( child_name == "." )
        child = _e;

    else if ( child_name == ".." )
        child = get_parent();

    else
        if ( _e->FirstChildElement(child_name.c_str()) == nullptr )
            return false;
        else
            child = _e->FirstChildElement(child_name.c_str());

    if (child._e->NextSiblingElement(child_name.c_str()) != nullptr ) 
    {
        std::ostringstream s_out;
        s_out << "[ERROR] Xml_element::has_child -> There are several children with name \"" + child_name + "\"" << std::endl;
        print(s_out);
        throw lion_exception(s_out.str());
    }

    if ( the_rest.size() == 0 )
        return true;

    else
        return child.has_child(the_rest);
}


inline std::string Xml_element::get_attribute(const std::string& attribute, std::string&&) const 
{ 
    if ( has_attribute(attribute) )
        return _e->Attribute(attribute.c_str()); 
    else
    {
        std::ostringstream s_out;
        s_out << "[ERROR] Xml_element::get_attribute -> Attribute \"" + attribute + "\" was not found" << std::endl;
        print(s_out);
        throw lion_exception(s_out.str());
    }
}


inline bool Xml_element::this_and_childs_have_attributes() const
{
    // check if this element has attributes
    if (_e -> FirstAttribute() != nullptr) {
        return true;
    }
    else {

        // check if children have attributes
        for (const auto& child : get_children()) {
            if (child.this_and_childs_have_attributes()) {
                return true;
            }
        }

    }

    return false;
}


inline bool Xml_element::this_and_childs_have_both_value_and_children() const
{
    if (!(get_value().empty()) && (get_children().size() > 0u)) {
        return true;
    }
    else {

        // check if children have attributes
        for (const auto& child : get_children()) {
            if (child.this_and_childs_have_both_value_and_children()) {
                return true;
            }
        }

    }

    return false;
}

#endif
