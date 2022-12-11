#ifndef __XML_ELEMENT_H__
#define __XML_ELEMENT_H__

#include "lion/thirdparty/include/tinyxml2.h"
#include "lion/foundation/utils.hpp"
#include "lion/math/vector3d.hpp"
#include "lion/math/matrix3x3.h"
#include "lion/math/matrix_extensions.h"
#include "lion/thirdparty/include/logger.hpp"

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
    std::string get_value()              { return (_e->GetText() == nullptr) ? "" : _e->GetText(); }
    std::string get_value(std::string&&) { return get_value(); }

    //! Get value as double
    //! @param[in] simply pass int() to overload this version
    int get_value(int&&) { return std::stoi(get_value()); }

    //! Get value as double
    //! @param[in] simply pass double() to overload this version
    double get_value(double&&) { return std::stod(get_value()); }

    //! Get value as double vector
    //! @param[in] simply pass std::vector<double>() to overload this version
    std::vector<double> get_value(std::vector<double>&&) { return string_to_double_vector(get_value()); }

    //! Get value as vector3d
    //! @param[in] simply pass Vector3d() to overload this function
    sVector3d get_value(sVector3d&&) { return sVector3d(string_to_double_vector(get_value())); }

    //! Get value as matrix3x3
    //! @param[in] simply pass Matrix3x3() to overload this function
    sMatrix3x3 get_value(sMatrix3x3&&) { return transpose(sMatrix3x3(string_to_double_vector(get_value()))); }

    //! Get value as bool
    //! @param[in] simply pass bool() to overload this function
    bool get_value(bool&&) { return to_bool(get_value()); }

    //! Set the value from string
    //! @param[in] val: new string value
    Xml_element& set_value(const std::string& val) { _e->SetText(val.c_str()); return *this;}

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
    std::string get_attribute(const std::string& attribute) const;

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


#include "Xml_element.hpp"

#endif
