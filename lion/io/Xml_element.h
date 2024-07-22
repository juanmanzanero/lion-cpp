#ifndef __XML_ELEMENT_H__
#define __XML_ELEMENT_H__

#include "tinyxml2.h"

#include "lion/thirdparty/include/logger.hpp"
#include "lion/foundation/utils.hpp"
#include "lion/math/vector3d.hpp"
#include "lion/math/matrix3x3.h"
#include "lion/math/matrix_extensions.h"
#include "document_element.h"

class Xml_element : public Document_element
{
 public:

    using base_type = Document_element;
    using base_type::get_value;
    using base_type::value_ptr;

    //! Constructor
    //! @param[in] e: pointer to a tinyxml2 XMLElement
    Xml_element(tinyxml2::XMLElement* e) : base_type(e, e->Name()) {}


    // cast the element void* in the parent class to the actual type for its usage
    // in this implementation
    tinyxml2::XMLElement& e_xml() { return *static_cast<tinyxml2::XMLElement*>(this->e()); }
    const tinyxml2::XMLElement& e_xml_const() const { return *static_cast<const tinyxml2::XMLElement*>(this->e()); }
    tinyxml2::XMLElement* e_xml_ptr() { return static_cast<tinyxml2::XMLElement*>(this->e()); }

    //! Set the element name
    //! @param[in] name: new name for the element
    void set_name(const std::string& name) { e_xml().SetName(name.c_str()); }

    //! Get value as string
    virtual std::string get_value() const override { return (e_xml_const().GetText() == nullptr) ? "" : e_xml_const().GetText(); }

    //! Get value as double
    //! @param[in] simply pass int() to overload this version
    virtual int get_value(int&&) const override { return std::stoi(get_value()); }

    //! Get value as double
    //! @param[in] simply pass double() to overload this version
    virtual double get_value(double&&) const override { return std::stod(get_value()); }

    //! Get value as double vector
    //! @param[in] simply pass std::vector<double>() to overload this version
    virtual std::vector<double> get_value(std::vector<double>&&) const override { return string_to_double_vector<double>(get_value()); }

    //! Get value as float vector
    //! @param[in] simply pass std::vector<float>() to overload this version
    virtual std::vector<float> get_value(std::vector<float>&&) const override { return string_to_double_vector<float>(get_value()); }

    //! Get value as vector3d
    //! @param[in] simply pass Vector3d() to overload this function
    virtual sVector3d get_value(sVector3d&&) const override { return sVector3d(string_to_double_vector<double>(get_value())); }

    //! Get value as matrix3x3
    //! @param[in] simply pass Matrix3x3() to overload this function
    virtual sMatrix3x3 get_value(sMatrix3x3&&) const override { return transpose(sMatrix3x3(string_to_double_vector<double>(get_value()))); }


    //! Set the value from different variable types
    virtual void set_value(const std::string& val) override { e_xml().SetText(val.c_str()); }

    virtual void set_value(const char* val) override { e_xml().SetText(val); }

    virtual void set_value(const double& val) override { set_value_generic(val); }

    virtual void set_value(const int& val) override { set_value_generic(val); }

    virtual void set_value(const float& val) override { set_value_generic(val); }

    virtual void set_value(const std::vector<double>& val) override { set_value_generic(val); }

    virtual void set_value(const std::vector<int>& val) override { set_value_generic(val); }

    virtual void set_value(const std::vector<float>& val) override { set_value_generic(val); }

    virtual void set_value(const sVector3d& val) override { set_value_generic(val); }

    virtual void set_value(const sMatrix3x3& val) override { set_value_generic(val); }


    //! Set the value from numbers (scalars, std::vector, std::array)
    template<typename T,
             typename std::enable_if_t<std::is_arithmetic_v<T> >* = nullptr>
    void set_value_generic(const T &num) { return set_value(num2str(num)); }

    template<typename T>
    void set_value_generic(const std::vector<T> &vec) { return set_value(vec2str(vec)); }

    template<typename T, std::size_t N>
    void set_value_generic(const std::array<T, N> &arr) { return set_value(vec2str(arr)); }



    //! Get all children to a std::vector
    virtual std::vector<value_ptr> get_children() override;

    //! Add children
    virtual value_ptr add_child(const std::string& name) override; 

    //! Add comment
    void add_comment(const std::string& text) { e_xml().InsertNewComment(text.c_str()); }

    //! See if the element has children
    bool has_children() const override {return (e_xml_const().FirstChildElement() ? true : false);}

    //! Get the first child
    Xml_element get_first_child()  {return e_xml().FirstChildElement();}

    //! Get child by name
    virtual value_ptr get_child(const std::string& name) override;

    //! See if a child exists
    virtual bool has_child(const std::string& name) override;


    //! Try to get the value of a child, or return a default one
    template<typename ValueType>
    ValueType try_get_child_value(const std::string &childname, ValueType default_value,
                                  bool warn_if_returning_default_value = true)
    {
        if (has_child(childname)) {
            return get_child(childname)->get_value(ValueType{});
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
    bool has_parent() const { return (e_xml_const().Parent()->ToElement() ? true : false); }

    //! Get parent
    Xml_element get_parent()  { return e_xml().Parent()->ToElement(); }

    //! See if sibling exists
    bool has_sibling() const { return (e_xml_const().NextSiblingElement() ? true : false); }

    //! Get Next Sibling
    Xml_element get_sibling() { return e_xml().NextSiblingElement()->ToElement(); }

    //! See if an attribute exists
    //! @param[in] attribute: name of the attribute
    bool has_attribute(const std::string& attribute) const {return (e_xml_const().FindAttribute(attribute.c_str()) ? true : false);}

    //! Get attribute by name
    //! @param[in] attribute: name of the attribute
    std::string get_attribute(const std::string& attribute, std::string&& = "") const;

    //! Get attribute as double
    //! @param[in] attribute: name of the attribute
    //! @param[in] simply pass double() to overload this version
    double get_attribute(const std::string& attribute, double&&) const { return std::stod(e_xml_const().Attribute(attribute.c_str())); }

    //! Sets a new attribute or modifies the value of an existing one
    void set_attribute(const std::string& attrib_name, const std::string& attrib_value) { e_xml().SetAttribute(attrib_name.c_str(), attrib_value.c_str()); }

    template<typename T>
    void set_attribute(const std::string& attrib_name, const T& attrib_value) { e_xml().SetAttribute(attrib_name.c_str(), attrib_value); }

    //! Print
    virtual void print(std::ostream& os) const override;
    
    //! Copy the contents of other into this node
    virtual void copy_contents(Document_element& other) override;

    //! Delete attribute
    void delete_attribute(const std::string& attribute_name)
    {

        if (has_attribute(attribute_name))
            e_xml().DeleteAttribute(attribute_name.c_str());
        else
            throw lion_exception("[ERROR] delete_attribute -> attribute does not exist");
    } 


    //! Check if this element and its childs have any attribute
    bool this_and_childs_have_attributes();


    //! Check if this element and its childs have value and children
    bool this_and_childs_have_both_value_and_children();


    inline Document_element::value_ptr to_value_ptr() { return Document_element_ptr(std::make_shared<Xml_element>(e_xml_ptr())); }
};


inline void Xml_element::print(std::ostream& os) const
{
    tinyxml2::XMLPrinter printer;
    e_xml_const().Accept(&printer);
    os << printer.CStr();
}


inline std::ostream& operator<<(std::ostream& os, const Xml_element& e)
{
    e.print(os);
    return os;
}


inline void Xml_element::copy_contents(Document_element& other_)
{
    auto& other = static_cast<Xml_element&>(other_);
    tinyxml2::XMLElement* element = other.e_xml().FirstChildElement(); 
    tinyxml2::XMLNode* previous_sibling_ptr = e_xml().LastChild();
    while (element != nullptr)
    {
        auto* element_copy = element->DeepClone(e_xml().GetDocument()); 

        if (previous_sibling_ptr != nullptr)
            e_xml().InsertAfterChild(previous_sibling_ptr, element_copy);
        else
            e_xml().InsertFirstChild(element_copy);

        element = element->NextSiblingElement();
        previous_sibling_ptr = element_copy;
    }
}


inline std::vector<Document_element::value_ptr> Xml_element::get_children()
{
    std::vector<value_ptr> output;

    tinyxml2::XMLElement* element = e_xml().FirstChildElement();

    while (element != nullptr)
    {
        output.push_back(std::make_shared<Xml_element>(element));
        element = element->NextSiblingElement();
    }

    return output;
}


inline Document_element::value_ptr Xml_element::get_child(const std::string& name)
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
    Xml_element child(&e_xml());

    if ( child_name == "." )
        child = &e_xml();

    else if ( child_name == ".." )
        child = get_parent();

    else
    {
        if ( e_xml().FirstChildElement(child_name.c_str()) == nullptr )
        {
            std::ostringstream s_out;
            s_out << "[ERROR] Xml_element::get_child -> No child found with name \"" + child_name + "\"" << std::endl;
            print(s_out);
            throw lion_exception(s_out.str());
        }

        child = e_xml().FirstChildElement(child_name.c_str());
    }

    if (child.e_xml().NextSiblingElement(child_name.c_str()) != nullptr ) 
    {
        std::ostringstream s_out;
        s_out << "[ERROR] Xml_element::get_child -> There are several children with name \"" + child_name + "\"" << std::endl;
        print(s_out);
        throw lion_exception("There are several children with name \"" + child_name + "\"");
    }

    if ( the_rest.size() == 0 )
        return std::make_shared<Xml_element>(child);

    else
        return child.get_child(the_rest);
}


inline Document_element::value_ptr Xml_element::add_child(const std::string& name) 
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
    Xml_element child(&e_xml());

    if ( the_rest.size() > 0 )
    {
        // The child exists and its me
        if ( child_name == "." )
            child = &e_xml();

        // The child exists and its my parent
        else if ( child_name == ".." )
            child = get_parent();

        else
        {
            // The children does not exist -> create it
            if ( e_xml().FirstChildElement(child_name.c_str()) == nullptr )
            {
                child = e_xml().InsertNewChildElement(child_name.c_str());
            }
            else
            {
                // The children does exist -> point to it
                child = e_xml().FirstChildElement(child_name.c_str());
        
                // Check uniqueness
                if (child.e_xml().NextSiblingElement(child_name.c_str()) != nullptr ) 
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
        if ( e_xml().FirstChildElement(child_name.c_str()) != nullptr )
        {
            std::ostringstream s_out;
            s_out << "[ERROR] Xml_element::add_child -> Node \"" + child_name + "\" already exists" << std::endl;
            print(s_out);
            throw lion_exception("Node \"" + child_name + "\" already exists");
        }

        return std::make_shared<Xml_element>(e_xml().InsertNewChildElement(name.c_str()));
    }
}


inline bool Xml_element::has_child(const std::string& name)
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
    Xml_element child(&e_xml());

    if ( child_name == "." )
        child = &e_xml();

    else if ( child_name == ".." )
        child = get_parent();

    else
        if ( e_xml().FirstChildElement(child_name.c_str()) == nullptr )
            return false;
        else
            child = e_xml().FirstChildElement(child_name.c_str());

    if (child.e_xml().NextSiblingElement(child_name.c_str()) != nullptr ) 
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
        return e_xml_const().Attribute(attribute.c_str()); 
    else
    {
        std::ostringstream s_out;
        s_out << "[ERROR] Xml_element::get_attribute -> Attribute \"" + attribute + "\" was not found" << std::endl;
        print(s_out);
        throw lion_exception(s_out.str());
    }
}


inline bool Xml_element::this_and_childs_have_attributes()
{
    // check if this element has attributes
    if (e_xml_const().FirstAttribute() != nullptr) {
        return true;
    }
    else {

        // check if children have attributes
        for (auto& child : get_children()) {
            if (child.cast<Xml_element>().this_and_childs_have_attributes()) {
            //if (static_cast<Xml_element&>(*child).this_and_childs_have_attributes()) {
                return true;
            }
        }

    }

    return false;
}


inline bool Xml_element::this_and_childs_have_both_value_and_children()
{
    if (!(get_value().empty()) && (get_children().size() > 0u)) {
        return true;
    }
    else {

        // check if children have attributes
        for (auto& child : get_children()) {
            if (child.cast<Xml_element>().this_and_childs_have_both_value_and_children()) {
                return true;
            }
        }
    }

    return false;
}

#endif
