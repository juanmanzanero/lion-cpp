#ifndef __XML_ELEMENT_HPP__
#define __XML_ELEMENT_HPP__


inline Xml_element::Xml_element(tinyxml2::XMLElement* e) : _e(e) 
{
    if ( e == nullptr )
        throw std::runtime_error("Attempted to construct Xml_element from nullptr"); 
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
 try
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
            throw std::runtime_error("No child found with name \"" + child_name + "\"");

        child = _e->FirstChildElement(child_name.c_str());
    }

    if (_e->NextSiblingElement(child_name.c_str()) != nullptr ) 
        throw std::runtime_error("There are several children with name \"" + child_name + "\"");

    if ( the_rest.size() == 0 )
        return child;

    else
        return child.get_child(the_rest);
 }
 catch (std::runtime_error& err)
 {
    out(2) << "Full name: " << name << std::endl;
    throw;
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

    if (_e->NextSiblingElement(child_name.c_str()) != nullptr ) 
        throw std::runtime_error("There are several children with name \"" + child_name + "\"");

    if ( the_rest.size() == 0 )
        return true;

    else
        return child.has_child(the_rest);
}



#endif
