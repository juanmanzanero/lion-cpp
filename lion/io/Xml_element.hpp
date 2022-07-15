#ifndef __XML_ELEMENT_HPP__
#define __XML_ELEMENT_HPP__


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


#endif
