#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include "Xml_document.h"
#include "lion/math/vector3d.h"
#include "lion/math/matrix3x3.h"

struct Database_parameter
{
    Database_parameter(const std::string& n, double& v) : name(n), type(DOUBLE), address(&v) {}
    
    Database_parameter(const std::string& n, int& v): name(n), type(INT), address(&v) {}

    Database_parameter(const std::string& n, std::vector<double>& v): name(n), type(STD_VECTOR_DOUBLE), address(&v) {}

    Database_parameter(const std::string& n, sVector3d& v): name(n), type(VECTOR3), address(&v) {}
    
    Database_parameter(const std::string& n, sMatrix3x3& v): name(n), type(MATRIX3X3), address(&v) {}
        
    enum Parameter_type { DOUBLE, INT, STD_VECTOR_DOUBLE, VECTOR3, MATRIX3X3 };
    std::string name;
    Parameter_type type;
    void* address;
};


inline void read_parameters(Xml_document& doc, const std::string& path, const std::vector<Database_parameter>& p)
{
    for (auto ip = p.cbegin(); ip != p.cend(); ++ip)
    {
        switch (ip->type)
        {
         case(Database_parameter::DOUBLE): 
            *static_cast<double*>(ip->address) = doc.get_element(path + ip->name).get_value(double());
            break;

         case(Database_parameter::INT): 
            *static_cast<int*>(ip->address) = doc.get_element(path + ip->name).get_value(int());
            break;

         case(Database_parameter::STD_VECTOR_DOUBLE): 
            *static_cast<std::vector<double>*>(ip->address) = doc.get_element(path + ip->name).get_value(std::vector<double>());
            break;

         case(Database_parameter::VECTOR3): 
            *static_cast<sVector3d*>(ip->address) = doc.get_element(path + ip->name).get_value(sVector3d());
            break;

         case(Database_parameter::MATRIX3X3): 
            *static_cast<sMatrix3x3*>(ip->address) = doc.get_element(path + ip->name).get_value(sMatrix3x3());
            break;

         default:
            break;
        }
    }
}


#endif
