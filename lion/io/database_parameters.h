#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include "Xml_document.h"
#include "lion/math/vector3d.hpp"
#include "lion/math/matrix3x3.h"

template<typename T>
struct Database_parameter
{
    static_assert(std::is_same_v<T,void> || std::is_same_v<T,const void>);

    enum Parameter_type { DOUBLE, INT, STD_VECTOR_DOUBLE, VECTOR3, MATRIX3X3 };

    Database_parameter(const std::string& name_, const Parameter_type type_, T* address_) : name(name_), type(type_), address(address_) {}

    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const double&, double&>::type v) : name(n), type(DOUBLE), address(&v) {}
    
    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const int&, int&>::type v): name(n), type(INT), address(&v) {}

    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const std::vector<double>&, std::vector<double>&>::type v): name(n), type(STD_VECTOR_DOUBLE), address(&v) {}

    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const sVector3d&, sVector3d&>::type v): name(n), type(VECTOR3), address(&v) {}
    
    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const sMatrix3x3&, sMatrix3x3&>::type v): name(n), type(MATRIX3X3), address(&v) {}
        
    std::string name;
    Parameter_type type;
    T* address;
};


using Database_parameter_mutable = Database_parameter<void>;
using Database_parameter_const   = Database_parameter<const void>;


#define DECLARE_PARAMS(...) \
    std::vector<Database_parameter_mutable> get_parameters() { return { __VA_ARGS__ }; } \
    std::vector<Database_parameter_const> get_parameters() const { return { __VA_ARGS__ }; }

inline void read_parameters(Xml_document& doc, const std::string& path, const std::vector<Database_parameter_mutable>& p)
{
    for (auto ip = p.cbegin(); ip != p.cend(); ++ip)
    {
        switch (ip->type)
        {
         case(Database_parameter_mutable::DOUBLE): 
            *static_cast<double*>(ip->address) = doc.get_element(path + ip->name).get_value(double());
            break;

         case(Database_parameter_mutable::INT): 
            *static_cast<int*>(ip->address) = doc.get_element(path + ip->name).get_value(int());
            break;

         case(Database_parameter_mutable::STD_VECTOR_DOUBLE): 
            *static_cast<std::vector<double>*>(ip->address) = doc.get_element(path + ip->name).get_value(std::vector<double>());
            break;

         case(Database_parameter_mutable::VECTOR3): 
            *static_cast<sVector3d*>(ip->address) = doc.get_element(path + ip->name).get_value(sVector3d());
            break;

         case(Database_parameter_mutable::MATRIX3X3): 
            *static_cast<sMatrix3x3*>(ip->address) = doc.get_element(path + ip->name).get_value(sMatrix3x3());
            break;

         default:
            break;
        }
    }
}


template<typename T>
inline bool set_parameter(const std::vector<Database_parameter_mutable>& p, const std::string& parameter, const std::string& path, const T value)
{
    for (auto ip = p.cbegin(); ip != p.cend(); ++ip)
    {
        if ( (path+ip->name) != parameter)
            continue;

        switch (ip->type)
        {
         case(Database_parameter_mutable::DOUBLE): 
            if constexpr (std::is_same<double,T>::value)
            {
                *static_cast<double*>(ip->address) = value;
                return true;
            }
            else
            {   
                throw std::runtime_error("Attempt to set double from non-double type");
            }
            break;

         case(Database_parameter_mutable::INT): 
            if constexpr (std::is_same<int,T>::value)
            {
                *static_cast<int*>(ip->address) = value;
                return true;
            }
            else
            {   
                throw std::runtime_error("Attempt to set int from non-int type");
            }
            break;

         case(Database_parameter_mutable::STD_VECTOR_DOUBLE): 
            if constexpr (std::is_same<std::vector<double>,T>::value)
            {
                *static_cast<std::vector<double>*>(ip->address) = value;
                return true;
            }
            else
            {   
                throw std::runtime_error("Attempt to set double vector from non-double vector type");
            }
            break;

         case(Database_parameter_mutable::VECTOR3): 
            if constexpr (std::is_same<sVector3d,T>::value)
            {
                *static_cast<sVector3d*>(ip->address) = value;
                return true;
            }
            else
            {
                throw std::runtime_error("Attempt to set vector3d from non-vector3d type");
            }
            break;

         case(Database_parameter_mutable::MATRIX3X3): 
            if constexpr (std::is_same<sMatrix3x3,T>::value)
            {
                *static_cast<sMatrix3x3*>(ip->address) = value;
                return true;
            }
            else
            {
                throw std::runtime_error("Attempt to set matrix3x3 from non-matrix3x3 type");
            }
            break;

         default:
            break;
        }
    }

    return false;
}

template<typename T>
inline void write_parameters(Xml_document& doc, const std::string& path, const std::vector<Database_parameter<T>>& p)
{
    for (auto ip = p.cbegin(); ip != p.cend(); ++ip)
    {
        std::ostringstream s_out;
        s_out << std::setprecision(17);
        switch (ip->type)
        {
         case(Database_parameter<T>::DOUBLE): 
         {
            s_out << *static_cast<const double*>(ip->address);
            doc.add_element(path + ip->name).set_value(s_out.str());
            break;
         }
         case(Database_parameter<T>::INT): 
         {
            doc.add_element(path + ip->name).set_value(std::to_string(*static_cast<const int*>(ip->address)));
            break;
         }
         case(Database_parameter<T>::STD_VECTOR_DOUBLE): 
         {
            const auto& v = *static_cast<const std::vector<double>*>(ip->address);
            
            for (auto it = v.cbegin(); it != v.cend()-1; ++it )
                s_out << *it << ", " ;
            s_out << v.back();
            doc.add_element(path + ip->name).set_value(s_out.str());
            break;
         }
         case(Database_parameter<T>::VECTOR3): 
         {
            const auto& v = *static_cast<const sVector3d*>(ip->address);
            s_out << v.x() << ", " << v.y() << ", " << v.z() ;
            doc.add_element(path + ip->name).set_value(s_out.str());
            break;
         }
         case(Database_parameter<T>::MATRIX3X3): 
         {
            const auto& m = *static_cast<const sMatrix3x3*>(ip->address);
            s_out << m.xx() << ", " << m.xy() << ", " << m.xz() << ", "; 
            s_out << m.yx() << ", " << m.yy() << ", " << m.yz() << ", "; 
            s_out << m.zx() << ", " << m.zy() << ", " << m.zz() ; 
            doc.add_element(path + ip->name).set_value(s_out.str());
            break;
         }
         default:
            break;
        }
    }
}

#endif
