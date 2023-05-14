#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include "cppad/cppad.hpp"
#include "Xml_document.h"
#include "lion/math/vector3d.hpp"
#include "lion/math/matrix3x3.h"

template<typename T>
struct Database_parameter
{
    static_assert(std::is_same_v<T,void> || std::is_same_v<T,const void>);

    enum Parameter_type { DOUBLE, INT, STD_VECTOR_DOUBLE, VECTOR3, MATRIX3X3, AD, BOOL };

    Database_parameter(const std::string& name_, const Parameter_type type_, T* address_) : name(name_), type(type_), address(address_) {}

    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const scalar&, scalar&>::type v) : name(n), type(DOUBLE), address(&v) {}
    
    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const int&, int&>::type v): name(n), type(INT), address(&v) {}

    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const std::vector<scalar>&, std::vector<scalar>&>::type v): name(n), type(STD_VECTOR_DOUBLE), address(&v) {}

    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const sVector3d&, sVector3d&>::type v): name(n), type(VECTOR3), address(&v) {}
    
    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const sMatrix3x3&, sMatrix3x3&>::type v): name(n), type(MATRIX3X3), address(&v) {}

    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const CppAD::AD<scalar>&, CppAD::AD<scalar>&>::type v)
    : name(n), type(AD), address(&v) {}

    Database_parameter(const std::string& n, typename std::conditional<std::is_const<T>::value, const bool&, bool&>::type v): name(n), type(BOOL), address(&v) {}
        
    std::string name;
    Parameter_type type;
    T* address;
};


using Database_parameter_mutable = Database_parameter<void>;
using Database_parameter_const   = Database_parameter<const void>;

#define DECLARE_PARAMS(...) \
    std::vector<Database_parameter_mutable> __get_parameters() { return { __VA_ARGS__ }; } \
    std::vector<Database_parameter_const> __get_parameters() const { return { __VA_ARGS__ }; } \
    std::vector<bool> __used_parameters = std::vector<bool>(__get_parameters().size(),false); 

inline void read_parameters(Xml_document& doc, const std::string& path, const std::vector<Database_parameter_mutable>& p, std::vector<bool>& used_parameters)
{
    assert(p.size() == used_parameters.size());

    for (size_t i = 0; i < p.size(); ++i)
    {
        Xml_element element = doc.get_element(path + p[i].name);

        switch (p[i].type)
        {
         case(Database_parameter_mutable::DOUBLE): 
            *static_cast<scalar*>(p[i].address) = element.get_value(scalar());
            break;

         case(Database_parameter_mutable::INT): 
            *static_cast<int*>(p[i].address) = element.get_value(int());
            break;

         case(Database_parameter_mutable::STD_VECTOR_DOUBLE): 
            *static_cast<std::vector<scalar>*>(p[i].address) = element.get_value(std::vector<scalar>());
            break;

         case(Database_parameter_mutable::VECTOR3): 
            *static_cast<sVector3d*>(p[i].address) = element.get_value(sVector3d());
            break;

         case(Database_parameter_mutable::MATRIX3X3): 
            *static_cast<sMatrix3x3*>(p[i].address) = element.get_value(sMatrix3x3());
            break;

         case(Database_parameter_mutable::AD): 
            *static_cast<CppAD::AD<scalar>*>(p[i].address) = element.get_value(scalar());
            break;

         case(Database_parameter_mutable::BOOL): 
            *static_cast<bool*>(p[i].address) = element.get_value(bool());
            break;

         default:
            throw lion_exception("[ERROR] read_parameters -> parameter type was not recognized");
            break;
        }
        used_parameters[i] = true;
        element.set_attribute<bool>("__unused__", false);

        while(element.has_parent())
        {
            element = element.get_parent();

            if(element.has_attribute("__unused__"))
                break;
            else
                element.set_attribute<bool>("__unused__", false);
        }
    }
}


template<typename T>
inline bool set_parameter(const std::vector<Database_parameter_mutable>& p, std::vector<bool>& used_parameters, const std::string& parameter, const std::string& path, const T value)
{
    assert(p.size() == used_parameters.size());
    for (size_t i = 0; i < p.size(); ++i)
    {
        if ( (path+p[i].name) != parameter)
            continue;

        used_parameters[i] = true;
        switch (p[i].type)
        {
         case(Database_parameter_mutable::DOUBLE): 
            if constexpr (std::is_same<scalar,T>::value)
            {
                *static_cast<scalar*>(p[i].address) = value;
                return true;
            }
            else if constexpr (std::is_same<CppAD::AD<scalar>,T>::value)
            {
                *static_cast<scalar*>(p[i].address) = Value(CppAD::Var2Par(value));
                return true;
            }
            else
            {   
                throw lion_exception("Attempt to set scalar from non-scalar type");
            }
            break;

         case(Database_parameter_mutable::INT): 
            if constexpr (std::is_same<int,T>::value)
            {
                *static_cast<int*>(p[i].address) = value;
                return true;
            }
            else
            {   
                throw lion_exception("Attempt to set int from non-int type");
            }
            break;

         case(Database_parameter_mutable::STD_VECTOR_DOUBLE): 
            if constexpr (std::is_same<std::vector<scalar>,T>::value)
            {
                *static_cast<std::vector<scalar>*>(p[i].address) = value;
                return true;
            }
            else
            {   
                throw lion_exception("Attempt to set scalar vector from non-scalar vector type");
            }
            break;

         case(Database_parameter_mutable::VECTOR3): 
            if constexpr (std::is_same<sVector3d,T>::value)
            {
                *static_cast<sVector3d*>(p[i].address) = value;
                return true;
            }
            else
            {
                throw lion_exception("Attempt to set vector3d from non-vector3d type");
            }
            break;

         case(Database_parameter_mutable::MATRIX3X3): 
            if constexpr (std::is_same<sMatrix3x3,T>::value)
            {
                *static_cast<sMatrix3x3*>(p[i].address) = value;
                return true;
            }
            else
            {
                throw lion_exception("Attempt to set matrix3x3 from non-matrix3x3 type");
            }
            break;

         case(Database_parameter_mutable::AD): 
            if constexpr (std::is_same<scalar,T>::value || std::is_same<CppAD::AD<scalar>,T>::value)
            {
                *static_cast<CppAD::AD<scalar>*>(p[i].address) = value;
                return true;
            }
            else
            {
                throw lion_exception("Attempt to set AD variable from non-scalar or CppAD::AD<scalar> types");
            }
            break;

         case(Database_parameter_mutable::BOOL): 
            if constexpr (std::is_same<bool,T>::value)
            {
                *static_cast<bool*>(p[i].address) = value;
                return true;
            }
            else
            {
                throw lion_exception("Attempt to set BOOL variable from non-bool type");
            }
            break;

         default:
            throw lion_exception("[ERROR] set_parameter -> type not recognized");
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
            s_out << *static_cast<const scalar*>(ip->address);
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
            const auto& v = *static_cast<const std::vector<scalar>*>(ip->address);
            
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
         case(Database_parameter<T>::AD): 
         {
            s_out << Value(*static_cast<const CppAD::AD<scalar>*>(ip->address));
            doc.add_element(path + ip->name).set_value(s_out.str());
            break;
         }
         case(Database_parameter<T>::BOOL): 
         {
            doc.add_element(path + ip->name).set_value(std::to_string(*static_cast<const bool*>(ip->address)));
            break;
         }

         default:
            break;
        }
    }
}


inline bool database_parameters_all_used(const Xml_element& element)
{
    if( !element.has_attribute("__unused__") )
        return false;

    if ( !element.has_children() )
        return true;

    for (auto child : element.get_children() )
        if ( !database_parameters_all_used(child) )
            return false;

    return true;
}

#endif
