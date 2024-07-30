#ifndef LION_IO_DOCUMENT_H
#define LION_IO_DOCUMENT_H

#include "document_element.h"

class Document
{
public:
    Document(const std::string& name = "") : _name(name) {}

    virtual void load() = 0;

    virtual bool save() = 0;

    bool save(const std::string& name)
    {
        _name = name;
        return save();
    }

    const std::string& get_file_name() const { return _name; }

    virtual void parse(const char* contents) = 0;

    virtual void parse(const std::string& contents) = 0;

    virtual Document_element::value_ptr get_root_element() = 0;

    virtual Document_element::value_ptr get_element(const std::string& name) = 0;


private:

    std::string _name;
};


using Document_ptr = std::shared_ptr<Document>;

#endif
